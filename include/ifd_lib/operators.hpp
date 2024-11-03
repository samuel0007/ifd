#ifndef IFD_OPERATORS_HPP
#define IFD_OPERATORS_HPP

#include <algorithm>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace operators {

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triplet;


template <typename OperatorT, typename SolverT>
class EvolutionProblem {
public:
    EvolutionProblem(std::function<double(std::array<double, 3>)> IC, double dt, std::shared_ptr<mesh::Mesh> mesh_p): _dt(dt), _mesh_p(mesh_p) {
        this->_N = _mesh_p->cells.size();
        x.resize(_N);
        const std::vector<std::array<double, 3>>& cells_centroid = mesh_p->cells.getCellsCentroid();
        for(size_t i = 0; i < _N; ++i) {
            x[i] = IC(cells_centroid[i]);
        }

        _A.resize(_N, _N);
        _A.setZero();
        _b.resize(_N);
    }

    // Assemble Jacobian and factorize matrix
    void initialize() {
        _assemble_jacobian(this->_A);
        solver.compute(this->_A);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Matrix is not PD" << std::endl;
        }
    }

    void evolve() {
        _assemble_RHS(_b);

        // Solve Ax = b;
        x = solver.solve(_b);
    }

    double getDt() {
        return this-> _dt;
    }

    Eigen::VectorXd x;

protected:
    virtual void _assemble_jacobian(OperatorT& A) = 0;
    virtual void _assemble_RHS(Eigen::VectorXd& b) = 0;

    double _dt;
    size_t _N;
    std::shared_ptr<mesh::Mesh> _mesh_p;
    OperatorT _A;
    Eigen::VectorXd _b;
    SolverT solver; // Eigen Solver
};


// For now just transient laplacian, will be abstracted in the future
template <typename SolverT>
class DenseTransientLaplacianProblem: public EvolutionProblem<Eigen::MatrixXd, SolverT> {
public:
    DenseTransientLaplacianProblem(double gamma, std::function<double(std::array<double, 3>)> IC, double dt, std::shared_ptr<mesh::Mesh> mesh_p): _gamma(gamma), EvolutionProblem<Eigen::MatrixXd, SolverT>(IC, dt, mesh_p){
        std::vector<double> gamma_field(this->_N);
        std::ranges::fill(gamma_field, gamma);
        mesh_p->addCellData(gamma_field, "gamma");        
    };

private:
    double _gamma;
    
    void _assemble_ddt_A(Eigen::MatrixXd& A) {
        const std::vector<double>& cells_volume = this->_mesh_p->cells.getCellsVolume();
        for(size_t cell_idx = 0; cell_idx < this->_N; ++cell_idx) {
            A(cell_idx, cell_idx) += cells_volume[cell_idx] / this->_dt;
        }
    }

    void _assemble_ddt_b(Eigen::VectorXd& b) {
        const std::vector<double>& cells_volume = this->_mesh_p->cells.getCellsVolume();
        for(size_t cell_idx = 0; cell_idx < this->_N; ++cell_idx) {
            b(cell_idx) = cells_volume[cell_idx] * this->x[cell_idx] /this->_dt;
        }
    }

    void _assemble_laplacian(Eigen::MatrixXd& A) {
        // a_n
        const std::vector<std::vector<size_t>>& faces_cellsID = this->_mesh_p->getFacesCellsID();
        const std::vector<bool>& faces_boundaryFlag = this->_mesh_p->getFacesBoundaryFlag();
        const std::vector<double>& faces_delta = this->_mesh_p->getFacesDelta();
        const std::vector<double>& faces_gamma = this->_mesh_p->getFaceData("gamma");
        const std::vector<std::array<double, 3>>& faces_normal_vectors = this->_mesh_p->faces.getFacesNormalVector();

        for(const auto& [cellsID, gamma, normal_vector, delta, bd_flag]: std::views::zip(faces_cellsID, faces_gamma, faces_normal_vectors, faces_delta, faces_boundaryFlag)) {
            if(bd_flag) {
                continue;
            }

            double mag_normal_vector = std::sqrt(normal_vector[0]*normal_vector[0] + normal_vector[1]*normal_vector[1] + normal_vector[2] * normal_vector[2]);
            double a_n = - mag_normal_vector * delta * gamma;

            A(cellsID[0], cellsID[1]) += a_n;
            A(cellsID[1], cellsID[0]) += a_n;
        }

        // a_p
        for(size_t cell_idx = 0; cell_idx < this->_N; ++cell_idx) {
            A(cell_idx, cell_idx) = -A.row(cell_idx).sum();
        }
    }
    
    void _assemble_jacobian(Eigen::MatrixXd& A) {
        _assemble_laplacian(A);
        _assemble_ddt_A(A);
    }

    void _assemble_RHS(Eigen::VectorXd& b) {
        _assemble_ddt_b(b);
    }
};

// TESTED with: Eigen::SimplicialCholesky<SpMat, Eigen::Upper> solver;
template <typename SolverT>
class SparseTransientLaplacianProblem: public EvolutionProblem<SpMat, SolverT> {
public:
    SparseTransientLaplacianProblem(double gamma, std::function<double(std::array<double, 3>)> IC, double dt, std::shared_ptr<mesh::Mesh> mesh_p): _gamma(gamma), EvolutionProblem<SpMat, SolverT>(IC, dt, mesh_p){
        std::vector<double> gamma_field(this->_N);
        std::ranges::fill(gamma_field, gamma);
        mesh_p->addCellData(gamma_field, "gamma");
    };

private:
    double _gamma;
    std::vector<Triplet> _triplets;

    void _assemble_ddt_A(std::vector<Triplet>& triplets) {
        const std::vector<double>& cells_volume = this->_mesh_p->cells.getCellsVolume();
        for(size_t cell_idx = 0; cell_idx < this->_N; ++cell_idx) {
            triplets.push_back(Triplet(cell_idx, cell_idx, cells_volume[cell_idx] / this->_dt));
        }
    }

    void _assemble_ddt_b(Eigen::VectorXd& b) {
        const std::vector<double>& cells_volume = this->_mesh_p->cells.getCellsVolume();
        for(size_t cell_idx = 0; cell_idx < this->_N; ++cell_idx) {
            b(cell_idx) = cells_volume[cell_idx] * this->x[cell_idx] /this-> _dt;
        }
    }

    void _assemble_laplacian(std::vector<Triplet>& triplets) {
        const std::vector<std::vector<size_t>>& faces_cellsID = this->_mesh_p->getFacesCellsID();
        const std::vector<bool>& faces_boundaryFlag = this->_mesh_p->getFacesBoundaryFlag();
        const std::vector<double>& faces_delta = this->_mesh_p->getFacesDelta();
        const std::vector<double>& faces_gamma = this->_mesh_p->getFaceData("gamma");
        const std::vector<std::array<double, 3>>& faces_normal_vectors = this->_mesh_p->faces.getFacesNormalVector();

        std::vector<double> APs(this->_N);
        std::ranges::fill(APs, 0.);

        for(const auto& [cellsID, gamma, normal_vector, delta, bd_flag]: std::views::zip(faces_cellsID, faces_gamma, faces_normal_vectors, faces_delta, faces_boundaryFlag)) {
            if(bd_flag) {
                continue; // TODO?!
            }

            double mag_normal_vector = std::sqrt(normal_vector[0]*normal_vector[0] + normal_vector[1]*normal_vector[1] + normal_vector[2] * normal_vector[2]);
            double a_n = -mag_normal_vector * delta * gamma;

            // a_n
            // We only store upper diag part
            triplets.push_back(Triplet(cellsID[0], cellsID[1], a_n));

            // a_p
            // Both left and right contribution
            triplets.push_back(Triplet(cellsID[0], cellsID[0], -a_n));
            triplets.push_back(Triplet(cellsID[1], cellsID[1], -a_n));
        }
    }

    void _assemble_jacobian(SpMat& A) {
        // _triplets.reserve(some_estimation)
        std::cout << "Assembling laplacian... \n";
        _assemble_laplacian(_triplets);
        std::cout << "Assembling ddt... \n";
        _assemble_ddt_A(_triplets);
        A.setFromTriplets(_triplets.begin(), _triplets.end());
        // Free memory
        _triplets.resize(0);
        std::cout << "--- Sparse Jacobian assembled --- \n";
    }

    void _assemble_RHS(Eigen::VectorXd& b) {
        _assemble_ddt_b(b);
    }
};

// TESTED with: 
template <typename SolverT>
class SparseConvectionDiffusionProblem: public EvolutionProblem<SpMat, SolverT> {
public:
    SparseConvectionDiffusionProblem(double gamma, std::function<double(std::array<double, 3>)> Ux, std::function<double(std::array<double, 3>)> Uy, std::function<double(std::array<double, 3>)> Uz, std::function<double(std::array<double, 3>)> IC, double dt, std::shared_ptr<mesh::Mesh> mesh_p): _gamma(gamma), EvolutionProblem<SpMat, SolverT>(IC, dt, mesh_p){
        std::vector<double> gamma_field(this->_N);
        std::ranges::fill(gamma_field, gamma);
        mesh_p->addCellData(gamma_field, "gamma");

        const std::vector<std::array<double, 3>>& cells_centroid = mesh_p->cells.getCellsCentroid();

        std::vector<double> ux_field(this->_N);
        std::vector<double> uy_field(this->_N);
        std::vector<double> uz_field(this->_N);
        
        for(size_t i = 0; i < this->_N; ++i) {
            ux_field[i] = Ux(cells_centroid[i]);
            uy_field[i] = Uy(cells_centroid[i]);
            uy_field[i] = Uz(cells_centroid[i]);
        }

        mesh_p->addCellData(ux_field, "ux");
        mesh_p->addCellData(uy_field, "uy");
        mesh_p->addCellData(uz_field, "uz");
    };

private:
    double _gamma;
    std::vector<Triplet> _triplets;
    
    void _assemble_ddt_A(std::vector<Triplet>& triplets) {
        const std::vector<double>& cells_volume = this->_mesh_p->cells.getCellsVolume();
        for(size_t cell_idx = 0; cell_idx < this->_N; ++cell_idx) {
            triplets.push_back(Triplet(cell_idx, cell_idx, cells_volume[cell_idx] / this->_dt));
        }
    }

    void _assemble_ddt_b(Eigen::VectorXd& b) {
        const std::vector<double>& cells_volume = this->_mesh_p->cells.getCellsVolume();
        for(size_t cell_idx = 0; cell_idx < this->_N; ++cell_idx) {
            b(cell_idx) = cells_volume[cell_idx] * this->x[cell_idx] / this->_dt;
        }
    }

    void _assemble_laplacian(std::vector<Triplet>& triplets) {
        const std::vector<std::vector<size_t>>& faces_cellsID = this->_mesh_p->getFacesCellsID();
        const std::vector<bool>& faces_boundaryFlag = this->_mesh_p->getFacesBoundaryFlag();
        const std::vector<double>& faces_delta = this->_mesh_p->getFacesDelta();
        const std::vector<double>& faces_gamma = this->_mesh_p->getFaceData("gamma");
        const std::vector<std::array<double, 3>>& faces_normal_vectors = this->_mesh_p->faces.getFacesNormalVector();

        for(const auto& [cellsID, gamma, normal_vector, delta, bd_flag]: std::views::zip(faces_cellsID, faces_gamma, faces_normal_vectors, faces_delta, faces_boundaryFlag)) {
            double mag_normal_vector = std::sqrt(normal_vector[0]*normal_vector[0] + normal_vector[1]*normal_vector[1] + normal_vector[2] * normal_vector[2]);
            double a_n = -mag_normal_vector * delta * gamma;

            // a_n
            // We only store upper diag part
            if(bd_flag) {
                continue; // TODO?!
            }

            triplets.push_back(Triplet(cellsID[0], cellsID[1], a_n));
            triplets.push_back(Triplet(cellsID[1], cellsID[0], a_n));


            // a_p
            // Both left and right contribution
            triplets.push_back(Triplet(cellsID[0], cellsID[0], -a_n));
            triplets.push_back(Triplet(cellsID[1], cellsID[1], -a_n));
        }
    }

    void _assemble_convection(std::vector<Triplet>& triplets) {
        const std::vector<std::vector<size_t>>& faces_cellsID = this->_mesh_p->getFacesCellsID();
        const std::vector<bool>& faces_boundaryFlag = this->_mesh_p->getFacesBoundaryFlag();
        const std::vector<double>& faces_fx = this->_mesh_p->getFacesFX();
        const std::vector<double>& faces_ux = this->_mesh_p->getFaceData("ux");
        const std::vector<double>& faces_uy = this->_mesh_p->getFaceData("uy");
        const std::vector<double>& faces_uz = this->_mesh_p->getFaceData("uz");
        const std::vector<std::array<double, 3>>& faces_normal_vectors = this->_mesh_p->faces.getFacesNormalVector();

        for(const auto& [cellsID, ux, uy, uz, fx, normal_vector, bd_flag]: std::views::zip(faces_cellsID, faces_ux, faces_uy, faces_uz, faces_fx, faces_normal_vectors, faces_boundaryFlag)) {
            const double F = ux*normal_vector[0] + uy * normal_vector[1] + uz * normal_vector[2];
            const double aPP = fx*F;
            if(bd_flag) {
                triplets.push_back(Triplet(cellsID[0], cellsID[0], aPP));
                continue;
            }

            const double aNN = - fx*F; // Flip flux to be oriented correctly
            const double aPN = (1. - fx)*F;
            const double aNP = (fx - 1.)*F;
    
            triplets.push_back(Triplet(cellsID[0], cellsID[0], aPP));
            triplets.push_back(Triplet(cellsID[1], cellsID[1], aNN));
            triplets.push_back(Triplet(cellsID[0], cellsID[1], aPN));
            triplets.push_back(Triplet(cellsID[1], cellsID[0], aNP));
        }
    }
    void _assemble_jacobian(SpMat& A) {
        // _triplets.reserve(some_estimation)
        std::cout << "Assembling laplacian... \n";
        _assemble_laplacian(_triplets);
        std::cout << "Assembling convection... \n";
        _assemble_convection(_triplets);
        std::cout << "Assembling ddt... \n";
        _assemble_ddt_A(_triplets);
        A.setFromTriplets(_triplets.begin(), _triplets.end());
        // Free memory
        _triplets.resize(0);
        std::cout << "--- Sparse Jacobian assembled --- \n";
    }

    void _assemble_RHS(Eigen::VectorXd& b) {
        _assemble_ddt_b(b);
    }
};

}

#endif //IFD_OPERATORS_HPP
