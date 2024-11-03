#include "include/ifd_lib/mesh.hpp"
#include "include/ifd_lib/operators.hpp"
// #include "include/ifd_lib/utils.hpp"
#include "include/ifd_lib/vtkwriter.hpp"

#include <algorithm>
#include <vector>
#include <string_view>
#include <stdexcept>
#include <filesystem>
#include <Eigen/Dense>

namespace fs = std::filesystem;

struct Global_Options {
    fs::path mesh_folder_path;
    fs::path output_folder_path;

};




Global_Options parse_args(int argc, char *argv[]) {
    Global_Options options;
    const std::vector<std::string_view> args(argv + 1, argv + argc);
    if (argc < 3) {
        throw std::invalid_argument("User has to pass mesh and output folder");
    } 
    options.mesh_folder_path = fs::path(args[0]);
    options.output_folder_path = fs::path(args[1]);
    return options;
};

int main(int argc, char *argv[]) {
    // Load data
    Global_Options options = parse_args(argc, argv);
    const size_t N = 30;
    const double L = 1.;

    const double dx = L/N;
    // const double dx = 1.;
    const double gamma = 1e-5;
    //double T = 0.2;
    const int steps = 1000;
    double dt = 0.1;

    auto IC = [](std::array<double, 3> x) -> double {
        return std::max<double>(0., 1. - 4.*sqrt((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.25)*(x[1]-0.25)));
    };

    auto Ux = [](std::array<double, 3> x) -> double {
        return 0.;
        // return 0.0002 * std::sin(M_PI*x[0])*std::cos(M_PI*x[1]);
    };

    auto Uy = [](std::array<double, 3> x) -> double {
        return 0.;
        // return - 0.0002 * std::sin(M_PI*x[1])*std::cos(M_PI*x[0]);
    };

    auto Uz = [](std::array<double, 3> x) -> double {
        return 0.;
    };

    std::array<size_t, 3> Ns({N, N, 1});
    std::array<double, 3> DXs({dx, dx, dx});

    std::shared_ptr<mesh::Mesh> mesh_p = mesh::create_cartesian_mesh(Ns, DXs);
    VTKWriter vtkwriter(options.output_folder_path, mesh_p);

    operators::SparseTransientLaplacianProblem<Eigen::SimplicialCholesky<operators::SpMat, Eigen::Upper>> problem(gamma, IC, dt, mesh_p);
    problem.initialize();

    for(int i = 0; i < steps; ++i) {
        std::cout << "#step "<< i << "\n";
        problem.evolve();
        if(i%10 != 0) continue;
        // Pack data
        std::vector<double> output(problem.x.size());
        for(int j = 0; j < problem.x.size(); ++j) {
            output[j] = problem.x(j);
        }

        std::ostringstream filename;
        filename << "test" << "_" << std::setw(4) << std::setfill('0') << i;

        vtkwriter.write(output, filename.str(), "temperature");
    }

    return 0;
}