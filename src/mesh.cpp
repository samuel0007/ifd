
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <vector>
#include <string_view>
#include <stdexcept>
#include <string>
#include <cassert>
#include <tuple>
#include <ranges>
#include <initializer_list>


#include "../include/ifd_lib/mesh.hpp"

size_t assert_equal_size(auto t) {
    return t.size();
}
size_t assert_equal_size(auto t, auto... args) {
    assert(t.size() == assert_equal_size(args...));
    return t.size();
};

namespace utils {
    std::array<double, 3> cross_product(const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return {
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0],
        };
    };

    double dot_product(const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    std::array<double, 3> createVec(const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return {
            b[0] - a[0],
            b[1] - a[1],
            b[2] - a[2],
        };
    }


    double magnitude(const std::array<double, 3>& a) {
        return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }
}

namespace mesh {

    namespace parser {
        void clean_line(std::string &line) {
            if (line[line.size() - 1] == '\r') line.erase(line.size() - 1);
        }

        void remove_parentheses(std::string &line) {
            if (line[0] != '(') throw std::invalid_argument("");
            line.erase(0, 1);
            if (line[line.size() - 1] != ')') throw std::invalid_argument("");
            line.erase(line.size() - 1);
        }

        std::tuple<size_t, size_t>
        parse_points_ascii(std::ifstream &points_file, std::vector<std::array<point_t, 3>> &data) {
            std::string line;

            // First line is number of points
            getline(points_file, line);
            clean_line(line);

            size_t n = boost::lexical_cast<int>(line);
            data.resize(n);

            // Second line should be an open parenthesis
            getline(points_file, line);
            clean_line(line);
            if (line != "(") throw std::invalid_argument("");

            // Find dimension from first point
            getline(points_file, line);
            clean_line(line);
            remove_parentheses(line);

            std::stringstream first_point_ss(line);
            point_t value;
            std::array<point_t, 3> first_point;
            int i = 0;
            while (first_point_ss >> value) {
                first_point[i] = value;
                ++i;
            }
            size_t dim = first_point.size();
            data[0] = first_point;

            int point_idx = 1;
            while (getline(points_file, line)) {
                // First character should be an open parenthesis, last should be a closing one.
                clean_line(line);
                if (line == ")") break;
                remove_parentheses(line);
                std::array<point_t, 3> point;
                std::stringstream point_ss(line);
                int i = 0;
                while (point_ss >> value) {
                    point[i] = value;
                    ++i;
                }
                data[point_idx] = point;
                if (i != dim) throw std::invalid_argument("Points don't have uniform dimension");
                ++point_idx;
            }
            if (point_idx != n) throw std::invalid_argument("Invalid point number.");

            return {n, dim};
        }

        template<typename T>
        size_t parse_object_ascii(std::ifstream &objects_file, std::vector<std::vector<T>> &data) {
            std::string line;

            // First line is number of objects
            getline(objects_file, line);
            clean_line(line);

            size_t n = boost::lexical_cast<int>(line);
            data.resize(n);

            // Second line should be an open parenthesis
            getline(objects_file, line);
            clean_line(line);
            if (line != "(") throw std::invalid_argument("");

            size_t object_idx = 0;
            while (getline(objects_file, line)) {
                clean_line(line);
                if (line == ")") break;
                int obj_dim = stoi(line.substr(0, line.find_first_of('(')));
                line.erase(0, line.find_first_of('('));
                remove_parentheses(line);
                std::vector<T> object(obj_dim);
                std::stringstream object_ss(line);
                T value;
                int i = 0;
                while (object_ss >> value) {
                    object[i] = value;
                    ++i;
                }
                data[object_idx] = object;
                if (i != obj_dim) throw std::invalid_argument("Points don't have uniform dimension");
                ++object_idx;
            }
            if (object_idx != n) throw std::invalid_argument("Invalid point number.");

            return n;
        }

        size_t parse_boundary_ascii(std::ifstream &boundary_file, std::vector<Wall> &data) {
            std::string line;

            // First line is number of walls
            getline(boundary_file, line);
            clean_line(line);
            size_t n = boost::lexical_cast<int>(line);
            data.resize(n);

            // Second line should be an open parenthesis
            getline(boundary_file, line);
            clean_line(line);
            if (line != "(") throw std::invalid_argument("");

            for (int i = 0; i < n; ++i) {
                getline(boundary_file, line);
                clean_line(line);
                // First line is label
                Wall wall;
                wall.label = line;

                // Second line is number of faces
                getline(boundary_file, line);
                clean_line(line);
                int wall_dim = stoi(line);
                wall.facesID.resize(wall_dim);

                // Third line is open parenthesis
                getline(boundary_file, line);
                clean_line(line);
                if (line != "(") throw std::invalid_argument("");

                // Fourth line is data
                getline(boundary_file, line);
                clean_line(line);
                std::stringstream faces_idx_ss(line);

                size_t value;
                int j = 0;
                while (faces_idx_ss >> value) {
                    wall.facesID[j] = value;
                    ++j;
                    if (j > wall_dim) throw std::invalid_argument("bad boundary formatting");
                }

                // Fifth line is closing parenthesis
                getline(boundary_file, line);
                clean_line(line);
                if (line != ")") throw std::invalid_argument("");

                data[i] = wall;
            }
            return n;
        }
    }

    Boundary::Boundary(std::ifstream &boundary_file, const Cells& cells, Mesh_file_type FILE_TYPE): _cells(cells) {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            this->_n = parser::parse_boundary_ascii(boundary_file, this->_walls);
        } else {
            assert(false && "File type not implemented");
        }

        this->register_walls_faces_cellID();

        for(const auto& wall_faces_cellID: this->_walls_faces_cellID) {
            for(const auto& wall_face_cellID: wall_faces_cellID) {
                std::cout << wall_face_cellID << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;

    }

    void Boundary::register_walls_faces_cellID() {
        this->_walls_faces_cellID.resize(this->size());
        const auto& faces_cells_id = this->_cells.getFacesCellsID();
        for(const auto& [wall, wall_faces_cellID]: std::views::zip(this->_walls, this->_walls_faces_cellID)) {
            wall_faces_cellID.resize(wall.facesID.size());
            for(const auto& [faceID, face_cellID]: std::views::zip(wall.facesID, wall_faces_cellID)){
                // boundary faces should be attached to only one cell
                assert(faces_cells_id[faceID].size() == 1);
                face_cellID = faces_cells_id[faceID][0];
            }
        }
    }

    void Boundary::print() const {
        for (int i = 0; i < this->_n; ++i) {
            std::cout << this->_walls[i].label << " ";
            size_t n_faces = this->_walls[i].facesID.size();
            for (int j = 0; j < n_faces; ++j) {
                std::cout << this->_walls[i].facesID[j] << " ";
            }
            std::cout << std::endl;
        }
    }

    Points::Points(std::ifstream &points_file, Mesh_file_type FILE_TYPE) {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            std::tie(_n, _dim) = parser::parse_points_ascii(points_file, this->_pointsPOS);
        } else
            assert(false && "File type not implemented");
    }

    void Points::print() const {
        for (int i = 0; i < this->_n; ++i) {
            for (int j = 0; j < this->_dim; ++j) {
                std::cout << this->_pointsPOS[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    const std::vector<std::array<point_t, 3>> &Points::getPointsPOS() const {
        return _pointsPOS;
    }


    void Faces::compute_facesPointPOS() {
        this->_facesPointsPOS.resize(this->size());
        this->_faces_nPoints.resize(this->size());
        assert((this->_facesPointsID.size() == this->_facesPointsPOS.size()) && (this->_facesPointsID.size() == this->_faces_nPoints.size()));
        for (const auto& [facePointsID, facePointsPOS, face_nPoint]: std::views::zip(
                this->_facesPointsID, this->_facesPointsPOS, this->_faces_nPoints))
        {
            face_nPoint = facePointsID.size();
            facePointsPOS.resize(face_nPoint);
            const auto& pointsPOS = this->_points.getPointsPOS();
            std::ranges::transform(facePointsID, facePointsPOS.begin(), [&pointsPOS](size_t pointID){
                return pointsPOS[pointID];
            });
        }
    }

    void Faces::compute_faces_centroid(){
        this->_faces_centroid.resize(this->size());
        std::ranges::fill(this->_faces_centroid, std::array<value_t, 3>{});
        for(const auto& [centroid, pointsPOS, nPoints]: std::views::zip(this->_faces_centroid, this->_facesPointsPOS, this->_faces_nPoints)) {
            for(const auto& pointPOS: pointsPOS){
                std::ranges::transform(centroid, pointPOS, centroid.begin(), std::plus<value_t>{});
            }
            std::ranges::for_each(centroid, [&nPoints](value_t& x){x /= (double)nPoints;});
        }
    }

    void Faces::compute_faces_areas()
    {
        this->_faces_area.resize(this->size());
        this->_faces_area_vector.resize(this->size());
        assert_equal_size(this->_faces_area, this->_faces_centroid, this->_facesPointsPOS, this->_faces_area_vector);
        for(const auto& [area, area_vector, centroid, pointsPOS]: std::views::zip(this->_faces_area, this->_faces_area_vector, this->_faces_centroid, this->_facesPointsPOS)) {
            size_t face_nPoints = pointsPOS.size();
            area = 0.;
            std::ranges::fill(area_vector, 0.);         
            for(int i = 0; i < face_nPoints; ++i){
                const std::array<value_t, 3>& p0 = pointsPOS[i];
                const std::array<value_t, 3>& p1 = pointsPOS[(i+1)%face_nPoints];
                const std::array<value_t, 3>& p2 = centroid;

                std::array<value_t, 3> p2p0 = utils::createVec(p2, p0);
                std::array<value_t, 3> p2p1 = utils::createVec(p2, p1);

                // compute half cross product P2P0 x P2P1 to get normal area vector of each triangle
                std::array<value_t, 3> triangle_normal = utils::cross_product(p2p0, p2p1);
                std::ranges::transform(triangle_normal, triangle_normal.begin(), [](double x){return 0.5*x;});

                // Total area is the sum of all triangles on the face
                area += utils::magnitude(triangle_normal);

                // area_vector = sum_i triangle_normal_i / area
                std::ranges::transform(area_vector, triangle_normal, area_vector.begin(), std::plus<value_t>{});
            }
            std::ranges::transform(area_vector, area_vector.begin(), [area](value_t x){return x / area;});
        }
    }

    Faces::Faces(std::ifstream &faces_file, const Points &points, Mesh_file_type FILE_TYPE) : _points(points)
    {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            this->_n = parser::parse_object_ascii(faces_file, this->_facesPointsID);
        } else
            assert(false && "File type not implemented");

        this->compute_facesPointPOS();
        this->compute_faces_centroid();
        this->compute_faces_areas();

        // std::cout << "NORMAL" << "\n";
        // for(const auto& centroid: this->_faces_area_vector) {
        //     for(const auto& point: centroid) {
        //         std::cout << point << " ";
        //     }
        //     std::cout << "\n";
        // }

        // std::cout << "AREAS" << "\n";

        // for(const auto& centroid: this->_faces_area) {
        //     std::cout << centroid << "\n";
        // }
        // std::cout << std::endl;


    }

    void Faces::printPOS() const {
        for(const auto& facePointsPOS: this->_facesPointsPOS) {
            for(const auto& pointsPOS: facePointsPOS) {
                std::cout << "( ";
                for(const auto& pointPOS: pointsPOS) {
                    std::cout << pointPOS << " ";
                }
                std::cout << ")";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    void Faces::printID() const {
        for (const auto& facePointsID: this->_facesPointsID) {
            for (const auto& pointID: facePointsID) {
                std::cout << pointID << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    void Faces::print() const
    {
        this->printID();
        this->printPOS();
    }

    // template <int DIM>
    Cells::Cells(std::ifstream &cells_file, const Faces &faces, Mesh_file_type FILE_TYPE)
        : _faces(faces)
    {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            this->_n = parser::parse_object_ascii(cells_file, this->_cells_facesID);
        } else
            assert(false && "File type not implemented");

        this->register_cells_faces_pointsPOS();
        this->compute_cells_nPoints();
        this->compute_cells_centroid();
        this->compute_cells_volume();
        this->register_faces_cellsID();
        this->register_cells_ncellsID();
    }

//    template<int DIM>
    void Cells::print() const noexcept {
        this->printID();
    }

//    template<int DIM>
    void Cells::printID() const noexcept {
        for (const auto& cell_facesID: this->_cells_facesID) {
            for (const auto& faceID: cell_facesID) {
                std::cout << faceID << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

//    template<int DIM>
    const std::vector<value_t>& Cells::getCellsVolume() const noexcept {
        return this->_cells_volume;
    }

    // If we were to always retrieve the point, we would not
    // load the correct series of points in the cache as we would always be
    // RA the _points DS
//    template<int DIM>
    void Cells::register_cells_faces_pointsPOS() {
        this->_cells_facesID_pointsPOS.resize(this->_cells_facesID.size());
        for (const auto& [cell_facesID, cell_facesID_pointsPOS]: std::views::zip(this->_cells_facesID,
                                                                              this->_cells_facesID_pointsPOS)) {
            cell_facesID_pointsPOS.resize(cell_facesID.size());
            const auto& facesID_pointsPOS = this->_faces.getFacesPointsPOS();
            std::ranges::transform(cell_facesID, cell_facesID_pointsPOS.begin(), [&](size_t faceID){return facesID_pointsPOS[faceID];});
        }
    }

    void Cells::compute_cells_nPoints() {
        this->_cells_nPoints.resize(this->size());
        std::ranges::fill(this->_cells_nPoints, 0.);
        for(const auto& [cell_facesID_pointsPOS, cell_nPoints]: std::views::zip(this->_cells_facesID_pointsPOS, this->_cells_nPoints)) {
            for(const auto& faceID_pointsPOS: cell_facesID_pointsPOS){
                cell_nPoints += faceID_pointsPOS.size();
            }
        }
    }
//    template<int DIM>
    void Cells::compute_cells_centroid() {
        this->_cells_centroid.resize(this->size());
        std::ranges::fill(this->_cells_centroid, std::array<value_t, 3>{});
        assert((this->_cells_facesID_pointsPOS.size() == this->_cells_centroid.size()) && (this->_cells_centroid.size() == this->_cells_nPoints.size()));
        for(const auto& [cell_facesID_pointsPOS, cell_centroid, cell_nPoints]: std::views::zip(this->_cells_facesID_pointsPOS, this->_cells_centroid, this->_cells_nPoints)){
            for(const auto& faceID_pointsPOS: cell_facesID_pointsPOS) {
                for(const auto& pointsPOS: faceID_pointsPOS) {
                    std::ranges::transform(cell_centroid, pointsPOS, cell_centroid.begin(), std::plus<value_t>{});
                }
            }
            std::ranges::for_each(cell_centroid, [&cell_nPoints](value_t& x){x /= (double)cell_nPoints;});
        }
    }

//    template<int DIM>
    void Cells::compute_cells_volume() {
        // for each cell find all their faces
        const auto& facesID_centroid = this->_faces.getFacesCentroids();
        this->_cells_volume.resize(this->size());
        for(const auto&[cell_facesID_pointsPOS, cell_centroid, cell_facesID, cell_volume]: std::views::zip(this->_cells_facesID_pointsPOS, this->_cells_centroid, this->_cells_facesID, this->_cells_volume)){
            for(const auto&[pointsPOS, faceID] : std::views::zip(cell_facesID_pointsPOS, cell_facesID)) {
                size_t face_nPoints = pointsPOS.size();
                std::array<value_t, 3> face_centroid = facesID_centroid[faceID];
                for(size_t i = 0; i < face_nPoints; ++i) {
                    const std::array<value_t, 3>& p0 = face_centroid;
                    const std::array<value_t, 3>& p1 = pointsPOS[i];
                    const std::array<value_t, 3>& p2 = pointsPOS[(i+1)%face_nPoints];
                    const std::array<value_t, 3>& p3 = cell_centroid;
                    
                    const std::array<value_t, 3> b1 = utils::createVec(p3, p0);
                    const std::array<value_t, 3> b2 = utils::createVec(p3, p1);
                    const std::array<value_t, 3> b3 = utils::createVec(p3, p2);

                    cell_volume += fabs(utils::dot_product(b1, utils::cross_product(b2, b3))) / 6.;
                }
            }
        }
    }

    void Cells::register_faces_cellsID() {
        // For each cell (in order), register one self on each of the faces.
        this->_faces_cellsID.resize(this->_faces.size());
        const size_t N = this->size();
        for(size_t cell_idx = 0; cell_idx < N; ++cell_idx) {
            for(const auto& face_ID: this->_cells_facesID[cell_idx]){
                this->_faces_cellsID[face_ID].push_back(cell_idx);
            }
        }
    }

    void Cells::register_cells_ncellsID() {
        this->_cells_ncellsID.resize(this->size());
        for(const auto& face_cellsID: this->_faces_cellsID) {
            std::vector<size_t> neighboursID = face_cellsID;
            for(const auto& cellID: face_cellsID) {
                for(const auto& neighbourID: neighboursID){
                    if(neighbourID != cellID)
                        _cells_ncellsID[cellID].push_back(neighbourID);
                }
            }
        }
    }
//    template<int DIM>
    std::shared_ptr<Mesh> load_mesh(fs::path &mesh_folder_path) {
        // Try to open points, faces, cells, and boundary file
        fs::path points_file_path("points");
        fs::path faces_file_path("faces");
        fs::path cells_file_path("cells");
        fs::path boundary_file_path("boundary");

        std::ifstream points_file(mesh_folder_path / points_file_path);
        std::ifstream faces_file(mesh_folder_path / faces_file_path);
        std::ifstream cells_file(mesh_folder_path / cells_file_path);
        std::ifstream boundary_file(mesh_folder_path / boundary_file_path);

        if (points_file.fail()) {
            std::cout << mesh_folder_path / points_file_path << std::endl;
            throw std::invalid_argument("Indicated point file does not exist.");
        }
        if (faces_file.fail()) {
            std::cout << mesh_folder_path / faces_file_path << std::endl;
            throw std::invalid_argument("Indicated faces file does not exist.");
        }
        if (cells_file.fail()) {
            std::cout << mesh_folder_path / cells_file_path << std::endl;
            throw std::invalid_argument("Indicated cells file does not exist.");
        }
        if (boundary_file.fail()) {
            std::cout << mesh_folder_path / boundary_file_path << std::endl;
            throw std::invalid_argument("Indicated boundary file does not exist.");
        }


        // Parse files and create geometry objects
        Points points(points_file);
        Faces faces(faces_file, points);
        Cells cells(cells_file, faces);
        Boundary boundary(boundary_file, cells);

        return std::make_shared<Mesh>(points, faces, cells, boundary);
    }
}