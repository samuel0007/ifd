
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <vector>
#include <string_view>
#include <stdexcept>
#include <string>
#include <cassert>
#include <tuple>

#include "../include/ifd_lib/mesh.hpp"

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
        parse_points_ascii(std::ifstream &points_file, std::vector<std::vector<point_t>> &data) {
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
            double value;
            std::vector<double> first_point;
            while (first_point_ss >> value) {
                first_point.push_back(value);
            }
            size_t dim = first_point.size();
            data[0] = first_point;

            int point_idx = 1;
            while (getline(points_file, line)) {
                // First character should be an open parenthesis, last should be a closing one.
                clean_line(line);
                if (line == ")") break;
                remove_parentheses(line);
                std::vector<point_t> point(dim);
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
                wall.faces_idx.resize(wall_dim);

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
                    wall.faces_idx[j] = value;
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

    Boundary::Boundary(std::ifstream &boundary_file, Mesh_file_type FILE_TYPE) {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            this->_n = parser::parse_boundary_ascii(boundary_file, this->_data);
        } else
            assert(false && "File type not implemented");
    }

    void Boundary::print() const {
        for (int i = 0; i < this->_n; ++i) {
            std::cout << this->_data[i].label << " ";
            size_t n_faces = this->_data[i].faces_idx.size();
            for (int j = 0; j < n_faces; ++j) {
                std::cout << this->_data[i].faces_idx[j] << " ";
            }
            std::cout << std::endl;
        }
    }

    Points::Points(std::ifstream &points_file, Mesh_file_type FILE_TYPE) {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            std::tie(_n, _dim) = parser::parse_points_ascii(points_file, this->_data);
        } else
            assert(false && "File type not implemented");
    }

    void Points::print() const {
        for (int i = 0; i < this->_n; ++i) {
            for (int j = 0; j < this->_dim; ++j) {
                std::cout << this->_data[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    Faces::Faces(std::ifstream &faces_file, Mesh_file_type FILE_TYPE) {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            this->_n = parser::parse_object_ascii(faces_file, this->_data);
        } else
            assert(false && "File type not implemented");
    }

    void Faces::print() const {
        for (int i = 0; i < this->_n; ++i) {
            size_t face_size = this->_data[i].size();
            for (int j = 0; j < face_size; ++j) {
                std::cout << this->_data[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    Cells::Cells(std::ifstream &cells_file, Mesh_file_type FILE_TYPE) {
        if (FILE_TYPE == Mesh_file_type::ASCII) {
            this->_n = parser::parse_object_ascii(cells_file, this->_data);
        } else
            assert(false && "File type not implemented");
    }

    void Cells::print() const {
        for (int i = 0; i < this->_n; ++i) {
            size_t cell_dim = this->_data[i].size();
            for (int j = 0; j < cell_dim; ++j) {
                std::cout << this->_data[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

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
        Faces faces(faces_file);
        Cells cells(cells_file);
        Boundary boundary(boundary_file);

        return std::make_shared<Mesh>(points, faces, cells, boundary);
    }
}