//
// Created by russo on 11.10.2024.
//

#ifndef IFD_MESH_HPP
#define IFD_MESH_HPP

#include <boost/lexical_cast.hpp>
#include <iostream>
#include <utility>
#include <vector>
#include <string_view>
#include <stdexcept>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <string>
#include <cassert>
#include <tuple>

namespace fs = std::filesystem;

namespace mesh {
    enum class Mesh_file_type {
        ASCII
    };

    struct Wall {
        std::vector<size_t> faces_idx;
        std::string label;
    };

    typedef double point_t;


    class Points {
    public:
        explicit Points(std::ifstream &points_file, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);;

        [[nodiscard]] size_t size() const {
            return this->_n;
        }

        void print() const;

    private:
        size_t _n{};
        size_t _dim{};
        std::vector<std::vector<point_t>> _data;
    };


    class Faces {
    public:
        explicit Faces(std::ifstream &faces_file, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);;

        void print() const;

        [[nodiscard]] size_t size() const {
            return this->_n;
        };

    private:
        size_t _n{};
        std::vector<std::vector<size_t>> _data;
    };

    class Cells {
    public:
        explicit Cells(std::ifstream &cells_file, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);;

        void print() const;

        [[nodiscard]] size_t size() const {
            return this->_n;
        };


    private:
        size_t _n{};
        std::vector<std::vector<size_t>> _data;
    };


    class Boundary {
    public:
        explicit Boundary(std::ifstream &boundary_file, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);

        void print() const;

        [[nodiscard]] size_t size() const {
            return this->_n;
        };

    private:
        size_t _n{};
        std::vector<Wall> _data;
    };


    class Mesh {
    public:
        Mesh(Points _points, Faces _faces, Cells _cells, Boundary _boundary) : cells(std::move(_cells)),
                                                                               boundary(std::move(_boundary)),
                                                                               faces(std::move(_faces)),
                                                                               points(std::move(_points)) {};
        Cells cells;
        Boundary boundary;
        Faces faces;
        Points points;
    };


    std::shared_ptr<Mesh> load_mesh(fs::path &mesh_folder_path);;
}

#endif //IFD_MESH_HPP
