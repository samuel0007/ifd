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
#include <ranges>

namespace fs = std::filesystem;


namespace mesh {
    enum class Mesh_file_type {
        ASCII
    };

    struct Wall {
        std::vector<size_t> facesID;
        std::string label;
    };

    typedef double point_t;
    typedef double value_t;


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
        std::vector<std::array<point_t, 3>> _pointsPOS;
    public:
        [[nodiscard]] const std::vector<std::array<point_t, 3>> &getPointsPOS() const;
    };

    class Faces {
    public:
        explicit Faces(std::ifstream &faces_file, const Points& points, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);

        void print() const;

        [[nodiscard]] size_t size() const {
            return this->_n;
        };

        [[nodiscard]] const std::vector<std::vector<std::array<value_t, 3>>>& getFacesPointsPOS() const {
            return this->_facesPointsPOS;
        }

        const std::vector<std::array<value_t, 3>>& getFacesCentroids() const {
            return this->_faces_centroid;
        }

    private:
        size_t _n{};
        std::vector<std::vector<size_t>> _facesPointsID;
        std::vector<std::vector<std::array<value_t, 3>>> _facesPointsPOS;
        const Points& _points;
        std::vector<std::array<value_t, 3>> _faces_centroid;
        std::vector<value_t> _faces_nPoints;
        std::vector<value_t> _faces_area;
        std::vector<std::array<value_t, 3>> _faces_area_vector;

        void compute_facesPointPOS();
        void compute_faces_centroid();
        void compute_faces_areas();
        void printID() const;
        void printPOS() const;
    };

//    template <int DIM>
    class Cells {
    public:
        explicit Cells(std::ifstream &cells_file, const Faces& faces, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);

        [[nodiscard]] size_t size() const noexcept {
            return this->_n;
        };

        [[nodiscard]] const std::vector<value_t>& getCellsVolume() const noexcept;
        [[nodiscard]] const std::vector<std::vector<size_t>>& getFacesCellsID() const noexcept {
            return this->_faces_cellsID;
        };


        void print() const noexcept;

        auto x_centroid() const {
            return this->_cells_centroid | std::views::transform([](const auto& centroid){return centroid[0];});
        }
        auto y_centroid() const {
            return this->_cells_centroid | std::views::transform([](const auto& centroid){return centroid[1];});
        }
        auto z_centroid() const {
            return this->_cells_centroid | std::views::transform([](const auto& centroid){return centroid[2];});
        }
        std::vector<std::array<value_t,3>> _cells_centroid;


    private:
        size_t _n{};

        const Faces& _faces;
        std::vector<std::vector<size_t>> _cells_facesID;
        std::vector<std::vector<std::vector<std::array<value_t, 3>>>> _cells_facesID_pointsPOS;

        std::vector<std::vector<size_t>> _faces_cellsID; // In the future, we might have faces that connect to 4 cells
        std::vector<std::vector<size_t>> _cells_ncellsID;
        std::vector<value_t> _cells_volume;
        std::vector<size_t> _cells_nPoints;
        
        void compute_cells_nPoints();
        void compute_cells_centroid();
        void compute_cells_volume();
        void register_cells_faces_pointsPOS();
        void register_faces_cellsID();
        void register_cells_ncellsID();
        void printID() const noexcept;
    };


    class Boundary {
    public:
        explicit Boundary(std::ifstream &boundary_file, const Cells& cells, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);

        void print() const;

        [[nodiscard]] size_t size() const {
            return this->_n;
        };

    private:
        void register_walls_faces_cellID();
        size_t _n{};
        std::vector<Wall> _walls;
        std::vector<std::vector<size_t>> _walls_faces_cellID;
        const Cells& _cells;
    };

//    template <int DIM>
    class Mesh {
    public:
        Mesh(Points _points, Faces _faces, Cells _cells, Boundary _boundary) : cells(std::move(_cells)),
                                                                               boundary(std::move(_boundary)),
                                                                               faces(std::move(_faces)),
                                                                               points(std::move(_points)) {};
//        Cells<DIM> cells;
        Cells cells;
        Boundary boundary;
        Faces faces;
        Points points;

        void summarize(std::ofstream& out) {

        };

        void plot() {
        };
    };

//    template <int DIM>
    std::shared_ptr<Mesh> load_mesh(fs::path &mesh_folder_path);;
}

#endif //IFD_MESH_HPP
