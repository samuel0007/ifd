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
#include <set>
#include <map>


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
        Points(std::ifstream &points_file, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);
        Points(std::vector<std::array<point_t, 3>> pointsPOS);


        [[nodiscard]] size_t size() const {
            return this->_n;
        }
        [[nodiscard]] const std::vector<std::array<point_t, 3>> &getPointsPOS() const {
            return _pointsPOS;
        }

        void print() const;

    private:
        size_t _n{};
        std::vector<std::array<point_t, 3>> _pointsPOS;
    };

    class Faces {
    public:
        Faces(std::ifstream &faces_file, const Points& points, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);
        Faces(const Points &points, std::vector<std::vector<size_t>> facesPointsID);

        void print() const;

        [[nodiscard]] size_t size() const {
            return this->_n;
        };

        [[nodiscard]] const std::vector<std::vector<std::array<value_t, 3>>>& getFacesPointsPOS() const {
            return this->_facesPointsPOS;
        }

        const std::vector<std::vector<size_t>>& getFacesPointsID() const {
            return this->_facesPointsID;
        }

        const std::vector<std::array<value_t, 3>>& getFacesCentroid() const {
            return this->_faces_centroid;
        }

        const std::vector<std::array<value_t, 3>> getFacesNormalVector() const {
            return this->_faces_area_vector;
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

        void _init_faces();
        void compute_facesPointPOS();
        void compute_faces_centroid();
        void compute_faces_areas();
        void printID() const;
        void printPOS() const;
    };

    class Cells {
    public:
        Cells(std::ifstream &cells_file, const Faces& faces, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);
        Cells(const std::vector<std::vector<size_t>>& cells_facesID, const Faces &faces);

        [[nodiscard]] size_t size() const noexcept {
            return this->_n;
        };

        [[nodiscard]] const std::vector<value_t>& getCellsVolume() const noexcept;
        [[nodiscard]] const std::vector<std::vector<size_t>>& getCellsFacesID() const noexcept {
            return this->_cells_facesID;
        };

        [[nodiscard]] const std::vector<std::set<size_t>>& getCellsPointsID() const noexcept {
            return this->_cells_pointsID;
        };

        [[nodiscard]] const std::vector<std::array<value_t, 3>>& getCellsCentroid() const noexcept {
            return this->_cells_centroid;
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


    private:
        size_t _n{};

        const Faces& _faces;
        std::vector<std::vector<size_t>> _cells_facesID;

        std::vector<std::set<size_t>> _cells_pointsID;

        std::vector<std::vector<std::vector<std::array<value_t, 3>>>> _cells_facesID_pointsPOS;

        std::vector<std::array<value_t, 3>> _cells_centroid;
        std::vector<value_t> _cells_volume;
        std::vector<size_t> _cells_nPoints;
        
        void _init_cells();
        void compute_cells_nPoints();
        void compute_cells_centroid();
        void compute_cells_volume();
        void register_cells_faces_pointsPOS();
        void register_cells_pointsID();
        void printID() const noexcept;
    };


    class Boundary {
    public:
        Boundary(std::ifstream &boundary_file, const Cells& cells, Mesh_file_type FILE_TYPE = Mesh_file_type::ASCII);
        Boundary(const std::vector<Wall>& walls, const Cells& cells);
        
        void print() const;

        [[nodiscard]] size_t size() const {
            return this->_n;
        };

        std::vector<Wall> getWalls() const noexcept {
            return this->_walls;
        }

    private:
        size_t _n{};
        std::vector<Wall> _walls;
        const Cells& _cells;
    };

    class Mesh {
    public:
        Mesh(Points _points, Faces _faces, Cells _cells, Boundary _boundary);
        Cells cells;
        Boundary boundary;
        Faces faces;
        Points points;

        void summarize(std::ofstream& out) {};

        void addCellData (const std::vector<value_t>& cell_data, std::string label);

        value_t interpolateFaceValue(value_t phiP, value_t phiN, double fx) {
            return phiP * fx + phiN * (1. - fx);
        }

        const std::vector<value_t>& getCellData(std::string label) {
            return this->_container_cells_data[_container_map[label]];
        }

        const std::vector<value_t>& getFaceData(std::string label) {
            return this->_container_faces_data[_container_map[label]];
        }

        const std::vector<std::vector<size_t>>& getFacesCellsID() const {
            return this->_faces_cellsID;
        }

        const std::vector<bool>& getFacesBoundaryFlag() const {
            return this->_faces_boundaryFlag;
        }

        const std::vector<bool>& getCellsBoundaryFlag() const {
            return this->_cells_boundaryFlag;
        }
        
        const std::vector<value_t>& getFacesDelta() const {
            return this->_faces_delta;
        }

        const std::vector<value_t>& getFacesFX() const {
            return this->_faces_fx;
        }

    private:
        void registerBoundaryFaces();
        void registerFacesCellsID();
        void register_cells_ncellsID();
        void registerWallsFacesCellID();

        void computeCellToFaceRatio();
        void computeFaceInterpolatedValue();

        // Faces data
        std::vector<value_t> _faces_fx;
        std::vector<value_t> _faces_delta;
        std::vector<std::vector<size_t>> _faces_cellsID;
        std::vector<std::vector<std::array<value_t, 3>>> _faces_cellsCentroid;
        std::vector<bool> _faces_boundaryFlag;

        // Cells data
        std::vector<std::vector<size_t>> _cells_ncellsID;
        std::vector<bool> _cells_boundaryFlag;

        // Boundary data
        std::vector<std::vector<size_t>> _walls_faces_cellID;

        // Mesh Data container
        std::map<std::string, size_t> _container_map;
        std::vector<std::vector<value_t>> _container_cells_data;
        std::vector<std::vector<value_t>> _container_faces_data;
        size_t field_idx = 0;
    };

    std::shared_ptr<Mesh> load_mesh(fs::path &mesh_folder_path);

    std::shared_ptr<Mesh> create_cartesian_mesh(std::array<size_t, 3> Ns, std::array<double, 3> DXs);

}

#endif //IFD_MESH_HPP
