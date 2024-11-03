
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
#include "../include/ifd_lib/utils.hpp"


size_t assert_equal_size(auto t)
{
    return t.size();
}
size_t assert_equal_size(auto t, auto... args)
{
    assert(t.size() == assert_equal_size(args...));
    return t.size();
};


namespace mesh
{

    namespace parser
    {
        void clean_line(std::string &line)
        {
            if (line[line.size() - 1] == '\r')
                line.erase(line.size() - 1);
        }

        void remove_parentheses(std::string &line)
        {
            if (line[0] != '(')
                throw std::invalid_argument("");
            line.erase(0, 1);
            if (line[line.size() - 1] != ')')
                throw std::invalid_argument("");
            line.erase(line.size() - 1);
        }

        size_t parse_points_ascii(std::ifstream &points_file, std::vector<std::array<point_t, 3>> &data)
        {
            std::string line;

            // First line is number of points
            getline(points_file, line);
            clean_line(line);

            size_t n = boost::lexical_cast<int>(line);
            data.resize(n);

            // Second line should be an open parenthesis
            getline(points_file, line);
            clean_line(line);
            if (line != "(")
                throw std::invalid_argument("");

            // Find dimension from first point
            getline(points_file, line);
            clean_line(line);
            remove_parentheses(line);

            std::stringstream first_point_ss(line);
            point_t value;
            std::array<point_t, 3> first_point;
            int i = 0;
            while (first_point_ss >> value)
            {
                first_point[i] = value;
                ++i;
            }
            size_t dim = first_point.size();
            assert(dim == 3);
            data[0] = first_point;

            int point_idx = 1;
            while (getline(points_file, line))
            {
                // First character should be an open parenthesis, last should be a closing one.
                clean_line(line);
                if (line == ")")
                    break;
                remove_parentheses(line);
                std::array<point_t, 3> point;
                std::stringstream point_ss(line);
                int i = 0;
                while (point_ss >> value)
                {
                    point[i] = value;
                    ++i;
                }
                data[point_idx] = point;
                if (i != dim)
                    throw std::invalid_argument("Points don't have uniform dimension");
                ++point_idx;
            }
            if (point_idx != n)
                throw std::invalid_argument("Invalid point number.");

            return n;
        }

        template <typename T>
        size_t parse_object_ascii(std::ifstream &objects_file, std::vector<std::vector<T>> &data)
        {
            std::string line;

            // First line is number of objects
            getline(objects_file, line);
            clean_line(line);

            size_t n = boost::lexical_cast<int>(line);
            data.resize(n);

            // Second line should be an open parenthesis
            getline(objects_file, line);
            clean_line(line);
            if (line != "(")
                throw std::invalid_argument("");

            size_t object_idx = 0;
            while (getline(objects_file, line))
            {
                clean_line(line);
                if (line == ")")
                    break;
                int obj_dim = stoi(line.substr(0, line.find_first_of('(')));
                line.erase(0, line.find_first_of('('));
                remove_parentheses(line);
                std::vector<T> object(obj_dim);
                std::stringstream object_ss(line);
                T value;
                int i = 0;
                while (object_ss >> value)
                {
                    object[i] = value;
                    ++i;
                }
                data[object_idx] = object;
                if (i != obj_dim)
                    throw std::invalid_argument("Points don't have uniform dimension");
                ++object_idx;
            }
            if (object_idx != n)
                throw std::invalid_argument("Invalid point number.");

            return n;
        }

        size_t parse_boundary_ascii(std::ifstream &boundary_file, std::vector<Wall> &data)
        {
            std::string line;

            // First line is number of walls
            getline(boundary_file, line);
            clean_line(line);
            size_t n = boost::lexical_cast<int>(line);
            data.resize(n);

            // Second line should be an open parenthesis
            getline(boundary_file, line);
            clean_line(line);
            if (line != "(")
                throw std::invalid_argument("");

            for (int i = 0; i < n; ++i)
            {
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
                if (line != "(")
                    throw std::invalid_argument("");

                // Fourth line is data
                getline(boundary_file, line);
                clean_line(line);
                std::stringstream faces_idx_ss(line);

                size_t value;
                int j = 0;
                while (faces_idx_ss >> value)
                {
                    wall.facesID[j] = value;
                    ++j;
                    if (j > wall_dim)
                        throw std::invalid_argument("bad boundary formatting");
                }

                // Fifth line is closing parenthesis
                getline(boundary_file, line);
                clean_line(line);
                if (line != ")")
                    throw std::invalid_argument("");

                data[i] = wall;
            }
            return n;
        }
    }

    Boundary::Boundary(std::ifstream &boundary_file, const Cells &cells, Mesh_file_type FILE_TYPE) : _cells(cells)
    {
        if (FILE_TYPE == Mesh_file_type::ASCII)
        {
            this->_n = parser::parse_boundary_ascii(boundary_file, this->_walls);
        }
        else
        {
            assert(false && "File type not implemented");
        }
    }

    Boundary::Boundary(const std::vector<Wall> &walls, const Cells &cells) : _cells(cells), _walls(walls)
    {
        this->_n = walls.size();
    }

    void Boundary::print() const
    {
        for (int i = 0; i < this->_n; ++i)
        {
            std::cout << this->_walls[i].label << " ";
            size_t n_faces = this->_walls[i].facesID.size();
            for (int j = 0; j < n_faces; ++j)
            {
                std::cout << this->_walls[i].facesID[j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // --- Points --- //
    Points::Points(std::ifstream &points_file, Mesh_file_type FILE_TYPE)
    {
        if (FILE_TYPE == Mesh_file_type::ASCII)
        {
            this->_n = parser::parse_points_ascii(points_file, this->_pointsPOS);
        }
        else
            assert(false && "File type not implemented");
    }

    Points::Points(std::vector<std::array<value_t, 3>> pointsPOS)
    {
        this->_pointsPOS = pointsPOS;
        this->_n = pointsPOS.size();
    }

    void Points::print() const
    {
        for (int i = 0; i < this->_n; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                std::cout << this->_pointsPOS[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    // --- Faces --- //
    void Faces::compute_facesPointPOS()
    {
        this->_facesPointsPOS.resize(this->size());
        this->_faces_nPoints.resize(this->size());
        assert((this->_facesPointsID.size() == this->_facesPointsPOS.size()) && (this->_facesPointsID.size() == this->_faces_nPoints.size()));
        for (const auto &[facePointsID, facePointsPOS, face_nPoint] : std::views::zip(
                 this->_facesPointsID, this->_facesPointsPOS, this->_faces_nPoints))
        {
            face_nPoint = facePointsID.size();
            facePointsPOS.resize(face_nPoint);
            const auto &pointsPOS = this->_points.getPointsPOS();
            std::ranges::transform(facePointsID, facePointsPOS.begin(), [&pointsPOS](size_t pointID)
                                   { return pointsPOS[pointID]; });
        }
    }

    void Faces::compute_faces_centroid()
    {
        this->_faces_centroid.resize(this->size());
        std::ranges::fill(this->_faces_centroid, std::array<value_t, 3>{});
        for (const auto &[centroid, pointsPOS, nPoints] : std::views::zip(this->_faces_centroid, this->_facesPointsPOS, this->_faces_nPoints))
        {
            for (const auto &pointPOS : pointsPOS)
            {
                std::ranges::transform(centroid, pointPOS, centroid.begin(), std::plus<value_t>{});
            }
            std::ranges::for_each(centroid, [&nPoints](value_t &x)
                                  { x /= (double)nPoints; });
        }
    }

    void Faces::compute_faces_areas()
    {
        this->_faces_area.resize(this->size());
        this->_faces_area_vector.resize(this->size());
        assert_equal_size(this->_faces_area, this->_faces_centroid, this->_facesPointsPOS, this->_faces_area_vector);
        for (const auto &[area, area_vector, centroid, pointsPOS] : std::views::zip(this->_faces_area, this->_faces_area_vector, this->_faces_centroid, this->_facesPointsPOS))
        {
            size_t face_nPoints = pointsPOS.size();
            area = 0.;
            std::ranges::fill(area_vector, 0.);
            for (int i = 0; i < face_nPoints; ++i)
            {
                const std::array<value_t, 3> &p0 = pointsPOS[i];
                const std::array<value_t, 3> &p1 = pointsPOS[(i + 1) % face_nPoints];
                const std::array<value_t, 3> &p2 = centroid;

                std::array<value_t, 3> p2p0 = utils::createVec(p2, p0);
                std::array<value_t, 3> p2p1 = utils::createVec(p2, p1);

                // compute half cross product P2P0 x P2P1 to get normal area vector of each triangle
                std::array<value_t, 3> triangle_normal = utils::cross_product(p2p0, p2p1);
                std::ranges::transform(triangle_normal, triangle_normal.begin(), [](double x)
                                       { return 0.5 * x; });

                // Total area is the sum of all triangles on the face
                area += utils::magnitude(triangle_normal);

                // area_vector = sum_i triangle_normal_i / area
                std::ranges::transform(area_vector, triangle_normal, area_vector.begin(), std::plus<value_t>{});
            }
            std::ranges::transform(area_vector, area_vector.begin(), [area](value_t x)
                                   { return x / area; });
        }
    }

    Faces::Faces(std::ifstream &faces_file, const Points &points, Mesh_file_type FILE_TYPE) : _points(points)
    {
        if (FILE_TYPE == Mesh_file_type::ASCII)
        {
            this->_n = parser::parse_object_ascii(faces_file, this->_facesPointsID);
        }
        else
            assert(false && "File type not implemented");

        this->_init_faces();
    }

    Faces::Faces(const Points &points, std::vector<std::vector<size_t>> facesPointsID) : _points(points), _facesPointsID(facesPointsID)
    {
        this->_n = facesPointsID.size();
        this->_init_faces();
    }

    void Faces::_init_faces()
    {
        this->compute_facesPointPOS();
        this->compute_faces_centroid();
        this->compute_faces_areas();
    }

    void Faces::printPOS() const
    {
        for (const auto &facePointsPOS : this->_facesPointsPOS)
        {
            for (const auto &pointsPOS : facePointsPOS)
            {
                std::cout << "( ";
                for (const auto &pointPOS : pointsPOS)
                {
                    std::cout << pointPOS << " ";
                }
                std::cout << ")";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    void Faces::printID() const
    {
        for (const auto &facePointsID : this->_facesPointsID)
        {
            for (const auto &pointID : facePointsID)
            {
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

    // --- Cells --- //
    Cells::Cells(std::ifstream &cells_file, const Faces &faces, Mesh_file_type FILE_TYPE)
        : _faces(faces)
    {
        if (FILE_TYPE == Mesh_file_type::ASCII)
        {
            this->_n = parser::parse_object_ascii(cells_file, this->_cells_facesID);
        }
        else
            assert(false && "File type not implemented");

        this->_init_cells();
    }

    Cells::Cells(const std::vector<std::vector<size_t>> &cells_facesID, const Faces &faces)
        : _faces(faces), _cells_facesID(cells_facesID)
    {
        this->_n = cells_facesID.size();

        this->_init_cells();
    }

    void Cells::_init_cells()
    {
        this->register_cells_faces_pointsPOS();
        this->register_cells_pointsID();
        this->compute_cells_nPoints();
        this->compute_cells_centroid();
        this->compute_cells_volume();
    }

    void Cells::print() const noexcept
    {
        this->printID();
    }

    void Cells::printID() const noexcept
    {
        for (const auto &cell_facesID : this->_cells_facesID)
        {
            for (const auto &faceID : cell_facesID)
            {
                std::cout << faceID << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    const std::vector<value_t> &Cells::getCellsVolume() const noexcept
    {
        return this->_cells_volume;
    }

    void Cells::register_cells_faces_pointsPOS()
    {
        this->_cells_facesID_pointsPOS.resize(this->_cells_facesID.size());
        for (const auto &[cell_facesID, cell_facesID_pointsPOS] : std::views::zip(this->_cells_facesID,
                                                                                  this->_cells_facesID_pointsPOS))
        {
            cell_facesID_pointsPOS.resize(cell_facesID.size());
            const auto &facesID_pointsPOS = this->_faces.getFacesPointsPOS();
            std::ranges::transform(cell_facesID, cell_facesID_pointsPOS.begin(), [&](size_t faceID)
                                   { return facesID_pointsPOS[faceID]; });
        }
    }

    void Cells::register_cells_pointsID()
    {
        this->_cells_pointsID.resize(this->size());
        std::vector<std::vector<size_t>> faces_pointsID = this->_faces.getFacesPointsID();
        for (const auto &[cell_pointsID, cell_facesID] : std::views::zip(_cells_pointsID, _cells_facesID))
        {
            for (const auto &faceID : cell_facesID)
            {
                const std::vector<size_t> &face_pointsID = faces_pointsID[faceID];
                for (const auto &pointID : face_pointsID)
                {
                    cell_pointsID.insert(pointID);
                }
            }
        }
    }

    void Cells::compute_cells_nPoints()
    {
        this->_cells_nPoints.resize(this->size());
        std::ranges::fill(this->_cells_nPoints, 0.);
        for (const auto &[cell_facesID_pointsPOS, cell_nPoints] : std::views::zip(this->_cells_facesID_pointsPOS, this->_cells_nPoints))
        {
            for (const auto &faceID_pointsPOS : cell_facesID_pointsPOS)
            {
                cell_nPoints += faceID_pointsPOS.size();
            }
        }
    }

    void Cells::compute_cells_centroid()
    {
        this->_cells_centroid.resize(this->size());
        std::ranges::fill(this->_cells_centroid, std::array<value_t, 3>{});
        assert((this->_cells_facesID_pointsPOS.size() == this->_cells_centroid.size()) && (this->_cells_centroid.size() == this->_cells_nPoints.size()));
        for (const auto &[cell_facesID_pointsPOS, cell_centroid, cell_nPoints] : std::views::zip(this->_cells_facesID_pointsPOS, this->_cells_centroid, this->_cells_nPoints))
        {
            for (const auto &faceID_pointsPOS : cell_facesID_pointsPOS)
            {
                for (const auto &pointsPOS : faceID_pointsPOS)
                {
                    std::ranges::transform(cell_centroid, pointsPOS, cell_centroid.begin(), std::plus<value_t>{});
                }
            }
            std::ranges::for_each(cell_centroid, [&cell_nPoints](value_t &x)
                                  { x /= (double)cell_nPoints; });
        }
    }

    void Cells::compute_cells_volume()
    {
        // for each cell find all their faces
        const auto &facesID_centroid = this->_faces.getFacesCentroid();
        this->_cells_volume.resize(this->size());
        for (const auto &[cell_facesID_pointsPOS, cell_centroid, cell_facesID, cell_volume] : std::views::zip(this->_cells_facesID_pointsPOS, this->_cells_centroid, this->_cells_facesID, this->_cells_volume))
        {
            for (const auto &[pointsPOS, faceID] : std::views::zip(cell_facesID_pointsPOS, cell_facesID))
            {
                size_t face_nPoints = pointsPOS.size();
                std::array<value_t, 3> face_centroid = facesID_centroid[faceID];
                for (size_t i = 0; i < face_nPoints; ++i)
                {
                    const std::array<value_t, 3> &p0 = face_centroid;
                    const std::array<value_t, 3> &p1 = pointsPOS[i];
                    const std::array<value_t, 3> &p2 = pointsPOS[(i + 1) % face_nPoints];
                    const std::array<value_t, 3> &p3 = cell_centroid;

                    const std::array<value_t, 3> b1 = utils::createVec(p3, p0);
                    const std::array<value_t, 3> b2 = utils::createVec(p3, p1);
                    const std::array<value_t, 3> b3 = utils::createVec(p3, p2);

                    cell_volume += fabs(utils::dot_product(b1, utils::cross_product(b2, b3))) / 6.;
                }
            }
        }
    }
   

    // --- Mesh --- //

    Mesh::Mesh(Points _points, Faces _faces, Cells _cells, Boundary _boundary) : cells(std::move(_cells)),
                                                                                 boundary(std::move(_boundary)),
                                                                                 faces(std::move(_faces)),
                                                                                 points(std::move(_points))
    {
        this->registerFacesCellsID();
        this->register_cells_ncellsID();
        this->registerBoundaryFaces();
        this->computeCellToFaceRatio();
    };

     void Mesh::register_cells_ncellsID()
    {
        this->_cells_ncellsID.resize(this->cells.size());
        for (const auto &face_cellsID : this->_faces_cellsID)
        {
            std::vector<size_t> neighboursID = face_cellsID;
            for (const auto &cellID : face_cellsID)
            {
                for (const auto &neighbourID : neighboursID)
                {
                    if (neighbourID != cellID)
                        _cells_ncellsID[cellID].push_back(neighbourID);
                }
            }
        }
    }

    void Mesh::registerFacesCellsID()
    {
        // For each cell (in order), register one self on each of the faces.
        this->_faces_cellsID.resize(this->faces.size());
        this->_faces_cellsCentroid.resize(this->faces.size());

        const size_t N = this->cells.size();
        const auto& cells_facesID = this->cells.getCellsFacesID();
        const auto& cells_centroid = this->cells.getCellsCentroid();
        for (size_t cell_idx = 0; cell_idx < N; ++cell_idx)
        {
            for (const auto &face_ID : cells_facesID[cell_idx])
            {
                this->_faces_cellsID[face_ID].push_back(cell_idx);
                this->_faces_cellsCentroid[face_ID].push_back(cells_centroid[cell_idx]);
            }
        }
    }

    void Mesh::addCellData (const std::vector<value_t>& cell_data, std::string label) {
            _container_map[label] = field_idx;
            ++field_idx;

            if(cell_data.size() != cells.size()) {
                throw std::runtime_error("Cell data doesn't match size of the mesh.");
            }
        
            _container_cells_data.push_back(cell_data);

            // Build the values at the faces as well
            std::vector<value_t> faces_data(faces.size());
            for(const auto& [face_data, face_cellsID, face_fx, face_boundaryFlag]: std::views::zip(faces_data, _faces_cellsID, _faces_fx, _faces_boundaryFlag)) {
                value_t phiP = cell_data[face_cellsID[0]];
                if(face_boundaryFlag) {
                    face_data = phiP;
                    continue;
                }
                value_t phiN = cell_data[face_cellsID[1]];
                face_data = interpolateFaceValue(phiP, phiN, face_fx);
            }
            _container_faces_data.push_back(faces_data);
        }

    void Mesh::registerWallsFacesCellID()
    {
        this->_walls_faces_cellID.resize(this->boundary.size());
        const auto& walls = this->boundary.getWalls();

        for (const auto &[wall, wall_faces_cellID] : std::views::zip(walls, this->_walls_faces_cellID))
        {
            wall_faces_cellID.resize(wall.facesID.size());
            for (const auto &[faceID, face_cellID] : std::views::zip(wall.facesID, wall_faces_cellID))
            {
                // boundary faces should be attached to only one cell
                assert(_faces_cellsID[faceID].size() == 1);
                face_cellID = this->_faces_cellsID[faceID][0];
            }
        }
    }

    void Mesh::registerBoundaryFaces()
    {
        _faces_boundaryFlag.resize(faces.size());
        _cells_boundaryFlag.resize(cells.size());

        std::ranges::fill(_faces_boundaryFlag, false);
        std::ranges::fill(_cells_boundaryFlag, false);


        for (const auto &[face_cellsID, face_boundaryFlag] : std::views::zip(_faces_cellsID, _faces_boundaryFlag))
        {
            // TODO: check if all boundary faces have a boundary condition
            if (face_cellsID.size() == 1)
            {
                face_boundaryFlag = true;
                _cells_boundaryFlag[face_cellsID[0]] = true;
            }
        }
    }

    void Mesh::computeCellToFaceRatio()
    {
        const std::vector<std::array<value_t, 3>> &faces_centroid = faces.getFacesCentroid();
        const std::vector<std::array<value_t, 3>> &cells_centroid = cells.getCellsCentroid();

        _faces_fx.resize(faces.size());
        _faces_delta.resize(faces.size());

        for (const auto &[face_cellsID, face_boundaryFlag, face_fx, face_delta, face_centroid] : std::views::zip(_faces_cellsID, _faces_boundaryFlag, _faces_fx, _faces_delta, faces_centroid))
        {
            const size_t P_cellID = face_cellsID[0];
            const std::array<value_t, 3> &P = cells_centroid[P_cellID];

            if (face_boundaryFlag)
            {
                const std::array<value_t, 3> Pf = utils::createVec(P, face_centroid);
                face_delta = 1. / utils::magnitude(Pf);
                face_fx = 1.;

                continue;
            }

            assert(face_cellsID.size() > 1);

            const size_t N_cellID = face_cellsID[1];

            const std::array<value_t, 3> &N = cells_centroid[N_cellID];

            const std::array<value_t, 3> FN = utils::createVec(face_centroid, N);
            const std::array<value_t, 3> PN = utils::createVec(P, N);
            const double PN_mag_inv = 1. / utils::magnitude(PN);

            face_fx = utils::magnitude(FN) * PN_mag_inv;
            face_delta = PN_mag_inv;
            assert(face_fx >= 0 && face_fx <= 1);
        }
    }

    std::shared_ptr<Mesh> load_mesh(fs::path &mesh_folder_path)
    {
        // Try to open points, faces, cells, and boundary file
        fs::path points_file_path("points");
        fs::path faces_file_path("faces");
        fs::path cells_file_path("cells");
        fs::path boundary_file_path("boundary");

        std::ifstream points_file(mesh_folder_path / points_file_path);
        std::ifstream faces_file(mesh_folder_path / faces_file_path);
        std::ifstream cells_file(mesh_folder_path / cells_file_path);
        std::ifstream boundary_file(mesh_folder_path / boundary_file_path);

        if (points_file.fail())
        {
            std::cout << mesh_folder_path / points_file_path << std::endl;
            throw std::invalid_argument("Indicated point file does not exist.");
        }
        if (faces_file.fail())
        {
            std::cout << mesh_folder_path / faces_file_path << std::endl;
            throw std::invalid_argument("Indicated faces file does not exist.");
        }
        if (cells_file.fail())
        {
            std::cout << mesh_folder_path / cells_file_path << std::endl;
            throw std::invalid_argument("Indicated cells file does not exist.");
        }
        if (boundary_file.fail())
        {
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

    namespace cartesian_mesh
    {

        class CartesianMeshFactory
        {
        public:
            CartesianMeshFactory(std::array<size_t, 3> Ns, std::array<double, 3> DXs) : _Ns(Ns), _DXs(DXs) {};

            std::vector<std::array<value_t, 3>> create_points()
            {
                size_t N = 1;
                for (const auto &Ni : _Ns)
                {
                    N *= (Ni + 1); // Number points on one axis is the number of cells + 1
                }
                std::vector<std::array<value_t, 3>> res(N);

                for (size_t k = 0; k <= _Ns[2]; ++k)
                {
                    for (size_t j = 0; j <= _Ns[1]; ++j)
                    {
                        for (size_t i = 0; i <= _Ns[0]; ++i)
                        {
                            res[ijk(i, j, k)] = std::array<double, 3>({
                                i * _DXs[0],
                                j * _DXs[1],
                                k * _DXs[2],
                            });
                        }
                    }
                }
                return res;
            }

            std::vector<std::vector<size_t>> create_faces();
            std::vector<std::vector<size_t>> create_cells();

        private:
            std::array<double, 3> _DXs;
            std::array<size_t, 3> _Ns;

            std::array<value_t, 3> computePointCoordinate(size_t i, size_t j, size_t k);

            std::function<size_t(size_t, size_t, size_t)> ijk = [&](size_t i, size_t j, size_t k)
            {
                return i + j * (_Ns[0] + 1) + k * (_Ns[0] + 1) * (_Ns[1] + 1);
            };

            std::function<size_t(size_t, size_t, size_t)> cijk = [&](size_t i, size_t j, size_t k)
            {
                return i + j * (_Ns[0]) + k * (_Ns[0]) * (_Ns[1]);
            };

            // Return face id at i-1/2, j, k
            std::function<size_t(size_t, size_t, size_t)> sx = [&](size_t i, size_t j, size_t k)
            {
                return i + j * (_Ns[0]) + k * (_Ns[0]) * (_Ns[1] + 1);
            };

            // Return face id at i, j-1/2, k
            std::function<size_t(size_t, size_t, size_t)> sy = [&](size_t i, size_t j, size_t k)
            {
                return i + j * (_Ns[0] + 1) + k * (_Ns[0] + 1) * (_Ns[1]) + _Ns[0] * (_Ns[1] + 1) * _Ns[2];
            };

            // Return face id at i, j, k-1/2
            std::function<size_t(size_t, size_t, size_t)> sz = [&](size_t i, size_t j, size_t k)
            {
                return i + j * (_Ns[0]) + k * (_Ns[0]) * (_Ns[1]) + _Ns[0] * (_Ns[1] + 1) * _Ns[2] + (_Ns[0] + 1) * _Ns[1] * _Ns[2];
            };
        };

        std::vector<std::vector<size_t>> CartesianMeshFactory::create_faces()
        {
            std::vector<std::vector<size_t>> res;

            // xz axis
            for (size_t k = 0; k < _Ns[2]; ++k)
            {
                size_t j = 0;
                for (size_t i = 0; i < _Ns[0]; ++i)
                    {
                        res.push_back({
                            ijk(i, j, k),
                            ijk(i + 1, j, k),
                            ijk(i + 1, j, k + 1),
                            ijk(i, j, k + 1),
                        });
                    }

                for (size_t j = 1; j <= _Ns[1]; ++j)
                {
                    for (size_t i = 0; i < _Ns[0]; ++i)
                    {
                        res.push_back({
                            ijk(i, j, k + 1),
                            ijk(i + 1, j, k + 1),
                            ijk(i + 1, j, k),
                            ijk(i, j, k),
                        });
                    }
                }
            }

            // yz axis
            for (size_t k = 0; k < _Ns[2]; ++k)
            {
                for (size_t j = 0; j < _Ns[1]; ++j)
                {
                    size_t i = 0;
                    res.push_back({
                            ijk(i, j, k + 1),
                            ijk(i, j + 1, k + 1),
                            ijk(i, j + 1, k),
                            ijk(i, j, k),
                        });

                    for (size_t i = 1; i <= _Ns[0]; ++i)
                    {
                        res.push_back({
                            ijk(i, j, k),
                            ijk(i, j + 1, k),
                            ijk(i, j + 1, k + 1),
                            ijk(i, j, k + 1),
                        });
                    }
                }
            }

            // xy axis
            size_t k = 0;
            for (size_t j = 0; j < _Ns[1]; ++j)
                {
                    for (size_t i = 0; i < _Ns[0]; ++i)
                    {
                        res.push_back({
                            ijk(i, j, k),
                            ijk(i + 1, j, k),
                            ijk(i + 1, j + 1, k),
                            ijk(i, j + 1, k),
                        });
                    }
                }

            for (size_t k = 1; k <= _Ns[2]; ++k)
            {
                for (size_t j = 0; j < _Ns[1]; ++j)
                {
                    for (size_t i = 0; i < _Ns[0]; ++i)
                    {
                        res.push_back({
                            ijk(i, j + 1, k),
                            ijk(i + 1, j + 1, k),
                            ijk(i + 1, j, k),
                            ijk(i, j, k),
                        });
                    }
                }
            }

            return res;
        }

        std::vector<std::vector<size_t>> CartesianMeshFactory::create_cells()
        {
            std::vector<std::vector<size_t>> res;

            for (size_t k = 0; k < this->_Ns[2]; ++k)
            {
                for (size_t j = 0; j < this->_Ns[1]; ++j)
                {
                    for (size_t i = 0; i < this->_Ns[0]; ++i)
                    {
                        res.push_back({sx(i, j, k),
                                       sx(i, j + 1, k),
                                       sy(i, j, k),
                                       sy(i + 1, j, k),
                                       sz(i, j, k),
                                       sz(i, j, k + 1)});
                    }
                }
            }
            return res;
        }
    }

    std::shared_ptr<Mesh> create_cartesian_mesh(std::array<size_t, 3> Ns, std::array<double, 3> DXs)
    {
        std::cout << "Creating cartesian mesh" << std::endl;

        cartesian_mesh::CartesianMeshFactory meshFactory(Ns, DXs);
        std::vector<std::array<value_t, 3>> pointsPOS = meshFactory.create_points();
        // for(const auto& pointPOS: pointsPOS) {
        //     for(const auto& x: pointPOS) std::cout << x << " ";
        //     std::cout << "\n";
        // }
        Points points(pointsPOS);
        std::cout << "Number of points: " << points.size() << std::endl;

        std::vector<std::vector<size_t>> facesPointsID = meshFactory.create_faces();
        // for(const auto& facePointsID: facesPointsID) {
        //     for(const auto& pointID: facePointsID) {
        //         std::cout << pointID << " ";
        //     }
        //     std::cout << "\n";
        // }

        Faces faces(points, facesPointsID);
        std::cout << "Number of faces: " << faces.size() << std::endl;

        // for(const auto& n: faces.getFacesNormalVector()) {
        //     for(const auto& x: n) {
        //         std::cout << x << " ";
        //     }
        //     std::cout << "\n";
        // }

        // for(const auto& faceCentroid: faces.getFacesCentroids()) {
        //     for(const auto& x: faceCentroid) std::cout << x << " ";
        //     std::cout << "\n";
        // }

        std::vector<std::vector<size_t>> cellsFacesID = meshFactory.create_cells();
        Cells cells(cellsFacesID, faces);
        std::cout << "Number of cells: " << cells.size() << std::endl;

        // for(const auto& cellsID: cells.getFacesCellsID()) {
        //     for(const auto& cellID: cellsID) {
        //         std::cout << cellID << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // for(const auto& cellsID: cells.getCellsVolume()) {
        //     std::cout << cellsID << "\n";
        // }
        // std::cout << std::endl;

        Boundary boundary({}, cells);

        return std::make_shared<Mesh>(points, faces, cells, boundary);
    };
}