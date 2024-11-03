#ifndef IFD_VTKWRITER_HPP
#define IFD_VTKWRITER_HPP

#include <filesystem>

namespace fs = std::filesystem;

class VTKWriter
{
public:
    VTKWriter(const fs::path &folderpath, std::shared_ptr<mesh::Mesh> mesh_p) : _filepath(folderpath), _mesh_p(mesh_p) {};

    void write_mesh(std::ofstream &of)
    {
        auto pointsPOS = _mesh_p->points.getPointsPOS();
        of << "# vtk DataFile Version 2.0 \n";
        of << "test\n";
        of << "ASCII\n";
        of << "DATASET UNSTRUCTURED_GRID\n";
        of << "POINTS " << pointsPOS.size() << " float\n";
        for (const auto &pointPOS : pointsPOS)
        {
            of << pointPOS[0] << " " << pointPOS[1] << " " << pointPOS[2] << "\n";
        }
        auto cellsPointsID = _mesh_p->cells.getCellsPointsID();

        size_t totalPoints = 0;
        for (const auto &pointsID : cellsPointsID)
        {
            totalPoints += pointsID.size() + 1; // number of points per cell + leading size integer
        }
        of << "CELLS " << cellsPointsID.size() << " " << totalPoints << "\n";
        for (const auto &pointsID : cellsPointsID)
        {
            of << pointsID.size() << " ";
            for (const auto &pointID : pointsID)
            {
                of << pointID << " ";
            }
            of << "\n";
        }
        of << "CELL_TYPES " << cellsPointsID.size() << "\n";
        for (const auto &pointsID : cellsPointsID)
        {
            of << 11 << " "; // ASSUME voxels...............
        }
        of << "\n";
    }

    void void_write_cell_data(std::ofstream &of, const std::vector<double> &cell_data, std::string label)
    {
        assert(cell_data.size() == _mesh_p->cells.size());
        of << "CELL_DATA " << cell_data.size() << "\n";
        of << "SCALARS " << label << " float 1\n";
        of << "LOOKUP_TABLE default \n";
        for (const auto &el : cell_data)
        {
            of << el << " ";
        }
        of << std::endl;
    }

    void write(const std::vector<double> &cell_data, std::string filename, std::string label = "")
    {
        fs::path labelpath(filename + _extension);
        std::ofstream of;

        of.open(_filepath / labelpath);
        write_mesh(of);
        void_write_cell_data(of, cell_data, label);
        of.close();
    }

private:
    fs::path _filepath;
    std::shared_ptr<mesh::Mesh> _mesh_p;
    std::string _extension = ".vtk";
};

#endif // IFD_VTKWRITER_HPP
