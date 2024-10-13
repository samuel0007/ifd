#include "include/ifd_lib/mesh.hpp"
#include <vector>
#include <string_view>
#include <stdexcept>
#include <filesystem>

namespace fs = std::filesystem;

struct Global_Options {
    fs::path mesh_folder_path;
};

Global_Options parse_args(int argc, char *argv[]) {
    Global_Options options;
    const std::vector<std::string_view> args(argv + 1, argv + argc);
    if (argc < 2) {
        throw std::invalid_argument("User has to pass mesh data.");
    }
    options.mesh_folder_path = fs::path(args[0]);
    return options;
};


int main(int argc, char *argv[]) {
    // Load data
    Global_Options options = parse_args(argc, argv);
    std::shared_ptr<mesh::Mesh> mesh_p = mesh::load_mesh(options.mesh_folder_path);
    mesh_p->points.print();
    mesh_p->faces.print();
    mesh_p->cells.print();
    mesh_p->boundary.print();

    return 0;
}