#include "mfem.hpp"
#include <string>
#include<fstream>

#include "io.h"

void write_solution_to_csv(const mfem::GridFunction &u, const mfem::Mesh &mesh, const std::string &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return;
    }

    file << "xi,u\n";  // CSV header

    for (int i = 0; i < u.Size(); i++) {
        double xi = mesh.GetVertex(i)[0];  // Get spatial coordinate
        file << xi << "," << u[i] << "\n";  // Write to CSV
    }

    file.close();
    std::cout << "Solution written to " << filename << std::endl;
}
