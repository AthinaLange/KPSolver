/// CSV File IO module.
/** @file
*/

#include <iostream>
#include <fstream>
#include <sstream>

#include "csv.h"


void    read_wave_csv (
                        const std::string&          fn,
                        std::vector<double>&        eta) {

    std::cout << "fn: " << fn << std::endl;

    std::ifstream input_file(fn);
    if (input_file.is_open()) {

        std::string line;
        // no header
        // data
        while (std::getline(input_file, line)) {
            // single value per line
            std::cout << "line: " << line << std::endl;
            double value;
            std::stringstream(line) >> value;
            eta.emplace_back(value);
        }
        input_file.close();
    }
    return;
}


void    write_soliton_csv(
                            const std::string&          fn,
                            int                         nx,
                            int                         ny,
                            const std::vector<double>&  eta) {

    std::ofstream output_file(fn);
    if (output_file.is_open()) {

        for (int j = 0; j < nx*ny; ++j) {
            output_file << eta[j] << std::endl;
        }
        output_file.close();
    }
}

void    write_wave_csv(
                        const std::string&          fn,
                        int                         nx,
                        int                         ny,
                        double                      delx,
                        double                      dely,
                        const std::vector<double>&  eta) {

    std::ofstream output_file(fn);
    if (output_file.is_open()) {

        for (int j = 0; j < nx * ny; ++j) {
            output_file << eta[j] << std::endl;
        }
        /*for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                output_file << i*delx << ","
                            << j*dely << ","
                            << eta[i + j * nx]
                            << std::endl;
            }
            output_file << std::endl;
        }*/
        output_file.close();
    }
}

