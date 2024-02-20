/*==============================================================================*
 * PARAMETERS
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2023-02-01
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#include <fstream>
#include <iostream>
#include <string>
#include <regex>

#include "parameters.h"

Parameters::Parameters(std::string filename) {
    std::ifstream params_file;
    params_file.open(filename);

    std::string line;
    while(std::getline(params_file, line)) {
        // Filter out the comments and skip blank lines
        std::regex_replace(line, std::regex("\\s*\\(#.*\\)\\?$"), "");
        if (line.length() == 0) continue;

        std::smatch matches;
        if (std::regex_match(line, matches, std::regex("(\\w+)\\s*=\\s*(.*)$"))) {
            if (matches.size() != 3) {
                std::cerr << "Invalid input line: '" << line << "'" << std::endl;
                exit(-1);
            }

            std::string  key = matches[1];
            std::string  value = matches[2];

            if (key == "plasma_den") plasma_den = std::stod(value);
            else if (key == "ion_velocity") ion_velocity = std::stod(value);
            else if (key == "electron_temperature") electron_temperature = std::stod(value);
            else if (key == "max_iter") max_iter = std::stoi(value);
            else if (key == "wall_potential") wall_potential = std::stod(value);
            else if (key == "dt") dt = std::stod(value);

            else if (key == "global_mesh" or key == "inlet_mesh" or key == "wall_mesh") {
                mesh_files.insert({key, value});
            }

            else if (key == "invert_normals") {
                if (std::regex_match(value, std::regex("^t(rue)?$", std::regex_constants::icase))
                    || value == "1") {
                    invert_normals = true;
                } else if (std::regex_match(value, std::regex("^f(alse)?$", std::regex_constants::icase))
                           || value == "0") {
                    invert_normals = false;
                }
            }

            else if (key == "plasma_species") {
                if (std::regex_match(value, std::regex("D(uterium)?"))) plasma_species = Species::Deuterium;
                else if (std::regex_match(value, std::regex("O(xygen)?"))) plasma_species = Species::Oxygen;
                else {
                    std::cerr << "Invalid plasma species: '" << value << "'" << std::endl;
                    exit(-1);
                }
            }
            
            else if (key == "fesolver_method") {
                if      (std::regex_match(value, std::regex("nonlinear", std::regex_constants::icase))) fesolver_method = FESolver::NonLinear;
                else if (std::regex_match(value, std::regex("gaussseidel", std::regex_constants::icase))) fesolver_method = FESolver::GaussSeidel;
                else if (std::regex_match(value, std::regex("lapack", std::regex_constants::icase))) fesolver_method = FESolver::Lapack;
                else {
                    std::cerr << "Invalid fesolver_method: '" << value << "'" << std::endl;
                    exit(-1);
                }
            }
            
            else {
                std::cerr << "Unrecognised parameter: '" << key << "'" << std::endl;
                exit(-1);
            }
        }

    }
}


void Parameters::write(std::ostream &out) {

    out << "SIMULATION PARAMETERS"  << std::endl << "---------------------" << std::endl;
    out << "  plasma_den = " << plasma_den  << std::endl;
    out << "  ion_velocity = " << ion_velocity << std::endl;
    out << "  electron_temperature = " << electron_temperature << std::endl;
    out << "  max_iter = " << max_iter << std::endl;
    out << "  wall_potential = " << wall_potential << std::endl;
    out << "  dt = " << dt << std::endl;
    out << "  mesh_files:" << std::endl;
    for (auto const &[key, val]: mesh_files) {
        out << "    " << key << " = " << val << std::endl;
    }
    out << "  invert_normals = ";
    out << (invert_normals ? "True" : "False") << std::endl;

    out << "  plasma_species = ";
    switch (plasma_species) {
        case Species::Deuterium:
            out << "Deuterium";
            break;
        case Species::Oxygen:
            out << "Oxygen";
            break;
    }
    out << std::endl;

    out << "  fesolver_method = ";
    switch (fesolver_method) {
        case FESolver::NonLinear:
            out << "NonLinear";
            break;
        case FESolver::GaussSeidel:
            out << "GaussSeidel";
            break;
        case FESolver::Lapack:
            out << "Lapack";
            break;
    }
    out << std::endl << std::endl;
}
