/*==============================================================================*
 * PARAMETERS
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2023-02-01
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <map>

#include "particles.h"
#include "FESolver.h"

class Parameters{
    public:
        double plasma_den           = 1e20;
        double ion_velocity         = 7000;
        double electron_temperature = 24000;
        int    max_iter             = 250;
        double wall_potential       = 100;
        double dt                   = 1e-7;
        std::map<std::string, std::string> mesh_files;

        bool invert_normals   = false;

        Parameters(std::string fileanme);

        void write(std::ostream &out);

        Species::Name plasma_species = Species::Deuterium;
        FESolver::Method fesolver_method = FESolver::NonLinear;

};

#endif /* !PARAMETERS_H */

