/*==============================================================================*
 * PARTICLES
 *------------------------------------------------------------------------------*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 * Based on `fem-pic.cpp` by Lubos Brieda 
 * See https://www.particleincell.com/2015/fem-pic/ for more information
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <string>

/*particle*/
struct Particle {
    double pos[3];
    double vel[3];
    double lc[4];    /*particle's weights*/
    int cell_index;    /*last cell known to contain this particle*/
};

/*species class*/
class Species {
public:
    enum Name {Oxygen, Deuterium};
    std::vector<Particle> particles;
    double *den;
    double spwt;
    double mass;
    double charge;
    double rem;        /*fractional particle reminder*/

    Species (int n_nodes);
    Species (int n_nodes, Name species);
    ~Species () {delete[] den;}
};


void OutputParticles(std::vector<Particle> &particles);

#endif /* !PARTICLES_H */
