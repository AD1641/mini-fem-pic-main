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

#include <vector>
#include <iostream> 

#include "trace.h"
#include "particles.h"
#include "FESolver.h"

/*saves particle data*/
void OutputParticles(std::vector<Particle> &particles) { TRACE_ME;
    std::ofstream out("particles.vtp");
    if (!out.is_open()) {std::cerr<<"Failed to open output file "<<std::endl;exit(-1);}

    /*header*/
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<PolyData>\n";
    out<<"<Piece NumberOfPoints=\""<<particles.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
    out<<"NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

    /*points*/
    out<<"<Points>\n";
    out<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (Particle &part: particles)
        out<<part.pos[0]<<" "<<part.pos[1]<<" "<<part.pos[2]<<"\n";
    out<<"</DataArray>\n";
    out<<"</Points>\n";

    out<<"</Piece>\n";
    out<<"</PolyData>\n";
    out<<"</VTKFile>\n";

    out.close();
}

Species::Species(int n_nodes) {
    den=new double[n_nodes];
    rem=0;
}

Species::Species(int n_nodes, Species::Name n) {
    den=new double[n_nodes];
    rem=0;
    switch (n) {
        case Species::Deuterium:
            charge=1*QE;
            mass = 2*AMU;
            spwt = 2e2;
            break;

        case Species::Oxygen:
            charge=1*QE;
            mass = 16*AMU;
            spwt = 2e2;
            break;

        default:
            break;

    }

}
