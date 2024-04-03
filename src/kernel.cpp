#include <CL/sycl.hpp>

#include"kernel.h"
#include "parameters.h"
#include "trace.h"
#include "maths.h"
#include "particles.h"
#include "meshes.h"
#include "FESolver.h"
bool d_XtoLtet(Particle &part, Volume &volume) {
    /*first try the current tetrahedron*/
    Tetra &tet = volume.elements[part.cell_index];

    bool inside = true;
    /*loop over vertices*/
    for (int i=0;i<4;i++) {
        part.lc[i] = (1.0/6.0)*(tet.alpha[i] - part.pos[0]*tet.beta[i] +
                      part.pos[1]*tet.gamma[i] - part.pos[2]*tet.delta[i])/tet.volume;
        if (part.lc[i]<0 || part.lc[i]>1.0) inside=false;
    }

    if (inside) return true;

    
    /*we are outside the last known tet, find most negative weight*/
    int min_i=0;
    double min_lc=part.lc[0];
    for (int i=1;i<4;i++)
        if (part.lc[i]<min_lc) {min_lc=part.lc[i];min_i=i;}

    /*is there a neighbor in this direction?*/
    if (tet.cell_con[min_i]>=0) {
        part.cell_index = tet.cell_con[min_i];
        return d_XtoLtet2(part,volume);
    }

    return false;
}

bool d_XtoLtet2(Particle &part, Volume &volume) {
    /*first try the current tetrahedron*/
    Tetra &tet = volume.elements[part.cell_index];

    bool inside = true;
    /*loop over vertices*/
    for (int i=0;i<4;i++) {
        part.lc[i] = (1.0/6.0)*(tet.alpha[i] - part.pos[0]*tet.beta[i] +
                      part.pos[1]*tet.gamma[i] - part.pos[2]*tet.delta[i])/tet.volume;
        if (part.lc[i]<0 || part.lc[i]>1.0) inside=false;
    }

    if (inside) return true;

    
    /*we are outside the last known tet, find most negative weight*/
    int min_i=0;
    double min_lc=part.lc[0];
    for (int i=1;i<4;i++)
        if (part.lc[i]<min_lc) {min_lc=part.lc[i];min_i=i;}

    /*is there a neighbor in this direction?*/
    if (tet.cell_con[min_i]>=0) {
        part.cell_index = tet.cell_con[min_i];
        return d_XtoLtet(part,volume);
    }

    return false;
}