/*==============================================================================*
 * MESHES
 *------------------------------------------------------------------------------*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 * Based on `fem-pic.cpp` by Lubos Brieda 
 * See https://www.particleincell.com/2015/fem-pic/ for more information
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef MESHES_H
#define MESHES_H

#include <vector>
#include <string>

#include "particles.h"

/*node type*/
enum NodeType {NORMAL,OPEN,INLET,FIXED};

/*definition of a node*/
struct Node {
    Node(double x, double y, double z) {pos[0]=x;pos[1]=y;pos[2]=z;type=NORMAL;}
    double pos[3];    /*node position*/
    NodeType type;
    double volume;    /*node volume*/
};

/*definition of a tetrahedron*/
struct Tetra {
    int con[4];
    double volume;
    Tetra (int n1, int n2, int n3, int n4) {con[0]=n1;con[1]=n2;con[2]=n3;con[3]=n4;}

    /*data structures to hold precomputed 3x3 determinants*/
    double alpha[4], beta[4], gamma[4], delta[4];

    /*cell connectivity*/
    int cell_con[4];    /*index corresponds to the face opposite the i-th node*/
};

/* Definition of a triangle for the intlet faces*/
struct Face {
    Face(int n1, int n2, int n3) {con[0]=n1, con[1]=n2, con[2]=n3;}
    int con[3];     // IDs of Nodes comprising the face
    double area;
    double u[3];
    double v[3];
    int cell_con;
    double normal[3];
};

/*definition of a volume*/
struct Volume {
    std::vector <Node> nodes;
    std::vector <Tetra> elements;
    std::vector <Face> inlet_faces;
    double avg_edge_len;

    void summarize(std::ostream &out);
};

bool LoadVolumeMesh(const std::string file_name, Volume &volume);
bool LoadSurfaceMesh(const std::string file_name, Volume &volume, NodeType node_type, bool invert_faces);
void OutputMesh(int ts, Volume &volume, double *phi, double (*ef)[3], double *ion_den);

bool XtoLtet(Particle &part, Volume &volume, bool search=true);

/*converts physical coordinate to logical, walking through neighbouring tets
until the one containing the particle is found. Returns true if particle matched to a tet.
Uses a loop rather than recursion and plain arrays so it can also be called from a SYCL kernel*/
inline bool XtoLtet(Particle &part, const Tetra *elements, bool search=true) {
    while (true) {
        const Tetra &tet = elements[part.cell_index];

        bool inside = true;
        /*loop over vertices*/
        for (int i=0;i<4;i++) {
            part.lc[i] = (1.0/6.0)*(tet.alpha[i] - part.pos[0]*tet.beta[i] +
                          part.pos[1]*tet.gamma[i] - part.pos[2]*tet.delta[i])/tet.volume;
            if (part.lc[i]<0 || part.lc[i]>1.0) inside=false;
        }

        if (inside) return true;

        if (!search) return false;
        /*we are outside the last known tet, find most negative weight*/
        int min_i=0;
        double min_lc=part.lc[0];
        for (int i=1;i<4;i++)
            if (part.lc[i]<min_lc) {min_lc=part.lc[i];min_i=i;}

        /*is there a neighbor in this direction? if not the particle has left the domain*/
        if (tet.cell_con[min_i]<0) return false;
        part.cell_index = tet.cell_con[min_i];
    }
}

#endif /* !MESHES_H */
