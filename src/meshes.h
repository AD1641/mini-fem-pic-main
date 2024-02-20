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
void OutputMesh(int ts, Volume &volume, double *phi, double **ef, double *ion_den);

bool XtoLtet(Particle &part, Volume &volume, bool search=true);

#endif /* !MESHES_H */
