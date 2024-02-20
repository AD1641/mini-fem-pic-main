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


#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <algorithm>


#include "trace.h"
#include "maths.h"
#include "particles.h"
#include "meshes.h"

/*loads and initializes volume mesh*/
bool LoadVolumeMesh(const std::string file_name, Volume &volume) { TRACE_ME;
    /*open file*/
    std::ifstream in(file_name);
    if (!in.is_open()) {std::cerr<<"Failed to open "<<file_name<<std::endl; return false;}

    /*read number of nodes and elements*/
    int n_nodes, n_elements;
    in>>n_nodes>>n_elements;

    /*read the nodes*/
    for (int n=0;n<n_nodes;n++) {
        int index;
        double x, y, z;

        in >> index >> x >> y >> z;
        if (index!=n+1) std::cerr<<"Inconsistent node numbering"<<std::endl;

        volume.nodes.emplace_back(x/1000.,y/1000.,z/1000.);
    }

    std::vector<double> edge_lengths;
    /*read elements, this will also contain edges and triangles*/
    for (int e=0;e<n_elements;e++) {
        int index, type;
        int n1, n2, n3, n4;

        in >> index >> type;

        if (type==304) { // Volume element
            in >> n1 >> n2 >> n3 >> n4;
            /*flipping nodes 2 & 3 to get positive volumes*/
            volume.elements.emplace_back(n1-1, n2-1, n3-1, n4-1);

        } else if (type == 102) { // Edge element
            in >> n1 >> n2;
            double len = 0;
            for (int d=0; d<3; d++) {
                len += pow(volume.nodes[n1].pos[d] - volume.nodes[n2].pos[d], 2);
            }
            len = sqrt(len);
            edge_lengths.emplace_back(len);

        } else {
            std::string s; getline(in,s);continue;
        }
    }
    volume.avg_edge_len = 0;
    for (auto l: edge_lengths) volume.avg_edge_len += l;
    volume.avg_edge_len = volume.avg_edge_len / edge_lengths.size();



    /*reset number of nodes and elements since we skipped bunch of lines and triangles*/
    n_nodes = volume.nodes.size();
    n_elements = volume.elements.size();

    std::map<int, std::vector<int>> node_con;

    /*compute element volumes*/
    for (int l=0; l<n_elements; l++) {
        Tetra &tet = volume.elements[l];
        double M[4][4];

        /*set first column to 1*/
        for (int i=0;i<4;i++) M[i][0] = 1;

        /*loop over vertices*/
        for (int v=0;v<4;v++) {
            for (int dim=0;dim<3;dim++) {
                M[v][dim+1] = volume.nodes[tet.con[v]].pos[dim];
            }

            node_con[tet.con[v]].emplace_back(l);
        }

        /*volume is (1/6)*det4(M)*/
        tet.volume = (1.0/6.0)*det4(M);

        /*flip ABCD to ADBC if negative volume*/
        if (tet.volume<0) {int t=tet.con[1];tet.con[1]=tet.con[3];tet.con[3]=t;tet.volume=-tet.volume;}
    }

    /*precompute 3x3 determinants for LC computation*/
    for (Tetra &tet:volume.elements) {
        double M[3][3];
        /*loop over vertices*/
        for (int v=0;v<4;v++) {
            int v2,v3,v4;

            switch (v) {
                case 0: v2=1;v3=2;v4=3;break;
                case 1: v2=3;v3=2;v4=0;break;
                case 2: v2=3;v3=0;v4=1;break;
                case 3: v2=1;v3=0;v4=2;break;
            }

            double *p2 = volume.nodes[tet.con[v2]].pos;
            double *p3 = volume.nodes[tet.con[v3]].pos;
            double *p4 = volume.nodes[tet.con[v4]].pos;

            /*alpha*/
            M[0][0] = p2[0];
            M[0][1] = p2[1];
            M[0][2] = p2[2];
            M[1][0] = p3[0];
            M[1][1] = p3[1];
            M[1][2] = p3[2];
            M[2][0] = p4[0];
            M[2][1] = p4[1];
            M[2][2] = p4[2];
            tet.alpha[v] = det3(M);

            /*beta*/
            M[0][0] =1;
            M[0][1] = p2[1];
            M[0][2] = p2[2];
            M[1][0] = 1;
            M[1][1] = p3[1];
            M[1][2] = p3[2];
            M[2][0] = 1;
            M[2][1] = p4[1];
            M[2][2] = p4[2];
            tet.beta[v] = det3(M);

            /*gamma*/
            M[0][0] =1;
            M[0][1] = p2[0];
            M[0][2] = p2[2];
            M[1][0] = 1;
            M[1][1] = p3[0];
            M[1][2] = p3[2];
            M[2][0] = 1;
            M[2][1] = p4[0];
            M[2][2] = p4[2];
            tet.gamma[v] = det3(M);

            /*delta*/
            M[0][0] =1;
            M[0][1] = p2[0];
            M[0][2] = p2[1];
            M[1][0] = 1;
            M[1][1] = p3[0];
            M[1][2] = p3[1];
            M[2][0] = 1;
            M[2][1] = p4[0];
            M[2][2] = p4[1];
            tet.delta[v] = det3(M);
        }
    }

    /*build cell connectivity, there is probably a faster way*/

    /*reset connectivities*/
    for (int l=0;l<n_elements;l++) {
        Tetra &tet = volume.elements[l];
        for (int v=0;v<4;v++) tet.cell_con[v] = -1;    /*no neighbor*/
    }


    for (int l=0;l<n_elements;l++) {
        Tetra &tet = volume.elements[l];
        int v1,v2,v3;
        for (int v=0;v<4;v++) {
            /*skip if already set*/
            if (tet.cell_con[v]>=0) continue;

            switch(v) {
                case 0: v1=1;v2=2;v3=3;break;
                case 1: v1=2;v2=3;v3=0;break;
                case 2: v1=3;v2=0;v3=1;break;
                case 3: v1=0;v2=1;v3=2;break;
            }

            /*loop over the tets again looking for one with these three vertices*/
            /*(ejh - only look at elements that share a node with this element)*/
            for (int n=0; n<4; n++) {
                for (int m:node_con[tet.con[n]]) {
                    if (l == m) continue;
                    Tetra &other = volume.elements[m];

                    bool matches[4] = {false,false,false,false};
                    int count = 0;
                    for (int k=0;k<4;k++) {
                        if (other.con[k]==tet.con[v1] ||
                            other.con[k]==tet.con[v2] ||
                            other.con[k]==tet.con[v3]) {count++;matches[k]=true;}
                    }

                    /*if three vertices match*/
                    if (count==3) {
                        tet.cell_con[v] = m;

                        /*set the cell connectivity for the index without a matching vertex to l*/
                        for (int k=0;k<4;k++)
                            if(!matches[k]) other.cell_con[k] = l;
                    }
                }
            }
        }
    }

    /*also compute node volumes by scattering cell volumes,this can only be done after 3x3 dets are computed*/

    /*first set all to zero*/
    for (Node &node:volume.nodes) {node.volume=0;}

    for (int i=0;i<n_elements;i++) {
        Particle dummy_part;
        Tetra &tet = volume.elements[i];
        dummy_part.cell_index = i;
        /*compute centroid position*/
        for (int dim=0;dim<3;dim++) {
            dummy_part.pos[dim]=0;
            for (int v=0;v<4;v++) dummy_part.pos[dim]+=0.25*volume.nodes[tet.con[v]].pos[dim];
        }

        bool found = XtoLtet(dummy_part,volume,false);
        if (!found) std::cerr<<"something is wrong"<<std::endl;

        for (int v=0;v<4;v++) {
            volume.nodes[tet.con[v]].volume += dummy_part.lc[v]*tet.volume;
        }

    }

    /*mark nodes on open faces as open*/
    for (size_t e=0;e<volume.elements.size();e++) {
        Tetra &tet = volume.elements[e];
        for (int v=0;v<4;v++)
            if (tet.cell_con[v]<0)    /*no neighbor*/ {
                for (int i=0;i<4;i++) {
                    if (i!=v) volume.nodes[tet.con[i]].type=OPEN;
                }
            }
    }

    return true;
}

/*loads nodes from a surface mesh file and sets them to the specified node type*/
bool LoadSurfaceMesh(const std::string file_name, Volume &volume, NodeType node_type, bool invert_normal) { TRACE_ME;
    /*open file*/
    std::ifstream in(file_name);
    if (!in.is_open()) {std::cerr<<"Failed to open "<<file_name<<std::endl; return false;}

    /*read number of nodes and elements*/
    int n_nodes, n_elements;
    in>>n_nodes>>n_elements;

    int nn = volume.nodes.size();

    /*read the nodes*/
    for (int n=0;n<n_nodes;n++) {
        int index;
        double x, y, z;

        in >> index >> x >> y >> z;

        if (index<1 || index>nn) {std::cerr<<"Incorrect node number "<<index<<std::endl;continue;}
        volume.nodes[index-1].type=node_type;
    }

    if (node_type == INLET) {
        for (int e=0;e<n_elements;e++) {
            int index, type;
            int n1, n2, n3;

            in >> index >> type;

            if (type!=203) {std::string s; getline(in,s);continue;}

            in >> n1 >> n2 >> n3;

            /*flipping nodes 2 & 3 to get positive volumes*/
            volume.inlet_faces.emplace_back(n1-1, n2-1, n3-1);

            volume.inlet_faces.back().cell_con = -1;

            // Find the volume element that attaches to the inlet surface
            for (size_t v=0;v<volume.elements.size(); v++) {
                int matching_nodes = 0;
                for (int element_node: volume.elements[v].con) {
                    if (element_node == n1-1 || element_node == n2-1 || element_node == n3-1) {
                        matching_nodes += 1;
                    }
                }
                if (matching_nodes == 3) {
                    if (volume.inlet_faces.back().cell_con == -1) {
                         volume.inlet_faces.back().cell_con = v;
                    } else {
                        std::cerr<<"Inlet surface attached to more than one volume element"<<index<<std::endl;
                        exit(-1);
                    }
                }
            }
            if ( volume.inlet_faces.back().cell_con == -1) {
                std::cerr<<"No volume element attached to inlet surface"<<index<<std::endl;
                exit(-1);
            }

            // Set the inlet velocity normal to the inlet surface
            double normal[3];
            for (int i=0; i<3; i++) {
                volume.inlet_faces.back().u[i] = volume.nodes[n2-1].pos[i] - volume.nodes[n1-1].pos[i];
                volume.inlet_faces.back().v[i] = volume.nodes[n3-1].pos[i] - volume.nodes[n1-1].pos[i];
            }

            normal[0] = volume.inlet_faces.back().u[2]*volume.inlet_faces.back().v[1] 
                          - volume.inlet_faces.back().u[1]*volume.inlet_faces.back().v[2];
            normal[1] = volume.inlet_faces.back().u[0]*volume.inlet_faces.back().v[2] 
                          - volume.inlet_faces.back().u[2]*volume.inlet_faces.back().v[0];
            normal[2] = volume.inlet_faces.back().u[1]*volume.inlet_faces.back().v[0] 
                          - volume.inlet_faces.back().u[0]*volume.inlet_faces.back().v[1];
            if (invert_normal) {
                normal[0] = -normal[0];
                normal[1] = -normal[1];
                normal[2] = -normal[2];
            }

            double normal_len = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

            for (int i=0; i<3; i++) {
                volume.inlet_faces.back().normal[i] = normal[i]/normal_len;
            }
            volume.inlet_faces.back().area = normal_len / 2;
        }
    }

    return true;
}

/*saves volume mesh*/
void OutputMesh(int ts, Volume &volume, double *phi, double **ef, double *ion_den) { TRACE_ME;
    std::stringstream ss;
    ss<<"mesh_"<<std::setfill('0')<<std::setw(4)<<ts+1<<".vtu";
    std::ofstream out(ss.str());
    if (!out.is_open()) {std::cerr<<"Failed to open file "<<ss.str()<<std::endl;exit(-1);}

    /*header*/
    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out<<"<UnstructuredGrid>\n";
    out<<"<Piece NumberOfPoints=\""<<volume.nodes.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
    out<<"NumberOfStrips=\"0\" NumberOfCells=\""<<volume.elements.size()<<"\">\n";

    /*points*/
    out<<"<Points>\n";
    out<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (Node &node: volume.nodes)
        out<<node.pos[0]<<" "<<node.pos[1]<<" "<<node.pos[2]<<"\n";
    out<<"</DataArray>\n";
    out<<"</Points>\n";

    /*Cells*/
    out<<"<Cells>\n";
    out<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (Tetra &tetra: volume.elements)
        out<<tetra.con[0]<<" "<<tetra.con[1]<<" "<<tetra.con[2]<<" "<<tetra.con[3]<<"\n";
    out<<"</DataArray>\n";

    out<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t e=0; e<volume.elements.size();e++)
        out<<(e+1)*4<<" ";
    out<<"\n";
    out<<"</DataArray>\n";

    out<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t e=0; e<volume.elements.size();e++)
        out<<"10 ";
    out<<"\n";
    out<<"</DataArray>\n";
    out<<"</Cells>\n";

    /*save point data*/
    out<<"<PointData Scalars=\"phi\">\n";
    out<<"<DataArray type=\"Int32\" Name=\"node_index\" format=\"ascii\">\n";
    for (size_t n=0; n<volume.nodes.size();n++)
        out<<n<<" ";
    out<<"\n";
    out<<"</DataArray>\n";
    out<<"<DataArray type=\"Int32\" Name=\"node_type\" format=\"ascii\">\n";
    for (size_t n=0; n<volume.nodes.size();n++)
        out<<volume.nodes[n].type<<" ";
    out<<"\n";
    out<<"</DataArray>\n";

    out<<"<DataArray type=\"Float32\" Name=\"phi\" format=\"ascii\">\n";
    for (size_t n=0; n<volume.nodes.size();n++)
        out<<phi[n]<<" ";
    out<<"\n";
    out<<"</DataArray>\n";

    out<<"<DataArray type=\"Float32\" Name=\"ion_den\" format=\"ascii\">\n";
    for (size_t n=0; n<volume.nodes.size();n++)
        out<<ion_den[n]<<" ";
    out<<"\n";
    out<<"</DataArray>\n";

    out<<"</PointData>\n";

    /*save cell data*/
    out<<"<CellData Vectors=\"ef\">\n";
    out<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"ef\" format=\"ascii\">\n";
    for (size_t e=0; e<volume.elements.size();e++)
        out<<ef[e][0]<<" "<<ef[e][1]<<" "<<ef[e][2]<<" ";
    out<<"\n";
    out<<"</DataArray>\n";

    out<<"<DataArray type=\"Float32\" Name=\"cell_volume\" format=\"ascii\">\n";
    for (Tetra &tet:volume.elements)
        out<<tet.volume<<" ";
    out<<"\n";
    out<<"</DataArray>\n";

    out<<"</CellData>\n";

    out<<"</Piece>\n";
    out<<"</UnstructuredGrid>\n";
    out<<"</VTKFile>\n";

    out.close();
}

/*converts physical coordinate to logical
Returns true if particle matched to a tet
*/
bool XtoLtet(Particle &part, Volume &volume, bool search) {
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

    if (!search) return false;
    /*we are outside the last known tet, find most negative weight*/
    int min_i=0;
    double min_lc=part.lc[0];
    for (int i=1;i<4;i++)
        if (part.lc[i]<min_lc) {min_lc=part.lc[i];min_i=i;}

    /*is there a neighbor in this direction?*/
    if (tet.cell_con[min_i]>=0) {
        part.cell_index = tet.cell_con[min_i];
        return XtoLtet(part,volume);
    }

    return false;
}

void Volume::summarize(std::ostream &out) {
    out << "MESH INFORMATION" << std::endl << "----------------" << std::endl;
    out << "  Number of nodes          = " << nodes.size() << std::endl;
    out << "  Number of elements       = " << elements.size() << std::endl;
    out << "  Number of inlet faces    = " << inlet_faces.size() << std::endl;
    out << std::endl;

    double min_pos[3] = {99999,99999,99999};
    double max_pos[3] = {0,0,0};
    for (auto node: nodes) {
        for (int d=0; d<3; d++) {
            min_pos[d] = std::min(node.pos[d], min_pos[d]);
            max_pos[d] = std::max(node.pos[d], max_pos[d]);
        }
    }

    out << "  Simulation region: "
        << max_pos[0] - min_pos[0] << "m x "
        << max_pos[1] - min_pos[1] << "m x "
        << max_pos[2] - min_pos[2] << "m" << std::endl;
    out << "  Average edge length: " << avg_edge_len <<"m"<< std::endl;

    out << std::endl;
}
