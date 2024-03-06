/* Example of a Finite Element ES-PIC Code

  Written by Lubos Brieda for Advanced PIC 2015 Lesson 8
  See https://www.particleincell.com/2015/fem-pic/ for more information

  To compile and run:
    g++ -std=c++10 -O2 fem-pic.cpp -o fem-pic
    ./fem-pic

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <CL/sycl.hpp>

#include "parameters.h"
#include "trace.h"
#include "maths.h"
#include "particles.h"
#include "meshes.h"
#include "FESolver.h"

/*constants*/

int trace::enabled = 1;
Trace trace::current = Trace("__TRACE_BASE__");

std::time_t start_time;

/*PROTOTYPES*/
void write_header(std::ostream& out);
void write_footer(std::ostream& out);
int InjectIons(Species &ions, Volume &volume, FESolver &solver, Parameters params);
void MoveParticles(Species &ions, Volume &volume, FESolver &solver, Parameters params);


/**************** MAIN **************************/
int main(int argc, char **argv) {

    // Read in the simulation parameters from the file provided
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        exit(0);
    }
    std::string params_file = argv[1];

    write_header(std::cout);


    Parameters params(params_file);
    params.write(std::cout);

    // Load the meshes defined in the parameters file
    Volume volume;
    if (!LoadVolumeMesh(params.mesh_files["global_mesh"],volume) ||
        !LoadSurfaceMesh(params.mesh_files["inlet_mesh"],volume,INLET, params.invert_normals) ||
        !LoadSurfaceMesh(params.mesh_files["wall_mesh"],volume,FIXED, params.invert_normals)) return -1;

    volume.summarize(std::cout);

    /*instantiate solver*/
    FESolver solver(volume);

    /*set reference paramaters*/
    solver.phi0 = 0;
    solver.n0 = params.plasma_den;
    solver.kTe = Kb * params.electron_temperature;

    int n_nodes = volume.nodes.size();


    /*initialize solver "g" array*/
    for (int n=0;n<n_nodes;n++) {
        if (volume.nodes[n].type==INLET) solver.g[n]=0;    /*phi_inlet*/
        else if (volume.nodes[n].type==FIXED) solver.g[n]=-params.wall_potential; /*fixed phi points*/
        else solver.g[n]=0;    /*default*/
    }

    // sycl::queue q(sycl::default_selector{});

    // sycl::buffer<double,1> d_solver_g {solver.g.data(), sycl::range<1>(solver.g.size())};
    // sycl::buffer<double,1> d_params_wall_pot {params.wall_potential, sycl::range<1>(params.wall_potential.size())};

    // q.submit ([&](sycl::handler& h)){
    //     auto G = d_solver_g.template get_access<sycl::access::mode::write>(h);
    //     auto P = d_params_wall_pot.template get_access<sycl::access::mode::read>(h);

    //     h.parallel_for(sycl::<1>{n_nodes,[=] (sycl::id<1> i)})
    //     {
    //         if(volume.nodes[i].type == INLET)
    //         {
    //             G[i] = 0;
    //         }
    //         else if (volume.nodes[n].type==FIXED)
    //         {
    //             G[i] =- P;
    //         }
    //         else
    //         {
    //             G[i] = 0;
    //         }
    //     }
    // }


    /*sample assembly code*/
    solver.startAssembly();
    solver.preAssembly();    /*this will form K and F0*/

    solver.summarize(std::cout);

    /*ions species*/
    Species ions(n_nodes, params.plasma_species);
    ions.charge=1*QE;
    ions.mass = 16*AMU;
    ions.spwt = 2e2;

    /*main loop*/
    int ts;
    for (ts=0;ts<params.max_iter;ts++) {
        /*sample new particles*/
        int n_new_particles = InjectIons(ions, volume, solver, params);

        int old_nparts = ions.particles.size();
        /*update velocity and move particles*/
        MoveParticles(ions, volume, solver, params);

        /*call potential solver*/
        solver.computePhi(ions.den, params.fesolver_method);

        solver.updateEf();

        /*check values*/
        double max_den=0;
        //for (int n=0;n<n_nodes;n++) if (ions.den[n]>max_den) max_den=ions.den[n];
        {
            sycl::queue q(sycl::default_selector{});

            sycl::buffer<double,1> d_ions_den {ions.den, sycl::range<1>(sizeof(ions.den))};
            sycl::buffer<double,1> d_max_den {max_den, sycl::range<1>(sizeof(max_den))};

            q.submit ([&](sycl::handler& h){
                auto ION_DEN = d_ions_den.template get_access<sycl::access::mode::read>(h);
                auto MAX_DEN = d_max_den.template get_access<sycl::access::mode::write>(h);

                h.parallel_for(sycl::range<1>{n_nodes},[=] (sycl::id<1> i)
                {
                    if(ION_DEN[i] > MAX_DEN)
                    {
                        MAX_DEN = ION_DEN[i];
                    }
                });
            }).wait();
        }
     
    
        double max_phi=0;


        for (int n=0;n<n_nodes;n++) if (abs(solver.uh[n])>max_phi) max_phi=abs(solver.uh[n]);


        if ((ts+1)%10==0) OutputMesh(ts,volume, solver.uh, solver.ef, ions.den);

        std::cout<<"ts: "<<ts
                 <<"\t np: "<<ions.particles.size()
                 <<" (" <<  n_new_particles << " added, "<< old_nparts - ions.particles.size() << " removed)" 
                 <<"\t max den: "<<max_den
                 <<"\t max |phi|: "<<max_phi
                 <<std::endl;
    }

    /*output mesh*/
    OutputMesh(ts,volume, solver.uh, solver.ef, ions.den);

    /*output particles*/
    OutputParticles(ions.particles);
    if (trace::enabled) trace::current.write_profile("trace.csv");

    write_footer(std::cout);
    return 0;
}


/*** FUNCTIONS ***/


void write_header(std::ostream &out) {
    out << "MINI-PIC" << std::endl << "========" << std::endl;

    start_time = std::time(nullptr);
    std::tm tm = *std::localtime(&start_time);
    out << "Simulation started at " << std::put_time(&tm, "%Y-%m-%d, %T") << std::endl;
#pragma omp master
    {
        out << "Running on " << omp_get_max_threads() << " threads" << std::endl;
    }
    out << std::endl;
}

void write_footer(std::ostream &out) {
    std::time_t end_time = std::time(nullptr);
    std::tm tm = *std::localtime(&end_time);
    out << "Simulation finished at " << std::put_time(&tm, "%Y-%m-%d, %T") << std::endl;
    out << "Calculation took " << trace::current.calculation_time() << " seconds" << std::endl;
}


/*** Samples particle on the inlet surface
here we just sample on a known plane. A generic code should instead sample
from the surface triangles making up the inlet face*/
int InjectIons(Species &ions, Volume &volume, FESolver &solver, Parameters params) { TRACE_ME;
    /*set area of the k=0 face, this should be the sum of triangle areas on the inlet*/

    int total_new_particles = 0;
    for (auto face: volume.inlet_faces) {

        /*number of real ions per sec, given prescribed density and velocity*/
        double num_per_sec = params.plasma_den*params.ion_velocity*face.area;

        /*number of ions to generate in this time step*/
        double num_real = num_per_sec*params.dt;

        /*fraction number of macroparticles*/
        double fnum_mp = num_real/ions.spwt + ions.rem;

        /*integer number of macroparticles*/
        int num_mp = (int)fnum_mp;

        /*update reminder*/
        ions.rem = fnum_mp-num_mp;

        /*sample particles*/
        for (int p=0;p<num_mp;p++) {
            /*new particle*/
            Particle part;

            /*sample random position on the inlet face*/
            double a = rnd();
            double b = rnd();
            if ((a+b) > 1)  {
                a = 1-a;
                b = 1-b;
            }

            for (int i=0; i<3; i++) {
                part.pos[i] = a*face.u[i] + b*face.v[i] + volume.nodes[face.con[0]].pos[i];

                /*injecting cold beam*/
                part.vel[i] = face.normal[i] * params.ion_velocity;
            }

            /*set initial tetrahedron*/
            part.cell_index = face.cell_con;

            /*rewind velocity*/
            double ef_part[3];
            solver.evalEf(ef_part, part.cell_index);

            for (int i=0;i<3;i++)
                part.vel[i] -= ions.charge/ions.mass*ef_part[i]*(0.5*params.dt);

            /*add to list*/
            ions.particles.push_back(part);
        }
        total_new_particles += num_mp;
    }
    return total_new_particles;
}

/*updates ion velocities and positions*/
void MoveParticles(Species &ions, Volume &volume, FESolver &solver, Parameters params) {    TRACE_ME;
    int n_nodes = (int) volume.nodes.size();

    /*reset ion density*/
    for (int i=0;i<n_nodes;i++) ions.den[i] = 0;

    #pragma omp parallel
    {
        std::vector<Particle> thread_newparts;
        #pragma omp for
        for (auto part_it = ions.particles.begin(); part_it != ions.particles.end(); part_it++) {
            Particle &part = *part_it;

            /*update particle velocity*/
            double ef_part[3];
            solver.evalEf(ef_part, part.cell_index);

            for (int i=0;i<3;i++)
                part.vel[i] += ions.charge/ions.mass*ef_part[i]*params.dt;

            /*update particle positions*/
            for (int i=0;i<3;i++) part.pos[i]+=part.vel[i]*params.dt;

            //trace::current.enter("XtoLtet");
            bool inside = XtoLtet(part,volume);
            //trace::current.exit("XtoLtet");

            if (inside) {
                Tetra &tet = volume.elements[part.cell_index];
                /*now we know that we are inside this tetrahedron, scatter*/
                double sum=0;
                for (int v=0;v<4;v++) {
                    #pragma omp atomic update
                    ions.den[tet.con[v]]+=part.lc[v];
                    sum+=part.lc[v];    /*for testing*/
                }

                /*testing*/
                if (std::abs(sum-1.0)>0.001) std::cout<<sum<<std::endl;

                thread_newparts.push_back(part);
            }
        }

        #pragma omp master
        {
            ions.particles.clear();
        }
        #pragma omp barrier
        #pragma omp critical
        {
            ions.particles.insert(ions.particles.end(), thread_newparts.begin(), thread_newparts.end());
        }
    }

    /*convert to ion density*/
    for (int n=0;n<n_nodes;n++) ions.den[n] *= ions.spwt/volume.nodes[n].volume;
}

