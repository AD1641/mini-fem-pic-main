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
#include <sycl/sycl.hpp>
#include <time.h>

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

/*Copies of the data that MoveParticles needs on the SYCL device. The mesh never
changes so it is copied over once, the particles, electric field and ion density
are copied to and from the device every time step*/
struct DeviceData {
    sycl::queue q;
    int n_nodes;
    int n_elements;
    Tetra *elements;            /*mesh elements*/
    double (*ef)[3];            /*electric field in each element*/
    double *den;                /*ion density on each node*/
    Particle *particles = nullptr;
    size_t max_particles = 0;   /*number of particles there is room for*/

    DeviceData(Volume &volume);
    ~DeviceData();
    DeviceData(const DeviceData &) = delete;
    void reserveParticles(size_t n);
};

/*PROTOTYPES*/
void write_header(std::ostream& out);
void write_footer(std::ostream& out);
int InjectIons(Species &ions, Volume &volume, FESolver &solver, Parameters params);
void MoveParticles(Species &ions, Volume &volume, FESolver &solver, const Parameters &params, DeviceData &device);


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



    /*sample assembly code*/
    solver.startAssembly();
    solver.preAssembly();    /*this will form K and F0*/

    solver.summarize(std::cout);

    /*ions species*/
    Species ions(n_nodes, params.plasma_species);
    ions.charge=1*QE;
    ions.mass = 16*AMU;
    ions.spwt = 2e2;

    /*set up the SYCL device used to move the particles*/
    DeviceData device(volume);
    std::cout << "SYCL device: " << device.q.get_device().get_info<sycl::info::device::name>() << std::endl << std::endl;

    /*main loop*/
    int ts;
    for (ts=0;ts<params.max_iter;ts++) {
        struct timespec start, stop; 
        double duration; 

        
        /*sample new particles*/
        int n_new_particles = InjectIons(ions, volume, solver, params);

        int old_nparts = ions.particles.size();
        /*update velocity and move particles*/
        clock_gettime(CLOCK_MONOTONIC, &start);
        MoveParticles(ions, volume, solver, params, device);
        clock_gettime(CLOCK_MONOTONIC, &stop); 
        duration = (double) (stop.tv_sec * 1000000000 + stop.tv_nsec) - (start.tv_sec * 1000000000 + start.tv_nsec);
        duration = duration/1000000000;

        /*call potential solver*/
        solver.computePhi(ions.den, params.fesolver_method);

        solver.updateEf();

        /*check values*/
        double max_den=0;
        for (int n=0;n<n_nodes;n++) if (ions.den[n]>max_den) max_den=ions.den[n];
        
    
        double max_phi=0;


        for (int n=0;n<n_nodes;n++) if (abs(solver.uh[n])>max_phi) max_phi=abs(solver.uh[n]);


        if ((ts+1)%10==0) OutputMesh(ts,volume, solver.uh, solver.ef, ions.den);

        

        std::cout<<"ts: "<<ts
                 <<"\t np: "<<ions.particles.size()
                 <<" (" <<  n_new_particles << " added, "<< old_nparts - ions.particles.size() << " removed)" 
                 <<"\t max den: "<<max_den
                 <<"\t max |phi|: "<<max_phi
                 <<"\t step duration: " << duration
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

/*updates ion velocities and positions on the SYCL device*/
void MoveParticles(Species &ions, Volume &volume, FESolver &solver, const Parameters &params, DeviceData &device) {    TRACE_ME;
    int n_nodes = (int) volume.nodes.size();
    size_t n_parts = ions.particles.size();
    sycl::queue &q = device.q;

    /*copy the particles and electric field to the device and reset the ion density*/
    device.reserveParticles(n_parts);
    q.memcpy(device.particles, ions.particles.data(), n_parts*sizeof(Particle));
    q.memcpy(device.ef, solver.ef, device.n_elements*sizeof(double[3]));
    q.fill(device.den, 0.0, n_nodes);

    /*the kernel can only use plain values and device pointers, not the host objects*/
    Particle *particles = device.particles;
    const Tetra *elements = device.elements;
    const double (*ef)[3] = device.ef;
    double *den = device.den;
    double q_over_m = ions.charge/ions.mass;
    double dt = params.dt;

    q.parallel_for(sycl::range<1>(n_parts), [=](sycl::id<1> p) {
        Particle &part = particles[p];

        /*update particle velocity*/
        for (int i=0;i<3;i++)
            part.vel[i] += q_over_m*ef[part.cell_index][i]*dt;

        /*update particle positions*/
        for (int i=0;i<3;i++) part.pos[i]+=part.vel[i]*dt;

        if (XtoLtet(part, elements)) {
            const Tetra &tet = elements[part.cell_index];
            /*now we know that we are inside this tetrahedron, scatter*/
            for (int v=0;v<4;v++) {
                sycl::atomic_ref<double, sycl::memory_order::relaxed, sycl::memory_scope::device,
                                 sycl::access::address_space::global_space> node_den(den[tet.con[v]]);
                node_den.fetch_add(part.lc[v]);
            }
        } else {
            part.cell_index = -1;    /*particle has left the domain, flag it for removal*/
        }
    });

    /*copy the particles and ion density back*/
    q.memcpy(ions.particles.data(), particles, n_parts*sizeof(Particle));
    q.memcpy(ions.den, den, n_nodes*sizeof(double));
    q.wait();

    /*remove the particles that left the domain*/
    std::erase_if(ions.particles, [](const Particle &part) {return part.cell_index<0;});

    /*convert to ion density*/
    for (int n=0;n<n_nodes;n++) ions.den[n] *= ions.spwt/volume.nodes[n].volume;
}

DeviceData::DeviceData(Volume &volume): q(sycl::property::queue::in_order{}) {
    sycl::device dev = q.get_device();
    if (!dev.has(sycl::aspect::fp64) || !dev.has(sycl::aspect::atomic64)) {
        std::cerr<<"SYCL device "<<dev.get_info<sycl::info::device::name>()<<" does not support double precision"<<std::endl;
        exit(-1);
    }

    n_nodes = volume.nodes.size();
    n_elements = volume.elements.size();
    elements = sycl::malloc_device<Tetra>(n_elements, q);
    ef = sycl::malloc_device<double[3]>(n_elements, q);
    den = sycl::malloc_device<double>(n_nodes, q);
    if (!elements || !ef || !den) {std::cerr<<"Failed to allocate the mesh on the SYCL device"<<std::endl;exit(-1);}

    /*the mesh doesn't change, so only needs copying once*/
    q.memcpy(elements, volume.elements.data(), n_elements*sizeof(Tetra)).wait();
}

DeviceData::~DeviceData() {
    sycl::free(elements, q);
    sycl::free(ef, q);
    sycl::free(den, q);
    if (particles) sycl::free(particles, q);
}

/*makes sure there is room for n particles on the device*/
void DeviceData::reserveParticles(size_t n) {
    if (n <= max_particles) return;
    if (particles) sycl::free(particles, q);

    max_particles = 2*n;    /*leave room to grow so we don't reallocate every time step*/
    particles = sycl::malloc_device<Particle>(max_particles, q);
    if (!particles) {std::cerr<<"Failed to allocate "<<max_particles<<" particles on the SYCL device"<<std::endl;exit(-1);}
}
