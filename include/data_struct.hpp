#ifndef DATA_STRUCT_HPP
#define DATA_STRUCT_HPP

#include<vector>
// -- NS Struct --

// NS - Physical Parameters Struct
struct NavierStokesPhysicalParameters {

    double viscosity;     // viscosity NS
    double gamma;         // gamma NS
    double p_over_rho;    // Boundary condition at fluid outlet

};

// NS - Numerical Parameters Struct
struct NavierStokesNumericalParameters {

    double delta_t;
    double end_time;
    double output_interval;
    double tolerance;
    double max_iterations;
    int FE_degree;

};

// NS - Final Struct
struct NavierStokes {

    NavierStokesPhysicalParameters physical_parameters;
    NavierStokesNumericalParameters numerical_parameters;

};


// -- DD Struct --

// DD - Physical Parameters Struct
struct DriftDiffusionPhysicalParameters { // verificare e finire di scirvere unita misura

    double eps_0;    //permittivity of free space[F/m]= [C^2 s^2 / kg / m^3]
    double eps_r;    //permittivity
    double q0;       //unit charge [C]
    double A;        //pn junction drogaggio
    double D;        //pn junction drogaggio
    double ni;       //pn junction qualcosa
    double V_TH;     //pn junction ion temperature [V]
    double kB;       //[J/K]
    double mu_p;     //pn junction mobility p
    double mu_n;     //pn junction mobility n
    double mu0;      //Moseley (?) 
    bool stratosphere; //bool, if false use atmospheric 0 km condition
    double E_ON;     // onset field threshold [V/m]
    double E_ref;    // maximum field value [V/m]
    double N_ref;    // maximum density value [m^-3]
    double Mm;       // average air molar mass [kg m^-3]
    double Avo;      // Avogadro's number
    double Ve;       // emitter voltage [V]

};

// DD - Numerical Parameters Struct
struct DriftDiffusionNumericalParameters {

    double tolerance;
    int max_iterations;
    int FE_degree;

};

// DD - Final Struct
struct DriftDiffusion {

    DriftDiffusionPhysicalParameters physical_parameters;
    DriftDiffusionNumericalParameters numerical_parameters;

};


// -- Geometrical infos Struct --

struct GeometricalParameters {

    // for NACA mesh

    int NACA_digits;                       // last two digits of the NACA airfoil
    double chord_length;                   // length of the chord [m]
    std::vector<double> emitter_radius;    // radius of the circular emitter [m]
    double distance_emitter_collector;     // distance circular emitter surface, airfoil collector surface
    double distance_trailing_edge_outlet;  // distance trailing edge outlet
    double distance_emitter_inlet;         // distance circular emitter surface and inlet
    double distance_emitter_up_bottom;     // distance circular emitter surface and up/bottom boundary

    // for PNjunction mesh

    double L;       // length pn junciton
    double H;       // height pn junction
    int num_elem_L; // number of elements along the length of the junction
    int num_elem_H; // number of elements along the height of the junction

};

// -- Simulations infos Struct --

struct SimulationSpecification {
   int ID_simulation;  // ID for changing the simulation: (1) NS,  (2) NS-DD,  (3) NL Poisson,  (4) DD

};


// Compleate struct
struct data_struct {

    NavierStokes navier_stokes;
    DriftDiffusion drift_diffusion;
    GeometricalParameters geometrical_parameters;
    SimulationSpecification simulation_specification;

};

#endif //DATA_STRUCT_HPP