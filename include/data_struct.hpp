#ifndef DATA_STRUCT_HPP
#define DATA_STRUCT_HPP

#include<string>
// Fluid Physical Parameters Struct

struct FluidParameters {

    double viscosity;     // viscosity NS
    double gamma;         // gamma NS
    double p_over_rho;    // Boundary condition at fluid outlet

};

// Electrical Physical Parameters Struct
struct ElectricalParameters { // verificare e finire di scirvere unita misura

    double eps_0;    //permittivity of free space[F/m]= [C^2 s^2 / kg / m^3]
    double eps_r;    //permittivity
    double q0;       //unit charge [C]
    double kB;       //[J/K]
    double mu0;      //Mobility
    bool stratosphere; //bool, if false use atmospheric 0 km condition
    double E_ON;     // onset field threshold [V/m]
    double E_ref;    // maximum field value [V/m]
    double N_ref;    // maximum density value [m^-3]
    double N_min;
    double Mm;       // average air molar mass [kg m^-3]
    double Avo;      // Avogadro's number
    double Ve;       // emitter voltage [V]
    double theta;

};


// -- Geometrical infos Struct --

struct NACA {
     
    double chord_length;                   // length of the chord [m]
    int naca_digits;                       // type of naca
    double emitter_radius;                 // radius of the circular emitter [m]
    double distance_emitter_collector;     // distance circular emitter surface, airfoil collector surface
    double distance_trailing_edge_outlet;  // distance trailing edge outlet
    double distance_emitter_inlet;         // distance circular emitter surface and inlet
    double distance_emitter_up_bottom;     // distance circular emitter surface and up/bottom boundary

};


struct WW {
     
    double emitter_radius;                 // radius of the circular emitter [m]
    double collector_radius;               // radius of the circular collector [m]
    double distance_emitter_collector;     // distance circular emitter surface, airfoil collector surface
    double distance_collector_outlet;      // distance trailing edge outlet
    double distance_emitter_inlet;         // distance circular emitter surface and inlet
    double distance_emitter_up_bottom;     // distance circular emitter surface and up/bottom boundary

};

struct CYL {

    double emitter_radius;                 // radius of the circular emitter [m]
    double collector_radius;               // radius of the circular collector [m]

};

struct GeometricalParameters {

   NACA naca;
   WW ww;
   CYL cyl;

};

// -- Simulations infos Struct --

struct SimulationSpecification {
   
   std::string mesh_name;
   std::string mesh_TAG;

};


// Compleate struct
struct data_struct {
    
    FluidParameters fluid_parameters;
    ElectricalParameters electrical_parameters;
    GeometricalParameters geometrical_parameters;
    SimulationSpecification simulation_specification;

};

#endif //DATA_STRUCT_HPP