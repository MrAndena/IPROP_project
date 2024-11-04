#ifndef DATA_STRUCT_HPP
#define DATA_STRUCT_HPP

#include<string>
// Fluid Physical Parameters Struct

struct FluidParameters {

    double viscosity;     // viscosity of the fluid 
    double gamma;         // stabilizaton parameter in NS

};

// Electrical Physical Parameters Struct
struct ElectricalParameters { 

    double eps_0;    //permittivity of free space [F/m]= [C^2 s^2 / kg / m^3]
    double eps_r;    //permittivity [F/m]= [C^2 s^2 / kg / m^3]
    double q0;       //unit charge [C]
    double kB;       //Boltzman constant[J/K]
    double mu0;      //Mobility of ions
    bool stratosphere; //bool, if false use atmospheric 0 km condition
    double E_ON;     // onset field threshold [V/m]
    double E_ref;    // maximum field value   [V/m]
    double N_ref;    // maximum density value [m^-3]
    double N_min;    
    double Mm;       // average air molar mass [kg m^-3]
    double Ve;       // emitter voltage [V]

};


// -- Geometrical infos Struct --

struct NACA {
     
    double chord_length;                   // length of the chord [m]
    int naca_digits;                       // type of naca
    double emitter_radius;                 // radius of the circular emitter [m]
    double distance_emitter_collector;     // distance circular emitter surface, airfoil collector surface [m]
    double distance_trailing_edge_outlet;  // distance trailing edge outlet [m]
    double distance_emitter_inlet;         // distance circular emitter surface and inlet [m]
    double distance_emitter_up_bottom;     // distance circular emitter surface and up/bottom boundary [m]

};


struct WW {
     
    double emitter_radius;                 // radius of the circular emitter [m]
    double collector_radius;               // radius of the circular collector [m]
    double distance_emitter_collector;     // distance circular emitter surface, airfoil collector surface [m]
    double distance_collector_outlet;      // distance trailing edge outlet [m]
    double distance_emitter_inlet;         // distance circular emitter surface and inlet [m]
    double distance_emitter_up_bottom;     // distance circular emitter surface and up/bottom boundary [m]

};

struct CYL {

    double emitter_radius;                 // radius of the circular emitter [m]
    double collector_radius;               // radius of the circular collector [m]

};

struct GeometricalParameters {
   
   double emitter_center_X;               // X coordinate of the emitter center [m]
   double emitter_center_Y;               // Y coordinate of the emitter center [m]
   NACA naca;
   WW ww;
   CYL cyl;

};

// -- Simulations infos Struct --

struct SimulationSpecification {
   
   std::string mesh_name;                // name of the file containing the mesh
   std::string mesh_TAG;                 // type of the mesh

};


// Compleate struct
struct data_struct {
    
    FluidParameters fluid_parameters;
    ElectricalParameters electrical_parameters;
    GeometricalParameters geometrical_parameters;
    SimulationSpecification simulation_specification;

};

#endif //DATA_STRUCT_HPP