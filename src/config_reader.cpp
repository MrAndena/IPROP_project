#include "config_reader.hpp"

void reader(data_struct& data){

    // Open a file stream for reading the json file
    std::ifstream inFile("../config/configuration.json");

    // Check if the file stream is open
    if (!inFile.is_open()) {
    std::cerr << "Failed to open the file for reading." << std::endl;
    return ;
    }

    // Read JSON data from the file
    json dati_json;       
    inFile >> dati_json;

    // Close the file stream
    inFile.close();

    // Access the data from the JSON object and store them in data_struct object

    // Load the physical fluid parameters
    dati_json.at("fluid_parameters").at("viscosity").get_to(data.fluid_parameters.viscosity);
    dati_json.at("fluid_parameters").at("gamma").get_to(data.fluid_parameters.gamma);
    dati_json.at("fluid_parameters").at("p_over_rho").get_to(data.fluid_parameters.p_over_rho);

    // Load the electrical parameters
    dati_json.at("electrical_parameters").at("eps_0").get_to(data.electrical_parameters.eps_0);
    dati_json.at("electrical_parameters").at("eps_r").get_to(data.electrical_parameters.eps_r);
    dati_json.at("electrical_parameters").at("q0").get_to(data.electrical_parameters.q0);
    dati_json.at("electrical_parameters").at("eps_0").get_to(data.electrical_parameters.eps_0);
    dati_json.at("electrical_parameters").at("kB").get_to(data.electrical_parameters.kB);
    dati_json.at("electrical_parameters").at("mu0").get_to(data.electrical_parameters.mu0);
    dati_json.at("electrical_parameters").at("stratosphere").get_to(data.electrical_parameters.stratosphere);
    dati_json.at("electrical_parameters").at("E_ON").get_to(data.electrical_parameters.E_ON);
    dati_json.at("electrical_parameters").at("E_ref").get_to(data.electrical_parameters.E_ref);
    dati_json.at("electrical_parameters").at("N_ref").get_to(data.electrical_parameters.N_ref);
    dati_json.at("electrical_parameters").at("N_min").get_to(data.electrical_parameters.N_min);
    dati_json.at("electrical_parameters").at("Mm").get_to(data.electrical_parameters.Mm);
    dati_json.at("electrical_parameters").at("Avo").get_to(data.electrical_parameters.Avo);
    dati_json.at("electrical_parameters").at("Ve").get_to(data.electrical_parameters.Ve);
    dati_json.at("electrical_parameters").at("theta").get_to(data.electrical_parameters.theta);

    // Load Geometrical parameters

    //NACA
    dati_json.at("geometrical_parameters").at("NACA").at("chord_length").get_to(data.geometrical_parameters.naca.chord_length);
    dati_json.at("geometrical_parameters").at("NACA").at("naca_digits").get_to(data.geometrical_parameters.naca.naca_digits);
    dati_json.at("geometrical_parameters").at("NACA").at("emitter_radius").get_to(data.geometrical_parameters.naca.emitter_radius);
    dati_json.at("geometrical_parameters").at("NACA").at("distance_emitter_collector").get_to(data.geometrical_parameters.naca.distance_emitter_collector);
    dati_json.at("geometrical_parameters").at("NACA").at("distance_trailing_edge_outlet").get_to(data.geometrical_parameters.naca.distance_trailing_edge_outlet);
    dati_json.at("geometrical_parameters").at("NACA").at("distance_emitter_inlet").get_to(data.geometrical_parameters.naca.distance_emitter_inlet);
    dati_json.at("geometrical_parameters").at("NACA").at("distance_emitter_up_bottom").get_to(data.geometrical_parameters.naca.distance_emitter_up_bottom);
    
    //WW
    dati_json.at("geometrical_parameters").at("WW").at("emitter_radius").get_to(data.geometrical_parameters.ww.emitter_radius);
    dati_json.at("geometrical_parameters").at("WW").at("collector_radius").get_to(data.geometrical_parameters.ww.collector_radius);
    dati_json.at("geometrical_parameters").at("WW").at("distance_emitter_collector").get_to(data.geometrical_parameters.ww.distance_emitter_collector);
    dati_json.at("geometrical_parameters").at("WW").at("distance_collector_outlet").get_to(data.geometrical_parameters.ww.distance_collector_outlet);
    dati_json.at("geometrical_parameters").at("WW").at("distance_emitter_inlet").get_to(data.geometrical_parameters.ww.distance_emitter_inlet);
    dati_json.at("geometrical_parameters").at("WW").at("distance_emitter_up_bottom").get_to(data.geometrical_parameters.ww.distance_emitter_up_bottom);

    //CYL
    dati_json.at("geometrical_parameters").at("CYL").at("emitter_radius").get_to(data.geometrical_parameters.cyl.emitter_radius);
    dati_json.at("geometrical_parameters").at("CYL").at("collector_radius").get_to(data.geometrical_parameters.cyl.collector_radius);
    

    // Load Simulation Specification
    dati_json.at("simulation_specification").at("mesh_name").get_to(data.simulation_specification.mesh_name);
    dati_json.at("simulation_specification").at("mesh_TAG").get_to(data.simulation_specification.mesh_TAG);



    return;
}