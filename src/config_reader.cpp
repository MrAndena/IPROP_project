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

    // Load navier_stokes data

    // Load Navier-Stokes physical parameters
    dati_json.at("navier_stokes").at("physical_parameters").at("viscosity").get_to(data.navier_stokes.physical_parameters.viscosity);
    dati_json.at("navier_stokes").at("physical_parameters").at("gamma").get_to(data.navier_stokes.physical_parameters.gamma);
    dati_json.at("navier_stokes").at("physical_parameters").at("p_over_rho").get_to(data.navier_stokes.physical_parameters.p_over_rho);

    // Load Navier-Stokes numerical parameters
    dati_json.at("navier_stokes").at("numerical_parameters").at("delta_t").get_to(data.navier_stokes.numerical_parameters.delta_t);
    dati_json.at("navier_stokes").at("numerical_parameters").at("end_time").get_to(data.navier_stokes.numerical_parameters.end_time);
    dati_json.at("navier_stokes").at("numerical_parameters").at("output_interval").get_to(data.navier_stokes.numerical_parameters.output_interval);
    dati_json.at("navier_stokes").at("numerical_parameters").at("tolerance").get_to(data.navier_stokes.numerical_parameters.tolerance);
    dati_json.at("navier_stokes").at("numerical_parameters").at("max_iterations").get_to(data.navier_stokes.numerical_parameters.max_iterations);
    dati_json.at("navier_stokes").at("numerical_parameters").at("FE_degree").get_to(data.navier_stokes.numerical_parameters.FE_degree);

    // Load Drift-Diffusion physical parameters
    dati_json.at("drift_diffusion").at("physical_parameters").at("eps_0").get_to(data.drift_diffusion.physical_parameters.eps_0);
    dati_json.at("drift_diffusion").at("physical_parameters").at("eps_r").get_to(data.drift_diffusion.physical_parameters.eps_r);
    dati_json.at("drift_diffusion").at("physical_parameters").at("q0").get_to(data.drift_diffusion.physical_parameters.q0);
    dati_json.at("drift_diffusion").at("physical_parameters").at("A").get_to(data.drift_diffusion.physical_parameters.A);
    dati_json.at("drift_diffusion").at("physical_parameters").at("D").get_to(data.drift_diffusion.physical_parameters.D);
    dati_json.at("drift_diffusion").at("physical_parameters").at("ni").get_to(data.drift_diffusion.physical_parameters.ni);
    dati_json.at("drift_diffusion").at("physical_parameters").at("V_TH").get_to(data.drift_diffusion.physical_parameters.V_TH);
    dati_json.at("drift_diffusion").at("physical_parameters").at("kB").get_to(data.drift_diffusion.physical_parameters.kB);
    dati_json.at("drift_diffusion").at("physical_parameters").at("mu_p").get_to(data.drift_diffusion.physical_parameters.mu_p);
    dati_json.at("drift_diffusion").at("physical_parameters").at("mu_n").get_to(data.drift_diffusion.physical_parameters.mu_n);
    dati_json.at("drift_diffusion").at("physical_parameters").at("mu0").get_to(data.drift_diffusion.physical_parameters.mu0);
    dati_json.at("drift_diffusion").at("physical_parameters").at("stratosphere").get_to(data.drift_diffusion.physical_parameters.stratosphere);
    dati_json.at("drift_diffusion").at("physical_parameters").at("E_ON").get_to(data.drift_diffusion.physical_parameters.E_ON);
    dati_json.at("drift_diffusion").at("physical_parameters").at("E_ref").get_to(data.drift_diffusion.physical_parameters.E_ref);
    dati_json.at("drift_diffusion").at("physical_parameters").at("N_ref").get_to(data.drift_diffusion.physical_parameters.N_ref);
    dati_json.at("drift_diffusion").at("physical_parameters").at("Mm").get_to(data.drift_diffusion.physical_parameters.Mm);
    dati_json.at("drift_diffusion").at("physical_parameters").at("Avo").get_to(data.drift_diffusion.physical_parameters.Avo);
    dati_json.at("drift_diffusion").at("physical_parameters").at("Ve").get_to(data.drift_diffusion.physical_parameters.Ve);

    // Load Drift-Diffusion numerical parameters
    dati_json.at("drift_diffusion").at("numerical_parameters").at("tolerance").get_to(data.drift_diffusion.numerical_parameters.tolerance);
    dati_json.at("drift_diffusion").at("numerical_parameters").at("max_iterations").get_to(data.drift_diffusion.numerical_parameters.max_iterations);
    dati_json.at("drift_diffusion").at("numerical_parameters").at("FE_degree").get_to(data.drift_diffusion.numerical_parameters.FE_degree);

    // Load Geometrical parameters
    dati_json.at("geometrical_parameters").at("NACA_digits").get_to(data.geometrical_parameters.NACA_digits);
    dati_json.at("geometrical_parameters").at("chord_length").get_to(data.geometrical_parameters.chord_length);
    dati_json.at("geometrical_parameters").at("emitter_radius").get_to(data.geometrical_parameters.emitter_radius);
    dati_json.at("geometrical_parameters").at("distance_emitter_collector").get_to(data.geometrical_parameters.distance_emitter_collector);
    dati_json.at("geometrical_parameters").at("distance_trailing_edge_outlet").get_to(data.geometrical_parameters.distance_trailing_edge_outlet);
    dati_json.at("geometrical_parameters").at("distance_emitter_inlet").get_to(data.geometrical_parameters.distance_emitter_inlet);
    dati_json.at("geometrical_parameters").at("distance_emitter_up_bottom").get_to(data.geometrical_parameters.distance_emitter_up_bottom);
    dati_json.at("geometrical_parameters").at("L").get_to(data.geometrical_parameters.L);
    dati_json.at("geometrical_parameters").at("H").get_to(data.geometrical_parameters.H);
    dati_json.at("geometrical_parameters").at("num_elem_L").get_to(data.geometrical_parameters.num_elem_L);
    dati_json.at("geometrical_parameters").at("num_elem_H").get_to(data.geometrical_parameters.num_elem_H);
    
    // Load Simulation Specification
    dati_json.at("simulation_specification").at("ID_simulation").get_to(data.simulation_specification.ID_simulation);



    return;
}