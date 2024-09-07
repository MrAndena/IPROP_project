#include "data4geo.hpp"

void write_data_4_geo(const data_struct& d, unsigned short int i){

    // Now we write in a .txt file all the geometry data
    std::ofstream file("Data4Geo.txt");

    if (!file.is_open()) {
        std::cerr << "Failed to open the file for writing."  << std::endl;
        return;
    }
    
    if(d.simulation_specification.ID_simulation <=2){

        file << "NACA_digits = "<<d.geometrical_parameters.NACA_digits<<";"<<std::endl;
        file << "chord_length = " << d.geometrical_parameters.chord_length <<";"<<std::endl;
        file << "emitter_radius = " << d.geometrical_parameters.emitter_radius[i] <<";"<<std::endl;
        file << "distance_emitter_collector = " << d.geometrical_parameters.distance_emitter_collector <<";"<<std::endl;
        file << "distance_trailing_edge_outlet = " << d.geometrical_parameters.distance_trailing_edge_outlet <<";"<<std::endl;
        file << "distance_emitter_inlet = " << d.geometrical_parameters.distance_emitter_inlet <<";"<<std::endl;
        file << "distance_emitter_up_bottom = " << d.geometrical_parameters.distance_emitter_up_bottom <<";"<<std::endl;

    } else{

        file << "L = "<<d.geometrical_parameters.L<<";"<<std::endl;
        file << "H = "<<d.geometrical_parameters.H<<";"<<std::endl;
        file << "num_elem_L = "<<d.geometrical_parameters.num_elem_L<<";"<<std::endl;
        file << "num_elem_H = "<<d.geometrical_parameters.num_elem_H<<";"<<std::endl;

    }


    file.close();

    std::cout << "All the data are written in Data4Geo.txt" << std::endl;


}