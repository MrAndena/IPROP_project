#include "gmsh_runner.hpp"


void gmsh_runner(const data_struct& d, unsigned short int i){
    
    std::string gmsh_command;


    if(d.simulation_specification.ID_simulation<=2){
    
        std::string name_mesh = "naca_" + std::to_string(i) + ".msh";
        //std::string name_mesh = "WireWire_" + std::to_string(i) + ".msh";
        //gmsh_command = "gmsh -2 -format msh2 ../config/GeoMeshFile.geo -o ../output/meshes/"+name_mesh;           
        gmsh_command = "gmsh -2 -format msh2 ../config/GeoMeshFile_original.geo -o ../output/meshes/"+name_mesh; 
        //gmsh_command = "gmsh -2 -format msh2 ../config/coarse_WW.geo -o ../output/meshes/"+name_mesh;

    } else{

        std::string name_mesh = "Structured_Square.msh";
        gmsh_command = "gmsh -2 -format msh2 ../config/PNjunction.geo -o ../output/meshes/"+name_mesh;

    } 
       
    int result = std::system(gmsh_command.c_str());

    if (result == 0) {
        std::cout << "Gmsh executed successfully."<<std::endl;
    } else {
        std::cerr << "Gmsh execution failed with error code. " << std::endl;
    }

   

    return;
}