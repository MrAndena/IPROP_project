#include "config_reader.hpp"
#include "data4geo.hpp"
#include "NACA_points_generator.hpp"
#include "gmsh_runner.hpp"
#include "InsIMEX.hpp"
// #include "CollectorGeometryNS.hpp"
#include "PoissonProblem.h"
#include "Complete_problem.hpp"

int main(int argc, char** argv){

//################### - Initialization of parallel utilities - ############################################################################

    using namespace dealii;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    const unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD); // store rank of the current processor

    data_struct my_data; //empty struct data
    reader(my_data); //all the processor know my_data, struct that contains all the user specification

    unsigned short int number_of_simulations = my_data.geometrical_parameters.emitter_radius.size(); //3

//################# - Create all the meshes for the simulations - ##########################################################################
    
    if(rank==0){ // only rank zero create the meshes
        
        if(my_data.simulation_specification.ID_simulation <= 2){

                NACA_points_generator my_NACA_object(my_data); //generate the points of the desired profile
                my_NACA_object.write_profile();                //write the points in a .txt file

                for(unsigned short int i = 0; i <number_of_simulations; ++i){
                
                        write_data_4_geo(my_data, i); // create the i-th "Data4Geo.txt"
                        
                        gmsh_runner(my_data,i);// run gmsh to create the i-th .msh mesh file for the simulation
                                
                }


        }else {
                
                //per adesso solo una simulazione
                write_data_4_geo(my_data,0);

                gmsh_runner(my_data,0);

        }
    
    }

//####################  Switch to choose the desired simulation (by ID_simulation) - ################################################################################

        switch (my_data.simulation_specification.ID_simulation) {

                case 1:{ // NS on airfoil profile
                        
                        MPI_Barrier(MPI_COMM_WORLD); // Needed in order to wait rank zero !!

                        for(unsigned short int i=0; i<number_of_simulations;++i){
                        
                                        //code for the simulation
                                        try
                                        {
                
                                        parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);

                                        create_triangulation(tria, my_data, i);

                                        InsIMEX<2> flow(tria, my_data.navier_stokes, i);

                                        flow.run();
                                        }
                                        catch (std::exception &exc)
                                        {
                                        std::cerr << std::endl
                                                << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;
                                        std::cerr << "Exception on processing: " << std::endl
                                                << exc.what() << std::endl
                                                << "Aborting!" << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;
                                        return 1;
                                        }
                                        catch (...)
                                        {
                                        std::cerr << std::endl
                                                << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;    
                                        std::cerr << "Unknown exception!" << std::endl
                                                << "Aborting!" << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;
                                        return 1;
                                        }
                        }
                        
                        break;

                }
        
                

                case 2:{ // NS-DD on naca profile

                        MPI_Barrier(MPI_COMM_WORLD); // Needed in order to wait rank zero !!

                        for(unsigned short int i=0; i<number_of_simulations;++i){
                        
                                        //code for the simulation
                                        try
                                        {
                
                                        parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);

                                        create_triangulation(tria, my_data, i);

                                        CompleteProblem<2> problem(tria, my_data, i);

                                        problem.run();
                                        }
                                        catch (std::exception &exc)
                                        {
                                        std::cerr << std::endl
                                                << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;
                                        std::cerr << "Exception on processing: " << std::endl
                                                << exc.what() << std::endl
                                                << "Aborting!" << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;
                                        return 1;
                                        }
                                        catch (...)
                                        {
                                        std::cerr << std::endl
                                                << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;    
                                        std::cerr << "Unknown exception!" << std::endl
                                                << "Aborting!" << std::endl
                                                << "----------------------------------------------------"
                                                << std::endl;
                                        return 1;
                                        }
                        }

                
                        break;
                }


                case 3:{ // NL Poisson over pn junction
                         
                        MPI_Barrier(MPI_COMM_WORLD); // Needed in order to wait rank zero !!

                        try
                        {

                        parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);
                        create_triangulation(tria, my_data,0);                         
                        

                        PoissonProblem<2> poisson_problem_2d(tria,my_data);
                        poisson_problem_2d.run(my_data.drift_diffusion.numerical_parameters.max_iterations,
                                               my_data.drift_diffusion.numerical_parameters.tolerance);
                        
                        }
                        catch (std::exception &exc)
                        {
                        std::cerr << std::endl
                                        << std::endl
                                        << "----------------------------------------------------"
                                        << std::endl;
                        std::cerr << "Exception on processing: " << std::endl
                                        << exc.what() << std::endl
                                        << "Aborting!" << std::endl
                                        << "----------------------------------------------------"
                                        << std::endl;

                        return 1;
                        }
                        catch (...)
                        {
                        std::cerr << std::endl
                                        << std::endl
                                        << "----------------------------------------------------"
                                        << std::endl;
                        std::cerr << "Unknown exception!" << std::endl
                                        << "Aborting!" << std::endl
                                        << "----------------------------------------------------"
                                        << std::endl;
                        return 1;
                        }

                        break;
                }


                case 4:{ // DD over pn junction

                        break;
                }

                default:{

                        std::cout << "Invalid option!" << std::endl;

                        break;
                }
        }
    


    return 0;
}