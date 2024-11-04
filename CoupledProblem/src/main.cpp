#include "config_reader.hpp"
#include "CompleteProblem.hpp"

int main(int argc, char** argv){ 

//################### - Initialization of parallel utilities - ############################################################################

    using namespace dealii;

    // Initialize the MPI environment for parallel computation.
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    
    // Define an empty data structure to hold the simulation input parameters.
    data_struct my_data; 

    std::cout<<"   Reading JSON file inputs ...";

    // Call the `reader` function to populate `my_data` with values from the input JSON file.
    // This structure is available to all processors, ensuring consistent input across ranks.
    reader(my_data); 

    std::cout<<"   Done !"<<std::endl;


//#################### - Code for the simulation - ##########################################################################################

        
        try
        {
                // Create a distributed triangulation for a 2D domain using MPI for parallel processing.
                // This triangulation object will manage the mesh and handle data distribution across processors.
                parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);

                // Initialize the triangulation by creating the mesh based on `my_data`.
                create_triangulation(tria, my_data);

                // Instantiate the main problem class, `CompleteProblem`, with dimension 2 and pass in
                // the triangulation and input data. This class will handle the simulation.
                CompleteProblem<2> problem(tria, my_data);

                // Run the simulation by invoking the `run` method on the `CompleteProblem` instance.
                // This method contains the main loop or procedure for solving the problem.
                problem.run();

        }
        // Catch any standard exceptions that might occur during simulation execution.
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
        // Catch any non-standard exceptions (unknown types) and print a generic error message.
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
        // Exit successfully after completing the simulation without errors.
        return 0;
}