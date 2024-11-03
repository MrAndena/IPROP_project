#include "config_reader.hpp"
#include "CompleteProblem.hpp"

int main(int argc, char** argv){ 

//################### - Initialization of parallel utilities - ############################################################################

    using namespace dealii;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    //const unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD); // store rank of the current processor
    
    data_struct my_data; //empty struct data

    std::cout<<"   Reading JSON file inputs ...";

    reader(my_data); //all the processor know my_data, struct that contains all the user specification

    std::cout<<"   Done !"<<std::endl;


//#################### - Code for the simulation - ##########################################################################################

        
        try
        {

                parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);

                create_triangulation(tria, my_data);

                drift_diffusion<2> problem(tria, my_data);

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

        
    return 0;
}