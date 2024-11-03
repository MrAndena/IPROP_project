#include "DriftDiffusion.hpp"



//-----------------------------------------------------------------------------------------------------------------------------------------------------
// HELPER FUNCTION FOR GRID IN 

template <int dim>
void create_triangulation(parallel::distributed::Triangulation<dim> &tria)
{
  
  const std::string filename = "../../meshes/Structured_Square.msh"; //remember: you run from build directory

  ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  std::ifstream input_file(filename);

  if (!input_file.is_open())
  {
    std::cerr << "Error: impossible to open mesh file." << std::endl;
    return;
  }

  pcout <<"Reading the mesh from " << filename << std::endl;
  GridIn<2>  grid_in; //This class implements an input mechanism for grid data. It allows to read a grid structure into a triangulation object
  grid_in.attach_triangulation(tria); //we pass to grid_in our (empty) triangulation
  grid_in.read_msh(input_file); // read the msh file
  pcout << " Grid read correctly "<< std::endl;

}
  

//######################### MAIN #############################################################################################################


int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);
      create_triangulation(tria);
      
      

      DriftDiffusion<2> dd_2d(tria);
      dd_2d.run(20, 1e-12);
    
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
