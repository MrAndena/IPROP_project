#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>

// This file includes UMFPACK: the direct solver:
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_in.h>             
#include <deal.II/base/geometry_info.h>        
#include <deal.II/grid/manifold_lib.h>           
#include <deal.II/base/mpi.h>                    

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "InsIMEX.hpp"
#include "MyDataStruct.hpp"
#include "json.hpp"

using json = nlohmann::json;

using namespace dealii;


// Collector Manifold - START

double get_collector_height(const double &p)
{
  
  double collector_length = 1.;
  const double x = p/collector_length;
	double y = 0;

	if ( abs(x-1.) > 1e-12 && abs(x) > 1e-12 ) {
		double a0 = 0.2969;
		double a1 = -0.126;
		double a2 = -0.3516;
		double a3 = 0.2843;
		double a4 = -0.1036; // or -0.1015 for an open trailing edge
		double t = 12;    // Last 2 digits of the NACA by 100

		y = 5*t*( a0 * sqrt(x) + a1 * x + a2 * pow(x,2.0) + a3 * pow(x,3.0) + a4 * pow(x,4.0) );
	}

	return y * collector_length;
}



template <int dim>
class CollectorGeometry : public ChartManifold<dim, dim, dim-1>     //ChartManifold is a class describes mappings that can be expressed in terms of charts.
  {
public:
  virtual Point<dim-1> pull_back(const Point<dim> &space_point) const override;        //Pull back the given point in spacedim to the Euclidean chartdim dimensional space

  virtual Point<dim> push_forward(const Point<dim-1> &chart_point) const override;     //Given a point in the chartdim dimensional Euclidean space, this method returns a point on the manifold embedded in the spacedim Euclidean space.
    
  virtual std::unique_ptr<Manifold<dim, dim>> clone() const override;                  //Return a copy of this manifold

  };

template <int dim>
std::unique_ptr<Manifold<dim, dim>> CollectorGeometry<dim>::clone() const
  {
  return std::make_unique<CollectorGeometry<dim>>();
  }


template <int dim>
Point<dim> CollectorGeometry<dim>::push_forward(const Point<dim-1>  &x) const          //Input: a chart point that in our case is a 1D point 
{
  const double y = get_collector_height(x[0]);

  Point<dim> p;
  p[0] = x[0]; p[1] = y;

  if (dim == 3) {
  p[2] = x[1];
  }

  return p;                                                                              //Output: a point of our collector in 2D 
}


template <int dim>
Point<dim-1>  CollectorGeometry<dim>::pull_back(const Point<dim> &p) const             //Input: a point in our 2D mesh
{
  Point<dim-1> x;
  x[0] = p[0];

  if (dim == 3) {
  x[1] = p[2];
  }

  return x;                                                                              //Output: a chart point that in our case is a 1D point
}  
// Collector Manifold - END

 


 void create_triangulation(parallel::distributed::Triangulation<2> &tria)
{ 
  // const std::string filename = "../../meshes/coarse_WW.msh";
  const std::string filename = "../../meshes/structured_naca_2.msh";
  
  ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  pcout << "Reading from " << filename << std::endl;
  std::ifstream input_file(filename);
  GridIn<2>       grid_in;
  grid_in.attach_triangulation(tria);            //Attach this triangulation to be fed with the grid data
  grid_in.read_msh(input_file);                           //Read grid data from an msh file

  const types::manifold_id emitter = 1;                   //The type used to denote manifold indicators associated with every object of the mesh

  // FOR WIRE WIRE SIMULATION
  // double r_col = 1e-3;
  // double r_emi = 30e-5;
  // double dist_emi_col = 0.025;
  // const double X = -r_emi-dist_emi_col;

  // FOR NACA SIMULATION
  double X = -2.53;

  const Point<2> center(X,0.0);
  SphericalManifold<2> emitter_manifold(center);

  // FOR NACA SIMULATION
  const types::manifold_id collector = 2;
  CollectorGeometry<2> collector_manifold; 

  // FOR WIRE WIRE SIMULATION
  // const Point<2> center2(r_col,0.0);
  // SphericalManifold<2> collector_manifold(center2);               

  tria.set_all_manifold_ids_on_boundary(1, emitter);
  tria.set_manifold(emitter, emitter_manifold);
  tria.set_all_manifold_ids_on_boundary(2, collector);
  tria.set_manifold(collector, collector_manifold);
  pcout  << "Active cells: " << tria.n_active_cells() << std::endl;
}




// Main function
int main(int argc, char *argv[])
{


try
  {
    using namespace dealii;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);            //Initialize MPI (and, if deal.II was configured to use it, PETSc) and set the maximum number of threads used by deal.II to the given parameter.
    parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);
    create_triangulation(tria);
    InsIMEX<2> flow(tria);
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
return 0;
}
