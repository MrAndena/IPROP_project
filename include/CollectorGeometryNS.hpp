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

#include <deal.II/grid/grid_in.h>                //Aggiunto da me
#include <deal.II/base/geometry_info.h>          //Aggiunto da me
#include <deal.II/grid/manifold_lib.h>           // To use manifolds
#include <deal.II/base/mpi.h>                    //Aggiunto da me 

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
#include "BlockSchurPreconditioner.hpp"
#include "Time.hpp"
#include "BoundaryValues.hpp"
#include "data_struct.hpp"


using namespace dealii;


// Collector Manifold - START

double get_collector_height(const double &p, const data_struct &s_data)
{
  
    
  double collector_length = s_data.geometrical_parameters.chord_length;
  const double x = p/collector_length;
	double y = 0;

	if ( abs(x-1.) > 1e-12 && abs(x) > 1e-12 ) {
		double a0 = 0.2969;
		double a1 = -0.126;
		double a2 = -0.3516;
		double a3 = 0.2843;
		double a4 = -0.1036; // or -0.1015 for an open trailing edge
		double t = s_data.geometrical_parameters.NACA_digits; // Last 2 digits of the NACA by 100

		y = 5*t*( a0 * sqrt(x) + a1 * x + a2 * pow(x,2.0) + a3 * pow(x,3.0) + a4 * pow(x,4.0) );
	}

	return y * collector_length;
}



template <int dim>
class CollectorGeometry : public ChartManifold<dim, dim, dim-1>     //ChartManifold is a class describes mappings that can be expressed in terms of charts.
  {
public:
  virtual Point<dim-1> pull_back(const Point<dim> &space_point) const override;        //Pull back the given point in spacedim to the Euclidean chartdim dimensional space

  Point<dim> push_forward(const Point<dim-1> &chart_point) const override;     //serve per non rendere la classe virtuale
  Point<dim> push_forward(const Point<dim-1> &chart_point, const data_struct& s_data) const; 
  
  virtual std::unique_ptr<Manifold<dim, dim>> clone() const override;                  //Return a copy of this manifold

  };

template <int dim>
std::unique_ptr<Manifold<dim, dim>> CollectorGeometry<dim>::clone() const
  {
return std::make_unique<CollectorGeometry<dim>>();
  }

// Serve se no CollectorGeometry è una virtual class
template <int dim>
Point<dim> CollectorGeometry<dim>::push_forward(const Point<dim-1>  &x) const          //Input: a chart point that in our case is a 1D point 
{
  //const double y = get_collector_height(x[0], s_data);

  Point<dim> p;
  p[0] = x[0];  //NON SERVE A NIENTE
  //p[0] = x[0]; p[1] = y;

  // if (dim == 3) {
  // p[2] = x[1];
  // }

  return p;                                                                              //Output: a point of our collector in 2D 
}

template <int dim>
Point<dim> CollectorGeometry<dim>::push_forward(const Point<dim-1>  &x, const data_struct& s_data) const          //Input: a chart point that in our case is a 1D point 
{
  const double y = get_collector_height(x[0], s_data);

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

  

void create_triangulation(parallel::distributed::Triangulation<2> &tria, const data_struct& s_data, unsigned short int i)
{ 
  if(s_data.simulation_specification.ID_simulation <=2){ //NACA NS 

        std::string name_mesh = "naca_" + std::to_string(i) + ".msh";
        //std::string name_mesh = "WireWire_" + std::to_string(i) + ".msh";
      
        std::string filename = "../output/meshes/"+name_mesh;
        std::cout <<"Reading the mesh from " << filename << std::endl;
        std::ifstream input_file(filename);
        GridIn<2>       grid_in;
        grid_in.attach_triangulation(tria);           
        grid_in.read_msh(input_file);                 

        const types::manifold_id emitter = 3; 
        const types::manifold_id collector = 4;       
       
        const double distance_emitter_collector = s_data.geometrical_parameters.distance_emitter_collector;
        const double r_emi = s_data.geometrical_parameters.emitter_radius[i];
        double X = -distance_emitter_collector - r_emi; 
        const Point<2> center(X,0.0);

        SphericalManifold<2> emitter_manifold(center);

        

        CollectorGeometry<2> collector_manifold; // QUA USI LA CLASSE PER IL PROFILO NACA


        // for wire wire simulation
        //const double r_col = 1e-3;
        //const Point<2> center2(r_col,0.0);

        //SphericalManifold<2> collector_manifold(center2);               

        tria.set_all_manifold_ids_on_boundary(3, emitter);
        tria.set_manifold(emitter, emitter_manifold);
        tria.set_all_manifold_ids_on_boundary(4, collector);
        tria.set_manifold(collector, collector_manifold);
        
        std::cout  << "Active cells: " << tria.n_active_cells() << std::endl;

  }else{ //NL Poisson

        const std::string filename = "../output/meshes/Structured_Square.msh"; //remember: you run from build directory

        ConditionalOStream pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

        std::ifstream input_file(filename);

        if (!input_file.is_open())
        {
          std::cerr << "Errore: impossibile aprire il file di mesh." << std::endl;
          return;
        }

        pcout <<"Reading the mesh from " << filename << std::endl;
        GridIn<2>  grid_in; //This class implements an input mechanism for grid data. It allows to read a grid structure into a triangulation object
        grid_in.attach_triangulation(tria); //we pass to grid_in our (empty) triangulation
        grid_in.read_msh(input_file); // read the msh file
        pcout << " Grid read correctly "<< std::endl;

  }
}