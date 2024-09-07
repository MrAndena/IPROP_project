// Constructor of the class

template <int dim>
CompleateProblem<dim>::CompleateProblem(parallel::distributed::Triangulation<dim> &tria, 
                                        const data_struct &d,
                                        unsigned short int i)
  : m_data(d)
  , simulation_tag(i)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , triangulation(tria)
  , fe(1) 
  , dof_handler(tria)
  , mapping()
  , timestep(1.e-5)
  ,	viscosity(d.navier_stokes.physical_parameters.viscosity)
  , gamma(d.navier_stokes.physical_parameters.gamma) 
  , degree(d.navier_stokes.numerical_parameters.FE_degree) 
  , NS_fe(FE_Q<dim>(degree+1), dim, FE_Q<dim>(degree), 1)
  , NS_dof_handler(tria)
  , volume_quad_formula(degree + 2)
  , face_quad_formula(degree + 2)
  , NS_mapping() // serve ?
{}

//---------------------------------------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void CompleateProblem<dim>::setup_poisson()
  {

    dof_handler.distribute_dofs(fe);

    // INDEX SETS INITIALIZATION
    locally_owned_dofs = dof_handler.locally_owned_dofs();                           //local dofs
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);    //local dofs + ghost dofs
    
    // PETSC VECTORS DECLARATIONS 
    potential.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);   //ghosted
    poisson_newton_update.reinit(locally_owned_dofs, mpi_communicator);              //non-ghosted
    poisson_rhs.reinit(locally_owned_dofs, mpi_communicator);                        //non-ghosted

    Field_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);//  ghosted o non ghosted?
    Field_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    // ZERO_CONSTRAINTS FOR NEWTON POISSON PROBLEM
    zero_constraints_poisson.clear();
    zero_constraints_poisson.reinit(locally_relevant_dofs);
    zero_constraints_poisson.close(); 


    // DYNAMIC SPARSITY PATTERN AND POISSON MATRICES 
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints, false); 
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);

  
    system_matrix_poisson.clear(); //store the matrix used to solve nlpoisson
    system_matrix_poisson.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    laplace_matrix_poisson.clear(); //store laplace matrix
    laplace_matrix_poisson.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    mass_matrix_poisson.clear();  //store mass matrix
    mass_matrix_poisson.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    density_matrix.clear(); // store the term: M(n+p)q0/V_TH
    density_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void CompleateProblem<dim>::assemble_laplace_matrix()
  {
    const QTrapezoid<dim> quadrature_formula;

    laplace_matrix_poisson = 0;

    FEValues<dim> fe_values(fe,
			    quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  
    unsigned int q_point = 0, idof = 0, jdof = 0;
    for (const auto &cell : dof_handler.active_cell_iterators()){
      if (cell->is_locally_owned()){
	      cell_matrix = 0.;
        fe_values.reinit(cell);
	
        for (q_point = 0; q_point < n_q_points; ++q_point) {
	  for (idof = 0; idof < dofs_per_cell; ++idof) {
	    for (jdof = 0; jdof < dofs_per_cell; ++jdof)
	      cell_matrix(idof, jdof) += fe_values.shape_grad(idof, q_point) *
		fe_values.shape_grad(jdof, q_point) * fe_values.JxW(q_point);
	  }
	}

        cell->get_dof_indices(local_dof_indices);
        zero_constraints_poisson.distribute_local_to_global(cell_matrix,       
                                               local_dof_indices,
                                               laplace_matrix_poisson);
      }
    }

    laplace_matrix_poisson.compress(VectorOperation::add);

    // pcout << " The L_INF norm of the laplace matrix is "<<laplace_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the laplace matrix is "<<laplace_matrix.frobenius_norm() <<std::endl<<std::endl;
    // pcout << "   End of Assembling Laplce matrix "<< std::endl<<std::endl;
  }

//------------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void CompleateProblem<dim>::assemble_mass_matrix()
  {
    const QTrapezoid<dim> quadrature_formula;

    mass_matrix_poisson = 0;

    FEValues<dim> fe_values(fe,
			    quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    unsigned int q_point = 0, idof = 0, jdof = 0;
    for (const auto &cell : dof_handler.active_cell_iterators()){
      if (cell->is_locally_owned())
	{
	  cell_matrix = 0.;

	  fe_values.reinit(cell);

	  for (q_point = 0; q_point < n_q_points; ++q_point) {
	    for (idof = 0; idof < dofs_per_cell; ++idof) {
	      for (jdof = 0; jdof < dofs_per_cell; ++jdof)
		cell_matrix(idof, jdof) += fe_values.shape_value(idof, q_point) *
		  fe_values.shape_value(jdof, q_point) * fe_values.JxW(q_point);
	    }
	  }

	  cell->get_dof_indices(local_dof_indices);
	  zero_constraints_poisson.distribute_local_to_global(cell_matrix,   
						                               local_dof_indices,
						                               mass_matrix_poisson );
	}
    }

    mass_matrix_poisson.compress(VectorOperation::add);

    // pcout << " The L_INF norm of the mass matrix is "<<mass_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the mass matrix is "<<mass_matrix.frobenius_norm() <<std::endl<<std::endl;
    // pcout << "   End of Assembling Mass matrix "<< std::endl<<std::endl;

  }
  
//----------------------------------------------------------------------------------------------------------------------------  
  








  //HELPER FUNCTION IMPLEMENTATION

  void bernoulli (double x, double &bp, double &bn)
  {
    const double xlim = 1.0e-2;
    double ax  = fabs(x);

    bp  = 0.0;
    bn  = 0.0;

    //  X=0
    if (x == 0.0)
      {
        bp = 1.0;
        bn = 1.0;
        return;
      }

    // ASYMPTOTICS
    if (ax > 80.0)
      {
        if (x > 0.0)
          {
            bp = 0.0;
            bn = x;
          }
        else
          {
            bp = -x;
            bn = 0.0;
          }
        return;
      }

    // INTERMEDIATE VALUES
    if (ax <= 80 &&  ax > xlim)
      {
        bp = x / (exp (x) - 1.0);
        bn = x + bp;
        return;
      }

    // SMALL VALUES
    if (ax <= xlim &&  ax != 0.0)
      {
        double jj = 1.0;
        double fp = 1.0;
        double fn = 1.0;
        double df = 1.0;
        double segno = 1.0;
        while (fabs (df) > 1.0e-16)
          {
            jj += 1.0;
            segno = -segno;
            df = df * x / jj;
            fp = fp + df;
            fn = fn + segno * df;
          }
        bp = 1 / fp;
        bn = 1 / fn;
        return;
      }

  };

//-------------------------------------------------------------------------------------------------------------------------------------------------

double side_length (const Point<2> a, const Point<2> b)
{
	double length = 0.;

	if (a[0] == b[0])
		length = std::abs(a[1] - b[1]);
	else if (a[1] == b[1])
		length = std::abs(a[0] - b[0]);
	else
		length = std::sqrt(a[0]*a[0] + b[0]*b[0] - 2.*a[0]*b[0] + a[1]*a[1] + b[1]*b[1] - 2.*a[1]*b[1]);

	return length;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------

double triangle_denom(const Point<2> a, const Point<2> b, const Point<2> c)
{
	const double x1 = a[0];
	const double y1 = a[1];

	const double x2 = b[0];
	const double y2 = b[1];

	const double x3 = c[0];
	const double y3 = c[1];

	const double denom = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);

	return denom;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------

Tensor<1,2> face_normal(const Point<2> a, const Point<2> b) {

	Tensor<1,2> tangent, normal;

	tangent[0] = b[0] - a[0];
	tangent[1] = b[1] - a[1];

	normal[0] = -tangent[1];
	normal[1] = tangent[0];

	return normal;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------

FullMatrix<double> compute_triangle_matrix(const Point<2> a, const Point<2> b, const Point<2> c, const double alpha12, const double alpha23, const double alpha31, const double D)
{
	const unsigned int size = 3;
	FullMatrix<double> tria_matrix(size,size);

	tria_matrix = 0;

	const double denom = triangle_denom(a,b,c);
	const double area = 0.5*std::abs(denom);

	const Tensor<1,2> grad_psi_1 = face_normal(b,c)/denom;
	const Tensor<1,2> grad_psi_2 = face_normal(c,a)/denom;
	const Tensor<1,2> grad_psi_3 = face_normal(a,b)/denom;

	const double l_12 = grad_psi_1 * grad_psi_2;
	const double l_23 = grad_psi_2 * grad_psi_3;
	const double l_31 = grad_psi_1 * grad_psi_3;

	double bp12, bn12, bp23, bn23, bp31, bn31;

	bernoulli(alpha12,bp12,bn12);
	bernoulli(alpha23,bp23,bn23);
	bernoulli(alpha31,bp31,bn31);

	tria_matrix(0,1) = D * area * bp12 * l_12;
	tria_matrix(0,2) = D * area * bn31 * l_31;

	tria_matrix(1,0) = D * area * bn12 * l_12;
	tria_matrix(1,2) = D * area * bp23 * l_23;

	tria_matrix(2,0) = D * area * bp31 * l_31;
	tria_matrix(2,1) = D * area * bn23 * l_23;
	
	tria_matrix(0,0) = - (tria_matrix(1,0)+tria_matrix(2,0));
	tria_matrix(1,1) = - (tria_matrix(0,1)+tria_matrix(2,1));
	tria_matrix(2,2) = - (tria_matrix(0,2)+tria_matrix(1,2));

	return tria_matrix;
}