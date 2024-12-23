


template <int dim>
InsIMEX<dim>::InsIMEX(parallel::distributed::Triangulation<dim> &tria)
    : viscosity(0.01),
    gamma(0.1),                             //Coefficient of the Grad-div stabilization
    degree(1),
    triangulation(tria),
    fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1),
    dof_handler(triangulation),            //Stores a list of degrees of freedom
    volume_quad_formula(degree + 2),
    face_quad_formula(degree + 2),
    mpi_communicator(MPI_COMM_WORLD),
    pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    time(1e0, 1e-3, 5e-3),
    timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
{
}


template <int dim>
void InsIMEX<dim>::setup_dofs()
{
    // Associate degrees of freedom (DoFs) with the mesh.
    dof_handler.distribute_dofs(fe);    
    DoFRenumbering::Cuthill_McKee(dof_handler);                      //Renumber the degrees of freedom according to the Cuthill-McKee method.

    std::vector<unsigned int> block_component(dim + 1, 0);           //Musk for the reording       //Block: in our case are two, one with dimension = dim and one with dimension = 1
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);    //We want to sort degree of freedom for a NS discretization so that we first get all velocities and then all the pressures so that the resulting matrix naturally decomposes into a 2×2 system.
    dofs_per_block = DoFTools::count_dofs_per_fe_block(dof_handler, block_component); //Returns a vector (of dim=2) of types::global_dof_index  
    
    // Partitioning.
    unsigned int dof_u = dofs_per_block[0];
    unsigned int dof_p = dofs_per_block[1];

    owned_partitioning.resize(2);
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, dof_u);      //Extract the set of locally owned DoF indices for each component within the mask that are owned by the current processor.
    owned_partitioning[1] = dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);

    DoFTools::extract_locally_relevant_dofs(dof_handler,locally_relevant_dofs);     //Extract the set of global DoF indices that are active on the current DoFHandler. This is the union of DoFHandler::locally_owned_dofs() and the DoF indices on all ghost cells.
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, dof_u);
    relevant_partitioning[1] = locally_relevant_dofs.get_view(dof_u, dof_u + dof_p);

    pcout << "   Number of active fluid cells: "
        << triangulation.n_global_active_cells() << std::endl
        << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " ("
        << dof_u << '+' << dof_p << ')' << std::endl;
}


template <int dim>
void InsIMEX<dim>::make_constraints_init()
{
    // Because the equation is written in incremental form, two constraints are needed: initial constraints and update constraints.
    
    initial_NS_constraints.clear();    //clear() the AffineConstraints object and supply an IndexSet with lines that may be constrained
    initial_NS_constraints.reinit(locally_relevant_dofs);      
    // DoFTools::make_hanging_node_constraints(dof_handler, nonzero_NS_constraints);      //Necessary when work with not homogeneous meshes (mesh non conformi) since in this cases we have the presence of hanging nodes. We have to impose constraints also on these.
    // DoFTools::make_hanging_node_constraints(dof_handler, zero_NS_constraints);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar vertical_velocity(1);
    const FEValuesExtractors::Vector vertical_velocity_and_pressure(1);    

    //Nonzero_NS_constraints //Assign B.C. different from zero in nonzero_constraints
    VectorTools::interpolate_boundary_values(dof_handler,                           
                                            0,                                      // Up and down 
                                            Functions::ZeroFunction<dim>(dim+1),    
                                            initial_NS_constraints,
                                            fe.component_mask(vertical_velocity));
    VectorTools::interpolate_boundary_values(dof_handler,
                                            3,                                      //Emitter                  
                                            Functions::ZeroFunction<dim>(dim+1),
                                            initial_NS_constraints,
                                            fe.component_mask(velocities));
    VectorTools::interpolate_boundary_values(dof_handler,
                                            4,                                      //Collector                                     
                                            Functions::ZeroFunction<dim>(dim+1),
                                            initial_NS_constraints,
                                            fe.component_mask(velocities));
    VectorTools::interpolate_boundary_values(dof_handler, 
                                            1,                                      // Inlet
                                            BoundaryValues<dim>(), 
                                            initial_NS_constraints, 
                                            fe.component_mask(velocities));
    

    VectorTools::interpolate_boundary_values(dof_handler,
                                            2,                                     // Outlet
                                            Functions::ZeroFunction<dim>(dim+1),
                                            initial_NS_constraints,
                                            fe.component_mask(vertical_velocity_and_pressure));  //Vertical velocity and pressure at outlet equal to 0

    initial_NS_constraints.close();         //After closing, no more entries are accepted

}


template <int dim>
void InsIMEX<dim>::make_constraints_update()   
{
    
    update_NS_constraints.clear();                 
    update_NS_constraints.reinit(locally_relevant_dofs);
    
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar vertical_velocity(1);
    const FEValuesExtractors::Vector vertical_velocity_and_pressure(1);  

    double inlet_value = 0.01; // We want that our velocity goes form 0.0 to 1.0 in 100 timesteps  

    //Zero_NS_constraints
    VectorTools::interpolate_boundary_values(dof_handler,
                                            0,
                                            Functions::ZeroFunction<dim>(dim+1),
                                            update_NS_constraints,
                                            fe.component_mask(vertical_velocity));
    VectorTools::interpolate_boundary_values(dof_handler,
                                            3,
                                            Functions::ZeroFunction<dim>(dim+1),
                                            update_NS_constraints,
                                            fe.component_mask(velocities));
    VectorTools::interpolate_boundary_values(dof_handler,
                                            4,
                                            Functions::ZeroFunction<dim>(dim+1),
                                            update_NS_constraints,
                                            fe.component_mask(velocities));
    VectorTools::interpolate_boundary_values(dof_handler,
                                            1,
                                            BoundaryValues<dim>(inlet_value),
                                            update_NS_constraints,
                                            fe.component_mask(velocities));
    VectorTools::interpolate_boundary_values(dof_handler,
                                            2, 
                                            Functions::ZeroFunction<dim>(dim+1),
                                            update_NS_constraints,
                                            fe.component_mask(vertical_velocity_and_pressure));

    update_NS_constraints.close();

}



template <int dim>
void InsIMEX<dim>::initialize_system()
{
    preconditioner.reset();
    system_matrix.clear();
    mass_matrix.clear();
    mass_schur.clear();

    BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);            //This class implements an array of compressed sparsity patterns (one for velocity and one for pressure in our case) that can be used to initialize objects of type BlockSparsityPattern.
    DoFTools::make_sparsity_pattern(dof_handler, dsp, initial_NS_constraints);     //Compute which entries of a matrix built on the given dof_handler may possibly be nonzero, and create a sparsity pattern (assigning it to dsp) object that represents these nonzero locations.
    sparsity_pattern.copy_from(dsp);                                            //Copy data from an object of type BlockDynamicSparsityPattern
    SparsityTools::distribute_sparsity_pattern(                    //Communicate rows in a dynamic sparsity pattern over MPI.
        dsp,                                                       //A dynamic sparsity pattern that has been built locally and for which we need to exchange entries with other processors to make sure that each processor knows all the elements of the rows of a matrix it stores and that may eventually be written to. This sparsity pattern will be changed as a result of this function: All entries in rows that belong to a different processor are sent to them and added there.
        dof_handler.locally_owned_dofs(),                          
        mpi_communicator,
        locally_relevant_dofs);                                    //The range of elements stored on the local MPI process. Only rows contained in this set are checked in dsp for transfer

    system_matrix.reinit(owned_partitioning, dsp, mpi_communicator);         //Efficiently reinit the block matrix for a parallel computation
    mass_matrix.reinit(owned_partitioning, dsp, mpi_communicator);

    // Only the $(1, 1)$ block in the mass schur matrix is used. (B M_u B.T)
    // Compute the sparsity pattern for mass schur in advance.
    // The only nonzero block has the same sparsity pattern as $BB^T$.
    BlockDynamicSparsityPattern schur_dsp(dofs_per_block, dofs_per_block);
    schur_dsp.block(1, 1).compute_mmult_pattern(sparsity_pattern.block(1, 0), sparsity_pattern.block(0, 1));
    mass_schur.reinit(owned_partitioning, schur_dsp, mpi_communicator);       //We want to reinitialize our matrix with new partitions of data

    // present_solution is ghosted because it is used in the output and mesh refinement functions.
    present_solution.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);      
    // solution_increment is non-ghosted (so doesn't contain info of "ghost zone") because the linear solver needs a completely distributed vector.
    solution_increment.reinit(owned_partitioning, mpi_communicator);   
    // system_rhs is non-ghosted because it is only used in the linear solver and residual evaluation.
    system_rhs.reinit(owned_partitioning, mpi_communicator);
}




// Assemble the system matrix, mass matrix, and the RHS.
// It can be used to assemble the entire system or only the RHS.
// An additional option is added to determine whether nonzero
// constraints or zero constraints should be used.
template <int dim>
void InsIMEX<dim>::assemble(bool use_initial_constraints,
                            bool assemble_system)
{
    TimerOutput::Scope timer_section(timer, "Assemble system");        //Enter the given section in the timer

    if (assemble_system)           //If assemble_system is false i assembly only rhs
    {
        system_matrix = 0;
        mass_matrix = 0;
    }
    system_rhs = 0;

    // Setup FEValues for cell integration with quadrature points
    FEValues<dim> fe_values(fe,                              //FEValues contains finite element evaluated in quadrature points of a cell.                            
                            volume_quad_formula,             //It implicitely uses a Q1 mapping
                            update_values | update_quadrature_points |
                            update_JxW_values | update_gradients);
    FEFaceValues<dim> fe_face_values(fe,
                                    face_quad_formula,
                                    update_values | update_normal_vectors |
                                    update_quadrature_points |
                                    update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = volume_quad_formula.size();

    const FEValuesExtractors::Vector velocities(0);   //Extractor calls velocities that takes the vector in position 0
    const FEValuesExtractors::Scalar pressure(dim);   //Extractor calls pressure that takes the scalar in position dim

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> local_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> current_velocity_values(n_q_points);     //vector of values of velocity in quadrature points              
    std::vector<Tensor<2, dim>> current_velocity_gradients(n_q_points);
    std::vector<double> current_velocity_divergences(n_q_points);
    std::vector<double> current_pressure_values(n_q_points);

    std::vector<double> div_phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
    std::vector<double> phi_p(dofs_per_cell);

    for (auto cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)     //Iterator from the first active cell to the last one
    {
        if (cell->is_locally_owned())
        {
            fe_values.reinit(cell);    //Reinitialize the gradients, Jacobi determinants, etc for the given cell of type "iterator into a Triangulation object", and the given finite element.

            if (assemble_system)
            {
                local_matrix = 0;
                local_mass_matrix = 0;
            }
            local_rhs = 0;

            fe_values[velocities].get_function_values(present_solution,
                                                    current_velocity_values);     //in current_velocity_values i stored values of present_solution in quadrature points

            fe_values[velocities].get_function_gradients(
                            present_solution, current_velocity_gradients);

            fe_values[velocities].get_function_divergences(
                            present_solution, current_velocity_divergences);

            fe_values[pressure].get_function_values(present_solution,
                                                    current_pressure_values);

            // Assemble the system matrix and mass matrix simultaneouly.
            // The mass matrix only uses the $(0, 0)$ and $(1, 1)$ blocks.
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    div_phi_u[k] = fe_values[velocities].divergence(k, q);           //Returns the value of the k-th shape function in q quadrature point
                    grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                    phi_u[k] = fe_values[velocities].value(k, q);
                    phi_p[k] = fe_values[pressure].value(k, q);
                }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                    if (assemble_system)
                    {
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                            local_matrix(i, j) +=
                            (viscosity *
                                scalar_product(grad_phi_u[j], grad_phi_u[i]) -
                            div_phi_u[i] * phi_p[j] -
                            phi_p[i] * div_phi_u[j] +
                            gamma * div_phi_u[j] * div_phi_u[i] +
                            phi_u[i] * phi_u[j] / time.get_delta_t() +

                            phi_u[i] * (current_velocity_gradients[q] * phi_u[j]) +
                            phi_u[i] * (grad_phi_u[j] * current_velocity_values[q]) 
                            ) *
                            fe_values.JxW(q);

                            local_mass_matrix(i, j) +=
                            (phi_u[i] * phi_u[j] + phi_p[i] * phi_p[j]) *
                            fe_values.JxW(q);
                        }
                    }
                    local_rhs(i) -=
                    (viscosity * scalar_product(current_velocity_gradients[q],
                                                grad_phi_u[i]) -
                    current_velocity_divergences[q] * phi_p[i] -
                    current_pressure_values[q] * div_phi_u[i] +
                    gamma * current_velocity_divergences[q] * div_phi_u[i] +
                    current_velocity_gradients[q] *
                        current_velocity_values[q] * phi_u[i]) *
                    fe_values.JxW(q);
                }
            }

            cell->get_dof_indices(local_dof_indices);

            const AffineConstraints<double> &constraints_used =
            use_initial_constraints ? initial_NS_constraints : update_NS_constraints;
            if (assemble_system)
            {
                constraints_used.distribute_local_to_global(local_matrix,             
                                                            local_rhs,
                                                            local_dof_indices,        
                                                            system_matrix,
                                                            system_rhs);
                constraints_used.distribute_local_to_global(
                local_mass_matrix, local_dof_indices, mass_matrix);
            }
            else
            {
                constraints_used.distribute_local_to_global(
                local_rhs, local_dof_indices, system_rhs);
            }
        }
    }

    if (assemble_system)
    {
        system_matrix.compress(VectorOperation::add);
        mass_matrix.compress(VectorOperation::add);
    }
    system_rhs.compress(VectorOperation::add);
}




// Solve the linear system using BICG solver with block preconditioner.
template <int dim>
std::pair<unsigned int, double>
InsIMEX<dim>::solve(bool use_initial_constraints, bool assemble_system, double time_step)
{
    if (assemble_system)
    {
        preconditioner.reset(new BlockSchurPreconditioner(timer,
                                                        gamma,
                                                        viscosity,
                                                        time.get_delta_t(),
                                                        owned_partitioning,
                                                        system_matrix,
                                                        mass_matrix,
                                                        mass_schur));
    }
    
    double coeff = 0.0 ;   // Adaptive tolerance coefficient to have a more precise solution in the initial steps
    if (time_step < 4) {
        coeff = 1e-4;
    } else if (time_step >= 4 && time_step < 10) {
        coeff = 1e-3;
    } else {
        coeff = 1e-2;
    }
    
    SolverControl solver_control(                                         //Used by iterative methods to determine whether the iteration should be continued
    system_matrix.m(), coeff * system_rhs.l2_norm(), true);


    SolverBicgstab<PETScWrappers::MPI::BlockVector> bicg(solver_control);

    pcout << "system_matrix frob norm is " << system_matrix.frobenius_norm() << std::endl;
    pcout << "system_rhs l2 norm is " << system_rhs.l2_norm() << std::endl;
    // The solution vector must be non-ghosted
    bicg.solve(system_matrix, solution_increment, system_rhs, *preconditioner);

    const AffineConstraints<double> &constraints_used =
    use_initial_constraints ? initial_NS_constraints : update_NS_constraints;
    constraints_used.distribute(solution_increment);

    return {solver_control.last_step(), solver_control.last_value()};
}




template <int dim>
void InsIMEX<dim>::run()
{
    pcout << "Running with PETSc on "
        << Utilities::MPI::n_mpi_processes(mpi_communicator)         //Return the number of MPI processes there exist in the given communicator object
        << " MPI rank(s)..." << std::endl;

    triangulation.refine_global(0);
    setup_dofs();
    make_constraints_init();
    initialize_system();

    //Initialize the present_solution
    MappingQ1<dim>  mapping; 
    PETScWrappers::MPI::BlockVector tmp1;
	tmp1.reinit(owned_partitioning, mpi_communicator);
	tmp1 = present_solution;
	VectorTools::interpolate(mapping, dof_handler, Functions::ZeroFunction<dim>(dim+1), tmp1);   //Zero initial solution
	present_solution = tmp1;

    // Time loop.
    while (time.end() - time.current() > 1e-12)            //We are not still at the end
    {
        if (time.get_timestep() == 0)
        {
            output_results(0);                             //Print result at 0 timestep
            pcout << "L2 norm of the initial condition: " << present_solution.l2_norm() << std::endl;
        }
        time.increment();
        std::cout.precision(6);                           //This means that when floating-point numbers are output to std::cout, they will be displayed with up to 6 digits after the decimal point.
        std::cout.width(12);                              //Sets the minimum field width to 12 characters
        pcout << std::string(96, '*') << std::endl
            << "Time step = " << time.get_timestep()
            << ", at t = " << std::scientific << time.current() << std::endl;

        // Resetting
        solution_increment = 0;

        bool apply_initial_constraints = false;
        
        bool assemble_system = true;

        make_constraints_update();

        assemble(apply_initial_constraints, assemble_system);

        auto state = solve(apply_initial_constraints, assemble_system, time.get_timestep());
        // Note we have to use a non-ghosted vector to do the addition.
        PETScWrappers::MPI::BlockVector tmp;
        tmp.reinit(owned_partitioning, mpi_communicator);   //We do this since present solution is a ghost vector so is read-only
        tmp = present_solution;
        tmp += solution_increment;
        present_solution = tmp;
        pcout << std::scientific << std::left << " GMRES_ITR = " << std::setw(3)
            << state.first << " GMRES_RES = " << state.second << std::endl;

        pcout << "L2 norm of the present solution: " << present_solution.l2_norm() << std::endl;

        // Output
        if (time.time_to_output())
        {
            output_results(time.get_timestep());         //Print result at current timestep
        }
    }
}




template <int dim>
void InsIMEX<dim>::output_results(const unsigned int output_index) const
{

    // Base directory for output
    std::string base_directory = "../output";

    // Directory to store the results of this simulation
    std::string output_directory = base_directory + "/NS_Simulation/";

    // Ensure the output directory is created (if it doesn't exist)
    if (!std::filesystem::exists(output_directory))
    {
        std::filesystem::create_directory(output_directory);
    }

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);                                              
    // vector to be output must be ghosted
    data_out.add_data_vector(present_solution,
                            solution_names,
                            DataOut<dim>::type_dof_data,
                            data_component_interpretation);


    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = triangulation.locally_owned_subdomain();

    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches();

    data_out.write_vtu_with_pvtu_record(output_directory, "solution", output_index, mpi_communicator, 2, 1);


}

