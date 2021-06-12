

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2000 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE.md at 
 * the top level directory of deal.II. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004, 2005, 
 * Timo Heister, 2013 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/multithread_info.h> 
#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/petsc_vector.h> 
#include <deal.II/lac/petsc_sparse_matrix.h> 
#include <deal.II/lac/petsc_solver.h> 
#include <deal.II/lac/petsc_precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/sparsity_tools.h> 
#include <deal.II/distributed/shared_tria.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 


#include <deal.II/base/symmetric_tensor.h> 


#include <deal.II/physics/transformations.h> 


#include <fstream> 
#include <iostream> 
#include <iomanip> 


namespace Step18 
{ 
  using namespace dealii; 



  template <int dim> 
  struct PointHistory 
  { 
    SymmetricTensor<2, dim> old_stress; 
  }; 


  template <int dim> 
  SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda, 
                                                   const double mu) 
  { 
    SymmetricTensor<4, dim> tmp; 
    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = 0; j < dim; ++j) 
        for (unsigned int k = 0; k < dim; ++k) 
          for (unsigned int l = 0; l < dim; ++l) 
            tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) + 
                               ((i == l) && (j == k) ? mu : 0.0) + 
                               ((i == j) && (k == l) ? lambda : 0.0)); 
    return tmp; 
  } 









  template <int dim> 
  inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values, 
                                            const unsigned int   shape_func, 
                                            const unsigned int   q_point) 
  { 


    SymmetricTensor<2, dim> tmp; 


    for (unsigned int i = 0; i < dim; ++i) 
      tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i]; 


    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = i + 1; j < dim; ++j) 
        tmp[i][j] = 
          (fe_values.shape_grad_component(shape_func, q_point, i)[j] + 
           fe_values.shape_grad_component(shape_func, q_point, j)[i]) / 
          2; 

    return tmp; 
  } 



  template <int dim> 
  inline SymmetricTensor<2, dim> 
  get_strain(const std::vector<Tensor<1, dim>> &grad) 
  { 
    Assert(grad.size() == dim, ExcInternalError()); 

    SymmetricTensor<2, dim> strain; 
    for (unsigned int i = 0; i < dim; ++i) 
      strain[i][i] = grad[i][i]; 

    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = i + 1; j < dim; ++j) 
        strain[i][j] = (grad[i][j] + grad[j][i]) / 2; 

    return strain; 
  } 




  Tensor<2, 2> get_rotation_matrix(const std::vector<Tensor<1, 2>> &grad_u) 
  { 


    const double curl = (grad_u[1][0] - grad_u[0][1]); 


    const double angle = std::atan(curl); 


    return Physics::Transformations::Rotations::rotation_matrix_2d(-angle); 
  } 


  Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>> &grad_u) 
  { 


    const Point<3> curl(grad_u[2][1] - grad_u[1][2], 
                        grad_u[0][2] - grad_u[2][0], 
                        grad_u[1][0] - grad_u[0][1]); 


    const double tan_angle = std::sqrt(curl * curl); 
    const double angle     = std::atan(tan_angle); 



    if (std::abs(angle) < 1e-9) 
      { 
        static const double rotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; 
        static const Tensor<2, 3> rot(rotation); 
        return rot; 
      } 


    const Point<3> axis = curl / tan_angle; 
    return Physics::Transformations::Rotations::rotation_matrix_3d(axis, 
                                                                   -angle); 
  } 





  template <int dim> 
  class TopLevel 
  { 
  public: 
    TopLevel(); 
    ~TopLevel(); 
    void run(); 

  private: 


    void create_coarse_grid(); 

    void setup_system(); 

    void assemble_system(); 

    void solve_timestep(); 

    unsigned int solve_linear_problem(); 

    void output_results() const; 


    void do_initial_timestep(); 

    void do_timestep(); 


    void refine_initial_grid(); 


    void move_mesh(); 


    void setup_quadrature_point_history(); 


    void update_quadrature_point_history(); 


    parallel::shared::Triangulation<dim> triangulation; 

    FESystem<dim> fe; 

    DoFHandler<dim> dof_handler; 

    AffineConstraints<double> hanging_node_constraints; 


    const QGauss<dim> quadrature_formula; 


    std::vector<PointHistory<dim>> quadrature_point_history; 




    PETScWrappers::MPI::SparseMatrix system_matrix; 

    PETScWrappers::MPI::Vector system_rhs; 

    Vector<double> incremental_displacement; 


    double       present_time; 
    double       present_timestep; 
    double       end_time; 
    unsigned int timestep_no; 


    MPI_Comm mpi_communicator; 

    const unsigned int n_mpi_processes; 

    const unsigned int this_mpi_process; 

    ConditionalOStream pcout; 


    IndexSet locally_owned_dofs; 
    IndexSet locally_relevant_dofs; 


    static const SymmetricTensor<4, dim> stress_strain_tensor; 
  }; 



  template <int dim> 
  class BodyForce : public Function<dim> 
  { 
  public: 
    BodyForce(); 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  values) const override; 

    virtual void 
    vector_value_list(const std::vector<Point<dim>> &points, 
                      std::vector<Vector<double>> &  value_list) const override; 
  }; 

  template <int dim> 
  BodyForce<dim>::BodyForce() 
    : Function<dim>(dim) 
  {} 

  template <int dim> 
  inline void BodyForce<dim>::vector_value(const Point<dim> & /*p*/, 
                                           Vector<double> &values) const 
  { 
    Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim)); 

    const double g   = 9.81; 
    const double rho = 7700; 

    values          = 0; 
    values(dim - 1) = -rho * g; 
  } 

  template <int dim> 
  void BodyForce<dim>::vector_value_list( 
    const std::vector<Point<dim>> &points, 
    std::vector<Vector<double>> &  value_list) const 
  { 
    const unsigned int n_points = points.size(); 

    Assert(value_list.size() == n_points, 
           ExcDimensionMismatch(value_list.size(), n_points)); 

    for (unsigned int p = 0; p < n_points; ++p) 
      BodyForce<dim>::vector_value(points[p], value_list[p]); 
  } 






  template <int dim> 
  class IncrementalBoundaryValues : public Function<dim> 
  { 
  public: 
    IncrementalBoundaryValues(const double present_time, 
                              const double present_timestep); 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  values) const override; 

    virtual void 
    vector_value_list(const std::vector<Point<dim>> &points, 
                      std::vector<Vector<double>> &  value_list) const override; 

  private: 
    const double velocity; 
    const double present_time; 
    const double present_timestep; 
  }; 

  template <int dim> 
  IncrementalBoundaryValues<dim>::IncrementalBoundaryValues( 
    const double present_time, 
    const double present_timestep) 
    : Function<dim>(dim) 
    , velocity(.08) 
    , present_time(present_time) 
    , present_timestep(present_timestep) 
  {} 

  template <int dim> 
  void 
  IncrementalBoundaryValues<dim>::vector_value(const Point<dim> & /*p*/, 
                                               Vector<double> &values) const 
  { 
    Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim)); 

    values    = 0; 
    values(2) = -present_timestep * velocity; 
  } 

  template <int dim> 
  void IncrementalBoundaryValues<dim>::vector_value_list( 
    const std::vector<Point<dim>> &points, 
    std::vector<Vector<double>> &  value_list) const 
  { 
    const unsigned int n_points = points.size(); 

    Assert(value_list.size() == n_points, 
           ExcDimensionMismatch(value_list.size(), n_points)); 

    for (unsigned int p = 0; p < n_points; ++p) 
      IncrementalBoundaryValues<dim>::vector_value(points[p], value_list[p]); 
  } 



  template <int dim> 
  const SymmetricTensor<4, dim> TopLevel<dim>::stress_strain_tensor = 
    get_stress_strain_tensor<dim>(
      /*lambda = */ 9.695e10, 
      /*mu =  */  7.617e10)



  template <int dim> 
  TopLevel<dim>::TopLevel() 
    : triangulation(MPI_COMM_WORLD) 
    , fe(FE_Q<dim>(1), dim) 
    , dof_handler(triangulation) 
    , quadrature_formula(fe.degree + 1) 
    , present_time(0.0) 
    , present_timestep(1.0) 
    , end_time(10.0) 
    , timestep_no(0) 
    , mpi_communicator(MPI_COMM_WORLD) 
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)) 
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)) 
    , pcout(std::cout, this_mpi_process == 0) 
  {} 

  template <int dim> 
  TopLevel<dim>::~TopLevel() 
  { 
    dof_handler.clear(); 
  } 


  template <int dim> 
  void TopLevel<dim>::run() 
  { 
    do_initial_timestep(); 

    while (present_time < end_time) 
      do_timestep(); 
  } 



  template <int dim> 
  void TopLevel<dim>::create_coarse_grid() 
  { 
    const double inner_radius = 0.8, outer_radius = 1; 
    GridGenerator::cylinder_shell(triangulation, 3, inner_radius, outer_radius); 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary()) 
          { 
            const Point<dim> face_center = face->center(); 

            if (face_center[2] == 0) 
              face->set_boundary_id(0); 
            else if (face_center[2] == 3) 
              face->set_boundary_id(1); 
            else if (std::sqrt(face_center[0] * face_center[0] + 
                               face_center[1] * face_center[1]) < 
                     (inner_radius + outer_radius) / 2) 
              face->set_boundary_id(2); 
            else 
              face->set_boundary_id(3); 
          } 


    triangulation.refine_global(1); 


    setup_quadrature_point_history(); 
  } 



  template <int dim> 
  void TopLevel<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    locally_owned_dofs = dof_handler.locally_owned_dofs(); 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 


    hanging_node_constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, 
                                            hanging_node_constraints); 
    hanging_node_constraints.close(); 



    DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    sparsity_pattern, 
                                    hanging_node_constraints, 
                                    /*保持约束性dofs  */ false)

    SparsityTools::distribute_sparsity_pattern(sparsity_pattern, 
                                               locally_owned_dofs, 
                                               mpi_communicator, 
                                               locally_relevant_dofs); 





    system_matrix.reinit(locally_owned_dofs, 
                         locally_owned_dofs, 
                         sparsity_pattern, 
                         mpi_communicator); 



    system_rhs.reinit(locally_owned_dofs, mpi_communicator); 
    incremental_displacement.reinit(dof_handler.n_dofs()); 
  } 




  template <int dim> 
  void TopLevel<dim>::assemble_system() 
  { 
    system_rhs    = 0; 
    system_matrix = 0; 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    BodyForce<dim>              body_force; 
    std::vector<Vector<double>> body_force_values(n_q_points, 
                                                  Vector<double>(dim)); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell_matrix = 0; 
          cell_rhs    = 0; 

          fe_values.reinit(cell); 


          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
                { 
                  const SymmetricTensor<2, dim> 
                    eps_phi_i = get_strain(fe_values, i, q_point), 
                    eps_phi_j = get_strain(fe_values, j, q_point); 

                  cell_matrix(i, j) += (eps_phi_i *            // 
                                        stress_strain_tensor * // 
                                        eps_phi_j              // 
                                        ) *                    // 
                                       fe_values.JxW(q_point); // 
                } 


          const PointHistory<dim> *local_quadrature_points_data = 
            reinterpret_cast<PointHistory<dim> *>(cell->user_pointer()); 


          body_force.vector_value_list(fe_values.get_quadrature_points(), 
                                       body_force_values); 


          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const unsigned int component_i = 
                fe.system_to_component_index(i).first; 

              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
                { 
                  const SymmetricTensor<2, dim> &old_stress = 
                    local_quadrature_points_data[q_point].old_stress; 

                  cell_rhs(i) += 
                    (body_force_values[q_point](component_i) * 
                       fe_values.shape_value(i, q_point) - 
                     old_stress * get_strain(fe_values, i, q_point)) * 
                    fe_values.JxW(q_point); 
                } 
            } 


          cell->get_dof_indices(local_dof_indices); 

          hanging_node_constraints.distribute_local_to_global(cell_matrix, 
                                                              cell_rhs, 
                                                              local_dof_indices, 
                                                              system_matrix, 
                                                              system_rhs); 
        } 


    system_matrix.compress(VectorOperation::add); 
    system_rhs.compress(VectorOperation::add); 




    FEValuesExtractors::Scalar                z_component(dim - 1); 
    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(dim), 
                                             boundary_values); 
    VectorTools::interpolate_boundary_values( 
      dof_handler, 
      1, 
      IncrementalBoundaryValues<dim>(present_time, present_timestep), 
      boundary_values, 
      fe.component_mask(z_component)); 

    PETScWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator); 
    MatrixTools::apply_boundary_values( 
      boundary_values, system_matrix, tmp, system_rhs, false); 
    incremental_displacement = tmp; 
  } 



  template <int dim> 
  void TopLevel<dim>::solve_timestep() 
  { 
    pcout << "    Assembling system..." << std::flush; 
    assemble_system(); 
    pcout << " norm of rhs is " << system_rhs.l2_norm() << std::endl; 

    const unsigned int n_iterations = solve_linear_problem(); 

    pcout << "    Solver converged in " << n_iterations << " iterations." 
          << std::endl; 

    pcout << "    Updating quadrature point data..." << std::flush; 
    update_quadrature_point_history(); 
    pcout << std::endl; 
  } 



  template <int dim> 
  unsigned int TopLevel<dim>::solve_linear_problem() 
  { 
    PETScWrappers::MPI::Vector distributed_incremental_displacement( 
      locally_owned_dofs, mpi_communicator); 
    distributed_incremental_displacement = incremental_displacement; 

    SolverControl solver_control(dof_handler.n_dofs(), 
                                 1e-16 * system_rhs.l2_norm()); 

    PETScWrappers::SolverCG cg(solver_control, mpi_communicator); 

    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix); 

    cg.solve(system_matrix, 
             distributed_incremental_displacement, 
             system_rhs, 
             preconditioner); 

    incremental_displacement = distributed_incremental_displacement; 

    hanging_node_constraints.distribute(incremental_displacement); 

    return solver_control.last_step(); 
  } 




  template <int dim> 
  void TopLevel<dim>::output_results() const 
  { 
    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 

    std::vector<std::string> solution_names; 
    switch (dim) 
      { 
        case 1: 
          solution_names.emplace_back("delta_x"); 
          break; 
        case 2: 
          solution_names.emplace_back("delta_x"); 
          solution_names.emplace_back("delta_y"); 
          break; 
        case 3: 
          solution_names.emplace_back("delta_x"); 
          solution_names.emplace_back("delta_y"); 
          solution_names.emplace_back("delta_z"); 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    data_out.add_data_vector(incremental_displacement, solution_names); 


    Vector<double> norm_of_stress(triangulation.n_active_cells()); 
    { 


      for (auto &cell : triangulation.active_cell_iterators()) 
        if (cell->is_locally_owned()) 
          { 


            SymmetricTensor<2, dim> accumulated_stress; 
            for (unsigned int q = 0; q < quadrature_formula.size(); ++q) 
              accumulated_stress += 
                reinterpret_cast<PointHistory<dim> *>(cell->user_pointer())[q] 
                  .old_stress; 


            norm_of_stress(cell->active_cell_index()) = 
              (accumulated_stress / quadrature_formula.size()).norm(); 
          } 


        else 
          norm_of_stress(cell->active_cell_index()) = -1e+20; 
    } 


    data_out.add_data_vector(norm_of_stress, "norm_of_stress"); 


    std::vector<types::subdomain_id> partition_int( 
      triangulation.n_active_cells()); 
    GridTools::get_subdomain_association(triangulation, partition_int); 
    const Vector<double> partitioning(partition_int.begin(), 
                                      partition_int.end()); 
    data_out.add_data_vector(partitioning, "partitioning"); 


    data_out.build_patches(); 



    const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
      "./", "solution", timestep_no, mpi_communicator, 4); 


    if (this_mpi_process == 0) 
      { 


        static std::vector<std::pair<double, std::string>> times_and_names; 
        times_and_names.push_back( 
          std::pair<double, std::string>(present_time, pvtu_filename)); 
        std::ofstream pvd_output("solution.pvd"); 
        DataOutBase::write_pvd_record(pvd_output, times_and_names); 
      } 
  } 




  template <int dim> 
  void TopLevel<dim>::do_initial_timestep() 
  { 
    present_time += present_timestep; 
    ++timestep_no; 
    pcout << "Timestep " << timestep_no << " at time " << present_time 
          << std::endl; 

    for (unsigned int cycle = 0; cycle < 2; ++cycle) 
      { 
        pcout << "  Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          create_coarse_grid(); 
        else 
          refine_initial_grid(); 

        pcout << "    Number of active cells:       " 
              << triangulation.n_active_cells() << " (by partition:"; 
        for (unsigned int p = 0; p < n_mpi_processes; ++p) 
          pcout << (p == 0 ? ' ' : '+') 
                << (GridTools::count_cells_with_subdomain_association( 
                     triangulation, p)); 
        pcout << ")" << std::endl; 

        setup_system(); 

        pcout << "    Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (by partition:"; 
        for (unsigned int p = 0; p < n_mpi_processes; ++p) 
          pcout << (p == 0 ? ' ' : '+') 
                << (DoFTools::count_dofs_with_subdomain_association(dof_handler, 
                                                                    p)); 
        pcout << ")" << std::endl; 

        solve_timestep(); 
      } 

    move_mesh(); 
    output_results(); 

    pcout << std::endl; 
  } 



  template <int dim> 
  void TopLevel<dim>::do_timestep() 
  { 
    present_time += present_timestep; 
    ++timestep_no; 
    pcout << "Timestep " << timestep_no << " at time " << present_time 
          << std::endl; 
    if (present_time > end_time) 
      { 
        present_timestep -= (present_time - end_time); 
        present_time = end_time; 
      } 

    solve_timestep(); 

    move_mesh(); 
    output_results(); 

    pcout << std::endl; 
  } 


  template <int dim> 
  void TopLevel<dim>::refine_initial_grid() 
  { 


    Vector<float> error_per_cell(triangulation.n_active_cells()); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      incremental_displacement, 
      error_per_cell, 
      ComponentMask(), 
      nullptr, 
      MultithreadInfo::n_threads(), 
      this_mpi_process); 


    const unsigned int n_local_cells = 
      triangulation.n_locally_owned_active_cells(); 

    PETScWrappers::MPI::Vector distributed_error_per_cell( 
      mpi_communicator, triangulation.n_active_cells(), n_local_cells); 

    for (unsigned int i = 0; i < error_per_cell.size(); ++i) 
      if (error_per_cell(i) != 0) 
        distributed_error_per_cell(i) = error_per_cell(i); 
    distributed_error_per_cell.compress(VectorOperation::insert); 


    error_per_cell = distributed_error_per_cell; 
    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    error_per_cell, 
                                                    0.35, 
                                                    0.03); 
    triangulation.execute_coarsening_and_refinement(); 


    setup_quadrature_point_history(); 
  } 








  template <int dim> 
  void TopLevel<dim>::move_mesh() 
  { 
    pcout << "    Moving mesh..." << std::endl; 

    std::vector<bool> vertex_touched(triangulation.n_vertices(), false); 
    for (auto &cell : dof_handler.active_cell_iterators()) 
      for (const auto v : cell->vertex_indices()) 
        if (vertex_touched[cell->vertex_index(v)] == false) 
          { 
            vertex_touched[cell->vertex_index(v)] = true; 

            Point<dim> vertex_displacement; 
            for (unsigned int d = 0; d < dim; ++d) 
              vertex_displacement[d] = 
                incremental_displacement(cell->vertex_dof_index(v, d)); 

            cell->vertex(v) += vertex_displacement; 
          } 
  } 



  template <int dim> 
  void TopLevel<dim>::setup_quadrature_point_history() 
  { 


    triangulation.clear_user_data(); 


    { 
      std::vector<PointHistory<dim>> tmp; 
      quadrature_point_history.swap(tmp); 
    } 
    quadrature_point_history.resize( 
      triangulation.n_locally_owned_active_cells() * quadrature_formula.size()); 


    unsigned int history_index = 0; 
    for (auto &cell : triangulation.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell->set_user_pointer(&quadrature_point_history[history_index]); 
          history_index += quadrature_formula.size(); 
        } 


    Assert(history_index == quadrature_point_history.size(), 
           ExcInternalError()); 
  } 






  template <int dim> 
  void TopLevel<dim>::update_quadrature_point_history() 
  { 


    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients); 

    std::vector<std::vector<Tensor<1, dim>>> displacement_increment_grads( 
      quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim)); 


    for (auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 


          PointHistory<dim> *local_quadrature_points_history = 
            reinterpret_cast<PointHistory<dim> *>(cell->user_pointer()); 
          Assert(local_quadrature_points_history >= 
                   &quadrature_point_history.front(), 
                 ExcInternalError()); 
          Assert(local_quadrature_points_history <= 
                   &quadrature_point_history.back(), 
                 ExcInternalError()); 


          fe_values.reinit(cell); 
          fe_values.get_function_gradients(incremental_displacement, 
                                           displacement_increment_grads); 


          for (unsigned int q = 0; q < quadrature_formula.size(); ++q) 
            { 


              const SymmetricTensor<2, dim> new_stress = 
                (local_quadrature_points_history[q].old_stress + 
                 (stress_strain_tensor * 
                  get_strain(displacement_increment_grads[q]))); 


              const Tensor<2, dim> rotation = 
                get_rotation_matrix(displacement_increment_grads[q]); 



              const SymmetricTensor<2, dim> rotated_new_stress = 
                symmetrize(transpose(rotation) * 
                           static_cast<Tensor<2, dim>>(new_stress) * rotation); 



              local_quadrature_points_history[q].old_stress = 
                rotated_new_stress; 
            } 
        } 
  } 


} // namespace Step18 

int main(int argc, char **argv) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step18; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

      TopLevel<3> elastic_problem; 
      elastic_problem.run(); 
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


