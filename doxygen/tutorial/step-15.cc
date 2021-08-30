

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2012 - 2021 by the deal.II authors 
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
 * Author: Sven Wetterauer, University of Heidelberg, 2012 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <fstream> 
#include <iostream> 


#include <deal.II/numerics/solution_transfer.h> 


namespace Step15 
{ 
  using namespace dealii; 





  template <int dim> 
  class MinimalSurfaceProblem 
  { 
  public: 
    MinimalSurfaceProblem(); 
    void run(); 

  private: 
    void   setup_system(const bool initial_step); 
    void   assemble_system(); 
    void   solve(); 
    void   refine_mesh(); 
    void   set_boundary_values(); 
    double compute_residual(const double alpha) const; 
    double determine_step_length() const; 
    void   output_results(const unsigned int refinement_cycle) const; 

    Triangulation<dim> triangulation; 

    DoFHandler<dim> dof_handler; 
    FE_Q<dim>       fe; 

    AffineConstraints<double> hanging_node_constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> current_solution; 
    Vector<double> newton_update; 
    Vector<double> system_rhs; 
  }; 


  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> &p, 
                                    const unsigned int /*component*/) const 
  { 
    return std::sin(2 * numbers::PI * (p[0] + p[1])); 
  } 


  template <int dim> 
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem() 
    : dof_handler(triangulation) 
    , fe(2) 
  {} 


  template <int dim> 
  void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step) 
  { 
    if (initial_step) 
      { 
        dof_handler.distribute_dofs(fe); 
        current_solution.reinit(dof_handler.n_dofs()); 

        hanging_node_constraints.clear(); 
        DoFTools::make_hanging_node_constraints(dof_handler, 
                                                hanging_node_constraints); 
        hanging_node_constraints.close(); 
      } 


    newton_update.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

    hanging_node_constraints.condense(dsp); 

    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 
  } 



  template <int dim> 
  void MinimalSurfaceProblem<dim>::assemble_system() 
  { 
    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    system_matrix = 0; 
    system_rhs    = 0; 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_gradients | update_quadrature_points | 
                              update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        fe_values.reinit(cell); 


        fe_values.get_function_gradients(current_solution, 
                                         old_solution_gradients); 


        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double coeff = 
              1.0 / std::sqrt(1 + old_solution_gradients[q] * 
                                    old_solution_gradients[q]); 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              { 
                for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                  cell_matrix(i, j) += 
                    (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i 
                       * coeff                         //   * a_n 
                       * fe_values.shape_grad(j, q))   //   * \nabla \phi_j) 
                      -                                //  - 
                      (fe_values.shape_grad(i, q)      //  (\nabla \phi_i 
                       * coeff * coeff * coeff         //   * a_n^3 
                       * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j 
                          * old_solution_gradients[q]) //      * \nabla u_n) 
                       * old_solution_gradients[q]))   //   * \nabla u_n))) 
                     * fe_values.JxW(q));              // * dx 

                cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i 
                                * coeff                     // * a_n 
                                * old_solution_gradients[q] // * u_n 
                                * fe_values.JxW(q));        // * dx 
              } 
          } 

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          { 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              system_matrix.add(local_dof_indices[i], 
                                local_dof_indices[j], 
                                cell_matrix(i, j)); 

            system_rhs(local_dof_indices[i]) += cell_rhs(i); 
          } 
      } 


    hanging_node_constraints.condense(system_matrix); 
    hanging_node_constraints.condense(system_rhs); 

    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(), 
                                             boundary_values); 
    MatrixTools::apply_boundary_values(boundary_values, 
                                       system_matrix, 
                                       newton_update, 
                                       system_rhs); 
  } 



  template <int dim> 
  void MinimalSurfaceProblem<dim>::solve() 
  { 
    SolverControl            solver_control(system_rhs.size(), 
                                 system_rhs.l2_norm() * 1e-6); 
    SolverCG<Vector<double>> solver(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    solver.solve(system_matrix, newton_update, system_rhs, preconditioner); 

    hanging_node_constraints.distribute(newton_update); 

    const double alpha = determine_step_length(); 
    current_solution.add(alpha, newton_update); 
  } 


  template <int dim> 
  void MinimalSurfaceProblem<dim>::refine_mesh() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      current_solution, 
      estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.03); 


    triangulation.prepare_coarsening_and_refinement(); 


    SolutionTransfer<dim> solution_transfer(dof_handler); 
    solution_transfer.prepare_for_coarsening_and_refinement(current_solution); 

    triangulation.execute_coarsening_and_refinement(); 

    dof_handler.distribute_dofs(fe); 


    Vector<double> tmp(dof_handler.n_dofs()); 
    solution_transfer.interpolate(current_solution, tmp); 
    current_solution = tmp; 


    hanging_node_constraints.clear(); 

    DoFTools::make_hanging_node_constraints(dof_handler, 
                                            hanging_node_constraints); 
    hanging_node_constraints.close(); 


    set_boundary_values(); 


    setup_system(false); 
  } 




  template <int dim> 
  void MinimalSurfaceProblem<dim>::set_boundary_values() 
  { 
    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             BoundaryValues<dim>(), 
                                             boundary_values); 
    for (auto &boundary_value : boundary_values) 
      current_solution(boundary_value.first) = boundary_value.second; 

    hanging_node_constraints.distribute(current_solution); 
  } 



  template <int dim> 
  double MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const 
  { 
    Vector<double> residual(dof_handler.n_dofs()); 

    Vector<double> evaluation_point(dof_handler.n_dofs()); 
    evaluation_point = current_solution; 
    evaluation_point.add(alpha, newton_update); 

    const QGauss<dim> quadrature_formula(fe.degree + 1); 
    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_gradients | update_quadrature_points | 
                              update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double>              cell_residual(dofs_per_cell); 
    std::vector<Tensor<1, dim>> gradients(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_residual = 0; 
        fe_values.reinit(cell); 


        fe_values.get_function_gradients(evaluation_point, gradients); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double coeff = 
              1. / std::sqrt(1 + gradients[q] * gradients[q]); 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i 
                                   * coeff                    // * a_n 
                                   * gradients[q]             // * u_n 
                                   * fe_values.JxW(q));       // * dx 
          } 

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          residual(local_dof_indices[i]) += cell_residual(i); 
      } 


    hanging_node_constraints.condense(residual); 

    for (types::global_dof_index i : 
         DoFTools::extract_boundary_dofs(dof_handler)) 
      residual(i) = 0; 


    return residual.l2_norm(); 
  } 




  template <int dim> 
  double MinimalSurfaceProblem<dim>::determine_step_length() const 
  { 
    return 0.1; 
  } 



  template <int dim> 
  void MinimalSurfaceProblem<dim>::output_results( 
    const unsigned int refinement_cycle) const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(current_solution, "solution"); 
    data_out.add_data_vector(newton_update, "update"); 
    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu"; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 




  template <int dim> 
  void MinimalSurfaceProblem<dim>::run() 
  { 
    GridGenerator::hyper_ball(triangulation); 
    triangulation.refine_global(2); 

    setup_system(/*first time=*/true); 
    set_boundary_values(); 


    double       last_residual_norm = std::numeric_limits<double>::max(); 
    unsigned int refinement_cycle   = 0; 
    do 
      { 
        std::cout << "Mesh refinement step " << refinement_cycle << std::endl; 

        if (refinement_cycle != 0) 
          refine_mesh(); 




        std::cout << "  Initial residual: " << compute_residual(0) << std::endl; 

        for (unsigned int inner_iteration = 0; inner_iteration < 5; 
             ++inner_iteration) 
          { 
            assemble_system(); 
            last_residual_norm = system_rhs.l2_norm(); 

            solve(); 

            std::cout << "  Residual: " << compute_residual(0) << std::endl; 
          } 

        output_results(refinement_cycle); 

        ++refinement_cycle; 
        std::cout << std::endl; 
      } 
    while (last_residual_norm > 1e-3); 
  } 
} // namespace Step15 


int main() 
{ 
  try 
    { 
      using namespace Step15; 

      MinimalSurfaceProblem<2> laplace_problem_2d; 
      laplace_problem_2d.run(); 
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



