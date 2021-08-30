

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2006 - 2021 by the deal.II authors 
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
 * Authors: Wolfgang Bangerth, Texas A&M University, 2006, 2007; 
 *          Denis Davydov, University of Erlangen-Nuremberg, 2016; 
 *          Marc Fehling, Colorado State University, 2020. 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 


#include <deal.II/hp/fe_collection.h> 
#include <deal.II/hp/fe_values.h> 
#include <deal.II/hp/refinement.h> 
#include <deal.II/fe/fe_series.h> 
#include <deal.II/numerics/smoothness_estimator.h> 


#include <fstream> 
#include <iostream> 


namespace Step27 
{ 
  using namespace dealii; 



  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(); 
    ~LaplaceProblem(); 

    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void create_coarse_grid(); 
    void postprocess(const unsigned int cycle); 

    Triangulation<dim> triangulation; 

 
    DoFHandler<dim>          dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection; 
    hp::QCollection<dim - 1> face_quadrature_collection; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    const unsigned int max_degree; 
  }; 



  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component) const override; 
  }; 

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> &p, 
                                   const unsigned int /*component*/) const 
  { 
    double product = 1; 
    for (unsigned int d = 0; d < dim; ++d) 
      product *= (p[d] + 1); 
    return product; 
  } 




  template <int dim> 
  LaplaceProblem<dim>::LaplaceProblem() 
    : dof_handler(triangulation) 
    , max_degree(dim <= 2 ? 7 : 5) 
  { 
    for (unsigned int degree = 2; degree <= max_degree; ++degree) 
      { 
        fe_collection.push_back(FE_Q<dim>(degree)); 
        quadrature_collection.push_back(QGauss<dim>(degree + 1)); 
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1)); 
      } 
  } 


  template <int dim> 
  LaplaceProblem<dim>::~LaplaceProblem() 
  { 
    dof_handler.clear(); 
  } 


  template <int dim> 
  void LaplaceProblem<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe_collection); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(), 
                                             constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
  } 




  template <int dim> 
  void LaplaceProblem<dim>::assemble_system() 
  { 
    hp::FEValues<dim> hp_fe_values(fe_collection, 
                                   quadrature_collection, 
                                   update_values | update_gradients | 
                                     update_quadrature_points | 
                                     update_JxW_values); 

    RightHandSide<dim> rhs_function; 

    FullMatrix<double> cell_matrix; 
    Vector<double>     cell_rhs; 

    std::vector<types::global_dof_index> local_dof_indices; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 

        cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
        cell_matrix = 0; 

        cell_rhs.reinit(dofs_per_cell); 
        cell_rhs = 0; 

        hp_fe_values.reinit(cell); 

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values(); 

        std::vector<double> rhs_values(fe_values.n_quadrature_points); 
        rhs_function.value_list(fe_values.get_quadrature_points(), rhs_values); 

        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points; 
             ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                cell_matrix(i, j) += 
                  (fe_values.shape_grad(i, q_point) * // grad phi_i(x_q) 
                   fe_values.shape_grad(j, q_point) * // grad phi_j(x_q) 
                   fe_values.JxW(q_point));           // dx 

              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q) 
                              rhs_values[q_point] *               // f(x_q) 
                              fe_values.JxW(q_point));            // dx 
            } 

        local_dof_indices.resize(dofs_per_cell); 
        cell->get_dof_indices(local_dof_indices); 

        constraints.distribute_local_to_global( 
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::solve() 
  { 
    SolverControl            solver_control(system_rhs.size(), 
                                 1e-12 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    constraints.distribute(solution); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::postprocess(const unsigned int cycle) 
  { 


    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      face_quadrature_collection, 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      estimated_error_per_cell); 


    Vector<float> smoothness_indicators(triangulation.n_active_cells()); 
    FESeries::Fourier<dim> fourier = 
      SmoothnessEstimator::Fourier::default_fe_series(fe_collection); 
    SmoothnessEstimator::Fourier::coefficient_decay(fourier, 
                                                    dof_handler, 
                                                    solution, 
                                                    smoothness_indicators); 



    { 
      Vector<float> fe_degrees(triangulation.n_active_cells()); 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        fe_degrees(cell->active_cell_index()) = 
          fe_collection[cell->active_fe_index()].degree; 


      DataOut<dim> data_out; 

      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution, "solution"); 
      data_out.add_data_vector(estimated_error_per_cell, "error"); 
      data_out.add_data_vector(smoothness_indicators, "smoothness"); 
      data_out.add_data_vector(fe_degrees, "fe_degree"); 
      data_out.build_patches(); 


      const std::string filename = 
        "solution-" + Utilities::int_to_string(cycle, 2) + ".vtk"; 
      std::ofstream output(filename); 
      data_out.write_vtk(output); 
    } 


    { 
      GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                      estimated_error_per_cell, 
                                                      0.3, 
                                                      0.03); 



      hp::Refinement::p_adaptivity_from_relative_threshold( 
        dof_handler, smoothness_indicators, 0.2, 0.2); 




      hp::Refinement::choose_p_over_h(dof_handler); 


      triangulation.prepare_coarsening_and_refinement(); 
      hp::Refinement::limit_p_level_difference(dof_handler); 

      triangulation.execute_coarsening_and_refinement(); 
    } 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::create_coarse_grid() 
  { 
    Triangulation<dim> cube; 
    GridGenerator::subdivided_hyper_cube(cube, 4, -1., 1.); 

    std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove; 
    for (const auto &cell : cube.active_cell_iterators()) 
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v) 
        if (cell->vertex(v).square() < .1) 
          cells_to_remove.insert(cell); 

    GridGenerator::create_triangulation_with_removed_cells(cube, 
                                                           cells_to_remove, 
                                                           triangulation); 

    triangulation.refine_global(3); 
  } 




  template <int dim> 
  void LaplaceProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          create_coarse_grid(); 

        setup_system(); 

        std::cout << "   Number of active cells      : " 
                  << triangulation.n_active_cells() << std::endl 
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl 
                  << "   Number of constraints       : " 
                  << constraints.n_constraints() << std::endl; 

        assemble_system(); 
        solve(); 
        postprocess(cycle); 
      } 
  } 
} // namespace Step27 


int main() 
{ 
  try 
    { 
      using namespace Step27; 

      LaplaceProblem<2> laplace_problem; 
      laplace_problem.run(); 
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

