

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2001 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 2001 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/table_handler.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 


#include <deal.II/lac/dynamic_sparsity_pattern.h> 


#include <algorithm> 
#include <iostream> 
#include <iomanip> 
#include <cmath> 


namespace Step11 
{ 
  using namespace dealii; 



  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(const unsigned int mapping_degree); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_and_solve(); 
    void solve(); 
    void write_high_order_mesh(const unsigned cycle); 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 
    MappingQ<dim>      mapping; 

    SparsityPattern           sparsity_pattern; 
    SparseMatrix<double>      system_matrix; 
    AffineConstraints<double> mean_value_constraints; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    TableHandler output_table; 
  }; 


  template <int dim> 
  LaplaceProblem<dim>::LaplaceProblem(const unsigned int mapping_degree) 
    : fe(1) 
    , dof_handler(triangulation) 
    , mapping(mapping_degree) 
  { 
    std::cout << "Using mapping with degree " << mapping_degree << ":" 
              << std::endl 
              << "============================" << std::endl; 
  } 


  template <int dim> 
  void LaplaceProblem<dim>::setup_system() 
  { 


    dof_handler.distribute_dofs(fe); 
    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 



    const IndexSet boundary_dofs = DoFTools::extract_boundary_dofs(dof_handler); 

    const types::global_dof_index first_boundary_dof = 
      boundary_dofs.nth_index_in_set(0); 


    mean_value_constraints.clear(); 
    mean_value_constraints.add_line(first_boundary_dof); 
    for (types::global_dof_index i : boundary_dofs) 
      if (i != first_boundary_dof) 
        mean_value_constraints.add_entry(first_boundary_dof, i, -1); 
    mean_value_constraints.close(); 



    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    mean_value_constraints.condense(dsp); 



    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 
  } 


  template <int dim> 
  void LaplaceProblem<dim>::assemble_and_solve() 
  { 




    const unsigned int gauss_degree = 
      std::max(static_cast<unsigned int>( 
                 std::ceil(1. * (mapping.get_degree() + 1) / 2)), 
               2U); 
    MatrixTools::create_laplace_matrix(mapping, 
                                       dof_handler, 
                                       QGauss<dim>(gauss_degree), 
                                       system_matrix); 
    VectorTools::create_right_hand_side(mapping, 
                                        dof_handler, 
                                        QGauss<dim>(gauss_degree), 
                                        Functions::ConstantFunction<dim>(-2), 
                                        system_rhs); 





    Vector<double> tmp(system_rhs.size()); 
    VectorTools::create_boundary_right_hand_side( 
      mapping, 
      dof_handler, 
      QGauss<dim - 1>(gauss_degree), 
      Functions::ConstantFunction<dim>(1), 
      tmp); 


    system_rhs += tmp; 



    mean_value_constraints.condense(system_matrix); 
    mean_value_constraints.condense(system_rhs); 

    solve(); 
    mean_value_constraints.distribute(solution); 


    Vector<float> norm_per_cell(triangulation.n_active_cells()); 
    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      solution, 
                                      Functions::ZeroFunction<dim>(), 
                                      norm_per_cell, 
                                      QGauss<dim>(gauss_degree + 1), 
                                      VectorTools::H1_seminorm); 


    const double norm = 
      VectorTools::compute_global_error(triangulation, 
                                        norm_per_cell, 
                                        VectorTools::H1_seminorm); 


    output_table.add_value("cells", triangulation.n_active_cells()); 
    output_table.add_value("|u|_1", norm); 
    output_table.add_value("error", 
                           std::fabs(norm - std::sqrt(3.14159265358 / 2))); 
  } 


  template <int dim> 
  void LaplaceProblem<dim>::solve() 
  { 
    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::write_high_order_mesh(const unsigned cycle) 
  { 
    DataOut<dim> data_out; 

    DataOutBase::VtkFlags flags; 
    flags.write_higher_order_cells = true; 
    data_out.set_flags(flags); 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 

    data_out.build_patches(mapping, 
                           mapping.get_degree(), 
                           DataOut<dim>::curved_inner_cells); 

    std::ofstream file("solution-c=" + std::to_string(cycle) + 
                       ".p=" + std::to_string(mapping.get_degree()) + ".vtu"); 

    data_out.write_vtu(file); 
  } 


  template <int dim> 
  void LaplaceProblem<dim>::run() 
  { 
    GridGenerator::hyper_ball(triangulation); 

    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      { 
        setup_system(); 
        assemble_and_solve(); 
        write_high_order_mesh(cycle); 

        triangulation.refine_global(); 
      } 


    output_table.set_precision("|u|_1", 6); 
    output_table.set_precision("error", 6); 
    output_table.write_text(std::cout); 
    std::cout << std::endl; 
  } 
} // namespace Step11 


int main() 
{ 
  try 
    { 
      std::cout.precision(5); 


      for (unsigned int mapping_degree = 1; mapping_degree <= 3; 
           ++mapping_degree) 
        Step11::LaplaceProblem<2>(mapping_degree).run(); 
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
    }; 

  return 0; 
} 


