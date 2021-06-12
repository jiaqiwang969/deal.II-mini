

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2021 by the deal.II authors 
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
 * Author: Tyler Anderson, Colorado State University, 2021 
 */ 




#include <deal.II/base/convergence_table.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_accessor.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/tria_accessor.h> 
#include <deal.II/grid/tria_iterator.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/data_out_stack.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/solution_transfer.h> 
#include <deal.II/numerics/vector_tools.h> 

#include <fstream> 
#include <iostream> 


namespace BlackScholesSolver 
{ 
  using namespace dealii; 

#define MMS 


  template <int dim> 
  class Solution : public Function<dim> 
  { 
  public: 
    Solution(const double maturity_time); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual Tensor<1, dim> 
    gradient(const Point<dim> & p, 
             const unsigned int component = 0) const override; 

  private: 
    const double maturity_time; 
  }; 

  template <int dim> 
  Solution<dim>::Solution(const double maturity_time) 
    : maturity_time(maturity_time) 
  { 
    Assert(dim == 1, ExcNotImplemented()); 
  } 

  template <int dim> 
  double Solution<dim>::value(const Point<dim> & p, 
                              const unsigned int component) const 
  { 
    return -Utilities::fixed_power<2, double>(p(component)) - 
           Utilities::fixed_power<2, double>(this->get_time()) + 6; 
  } 

  template <int dim> 
  Tensor<1, dim> Solution<dim>::gradient(const Point<dim> & p, 
                                         const unsigned int component) const 
  { 
    return Point<dim>(-2 * p(component)); 
  } 




  template <int dim> 
  class InitialConditions : public Function<dim> 
  { 
  public: 
    InitialConditions(const double strike_price); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

  private: 
    const double strike_price; 
  }; 

  template <int dim> 
  InitialConditions<dim>::InitialConditions(const double strike_price) 
    : strike_price(strike_price) 
  {} 

  template <int dim> 
  double InitialConditions<dim>::value(const Point<dim> & p, 
                                       const unsigned int component) const 
  { 
#ifdef MMS 
    return -Utilities::fixed_power<2, double>(p(component)) + 6; 
#else 
    return std::max(p(component) - strike_price, 0.); 
#endif 
  } 


  template <int dim> 
  class LeftBoundaryValues : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double LeftBoundaryValues<dim>::value(const Point<dim> &, 
                                        const unsigned int /*component*/) const 
  { 
#ifdef MMS 
    return -Utilities::fixed_power<2, double>(this->get_time()) + 6; 
#else 
    return 0.; 
#endif 
  } 


  template <int dim> 
  class RightBoundaryValues : public Function<dim> 
  { 
  public: 
    RightBoundaryValues(const double strike_price, const double interest_rate); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

  private: 
    const double strike_price; 
    const double interest_rate; 
  }; 

  template <int dim> 
  RightBoundaryValues<dim>::RightBoundaryValues(const double strike_price, 
                                                const double interest_rate) 
    : strike_price(strike_price) 
    , interest_rate(interest_rate) 
  {} 

  template <int dim> 
  double RightBoundaryValues<dim>::value(const Point<dim> & p, 
                                         const unsigned int component) const 
  { 
#ifdef MMS 
    return -Utilities::fixed_power<2, double>(p(component)) - 
           Utilities::fixed_power<2, double>(this->get_time()) + 6; 
#else 
    return (p(component) - strike_price) * 
           exp((-interest_rate) * (this->get_time())); 
#endif 
  } 


  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    RightHandSide(const double asset_volatility, const double interest_rate); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

  private: 
    const double asset_volatility; 
    const double interest_rate; 
  }; 

  template <int dim> 
  RightHandSide<dim>::RightHandSide(const double asset_volatility, 
                                    const double interest_rate) 
    : asset_volatility(asset_volatility) 
    , interest_rate(interest_rate) 
  {} 

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> & p, 
                                   const unsigned int component) const 
  { 
#ifdef MMS 
    return 2 * (this->get_time()) - 
           Utilities::fixed_power<2, double>(asset_volatility * p(component)) - 
           2 * interest_rate * Utilities::fixed_power<2, double>(p(component)) - 
           interest_rate * 
             (-Utilities::fixed_power<2, double>(p(component)) - 
              Utilities::fixed_power<2, double>(this->get_time()) + 6); 
#else 
    (void)p; 
    (void)component; 
    return 0.0; 
#endif 
  } 












  template <int dim> 
  class BlackScholes 
  { 
  public: 
    BlackScholes(); 

    void run(); 

  private: 
    void setup_system(); 
    void solve_time_step(); 
    void refine_grid(); 
    void process_solution(); 
    void add_results_for_output(); 
    void write_convergence_table(); 

    const double maximum_stock_price; 
    const double maturity_time; 
    const double asset_volatility; 
    const double interest_rate; 
    const double strike_price; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> mass_matrix; 
    SparseMatrix<double> laplace_matrix; 
    SparseMatrix<double> a_matrix; 
    SparseMatrix<double> b_matrix; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    double time; 
    double time_step; 
    
    const double       theta;
    const unsigned int n_cycles;
    const unsigned int n_time_steps;

    DataOutStack<dim>        data_out_stack; 
    std::vector<std::string> solution_names; 

    ConvergenceTable convergence_table; 
  }; 



  template <int dim> 
  BlackScholes<dim>::BlackScholes() 
    : maximum_stock_price(1.) 
    , maturity_time(1.) 
    , asset_volatility(.2) 
    , interest_rate(0.05) 
    , strike_price(0.5) 
    , fe(1) 
    , dof_handler(triangulation) 
    , time(0.0) 
    , theta(0.5) 
    , n_cycles(4) 
    , n_time_steps(5000) 
  { 
    Assert(dim == 1, ExcNotImplemented()); 
  } 



  template <int dim> 
  void BlackScholes<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 

    time_step = maturity_time / n_time_steps; 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    constraints.close(); 
    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    dsp, 
                                    constraints, 
                                    /*keep_constrained_dofs =  */ true);
    sparsity_pattern.copy_from(dsp); 

    mass_matrix.reinit(sparsity_pattern); 
    laplace_matrix.reinit(sparsity_pattern); 
    a_matrix.reinit(sparsity_pattern); 
    b_matrix.reinit(sparsity_pattern); 
    system_matrix.reinit(sparsity_pattern); 

    MatrixCreator::create_mass_matrix(dof_handler, 
                                      QGauss<dim>(fe.degree + 1), 
                                      mass_matrix); 


    const unsigned int dofs_per_cell = fe.dofs_per_cell; 
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    QGauss<dim>        quadrature_formula(fe.degree + 1); 
    FEValues<dim>      fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0.; 
        fe_values.reinit(cell); 
        for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
          { 
            const double current_coefficient = 
              fe_values.quadrature_point(q_index).square(); 
            for (const unsigned int i : fe_values.dof_indices()) 
              { 
                for (const unsigned int j : fe_values.dof_indices()) 
                  cell_matrix(i, j) += 
                    (current_coefficient *              // (x_q)^2 
                     fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
                     fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
                     fe_values.JxW(q_index));           // dx 
              } 
          } 
        cell->get_dof_indices(local_dof_indices); 
        for (const unsigned int i : fe_values.dof_indices()) 
          { 
            for (const unsigned int j : fe_values.dof_indices()) 
              laplace_matrix.add(local_dof_indices[i], 
                                 local_dof_indices[j], 
                                 cell_matrix(i, j)); 
          } 
      } 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0.; 
        fe_values.reinit(cell); 
        for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
          { 
            const Tensor<1, dim> current_coefficient = 
              fe_values.quadrature_point(q_index); 
            for (const unsigned int i : fe_values.dof_indices()) 
              { 
                for (const unsigned int j : fe_values.dof_indices()) 
                  { 
                    cell_matrix(i, j) += 
                      (current_coefficient *               // x_q 
                       fe_values.shape_grad(i, q_index) *  // grad phi_i(x_q) 
                       fe_values.shape_value(j, q_index) * // phi_j(x_q) 
                       fe_values.JxW(q_index));            // dx 
                  } 
              } 
          } 
        cell->get_dof_indices(local_dof_indices); 
        for (const unsigned int i : fe_values.dof_indices()) 
          { 
            for (const unsigned int j : fe_values.dof_indices()) 
              a_matrix.add(local_dof_indices[i], 
                           local_dof_indices[j], 
                           cell_matrix(i, j)); 
          } 
      } 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0.; 
        fe_values.reinit(cell); 
        for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
          { 
            const Tensor<1, dim> current_coefficient = 
              fe_values.quadrature_point(q_index); 
            for (const unsigned int i : fe_values.dof_indices()) 
              { 
                for (const unsigned int j : fe_values.dof_indices()) 
                  cell_matrix(i, j) += 
                    (current_coefficient *               // x_q 
                     fe_values.shape_value(i, q_index) * // phi_i(x_q) 
                     fe_values.shape_grad(j, q_index) *  // grad phi_j(x_q) 
                     fe_values.JxW(q_index));            // dx 
              } 
          } 
        cell->get_dof_indices(local_dof_indices); 
        for (const unsigned int i : fe_values.dof_indices()) 
          { 
            for (const unsigned int j : fe_values.dof_indices()) 
              b_matrix.add(local_dof_indices[i], 
                           local_dof_indices[j], 
                           cell_matrix(i, j)); 
          } 
      } 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 


  template <int dim> 
  void BlackScholes<dim>::solve_time_step() 
  { 
    SolverControl                          solver_control(1000, 1e-12); 
    SolverCG<Vector<double>>               cg(solver_control); 
    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.0); 
    cg.solve(system_matrix, solution, system_rhs, preconditioner); 
    constraints.distribute(solution); 
  } 


  template <int dim> 
  void BlackScholes<dim>::add_results_for_output() 
  { 
    data_out_stack.new_parameter_value(time, time_step); 
    data_out_stack.attach_dof_handler(dof_handler); 
    data_out_stack.add_data_vector(solution, solution_names); 
    data_out_stack.build_patches(2); 
    data_out_stack.finish_parameter_value(); 
  } 


  template <int dim> 
  void BlackScholes<dim>::refine_grid() 
  { 
    triangulation.refine_global(1); 
  } 


  template <int dim> 
  void BlackScholes<dim>::process_solution() 
  { 
    Solution<dim> sol(maturity_time); 
    sol.set_time(time); 
    Vector<float> difference_per_cell(triangulation.n_active_cells()); 
    VectorTools::integrate_difference(dof_handler, 
                                      solution, 
                                      sol, 
                                      difference_per_cell, 
                                      QGauss<dim>(fe.degree + 1), 
                                      VectorTools::L2_norm); 
    const double L2_error = 
      VectorTools::compute_global_error(triangulation, 
                                        difference_per_cell, 
                                        VectorTools::L2_norm); 
    VectorTools::integrate_difference(dof_handler, 
                                      solution, 
                                      sol, 
                                      difference_per_cell, 
                                      QGauss<dim>(fe.degree + 1), 
                                      VectorTools::H1_seminorm); 
    const double H1_error = 
      VectorTools::compute_global_error(triangulation, 
                                        difference_per_cell, 
                                        VectorTools::H1_seminorm); 
    const QTrapezoid<1>  q_trapezoid; 
    const QIterated<dim> q_iterated(q_trapezoid, fe.degree * 2 + 1); 
    VectorTools::integrate_difference(dof_handler, 
                                      solution, 
                                      sol, 
                                      difference_per_cell, 
                                      q_iterated, 
                                      VectorTools::Linfty_norm); 
    const double Linfty_error = 
      VectorTools::compute_global_error(triangulation, 
                                        difference_per_cell, 
                                        VectorTools::Linfty_norm); 
    const unsigned int n_active_cells = triangulation.n_active_cells(); 
    const unsigned int n_dofs         = dof_handler.n_dofs(); 
    convergence_table.add_value("cells", n_active_cells); 
    convergence_table.add_value("dofs", n_dofs); 
    convergence_table.add_value("L2", L2_error); 
    convergence_table.add_value("H1", H1_error); 
    convergence_table.add_value("Linfty", Linfty_error); 
  } 


  template <int dim> 
  void BlackScholes<dim>::write_convergence_table() 
  { 
    convergence_table.set_precision("L2", 3); 
    convergence_table.set_precision("H1", 3); 
    convergence_table.set_precision("Linfty", 3); 
    convergence_table.set_scientific("L2", true); 
    convergence_table.set_scientific("H1", true); 
    convergence_table.set_scientific("Linfty", true); 
    convergence_table.set_tex_caption("cells", "\\# cells"); 
    convergence_table.set_tex_caption("dofs", "\\# dofs"); 
    convergence_table.set_tex_caption("L2", "@f$L^2@f$-error"); 
    convergence_table.set_tex_caption("H1", "@f$H^1@f$-error"); 
    convergence_table.set_tex_caption("Linfty", "@f$L^\\infty@f$-error"); 
    convergence_table.set_tex_format("cells", "r"); 
    convergence_table.set_tex_format("dofs", "r"); 
    std::cout << std::endl; 
    convergence_table.write_text(std::cout); 
    std::string error_filename = "error"; 
    error_filename += "-global"; 
    error_filename += ".tex"; 
    std::ofstream error_table_file(error_filename); 
    convergence_table.write_tex(error_table_file); 


    convergence_table.add_column_to_supercolumn("cells", "n cells"); 
    std::vector<std::string> new_order; 
    new_order.emplace_back("n cells"); 
    new_order.emplace_back("H1"); 
    new_order.emplace_back("L2"); 
    convergence_table.set_column_order(new_order); 
    convergence_table.evaluate_convergence_rates( 
      "L2", ConvergenceTable::reduction_rate); 
    convergence_table.evaluate_convergence_rates( 
      "L2", ConvergenceTable::reduction_rate_log2); 
    convergence_table.evaluate_convergence_rates( 
      "H1", ConvergenceTable::reduction_rate); 
    convergence_table.evaluate_convergence_rates( 
      "H1", ConvergenceTable::reduction_rate_log2); 
    std::cout << std::endl; 
    convergence_table.write_text(std::cout); 
    std::string conv_filename = "convergence"; 
    conv_filename += "-global"; 
    switch (fe.degree) 
      { 
        case 1: 
          conv_filename += "-q1"; 
          break; 
        case 2: 
          conv_filename += "-q2"; 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
      } 
    conv_filename += ".tex"; 
    std::ofstream table_file(conv_filename); 
    convergence_table.write_tex(table_file); 
  } 


  template <int dim> 
  void BlackScholes<dim>::run() 
  { 
    GridGenerator::hyper_cube(triangulation, 0.0, maximum_stock_price, true); 
    triangulation.refine_global(0); 

    solution_names.emplace_back("u"); 
    data_out_stack.declare_data_vector(solution_names, 
                                       DataOutStack<dim>::dof_vector); 

    Vector<double> vmult_result; 
    Vector<double> forcing_terms; 

    for (unsigned int cycle = 0; cycle < n_cycles; cycle++) 
      { 
        if (cycle != 0) 
          { 
            refine_grid(); 
            time = 0.0; 
          } 

        setup_system(); 

        std::cout << std::endl 
                  << "===========================================" << std::endl 
                  << "Cycle " << cycle << ':' << std::endl 
                  << "Number of active cells: " 
                  << triangulation.n_active_cells() << std::endl 
                  << "Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl 
                  << std::endl; 

        VectorTools::interpolate(dof_handler, 
                                 InitialConditions<dim>(strike_price), 
                                 solution); 

        if (cycle == (n_cycles - 1)) 
          { 
            add_results_for_output(); 
          } 


        vmult_result.reinit(dof_handler.n_dofs()); 
        forcing_terms.reinit(dof_handler.n_dofs()); 
        for (unsigned int timestep_number = 0; timestep_number < n_time_steps; 
             ++timestep_number) 
          { 
            time += time_step; 

            if (timestep_number % 1000 == 0) 
              std::cout << "Time step " << timestep_number << " at t=" << time 
                        << std::endl; 

            mass_matrix.vmult(system_rhs, solution); 

            laplace_matrix.vmult(vmult_result, solution); 
            system_rhs.add( 
              (-1) * (1 - theta) * time_step * 
                Utilities::fixed_power<2, double>(asset_volatility) * 0.5, 
              vmult_result); 
            mass_matrix.vmult(vmult_result, solution); 

            system_rhs.add((-1) * (1 - theta) * time_step * interest_rate * 2, 
                           vmult_result); 

            a_matrix.vmult(vmult_result, solution); 
            system_rhs.add((-1) * time_step * interest_rate, vmult_result); 

            b_matrix.vmult(vmult_result, solution); 
            system_rhs.add( 
              (-1) * Utilities::fixed_power<2, double>(asset_volatility) * 
                time_step * 1, 
              vmult_result); 


            RightHandSide<dim> rhs_function(asset_volatility, interest_rate); 
            rhs_function.set_time(time); 
            VectorTools::create_right_hand_side(dof_handler, 
                                                QGauss<dim>(fe.degree + 1), 
                                                rhs_function, 
                                                forcing_terms); 
            forcing_terms *= time_step * theta; 
            system_rhs -= forcing_terms; 

            rhs_function.set_time(time - time_step); 
            VectorTools::create_right_hand_side(dof_handler, 
                                                QGauss<dim>(fe.degree + 1), 
                                                rhs_function, 
                                                forcing_terms); 
            forcing_terms *= time_step * (1 - theta); 
            system_rhs -= forcing_terms; 


            system_matrix.copy_from(mass_matrix); 
            system_matrix.add( 
              (theta)*time_step * 
                Utilities::fixed_power<2, double>(asset_volatility) * 0.5, 
              laplace_matrix); 
            system_matrix.add((time_step)*interest_rate * theta * (1 + 1), 
                              mass_matrix); 

            constraints.condense(system_matrix, system_rhs); 


            { 
              RightBoundaryValues<dim> right_boundary_function(strike_price, 
                                                               interest_rate); 
              LeftBoundaryValues<dim>  left_boundary_function; 
              right_boundary_function.set_time(time); 
              left_boundary_function.set_time(time); 
              std::map<types::global_dof_index, double> boundary_values; 
              VectorTools::interpolate_boundary_values(dof_handler, 
                                                       0, 
                                                       left_boundary_function, 
                                                       boundary_values); 
              VectorTools::interpolate_boundary_values(dof_handler, 
                                                       1, 
                                                       right_boundary_function, 
                                                       boundary_values); 
              MatrixTools::apply_boundary_values(boundary_values, 
                                                 system_matrix, 
                                                 solution, 
                                                 system_rhs); 
            } 


            solve_time_step(); 

            if (cycle == (n_cycles - 1)) 
              { 
                add_results_for_output(); 
              } 
          } 
#ifdef MMS 
        process_solution(); 
#endif 
      } 

    const std::string filename = "solution.vtk"; 
    std::ofstream     output(filename); 
    data_out_stack.write_vtk(output); 

#ifdef MMS 
    write_convergence_table(); 
#endif 
  } 

} // namespace BlackScholesSolver 


int main() 
{ 
  try 
    { 
      using namespace BlackScholesSolver; 

      BlackScholes<1> black_scholes_solver; 
      black_scholes_solver.run(); 
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

