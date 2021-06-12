

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
 * Author: Xing Jin, Wolfgang Bangerth, Texas A&M University, 2006 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
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

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/vector_tools.h> 

#include <fstream> 
#include <iostream> 


#include <deal.II/grid/grid_tools.h> 


namespace Step24 
{ 
  using namespace dealii; 


  template <int dim> 
  class TATForwardProblem 
  { 
  public: 
    TATForwardProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void solve_p(); 
    void solve_v(); 
    void output_results() const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 
    SparseMatrix<double> mass_matrix; 
    SparseMatrix<double> laplace_matrix; 

    Vector<double> solution_p, solution_v; 
    Vector<double> old_solution_p, old_solution_v; 
    Vector<double> system_rhs_p, system_rhs_v; 

    double       time_step, time; 
    unsigned int timestep_number; 
    const double theta; 


    SparseMatrix<double> boundary_matrix; 
    const double         wave_speed; 


    std::vector<Point<dim>> detector_locations; 
  }; 



  template <int dim> 
  class InitialValuesP : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      static const std::array<Source, 5> sources{ 
        {Source(Point<dim>(0, 0), 0.025), 
         Source(Point<dim>(-0.135, 0), 0.05), 
         Source(Point<dim>(0.17, 0), 0.03), 
         Source(Point<dim>(-0.25, 0), 0.02), 
         Source(Point<dim>(-0.05, -0.15), 0.015)}}; 

      for (const auto &source : sources) 
        if (p.distance(source.location) < source.radius) 
          return 1; 

      return 0; 
    } 

 
    struct Source 
    { 
      Source(const Point<dim> &l, const double r) 
        : location(l) 
        , radius(r) 
      {} 

      const Point<dim> location; 
      const double     radius; 
    }; 
  }; 


  template <int dim> 
  TATForwardProblem<dim>::TATForwardProblem() 
    : fe(1) 
    , dof_handler(triangulation) 
    , time_step(std::numeric_limits<double>::quiet_NaN()) 
    , time(time_step) 
    , timestep_number(1) 
    , theta(0.5) 
    , wave_speed(1.437) 
  { 



    Assert(dim == 2, ExcNotImplemented()); 

    const double detector_step_angle = 2.25; 
    const double detector_radius     = 0.5; 

    for (double detector_angle = 2 * numbers::PI; detector_angle >= 0; 
         detector_angle -= detector_step_angle / 360 * 2 * numbers::PI) 
      detector_locations.push_back( 
        Point<dim>(std::cos(detector_angle), std::sin(detector_angle)) * 
        detector_radius); 
  } 






  template <int dim> 
  void TATForwardProblem<dim>::setup_system() 
  { 
    const Point<dim> center; 
    GridGenerator::hyper_ball(triangulation, center, 1.); 
    triangulation.refine_global(7); 

    time_step = GridTools::minimal_cell_diameter(triangulation) / wave_speed / 
                std::sqrt(1. * dim); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl; 

    dof_handler.distribute_dofs(fe); 

    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl 
              << std::endl; 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
    mass_matrix.reinit(sparsity_pattern); 
    laplace_matrix.reinit(sparsity_pattern); 

    MatrixCreator::create_mass_matrix(dof_handler, 
                                      QGauss<dim>(fe.degree + 1), 
                                      mass_matrix); 
    MatrixCreator::create_laplace_matrix(dof_handler, 
                                         QGauss<dim>(fe.degree + 1), 
                                         laplace_matrix); 




    boundary_matrix.reinit(sparsity_pattern); 


    { 
      const QGauss<dim - 1> quadrature_formula(fe.degree + 1); 
      FEFaceValues<dim>     fe_values(fe, 
                                  quadrature_formula, 
                                  update_values | update_JxW_values); 

      const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
      const unsigned int n_q_points    = quadrature_formula.size(); 

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary()) 
            { 
              cell_matrix = 0; 

              fe_values.reinit(cell, face); 

              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    cell_matrix(i, j) += (fe_values.shape_value(i, q_point) * 
                                          fe_values.shape_value(j, q_point) * 
                                          fe_values.JxW(q_point)); 

              cell->get_dof_indices(local_dof_indices); 
              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                  boundary_matrix.add(local_dof_indices[i], 
                                      local_dof_indices[j], 
                                      cell_matrix(i, j)); 
            } 
    } 

    system_matrix.copy_from(mass_matrix); 
    system_matrix.add(time_step * time_step * theta * theta * wave_speed * 
                        wave_speed, 
                      laplace_matrix); 
    system_matrix.add(wave_speed * theta * time_step, boundary_matrix); 

    solution_p.reinit(dof_handler.n_dofs()); 
    old_solution_p.reinit(dof_handler.n_dofs()); 
    system_rhs_p.reinit(dof_handler.n_dofs()); 

    solution_v.reinit(dof_handler.n_dofs()); 
    old_solution_v.reinit(dof_handler.n_dofs()); 
    system_rhs_v.reinit(dof_handler.n_dofs()); 

    constraints.close(); 
  } 


  template <int dim> 
  void TATForwardProblem<dim>::solve_p() 
  { 
    SolverControl solver_control(1000, 1e-8 * system_rhs_p.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    cg.solve(system_matrix, solution_p, system_rhs_p, PreconditionIdentity()); 

    std::cout << "   p-equation: " << solver_control.last_step() 
              << " CG iterations." << std::endl; 
  } 

  template <int dim> 
  void TATForwardProblem<dim>::solve_v() 
  { 
    SolverControl solver_control(1000, 1e-8 * system_rhs_v.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    cg.solve(mass_matrix, solution_v, system_rhs_v, PreconditionIdentity()); 

    std::cout << "   v-equation: " << solver_control.last_step() 
              << " CG iterations." << std::endl; 
  } 



  template <int dim> 
  void TATForwardProblem<dim>::output_results() const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution_p, "P"); 
    data_out.add_data_vector(solution_v, "V"); 

    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 
    DataOutBase::VtkFlags vtk_flags; 
    vtk_flags.compression_level = 
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 




  template <int dim> 
  void TATForwardProblem<dim>::run() 
  { 
    setup_system(); 

    VectorTools::project(dof_handler, 
                         constraints, 
                         QGauss<dim>(fe.degree + 1), 
                         InitialValuesP<dim>(), 
                         old_solution_p); 
    old_solution_v = 0; 

    std::ofstream detector_data("detectors.dat"); 

    Vector<double> tmp(solution_p.size()); 
    Vector<double> G1(solution_p.size()); 
    Vector<double> G2(solution_v.size()); 

    const double end_time = 0.7; 
    for (time = time_step; time <= end_time; 
         time += time_step, ++timestep_number) 
      { 
        std::cout << std::endl; 
        std::cout << "time_step " << timestep_number << " @ t=" << time 
                  << std::endl; 

        mass_matrix.vmult(G1, old_solution_p); 
        mass_matrix.vmult(tmp, old_solution_v); 
        G1.add(time_step * (1 - theta), tmp); 

        mass_matrix.vmult(G2, old_solution_v); 
        laplace_matrix.vmult(tmp, old_solution_p); 
        G2.add(-wave_speed * wave_speed * time_step * (1 - theta), tmp); 

        boundary_matrix.vmult(tmp, old_solution_p); 
        G2.add(wave_speed, tmp); 

        system_rhs_p = G1; 
        system_rhs_p.add(time_step * theta, G2); 

        solve_p(); 

        system_rhs_v = G2; 
        laplace_matrix.vmult(tmp, solution_p); 
        system_rhs_v.add(-time_step * theta * wave_speed * wave_speed, tmp); 

        boundary_matrix.vmult(tmp, solution_p); 
        system_rhs_v.add(-wave_speed, tmp); 

        solve_v(); 

        output_results(); 

        detector_data << time; 
        for (unsigned int i = 0; i < detector_locations.size(); ++i) 
          detector_data << " " 
                        << VectorTools::point_value(dof_handler, 
                                                    solution_p, 
                                                    detector_locations[i]) 
                        << " "; 
        detector_data << std::endl; 

        old_solution_p = solution_p; 
        old_solution_v = solution_v; 
      } 
  } 
} // namespace Step24 



int main() 
{ 
  try 
    { 
      using namespace Step24; 

      TATForwardProblem<2> forward_problem_solver; 
      forward_problem_solver.run(); 
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


