

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2009 - 2021 by the deal.II authors 
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
 * Author: Guido Kanschat, Texas A&M University, 2009 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/fe/mapping_q1.h> 


#include <deal.II/fe/fe_dgq.h> 


#include <deal.II/lac/solver_richardson.h> 
#include <deal.II/lac/precondition_block.h> 


#include <deal.II/numerics/derivative_approximation.h> 


#include <deal.II/meshworker/dof_info.h> 
#include <deal.II/meshworker/integration_info.h> 
#include <deal.II/meshworker/simple.h> 
#include <deal.II/meshworker/loop.h> 


#include <iostream> 
#include <fstream> 

namespace Step12 
{ 
  using namespace dealii; 


  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    BoundaryValues() = default; 
    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int component = 0) const override; 
  }; 


  template <int dim> 
  void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points, 
                                       std::vector<double> &          values, 
                                       const unsigned int component) const 
  { 
    (void)component; 
    AssertIndexRange(component, 1); 
    Assert(values.size() == points.size(), 
           ExcDimensionMismatch(values.size(), points.size())); 

    for (unsigned int i = 0; i < values.size(); ++i) 
      { 
        if (points[i](0) < 0.5) 
          values[i] = 1.; 
        else 
          values[i] = 0.; 
      } 
  } 


  template <int dim> 
  Tensor<1, dim> beta(const Point<dim> &p) 
  { 
    Assert(dim >= 2, ExcNotImplemented()); 

    Tensor<1, dim> wind_field; 
    wind_field[0] = -p[1]; 
    wind_field[1] = p[0]; 
    wind_field /= wind_field.norm(); 

    return wind_field; 
  } 



  template <int dim> 
  class AdvectionProblem 
  { 
  public: 
    AdvectionProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(Vector<double> &solution); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim>   triangulation; 
    const MappingQ1<dim> mapping; 


    FE_DGQ<dim>     fe; 
    DoFHandler<dim> dof_handler; 


    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> right_hand_side; 


    using DoFInfo  = MeshWorker::DoFInfo<dim>; 
    using CellInfo = MeshWorker::IntegrationInfo<dim>; 



    static void integrate_cell_term(DoFInfo &dinfo, CellInfo &info); 
    static void integrate_boundary_term(DoFInfo &dinfo, CellInfo &info); 
    static void integrate_face_term(DoFInfo & dinfo1, 
                                    DoFInfo & dinfo2, 
                                    CellInfo &info1, 
                                    CellInfo &info2); 
  }; 


  template <int dim> 
  AdvectionProblem<dim>::AdvectionProblem() 
    : mapping() 
    , fe(1) 
    , dof_handler(triangulation) 
  {} 

  template <int dim> 
  void AdvectionProblem<dim>::setup_system() 
  { 


    dof_handler.distribute_dofs(fe); 


    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 


    system_matrix.reinit(sparsity_pattern); 
    solution.reinit(dof_handler.n_dofs()); 
    right_hand_side.reinit(dof_handler.n_dofs()); 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::assemble_system() 
  { 


    MeshWorker::IntegrationInfoBox<dim> info_box; 


    const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1; 
    info_box.initialize_gauss_quadrature(n_gauss_points, 
                                         n_gauss_points, 
                                         n_gauss_points); 


    info_box.initialize_update_flags(); 
    UpdateFlags update_flags = 
      update_quadrature_points | update_values | update_gradients; 
    info_box.add_update_flags(update_flags, true, true, true, true); 


    info_box.initialize(fe, mapping); 


    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 


    MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double>> 
      assembler; 
    assembler.initialize(system_matrix, right_hand_side); 



    MeshWorker::loop<dim, 
                     dim, 
                     MeshWorker::DoFInfo<dim>, 
                     MeshWorker::IntegrationInfoBox<dim>>( 
      dof_handler.begin_active(), 
      dof_handler.end(), 
      dof_info, 
      info_box, 
      &AdvectionProblem<dim>::integrate_cell_term, 
      &AdvectionProblem<dim>::integrate_boundary_term, 
      &AdvectionProblem<dim>::integrate_face_term, 
      assembler); 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::integrate_cell_term(DoFInfo & dinfo, 
                                                  CellInfo &info) 
  { 


    const FEValuesBase<dim> &  fe_values    = info.fe_values(); 
    FullMatrix<double> &       local_matrix = dinfo.matrix(0).matrix; 
    const std::vector<double> &JxW          = fe_values.get_JxW_values(); 


    for (unsigned int point = 0; point < fe_values.n_quadrature_points; ++point) 
      { 
        const Tensor<1, dim> beta_at_q_point = 
          beta(fe_values.quadrature_point(point)); 


        for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < fe_values.dofs_per_cell; ++j) 
            local_matrix(i, j) += -beta_at_q_point *                // 
                                  fe_values.shape_grad(i, point) *  // 
                                  fe_values.shape_value(j, point) * // 
                                  JxW[point]; 
      } 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::integrate_boundary_term(DoFInfo & dinfo, 
                                                      CellInfo &info) 
  { 
    const FEValuesBase<dim> &fe_face_values = info.fe_values(); 
    FullMatrix<double> &     local_matrix   = dinfo.matrix(0).matrix; 
    Vector<double> &         local_vector   = dinfo.vector(0).block(0); 

    const std::vector<double> &        JxW = fe_face_values.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = 
      fe_face_values.get_normal_vectors(); 

    std::vector<double> g(fe_face_values.n_quadrature_points); 

    static BoundaryValues<dim> boundary_function; 
    boundary_function.value_list(fe_face_values.get_quadrature_points(), g); 

    for (unsigned int point = 0; point < fe_face_values.n_quadrature_points; 
         ++point) 
      { 
        const double beta_dot_n = 
          beta(fe_face_values.quadrature_point(point)) * normals[point]; 
        if (beta_dot_n > 0) 
          for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < fe_face_values.dofs_per_cell; ++j) 
              local_matrix(i, j) += beta_dot_n *                           // 
                                    fe_face_values.shape_value(j, point) * // 
                                    fe_face_values.shape_value(i, point) * // 
                                    JxW[point]; 
        else 
          for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i) 
            local_vector(i) += -beta_dot_n *                          // 
                               g[point] *                             // 
                               fe_face_values.shape_value(i, point) * // 
                               JxW[point]; 
      } 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::integrate_face_term(DoFInfo & dinfo1, 
                                                  DoFInfo & dinfo2, 
                                                  CellInfo &info1, 
                                                  CellInfo &info2) 
  { 


    const FEValuesBase<dim> &fe_face_values = info1.fe_values(); 
    const unsigned int       dofs_per_cell  = fe_face_values.dofs_per_cell; 


    const FEValuesBase<dim> &fe_face_values_neighbor = info2.fe_values(); 
    const unsigned int       neighbor_dofs_per_cell = 
      fe_face_values_neighbor.dofs_per_cell; 


    FullMatrix<double> &u1_v1_matrix = dinfo1.matrix(0, false).matrix; 
    FullMatrix<double> &u2_v1_matrix = dinfo1.matrix(0, true).matrix; 
    FullMatrix<double> &u1_v2_matrix = dinfo2.matrix(0, true).matrix; 
    FullMatrix<double> &u2_v2_matrix = dinfo2.matrix(0, false).matrix; 


    const std::vector<double> &        JxW = fe_face_values.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = 
      fe_face_values.get_normal_vectors(); 

    for (unsigned int point = 0; point < fe_face_values.n_quadrature_points; 
         ++point) 
      { 
        const double beta_dot_n = 
          beta(fe_face_values.quadrature_point(point)) * normals[point]; 
        if (beta_dot_n > 0) 
          { 


            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                u1_v1_matrix(i, j) += beta_dot_n *                           // 
                                      fe_face_values.shape_value(j, point) * // 
                                      fe_face_values.shape_value(i, point) * // 
                                      JxW[point]; 


            for (unsigned int k = 0; k < neighbor_dofs_per_cell; ++k) 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                u1_v2_matrix(k, j) += 
                  -beta_dot_n *                                   // 
                  fe_face_values.shape_value(j, point) *          // 
                  fe_face_values_neighbor.shape_value(k, point) * // 
                  JxW[point]; 
          } 
        else 
          { 


            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              for (unsigned int l = 0; l < neighbor_dofs_per_cell; ++l) 
                u2_v1_matrix(i, l) += 
                  beta_dot_n *                                    // 
                  fe_face_values_neighbor.shape_value(l, point) * // 
                  fe_face_values.shape_value(i, point) *          // 
                  JxW[point]; 


            for (unsigned int k = 0; k < neighbor_dofs_per_cell; ++k) 
              for (unsigned int l = 0; l < neighbor_dofs_per_cell; ++l) 
                u2_v2_matrix(k, l) += 
                  -beta_dot_n *                                   // 
                  fe_face_values_neighbor.shape_value(l, point) * // 
                  fe_face_values_neighbor.shape_value(k, point) * // 
                  JxW[point]; 
          } 
      } 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::solve(Vector<double> &solution) 
  { 
    SolverControl                    solver_control(1000, 1e-12); 
    SolverRichardson<Vector<double>> solver(solver_control); 


    PreconditionBlockSSOR<SparseMatrix<double>> preconditioner; 


    preconditioner.initialize(system_matrix, fe.n_dofs_per_cell()); 


    solver.solve(system_matrix, solution, right_hand_side, preconditioner); 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::refine_grid() 
  { 


    Vector<float> gradient_indicator(triangulation.n_active_cells()); 


    DerivativeApproximation::approximate_gradient(mapping, 
                                                  dof_handler, 
                                                  solution, 
                                                  gradient_indicator); 

    unsigned int cell_no = 0; 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      gradient_indicator(cell_no++) *= 
        std::pow(cell->diameter(), 1 + 1.0 * dim / 2); 


    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    gradient_indicator, 
                                                    0.3, 
                                                    0.1); 

    triangulation.execute_coarsening_and_refinement(); 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::output_results(const unsigned int cycle) const 
  { 


    { 
      const std::string filename = "grid-" + std::to_string(cycle) + ".eps"; 
      deallog << "Writing grid to <" << filename << ">" << std::endl; 
      std::ofstream eps_output(filename); 

      GridOut grid_out; 
      grid_out.write_eps(triangulation, eps_output); 
    } 


    { 
      const std::string filename = "sol-" + std::to_string(cycle) + ".gnuplot"; 
      deallog << "Writing solution to <" << filename << ">" << std::endl; 
      std::ofstream gnuplot_output(filename); 

      DataOut<dim> data_out; 
      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution, "u"); 

      data_out.build_patches(); 

      data_out.write_gnuplot(gnuplot_output); 
    } 
  } 


  template <int dim> 
  void AdvectionProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      { 
        deallog << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_cube(triangulation); 

            triangulation.refine_global(3); 
          } 
        else 
          refine_grid(); 

        deallog << "Number of active cells:       " 
                << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        deallog << "Number of degrees of freedom: " << dof_handler.n_dofs() 
                << std::endl; 

        assemble_system(); 
        solve(solution); 

        output_results(cycle); 
      } 
  } 
} // namespace Step12 


int main() 
{ 
  try 
    { 
      dealii::deallog.depth_console(5); 

      Step12::AdvectionProblem<2> dgmethod; 
      dgmethod.run(); 
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


