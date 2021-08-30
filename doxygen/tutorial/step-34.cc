

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
 * Author: Luca Heltai, Cataldo Manigrasso, 2009 
 */ 




#include <deal.II/base/smartpointer.h> 
#include <deal.II/base/convergence_table.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/quadrature_selector.h> 
#include <deal.II/base/parsed_function.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/solver_control.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/precondition.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_in.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 


#include <cmath> 
#include <iostream> 
#include <fstream> 
#include <string> 


namespace Step34 
{ 
  using namespace dealii; 



  namespace LaplaceKernel 
  { 
    template <int dim> 
    double single_layer(const Tensor<1, dim> &R) 
    { 
      switch (dim) 
        { 
          case 2: 
            return (-std::log(R.norm()) / (2 * numbers::PI)); 

          case 3: 
            return (1. / (R.norm() * 4 * numbers::PI)); 

          default: 
            Assert(false, ExcInternalError()); 
            return 0.; 
        } 
    } 

    template <int dim> 
    Tensor<1, dim> double_layer(const Tensor<1, dim> &R) 
    { 
      switch (dim) 
        { 
          case 2: 
            return R / (-2 * numbers::PI * R.norm_square()); 
          case 3: 
            return R / (-4 * numbers::PI * R.norm_square() * R.norm()); 

          default: 
            Assert(false, ExcInternalError()); 
            return Tensor<1, dim>(); 
        } 
    } 
  } // namespace LaplaceKernel 


  template <int dim> 
  class BEMProblem 
  { 
  public: 
    BEMProblem(const unsigned int fe_degree      = 1, 
               const unsigned int mapping_degree = 1); 

    void run(); 

  private: 
    void read_parameters(const std::string &filename); 

    void read_domain(); 

    void refine_and_resize(); 




    void assemble_system(); 



    void solve_system(); 



    void compute_errors(const unsigned int cycle); 




    void compute_exterior_solution(); 

    void output_results(const unsigned int cycle); 


    const Quadrature<dim - 1> &get_singular_quadrature( 
      const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, 
      const unsigned int index) const; 




    Triangulation<dim - 1, dim> tria; 
    FE_Q<dim - 1, dim>          fe; 
    DoFHandler<dim - 1, dim>    dof_handler; 
    MappingQ<dim - 1, dim>      mapping; 


    FullMatrix<double> system_matrix; 
    Vector<double>     system_rhs; 


    Vector<double> phi; 
    Vector<double> alpha; 


    ConvergenceTable convergence_table; 





    Functions::ParsedFunction<dim> wind; 
    Functions::ParsedFunction<dim> exact_solution; 

    unsigned int                         singular_quadrature_order; 
    std::shared_ptr<Quadrature<dim - 1>> quadrature; 

    SolverControl solver_control; 

    unsigned int n_cycles; 
    unsigned int external_refinement; 

    bool run_in_this_dimension; 
    bool extend_solution; 
  }; 



  template <int dim> 
  BEMProblem<dim>::BEMProblem(const unsigned int fe_degree, 
                              const unsigned int mapping_degree) 
    : fe(fe_degree) 
    , dof_handler(tria) 
    , mapping(mapping_degree, true) 
    , wind(dim) 
    , singular_quadrature_order(5) 
    , n_cycles(4) 
    , external_refinement(5) 
    , run_in_this_dimension(true) 
    , extend_solution(true) 
  {} 

  template <int dim> 
  void BEMProblem<dim>::read_parameters(const std::string &filename) 
  { 
    deallog << std::endl 
            << "Parsing parameter file " << filename << std::endl 
            << "for a " << dim << " dimensional simulation. " << std::endl; 

    ParameterHandler prm; 

    prm.declare_entry("Number of cycles", "4", Patterns::Integer()); 
    prm.declare_entry("External refinement", "5", Patterns::Integer()); 
    prm.declare_entry("Extend solution on the -2,2 box", 
                      "true", 
                      Patterns::Bool()); 
    prm.declare_entry("Run 2d simulation", "true", Patterns::Bool()); 
    prm.declare_entry("Run 3d simulation", "true", Patterns::Bool()); 

    prm.enter_subsection("Quadrature rules"); 
    { 
      prm.declare_entry( 
        "Quadrature type", 
        "gauss", 
        Patterns::Selection( 
          QuadratureSelector<(dim - 1)>::get_quadrature_names())); 
      prm.declare_entry("Quadrature order", "4", Patterns::Integer()); 
      prm.declare_entry("Singular quadrature order", "5", Patterns::Integer()); 
    } 
    prm.leave_subsection(); 





    prm.enter_subsection("Wind function 2d"); 
    { 
      Functions::ParsedFunction<2>::declare_parameters(prm, 2); 
      prm.set("Function expression", "1; 1"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Wind function 3d"); 
    { 
      Functions::ParsedFunction<3>::declare_parameters(prm, 3); 
      prm.set("Function expression", "1; 1; 1"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Exact solution 2d"); 
    { 
      Functions::ParsedFunction<2>::declare_parameters(prm); 
      prm.set("Function expression", "x+y"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Exact solution 3d"); 
    { 
      Functions::ParsedFunction<3>::declare_parameters(prm); 
      prm.set("Function expression", "x+y+z"); 
    } 
    prm.leave_subsection(); 


    prm.enter_subsection("Solver"); 
    SolverControl::declare_parameters(prm); 
    prm.leave_subsection(); 


    prm.parse_input(filename); 

    n_cycles            = prm.get_integer("Number of cycles"); 
    external_refinement = prm.get_integer("External refinement"); 
    extend_solution     = prm.get_bool("Extend solution on the -2,2 box"); 

    prm.enter_subsection("Quadrature rules"); 
    { 
      quadrature = std::shared_ptr<Quadrature<dim - 1>>( 
        new QuadratureSelector<dim - 1>(prm.get("Quadrature type"), 
                                        prm.get_integer("Quadrature order"))); 
      singular_quadrature_order = prm.get_integer("Singular quadrature order"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Wind function " + std::to_string(dim) + "d"); 
    { 
      wind.parse_parameters(prm); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Exact solution " + std::to_string(dim) + "d"); 
    { 
      exact_solution.parse_parameters(prm); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Solver"); 
    solver_control.parse_parameters(prm); 
    prm.leave_subsection(); 


    run_in_this_dimension = 
      prm.get_bool("Run " + std::to_string(dim) + "d simulation"); 
  } 




  template <int dim> 
  void BEMProblem<dim>::read_domain() 
  { 
    const Point<dim>                      center = Point<dim>(); 
    const SphericalManifold<dim - 1, dim> manifold(center); 

    std::ifstream in; 
    switch (dim) 
      { 
        case 2: 
          in.open("coarse_circle.inp"); 
          break; 

        case 3: 
          in.open("coarse_sphere.inp"); 
          break; 

        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    GridIn<dim - 1, dim> gi; 
    gi.attach_triangulation(tria); 
    gi.read_ucd(in); 

    tria.set_all_manifold_ids(1); 


    tria.set_manifold(1, manifold); 
  } 


  template <int dim> 
  void BEMProblem<dim>::refine_and_resize() 
  { 
    tria.refine_global(1); 

    dof_handler.distribute_dofs(fe); 

    const unsigned int n_dofs = dof_handler.n_dofs(); 

    system_matrix.reinit(n_dofs, n_dofs); 

    system_rhs.reinit(n_dofs); 
    phi.reinit(n_dofs); 
    alpha.reinit(n_dofs); 
  } 


  template <int dim> 
  void BEMProblem<dim>::assemble_system() 
  { 


    FEValues<dim - 1, dim> fe_v(mapping, 
                                fe, 
                                *quadrature, 
                                update_values | update_normal_vectors | 
                                  update_quadrature_points | update_JxW_values); 

    const unsigned int n_q_points = fe_v.n_quadrature_points; 

    std::vector<types::global_dof_index> local_dof_indices( 
      fe.n_dofs_per_cell()); 

    std::vector<Vector<double>> cell_wind(n_q_points, Vector<double>(dim)); 
    double                      normal_wind; 


    Vector<double> local_matrix_row_i(fe.n_dofs_per_cell()); 



    std::vector<Point<dim>> support_points(dof_handler.n_dofs()); 
    DoFTools::map_dofs_to_support_points<dim - 1, dim>(mapping, 
                                                       dof_handler, 
                                                       support_points); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_v.reinit(cell); 
        cell->get_dof_indices(local_dof_indices); 

        const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points(); 
        const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 
        wind.vector_value_list(q_points, cell_wind); 


        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
          { 
            local_matrix_row_i = 0; 

            bool         is_singular    = false; 
            unsigned int singular_index = numbers::invalid_unsigned_int; 

            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
              if (local_dof_indices[j] == i) 
                { 
                  singular_index = j; 
                  is_singular    = true; 
                  break; 
                } 


            if (is_singular == false) 
              { 
                for (unsigned int q = 0; q < n_q_points; ++q) 
                  { 
                    normal_wind = 0; 
                    for (unsigned int d = 0; d < dim; ++d) 
                      normal_wind += normals[q][d] * cell_wind[q](d); 

                    const Tensor<1, dim> R = q_points[q] - support_points[i]; 

                    system_rhs(i) += (LaplaceKernel::single_layer(R) * 
                                      normal_wind * fe_v.JxW(q)); 

                    for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 

                      local_matrix_row_i(j) -= 
                        ((LaplaceKernel::double_layer(R) * normals[q]) * 
                         fe_v.shape_value(j, q) * fe_v.JxW(q)); 
                  } 
              } 
            else 
              { 


                Assert(singular_index != numbers::invalid_unsigned_int, 
                       ExcInternalError()); 

                const Quadrature<dim - 1> &singular_quadrature = 
                  get_singular_quadrature(cell, singular_index); 

                FEValues<dim - 1, dim> fe_v_singular( 
                  mapping, 
                  fe, 
                  singular_quadrature, 
                  update_jacobians | update_values | update_normal_vectors | 
                    update_quadrature_points); 

                fe_v_singular.reinit(cell); 

                std::vector<Vector<double>> singular_cell_wind( 
                  singular_quadrature.size(), Vector<double>(dim)); 

                const std::vector<Tensor<1, dim>> &singular_normals = 
                  fe_v_singular.get_normal_vectors(); 
                const std::vector<Point<dim>> &singular_q_points = 
                  fe_v_singular.get_quadrature_points(); 

                wind.vector_value_list(singular_q_points, singular_cell_wind); 

                for (unsigned int q = 0; q < singular_quadrature.size(); ++q) 
                  { 
                    const Tensor<1, dim> R = 
                      singular_q_points[q] - support_points[i]; 
                    double normal_wind = 0; 
                    for (unsigned int d = 0; d < dim; ++d) 
                      normal_wind += 
                        (singular_cell_wind[q](d) * singular_normals[q][d]); 

                    system_rhs(i) += (LaplaceKernel::single_layer(R) * 
                                      normal_wind * fe_v_singular.JxW(q)); 

                    for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
                      { 
                        local_matrix_row_i(j) -= 
                          ((LaplaceKernel::double_layer(R) * 
                            singular_normals[q]) * 
                           fe_v_singular.shape_value(j, q) * 
                           fe_v_singular.JxW(q)); 
                      } 
                  } 
              } 


            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
              system_matrix(i, local_dof_indices[j]) += local_matrix_row_i(j); 
          } 
      } 



    Vector<double> ones(dof_handler.n_dofs()); 
    ones.add(-1.); 

    system_matrix.vmult(alpha, ones); 
    alpha.add(1); 
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
      system_matrix(i, i) += alpha(i); 
  } 


  template <int dim> 
  void BEMProblem<dim>::solve_system() 
  { 
    SolverGMRES<Vector<double>> solver(solver_control); 
    solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity()); 
  } 


  template <int dim> 
  void BEMProblem<dim>::compute_errors(const unsigned int cycle) 
  { 
    Vector<float> difference_per_cell(tria.n_active_cells()); 
    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      phi, 
                                      exact_solution, 
                                      difference_per_cell, 
                                      QGauss<(dim - 1)>(2 * fe.degree + 1), 
                                      VectorTools::L2_norm); 
    const double L2_error = 
      VectorTools::compute_global_error(tria, 
                                        difference_per_cell, 
                                        VectorTools::L2_norm); 


    Vector<double> difference_per_node(alpha); 
    difference_per_node.add(-.5); 

    const double       alpha_error    = difference_per_node.linfty_norm(); 
    const unsigned int n_active_cells = tria.n_active_cells(); 
    const unsigned int n_dofs         = dof_handler.n_dofs(); 

    deallog << "Cycle " << cycle << ':' << std::endl 
            << "   Number of active cells:       " << n_active_cells 
            << std::endl 
            << "   Number of degrees of freedom: " << n_dofs << std::endl; 

    convergence_table.add_value("cycle", cycle); 
    convergence_table.add_value("cells", n_active_cells); 
    convergence_table.add_value("dofs", n_dofs); 
    convergence_table.add_value("L2(phi)", L2_error); 
    convergence_table.add_value("Linfty(alpha)", alpha_error); 
  } 















  template <> 
  const Quadrature<2> &BEMProblem<3>::get_singular_quadrature( 
    const DoFHandler<2, 3>::active_cell_iterator &, 
    const unsigned int index) const 
  { 
    Assert(index < fe.n_dofs_per_cell(), 
           ExcIndexRange(0, fe.n_dofs_per_cell(), index)); 

    static std::vector<QGaussOneOverR<2>> quadratures; 
    if (quadratures.size() == 0) 
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
        quadratures.emplace_back(singular_quadrature_order, 
                                 fe.get_unit_support_points()[i], 
                                 true); 
    return quadratures[index]; 
  } 

  template <> 
  const Quadrature<1> &BEMProblem<2>::get_singular_quadrature( 
    const DoFHandler<1, 2>::active_cell_iterator &cell, 
    const unsigned int                            index) const 
  { 
    Assert(index < fe.n_dofs_per_cell(), 
           ExcIndexRange(0, fe.n_dofs_per_cell(), index)); 

    static Quadrature<1> *q_pointer = nullptr; 
    if (q_pointer) 
      delete q_pointer; 

    q_pointer = new QGaussLogR<1>(singular_quadrature_order, 
                                  fe.get_unit_support_points()[index], 
                                  1. / cell->measure(), 
                                  true); 
    return (*q_pointer); 
  } 





  template <int dim> 
  void BEMProblem<dim>::compute_exterior_solution() 
  { 
    Triangulation<dim> external_tria; 
    GridGenerator::hyper_cube(external_tria, -2, 2); 

    FE_Q<dim>       external_fe(1); 
    DoFHandler<dim> external_dh(external_tria); 
    Vector<double>  external_phi; 

    external_tria.refine_global(external_refinement); 
    external_dh.distribute_dofs(external_fe); 
    external_phi.reinit(external_dh.n_dofs()); 

    FEValues<dim - 1, dim> fe_v(mapping, 
                                fe, 
                                *quadrature, 
                                update_values | update_normal_vectors | 
                                  update_quadrature_points | update_JxW_values); 

    const unsigned int n_q_points = fe_v.n_quadrature_points; 

    std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell()); 

    std::vector<double>         local_phi(n_q_points); 
    std::vector<double>         normal_wind(n_q_points); 
    std::vector<Vector<double>> local_wind(n_q_points, Vector<double>(dim)); 

    std::vector<Point<dim>> external_support_points(external_dh.n_dofs()); 
    DoFTools::map_dofs_to_support_points<dim>(StaticMappingQ1<dim>::mapping, 
                                              external_dh, 
                                              external_support_points); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_v.reinit(cell); 

        const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points(); 
        const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 

        cell->get_dof_indices(dofs); 
        fe_v.get_function_values(phi, local_phi); 

        wind.vector_value_list(q_points, local_wind); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            normal_wind[q] = 0; 
            for (unsigned int d = 0; d < dim; ++d) 
              normal_wind[q] += normals[q][d] * local_wind[q](d); 
          } 

        for (unsigned int i = 0; i < external_dh.n_dofs(); ++i) 
          for (unsigned int q = 0; q < n_q_points; ++q) 
            { 
              const Tensor<1, dim> R = q_points[q] - external_support_points[i]; 

              external_phi(i) += 
                ((LaplaceKernel::single_layer(R) * normal_wind[q] + 
                  (LaplaceKernel::double_layer(R) * normals[q]) * 
                    local_phi[q]) * 
                 fe_v.JxW(q)); 
            } 
      } 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(external_dh); 
    data_out.add_data_vector(external_phi, "external_phi"); 
    data_out.build_patches(); 

    const std::string filename = std::to_string(dim) + "d_external.vtk"; 
    std::ofstream     file(filename); 

 
  } 


  template <int dim> 
  void BEMProblem<dim>::output_results(const unsigned int cycle) 
  { 
    DataOut<dim - 1, dim> dataout; 

    dataout.attach_dof_handler(dof_handler); 
    dataout.add_data_vector(phi, "phi", DataOut<dim - 1, dim>::type_dof_data); 
    dataout.add_data_vector(alpha, 
                            "alpha", 
                            DataOut<dim - 1, dim>::type_dof_data); 
    dataout.build_patches(mapping, 
                          mapping.get_degree(), 
                          DataOut<dim - 1, dim>::curved_inner_cells); 

    const std::string filename = std::to_string(dim) + "d_boundary_solution_" + 
                                 std::to_string(cycle) + ".vtk"; 
    std::ofstream file(filename); 

    dataout.write_vtk(file); 

    if (cycle == n_cycles - 1) 
      { 
        convergence_table.set_precision("L2(phi)", 3); 
        convergence_table.set_precision("Linfty(alpha)", 3); 

        convergence_table.set_scientific("L2(phi)", true); 
        convergence_table.set_scientific("Linfty(alpha)", true); 

        convergence_table.evaluate_convergence_rates(
          "L2(phi)", ConvergenceTable::reduction_rate_log2); 
        convergence_table.evaluate_convergence_rates( 
          "Linfty(alpha)", ConvergenceTable::reduction_rate_log2); 
        deallog << std::endl; 
        convergence_table.write_text(std::cout); 
      } 
  } 


  template <int dim> 
  void BEMProblem<dim>::run() 
  { 
    read_parameters("parameters.prm"); 

    if (run_in_this_dimension == false) 
      { 
        deallog << "Run in dimension " << dim 
                << " explicitly disabled in parameter file. " << std::endl; 
        return; 
      } 

    read_domain(); 

    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
      { 
        refine_and_resize(); 
        assemble_system(); 
        solve_system(); 
        compute_errors(cycle); 
        output_results(cycle); 
      } 

    if (extend_solution == true) 
      compute_exterior_solution(); 
  } 
} // namespace Step34 


int main() 
{ 
  try 
    { 
      using namespace Step34; 

      const unsigned int degree         = 1; 
      const unsigned int mapping_degree = 1; 

      deallog.depth_console(3); 
      BEMProblem<2> laplace_problem_2d(degree, mapping_degree); 
      laplace_problem_2d.run(); 

      BEMProblem<3> laplace_problem_3d(degree, mapping_degree); 
      laplace_problem_3d.run(); 
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


