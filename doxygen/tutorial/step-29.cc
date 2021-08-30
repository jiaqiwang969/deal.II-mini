

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2007 - 2021 by the deal.II authors 
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
 * Author: Moritz Allmaras, Texas A&M University, 2007 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 

#include <iostream> 
#include <fstream> 


#include <deal.II/base/parameter_handler.h> 


#include <deal.II/lac/sparse_direct.h> 


#include <deal.II/fe/fe_system.h> 


#include <deal.II/base/timer.h> 


namespace Step29 
{ 
  using namespace dealii; 



  template <int dim> 
  class DirichletBoundaryValues : public Function<dim> 
  { 
  public: 
    DirichletBoundaryValues() 
      : Function<dim>(2) 
    {} 

    virtual void vector_value(const Point<dim> & /*p*/, 
                              Vector<double> &values) const override 
    { 
      Assert(values.size() == 2, ExcDimensionMismatch(values.size(), 2)); 

      values(0) = 1; 
      values(1) = 0; 
    } 

    virtual void 
    vector_value_list(const std::vector<Point<dim>> &points, 
                      std::vector<Vector<double>> &  value_list) const override 
    { 
      Assert(value_list.size() == points.size(), 
             ExcDimensionMismatch(value_list.size(), points.size())); 

      for (unsigned int p = 0; p < points.size(); ++p) 
        DirichletBoundaryValues<dim>::vector_value(points[p], value_list[p]); 
    } 
  }; 


  class ParameterReader : public Subscriptor 
  { 
  public: 
    ParameterReader(ParameterHandler &); 
    void read_parameters(const std::string &); 

  private: 
    void              declare_parameters(); 
    ParameterHandler &prm; 
  }; 


  ParameterReader::ParameterReader(ParameterHandler &paramhandler) 
    : prm(paramhandler) 
  {} 


  void ParameterReader::declare_parameters() 
  { 


    prm.enter_subsection("Mesh & geometry parameters"); 
    { 
      prm.declare_entry("Number of refinements", 
                        "6", 
                        Patterns::Integer(0), 
                        "Number of global mesh refinement steps " 
                        "applied to initial coarse grid"); 

      prm.declare_entry("Focal distance", 
                        "0.3", 
                        Patterns::Double(0), 
                        "Distance of the focal point of the lens " 
                        "to the x-axis"); 
    } 
    prm.leave_subsection(); 


    prm.enter_subsection("Physical constants"); 
    { 
      prm.declare_entry("c", "1.5e5", Patterns::Double(0), "Wave speed"); 

      prm.declare_entry("omega", "5.0e7", Patterns::Double(0), "Frequency"); 
    } 
    prm.leave_subsection(); 


    prm.enter_subsection("Output parameters"); 
    { 
      prm.declare_entry("Output filename", 
                        "solution", 
                        Patterns::Anything(), 
                        "Name of the output file (without extension)"); 


      DataOutInterface<1>::declare_parameters(prm); 
    } 
    prm.leave_subsection(); 
  } 


  void ParameterReader::read_parameters(const std::string &parameter_file) 
  { 
    declare_parameters(); 

    prm.parse_input(parameter_file); 
  } 





  template <int dim> 
  class ComputeIntensity : public DataPostprocessorScalar<dim> 
  { 
  public: 
    ComputeIntensity(); 

    virtual void evaluate_vector_field( 
      const DataPostprocessorInputs::Vector<dim> &inputs, 
      std::vector<Vector<double>> &computed_quantities) const override; 
  }; 



  template <int dim> 
  ComputeIntensity<dim>::ComputeIntensity() 
    : DataPostprocessorScalar<dim>("Intensity", update_values) 
  {} 


  template <int dim> 
  void ComputeIntensity<dim>::evaluate_vector_field( 
    const DataPostprocessorInputs::Vector<dim> &inputs, 
    std::vector<Vector<double>> &               computed_quantities) const 
  { 
    Assert(computed_quantities.size() == inputs.solution_values.size(), 
           ExcDimensionMismatch(computed_quantities.size(), 
                                inputs.solution_values.size())); 


    for (unsigned int i = 0; i < computed_quantities.size(); i++) 
      { 
        Assert(computed_quantities[i].size() == 1, 
               ExcDimensionMismatch(computed_quantities[i].size(), 1)); 
        Assert(inputs.solution_values[i].size() == 2, 
               ExcDimensionMismatch(inputs.solution_values[i].size(), 2)); 

        const std::complex<double> u(inputs.solution_values[i](0), 
                                     inputs.solution_values[i](1)); 

        computed_quantities[i](0) = std::abs(u); 
      } 
  } 


  template <int dim> 
  class UltrasoundProblem 
  { 
  public: 
    UltrasoundProblem(ParameterHandler &); 
    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void output_results() const; 

    ParameterHandler &prm; 

    Triangulation<dim> triangulation; 
    DoFHandler<dim>    dof_handler; 
    FESystem<dim>      fe; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 
    Vector<double>       solution, system_rhs; 
  }; 


  template <int dim> 
  UltrasoundProblem<dim>::UltrasoundProblem(ParameterHandler &param) 
    : prm(param) 
    , dof_handler(triangulation) 
    , fe(FE_Q<dim>(1), 2) 
  {} 


  template <int dim> 
  void UltrasoundProblem<dim>::make_grid() 
  { 


    std::cout << "Generating grid... "; 
    Timer timer; 


    prm.enter_subsection("Mesh & geometry parameters"); 

    const double       focal_distance = prm.get_double("Focal distance"); 
    const unsigned int n_refinements = prm.get_integer("Number of refinements"); 

    prm.leave_subsection(); 


    const Point<dim> transducer = 
      (dim == 2) ? Point<dim>(0.5, 0.0) : Point<dim>(0.5, 0.5, 0.0); 
    const Point<dim> focal_point = (dim == 2) ? 
                                     Point<dim>(0.5, focal_distance) : 
                                     Point<dim>(0.5, 0.5, focal_distance); 


    GridGenerator::subdivided_hyper_cube(triangulation, 5, 0, 1); 

    for (auto &cell : triangulation.cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary() && 
            ((face->center() - transducer).norm_square() < 0.01)) 
          { 
            face->set_boundary_id(1); 
            face->set_manifold_id(1); 
          } 


    triangulation.set_manifold(1, SphericalManifold<dim>(focal_point)); 


    triangulation.refine_global(n_refinements); 


    timer.stop(); 
    std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 

    std::cout << "  Number of active cells:  " << triangulation.n_active_cells() 
              << std::endl; 
  } 


  template <int dim> 
  void UltrasoundProblem<dim>::setup_system() 
  { 
    std::cout << "Setting up system... "; 
    Timer timer; 

    dof_handler.distribute_dofs(fe); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
    system_rhs.reinit(dof_handler.n_dofs()); 
    solution.reinit(dof_handler.n_dofs()); 

    timer.stop(); 
    std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 

    std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 
  } 


  template <int dim> 
  void UltrasoundProblem<dim>::assemble_system() 
  { 
    std::cout << "Assembling system matrix... "; 
    Timer timer; 


    prm.enter_subsection("Physical constants"); 

    const double omega = prm.get_double("omega"), c = prm.get_double("c"); 

    prm.leave_subsection(); 


    QGauss<dim>     quadrature_formula(fe.degree + 1); 
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 

    const unsigned int n_q_points      = quadrature_formula.size(), 
                       n_face_q_points = face_quadrature_formula.size(), 
                       dofs_per_cell   = fe.n_dofs_per_cell(); 


    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_JxW_values); 

    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_JxW_values); 


    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 


        cell_matrix = 0; 
        fe_values.reinit(cell); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          { 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              { 


                if (fe.system_to_component_index(i).first == 
                    fe.system_to_component_index(j).first) 
                  { 


                    for (unsigned int q_point = 0; q_point < n_q_points; 
                         ++q_point) 
                      cell_matrix(i, j) += 
                        (((fe_values.shape_value(i, q_point) * 
                           fe_values.shape_value(j, q_point)) * 
                            (-omega * omega) + 
                          (fe_values.shape_grad(i, q_point) * 
                           fe_values.shape_grad(j, q_point)) * 
                            c * c) * 
                         fe_values.JxW(q_point)); 


                  } 
              } 
          } 


        for (const auto face_no : cell->face_indices()) 
          if (cell->face(face_no)->at_boundary() && 
              (cell->face(face_no)->boundary_id() == 0)) 
            { 


              fe_face_values.reinit(cell, face_no); 


              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                  if ((fe.system_to_component_index(i).first != 
                       fe.system_to_component_index(j).first) && 
                      fe.has_support_on_face(i, face_no) && 
                      fe.has_support_on_face(j, face_no)) 



                    for (unsigned int q_point = 0; q_point < n_face_q_points; 
                         ++q_point) 
                      cell_matrix(i, j) += 
                        ((fe.system_to_component_index(i).first == 0) ? -1 : 
                                                                        1) * 
                        fe_face_values.shape_value(i, q_point) * 
                        fe_face_values.shape_value(j, q_point) * c * omega * 
                        fe_face_values.JxW(q_point); 
            } 


        cell->get_dof_indices(local_dof_indices); 


        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(local_dof_indices[i], 
                              local_dof_indices[j], 
                              cell_matrix(i, j)); 
      } 


    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             1, 
                                             DirichletBoundaryValues<dim>(), 
                                             boundary_values); 

    MatrixTools::apply_boundary_values(boundary_values, 
                                       system_matrix, 
                                       solution, 
                                       system_rhs); 

    timer.stop(); 
    std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
  } 




  template <int dim> 
  void UltrasoundProblem<dim>::solve() 
  { 
    std::cout << "Solving linear system... "; 
    Timer timer; 


    SparseDirectUMFPACK A_direct; 
    A_direct.initialize(system_matrix); 


    A_direct.vmult(solution, system_rhs); 

    timer.stop(); 
    std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
  } 



  template <int dim> 
  void UltrasoundProblem<dim>::output_results() const 
  { 
    std::cout << "Generating output... "; 
    Timer timer; 


    ComputeIntensity<dim> intensities; 
    DataOut<dim>          data_out; 

    data_out.attach_dof_handler(dof_handler); 


    prm.enter_subsection("Output parameters"); 

    const std::string output_filename = prm.get("Output filename"); 
    data_out.parse_parameters(prm); 

    prm.leave_subsection(); 


    const std::string filename = output_filename + data_out.default_suffix(); 

    std::ofstream output(filename); 


    std::vector<std::string> solution_names; 
    solution_names.emplace_back("Re_u"); 
    solution_names.emplace_back("Im_u"); 

    data_out.add_data_vector(solution, solution_names); 


    data_out.add_data_vector(solution, intensities); 


    data_out.build_patches(); 
    data_out.write(output); 

    timer.stop(); 
    std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
  } 



  template <int dim> 
  void UltrasoundProblem<dim>::run() 
  { 
    make_grid(); 
    setup_system(); 
    assemble_system(); 
    solve(); 
    output_results(); 
  } 
} // namespace Step29 


int main() 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step29; 

      ParameterHandler prm; 
      ParameterReader  param(prm); 
      param.read_parameters("step-29.prm"); 

      UltrasoundProblem<2> ultrasound_problem(prm); 
      ultrasound_problem.run(); 
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


