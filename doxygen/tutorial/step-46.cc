

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2011 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2011 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/sparse_direct.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_nothing.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/hp/fe_collection.h> 
#include <deal.II/hp/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <iostream> 
#include <fstream> 

namespace Step46 
{ 
  using namespace dealii; 




  template <int dim> 
  class FluidStructureProblem 
  { 
  public: 
    FluidStructureProblem(const unsigned int stokes_degree, 
                          const unsigned int elasticity_degree); 
    void run(); 

  private: 
    enum 
    { 
      fluid_domain_id, 
      solid_domain_id 
    }; 

    static bool cell_is_in_fluid_domain( 
      const typename DoFHandler<dim>::cell_iterator &cell); 

    static bool cell_is_in_solid_domain( 
      const typename DoFHandler<dim>::cell_iterator &cell); 

    void make_grid(); 
    void set_active_fe_indices(); 
    void setup_dofs(); 
    void assemble_system(); 
    void assemble_interface_term( 
      const FEFaceValuesBase<dim> &         elasticity_fe_face_values, 
      const FEFaceValuesBase<dim> &         stokes_fe_face_values, 
      std::vector<Tensor<1, dim>> &         elasticity_phi, 
      std::vector<SymmetricTensor<2, dim>> &stokes_symgrad_phi_u, 
      std::vector<double> &                 stokes_phi_p, 
      FullMatrix<double> &                  local_interface_matrix) const; 
    void solve(); 
    void output_results(const unsigned int refinement_cycle) const; 
    void refine_mesh(); 

    const unsigned int stokes_degree; 
    const unsigned int elasticity_degree; 

    Triangulation<dim>    triangulation; 
    FESystem<dim>         stokes_fe; 
    FESystem<dim>         elasticity_fe; 
    hp::FECollection<dim> fe_collection; 
    DoFHandler<dim>       dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    const double viscosity; 
    const double lambda; 
    const double mu; 
  }; 


  template <int dim> 
  class StokesBoundaryValues : public Function<dim> 
  { 
  public: 
    StokesBoundaryValues() 
      : Function<dim>(dim + 1 + dim) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  value) const override; 
  }; 

  template <int dim> 
  double StokesBoundaryValues<dim>::value(const Point<dim> & p, 
                                          const unsigned int component) const 
  { 
    Assert(component < this->n_components, 
           ExcIndexRange(component, 0, this->n_components)); 

    if (component == dim - 1) 
      switch (dim) 
        { 
          case 2: 
            return std::sin(numbers::PI * p[0]); 
          case 3: 
            return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]); 
          default: 
            Assert(false, ExcNotImplemented()); 
        } 

    return 0; 
  } 

  template <int dim> 
  void StokesBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                                               Vector<double> &  values) const 
  { 
    for (unsigned int c = 0; c < this->n_components; ++c) 
      values(c) = StokesBoundaryValues<dim>::value(p, c); 
  } 



  template <int dim> 
  FluidStructureProblem<dim>::FluidStructureProblem( 
    const unsigned int stokes_degree, 
    const unsigned int elasticity_degree) 
    : stokes_degree(stokes_degree) 
    , elasticity_degree(elasticity_degree) 
    , triangulation(Triangulation<dim>::maximum_smoothing) 
    , stokes_fe(FE_Q<dim>(stokes_degree + 1), 
                dim, 
                FE_Q<dim>(stokes_degree), 
                1, 
                FE_Nothing<dim>(), 
                dim) 
    , elasticity_fe(FE_Nothing<dim>(), 
                    dim, 
                    FE_Nothing<dim>(), 
                    1, 
                    FE_Q<dim>(elasticity_degree), 
                    dim) 
    , dof_handler(triangulation) 
    , viscosity(2) 
    , lambda(1) 
    , mu(1) 
  { 
    fe_collection.push_back(stokes_fe); 
    fe_collection.push_back(elasticity_fe); 
  } 

  template <int dim> 
  bool FluidStructureProblem<dim>::cell_is_in_fluid_domain( 
    const typename DoFHandler<dim>::cell_iterator &cell) 
  { 
    return (cell->material_id() == fluid_domain_id); 
  } 

  template <int dim> 
  bool FluidStructureProblem<dim>::cell_is_in_solid_domain( 
    const typename DoFHandler<dim>::cell_iterator &cell) 
  { 
    return (cell->material_id() == solid_domain_id); 
  } 


  template <int dim> 
  void FluidStructureProblem<dim>::make_grid() 
  { 
    GridGenerator::subdivided_hyper_cube(triangulation, 8, -1, 1); 

    for (const auto &cell : triangulation.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary() && (face->center()[dim - 1] == 1)) 
          face->set_all_boundary_ids(1); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (((std::fabs(cell->center()[0]) < 0.25) && 
           (cell->center()[dim - 1] > 0.5)) || 
          ((std::fabs(cell->center()[0]) >= 0.25) && 
           (cell->center()[dim - 1] > -0.5))) 
        cell->set_material_id(fluid_domain_id); 
      else 
        cell->set_material_id(solid_domain_id); 
  } 



  template <int dim> 
  void FluidStructureProblem<dim>::set_active_fe_indices() 
  { 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        if (cell_is_in_fluid_domain(cell)) 
          cell->set_active_fe_index(0); 
        else if (cell_is_in_solid_domain(cell)) 
          cell->set_active_fe_index(1); 
        else 
          Assert(false, ExcNotImplemented()); 
      } 
  } 


  template <int dim> 
  void FluidStructureProblem<dim>::setup_dofs() 
  { 
    set_active_fe_indices(); 
    dof_handler.distribute_dofs(fe_collection); 

    { 
      constraints.clear(); 
      DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

      const FEValuesExtractors::Vector velocities(0); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               1, 
                                               StokesBoundaryValues<dim>(), 
                                               constraints, 
                                               fe_collection.component_mask( 
                                                 velocities)); 

      const FEValuesExtractors::Vector displacements(dim + 1); 
      VectorTools::interpolate_boundary_values( 
        dof_handler, 
        0, 
        Functions::ZeroFunction<dim>(dim + 1 + dim), 
        constraints, 
        fe_collection.component_mask(displacements)); 
    } 


    { 
      std::vector<types::global_dof_index> local_face_dof_indices( 
        stokes_fe.n_dofs_per_face()); 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        if (cell_is_in_fluid_domain(cell)) 
          for (const auto face_no : cell->face_indices()) 
            if (cell->face(face_no)->at_boundary() == false) 
              { 
                bool face_is_on_interface = false; 

                if ((cell->neighbor(face_no)->has_children() == false) && 
                    (cell_is_in_solid_domain(cell->neighbor(face_no)))) 
                  face_is_on_interface = true; 
                else if (cell->neighbor(face_no)->has_children() == true) 
                  { 
                    for (unsigned int sf = 0; 
                         sf < cell->face(face_no)->n_children(); 
                         ++sf) 
                      if (cell_is_in_solid_domain( 
                            cell->neighbor_child_on_subface(face_no, sf))) 
                        { 
                          face_is_on_interface = true; 
                          break; 
                        } 
                  } 

                if (face_is_on_interface) 
                  { 
                    cell->face(face_no)->get_dof_indices(local_face_dof_indices, 
                                                         0); 
                    for (unsigned int i = 0; i < local_face_dof_indices.size(); 
                         ++i) 
                      if (stokes_fe.face_system_to_component_index(i).first < 
                          dim) 
                        constraints.add_line(local_face_dof_indices[i]); 
                  } 
              } 
    } 


    constraints.close(); 

    std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 


    { 
      DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 

      Table<2, DoFTools::Coupling> cell_coupling(fe_collection.n_components(), 
                                                 fe_collection.n_components()); 
      Table<2, DoFTools::Coupling> face_coupling(fe_collection.n_components(), 
                                                 fe_collection.n_components()); 

      for (unsigned int c = 0; c < fe_collection.n_components(); ++c) 
        for (unsigned int d = 0; d < fe_collection.n_components(); ++d) 
          { 
            if (((c < dim + 1) && (d < dim + 1) && 
                 !((c == dim) && (d == dim))) || 
                ((c >= dim + 1) && (d >= dim + 1))) 
              cell_coupling[c][d] = DoFTools::always; 

            if ((c >= dim + 1) && (d < dim + 1)) 
              face_coupling[c][d] = DoFTools::always; 
          } 

      DoFTools::make_flux_sparsity_pattern(dof_handler, 
                                           dsp, 
                                           cell_coupling, 
                                           face_coupling); 
      constraints.condense(dsp); 
      sparsity_pattern.copy_from(dsp); 
    } 

    system_matrix.reinit(sparsity_pattern); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 



  template <int dim> 
  void FluidStructureProblem<dim>::assemble_system() 
  { 
    system_matrix = 0; 
    system_rhs    = 0; 

    const QGauss<dim> stokes_quadrature(stokes_degree + 2); 
    const QGauss<dim> elasticity_quadrature(elasticity_degree + 2); 

    hp::QCollection<dim> q_collection; 
    q_collection.push_back(stokes_quadrature); 
    q_collection.push_back(elasticity_quadrature); 

    hp::FEValues<dim> hp_fe_values(fe_collection, 
                                   q_collection, 
                                   update_values | update_quadrature_points | 
                                     update_JxW_values | update_gradients); 

    const QGauss<dim - 1> common_face_quadrature( 
      std::max(stokes_degree + 2, elasticity_degree + 2)); 

    FEFaceValues<dim>    stokes_fe_face_values(stokes_fe, 
                                            common_face_quadrature, 
                                            update_JxW_values | 
                                              update_gradients | update_values); 
    FEFaceValues<dim>    elasticity_fe_face_values(elasticity_fe, 
                                                common_face_quadrature, 
                                                update_normal_vectors | 
                                                  update_values); 
    FESubfaceValues<dim> stokes_fe_subface_values(stokes_fe, 
                                                  common_face_quadrature, 
                                                  update_JxW_values | 
                                                    update_gradients | 
                                                    update_values); 
    FESubfaceValues<dim> elasticity_fe_subface_values(elasticity_fe, 
                                                      common_face_quadrature, 
                                                      update_normal_vectors | 
                                                        update_values); 


    const unsigned int stokes_dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
    const unsigned int elasticity_dofs_per_cell = 
      elasticity_fe.n_dofs_per_cell(); 

    FullMatrix<double> local_matrix; 
    FullMatrix<double> local_interface_matrix(elasticity_dofs_per_cell, 
                                              stokes_dofs_per_cell); 
    Vector<double>     local_rhs; 

    std::vector<types::global_dof_index> local_dof_indices; 
    std::vector<types::global_dof_index> neighbor_dof_indices( 
      stokes_dofs_per_cell); 

    const Functions::ZeroFunction<dim> right_hand_side(dim + 1); 


    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 
    const FEValuesExtractors::Vector displacements(dim + 1); 

    std::vector<SymmetricTensor<2, dim>> stokes_symgrad_phi_u( 
      stokes_dofs_per_cell); 
    std::vector<double> stokes_div_phi_u(stokes_dofs_per_cell); 
    std::vector<double> stokes_phi_p(stokes_dofs_per_cell); 

    std::vector<Tensor<2, dim>> elasticity_grad_phi(elasticity_dofs_per_cell); 
    std::vector<double>         elasticity_div_phi(elasticity_dofs_per_cell); 
    std::vector<Tensor<1, dim>> elasticity_phi(elasticity_dofs_per_cell); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        hp_fe_values.reinit(cell); 

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values(); 

        local_matrix.reinit(cell->get_fe().n_dofs_per_cell(), 
                            cell->get_fe().n_dofs_per_cell()); 
        local_rhs.reinit(cell->get_fe().n_dofs_per_cell()); 



        if (cell_is_in_fluid_domain(cell)) 
          { 
            const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 
            Assert(dofs_per_cell == stokes_dofs_per_cell, ExcInternalError()); 

            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) 
              { 
                for (unsigned int k = 0; k < dofs_per_cell; ++k) 
                  { 
                    stokes_symgrad_phi_u[k] = 
                      fe_values[velocities].symmetric_gradient(k, q); 
                    stokes_div_phi_u[k] = 
                      fe_values[velocities].divergence(k, q); 
                    stokes_phi_p[k] = fe_values[pressure].value(k, q); 
                  } 

                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    local_matrix(i, j) += 
                      (2 * viscosity * stokes_symgrad_phi_u[i] * 
                         stokes_symgrad_phi_u[j] - 
                       stokes_div_phi_u[i] * stokes_phi_p[j] - 
                       stokes_phi_p[i] * stokes_div_phi_u[j]) * 
                      fe_values.JxW(q); 
              } 
          } 
        else 
          { 
            const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 
            Assert(dofs_per_cell == elasticity_dofs_per_cell, 
                   ExcInternalError()); 

            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) 
              { 
                for (unsigned int k = 0; k < dofs_per_cell; ++k) 
                  { 
                    elasticity_grad_phi[k] = 
                      fe_values[displacements].gradient(k, q); 
                    elasticity_div_phi[k] = 
                      fe_values[displacements].divergence(k, q); 
                  } 

                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    { 
                      local_matrix(i, j) += 
                        (lambda * elasticity_div_phi[i] * 
                           elasticity_div_phi[j] + 
                         mu * scalar_product(elasticity_grad_phi[i], 
                                             elasticity_grad_phi[j]) + 
                         mu * 
                           scalar_product(elasticity_grad_phi[i], 
                                          transpose(elasticity_grad_phi[j]))) * 
                        fe_values.JxW(q); 
                    } 
              } 
          } 


        local_dof_indices.resize(cell->get_fe().n_dofs_per_cell()); 
        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global(local_matrix, 
                                               local_rhs, 
                                               local_dof_indices, 
                                               system_matrix, 
                                               system_rhs); 


        if (cell_is_in_solid_domain(cell)) 
          for (const auto f : cell->face_indices()) 
            if (cell->face(f)->at_boundary() == false) 
              { 





                if ((cell->neighbor(f)->level() == cell->level()) && 
                    (cell->neighbor(f)->has_children() == false) && 
                    cell_is_in_fluid_domain(cell->neighbor(f))) 
                  { 
                    elasticity_fe_face_values.reinit(cell, f); 
                    stokes_fe_face_values.reinit(cell->neighbor(f), 
                                                 cell->neighbor_of_neighbor(f)); 

                    assemble_interface_term(elasticity_fe_face_values, 
                                            stokes_fe_face_values, 
                                            elasticity_phi, 
                                            stokes_symgrad_phi_u, 
                                            stokes_phi_p, 
                                            local_interface_matrix); 

                    cell->neighbor(f)->get_dof_indices(neighbor_dof_indices); 
                    constraints.distribute_local_to_global( 
                      local_interface_matrix, 
                      local_dof_indices, 
                      neighbor_dof_indices, 
                      system_matrix); 
                  } 


                else if ((cell->neighbor(f)->level() == cell->level()) && 
                         (cell->neighbor(f)->has_children() == true)) 
                  { 
                    for (unsigned int subface = 0; 
                         subface < cell->face(f)->n_children(); 
                         ++subface) 
                      if (cell_is_in_fluid_domain( 
                            cell->neighbor_child_on_subface(f, subface))) 
                        { 
                          elasticity_fe_subface_values.reinit(cell, f, subface); 
                          stokes_fe_face_values.reinit( 
                            cell->neighbor_child_on_subface(f, subface), 
                            cell->neighbor_of_neighbor(f)); 

                          assemble_interface_term(elasticity_fe_subface_values, 
                                                  stokes_fe_face_values, 
                                                  elasticity_phi, 
                                                  stokes_symgrad_phi_u, 
                                                  stokes_phi_p, 
                                                  local_interface_matrix); 

                          cell->neighbor_child_on_subface(f, subface) 
                            ->get_dof_indices(neighbor_dof_indices); 
                          constraints.distribute_local_to_global( 
                            local_interface_matrix, 
                            local_dof_indices, 
                            neighbor_dof_indices, 
                            system_matrix); 
                        } 
                  } 


                else if (cell->neighbor_is_coarser(f) && 
                         cell_is_in_fluid_domain(cell->neighbor(f))) 
                  { 
                    elasticity_fe_face_values.reinit(cell, f); 
                    stokes_fe_subface_values.reinit( 
                      cell->neighbor(f), 
                      cell->neighbor_of_coarser_neighbor(f).first, 
                      cell->neighbor_of_coarser_neighbor(f).second); 

                    assemble_interface_term(elasticity_fe_face_values, 
                                            stokes_fe_subface_values, 
                                            elasticity_phi, 
                                            stokes_symgrad_phi_u, 
                                            stokes_phi_p, 
                                            local_interface_matrix); 

                    cell->neighbor(f)->get_dof_indices(neighbor_dof_indices); 
                    constraints.distribute_local_to_global( 
                      local_interface_matrix, 
                      local_dof_indices, 
                      neighbor_dof_indices, 
                      system_matrix); 
                  } 
              } 
      } 
  } 


  template <int dim> 
  void FluidStructureProblem<dim>::assemble_interface_term( 
    const FEFaceValuesBase<dim> &         elasticity_fe_face_values, 
    const FEFaceValuesBase<dim> &         stokes_fe_face_values, 
    std::vector<Tensor<1, dim>> &         elasticity_phi, 
    std::vector<SymmetricTensor<2, dim>> &stokes_symgrad_phi_u, 
    std::vector<double> &                 stokes_phi_p, 
    FullMatrix<double> &                  local_interface_matrix) const 
  { 
    Assert(stokes_fe_face_values.n_quadrature_points == 
             elasticity_fe_face_values.n_quadrature_points, 
           ExcInternalError()); 
    const unsigned int n_face_quadrature_points = 
      elasticity_fe_face_values.n_quadrature_points; 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 
    const FEValuesExtractors::Vector displacements(dim + 1); 

    local_interface_matrix = 0; 
    for (unsigned int q = 0; q < n_face_quadrature_points; ++q) 
      { 
        const Tensor<1, dim> normal_vector = 
          elasticity_fe_face_values.normal_vector(q); 

        for (unsigned int k = 0; k < stokes_fe_face_values.dofs_per_cell; ++k) 
          { 
            stokes_symgrad_phi_u[k] = 
              stokes_fe_face_values[velocities].symmetric_gradient(k, q); 
            stokes_phi_p[k] = stokes_fe_face_values[pressure].value(k, q); 
          } 
        for (unsigned int k = 0; k < elasticity_fe_face_values.dofs_per_cell; 
             ++k) 
          elasticity_phi[k] = 
            elasticity_fe_face_values[displacements].value(k, q); 

        for (unsigned int i = 0; i < elasticity_fe_face_values.dofs_per_cell; 
             ++i) 
          for (unsigned int j = 0; j < stokes_fe_face_values.dofs_per_cell; ++j) 
            local_interface_matrix(i, j) += 
              -((2 * viscosity * (stokes_symgrad_phi_u[j] * normal_vector) - 
                 stokes_phi_p[j] * normal_vector) * 
                elasticity_phi[i] * stokes_fe_face_values.JxW(q)); 
      } 
  } 


  template <int dim> 
  void FluidStructureProblem<dim>::solve() 
  { 
    SparseDirectUMFPACK direct_solver; 
    direct_solver.initialize(system_matrix); 
    direct_solver.vmult(solution, system_rhs); 

    constraints.distribute(solution); 
  } 



  template <int dim> 
  void FluidStructureProblem<dim>::output_results( 
    const unsigned int refinement_cycle) const 
  { 
    std::vector<std::string> solution_names(dim, "velocity"); 
    solution_names.emplace_back("pressure"); 
    for (unsigned int d = 0; d < dim; ++d) 
      solution_names.emplace_back("displacement"); 

    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    for (unsigned int d = 0; d < dim; ++d) 
      data_component_interpretation.push_back( 
        DataComponentInterpretation::component_is_part_of_vector); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 

    data_out.add_data_vector(solution, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.build_patches(); 

    std::ofstream output( 
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk"); 
    data_out.write_vtk(output); 
  } 


  template <int dim> 
  void FluidStructureProblem<dim>::refine_mesh() 
  { 
    Vector<float> stokes_estimated_error_per_cell( 
      triangulation.n_active_cells()); 
    Vector<float> elasticity_estimated_error_per_cell( 
      triangulation.n_active_cells()); 

    const QGauss<dim - 1> stokes_face_quadrature(stokes_degree + 2); 
    const QGauss<dim - 1> elasticity_face_quadrature(elasticity_degree + 2); 

    hp::QCollection<dim - 1> face_q_collection; 
    face_q_collection.push_back(stokes_face_quadrature); 
    face_q_collection.push_back(elasticity_face_quadrature); 

    const FEValuesExtractors::Vector velocities(0); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      face_q_collection, 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      stokes_estimated_error_per_cell, 
      fe_collection.component_mask(velocities)); 

    const FEValuesExtractors::Vector displacements(dim + 1); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      face_q_collection, 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      elasticity_estimated_error_per_cell, 
      fe_collection.component_mask(displacements)); 


    stokes_estimated_error_per_cell *= 
      4. / stokes_estimated_error_per_cell.l2_norm(); 
    elasticity_estimated_error_per_cell *= 
      1. / elasticity_estimated_error_per_cell.l2_norm(); 

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    estimated_error_per_cell += stokes_estimated_error_per_cell; 
    estimated_error_per_cell += elasticity_estimated_error_per_cell; 



    for (const auto &cell : dof_handler.active_cell_iterators()) 
      for (const auto f : cell->face_indices()) 
        if (cell_is_in_solid_domain(cell)) 
          { 
            if ((cell->at_boundary(f) == false) && 
                (((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == false) && 
                  cell_is_in_fluid_domain(cell->neighbor(f))) || 
                 ((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == true) && 
                  (cell_is_in_fluid_domain( 
                    cell->neighbor_child_on_subface(f, 0)))) || 
                 (cell->neighbor_is_coarser(f) && 
                  cell_is_in_fluid_domain(cell->neighbor(f))))) 
              estimated_error_per_cell(cell->active_cell_index()) = 0; 
          } 
        else 
          { 
            if ((cell->at_boundary(f) == false) && 
                (((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == false) && 
                  cell_is_in_solid_domain(cell->neighbor(f))) || 
                 ((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == true) && 
                  (cell_is_in_solid_domain( 
                    cell->neighbor_child_on_subface(f, 0)))) || 
                 (cell->neighbor_is_coarser(f) && 
                  cell_is_in_solid_domain(cell->neighbor(f))))) 
              estimated_error_per_cell(cell->active_cell_index()) = 0; 
          } 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.0); 
    triangulation.execute_coarsening_and_refinement(); 
  } 



  template <int dim> 
  void FluidStructureProblem<dim>::run() 
  { 
    make_grid(); 

    for (unsigned int refinement_cycle = 0; refinement_cycle < 10 - 2 * dim; 
         ++refinement_cycle) 
      { 
        std::cout << "Refinement cycle " << refinement_cycle << std::endl; 

        if (refinement_cycle > 0) 
          refine_mesh(); 

        setup_dofs(); 

        std::cout << "   Assembling..." << std::endl; 
        assemble_system(); 

        std::cout << "   Solving..." << std::endl; 
        solve(); 

        std::cout << "   Writing output..." << std::endl; 
        output_results(refinement_cycle); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step46 



int main() 
{ 
  try 
    { 
      using namespace Step46; 

      FluidStructureProblem<2> flow_problem(1, 1); 
      flow_problem.run(); 
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


