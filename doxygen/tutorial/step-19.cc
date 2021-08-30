

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2020 - 2021 by the deal.II authors 
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
 * Authors: Wolfgang Bangerth, Rene Gassmoeller, Peter Munch, 2020. 
 */ 




#include <deal.II/base/quadrature_lib.h> 

#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/fe/mapping_q.h> 
#include <deal.II/matrix_free/fe_point_evaluation.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/error_estimator.h> 


#include <deal.II/base/discrete_time.h> 
#include <deal.II/particles/particle_handler.h> 
#include <deal.II/particles/data_out.h> 

#include <fstream> 

using namespace dealii; 




namespace Step19 
{ 
  namespace BoundaryIds 
  { 
    constexpr types::boundary_id open          = 101; 
    constexpr types::boundary_id cathode       = 102; 
    constexpr types::boundary_id focus_element = 103; 
    constexpr types::boundary_id anode         = 104; 
  } // namespace BoundaryIds 

  namespace Constants 
  { 
    constexpr double electron_mass   = 9.1093837015e-31; 
    constexpr double electron_charge = 1.602176634e-19; 

    constexpr double V0 = 1; 

    constexpr double E_threshold = 0.05; 

    constexpr double electrons_per_particle = 3e15; 
  } // namespace Constants 


  template <int dim> 
  class CathodeRaySimulator 
  { 
  public: 
    CathodeRaySimulator(); 

    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void solve_field(); 
    void refine_grid(); 

    void create_particles(); 
    void move_particles(); 
    void track_lost_particle( 
      const typename Particles::ParticleIterator<dim> &        particle, 
      const typename Triangulation<dim>::active_cell_iterator &cell); 

    void update_timestep_size(); 
    void output_results() const; 

    Triangulation<dim>        triangulation; 
    MappingQGeneric<dim>      mapping; 
    FE_Q<dim>                 fe; 
    DoFHandler<dim>           dof_handler; 
    AffineConstraints<double> constraints; 

    SparseMatrix<double> system_matrix; 
    SparsityPattern      sparsity_pattern; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    Particles::ParticleHandler<dim> particle_handler; 
    types::particle_index           next_unused_particle_id; 
    types::particle_index           n_recently_lost_particles; 
    types::particle_index           n_total_lost_particles; 
    types::particle_index           n_particles_lost_through_anode; 

    DiscreteTime time; 
  }; 





  template <int dim> 
  CathodeRaySimulator<dim>::CathodeRaySimulator() 
    : mapping(1) 
    , fe(2) 
    , dof_handler(triangulation) 
    , particle_handler(triangulation, mapping, /*n_properties=*/dim) 
    , next_unused_particle_id(0) 
    , n_recently_lost_particles(0) 
    , n_total_lost_particles(0) 
    , n_particles_lost_through_anode(0) 
    , time(0, 1e-4) 
  { 
    particle_handler.signals.particle_lost.connect( 
      [this](const typename Particles::ParticleIterator<dim> &        particle, 
             const typename Triangulation<dim>::active_cell_iterator &cell) { 
        this->track_lost_particle(particle, cell); 
      }); 
  } 




  template <int dim> 
  void CathodeRaySimulator<dim>::make_grid() 
  { 
    static_assert(dim == 2, 
                  "This function is currently only implemented for 2d."); 

    const double       delta = 0.5; 
    const unsigned int nx    = 5; 
    const unsigned int ny    = 3; 

    const std::vector<Point<dim>> vertices // 
      = {{0, 0}, 
         {1, 0}, 
         {2, 0}, 
         {3, 0}, 
         {4, 0}, 
         {delta, 1}, 
         {1, 1}, 
         {2, 1}, 
         {3, 1}, 
         {4, 1}, 
         {0, 2}, 
         {1, 2}, 
         {2, 2}, 
         {3, 2}, 
         {4, 2}}; 
    AssertDimension(vertices.size(), nx * ny); 

    const std::vector<unsigned int> cell_vertices[(nx - 1) * (ny - 1)] = { 
      {0, 1, nx + 0, nx + 1}, 
      {1, 2, nx + 1, nx + 2}, 
      {2, 3, nx + 2, nx + 3}, 
      {3, 4, nx + 3, nx + 4}, 

      {5, nx + 1, 2 * nx + 0, 2 * nx + 1}, 
      {nx + 1, nx + 2, 2 * nx + 1, 2 * nx + 2}, 
      {nx + 2, nx + 3, 2 * nx + 2, 2 * nx + 3}, 
      {nx + 3, nx + 4, 2 * nx + 3, 2 * nx + 4}}; 



    std::vector<CellData<dim>> cells((nx - 1) * (ny - 1), CellData<dim>()); 
    for (unsigned int i = 0; i < cells.size(); ++i) 
      { 
        cells[i].vertices    = cell_vertices[i]; 
        cells[i].material_id = 0; 
      } 

    triangulation.create_triangulation( 
      vertices, 
      cells, 
      SubCellData()); // No boundary information 

    triangulation.refine_global(2); 



    for (auto &cell : triangulation.active_cell_iterators()) 
      for (auto &face : cell->face_iterators()) 
        if (face->at_boundary()) 
          { 
            if ((face->center()[0] > 0) && (face->center()[0] < 0.5) && 
                (face->center()[1] > 0) && (face->center()[1] < 2)) 
              face->set_boundary_id(BoundaryIds::cathode); 
            else if ((face->center()[0] > 0) && (face->center()[0] < 2)) 
              face->set_boundary_id(BoundaryIds::focus_element); 
            else if ((face->center()[0] > 4 - 1e-12) && 
                     ((face->center()[1] > 1.5) || (face->center()[1] < 0.5))) 
              face->set_boundary_id(BoundaryIds::anode); 
            else 
              face->set_boundary_id(BoundaryIds::open); 
          } 

    triangulation.refine_global(1); 
  } 


  template <int dim> 
  void CathodeRaySimulator<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    VectorTools::interpolate_boundary_values(dof_handler, 
                                             BoundaryIds::cathode, 
                                             Functions::ConstantFunction<dim>( 
                                               -Constants::V0), 
                                             constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             BoundaryIds::focus_element, 
                                             Functions::ConstantFunction<dim>( 
                                               -Constants::V0), 
                                             constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             BoundaryIds::anode, 
                                             Functions::ConstantFunction<dim>( 
                                               +Constants::V0), 
                                             constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    dsp, 
                                    constraints, 
                                    /*keep_constrained_dofs =  */ false);

    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
  } 


  template <int dim> 
  void CathodeRaySimulator<dim>::assemble_system() 
  { 
    system_matrix = 0; 
    system_rhs    = 0; 

    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.dofs_per_cell; 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        fe_values.reinit(cell); 

        for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
          for (const unsigned int i : fe_values.dof_indices()) 
            { 
              for (const unsigned int j : fe_values.dof_indices()) 
                cell_matrix(i, j) += 
                  (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
                   fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
                   fe_values.JxW(q_index));           // dx 
            } 




        if (particle_handler.n_particles_in_cell(cell) > 0) 
          for (const auto &particle : particle_handler.particles_in_cell(cell)) 
            { 
              const Point<dim> &reference_location = 
                particle.get_reference_location(); 
              for (const unsigned int i : fe_values.dof_indices()) 
                cell_rhs(i) += 
                  (fe.shape_value(i, reference_location) * // phi_i(x_p) 
                   (-Constants::electrons_per_particle *   // N 
                    Constants::electron_charge));          // e 
            } 


        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global( 
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 


  template <int dim> 
  void CathodeRaySimulator<dim>::solve_field() 
  { 
    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> solver(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    solver.solve(system_matrix, solution, system_rhs, preconditioner); 

    constraints.distribute(solution); 
  } 


  template <int dim> 
  void CathodeRaySimulator<dim>::refine_grid() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate(dof_handler, 
                                       QGauss<dim - 1>(fe.degree + 1), 
                                       {}, 
                                       solution, 
                                       estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.1, 
                                                    0.03); 

    triangulation.execute_coarsening_and_refinement(); 
  } 



  template <int dim> 
  void CathodeRaySimulator<dim>::create_particles() 
  { 
    FEFaceValues<dim> fe_face_values(fe, 
                                     QMidpoint<dim - 1>(), 
                                     update_quadrature_points | 
                                       update_gradients | 
                                       update_normal_vectors); 

    std::vector<Tensor<1, dim>> solution_gradients( 
      fe_face_values.n_quadrature_points); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary() && 
            (face->boundary_id() == BoundaryIds::cathode)) 
          { 
            fe_face_values.reinit(cell, face); 


            const FEValuesExtractors::Scalar electric_potential(0); 
            fe_face_values[electric_potential].get_function_gradients( 
              solution, solution_gradients); 
            for (const unsigned int q_point : 
                 fe_face_values.quadrature_point_indices()) 
              { 
                const Tensor<1, dim> E = solution_gradients[q_point]; 


                if ((E * fe_face_values.normal_vector(q_point) < 0) && 
                    (E.norm() > Constants::E_threshold)) 
                  { 
                    const Point<dim> &location = 
                      fe_face_values.quadrature_point(q_point); 

                    Particles::Particle<dim> new_particle; 
                    new_particle.set_location(location); 
                    new_particle.set_reference_location( 
                      mapping.transform_real_to_unit_cell(cell, location)); 
                    new_particle.set_id(next_unused_particle_id); 
                    particle_handler.insert_particle(new_particle, cell); 

                    ++next_unused_particle_id; 
                  } 
              } 
          } 


    particle_handler.update_cached_numbers(); 
  } 



  template <int dim> 
  void CathodeRaySimulator<dim>::move_particles() 
  { 
    const double dt = time.get_next_step_size(); 

    Vector<double>            solution_values(fe.n_dofs_per_cell()); 
    FEPointEvaluation<1, dim> evaluator(mapping, fe, update_gradients); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (particle_handler.n_particles_in_cell(cell) > 0) 
        { 
          const typename Particles::ParticleHandler< 
            dim>::particle_iterator_range particles_in_cell = 
            particle_handler.particles_in_cell(cell); 

          std::vector<Point<dim>> particle_positions; 
          for (const auto &particle : particles_in_cell) 
            particle_positions.push_back(particle.get_reference_location()); 

          cell->get_dof_values(solution, solution_values); 


          evaluator.reinit(cell, particle_positions); 
          evaluator.evaluate(make_array_view(solution_values), 
                             EvaluationFlags::gradients); 

          { 
            typename Particles::ParticleHandler<dim>::particle_iterator 
              particle = particles_in_cell.begin(); 
            for (unsigned int particle_index = 0; 
                 particle != particles_in_cell.end(); 
                 ++particle, ++particle_index) 
              { 
                const Tensor<1, dim> &E = 
                  evaluator.get_gradient(particle_index); 


                const Tensor<1, dim> old_velocity(particle->get_properties()); 

                const Tensor<1, dim> acceleration = 
                  Constants::electron_charge / Constants::electron_mass * E; 

                const Tensor<1, dim> new_velocity = 
                  old_velocity + acceleration * dt; 

                particle->set_properties(make_array_view(new_velocity)); 


                const Point<dim> new_location = 
                  particle->get_location() + dt * new_velocity; 
                particle->set_location(new_location); 
              } 
          } 
        } 


    particle_handler.sort_particles_into_subdomains_and_cells(); 
  } 


  template <int dim> 
  void CathodeRaySimulator<dim>::track_lost_particle( 
    const typename Particles::ParticleIterator<dim> &        particle, 
    const typename Triangulation<dim>::active_cell_iterator &cell) 
  { 
    ++n_recently_lost_particles; 
    ++n_total_lost_particles; 

    const Point<dim> current_location              = particle->get_location(); 
    const Point<dim> approximate_previous_location = cell->center(); 

    if ((approximate_previous_location[0] < 4) && (current_location[0] > 4)) 
      { 
        const Tensor<1, dim> direction = 
          (current_location - approximate_previous_location) / 
          (current_location[0] - approximate_previous_location[0]); 

        const double right_boundary_intercept = 
          approximate_previous_location[1] + 
          (4 - approximate_previous_location[0]) * direction[1]; 
        if ((right_boundary_intercept > 0.5) && 
            (right_boundary_intercept < 1.5)) 
          ++n_particles_lost_through_anode; 
      } 
  } 



  template <int dim> 
  void CathodeRaySimulator<dim>::update_timestep_size() 
  { 
    if (time.get_step_number() > 0) 
      { 
        double min_cell_size_over_velocity = std::numeric_limits<double>::max(); 

        for (const auto &cell : dof_handler.active_cell_iterators()) 
          if (particle_handler.n_particles_in_cell(cell) > 0) 
            { 
              const double cell_size = cell->minimum_vertex_distance(); 

              double max_particle_velocity(0.0); 

              for (const auto &particle : 
                   particle_handler.particles_in_cell(cell)) 
                { 
                  const Tensor<1, dim> velocity(particle.get_properties()); 
                  max_particle_velocity = 
                    std::max(max_particle_velocity, velocity.norm()); 
                } 

              if (max_particle_velocity > 0) 
                min_cell_size_over_velocity = 
                  std::min(min_cell_size_over_velocity, 
                           cell_size / max_particle_velocity); 
            } 

        constexpr double c_safety = 0.5; 
        time.set_desired_next_step_size(c_safety * 0.5 * 
                                        min_cell_size_over_velocity); 
      } 


    else 
      { 
        const QTrapezoid<dim> vertex_quadrature; 
        FEValues<dim> fe_values(fe, vertex_quadrature, update_gradients); 

        std::vector<Tensor<1, dim>> field_gradients(vertex_quadrature.size()); 

        double min_timestep = std::numeric_limits<double>::max(); 

        for (const auto &cell : dof_handler.active_cell_iterators()) 
          if (particle_handler.n_particles_in_cell(cell) > 0) 
            { 
              const double cell_size = cell->minimum_vertex_distance(); 

              fe_values.reinit(cell); 
              fe_values.get_function_gradients(solution, field_gradients); 

              double max_E = 0; 
              for (const auto q_point : fe_values.quadrature_point_indices()) 
                max_E = std::max(max_E, field_gradients[q_point].norm()); 

              if (max_E > 0) 
                min_timestep = 
                  std::min(min_timestep, 
                           std::sqrt(0.5 * cell_size * 
                                     Constants::electron_mass / 
                                     Constants::electron_charge / max_E)); 
            } 

        time.set_desired_next_step_size(min_timestep); 
      } 
  } 





  template <int dim> 
  class ElectricFieldPostprocessor : public DataPostprocessorVector<dim> 
  { 
  public: 
    ElectricFieldPostprocessor() 
      : DataPostprocessorVector<dim>("electric_field", update_gradients) 
    {} 

    virtual void evaluate_scalar_field( 
      const DataPostprocessorInputs::Scalar<dim> &input_data, 
      std::vector<Vector<double>> &computed_quantities) const override 
    { 
      AssertDimension(input_data.solution_gradients.size(), 
                      computed_quantities.size()); 

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p) 
        { 
          AssertDimension(computed_quantities[p].size(), dim); 
          for (unsigned int d = 0; d < dim; ++d) 
            computed_quantities[p][d] = input_data.solution_gradients[p][d]; 
        } 
    } 
  }; 


  template <int dim> 
  void CathodeRaySimulator<dim>::output_results() const 
  { 
    { 
      ElectricFieldPostprocessor<dim> electric_field; 
      DataOut<dim>                    data_out; 
      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution, "electric_potential"); 
      data_out.add_data_vector(solution, electric_field); 
      data_out.build_patches(); 

      data_out.set_flags( 
        DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number())); 

      std::ofstream output("solution-" + 
                           Utilities::int_to_string(time.get_step_number(), 4) + 
                           ".vtu"); 
      data_out.write_vtu(output); 
    } 


    { 
      Particles::DataOut<dim, dim> particle_out; 
      particle_out.build_patches( 
        particle_handler, 
        std::vector<std::string>(dim, "velocity"), 
        std::vector<DataComponentInterpretation::DataComponentInterpretation>( 
          dim, DataComponentInterpretation::component_is_part_of_vector)); 

      particle_out.set_flags( 
        DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number())); 

      std::ofstream output("particles-" + 
                           Utilities::int_to_string(time.get_step_number(), 4) + 
                           ".vtu"); 
      particle_out.write_vtu(output); 
    } 
  } 


  template <int dim> 
  void CathodeRaySimulator<dim>::run() 
  { 
    make_grid(); 


    const unsigned int n_pre_refinement_cycles = 3; 
    for (unsigned int refinement_cycle = 0; 
         refinement_cycle < n_pre_refinement_cycles; 
         ++refinement_cycle) 
      { 
        setup_system(); 
        assemble_system(); 
        solve_field(); 
        refine_grid(); 
      } 


    setup_system(); 
    do 
      { 
        std::cout << "Timestep " << time.get_step_number() + 1 << std::endl; 
        std::cout << "  Field degrees of freedom:                 " 
                  << dof_handler.n_dofs() << std::endl; 

        assemble_system(); 
        solve_field(); 

        create_particles(); 
        std::cout << "  Total number of particles in simulation:  " 
                  << particle_handler.n_global_particles() << std::endl; 

        n_recently_lost_particles = 0; 
        update_timestep_size(); 
        move_particles(); 

        time.advance_time(); 

        output_results(); 

        std::cout << "  Number of particles lost this time step:  " 
                  << n_recently_lost_particles << std::endl; 
        if (n_total_lost_particles > 0) 
          std::cout << "  Fraction of particles lost through anode: " 
                    << 1. * n_particles_lost_through_anode / 
                         n_total_lost_particles 
                    << std::endl; 

        std::cout << std::endl 
                  << "  Now at t=" << time.get_current_time() 
                  << ", dt=" << time.get_previous_step_size() << '.' 
                  << std::endl 
                  << std::endl; 
      } 
    while (time.is_at_end() == false); 
  } 
} // namespace Step19 



int main() 
{ 
  try 
    { 
      Step19::CathodeRaySimulator<2> cathode_ray_simulator_2d; 
      cathode_ray_simulator_2d.run(); 
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


