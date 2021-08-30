

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
 * Authors: Bruno Blais, Toni El Geitani Nehme, Rene Gassmoeller, Peter Munch 
 */ 


#include <deal.II/base/bounding_box.h> 
#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/discrete_time.h> 
#include <deal.II/base/mpi.h> 
#include <deal.II/base/parameter_acceptor.h> 
#include <deal.II/base/timer.h> 

#include <deal.II/distributed/cell_weights.h> 
#include <deal.II/distributed/solution_transfer.h> 
#include <deal.II/distributed/tria.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 

#include <deal.II/lac/la_parallel_vector.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 


#include <deal.II/particles/particle_handler.h> 


#include <deal.II/particles/generators.h> 


#include <deal.II/particles/data_out.h> 

#include <cmath> 
#include <iostream> 

namespace Step68 
{ 
  using namespace dealii; 




  class ParticleTrackingParameters : public ParameterAcceptor 
  { 
  public: 
    ParticleTrackingParameters(); 


    std::string output_directory = "./"; 

    unsigned int velocity_degree       = 1; 
    double       time_step             = 0.002; 
    double       final_time            = 4.0; 
    unsigned int output_frequency      = 10; 
    unsigned int repartition_frequency = 5; 


    unsigned int fluid_refinement              = 4; 
    unsigned int particle_insertion_refinement = 3; 
  }; 


  ParticleTrackingParameters::ParticleTrackingParameters() 
    : ParameterAcceptor("Particle Tracking Problem/") 
  { 
    add_parameter( 
      "Velocity degree", velocity_degree, "", prm, Patterns::Integer(1)); 

    add_parameter("Output frequency", 
                  output_frequency, 
                  "Iteration frequency at which output results are written", 
                  prm, 
                  Patterns::Integer(1)); 

    add_parameter("Repartition frequency", 
                  repartition_frequency, 
                  "Iteration frequency at which the mesh is load balanced", 
                  prm, 
                  Patterns::Integer(1)); 

    add_parameter("Output directory", output_directory); 

    add_parameter("Time step", time_step, "", prm, Patterns::Double()); 

    add_parameter("Final time", 
                  final_time, 
                  "End time of the simulation", 
                  prm, 
                  Patterns::Double()); 

    add_parameter("Fluid refinement", 
                  fluid_refinement, 
                  "Refinement level of the fluid domain", 
                  prm, 
                  Patterns::Integer(0)); 

    add_parameter( 
      "Particle insertion refinement", 
      particle_insertion_refinement, 
      "Refinement of the volumetric mesh used to insert the particles", 
      prm, 
      Patterns::Integer(0)); 
  } 



  template <int dim> 
  class Vortex : public Function<dim> 
  { 
  public: 
    Vortex() 
      : Function<dim>(dim) 
    {} 

    virtual void vector_value(const Point<dim> &point, 
                              Vector<double> &  values) const override; 
  }; 


  template <int dim> 
  void Vortex<dim>::vector_value(const Point<dim> &point, 
                                 Vector<double> &  values) const 
  { 
    const double T = 4; 
    const double t = this->get_time(); 

    const double px = numbers::PI * point(0); 
    const double py = numbers::PI * point(1); 
    const double pt = numbers::PI / T * t; 

    values[0] = -2 * cos(pt) * pow(sin(px), 2) * sin(py) * cos(py); 
    values[1] = 2 * cos(pt) * pow(sin(py), 2) * sin(px) * cos(px); 
    if (dim == 3) 
      { 
        values[2] = 0; 
      } 
  } 



  template <int dim> 
  class ParticleTracking 
  { 
  public: 
    ParticleTracking(const ParticleTrackingParameters &par, 
                     const bool                        interpolated_velocity); 
    void run(); 

  private: 


    void generate_particles(); 


    void setup_background_dofs(); 


    void interpolate_function_to_field(); 


    void euler_step_interpolated(const double dt); 
    void euler_step_analytical(const double dt); 


    unsigned int cell_weight( 
      const typename parallel::distributed::Triangulation<dim>::cell_iterator 
        &cell, 
      const typename parallel::distributed::Triangulation<dim>::CellStatus 
        status) const; 


    void output_particles(const unsigned int it); 
    void output_background(const unsigned int it); 


    const ParticleTrackingParameters &par; 

    MPI_Comm                                  mpi_communicator; 
    parallel::distributed::Triangulation<dim> background_triangulation; 
    Particles::ParticleHandler<dim>           particle_handler; 

    DoFHandler<dim>                            fluid_dh; 
    FESystem<dim>                              fluid_fe; 
    MappingQ1<dim>                             mapping; 
    LinearAlgebra::distributed::Vector<double> velocity_field; 

    Vortex<dim> velocity; 

    ConditionalOStream pcout; 

    bool interpolated_velocity; 
  }; 



  template <int dim> 
  ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par, 
                                          const bool interpolated_velocity) 
    : par(par) 
    , mpi_communicator(MPI_COMM_WORLD) 
    , background_triangulation(mpi_communicator) 
    , fluid_dh(background_triangulation) 
    , fluid_fe(FE_Q<dim>(par.velocity_degree), dim) 
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
    , interpolated_velocity(interpolated_velocity) 

  {} 



  template <int dim> 
  unsigned int ParticleTracking<dim>::cell_weight( 
    const typename parallel::distributed::Triangulation<dim>::cell_iterator 
      &                                                                  cell, 
    const typename parallel::distributed::Triangulation<dim>::CellStatus status) 
    const 
  { 


    if (!cell->is_locally_owned()) 
      return 0; 


    const unsigned int particle_weight = 10000; 


    if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST || 
        status == parallel::distributed::Triangulation<dim>::CELL_REFINE) 
      { 
        const unsigned int n_particles_in_cell = 
          particle_handler.n_particles_in_cell(cell); 
        return n_particles_in_cell * particle_weight; 
      } 
    else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN) 
      { 
        unsigned int n_particles_in_cell = 0; 

        for (unsigned int child_index = 0; child_index < cell->n_children(); 
             ++child_index) 
          n_particles_in_cell += 
            particle_handler.n_particles_in_cell(cell->child(child_index)); 

        return n_particles_in_cell * particle_weight; 
      } 

    Assert(false, ExcInternalError()); 
    return 0; 
  } 



  template <int dim> 
  void ParticleTracking<dim>::generate_particles() 
  { 


    GridGenerator::hyper_cube(background_triangulation, 0, 1); 
    background_triangulation.refine_global(par.fluid_refinement); 




    background_triangulation.signals.cell_weight.connect( 
      [&]( 
        const typename parallel::distributed::Triangulation<dim>::cell_iterator 
          &cell, 
        const typename parallel::distributed::Triangulation<dim>::CellStatus 
          status) -> unsigned int { return this->cell_weight(cell, status); }); 

    background_triangulation.signals.pre_distributed_repartition.connect( 
      [this]() { this->particle_handler.register_store_callback_function(); }); 

    background_triangulation.signals.post_distributed_repartition.connect( 
      [&]() { this->particle_handler.register_load_callback_function(false); }); 


    particle_handler.initialize(background_triangulation, mapping, 1 + dim); 


    Point<dim> center; 
    center[0] = 0.5; 
    center[1] = 0.75; 
    if (dim == 3) 
      center[2] = 0.5; 

    const double outer_radius = 0.15; 
    const double inner_radius = 0.01; 

    parallel::distributed::Triangulation<dim> particle_triangulation( 
      MPI_COMM_WORLD); 

    GridGenerator::hyper_shell( 
      particle_triangulation, center, inner_radius, outer_radius, 6); 
    particle_triangulation.refine_global(par.particle_insertion_refinement); 


    const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box( 
      background_triangulation, IteratorFilters::LocallyOwnedCell()); 
    const auto global_bounding_boxes = 
      Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box); 


    std::vector<std::vector<double>> properties( 
      particle_triangulation.n_locally_owned_active_cells(), 
      std::vector<double>(dim + 1, 0.)); 


    Particles::Generators::quadrature_points(particle_triangulation, 
                                             QMidpoint<dim>(), 
                                             global_bounding_boxes, 
                                             particle_handler, 
                                             mapping, 
                                             properties); 

    pcout << "Number of particles inserted: " 
          << particle_handler.n_global_particles() << std::endl; 
  } 



  template <int dim> 
  void ParticleTracking<dim>::setup_background_dofs() 
  { 
    fluid_dh.distribute_dofs(fluid_fe); 
    const IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs(); 
    IndexSet       locally_relevant_dofs; 
    DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs); 

    velocity_field.reinit(locally_owned_dofs, 
                          locally_relevant_dofs, 
                          mpi_communicator); 
  } 


  template <int dim> 
  void ParticleTracking<dim>::interpolate_function_to_field() 
  { 
    velocity_field.zero_out_ghost_values(); 
    VectorTools::interpolate(mapping, fluid_dh, velocity, velocity_field); 
    velocity_field.update_ghost_values(); 
  } 



  template <int dim> 
  void ParticleTracking<dim>::euler_step_analytical(const double dt) 
  { 
    const unsigned int this_mpi_rank = 
      Utilities::MPI::this_mpi_process(mpi_communicator); 
    Vector<double> particle_velocity(dim); 


    for (auto &particle : particle_handler) 
      { 


        Point<dim> particle_location = particle.get_location(); 
        velocity.vector_value(particle_location, particle_velocity); 


        for (int d = 0; d < dim; ++d) 
          particle_location[d] += particle_velocity[d] * dt; 

        particle.set_location(particle_location); 


        ArrayView<double> properties = particle.get_properties(); 
        for (int d = 0; d < dim; ++d) 
          properties[d] = particle_velocity[d]; 
        properties[dim] = this_mpi_rank; 
      } 
  } 


  template <int dim> 
  void ParticleTracking<dim>::euler_step_interpolated(const double dt) 
  { 
    Vector<double> local_dof_values(fluid_fe.dofs_per_cell); 


    auto particle = particle_handler.begin(); 
    while (particle != particle_handler.end()) 
      { 
        const auto cell = 
          particle->get_surrounding_cell(background_triangulation); 
        const auto dh_cell = 
          typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh); 

        dh_cell->get_dof_values(velocity_field, local_dof_values); 


        const auto pic = particle_handler.particles_in_cell(cell); 
        Assert(pic.begin() == particle, ExcInternalError()); 
        for (auto &p : pic) 
          { 
            const Point<dim> reference_location = p.get_reference_location(); 
            Tensor<1, dim>   particle_velocity; 
            for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j) 
              { 
                const auto comp_j = fluid_fe.system_to_component_index(j); 

                particle_velocity[comp_j.first] += 
                  fluid_fe.shape_value(j, reference_location) * 
                  local_dof_values[j]; 
              } 

            Point<dim> particle_location = particle->get_location(); 
            for (int d = 0; d < dim; ++d) 
              particle_location[d] += particle_velocity[d] * dt; 
            p.set_location(particle_location); 


            ArrayView<double> properties = p.get_properties(); 
            for (int d = 0; d < dim; ++d) 
              properties[d] = particle_velocity[d]; 

            properties[dim] = 
              Utilities::MPI::this_mpi_process(mpi_communicator); 

            ++particle; 
          } 
      } 
  } 



  template <int dim> 
  void ParticleTracking<dim>::output_particles(const unsigned int it) 
  { 
    Particles::DataOut<dim, dim> particle_output; 

    std::vector<std::string> solution_names(dim, "velocity"); 
    solution_names.push_back("process_id"); 

    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 

    particle_output.build_patches(particle_handler, 
                                  solution_names, 
                                  data_component_interpretation); 
    const std::string output_folder(par.output_directory); 
    const std::string file_name(interpolated_velocity ? 
                                  "interpolated-particles" : 
                                  "analytical-particles"); 

    pcout << "Writing particle output file: " << file_name << "-" << it 
          << std::endl; 

    particle_output.write_vtu_with_pvtu_record( 
      output_folder, file_name, it, mpi_communicator, 6); 
  } 

  template <int dim> 
  void ParticleTracking<dim>::output_background(const unsigned int it) 
  { 
    std::vector<std::string> solution_names(dim, "velocity"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 

    DataOut<dim> data_out; 


    data_out.attach_dof_handler(fluid_dh); 
    data_out.add_data_vector(velocity_field, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    Vector<float> subdomain(background_triangulation.n_active_cells()); 
    for (unsigned int i = 0; i < subdomain.size(); ++i) 
      subdomain(i) = background_triangulation.locally_owned_subdomain(); 
    data_out.add_data_vector(subdomain, "subdomain"); 

    data_out.build_patches(mapping); 

    const std::string output_folder(par.output_directory); 
    const std::string file_name("background"); 

    pcout << "Writing background field file: " << file_name << "-" << it 
          << std::endl; 

    data_out.write_vtu_with_pvtu_record( 
      output_folder, file_name, it, mpi_communicator, 6); 
  } 


  template <int dim> 
  void ParticleTracking<dim>::run() 
  { 
    DiscreteTime discrete_time(0, par.final_time, par.time_step); 

    generate_particles(); 

    pcout << "Repartitioning triangulation after particle generation" 
          << std::endl; 
    background_triangulation.repartition(); 


    if (interpolated_velocity) 
      { 
        setup_background_dofs(); 
        interpolate_function_to_field(); 
        euler_step_interpolated(0.); 
      } 
    else 
      euler_step_analytical(0.); 

    output_particles(discrete_time.get_step_number()); 
    if (interpolated_velocity) 
      output_background(discrete_time.get_step_number()); 


    while (!discrete_time.is_at_end()) 
      { 
        discrete_time.advance_time(); 
        velocity.set_time(discrete_time.get_previous_time()); 

        if ((discrete_time.get_step_number() % par.repartition_frequency) == 0) 
          { 
            background_triangulation.repartition(); 
            if (interpolated_velocity) 
              setup_background_dofs(); 
          } 

        if (interpolated_velocity) 
          { 
            interpolate_function_to_field(); 
            euler_step_interpolated(discrete_time.get_previous_step_size()); 
          } 
        else 
          euler_step_analytical(discrete_time.get_previous_step_size()); 

        particle_handler.sort_particles_into_subdomains_and_cells(); 

        if ((discrete_time.get_step_number() % par.output_frequency) == 0) 
          { 
            output_particles(discrete_time.get_step_number()); 
            if (interpolated_velocity) 
              output_background(discrete_time.get_step_number()); 
          } 
      } 
  } 

} // namespace Step68 



int main(int argc, char *argv[]) 
{ 
  using namespace Step68; 
  using namespace dealii; 
  deallog.depth_console(1); 

  try 
    { 
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

      std::string prm_file; 
      if (argc > 1) 
        prm_file = argv[1]; 
      else 
        prm_file = "parameters.prm"; 

      ParticleTrackingParameters par; 
      ParameterAcceptor::initialize(prm_file); 
      { 
        Step68::ParticleTracking<2> particle_tracking(par, false); 
        particle_tracking.run(); 
      } 
      { 
        Step68::ParticleTracking<2> particle_tracking(par, true); 
        particle_tracking.run(); 
      } 
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


