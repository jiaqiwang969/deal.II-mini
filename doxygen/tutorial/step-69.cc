

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2019 - 2021 by the deal.II authors 
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
 * Authors: Matthias Maier, Texas A&M University; 
 *          Ignacio Tomas, Texas A&M University, Sandia National Laboratories 
 * 
 * Sandia National Laboratories is a multimission laboratory managed and 
 * operated by National Technology & Engineering Solutions of Sandia, LLC, a 
 * wholly owned subsidiary of Honeywell International Inc., for the U.S. 
 * Department of Energy's National Nuclear Security Administration under 
 * contract DE-NA0003525. This document describes objective technical results 
 * and analysis. Any subjective views or opinions that might be expressed in 
 * the paper do not necessarily represent the views of the U.S. Department of 
 * Energy or the United States Government. 
 */ 




#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/parallel.h> 
#include <deal.II/base/parameter_acceptor.h> 
#include <deal.II/base/partitioner.h> 
#include <deal.II/base/quadrature.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/work_stream.h> 

#include <deal.II/distributed/tria.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/la_parallel_vector.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/sparse_matrix.templates.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/meshworker/scratch_data.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 


#include <boost/archive/binary_iarchive.hpp> 
#include <boost/archive/binary_oarchive.hpp> 


#include <deal.II/base/std_cxx20/iota_view.h> 
#include <boost/range/iterator_range.hpp> 

#include <cmath> 
#include <fstream> 
#include <future> 




namespace Step69 
{ 
  using namespace dealii; 


  namespace Boundaries 
  { 
    constexpr types::boundary_id do_nothing = 0; 
    constexpr types::boundary_id free_slip  = 1; 
    constexpr types::boundary_id dirichlet  = 2; 
  } // namespace Boundaries 

  template <int dim> 
  class Discretization : public ParameterAcceptor 
  { 
  public: 
    Discretization(const MPI_Comm     mpi_communicator, 
                   TimerOutput &      computing_timer, 
                   const std::string &subsection = "Discretization"); 

    void setup(); 

    const MPI_Comm mpi_communicator; 

    parallel::distributed::Triangulation<dim> triangulation; 

    const MappingQ<dim>   mapping; 
    const FE_Q<dim>       finite_element; 
    const QGauss<dim>     quadrature; 
    const QGauss<dim - 1> face_quadrature; 

  private: 
    TimerOutput &computing_timer; 

    double length; 
    double height; 
    double disk_position; 
    double disk_diameter; 

    unsigned int refinement; 
  }; 




  template <int dim> 
  class OfflineData : public ParameterAcceptor 
  { 
  public: 
    using BoundaryNormalMap = 
      std::map<types::global_dof_index, 
               std::tuple<Tensor<1, dim>, types::boundary_id, Point<dim>>>; 

    OfflineData(const MPI_Comm             mpi_communicator, 
                TimerOutput &              computing_timer, 
                const Discretization<dim> &discretization, 
                const std::string &        subsection = "OfflineData"); 

    void setup(); 
    void assemble(); 

    DoFHandler<dim> dof_handler; 

    std::shared_ptr<const Utilities::MPI::Partitioner> partitioner; 

    unsigned int n_locally_owned; 
    unsigned int n_locally_relevant; 

    SparsityPattern sparsity_pattern; 

    BoundaryNormalMap boundary_normal_map; 

    SparseMatrix<double>                  lumped_mass_matrix; 
    std::array<SparseMatrix<double>, dim> cij_matrix; 
    std::array<SparseMatrix<double>, dim> nij_matrix; 
    SparseMatrix<double>                  norm_matrix; 

  private: 
    const MPI_Comm mpi_communicator; 
    TimerOutput &  computing_timer; 

    SmartPointer<const Discretization<dim>> discretization; 
  }; 









  template <int dim> 
  class ProblemDescription 
  { 
  public: 
    static constexpr unsigned int problem_dimension = 2 + dim; 

    using state_type = Tensor<1, problem_dimension>; 
    using flux_type  = Tensor<1, problem_dimension, Tensor<1, dim>>; 

    const static std::array<std::string, problem_dimension> component_names; 

    static constexpr double gamma = 7. / 5.; 

    static DEAL_II_ALWAYS_INLINE inline Tensor<1, dim> 
    momentum(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double 
    internal_energy(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double pressure(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double 
    speed_of_sound(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline flux_type flux(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double 
    compute_lambda_max(const state_type &    U_i, 
                       const state_type &    U_j, 
                       const Tensor<1, dim> &n_ij); 
  }; 




  template <int dim> 
  class InitialValues : public ParameterAcceptor 
  { 
  public: 
    using state_type = typename ProblemDescription<dim>::state_type; 

    InitialValues(const std::string &subsection = "InitialValues"); 

    std::function<state_type(const Point<dim> &point, double t)> initial_state; 

  private: 


    void parse_parameters_callback(); 

    Tensor<1, dim> initial_direction; 
    Tensor<1, 3>   initial_1d_state; 
  }; 



  template <int dim> 
  class TimeStepping : public ParameterAcceptor 
  { 
  public: 
    static constexpr unsigned int problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 

    using state_type = typename ProblemDescription<dim>::state_type; 
    using flux_type  = typename ProblemDescription<dim>::flux_type; 

    using vector_type = 
      std::array<LinearAlgebra::distributed::Vector<double>, problem_dimension>; 

    TimeStepping(const MPI_Comm            mpi_communicator, 
                 TimerOutput &             computing_timer, 
                 const OfflineData<dim> &  offline_data, 
                 const InitialValues<dim> &initial_values, 
                 const std::string &       subsection = "TimeStepping"); 

    void prepare(); 

    double make_one_step(vector_type &U, double t); 

  private: 
    const MPI_Comm mpi_communicator; 
    TimerOutput &  computing_timer; 

    SmartPointer<const OfflineData<dim>>   offline_data; 
    SmartPointer<const InitialValues<dim>> initial_values; 

    SparseMatrix<double> dij_matrix; 

    vector_type temporary_vector; 

    double cfl_update; 
  }; 


  template <int dim> 
  class SchlierenPostprocessor : public ParameterAcceptor 
  { 
  public: 
    static constexpr unsigned int problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 

    using state_type = typename ProblemDescription<dim>::state_type; 

    using vector_type = 
      std::array<LinearAlgebra::distributed::Vector<double>, problem_dimension>; 

    SchlierenPostprocessor( 
      const MPI_Comm          mpi_communicator, 
      TimerOutput &           computing_timer, 
      const OfflineData<dim> &offline_data, 
      const std::string &     subsection = "SchlierenPostprocessor"); 

    void prepare(); 

    void compute_schlieren(const vector_type &U); 

    LinearAlgebra::distributed::Vector<double> schlieren; 

  private: 
    const MPI_Comm mpi_communicator; 
    TimerOutput &  computing_timer; 

    SmartPointer<const OfflineData<dim>> offline_data; 

    Vector<double> r; 

    unsigned int schlieren_index; 
    double       schlieren_beta; 
  }; 


  template <int dim> 
  class MainLoop : public ParameterAcceptor 
  { 
  public: 
    using vector_type = typename TimeStepping<dim>::vector_type; 

    MainLoop(const MPI_Comm mpi_communnicator); 

    void run(); 

  private: 
    vector_type interpolate_initial_values(const double t = 0); 

    void output(const vector_type &U, 
                const std::string &name, 
                double             t, 
                unsigned int       cycle, 
                bool               checkpoint = false); 

    const MPI_Comm     mpi_communicator; 
    std::ostringstream timer_output; 
    TimerOutput        computing_timer; 

    ConditionalOStream pcout; 

    std::string base_name; 
    double      t_final; 
    double      output_granularity; 

    bool asynchronous_writeback; 

    bool resume; 

    Discretization<dim>         discretization; 
    OfflineData<dim>            offline_data; 
    InitialValues<dim>          initial_values; 
    TimeStepping<dim>           time_stepping; 
    SchlierenPostprocessor<dim> schlieren_postprocessor; 

    vector_type output_vector; 

    std::future<void> background_thread_state; 
  }; 

  template <int dim> 
  Discretization<dim>::Discretization(const MPI_Comm     mpi_communicator, 
                                      TimerOutput &      computing_timer, 
                                      const std::string &subsection) 
    : ParameterAcceptor(subsection) 
    , mpi_communicator(mpi_communicator) 
    , triangulation(mpi_communicator) 
    , mapping(1) 
    , finite_element(1) 
    , quadrature(3) 
    , face_quadrature(3) 
    , computing_timer(computing_timer) 
  { 
    length = 4.; 
    add_parameter("length", length, "Length of computational domain"); 

    height = 2.; 
    add_parameter("height", height, "Height of computational domain"); 

    disk_position = 0.6; 
    add_parameter("object position", 
                  disk_position, 
                  "x position of immersed disk center point"); 

    disk_diameter = 0.5; 
    add_parameter("object diameter", 
                  disk_diameter, 
                  "Diameter of immersed disk"); 

    refinement = 5; 
    add_parameter("refinement", 
                  refinement, 
                  "Number of refinement steps of the geometry"); 
  } 



  template <int dim> 
  void Discretization<dim>::setup() 
  { 
    TimerOutput::Scope scope(computing_timer, "discretization - setup"); 

    triangulation.clear(); 

    Triangulation<dim> tria1, tria2, tria3, tria4, tria5, tria6; 

    GridGenerator::hyper_cube_with_cylindrical_hole( 
      tria1, disk_diameter / 2., disk_diameter, 0.5, 1, false); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria2, 
      {2, 1}, 
      Point<2>(-disk_diameter, disk_diameter), 
      Point<2>(disk_diameter, height / 2.)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria3, 
      {2, 1}, 
      Point<2>(-disk_diameter, -disk_diameter), 
      Point<2>(disk_diameter, -height / 2.)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria4, 
      {6, 2}, 
      Point<2>(disk_diameter, -disk_diameter), 
      Point<2>(length - disk_position, disk_diameter)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria5, 
      {6, 1}, 
      Point<2>(disk_diameter, disk_diameter), 
      Point<2>(length - disk_position, height / 2.)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria6, 
      {6, 1}, 
      Point<2>(disk_diameter, -height / 2.), 
      Point<2>(length - disk_position, -disk_diameter)); 

    GridGenerator::merge_triangulations( 
      {&tria1, &tria2, &tria3, &tria4, &tria5, &tria6}, 
      triangulation, 
      1.e-12, 
      true); 

    triangulation.set_manifold(0, PolarManifold<2>(Point<2>())); 


    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        for (const auto v : cell->vertex_indices()) 
          { 
            if (cell->vertex(v)[0] <= -disk_diameter + 1.e-6) 
              cell->vertex(v)[0] = -disk_position; 
          } 
      } 

    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        for (const auto f : cell->face_indices()) 
          { 
            const auto face = cell->face(f); 

            if (face->at_boundary()) 
              { 
                const auto center = face->center(); 

                if (center[0] > length - disk_position - 1.e-6) 
                  face->set_boundary_id(Boundaries::do_nothing); 
                else if (center[0] < -disk_position + 1.e-6) 
                  face->set_boundary_id(Boundaries::dirichlet); 
                else 
                  face->set_boundary_id(Boundaries::free_slip); 
              } 
          } 
      } 

    triangulation.refine_global(refinement); 
  } 


  template <int dim> 
  OfflineData<dim>::OfflineData(const MPI_Comm             mpi_communicator, 
                                TimerOutput &              computing_timer, 
                                const Discretization<dim> &discretization, 
                                const std::string &        subsection) 
    : ParameterAcceptor(subsection) 
    , dof_handler(discretization.triangulation) 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(computing_timer) 
    , discretization(&discretization) 
  {} 


  template <int dim> 
  void OfflineData<dim>::setup() 
  { 
    IndexSet locally_owned; 
    IndexSet locally_relevant; 

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "offline_data - distribute dofs"); 

      dof_handler.distribute_dofs(discretization->finite_element); 

      locally_owned   = dof_handler.locally_owned_dofs(); 
      n_locally_owned = locally_owned.n_elements(); 

      DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant); 
      n_locally_relevant = locally_relevant.n_elements(); 

      partitioner = 
        std::make_shared<Utilities::MPI::Partitioner>(locally_owned, 
                                                      locally_relevant, 
                                                      mpi_communicator); 
    } 




    { 
      TimerOutput::Scope scope( 
        computing_timer, 
        "offline_data - create sparsity pattern and set up matrices"); 



      DynamicSparsityPattern dsp(n_locally_relevant, n_locally_relevant); 

      const auto dofs_per_cell = 
        discretization->finite_element.n_dofs_per_cell(); 
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        { 
          if (cell->is_artificial()) 
            continue; 

          /* We transform the set of global dof indices on the cell to the
           * corresponding "local" index range on the MPI process: */
          cell->get_dof_indices(dof_indices); 
          std::transform(dof_indices.begin(), 
                         dof_indices.end(), 
                         dof_indices.begin(), 
                         [&](types::global_dof_index index) { 
                           return partitioner->global_to_local(index); 
                         }); 

/* 为每个dof简单地添加一个与所有其他 "本地 "的联接。 */

           /* dofs on the cell: */ 


          for (const auto dof : dof_indices) 
            dsp.add_entries(dof, dof_indices.begin(), dof_indices.end()); 
        } 

      sparsity_pattern.copy_from(dsp); 

      lumped_mass_matrix.reinit(sparsity_pattern); 
      norm_matrix.reinit(sparsity_pattern); 
      for (auto &matrix : cij_matrix) 
        matrix.reinit(sparsity_pattern); 
      for (auto &matrix : nij_matrix) 
        matrix.reinit(sparsity_pattern); 
    } 
  } 


  namespace 
  { 

    template <int dim> 
    struct CopyData 
    { 
      bool                                         is_artificial; 
      std::vector<types::global_dof_index>         local_dof_indices; 
      typename OfflineData<dim>::BoundaryNormalMap local_boundary_normal_map; 
      FullMatrix<double>                           cell_lumped_mass_matrix; 
      std::array<FullMatrix<double>, dim>          cell_cij_matrix; 
    }; 



    template <typename IteratorType> 
    DEAL_II_ALWAYS_INLINE inline SparseMatrix<double>::value_type 
    get_entry(const SparseMatrix<double> &matrix, const IteratorType &it) 
    { 
      const SparseMatrix<double>::const_iterator matrix_iterator( 
        &matrix, it->global_index()); 
      return matrix_iterator->value(); 
    } 


    template <typename IteratorType> 
    DEAL_II_ALWAYS_INLINE inline void 
    set_entry(SparseMatrix<double> &           matrix, 
              const IteratorType &             it, 
              SparseMatrix<double>::value_type value) 
    { 
      SparseMatrix<double>::iterator matrix_iterator(&matrix, 
                                                     it->global_index()); 
      matrix_iterator->value() = value; 
    } 

    template <std::size_t k, typename IteratorType> 
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k> 
    gather_get_entry(const std::array<SparseMatrix<double>, k> &c_ij, 
                     const IteratorType                         it) 
    { 
      Tensor<1, k> result; 
      for (unsigned int j = 0; j < k; ++j) 
        result[j] = get_entry(c_ij[j], it); 
      return result; 
    } 


    template <std::size_t k> 
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k> 
    gather(const std::array<SparseMatrix<double>, k> &n_ij, 
           const unsigned int                         i, 
           const unsigned int                         j) 
    { 
      Tensor<1, k> result; 
      for (unsigned int l = 0; l < k; ++l) 
        result[l] = n_ij[l](i, j); 
      return result; 
    } 

    template <std::size_t k> 
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k> 
    gather(const std::array<LinearAlgebra::distributed::Vector<double>, k> &U, 
           const unsigned int                                               i) 
    { 
      Tensor<1, k> result; 
      for (unsigned int j = 0; j < k; ++j) 
        result[j] = U[j].local_element(i); 
      return result; 
    } 

    template <std::size_t k, int k2> 
    DEAL_II_ALWAYS_INLINE inline void 
    scatter(std::array<LinearAlgebra::distributed::Vector<double>, k> &U, 
            const Tensor<1, k2> &                                      tensor, 
            const unsigned int                                         i) 
    { 
      static_assert(k == k2, 
                    "The dimensions of the input arguments must agree"); 
      for (unsigned int j = 0; j < k; ++j) 
        U[j].local_element(i) = tensor[j]; 
    } 
  } // namespace 











  template <int dim> 
  void OfflineData<dim>::assemble() 
  { 
    lumped_mass_matrix = 0.; 
    norm_matrix        = 0.; 
    for (auto &matrix : cij_matrix) 
      matrix = 0.; 
    for (auto &matrix : nij_matrix) 
      matrix = 0.; 

    unsigned int dofs_per_cell = 
      discretization->finite_element.n_dofs_per_cell(); 
    unsigned int n_q_points = discretization->quadrature.size(); 


    MeshWorker::ScratchData<dim> scratch_data( 
      discretization->mapping, 
      discretization->finite_element, 
      discretization->quadrature, 
      update_values | update_gradients | update_quadrature_points | 
        update_JxW_values, 
      discretization->face_quadrature, 
      update_normal_vectors | update_values | update_JxW_values); 

    { 
      TimerOutput::Scope scope( 
        computing_timer, 
        "offline_data - assemble lumped mass matrix, and c_ij"); 

      const auto local_assemble_system = // 
        [&](const typename DoFHandler<dim>::cell_iterator &cell, 
            MeshWorker::ScratchData<dim> &                 scratch, 
            CopyData<dim> &                                copy) { 
          copy.is_artificial = cell->is_artificial(); 
          if (copy.is_artificial) 
            return; 

          copy.local_boundary_normal_map.clear(); 
          copy.cell_lumped_mass_matrix.reinit(dofs_per_cell, dofs_per_cell); 
          for (auto &matrix : copy.cell_cij_matrix) 
            matrix.reinit(dofs_per_cell, dofs_per_cell); 

          const auto &fe_values = scratch.reinit(cell); 

          copy.local_dof_indices.resize(dofs_per_cell); 
          cell->get_dof_indices(copy.local_dof_indices); 

          std::transform(copy.local_dof_indices.begin(), 
                         copy.local_dof_indices.end(), 
                         copy.local_dof_indices.begin(), 
                         [&](types::global_dof_index index) { 
                           return partitioner->global_to_local(index); 
                         }); 


          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            { 
              const auto JxW = fe_values.JxW(q_point); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const auto value_JxW = 
                    fe_values.shape_value(j, q_point) * JxW; 
                  const auto grad_JxW = fe_values.shape_grad(j, q_point) * JxW; 

                  copy.cell_lumped_mass_matrix(j, j) += value_JxW; 

                  for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                    { 
                      const auto value = fe_values.shape_value(i, q_point); 
                      for (unsigned int d = 0; d < dim; ++d) 
                        copy.cell_cij_matrix[d](i, j) += value * grad_JxW[d]; 

                    } /* i */ 


                }     /* j */ 


            }         /* q */ 




          for (const auto f : cell->face_indices()) 
            { 
              const auto face = cell->face(f); 
              const auto id   = face->boundary_id(); 

              if (!face->at_boundary()) 
                continue; 

              const auto &fe_face_values = scratch.reinit(cell, f); 

              const unsigned int n_face_q_points = 
                fe_face_values.get_quadrature().size(); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  if (!discretization->finite_element.has_support_on_face(j, f)) 
                    continue; 


                  Tensor<1, dim> normal; 
                  if (id == Boundaries::free_slip) 
                    { 
                      for (unsigned int q = 0; q < n_face_q_points; ++q) 
                        normal += fe_face_values.normal_vector(q) * 
                                  fe_face_values.shape_value(j, q); 
                    } 

                  const auto index = copy.local_dof_indices[j]; 

                  Point<dim> position; 
                  for (const auto v : cell->vertex_indices()) 
                    if (cell->vertex_dof_index(v, 0) == 
                        partitioner->local_to_global(index)) 
                      { 
                        position = cell->vertex(v); 
                        break; 
                      } 

                  const auto old_id = 
                    std::get<1>(copy.local_boundary_normal_map[index]); 
                  copy.local_boundary_normal_map[index] = 
                    std::make_tuple(normal, std::max(old_id, id), position); 
                } 
            } 
        }; 


      const auto copy_local_to_global = [&](const CopyData<dim> &copy) { 
        if (copy.is_artificial) 
          return; 

        for (const auto &it : copy.local_boundary_normal_map) 
          { 
            std::get<0>(boundary_normal_map[it.first]) += 
              std::get<0>(it.second); 
            std::get<1>(boundary_normal_map[it.first]) = 
              std::max(std::get<1>(boundary_normal_map[it.first]), 
                       std::get<1>(it.second)); 
            std::get<2>(boundary_normal_map[it.first]) = std::get<2>(it.second); 
          } 

        lumped_mass_matrix.add(copy.local_dof_indices, 
                               copy.cell_lumped_mass_matrix); 

        for (int k = 0; k < dim; ++k) 
          { 
            cij_matrix[k].add(copy.local_dof_indices, copy.cell_cij_matrix[k]); 
            nij_matrix[k].add(copy.local_dof_indices, copy.cell_cij_matrix[k]); 
          } 
      }; 

      WorkStream::run(dof_handler.begin_active(), 
                      dof_handler.end(), 
                      local_assemble_system, 
                      copy_local_to_global, 
                      scratch_data, 
                      CopyData<dim>()); 
    } 

















    { 
      TimerOutput::Scope scope(computing_timer, 
                               "offline_data - compute |c_ij|, and n_ij"); 

      const std_cxx20::ranges::iota_view<unsigned int, unsigned int> indices( 
        0, n_locally_relevant); 

      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto row_index : 
               std_cxx20::ranges::iota_view<unsigned int, unsigned int>(*i1, 
                                                                        *i2)) 
            { 


              std::for_each( 
                sparsity_pattern.begin(row_index), 
                sparsity_pattern.end(row_index), 
                [&](const dealii::SparsityPatternIterators::Accessor &jt) { 
                  const auto   c_ij = gather_get_entry(cij_matrix, &jt); 
                  const double norm = c_ij.norm(); 

                  set_entry(norm_matrix, &jt, norm); 
                  for (unsigned int j = 0; j < dim; ++j) 
                    set_entry(nij_matrix[j], &jt, c_ij[j] / norm); 
                }); 
            } 
        }; 

      parallel::apply_to_subranges(indices.begin(), 
                                   indices.end(), 
                                   on_subranges, 
                                   4096); 


      for (auto &it : boundary_normal_map) 
        { 
          auto &normal = std::get<0>(it.second); 
          normal /= (normal.norm() + std::numeric_limits<double>::epsilon()); 
        } 
    } 
  } 





  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline Tensor<1, dim> 
  ProblemDescription<dim>::momentum(const state_type &U) 
  { 
    Tensor<1, dim> result; 
    std::copy_n(&U[1], dim, &result[0]); 
    return result; 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::internal_energy(const state_type &U) 
  { 
    const double &rho = U[0]; 
    const auto    m   = momentum(U); 
    const double &E   = U[dim + 1]; 
    return E - 0.5 * m.norm_square() / rho; 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::pressure(const state_type &U) 
  { 
    return (gamma - 1.) * internal_energy(U); 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::speed_of_sound(const state_type &U) 
  { 
    const double &rho = U[0]; 
    const double  p   = pressure(U); 

    return std::sqrt(gamma * p / rho); 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline typename ProblemDescription<dim>::flux_type 
  ProblemDescription<dim>::flux(const state_type &U) 
  { 
    const double &rho = U[0]; 
    const auto    m   = momentum(U); 
    const auto    p   = pressure(U); 
    const double &E   = U[dim + 1]; 

    flux_type result; 

    result[0] = m; 
    for (unsigned int i = 0; i < dim; ++i) 
      { 
        result[1 + i] = m * m[i] / rho; 
        result[1 + i][i] += p; 
      } 
    result[dim + 1] = m / rho * (E + p); 

    return result; 
  } 






  namespace 
  { 
    template <int dim> 
    DEAL_II_ALWAYS_INLINE inline std::array<double, 4> riemann_data_from_state( 
      const typename ProblemDescription<dim>::state_type U, 
      const Tensor<1, dim> &                             n_ij) 
    { 
      Tensor<1, 3> projected_U; 
      projected_U[0] = U[0]; 


      const auto m   = ProblemDescription<dim>::momentum(U); 
      projected_U[1] = n_ij * m; 

      const auto perpendicular_m = m - projected_U[1] * n_ij; 
      projected_U[2] = U[1 + dim] - 0.5 * perpendicular_m.norm_square() / U[0]; 


      return {{projected_U[0], 
               projected_U[1] / projected_U[0], 
               ProblemDescription<1>::pressure(projected_U), 
               ProblemDescription<1>::speed_of_sound(projected_U)}}; 
    } 


    DEAL_II_ALWAYS_INLINE inline double positive_part(const double number) 
    { 
      return std::max(number, 0.); 
    } 

    DEAL_II_ALWAYS_INLINE inline double negative_part(const double number) 
    { 
      return -std::min(number, 0.); 
    } 


    DEAL_II_ALWAYS_INLINE inline double 
    lambda1_minus(const std::array<double, 4> &riemann_data, 
                  const double                 p_star) 
    { 

      /* Implements formula (3.7) in Guermond-Popov-2016 */


      constexpr double gamma = ProblemDescription<1>::gamma; 
      const auto       u     = riemann_data[1]; 
      const auto       p     = riemann_data[2]; 
      const auto       a     = riemann_data[3]; 

      const double factor = (gamma + 1.0) / 2.0 / gamma; 
      const double tmp    = positive_part((p_star - p) / p); 
      return u - a * std::sqrt(1.0 + factor * tmp); 
    } 


    DEAL_II_ALWAYS_INLINE inline double 
    lambda3_plus(const std::array<double, 4> &riemann_data, const double p_star) 
    { 

      /* Implements formula (3.8) in Guermond-Popov-2016 */

      constexpr double gamma = ProblemDescription<1>::gamma; 
      const auto       u     = riemann_data[1]; 
      const auto       p     = riemann_data[2]; 
      const auto       a     = riemann_data[3]; 

      const double factor = (gamma + 1.0) / 2.0 / gamma; 
      const double tmp    = positive_part((p_star - p) / p); 
      return u + a * std::sqrt(1.0 + factor * tmp); 
    } 


    DEAL_II_ALWAYS_INLINE inline double 
    lambda_max_two_rarefaction(const std::array<double, 4> &riemann_data_i, 
                               const std::array<double, 4> &riemann_data_j) 
    { 
      constexpr double gamma = ProblemDescription<1>::gamma; 
      const auto       u_i   = riemann_data_i[1]; 
      const auto       p_i   = riemann_data_i[2]; 
      const auto       a_i   = riemann_data_i[3]; 
      const auto       u_j   = riemann_data_j[1]; 
      const auto       p_j   = riemann_data_j[2]; 
      const auto       a_j   = riemann_data_j[3]; 

      const double numerator = a_i + a_j - (gamma - 1.) / 2. * (u_j - u_i); 

      const double denominator = 
        a_i * std::pow(p_i / p_j, -1. * (gamma - 1.) / 2. / gamma) + a_j * 1.; 

/* Guermond-Popov-2016 */ 

中的公式（4.3）。
      const double p_star = 
        p_j * std::pow(numerator / denominator, 2. * gamma / (gamma - 1)); 

      const double lambda1 = lambda1_minus(riemann_data_i, p_star); 
      const double lambda3 = lambda3_plus(riemann_data_j, p_star); 

/* Guermond-Popov-2016中的公式（2.11）  */ 

      return std::max(positive_part(lambda3), negative_part(lambda1)); 
    } 


    DEAL_II_ALWAYS_INLINE inline double 
    lambda_max_expansion(const std::array<double, 4> &riemann_data_i, 
                         const std::array<double, 4> &riemann_data_j) 
    { 
      const auto u_i = riemann_data_i[1]; 
      const auto a_i = riemann_data_i[3]; 
      const auto u_j = riemann_data_j[1]; 
      const auto a_j = riemann_data_j[3]; 

      return std::max(std::abs(u_i), std::abs(u_j)) + 5. * std::max(a_i, a_j); 
    } 
  } // namespace 


  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::compute_lambda_max(const state_type &    U_i, 
                                              const state_type &    U_j, 
                                              const Tensor<1, dim> &n_ij) 
  { 
    const auto riemann_data_i = riemann_data_from_state(U_i, n_ij); 
    const auto riemann_data_j = riemann_data_from_state(U_j, n_ij); 

    const double lambda_1 = 
      lambda_max_two_rarefaction(riemann_data_i, riemann_data_j); 

    const double lambda_2 = 
      lambda_max_expansion(riemann_data_i, riemann_data_j); 

    return std::min(lambda_1, lambda_2); 
  } 


  template <> 
  const std::array<std::string, 3> ProblemDescription<1>::component_names{ 
    {"rho", "m", "E"}}; 

  template <> 
  const std::array<std::string, 4> ProblemDescription<2>::component_names{ 
    {"rho", "m_1", "m_2", "E"}}; 

  template <> 
  const std::array<std::string, 5> ProblemDescription<3>::component_names{ 
    {"rho", "m_1", "m_2", "m_3", "E"}}; 




  template <int dim> 
  InitialValues<dim>::InitialValues(const std::string &subsection) 
    : ParameterAcceptor(subsection) 
  { 

    /* We wire up the slot InitialValues<dim>::parse_parameters_callback to
       the ParameterAcceptor::parse_parameters_call_back signal: */


    ParameterAcceptor::parse_parameters_call_back.connect( 
      std::bind(&InitialValues<dim>::parse_parameters_callback, this)); 

    initial_direction[0] = 1.; 
    add_parameter("initial direction", 
                  initial_direction, 
                  "Initial direction of the uniform flow field"); 

    initial_1d_state[0] = ProblemDescription<dim>::gamma; 
    initial_1d_state[1] = 3.; 
    initial_1d_state[2] = 1.; 
    add_parameter("initial 1d state", 
                  initial_1d_state, 
                  "Initial 1d state (rho, u, p) of the uniform flow field"); 
  } 



  template <int dim> 
  void InitialValues<dim>::parse_parameters_callback() 
  { 
    AssertThrow(initial_direction.norm() != 0., 
                ExcMessage( 
                  "Initial shock front direction is set to the zero vector.")); 
    initial_direction /= initial_direction.norm(); 


    initial_state = [this](const Point<dim> & /*point*/, double /*t*/) { 
      const double            rho   = initial_1d_state[0]; 
      const double            u     = initial_1d_state[1]; 
      const double            p     = initial_1d_state[2]; 
      static constexpr double gamma = ProblemDescription<dim>::gamma; 

      state_type state; 

      state[0] = rho; 
      for (unsigned int i = 0; i < dim; ++i) 
        state[1 + i] = rho * u * initial_direction[i]; 

      state[dim + 1] = p / (gamma - 1.) + 0.5 * rho * u * u; 

      return state; 
    }; 
  } 


  template <int dim> 
  TimeStepping<dim>::TimeStepping( 
    const MPI_Comm            mpi_communicator, 
    TimerOutput &             computing_timer, 
    const OfflineData<dim> &  offline_data, 
    const InitialValues<dim> &initial_values, 
    const std::string &       subsection /*= "TimeStepping"*/) 
    : ParameterAcceptor(subsection) 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(computing_timer) 
    , offline_data(&offline_data) 
    , initial_values(&initial_values) 
  { 
    cfl_update = 0.80; 
    add_parameter("cfl update", 
                  cfl_update, 
                  "Relative CFL constant used for update"); 
  } 


  template <int dim> 
  void TimeStepping<dim>::prepare() 
  { 
    TimerOutput::Scope scope(computing_timer, 
                             "time_stepping - prepare scratch space"); 

    for (auto &it : temporary_vector) 
      it.reinit(offline_data->partitioner); 

    dij_matrix.reinit(offline_data->sparsity_pattern); 
  } 


  template <int dim> 
  double TimeStepping<dim>::make_one_step(vector_type &U, double t) 
  { 
    const auto &n_locally_owned    = offline_data->n_locally_owned; 
    const auto &n_locally_relevant = offline_data->n_locally_relevant; 

    const std_cxx20::ranges::iota_view<unsigned int, unsigned int> 
      indices_owned(0, n_locally_owned); 
    const std_cxx20::ranges::iota_view<unsigned int, unsigned int> 
      indices_relevant(0, n_locally_relevant); 

    const auto &sparsity = offline_data->sparsity_pattern; 

    const auto &lumped_mass_matrix = offline_data->lumped_mass_matrix; 
    const auto &norm_matrix        = offline_data->norm_matrix; 
    const auto &nij_matrix         = offline_data->nij_matrix; 
    const auto &cij_matrix         = offline_data->cij_matrix; 

    const auto &boundary_normal_map = offline_data->boundary_normal_map; 






    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 1 compute d_ij"); 

      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto i : 
               std_cxx20::ranges::iota_view<unsigned int, unsigned int>(*i1, 
                                                                        *i2)) 
            { 
              const auto U_i = gather(U, i); 


              for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt) 
                { 
                  const auto j = jt->column(); 


                  if (j >= i) 
                    continue; 

                  const auto U_j = gather(U, j); 

                  const auto   n_ij = gather_get_entry(nij_matrix, jt); 
                  const double norm = get_entry(norm_matrix, jt); 

                  const auto lambda_max = 
                    ProblemDescription<dim>::compute_lambda_max(U_i, U_j, n_ij); 

                  double d = norm * lambda_max; 


                  if (boundary_normal_map.count(i) != 0 && 
                      boundary_normal_map.count(j) != 0) 
                    { 
                      const auto n_ji = gather(nij_matrix, j, i); 
                      const auto lambda_max_2 = 
                        ProblemDescription<dim>::compute_lambda_max(U_j, 
                                                                    U_i, 
                                                                    n_ji); 
                      const double norm_2 = norm_matrix(j, i); 

                      d = std::max(d, norm_2 * lambda_max_2); 
                    } 

                  set_entry(dij_matrix, jt, d); 
                  dij_matrix(j, i) = d; 
                } 
            } 
        }; 

      parallel::apply_to_subranges(indices_relevant.begin(), 
                                   indices_relevant.end(), 
                                   on_subranges, 
                                   4096); 
    } 




    std::atomic<double> tau_max{std::numeric_limits<double>::infinity()}; 

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 2 compute d_ii, and tau_max"); 


      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          double tau_max_on_subrange = std::numeric_limits<double>::infinity(); 

          for (const auto i : 
               std_cxx20::ranges::iota_view<unsigned int, unsigned int>(*i1, 
                                                                        *i2)) 
            { 
              double d_sum = 0.; 

              for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt) 
                { 
                  const auto j = jt->column(); 

                  if (j == i) 
                    continue; 

                  d_sum -= get_entry(dij_matrix, jt); 
                } 


              dij_matrix.diag_element(i) = d_sum; 


              const double mass   = lumped_mass_matrix.diag_element(i); 
              const double tau    = cfl_update * mass / (-2. * d_sum); 
              tau_max_on_subrange = std::min(tau_max_on_subrange, tau); 
            } 

          double current_tau_max = tau_max.load(); 
          while (current_tau_max > tau_max_on_subrange && 
                 !tau_max.compare_exchange_weak(current_tau_max, 
                                                tau_max_on_subrange)) 
            ; 
        }; 

      parallel::apply_to_subranges(indices_relevant.begin(), 
                                   indices_relevant.end(), 
                                   on_subranges, 
                                   4096); 


      tau_max.store(Utilities::MPI::min(tau_max.load(), mpi_communicator)); 


      AssertThrow( 
        !std::isnan(tau_max.load()) && !std::isinf(tau_max.load()) && 
          tau_max.load() > 0., 
        ExcMessage( 
          "I'm sorry, Dave. I'm afraid I can't do that. - We crashed.")); 
    } 


   \mathbf{U}_i^{n+1} = \mathbf{U}_i^{n} - \frac{\tau_{\text{max}} }{m_i}
   \sum_{j \in \mathcal{I}(i)} (\mathbb{f}(\mathbf{U}_j^{n}) -
   \mathbb{f}(\mathbf{U}_i^{n})) \cdot \mathbf{c}_{ij} - d_{ij}
   (\mathbf{U}_j^{n} - \mathbf{U}_i^{n})
 \f]


    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 3 perform update"); 

      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto i : boost::make_iterator_range(i1, i2)) 
            { 
              Assert(i < n_locally_owned, ExcInternalError()); 

              const auto U_i = gather(U, i); 

              const auto   f_i = ProblemDescription<dim>::flux(U_i); 
              const double m_i = lumped_mass_matrix.diag_element(i); 

              auto U_i_new = U_i; 

              for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt) 
                { 
                  const auto j = jt->column(); 

                  const auto U_j = gather(U, j); 
                  const auto f_j = ProblemDescription<dim>::flux(U_j); 

                  const auto c_ij = gather_get_entry(cij_matrix, jt); 
                  const auto d_ij = get_entry(dij_matrix, jt); 

                  for (unsigned int k = 0; k < problem_dimension; ++k) 
                    { 
                      U_i_new[k] += 
                        tau_max / m_i * 
                        (-(f_j[k] - f_i[k]) * c_ij + d_ij * (U_j[k] - U_i[k])); 
                    } 
                } 

              scatter(temporary_vector, U_i_new, i); 
            } 
        }; 

      parallel::apply_to_subranges(indices_owned.begin(), 
                                   indices_owned.end(), 
                                   on_subranges, 
                                   4096); 
    } 





    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 4 fix boundary states"); 

      for (auto it : boundary_normal_map) 
        { 
          const auto i = it.first; 


          if (i >= n_locally_owned) 
            continue; 

          const auto &normal   = std::get<0>(it.second); 
          const auto &id       = std::get<1>(it.second); 
          const auto &position = std::get<2>(it.second); 

          auto U_i = gather(temporary_vector, i); 


          if (id == Boundaries::free_slip) 
            { 
              auto m = ProblemDescription<dim>::momentum(U_i); 
              m -= (m * normal) * normal; 
              for (unsigned int k = 0; k < dim; ++k) 
                U_i[k + 1] = m[k]; 
            } 


          else if (id == Boundaries::dirichlet) 
            { 
              U_i = initial_values->initial_state(position, t + tau_max); 
            } 

          scatter(temporary_vector, U_i, i); 
        } 
    } 

    for (auto &it : temporary_vector) 
      it.update_ghost_values(); 

    U.swap(temporary_vector); 

    return tau_max; 
  } 




  template <int dim> 
  SchlierenPostprocessor<dim>::SchlierenPostprocessor( 
    const MPI_Comm          mpi_communicator, 
    TimerOutput &           computing_timer, 
    const OfflineData<dim> &offline_data, 
    const std::string &     subsection /*= "SchlierenPostprocessor"*/) 
    : ParameterAcceptor(subsection) 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(computing_timer) 
    , offline_data(&offline_data) 
  { 
    schlieren_beta = 10.; 
    add_parameter("schlieren beta", 
                  schlieren_beta, 
                  "Beta factor used in Schlieren-type postprocessor"); 

    schlieren_index = 0; 
    add_parameter("schlieren index", 
                  schlieren_index, 
                  "Use the corresponding component of the state vector for the " 
                  "schlieren plot"); 
  } 


  template <int dim> 
  void SchlierenPostprocessor<dim>::prepare() 
  { 
    TimerOutput::Scope scope(computing_timer, 
                             "schlieren_postprocessor - prepare scratch space"); 

    r.reinit(offline_data->n_locally_relevant); 
    schlieren.reinit(offline_data->partitioner); 
  } 














  template <int dim> 
  void SchlierenPostprocessor<dim>::compute_schlieren(const vector_type &U) 
  { 
    TimerOutput::Scope scope( 
      computing_timer, "schlieren_postprocessor - compute schlieren plot"); 

    const auto &sparsity            = offline_data->sparsity_pattern; 
    const auto &lumped_mass_matrix  = offline_data->lumped_mass_matrix; 
    const auto &cij_matrix          = offline_data->cij_matrix; 
    const auto &boundary_normal_map = offline_data->boundary_normal_map; 
    const auto &n_locally_owned     = offline_data->n_locally_owned; 

    const auto indices = 
      std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0, 
                                                               n_locally_owned); 


    std::atomic<double> r_i_max{0.}; 
    std::atomic<double> r_i_min{std::numeric_limits<double>::infinity()}; 


    { 
      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          double r_i_max_on_subrange = 0.; 
          double r_i_min_on_subrange = std::numeric_limits<double>::infinity(); 

          for (const auto i : boost::make_iterator_range(i1, i2)) 
            { 
              Assert(i < n_locally_owned, ExcInternalError()); 

              Tensor<1, dim> r_i; 

 
                { 
                  const auto j = jt->column(); 

                  if (i == j) 
                    continue; 

                  const auto U_js = U[schlieren_index].local_element(j); 
                  const auto c_ij = gather_get_entry(cij_matrix, jt); 
                  r_i += c_ij * U_js; 
                } 


              const auto bnm_it = boundary_normal_map.find(i); 
              if (bnm_it != boundary_normal_map.end()) 
                { 
                  const auto &normal = std::get<0>(bnm_it->second); 
                  const auto &id     = std::get<1>(bnm_it->second); 

                  if (id == Boundaries::free_slip) 
                    r_i -= 1. * (r_i * normal) * normal; 
                  else 
                    r_i = 0.; 
                } 


              const double m_i    = lumped_mass_matrix.diag_element(i); 
              r[i]                = r_i.norm() / m_i; 
              r_i_max_on_subrange = std::max(r_i_max_on_subrange, r[i]); 
              r_i_min_on_subrange = std::min(r_i_min_on_subrange, r[i]); 
            } 


          double current_r_i_max = r_i_max.load(); 
          while (current_r_i_max < r_i_max_on_subrange && 
                 !r_i_max.compare_exchange_weak(current_r_i_max, 
                                                r_i_max_on_subrange)) 
            ; 

          double current_r_i_min = r_i_min.load(); 
          while (current_r_i_min > r_i_min_on_subrange && 
                 !r_i_min.compare_exchange_weak(current_r_i_min, 
                                                r_i_min_on_subrange)) 
            ; 
        }; 

      parallel::apply_to_subranges(indices.begin(), 
                                   indices.end(), 
                                   on_subranges, 
                                   4096); 
    } 


    r_i_max.store(Utilities::MPI::max(r_i_max.load(), mpi_communicator)); 
    r_i_min.store(Utilities::MPI::min(r_i_min.load(), mpi_communicator)); 


    { 
      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto i : boost::make_iterator_range(i1, i2)) 
            { 
              Assert(i < n_locally_owned, ExcInternalError()); 

              schlieren.local_element(i) = 
                1. - std::exp(-schlieren_beta * (r[i] - r_i_min) / 
                              (r_i_max - r_i_min)); 
            } 
        }; 

      parallel::apply_to_subranges(indices.begin(), 
                                   indices.end(), 
                                   on_subranges, 
                                   4096); 
    } 


    schlieren.update_ghost_values(); 
  } 



  template <int dim> 
  MainLoop<dim>::MainLoop(const MPI_Comm mpi_communicator) 
    : ParameterAcceptor("A - MainLoop") 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(mpi_communicator, 
                      timer_output, 
                      TimerOutput::never, 
                      TimerOutput::cpu_and_wall_times) 
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
    , discretization(mpi_communicator, computing_timer, "B - Discretization") 
    , offline_data(mpi_communicator, 
                   computing_timer, 
                   discretization, 
                   "C - OfflineData") 
    , initial_values("D - InitialValues") 
    , time_stepping(mpi_communicator, 
                    computing_timer, 
                    offline_data, 
                    initial_values, 
                    "E - TimeStepping") 
    , schlieren_postprocessor(mpi_communicator, 
                              computing_timer, 
                              offline_data, 
                              "F - SchlierenPostprocessor") 
  { 
    base_name = "test"; 
    add_parameter("basename", base_name, "Base name for all output files"); 

    t_final = 4.; 
    add_parameter("final time", t_final, "Final time"); 

    output_granularity = 0.02; 
    add_parameter("output granularity", 
                  output_granularity, 
                  "time interval for output"); 

    asynchronous_writeback = true; 
    add_parameter("asynchronous writeback", 
                  asynchronous_writeback, 
                  "Write out solution in a background thread performing IO"); 

    resume = false; 
    add_parameter("resume", resume, "Resume an interrupted computation."); 
  } 


  namespace 
  { 
    void print_head(ConditionalOStream &pcout, 
                    const std::string & header, 
                    const std::string & secondary = "") 
    { 
      const auto header_size   = header.size(); 
      const auto padded_header = std::string((34 - header_size) / 2, ' ') + 
                                 header + 
                                 std::string((35 - header_size) / 2, ' '); 

      const auto secondary_size = secondary.size(); 
      const auto padded_secondary = 
        std::string((34 - secondary_size) / 2, ' ') + secondary + 
        std::string((35 - secondary_size) / 2, ' '); 

 /* 关闭clang-format  */ 

      pcout << std::endl; 
      pcout << "    ####################################################" << std::endl; 
      pcout << "    #########                                  #########" << std::endl; 
      pcout << "    #########"     <<  padded_header   <<     "#########" << std::endl; 
      pcout << "    #########"     << padded_secondary <<     "#########" << std::endl; 
      pcout << "    #########                                  #########" << std::endl; 
      pcout << "    ####################################################" << std::endl; 
      pcout << std::endl; 
    /* clang-format on  */ 
    } 
  } // namespace 


  template <int dim> 
  void MainLoop<dim>::run() 
  { 

    pcout << "Reading parameters and allocating objects... " << std::flush; 

    ParameterAcceptor::initialize("step-69.prm"); 
    pcout << "done" << std::endl; 


    { 
      print_head(pcout, "create triangulation"); 
      discretization.setup(); 

      pcout << "Number of active cells:       " 
            << discretization.triangulation.n_global_active_cells() 
            << std::endl; 

      print_head(pcout, "compute offline data"); 
      offline_data.setup(); 
      offline_data.assemble(); 

      pcout << "Number of degrees of freedom: " 
            << offline_data.dof_handler.n_dofs() << std::endl; 

      print_head(pcout, "set up time step"); 
      time_stepping.prepare(); 
      schlieren_postprocessor.prepare(); 
    } 


    double       t            = 0.; 
    unsigned int output_cycle = 0; 

    print_head(pcout, "interpolate initial values"); 
    vector_type U = interpolate_initial_values(); 


    if (resume) 
      { 
        print_head(pcout, "restore interrupted computation"); 

        const unsigned int i = 
          discretization.triangulation.locally_owned_subdomain(); 

        const std::string name = base_name + "-checkpoint-" + 
                                 Utilities::int_to_string(i, 4) + ".archive"; 
        std::ifstream file(name, std::ios::binary); 


        boost::archive::binary_iarchive ia(file); 
        ia >> t >> output_cycle; 

        for (auto &it1 : U) 
          { 

            for (auto &it2 : it1) 
              ia >> it2; 
            it1.update_ghost_values(); 
          } 
      } 


    output(U, base_name, t, output_cycle++); 

    print_head(pcout, "enter main loop"); 

    for (unsigned int cycle = 1; t < t_final; ++cycle) 
      { 


        std::ostringstream head; 
        std::ostringstream secondary; 

        head << "Cycle  " << Utilities::int_to_string(cycle, 6) << "  (" // 
             << std::fixed << std::setprecision(1) << t / t_final * 100  // 
             << "%)"; 
        secondary << "at time t = " << std::setprecision(8) << std::fixed << t; 

        print_head(pcout, head.str(), secondary.str()); 


        t += time_stepping.make_one_step(U, t); 



        if (t > output_cycle * output_granularity) 
          { 
            output(U, base_name, t, output_cycle, true); 
            ++output_cycle; 
          } 
      } 


    if (background_thread_state.valid()) 
      background_thread_state.wait(); 

    computing_timer.print_summary(); 
    pcout << timer_output.str() << std::endl; 
  } 


  template <int dim> 
  typename MainLoop<dim>::vector_type 
  MainLoop<dim>::interpolate_initial_values(const double t) 
  { 
    pcout << "MainLoop<dim>::interpolate_initial_values(t = " << t << ")" 
          << std::endl; 
    TimerOutput::Scope scope(computing_timer, 
                             "main_loop - setup scratch space"); 

    vector_type U; 

    for (auto &it : U) 
      it.reinit(offline_data.partitioner); 

    constexpr auto problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 


    for (unsigned int i = 0; i < problem_dimension; ++i) 
      VectorTools::interpolate(offline_data.dof_handler, 
                               ScalarFunctionFromFunctionObject<dim, double>( 
                                 [&](const Point<dim> &x) { 
                                   return initial_values.initial_state(x, t)[i]; 
                                 }), 
                               U[i]); 

    for (auto &it : U) 
      it.update_ghost_values(); 

    return U; 
  } 




  template <int dim> 
  void MainLoop<dim>::output(const typename MainLoop<dim>::vector_type &U, 
                             const std::string &                        name, 
                             const double                               t, 
                             const unsigned int                         cycle, 
                             const bool checkpoint) 
  { 
    pcout << "MainLoop<dim>::output(t = " << t 
          << ", checkpoint = " << checkpoint << ")" << std::endl; 


    if (background_thread_state.valid()) 
      { 
        TimerOutput::Scope timer(computing_timer, "main_loop - stalled output"); 
        background_thread_state.wait(); 
      } 

    constexpr auto problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 


    for (unsigned int i = 0; i < problem_dimension; ++i) 
      { 
        output_vector[i] = U[i]; 
        output_vector[i].update_ghost_values(); 
      } 

    schlieren_postprocessor.compute_schlieren(output_vector); 

    auto data_out = std::make_shared<DataOut<dim>>(); 

    data_out->attach_dof_handler(offline_data.dof_handler); 

    const auto &component_names = ProblemDescription<dim>::component_names; 

    for (unsigned int i = 0; i < problem_dimension; ++i) 
      data_out->add_data_vector(output_vector[i], component_names[i]); 

    data_out->add_data_vector(schlieren_postprocessor.schlieren, 
                              "schlieren_plot"); 

    data_out->build_patches(discretization.mapping, 
                            discretization.finite_element.degree - 1); 


    const auto output_worker = [this, name, t, cycle, checkpoint, data_out]() { 
      if (checkpoint) 
        { 


          const unsigned int i = 
            discretization.triangulation.locally_owned_subdomain(); 
          std::string filename = 
            name + "-checkpoint-" + Utilities::int_to_string(i, 4) + ".archive"; 

          std::ofstream file(filename, std::ios::binary | std::ios::trunc); 

          boost::archive::binary_oarchive oa(file); 
          oa << t << cycle; 
          for (const auto &it1 : output_vector) 
            for (const auto &it2 : it1) 
              oa << it2; 
        } 

      DataOutBase::VtkFlags flags(t, 
                                  cycle, 
                                  true, 
                                  DataOutBase::VtkFlags::best_speed); 
      data_out->set_flags(flags); 

      data_out->write_vtu_with_pvtu_record( 
        "", name + "-solution", cycle, mpi_communicator, 6); 
    }; 



    if (asynchronous_writeback) 
      { 
        background_thread_state = std::async(std::launch::async, output_worker); 
      } 
    else 
      { 
        output_worker(); 
      } 
  } 

} // namespace Step69 


int main(int argc, char *argv[]) 
{ 
  try 
    { 
      constexpr int dim = 2; 

      using namespace dealii; 
      using namespace Step69; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv); 

      MPI_Comm      mpi_communicator(MPI_COMM_WORLD); 
      MainLoop<dim> main_loop(mpi_communicator); 

      main_loop.run(); 
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
} 


