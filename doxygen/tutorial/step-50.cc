

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
 * Author: Thomas C. Clevenger, Clemson University 
 *         Timo Heister, Clemson University 
 *         Guido Kanschat, Heidelberg University 
 *         Martin Kronbichler, Technical University of Munich 
 */ 




#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/data_out_base.h> 
#include <deal.II/base/index_set.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/distributed/grid_refinement.h> 
#include <deal.II/distributed/tria.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 


#include <deal.II/lac/generic_linear_algebra.h> 


#define FORCE_USE_OF_TRILINOS 

namespace LA 
{ 
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \ 
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS)) 
  using namespace dealii::LinearAlgebraPETSc; 
#  define USE_PETSC_LA 
#elif defined(DEAL_II_WITH_TRILINOS) 
  using namespace dealii::LinearAlgebraTrilinos; 
#else 
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required 
#endif 
} // namespace LA 

#include <deal.II/matrix_free/matrix_free.h> 
#include <deal.II/matrix_free/operators.h> 
#include <deal.II/matrix_free/fe_evaluation.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_constrained_dofs.h> 
#include <deal.II/multigrid/mg_matrix.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_transfer.h> 
#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer_matrix_free.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 


#include <deal.II/fe/fe_interface_values.h> 
#include <deal.II/meshworker/mesh_loop.h> 

using namespace dealii; 


namespace ChangeVectorTypes 
{ 
  template <typename number> 
  void copy(LA::MPI::Vector &                                         out, 
            const dealii::LinearAlgebra::distributed::Vector<number> &in) 
  { 
    dealii::LinearAlgebra::ReadWriteVector<double> rwv( 
      out.locally_owned_elements()); 
    rwv.import(in, VectorOperation::insert); 
#ifdef USE_PETSC_LA 
    AssertThrow(false, 
                ExcMessage("CopyVectorTypes::copy() not implemented for " 
                           "PETSc vector types.")); 
#else 
    out.import(rwv, VectorOperation::insert); 
#endif 
  } 

  template <typename number> 
  void copy(dealii::LinearAlgebra::distributed::Vector<number> &out, 
            const LA::MPI::Vector &                             in) 
  { 
    dealii::LinearAlgebra::ReadWriteVector<double> rwv; 
#ifdef USE_PETSC_LA 
    (void)in; 
    AssertThrow(false, 
                ExcMessage("CopyVectorTypes::copy() not implemented for " 
                           "PETSc vector types.")); 
#else 
    rwv.reinit(in); 
#endif 
    out.import(rwv, VectorOperation::insert); 
  } 
} // namespace ChangeVectorTypes 


template <int dim> 
class RightHandSide : public Function<dim> 
{ 
public: 
  virtual double value(const Point<dim> & /*p*/, 
                       const unsigned int /*component*/ = 0) const override 
  { 
    return 1.0; 
  } 

  template <typename number> 
  VectorizedArray<number> 
  value(const Point<dim, VectorizedArray<number>> & /*p*/, 
        const unsigned int /*component*/ = 0) const 
  { 
    return VectorizedArray<number>(1.0); 
  } 
}; 


template <int dim> 
class Coefficient : public Function<dim> 
{ 
public: 
  virtual double value(const Point<dim> &p, 
                       const unsigned int /*component*/ = 0) const override; 

  template <typename number> 
  VectorizedArray<number> value(const Point<dim, VectorizedArray<number>> &p, 
                                const unsigned int /*component*/ = 0) const; 

  template <typename number> 
  number average_value(const std::vector<Point<dim, number>> &points) const; 


  template <typename number> 
  std::shared_ptr<Table<2, VectorizedArray<number>>> make_coefficient_table( 
    const MatrixFree<dim, number, VectorizedArray<number>> &mf_storage) const; 
}; 

template <int dim> 
double Coefficient<dim>::value(const Point<dim> &p, const unsigned int) const 
{ 
  for (int d = 0; d < dim; ++d) 
    { 
      if (p[d] < -0.5) 
        return 100.0; 
    } 
  return 1.0; 
} 

template <int dim> 
template <typename number> 
VectorizedArray<number> 
Coefficient<dim>::value(const Point<dim, VectorizedArray<number>> &p, 
                        const unsigned int) const 
{ 
  VectorizedArray<number> return_value = VectorizedArray<number>(1.0); 
  for (unsigned int i = 0; i < VectorizedArray<number>::size(); ++i) 
    { 
      for (int d = 0; d < dim; ++d) 
        if (p[d][i] < -0.5) 
          { 
            return_value[i] = 100.0; 
            break; 
          } 
    } 

  return return_value; 
} 

template <int dim> 
template <typename number> 
number Coefficient<dim>::average_value( 
  const std::vector<Point<dim, number>> &points) const 
{ 
  number average(0); 
  for (unsigned int i = 0; i < points.size(); ++i) 
    average += value(points[i]); 
  average /= points.size(); 

  return average; 
} 

template <int dim> 
template <typename number> 
std::shared_ptr<Table<2, VectorizedArray<number>>> 
Coefficient<dim>::make_coefficient_table( 
  const MatrixFree<dim, number, VectorizedArray<number>> &mf_storage) const 
{ 
  auto coefficient_table = 
    std::make_shared<Table<2, VectorizedArray<number>>>(); 

  FEEvaluation<dim, -1, 0, 1, number> fe_eval(mf_storage); 

  const unsigned int n_cells    = mf_storage.n_cell_batches(); 
  const unsigned int n_q_points = fe_eval.n_q_points; 

  coefficient_table->reinit(n_cells, 1); 

  for (unsigned int cell = 0; cell < n_cells; ++cell) 
    { 
      fe_eval.reinit(cell); 

      VectorizedArray<number> average_value = 0.; 
      for (unsigned int q = 0; q < n_q_points; ++q) 
        average_value += value(fe_eval.quadrature_point(q)); 
      average_value /= n_q_points; 

      (*coefficient_table)(cell, 0) = average_value; 
    } 

  return coefficient_table; 
} 



struct Settings 
{ 
  bool try_parse(const std::string &prm_filename); 

  enum SolverType 
  { 
    gmg_mb, 
    gmg_mf, 
    amg 
  }; 

  SolverType solver; 

  int          dimension; 
  double       smoother_dampen; 
  unsigned int smoother_steps; 
  unsigned int n_steps; 
  bool         output; 
}; 

bool Settings::try_parse(const std::string &prm_filename) 
{ 
  ParameterHandler prm; 
  prm.declare_entry("dim", "2", Patterns::Integer(), "The problem dimension."); 
  prm.declare_entry("n_steps", 
                    "10", 
                    Patterns::Integer(0), 
                    "Number of adaptive refinement steps."); 
  prm.declare_entry("smoother dampen", 
                    "1.0", 
                    Patterns::Double(0.0), 
                    "Dampen factor for the smoother."); 
  prm.declare_entry("smoother steps", 
                    "1", 
                    Patterns::Integer(1), 
                    "Number of smoother steps."); 
  prm.declare_entry("solver", 
                    "MF", 
                    Patterns::Selection("MF|MB|AMG"), 
                    "Switch between matrix-free GMG, " 
                    "matrix-based GMG, and AMG."); 
  prm.declare_entry("output", 
                    "false", 
                    Patterns::Bool(), 
                    "Output graphical results."); 

  if (prm_filename.size() == 0) 
    { 
      std::cout << "****  Error: No input file provided!\n" 
                << "****  Error: Call this program as './step-50 input.prm\n" 
                << "\n" 
                << "****  You may want to use one of the input files in this\n" 
                << "****  directory, or use the following default values\n" 
                << "****  to create an input file:\n"; 
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
        prm.print_parameters(std::cout, ParameterHandler::Text); 
      return false; 
    } 

  try 
    { 
      prm.parse_input(prm_filename); 
    } 
  catch (std::exception &e) 
    { 
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
        std::cerr << e.what() << std::endl; 
      return false; 
    } 

  if (prm.get("solver") == "MF") 
    this->solver = gmg_mf; 
  else if (prm.get("solver") == "MB") 
    this->solver = gmg_mb; 
  else if (prm.get("solver") == "AMG") 
    this->solver = amg; 
  else 
    AssertThrow(false, ExcNotImplemented()); 

  this->dimension       = prm.get_integer("dim"); 
  this->n_steps         = prm.get_integer("n_steps"); 
  this->smoother_dampen = prm.get_double("smoother dampen"); 
  this->smoother_steps  = prm.get_integer("smoother steps"); 
  this->output          = prm.get_bool("output"); 

  return true; 
} 



template <int dim, int degree> 
class LaplaceProblem 
{ 
public: 
  LaplaceProblem(const Settings &settings); 
  void run(); 

private: 


  using MatrixType         = LA::MPI::SparseMatrix; 
  using VectorType         = LA::MPI::Vector; 
  using PreconditionAMG    = LA::MPI::PreconditionAMG; 
  using PreconditionJacobi = LA::MPI::PreconditionJacobi; 

  using MatrixFreeLevelMatrix = MatrixFreeOperators::LaplaceOperator< 
    dim, 
    degree, 
    degree + 1, 
    1, 
    LinearAlgebra::distributed::Vector<float>>; 
  using MatrixFreeActiveMatrix = MatrixFreeOperators::LaplaceOperator< 
    dim, 
    degree, 
    degree + 1, 
    1, 
    LinearAlgebra::distributed::Vector<double>>; 

  using MatrixFreeLevelVector  = LinearAlgebra::distributed::Vector<float>; 
  using MatrixFreeActiveVector = LinearAlgebra::distributed::Vector<double>; 

  void setup_system(); 
  void setup_multigrid(); 
  void assemble_system(); 
  void assemble_multigrid(); 
  void assemble_rhs(); 
  void solve(); 
  void estimate(); 
  void refine_grid(); 
  void output_results(const unsigned int cycle); 

  Settings settings; 

  MPI_Comm           mpi_communicator; 
  ConditionalOStream pcout; 

  parallel::distributed::Triangulation<dim> triangulation; 
  const MappingQ1<dim>                      mapping; 
  FE_Q<dim>                                 fe; 

  DoFHandler<dim> dof_handler; 

 
  IndexSet                  locally_relevant_dofs; 
  AffineConstraints<double> constraints; 

  MatrixType             system_matrix; 
  MatrixFreeActiveMatrix mf_system_matrix; 
  VectorType             solution; 
  VectorType             right_hand_side; 
  Vector<double>         estimated_error_square_per_cell; 

  MGLevelObject<MatrixType> mg_matrix; 
  MGLevelObject<MatrixType> mg_interface_in; 
  MGConstrainedDoFs         mg_constrained_dofs; 

  MGLevelObject<MatrixFreeLevelMatrix> mf_mg_matrix; 

  TimerOutput computing_timer; 
}; 


template <int dim, int degree> 
LaplaceProblem<dim, degree>::LaplaceProblem(const Settings &settings) 
  : settings(settings) 
  , mpi_communicator(MPI_COMM_WORLD) 
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
  , triangulation(mpi_communicator, 
                  Triangulation<dim>::limit_level_difference_at_vertices, 
                  (settings.solver == Settings::amg) ? 
                    parallel::distributed::Triangulation<dim>::default_setting : 
                    parallel::distributed::Triangulation< 
                      dim>::construct_multigrid_hierarchy) 
  , mapping() 
  , fe(degree) 
  , dof_handler(triangulation) 
  , computing_timer(pcout, TimerOutput::never, TimerOutput::wall_times) 
{ 
  GridGenerator::hyper_L(triangulation, -1., 1., /*colorize*/ false); 
  triangulation.refine_global(1); 
} 



template <int dim, int degree> 
void LaplaceProblem<dim, degree>::setup_system() 
{ 
  TimerOutput::Scope timing(computing_timer, "Setup"); 

  dof_handler.distribute_dofs(fe); 

  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
  locally_owned_dofs = dof_handler.locally_owned_dofs(); 

  solution.reinit(locally_owned_dofs, mpi_communicator); 
  right_hand_side.reinit(locally_owned_dofs, mpi_communicator); 
  constraints.reinit(locally_relevant_dofs); 
  DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

  VectorTools::interpolate_boundary_values( 
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints); 
  constraints.close(); 

  switch (settings.solver) 
    { 
      case Settings::gmg_mf: 
        { 
          typename MatrixFree<dim, double>::AdditionalData additional_data; 
          additional_data.tasks_parallel_scheme = 
            MatrixFree<dim, double>::AdditionalData::none; 
          additional_data.mapping_update_flags = 
            (update_gradients | update_JxW_values | update_quadrature_points); 
          std::shared_ptr<MatrixFree<dim, double>> mf_storage = 
            std::make_shared<MatrixFree<dim, double>>(); 
          mf_storage->reinit(mapping, 
                             dof_handler, 
                             constraints, 
                             QGauss<1>(degree + 1), 
                             additional_data); 

          mf_system_matrix.initialize(mf_storage); 

          const Coefficient<dim> coefficient; 
          mf_system_matrix.set_coefficient( 
            coefficient.make_coefficient_table(*mf_storage)); 

          break; 
        } 

      case Settings::gmg_mb: 
      case Settings::amg: 
        { 
#ifdef USE_PETSC_LA 
          DynamicSparsityPattern dsp(locally_relevant_dofs); 
          DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 

          SparsityTools::distribute_sparsity_pattern(dsp, 
                                                     locally_owned_dofs, 
                                                     mpi_communicator, 
                                                     locally_relevant_dofs); 

          system_matrix.reinit(locally_owned_dofs, 
                               locally_owned_dofs, 
                               dsp, 
                               mpi_communicator); 
#else 
          TrilinosWrappers::SparsityPattern dsp(locally_owned_dofs, 
                                                locally_owned_dofs, 
                                                locally_relevant_dofs, 
                                                mpi_communicator); 
          DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
          dsp.compress(); 
          system_matrix.reinit(dsp); 
#endif 

          break; 
        } 

      default: 
        Assert(false, ExcNotImplemented()); 
    } 
} 



template <int dim, int degree> 
void LaplaceProblem<dim, degree>::setup_multigrid() 
{ 
  TimerOutput::Scope timing(computing_timer, "Setup multigrid"); 

  dof_handler.distribute_mg_dofs(); 

  mg_constrained_dofs.clear(); 
  mg_constrained_dofs.initialize(dof_handler); 

  const std::set<types::boundary_id> boundary_ids = {types::boundary_id(0)}; 
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, boundary_ids); 

  const unsigned int n_levels = triangulation.n_global_levels(); 

  switch (settings.solver) 
    { 
      case Settings::gmg_mf: 
        { 
          mf_mg_matrix.resize(0, n_levels - 1); 

          for (unsigned int level = 0; level < n_levels; ++level) 
            { 
              IndexSet relevant_dofs; 
              DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
                                                            level, 
                                                            relevant_dofs); 
              AffineConstraints<double> level_constraints; 
              level_constraints.reinit(relevant_dofs); 
              level_constraints.add_lines( 
                mg_constrained_dofs.get_boundary_indices(level)); 
              level_constraints.close(); 

              typename MatrixFree<dim, float>::AdditionalData additional_data; 
              additional_data.tasks_parallel_scheme = 
                MatrixFree<dim, float>::AdditionalData::none; 
              additional_data.mapping_update_flags = 
                (update_gradients | update_JxW_values | 
                 update_quadrature_points); 
              additional_data.mg_level = level; 
              std::shared_ptr<MatrixFree<dim, float>> mf_storage_level( 
                new MatrixFree<dim, float>()); 
              mf_storage_level->reinit(mapping, 
                                       dof_handler, 
                                       level_constraints, 
                                       QGauss<1>(degree + 1), 
                                       additional_data); 

              mf_mg_matrix[level].initialize(mf_storage_level, 
                                             mg_constrained_dofs, 
                                             level); 

              const Coefficient<dim> coefficient; 
              mf_mg_matrix[level].set_coefficient( 
                coefficient.make_coefficient_table(*mf_storage_level)); 

              mf_mg_matrix[level].compute_diagonal(); 
            } 

          break; 
        } 

      case Settings::gmg_mb: 
        { 
          mg_matrix.resize(0, n_levels - 1); 
          mg_matrix.clear_elements(); 
          mg_interface_in.resize(0, n_levels - 1); 
          mg_interface_in.clear_elements(); 

          for (unsigned int level = 0; level < n_levels; ++level) 
            { 
              IndexSet dof_set; 
              DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
                                                            level, 
                                                            dof_set); 

              { 
#ifdef USE_PETSC_LA 
                DynamicSparsityPattern dsp(dof_set); 
                MGTools::make_sparsity_pattern(dof_handler, dsp, level); 
                dsp.compress(); 
                SparsityTools::distribute_sparsity_pattern( 
                  dsp, 
                  dof_handler.locally_owned_mg_dofs(level), 
                  mpi_communicator, 
                  dof_set); 

                mg_matrix[level].reinit( 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dsp, 
                  mpi_communicator); 
#else 
                TrilinosWrappers::SparsityPattern dsp( 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dof_set, 
                  mpi_communicator); 
                MGTools::make_sparsity_pattern(dof_handler, dsp, level); 

                dsp.compress(); 
                mg_matrix[level].reinit(dsp); 
#endif 
              } 

              { 
#ifdef USE_PETSC_LA 
                DynamicSparsityPattern dsp(dof_set); 
                MGTools::make_interface_sparsity_pattern(dof_handler, 
                                                         mg_constrained_dofs, 
                                                         dsp, 
                                                         level); 
                dsp.compress(); 
                SparsityTools::distribute_sparsity_pattern( 
                  dsp, 
                  dof_handler.locally_owned_mg_dofs(level), 
                  mpi_communicator, 
                  dof_set); 

                mg_interface_in[level].reinit( 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dsp, 
                  mpi_communicator); 
#else 
                TrilinosWrappers::SparsityPattern dsp( 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dof_handler.locally_owned_mg_dofs(level), 
                  dof_set, 
                  mpi_communicator); 

                MGTools::make_interface_sparsity_pattern(dof_handler, 
                                                         mg_constrained_dofs, 
                                                         dsp, 
                                                         level); 
                dsp.compress(); 
                mg_interface_in[level].reinit(dsp); 
#endif 
              } 
            } 
          break; 
        } 

      default: 
        Assert(false, ExcNotImplemented()); 
    } 
} 



template <int dim, int degree> 
void LaplaceProblem<dim, degree>::assemble_system() 
{ 
  TimerOutput::Scope timing(computing_timer, "Assemble"); 

  const QGauss<dim> quadrature_formula(degree + 1); 

  FEValues<dim> fe_values(fe, 
                          quadrature_formula, 
                          update_values | update_gradients | 
                            update_quadrature_points | update_JxW_values); 

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
  const unsigned int n_q_points    = quadrature_formula.size(); 

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
  Vector<double>     cell_rhs(dofs_per_cell); 

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

  const Coefficient<dim> coefficient; 
  RightHandSide<dim>     rhs; 
  std::vector<double>    rhs_values(n_q_points); 

  for (const auto &cell : dof_handler.active_cell_iterators()) 
    if (cell->is_locally_owned()) 
      { 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        fe_values.reinit(cell); 

        const double coefficient_value = 
          coefficient.average_value(fe_values.get_quadrature_points()); 
        rhs.value_list(fe_values.get_quadrature_points(), rhs_values); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                cell_matrix(i, j) += 
                  coefficient_value *                // epsilon(x) 
                  fe_values.shape_grad(i, q_point) * // * grad phi_i(x) 
                  fe_values.shape_grad(j, q_point) * // * grad phi_j(x) 
                  fe_values.JxW(q_point);            // * dx 

              cell_rhs(i) += 
                fe_values.shape_value(i, q_point) * // grad phi_i(x) 
                rhs_values[q_point] *               // * f(x) 
                fe_values.JxW(q_point);             // * dx 
            } 

        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global(cell_matrix, 
                                               cell_rhs, 
                                               local_dof_indices, 
                                               system_matrix, 
                                               right_hand_side); 
      } 

  system_matrix.compress(VectorOperation::add); 
  right_hand_side.compress(VectorOperation::add); 
} 


template <int dim, int degree> 
void LaplaceProblem<dim, degree>::assemble_multigrid() 
{ 
  TimerOutput::Scope timing(computing_timer, "Assemble multigrid"); 

  QGauss<dim> quadrature_formula(degree + 1); 

  FEValues<dim> fe_values(fe, 
                          quadrature_formula, 
                          update_values | update_gradients | 
                            update_quadrature_points | update_JxW_values); 

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
  const unsigned int n_q_points    = quadrature_formula.size(); 

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

  const Coefficient<dim> coefficient; 

  std::vector<AffineConstraints<double>> boundary_constraints( 
    triangulation.n_global_levels()); 
  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level) 
    { 
      IndexSet dof_set; 
      DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
                                                    level, 
                                                    dof_set); 
      boundary_constraints[level].reinit(dof_set); 
      boundary_constraints[level].add_lines( 
        mg_constrained_dofs.get_refinement_edge_indices(level)); 
      boundary_constraints[level].add_lines( 
        mg_constrained_dofs.get_boundary_indices(level)); 

      boundary_constraints[level].close(); 
    } 

  for (const auto &cell : dof_handler.cell_iterators()) 
    if (cell->level_subdomain_id() == triangulation.locally_owned_subdomain()) 
      { 
        cell_matrix = 0; 
        fe_values.reinit(cell); 

        const double coefficient_value = 
          coefficient.average_value(fe_values.get_quadrature_points()); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              cell_matrix(i, j) += 
                coefficient_value * fe_values.shape_grad(i, q_point) * 
                fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point); 

        cell->get_mg_dof_indices(local_dof_indices); 

        boundary_constraints[cell->level()].distribute_local_to_global( 
          cell_matrix, local_dof_indices, mg_matrix[cell->level()]); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            if (mg_constrained_dofs.is_interface_matrix_entry( 
                  cell->level(), local_dof_indices[i], local_dof_indices[j])) 
              mg_interface_in[cell->level()].add(local_dof_indices[i], 
                                                 local_dof_indices[j], 
                                                 cell_matrix(i, j)); 
      } 

  for (unsigned int i = 0; i < triangulation.n_global_levels(); ++i) 
    { 
      mg_matrix[i].compress(VectorOperation::add); 
      mg_interface_in[i].compress(VectorOperation::add); 
    } 
} 








template <int dim, int degree> 
void LaplaceProblem<dim, degree>::assemble_rhs() 
{ 
  TimerOutput::Scope timing(computing_timer, "Assemble right-hand side"); 

  MatrixFreeActiveVector solution_copy; 
  MatrixFreeActiveVector right_hand_side_copy; 
  mf_system_matrix.initialize_dof_vector(solution_copy); 
  mf_system_matrix.initialize_dof_vector(right_hand_side_copy); 

  solution_copy = 0.; 
  constraints.distribute(solution_copy); 
  solution_copy.update_ghost_values(); 
  right_hand_side_copy = 0; 
  const Table<2, VectorizedArray<double>> &coefficient = 
    *(mf_system_matrix.get_coefficient()); 

  RightHandSide<dim> right_hand_side_function; 

  FEEvaluation<dim, degree, degree + 1, 1, double> phi( 
    *mf_system_matrix.get_matrix_free()); 

  for (unsigned int cell = 0; 
       cell < mf_system_matrix.get_matrix_free()->n_cell_batches(); 
       ++cell) 
    { 
      phi.reinit(cell); 
      phi.read_dof_values_plain(solution_copy); 
      phi.evaluate(EvaluationFlags::gradients); 

      for (unsigned int q = 0; q < phi.n_q_points; ++q) 
        { 
          phi.submit_gradient(-1.0 * 
                                (coefficient(cell, 0) * phi.get_gradient(q)), 
                              q); 
          phi.submit_value( 
            right_hand_side_function.value(phi.quadrature_point(q)), q); 
        } 

      phi.integrate_scatter(EvaluationFlags::values | 
                              EvaluationFlags::gradients, 
                            right_hand_side_copy); 
    } 

  right_hand_side_copy.compress(VectorOperation::add); 

  ChangeVectorTypes::copy(right_hand_side, right_hand_side_copy); 
} 



template <int dim, int degree> 
void LaplaceProblem<dim, degree>::solve() 
{ 
  TimerOutput::Scope timing(computing_timer, "Solve"); 

  SolverControl solver_control(1000, 1.e-10 * right_hand_side.l2_norm()); 
  solver_control.enable_history_data(); 

  solution = 0.; 


  switch (settings.solver) 
    { 
      case Settings::gmg_mf: 
        { 
          computing_timer.enter_subsection("Solve: Preconditioner setup"); 

          MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs); 
          mg_transfer.build(dof_handler); 

          SolverControl coarse_solver_control(1000, 1e-12, false, false); 
          SolverCG<MatrixFreeLevelVector> coarse_solver(coarse_solver_control); 
          PreconditionIdentity            identity; 
          MGCoarseGridIterativeSolver<MatrixFreeLevelVector, 
                                      SolverCG<MatrixFreeLevelVector>, 
                                      MatrixFreeLevelMatrix, 
                                      PreconditionIdentity> 
            coarse_grid_solver(coarse_solver, mf_mg_matrix[0], identity); 

          using Smoother = dealii::PreconditionJacobi<MatrixFreeLevelMatrix>; 
          MGSmootherPrecondition<MatrixFreeLevelMatrix, 
                                 Smoother, 
                                 MatrixFreeLevelVector> 
            smoother; 
          smoother.initialize(mf_mg_matrix, 
                              typename Smoother::AdditionalData( 
                                settings.smoother_dampen)); 
          smoother.set_steps(settings.smoother_steps); 

          mg::Matrix<MatrixFreeLevelVector> mg_m(mf_mg_matrix); 

          MGLevelObject< 
            MatrixFreeOperators::MGInterfaceOperator<MatrixFreeLevelMatrix>> 
            mg_interface_matrices; 
          mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1); 
          for (unsigned int level = 0; level < triangulation.n_global_levels(); 
               ++level) 
            mg_interface_matrices[level].initialize(mf_mg_matrix[level]); 
          mg::Matrix<MatrixFreeLevelVector> mg_interface(mg_interface_matrices); 

          Multigrid<MatrixFreeLevelVector> mg( 
            mg_m, coarse_grid_solver, mg_transfer, smoother, smoother); 
          mg.set_edge_matrices(mg_interface, mg_interface); 

          PreconditionMG<dim, 
                         MatrixFreeLevelVector, 
                         MGTransferMatrixFree<dim, float>> 
            preconditioner(dof_handler, mg, mg_transfer); 


          MatrixFreeActiveVector solution_copy; 
          MatrixFreeActiveVector right_hand_side_copy; 
          mf_system_matrix.initialize_dof_vector(solution_copy); 
          mf_system_matrix.initialize_dof_vector(right_hand_side_copy); 

          ChangeVectorTypes::copy(solution_copy, solution); 
          ChangeVectorTypes::copy(right_hand_side_copy, right_hand_side); 
          computing_timer.leave_subsection("Solve: Preconditioner setup"); 


          { 
            TimerOutput::Scope timing(computing_timer, 
                                      "Solve: 1 multigrid V-cycle"); 
            preconditioner.vmult(solution_copy, right_hand_side_copy); 
          } 
          solution_copy = 0.; 


          { 
            SolverCG<MatrixFreeActiveVector> solver(solver_control); 

            TimerOutput::Scope timing(computing_timer, "Solve: CG"); 
            solver.solve(mf_system_matrix, 
                         solution_copy, 
                         right_hand_side_copy, 
                         preconditioner); 
          } 

          solution_copy.update_ghost_values(); 
          ChangeVectorTypes::copy(solution, solution_copy); 
          constraints.distribute(solution); 

          break; 
        } 


      case Settings::gmg_mb: 
        { 
          computing_timer.enter_subsection("Solve: Preconditioner setup"); 

          MGTransferPrebuilt<VectorType> mg_transfer(mg_constrained_dofs); 
          mg_transfer.build(dof_handler); 

          SolverControl        coarse_solver_control(1000, 1e-12, false, false); 
          SolverCG<VectorType> coarse_solver(coarse_solver_control); 
          PreconditionIdentity identity; 
          MGCoarseGridIterativeSolver<VectorType, 
                                      SolverCG<VectorType>, 
                                      MatrixType, 
                                      PreconditionIdentity> 
            coarse_grid_solver(coarse_solver, mg_matrix[0], identity); 

          using Smoother = LA::MPI::PreconditionJacobi; 
          MGSmootherPrecondition<MatrixType, Smoother, VectorType> smoother; 

#ifdef USE_PETSC_LA 
          smoother.initialize(mg_matrix); 
          Assert( 
            settings.smoother_dampen == 1.0, 
            ExcNotImplemented( 
              "PETSc's PreconditionJacobi has no support for a damping parameter.")); 
#else 
          smoother.initialize(mg_matrix, settings.smoother_dampen); 
#endif 

          smoother.set_steps(settings.smoother_steps); 

          mg::Matrix<VectorType> mg_m(mg_matrix); 
          mg::Matrix<VectorType> mg_in(mg_interface_in); 
          mg::Matrix<VectorType> mg_out(mg_interface_in); 

          Multigrid<VectorType> mg( 
            mg_m, coarse_grid_solver, mg_transfer, smoother, smoother); 
          mg.set_edge_matrices(mg_out, mg_in); 

          PreconditionMG<dim, VectorType, MGTransferPrebuilt<VectorType>> 
            preconditioner(dof_handler, mg, mg_transfer); 

          computing_timer.leave_subsection("Solve: Preconditioner setup"); 


          { 
            TimerOutput::Scope timing(computing_timer, 
                                      "Solve: 1 multigrid V-cycle"); 
            preconditioner.vmult(solution, right_hand_side); 
          } 
          solution = 0.; 


          { 
            SolverCG<VectorType> solver(solver_control); 

            TimerOutput::Scope timing(computing_timer, "Solve: CG"); 
            solver.solve(system_matrix, 
                         solution, 
                         right_hand_side, 
                         preconditioner); 
          } 

          constraints.distribute(solution); 

 
        } 


      case Settings::amg: 
        { 
          computing_timer.enter_subsection("Solve: Preconditioner setup"); 

          PreconditionAMG                 preconditioner; 
          PreconditionAMG::AdditionalData Amg_data; 

#ifdef USE_PETSC_LA 
          Amg_data.symmetric_operator = true; 
#else 
          Amg_data.elliptic              = true; 
          Amg_data.smoother_type         = "Jacobi"; 
          Amg_data.higher_order_elements = true; 
          Amg_data.smoother_sweeps       = settings.smoother_steps; 
          Amg_data.aggregation_threshold = 0.02; 
#endif 

          Amg_data.output_details = false; 

          preconditioner.initialize(system_matrix, Amg_data); 
          computing_timer.leave_subsection("Solve: Preconditioner setup"); 


          { 
            TimerOutput::Scope timing(computing_timer, 
                                      "Solve: 1 multigrid V-cycle"); 
            preconditioner.vmult(solution, right_hand_side); 
          } 
          solution = 0.; 


          { 
            SolverCG<VectorType> solver(solver_control); 

            TimerOutput::Scope timing(computing_timer, "Solve: CG"); 
            solver.solve(system_matrix, 
                         solution, 
                         right_hand_side, 
                         preconditioner); 
          } 
          constraints.distribute(solution); 

          break; 
        } 

      default: 
        Assert(false, ExcInternalError()); 
    } 

  pcout << "   Number of CG iterations:      " << solver_control.last_step() 
        << std::endl; 
} 


template <int dim> 
struct ScratchData 
{ 
  ScratchData(const Mapping<dim> &      mapping, 
              const FiniteElement<dim> &fe, 
              const unsigned int        quadrature_degree, 
              const UpdateFlags         update_flags, 
              const UpdateFlags         interface_update_flags) 
    : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags) 
    , fe_interface_values(mapping, 
                          fe, 
                          QGauss<dim - 1>(quadrature_degree), 
                          interface_update_flags) 
  {} 

  ScratchData(const ScratchData<dim> &scratch_data) 
    : fe_values(scratch_data.fe_values.get_mapping(), 
                scratch_data.fe_values.get_fe(), 
                scratch_data.fe_values.get_quadrature(), 
                scratch_data.fe_values.get_update_flags()) 
    , fe_interface_values(scratch_data.fe_values.get_mapping(), 
                          scratch_data.fe_values.get_fe(), 
                          scratch_data.fe_interface_values.get_quadrature(), 
                          scratch_data.fe_interface_values.get_update_flags()) 
  {} 

  FEValues<dim>          fe_values; 
  FEInterfaceValues<dim> fe_interface_values; 
}; 

struct CopyData 
{ 
  CopyData() 
    : cell_index(numbers::invalid_unsigned_int) 
    , value(0.) 
  {} 

  CopyData(const CopyData &) = default; 

  struct FaceData 
  { 
    unsigned int cell_indices[2]; 
    double       values[2]; 
  }; 

  unsigned int          cell_index; 
  double                value; 
  std::vector<FaceData> face_data; 
}; 

template <int dim, int degree> 
void LaplaceProblem<dim, degree>::estimate() 
{ 
  TimerOutput::Scope timing(computing_timer, "Estimate"); 

  VectorType temp_solution; 
  temp_solution.reinit(locally_owned_dofs, 
                       locally_relevant_dofs, 
                       mpi_communicator); 
  temp_solution = solution; 

  const Coefficient<dim> coefficient; 

  estimated_error_square_per_cell.reinit(triangulation.n_active_cells()); 

  using Iterator = typename DoFHandler<dim>::active_cell_iterator; 

  auto cell_worker = [&](const Iterator &  cell, 
                         ScratchData<dim> &scratch_data, 
                         CopyData &        copy_data) { 
    FEValues<dim> &fe_values = scratch_data.fe_values; 
    fe_values.reinit(cell); 

    RightHandSide<dim> rhs; 
    const double       rhs_value = rhs.value(cell->center()); 

    const double nu = coefficient.value(cell->center()); 

 
    fe_values.get_function_hessians(temp_solution, hessians); 

    copy_data.cell_index = cell->active_cell_index(); 

    double residual_norm_square = 0.; 
    for (unsigned k = 0; k < fe_values.n_quadrature_points; ++k) 
      { 
        const double residual = (rhs_value + nu * trace(hessians[k])); 
        residual_norm_square += residual * residual * fe_values.JxW(k); 
      } 

    copy_data.value = 
      cell->diameter() * cell->diameter() * residual_norm_square; 
  }; 

  auto face_worker = [&](const Iterator &    cell, 
                         const unsigned int &f, 
                         const unsigned int &sf, 
                         const Iterator &    ncell, 
                         const unsigned int &nf, 
                         const unsigned int &nsf, 
                         ScratchData<dim> &  scratch_data, 
                         CopyData &          copy_data) { 
    FEInterfaceValues<dim> &fe_interface_values = 
      scratch_data.fe_interface_values; 
    fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf); 

    copy_data.face_data.emplace_back(); 
    CopyData::FaceData &copy_data_face = copy_data.face_data.back(); 

    copy_data_face.cell_indices[0] = cell->active_cell_index(); 
    copy_data_face.cell_indices[1] = ncell->active_cell_index(); 

    const double coeff1 = coefficient.value(cell->center()); 
    const double coeff2 = coefficient.value(ncell->center()); 

    std::vector<Tensor<1, dim>> grad_u[2]; 

    for (unsigned int i = 0; i < 2; ++i) 
      { 
        grad_u[i].resize(fe_interface_values.n_quadrature_points); 
        fe_interface_values.get_fe_face_values(i).get_function_gradients( 
          temp_solution, grad_u[i]); 
      } 

    double jump_norm_square = 0.; 

    for (unsigned int qpoint = 0; 
         qpoint < fe_interface_values.n_quadrature_points; 
         ++qpoint) 
      { 
        const double jump = 
          coeff1 * grad_u[0][qpoint] * fe_interface_values.normal(qpoint) - 
          coeff2 * grad_u[1][qpoint] * fe_interface_values.normal(qpoint); 

        jump_norm_square += jump * jump * fe_interface_values.JxW(qpoint); 
      } 

    const double h           = cell->face(f)->measure(); 
    copy_data_face.values[0] = 0.5 * h * jump_norm_square; 
    copy_data_face.values[1] = copy_data_face.values[0]; 
  }; 

  auto copier = [&](const CopyData &copy_data) { 
    if (copy_data.cell_index != numbers::invalid_unsigned_int) 
      estimated_error_square_per_cell[copy_data.cell_index] += copy_data.value; 

    for (auto &cdf : copy_data.face_data) 
      for (unsigned int j = 0; j < 2; ++j) 
        estimated_error_square_per_cell[cdf.cell_indices[j]] += cdf.values[j]; 
  }; 

  const unsigned int n_gauss_points = degree + 1; 
  ScratchData<dim>   scratch_data(mapping, 
                                fe, 
                                n_gauss_points, 
                                update_hessians | update_quadrature_points | 
                                  update_JxW_values, 
                                update_values | update_gradients | 
                                  update_JxW_values | update_normal_vectors); 
  CopyData           copy_data; 


  MeshWorker::mesh_loop(dof_handler.begin_active(), 
                        dof_handler.end(), 
                        cell_worker, 
                        copier, 
                        scratch_data, 
                        copy_data, 
                        MeshWorker::assemble_own_cells | 
                          MeshWorker::assemble_ghost_faces_both | 
                          MeshWorker::assemble_own_interior_faces_once, 
                          /*boundary_worker=*/nullptr, face_worker);


 

  const double global_error_estimate = 
    std::sqrt(Utilities::MPI::sum(estimated_error_square_per_cell.l1_norm(), 
                                  mpi_communicator)); 
  pcout << "   Global error estimate:        " << global_error_estimate 
        << std::endl; 
} 


template <int dim, int degree> 
void LaplaceProblem<dim, degree>::refine_grid() 
{ 
  TimerOutput::Scope timing(computing_timer, "Refine grid"); 

  const double refinement_fraction = 1. / (std::pow(2.0, dim) - 1.); 
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number( 
    triangulation, estimated_error_square_per_cell, refinement_fraction, 0.0); 

  triangulation.execute_coarsening_and_refinement(); 
} 


template <int dim, int degree> 
void LaplaceProblem<dim, degree>::output_results(const unsigned int cycle) 
{ 
  TimerOutput::Scope timing(computing_timer, "Output results"); 

  VectorType temp_solution; 
  temp_solution.reinit(locally_owned_dofs, 
                       locally_relevant_dofs, 
                       mpi_communicator); 
  temp_solution = solution; 

  DataOut<dim> data_out; 
  data_out.attach_dof_handler(dof_handler); 
  data_out.add_data_vector(temp_solution, "solution"); 

  Vector<float> subdomain(triangulation.n_active_cells()); 
  for (unsigned int i = 0; i < subdomain.size(); ++i) 
    subdomain(i) = triangulation.locally_owned_subdomain(); 
  data_out.add_data_vector(subdomain, "subdomain"); 

  Vector<float> level(triangulation.n_active_cells()); 
  for (const auto &cell : triangulation.active_cell_iterators()) 
    level(cell->active_cell_index()) = cell->level(); 
  data_out.add_data_vector(level, "level"); 

  if (estimated_error_square_per_cell.size() > 0) 
    data_out.add_data_vector(estimated_error_square_per_cell, 
                             "estimated_error_square_per_cell"); 

  data_out.build_patches(); 

  const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
    "", "solution", cycle, mpi_communicator, 2 /*n_digits*/, 1 /*n_groups*/); 

  pcout << "   Wrote " << pvtu_filename << std::endl; 
} 


template <int dim, int degree> 
void LaplaceProblem<dim, degree>::run() 
{ 
  for (unsigned int cycle = 0; cycle < settings.n_steps; ++cycle) 
    { 
      pcout << "Cycle " << cycle << ':' << std::endl; 
      if (cycle > 0) 
        refine_grid(); 

      pcout << "   Number of active cells:       " 
            << triangulation.n_global_active_cells(); 


      if (settings.solver == Settings::gmg_mf || 
          settings.solver == Settings::gmg_mb) 
        pcout << " (" << triangulation.n_global_levels() << " global levels)" 
              << std::endl 
              << "   Partition efficiency:         " 
              << 1.0 / MGTools::workload_imbalance(triangulation); 
      pcout << std::endl; 

      setup_system(); 


      if (settings.solver == Settings::gmg_mf || 
          settings.solver == Settings::gmg_mb) 
        setup_multigrid(); 

      pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs(); 
      if (settings.solver == Settings::gmg_mf || 
          settings.solver == Settings::gmg_mb) 
        { 
          pcout << " (by level: "; 
          for (unsigned int level = 0; level < triangulation.n_global_levels(); 
               ++level) 
            pcout << dof_handler.n_dofs(level) 
                  << (level == triangulation.n_global_levels() - 1 ? ")" : 
                                                                     ", "); 
        } 
      pcout << std::endl; 


      if (settings.solver == Settings::gmg_mf) 
        assemble_rhs(); 
      else /*gmg_mb or amg*/ 
        { 
          assemble_system(); 
          if (settings.solver == Settings::gmg_mb) 
            assemble_multigrid(); 
        } 

      solve(); 
      estimate(); 

      if (settings.output) 
        output_results(cycle); 

      computing_timer.print_summary(); 
      computing_timer.reset(); 
    } 
} 


int main(int argc, char *argv[]) 
{ 
  using namespace dealii; 
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

  Settings settings; 
  if (!settings.try_parse((argc > 1) ? (argv[1]) : "")) 
    return 0; 

  try 
    { 
      constexpr unsigned int fe_degree = 2; 

      switch (settings.dimension) 
        { 
          case 2: 
            { 
              LaplaceProblem<2, fe_degree> test(settings); 
              test.run(); 

              break; 
            } 

          case 3: 
            { 
              LaplaceProblem<3, fe_degree> test(settings); 
              test.run(); 

              break; 
            } 

          default: 
            Assert(false, ExcMessage("This program only works in 2d and 3d.")); 
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
      MPI_Abort(MPI_COMM_WORLD, 1); 
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
      MPI_Abort(MPI_COMM_WORLD, 2); 
      return 1; 
    } 

  return 0; 
} 

