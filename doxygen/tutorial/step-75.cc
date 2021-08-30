

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
 * Author: Marc Fehling, Colorado State University, 2021 
 *         Peter Munch, Technical University of Munich and Helmholtz-Zentrum 
 *                      hereon, 2021 
 *         Wolfgang Bangerth, Colorado State University, 2021 
 */ 




#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/index_set.h> 
#include <deal.II/base/mpi.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/timer.h> 

#include <deal.II/distributed/grid_refinement.h> 
#include <deal.II/distributed/tria.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/grid/grid_generator.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_series.h> 

#include <deal.II/hp/fe_collection.h> 
#include <deal.II/hp/refinement.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/trilinos_precondition.h> 
#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/smoothness_estimator.h> 
#include <deal.II/numerics/vector_tools.h> 

#include <algorithm> 
#include <fstream> 
#include <iostream> 

#include <deal.II/distributed/cell_weights.h> 


#include <deal.II/base/function.h> 
#include <deal.II/base/geometric_utilities.h> 


#include <deal.II/matrix_free/matrix_free.h> 
#include <deal.II/matrix_free/fe_evaluation.h> 
#include <deal.II/matrix_free/tools.h> 


#include <deal.II/lac/la_parallel_vector.h> 


#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_constrained_dofs.h> 
#include <deal.II/multigrid/mg_matrix.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_transfer_global_coarsening.h> 
#include <deal.II/multigrid/multigrid.h> 

namespace Step75 
{ 
  using namespace dealii; 


  template <int dim> 
  class Solution : public Function<dim> 
  { 
  public: 
    Solution() 
      : Function<dim>() 
    {} 

    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/) const override 
    { 
      const std::array<double, dim> p_sphere = 
        GeometricUtilities::Coordinates::to_spherical(p); 

      constexpr const double alpha = 2. / 3.; 
      return std::pow(p_sphere[0], alpha) * std::sin(alpha * p_sphere[1]); 
    } 
  }; 




  struct MultigridParameters 
  { 
    struct 
    { 
      std::string  type            = "cg_with_amg"; 
      unsigned int maxiter         = 10000; 
      double       abstol          = 1e-20; 
      double       reltol          = 1e-4; 
      unsigned int smoother_sweeps = 1; 
      unsigned int n_cycles        = 1; 
      std::string  smoother_type   = "ILU"; 
    } coarse_solver; 

    struct 
    { 
      std::string  type                = "chebyshev"; 
      double       smoothing_range     = 20; 
      unsigned int degree              = 5; 
      unsigned int eig_cg_n_iterations = 20; 
    } smoother; 

    struct 
    { 
      MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType 
        p_sequence = MGTransferGlobalCoarseningTools:: 
          PolynomialCoarseningSequenceType::decrease_by_one; 
      bool perform_h_transfer = true; 
    } transfer; 
  }; 


  struct Parameters 
  { 
    unsigned int n_cycles         = 8; 
    double       tolerance_factor = 1e-12; 

    MultigridParameters mg_data; 

    unsigned int min_h_level            = 5; 
    unsigned int max_h_level            = 12; 
    unsigned int min_p_degree           = 2; 
    unsigned int max_p_degree           = 6; 
    unsigned int max_p_level_difference = 1; 

    double refine_fraction    = 0.3; 
    double coarsen_fraction   = 0.03; 
    double p_refine_fraction  = 0.9; 
    double p_coarsen_fraction = 0.9; 

    double weighting_factor   = 1e6; 
    double weighting_exponent = 1.; 
  }; 




  template <int dim, typename number> 
  class LaplaceOperator : public Subscriptor 
  { 
  public: 
    using VectorType = LinearAlgebra::distributed::Vector<number>; 

    using FECellIntegrator = FEEvaluation<dim, -1, 0, 1, number>; 

    LaplaceOperator() = default; 

    LaplaceOperator(const hp::MappingCollection<dim> &mapping, 
                    const DoFHandler<dim> &           dof_handler, 
                    const hp::QCollection<dim> &      quad, 
                    const AffineConstraints<number> & constraints, 
                    VectorType &                      system_rhs); 

    void reinit(const hp::MappingCollection<dim> &mapping, 
                const DoFHandler<dim> &           dof_handler, 
                const hp::QCollection<dim> &      quad, 
                const AffineConstraints<number> & constraints, 
                VectorType &                      system_rhs); 

    types::global_dof_index m() const; 

    number el(unsigned int, unsigned int) const; 

    void initialize_dof_vector(VectorType &vec) const; 

    void vmult(VectorType &dst, const VectorType &src) const; 

    void Tvmult(VectorType &dst, const VectorType &src) const; 

    const TrilinosWrappers::SparseMatrix &get_system_matrix() const; 

    void compute_inverse_diagonal(VectorType &diagonal) const; 

  private: 
    void do_cell_integral_local(FECellIntegrator &integrator) const; 

    void do_cell_integral_global(FECellIntegrator &integrator, 
                                 VectorType &      dst, 
                                 const VectorType &src) const; 

    void do_cell_integral_range( 
      const MatrixFree<dim, number> &              matrix_free, 
      VectorType &                                 dst, 
      const VectorType &                           src, 
      const std::pair<unsigned int, unsigned int> &range) const; 

    MatrixFree<dim, number> matrix_free; 


    AffineConstraints<number>              constraints; 
    mutable TrilinosWrappers::SparseMatrix system_matrix; 
  }; 


  template <int dim, typename number> 
  LaplaceOperator<dim, number>::LaplaceOperator( 
    const hp::MappingCollection<dim> &mapping, 
    const DoFHandler<dim> &           dof_handler, 
    const hp::QCollection<dim> &      quad, 
    const AffineConstraints<number> & constraints, 
    VectorType &                      system_rhs) 
  { 
    this->reinit(mapping, dof_handler, quad, constraints, system_rhs); 
  } 

  template <int dim, typename number> 
  void LaplaceOperator<dim, number>::reinit( 
    const hp::MappingCollection<dim> &mapping, 
    const DoFHandler<dim> &           dof_handler, 
    const hp::QCollection<dim> &      quad, 
    const AffineConstraints<number> & constraints, 
    VectorType &                      system_rhs) 
  { 


    this->system_matrix.clear(); 


    this->constraints.copy_from(constraints); 


    typename MatrixFree<dim, number>::AdditionalData data; 
    data.mapping_update_flags = update_gradients; 

    matrix_free.reinit(mapping, dof_handler, constraints, quad, data); 


    { 
      AffineConstraints<number> constraints_without_dbc; 

      IndexSet locally_relevant_dofs; 
      DoFTools::extract_locally_relevant_dofs(dof_handler, 
                                              locally_relevant_dofs); 
      constraints_without_dbc.reinit(locally_relevant_dofs); 

      DoFTools::make_hanging_node_constraints(dof_handler, 
                                              constraints_without_dbc); 
      constraints_without_dbc.close(); 

      VectorType b, x; 

      this->initialize_dof_vector(system_rhs); 

      MatrixFree<dim, number> matrix_free; 
      matrix_free.reinit( 
        mapping, dof_handler, constraints_without_dbc, quad, data); 

      matrix_free.initialize_dof_vector(b); 
      matrix_free.initialize_dof_vector(x); 

      constraints.distribute(x); 

      matrix_free.cell_loop(&LaplaceOperator::do_cell_integral_range, 
                            this, 
                            b, 
                            x); 

      constraints.set_zero(b); 

      system_rhs -= b; 
    } 
  } 



  template <int dim, typename number> 
  types::global_dof_index LaplaceOperator<dim, number>::m() const 
  { 
    return matrix_free.get_dof_handler().n_dofs(); 
  } 


  template <int dim, typename number> 
  number LaplaceOperator<dim, number>::el(unsigned int, unsigned int) const 
  { 
    Assert(false, ExcNotImplemented()); 
    return 0; 
  } 


  template <int dim, typename number> 
  void 
  LaplaceOperator<dim, number>::initialize_dof_vector(VectorType &vec) const 
  { 
    matrix_free.initialize_dof_vector(vec); 
  } 


  template <int dim, typename number> 
  void LaplaceOperator<dim, number>::vmult(VectorType &      dst, 
                                           const VectorType &src) const 
  { 
    this->matrix_free.cell_loop( 
      &LaplaceOperator::do_cell_integral_range, this, dst, src, true); 
  } 


  template <int dim, typename number> 
  void LaplaceOperator<dim, number>::Tvmult(VectorType &      dst, 
                                            const VectorType &src) const 
  { 
    this->vmult(dst, src); 
  } 


  template <int dim, typename number> 
  void LaplaceOperator<dim, number>::compute_inverse_diagonal( 
    VectorType &diagonal) const 
  { 
    MatrixFreeTools::compute_diagonal(matrix_free, 
                                      diagonal, 
                                      &LaplaceOperator::do_cell_integral_local, 
                                      this); 

    for (auto &i : diagonal) 
      i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0; 
  } 


  template <int dim, typename number> 
  const TrilinosWrappers::SparseMatrix & 
  LaplaceOperator<dim, number>::get_system_matrix() const 
  { 
    if (system_matrix.m() == 0 && system_matrix.n() == 0) 
      { 
        const auto &dof_handler = this->matrix_free.get_dof_handler(); 

        TrilinosWrappers::SparsityPattern dsp( 
          dof_handler.locally_owned_dofs(), 
          dof_handler.get_triangulation().get_communicator()); 

        DoFTools::make_sparsity_pattern(dof_handler, dsp, this->constraints); 

        dsp.compress(); 
        system_matrix.reinit(dsp); 

        MatrixFreeTools::compute_matrix( 
          matrix_free, 
          constraints, 
          system_matrix, 
          &LaplaceOperator::do_cell_integral_local, 
          this); 
      } 

    return this->system_matrix; 
  } 


  template <int dim, typename number> 
  void LaplaceOperator<dim, number>::do_cell_integral_local( 
    FECellIntegrator &integrator) const 
  { 
    integrator.evaluate(EvaluationFlags::gradients); 

    for (unsigned int q = 0; q < integrator.n_q_points; ++q) 
      integrator.submit_gradient(integrator.get_gradient(q), q); 

    integrator.integrate(EvaluationFlags::gradients); 
  } 


  template <int dim, typename number> 
  void LaplaceOperator<dim, number>::do_cell_integral_global( 
    FECellIntegrator &integrator, 
    VectorType &      dst, 
    const VectorType &src) const 
  { 
    integrator.gather_evaluate(src, EvaluationFlags::gradients); 

    for (unsigned int q = 0; q < integrator.n_q_points; ++q) 
      integrator.submit_gradient(integrator.get_gradient(q), q); 

    integrator.integrate_scatter(EvaluationFlags::gradients, dst); 
  } 


  template <int dim, typename number> 
  void LaplaceOperator<dim, number>::do_cell_integral_range( 
    const MatrixFree<dim, number> &              matrix_free, 
    VectorType &                                 dst, 
    const VectorType &                           src, 
    const std::pair<unsigned int, unsigned int> &range) const 
  { 
    FECellIntegrator integrator(matrix_free, range); 

    for (unsigned cell = range.first; cell < range.second; ++cell) 
      { 
        integrator.reinit(cell); 

        do_cell_integral_global(integrator, dst, src); 
      } 
  } 



  template <typename VectorType, 
            int dim, 
            typename SystemMatrixType, 
            typename LevelMatrixType, 
            typename MGTransferType> 
  static void 
  mg_solve(SolverControl &            solver_control, 
           VectorType &               dst, 
           const VectorType &         src, 
           const MultigridParameters &mg_data, 
           const DoFHandler<dim> &    dof, 
           const SystemMatrixType &   fine_matrix, 
           const MGLevelObject<std::unique_ptr<LevelMatrixType>> &mg_matrices, 
           const MGTransferType &                                 mg_transfer) 
  { 
    AssertThrow(mg_data.coarse_solver.type == "cg_with_amg", 
                ExcNotImplemented()); 
    AssertThrow(mg_data.smoother.type == "chebyshev", ExcNotImplemented()); 

    const unsigned int min_level = mg_matrices.min_level(); 
    const unsigned int max_level = mg_matrices.max_level(); 

    using SmootherPreconditionerType = DiagonalMatrix<VectorType>; 
    using SmootherType               = PreconditionChebyshev<LevelMatrixType, 
                                               VectorType, 
                                               SmootherPreconditionerType>; 
    using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>; 


    mg::Matrix<VectorType> mg_matrix(mg_matrices); 

    MGLevelObject<typename SmootherType::AdditionalData> smoother_data( 
      min_level, max_level); 

    for (unsigned int level = min_level; level <= max_level; level++) 
      { 
        smoother_data[level].preconditioner = 
          std::make_shared<SmootherPreconditionerType>(); 
        mg_matrices[level]->compute_inverse_diagonal( 
          smoother_data[level].preconditioner->get_vector()); 
        smoother_data[level].smoothing_range = mg_data.smoother.smoothing_range; 
        smoother_data[level].degree          = mg_data.smoother.degree; 
        smoother_data[level].eig_cg_n_iterations = 
          mg_data.smoother.eig_cg_n_iterations; 
      } 

    MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType> 
      mg_smoother; 
    mg_smoother.initialize(mg_matrices, smoother_data); 


    ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter, 
                                                mg_data.coarse_solver.abstol, 
                                                mg_data.coarse_solver.reltol, 
                                                false, 
                                                false); 
    SolverCG<VectorType> coarse_grid_solver(coarse_grid_solver_control); 

    std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse; 

    TrilinosWrappers::PreconditionAMG                 precondition_amg; 
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_data; 
    amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps; 
    amg_data.n_cycles        = mg_data.coarse_solver.n_cycles; 
    amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str(); 

    precondition_amg.initialize(mg_matrices[min_level]->get_system_matrix(), 
                                amg_data); 

    mg_coarse = 
      std::make_unique<MGCoarseGridIterativeSolver<VectorType, 
                                                   SolverCG<VectorType>, 
                                                   LevelMatrixType, 
                                                   decltype(precondition_amg)>>( 
        coarse_grid_solver, *mg_matrices[min_level], precondition_amg); 


    Multigrid<VectorType> mg( 
      mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother); 

    PreconditionerType preconditioner(dof, mg, mg_transfer); 

    SolverCG<VectorType>(solver_control) 
      .solve(fine_matrix, dst, src, preconditioner); 
  } 



  template <typename VectorType, typename OperatorType, int dim> 
  void solve_with_gmg(SolverControl &                  solver_control, 
                      const OperatorType &             system_matrix, 
                      VectorType &                     dst, 
                      const VectorType &               src, 
                      const MultigridParameters &      mg_data, 
                      const hp::MappingCollection<dim> mapping_collection, 
                      const DoFHandler<dim> &          dof_handler, 
                      const hp::QCollection<dim> &     quadrature_collection) 
  { 



    MGLevelObject<DoFHandler<dim>>                     dof_handlers; 
    MGLevelObject<std::unique_ptr<OperatorType>>       operators; 
    MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers; 

    std::vector<std::shared_ptr<const Triangulation<dim>>> 
      coarse_grid_triangulations; 
    if (mg_data.transfer.perform_h_transfer) 
      coarse_grid_triangulations = 
        MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence( 
          dof_handler.get_triangulation()); 
    else 
      coarse_grid_triangulations.emplace_back( 
        const_cast<Triangulation<dim> *>(&(dof_handler.get_triangulation())), 
        [](auto &) {}); 


    const unsigned int n_h_levels = coarse_grid_triangulations.size() - 1; 

    const auto get_max_active_fe_degree = [&](const auto &dof_handler) { 
      unsigned int max = 0; 

      for (auto &cell : dof_handler.active_cell_iterators()) 
        if (cell->is_locally_owned()) 
          max = 
            std::max(max, dof_handler.get_fe(cell->active_fe_index()).degree); 

      return Utilities::MPI::max(max, MPI_COMM_WORLD); 
    }; 

    const unsigned int n_p_levels = 
      MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence( 
        get_max_active_fe_degree(dof_handler), mg_data.transfer.p_sequence) 
        .size(); 

    std::map<unsigned int, unsigned int> fe_index_for_degree; 
    for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i) 
      { 
        const unsigned int degree = dof_handler.get_fe(i).degree; 
        Assert(fe_index_for_degree.find(degree) == fe_index_for_degree.end(), 
               ExcMessage("FECollection does not contain unique degrees.")); 
        fe_index_for_degree[degree] = i; 
      } 

    unsigned int minlevel   = 0; 
    unsigned int minlevel_p = n_h_levels; 
    unsigned int maxlevel   = n_h_levels + n_p_levels - 1; 

    dof_handlers.resize(minlevel, maxlevel); 
    operators.resize(minlevel, maxlevel); 
    transfers.resize(minlevel, maxlevel); 


    for (unsigned int l = 0; l < n_h_levels; ++l) 
      { 
        dof_handlers[l].reinit(*coarse_grid_triangulations[l]); 
        dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection()); 
      } 


    for (unsigned int i = 0, l = maxlevel; i < n_p_levels; ++i, --l) 
      { 
        dof_handlers[l].reinit(dof_handler.get_triangulation()); 

        if (l == maxlevel) // finest level 
          { 
            auto &dof_handler_mg = dof_handlers[l]; 

            auto cell_other = dof_handler.begin_active(); 
            for (auto &cell : dof_handler_mg.active_cell_iterators()) 
              { 
                if (cell->is_locally_owned()) 
                  cell->set_active_fe_index(cell_other->active_fe_index()); 
                cell_other++; 
              } 
          } 
        else // coarse level 
          { 
            auto &dof_handler_fine   = dof_handlers[l + 1]; 
            auto &dof_handler_coarse = dof_handlers[l + 0]; 

            auto cell_other = dof_handler_fine.begin_active(); 
            for (auto &cell : dof_handler_coarse.active_cell_iterators()) 
              { 
                if (cell->is_locally_owned()) 
                  { 
                    const unsigned int next_degree = 
                      MGTransferGlobalCoarseningTools:: 
                        create_next_polynomial_coarsening_degree( 
                          cell_other->get_fe().degree, 
                          mg_data.transfer.p_sequence); 
                    Assert(fe_index_for_degree.find(next_degree) != 
                             fe_index_for_degree.end(), 
                           ExcMessage("Next polynomial degree in sequence " 
                                      "does not exist in FECollection.")); 

                    cell->set_active_fe_index(fe_index_for_degree[next_degree]); 
                  } 
                cell_other++; 
              } 
          } 

        dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection()); 
      } 


    MGLevelObject<AffineConstraints<typename VectorType::value_type>> 
      constraints(minlevel, maxlevel); 

    for (unsigned int level = minlevel; level <= maxlevel; ++level) 
      { 
        const auto &dof_handler = dof_handlers[level]; 
        auto &      constraint  = constraints[level]; 

        IndexSet locally_relevant_dofs; 
        DoFTools::extract_locally_relevant_dofs(dof_handler, 
                                                locally_relevant_dofs); 
        constraint.reinit(locally_relevant_dofs); 

        DoFTools::make_hanging_node_constraints(dof_handler, constraint); 
        VectorTools::interpolate_boundary_values(mapping_collection, 
                                                 dof_handler, 
                                                 0, 
                                                 Functions::ZeroFunction<dim>(), 
                                                 constraint); 
        constraint.close(); 

        VectorType dummy; 

        operators[level] = std::make_unique<OperatorType>(mapping_collection, 
                                                          dof_handler, 
                                                          quadrature_collection, 
                                                          constraint, 
                                                          dummy); 
      } 


    for (unsigned int level = minlevel; level < minlevel_p; ++level) 
      transfers[level + 1].reinit_geometric_transfer(dof_handlers[level + 1], 
                                                     dof_handlers[level], 
                                                     constraints[level + 1], 
                                                     constraints[level]); 

    for (unsigned int level = minlevel_p; level < maxlevel; ++level) 
      transfers[level + 1].reinit_polynomial_transfer(dof_handlers[level + 1], 
                                                      dof_handlers[level], 
                                                      constraints[level + 1], 
                                                      constraints[level]); 

    MGTransferGlobalCoarsening<dim, VectorType> transfer( 
      transfers, [&](const auto l, auto &vec) { 
        operators[l]->initialize_dof_vector(vec); 
      }); 


    mg_solve(solver_control, 
             dst, 
             src, 
             mg_data, 
             dof_handler, 
             system_matrix, 
             operators, 
             transfer); 
  } 





  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(const Parameters &parameters); 

    void run(); 

  private: 
    void initialize_grid(); 
    void setup_system(); 
    void print_diagnostics(); 
    void solve_system(); 
    void compute_indicators(); 
    void adapt_resolution(); 
    void output_results(const unsigned int cycle); 

    MPI_Comm mpi_communicator; 

    const Parameters prm; 

    parallel::distributed::Triangulation<dim> triangulation; 
    DoFHandler<dim>                           dof_handler; 

    hp::MappingCollection<dim> mapping_collection; 
    hp::FECollection<dim>      fe_collection; 
    hp::QCollection<dim>       quadrature_collection; 
    hp::QCollection<dim - 1>   face_quadrature_collection; 

    IndexSet locally_owned_dofs; 
    IndexSet locally_relevant_dofs; 

    AffineConstraints<double> constraints; 

    LaplaceOperator<dim, double>               laplace_operator; 
    LinearAlgebra::distributed::Vector<double> locally_relevant_solution; 
    LinearAlgebra::distributed::Vector<double> system_rhs; 

    std::unique_ptr<FESeries::Legendre<dim>>    legendre; 
    std::unique_ptr<parallel::CellWeights<dim>> cell_weights; 

    Vector<float> estimated_error_per_cell; 
    Vector<float> hp_decision_indicators; 

    ConditionalOStream pcout; 
    TimerOutput        computing_timer; 
  }; 



  template <int dim> 
  LaplaceProblem<dim>::LaplaceProblem(const Parameters &parameters) 
    : mpi_communicator(MPI_COMM_WORLD) 
    , prm(parameters) 
    , triangulation(mpi_communicator) 
    , dof_handler(triangulation) 
    , pcout(std::cout, 
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
    , computing_timer(mpi_communicator, 
                      pcout, 
                      TimerOutput::summary, 
                      TimerOutput::wall_times) 
  { 
    Assert(prm.min_h_level <= prm.max_h_level, 
           ExcMessage( 
             "Triangulation level limits have been incorrectly set up.")); 
    Assert(prm.min_p_degree <= prm.max_p_degree, 
           ExcMessage("FECollection degrees have been incorrectly set up.")); 



    mapping_collection.push_back(MappingQ1<dim>()); 

    for (unsigned int degree = 1; degree <= prm.max_p_degree; ++degree) 
      { 
        fe_collection.push_back(FE_Q<dim>(degree)); 
        quadrature_collection.push_back(QGauss<dim>(degree + 1)); 
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1)); 
      } 


    const unsigned int min_fe_index = prm.min_p_degree - 1; 
    fe_collection.set_hierarchy( 

   /*下一个_index=  */ 
      [](const typename hp::FECollection<dim> &fe_collection, 
         const unsigned int                    fe_index) -> unsigned int { 
        return ((fe_index + 1) < fe_collection.size()) ? fe_index + 1 : 
                                                         fe_index; 
      }, 
    /*上一页_index=  */ 
      [min_fe_index](const typename hp::FECollection<dim> &, 
                     const unsigned int fe_index) -> unsigned int { 
        Assert(fe_index >= min_fe_index, 
               ExcMessage("Finite element is not part of hierarchy!")); 
        return (fe_index > min_fe_index) ? fe_index - 1 : fe_index; 
      }); 


    legendre = std::make_unique<FESeries::Legendre<dim>>( 
      SmoothnessEstimator::Legendre::default_fe_series(fe_collection)); 




    cell_weights = std::make_unique<parallel::CellWeights<dim>>( 
      dof_handler, 
      parallel::CellWeights<dim>::ndofs_weighting( 
        {prm.weighting_factor, prm.weighting_exponent})); 


    triangulation.signals.post_p4est_refinement.connect( 
      [&, min_fe_index]() { 
        const parallel::distributed::TemporarilyMatchRefineFlags<dim> 
          refine_modifier(triangulation); 
        hp::Refinement::limit_p_level_difference(dof_handler, 
                                                 prm.max_p_level_difference, 
                                                 /*包含=  */ min_fe_index);
      }, 
      boost::signals2::at_front); 
  } 






  template <int dim> 
  void LaplaceProblem<dim>::initialize_grid() 
  { 
    TimerOutput::Scope t(computing_timer, "initialize grid"); 

    std::vector<unsigned int> repetitions(dim); 
    Point<dim>                bottom_left, top_right; 
    for (unsigned int d = 0; d < dim; ++d) 
      if (d < 2) 
        { 
          repetitions[d] = 2; 
          bottom_left[d] = -1.; 
          top_right[d]   = 1.; 
        } 
      else 
        { 
          repetitions[d] = 1; 
          bottom_left[d] = 0.; 
          top_right[d]   = 1.; 
        } 

    std::vector<int> cells_to_remove(dim, 1); 
    cells_to_remove[0] = -1; 

    GridGenerator::subdivided_hyper_L( 
      triangulation, repetitions, bottom_left, top_right, cells_to_remove); 

    triangulation.refine_global(prm.min_h_level); 

    const unsigned int min_fe_index = prm.min_p_degree - 1; 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        cell->set_active_fe_index(min_fe_index); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::setup_system() 
  { 
    TimerOutput::Scope t(computing_timer, "setup system"); 

    dof_handler.distribute_dofs(fe_collection); 

    locally_owned_dofs = dof_handler.locally_owned_dofs(); 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 

    locally_relevant_solution.reinit(locally_owned_dofs, 
                                     locally_relevant_dofs, 
                                     mpi_communicator); 
    system_rhs.reinit(locally_owned_dofs, mpi_communicator); 

    constraints.clear(); 
    constraints.reinit(locally_relevant_dofs); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    VectorTools::interpolate_boundary_values( 
      mapping_collection, dof_handler, 0, Solution<dim>(), constraints); 
    constraints.close(); 

    laplace_operator.reinit(mapping_collection, 
                            dof_handler, 
                            quadrature_collection, 
                            constraints, 
                            system_rhs); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::print_diagnostics() 
  { 
    const unsigned int first_n_processes = 
      std::min<unsigned int>(8, 
                             Utilities::MPI::n_mpi_processes(mpi_communicator)); 
    const bool output_cropped = 
      first_n_processes < Utilities::MPI::n_mpi_processes(mpi_communicator); 

    { 
      pcout << "   Number of active cells:       " 
            << triangulation.n_global_active_cells() << std::endl 
            << "     by partition:              "; 

      std::vector<unsigned int> n_active_cells_per_subdomain = 
        Utilities::MPI::gather(mpi_communicator, 
                               triangulation.n_locally_owned_active_cells()); 
      for (unsigned int i = 0; i < first_n_processes; ++i) 
        pcout << ' ' << n_active_cells_per_subdomain[i]; 
      if (output_cropped) 
        pcout << " ..."; 
      pcout << std::endl; 
    } 

    { 
      pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
            << std::endl 
            << "     by partition:              "; 

      std::vector<types::global_dof_index> n_dofs_per_subdomain = 
        Utilities::MPI::gather(mpi_communicator, 
                               dof_handler.n_locally_owned_dofs()); 
      for (unsigned int i = 0; i < first_n_processes; ++i) 
        pcout << ' ' << n_dofs_per_subdomain[i]; 
      if (output_cropped) 
        pcout << " ..."; 
      pcout << std::endl; 
    } 

    { 
      std::vector<types::global_dof_index> n_constraints_per_subdomain = 
        Utilities::MPI::gather(mpi_communicator, constraints.n_constraints()); 

      pcout << "   Number of constraints:        " 
            << std::accumulate(n_constraints_per_subdomain.begin(), 
                               n_constraints_per_subdomain.end(), 
                               0) 
            << std::endl 
            << "     by partition:              "; 
      for (unsigned int i = 0; i < first_n_processes; ++i) 
        pcout << ' ' << n_constraints_per_subdomain[i]; 
      if (output_cropped) 
        pcout << " ..."; 
      pcout << std::endl; 
    } 

    { 
      std::vector<unsigned int> n_fe_indices(fe_collection.size(), 0); 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        if (cell->is_locally_owned()) 
          n_fe_indices[cell->active_fe_index()]++; 

      Utilities::MPI::sum(n_fe_indices, mpi_communicator, n_fe_indices); 

      pcout << "   Frequencies of poly. degrees:"; 
      for (unsigned int i = 0; i < fe_collection.size(); ++i) 
        if (n_fe_indices[i] > 0) 
          pcout << ' ' << fe_collection[i].degree << ":" << n_fe_indices[i]; 
      pcout << std::endl; 
    } 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::solve_system() 
  { 
    TimerOutput::Scope t(computing_timer, "solve system"); 

    LinearAlgebra::distributed::Vector<double> completely_distributed_solution; 
    laplace_operator.initialize_dof_vector(completely_distributed_solution); 

    SolverControl solver_control(system_rhs.size(), 
                                 prm.tolerance_factor * system_rhs.l2_norm()); 

    solve_with_gmg(solver_control, 
                   laplace_operator, 
                   completely_distributed_solution, 
                   system_rhs, 
                   prm.mg_data, 
                   mapping_collection, 
                   dof_handler, 
                   quadrature_collection); 

    pcout << "   Solved in " << solver_control.last_step() << " iterations." 
          << std::endl; 

    constraints.distribute(completely_distributed_solution); 

    locally_relevant_solution.copy_locally_owned_data_from( 
      completely_distributed_solution); 
    locally_relevant_solution.update_ghost_values(); 
  } 





  template <int dim> 
  void LaplaceProblem<dim>::compute_indicators() 
  { 
    TimerOutput::Scope t(computing_timer, "compute indicators"); 

    estimated_error_per_cell.grow_or_shrink(triangulation.n_active_cells()); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      face_quadrature_collection, 
      std::map<types::boundary_id, const Function<dim> *>(), 
      locally_relevant_solution, 
      estimated_error_per_cell, 
      /*component_mask=  */ 
      ComponentMask(), 
      /*coefficients=  */ 
      nullptr, 
      /*n_threads=*/
      numbers::invalid_unsigned_int,  
      /*subdomain_id=*/ 
      numbers::invalid_subdomain_id,  
      /*material_id=*/ 
      numbers::invalid_material_id,  
      /*策略=  */ 
      KellyErrorEstimator<dim>::Strategy::face_diameter_over_twice_max_degree); 

    hp_decision_indicators.grow_or_shrink(triangulation.n_active_cells()); 
    SmoothnessEstimator::Legendre::coefficient_decay(*legendre, 
                                                     dof_handler, 
                                                     locally_relevant_solution, 
                                                     hp_decision_indicators); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::adapt_resolution() 
  { 
    TimerOutput::Scope t(computing_timer, "adapt resolution"); 



    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number( 
      triangulation, 
      estimated_error_per_cell, 
      prm.refine_fraction, 
      prm.coarsen_fraction); 




    hp::Refinement::p_adaptivity_fixed_number(dof_handler, 
                                              hp_decision_indicators, 
                                              prm.p_refine_fraction, 
                                              prm.p_coarsen_fraction); 


    hp::Refinement::choose_p_over_h(dof_handler); 



    Assert(triangulation.n_levels() >= prm.min_h_level + 1 && 
             triangulation.n_levels() <= prm.max_h_level + 1, 
           ExcInternalError()); 

    if (triangulation.n_levels() > prm.max_h_level) 
      for (const auto &cell : 
           triangulation.active_cell_iterators_on_level(prm.max_h_level)) 
        cell->clear_refine_flag(); 

    for (const auto &cell : 
         triangulation.active_cell_iterators_on_level(prm.min_h_level)) 
      cell->clear_coarsen_flag(); 



    triangulation.execute_coarsening_and_refinement(); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::output_results(const unsigned int cycle) 
  { 
    TimerOutput::Scope t(computing_timer, "output results"); 

    Vector<float> fe_degrees(triangulation.n_active_cells()); 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        fe_degrees(cell->active_cell_index()) = cell->get_fe().degree; 

    Vector<float> subdomain(triangulation.n_active_cells()); 
    for (auto &subd : subdomain) 
      subd = triangulation.locally_owned_subdomain(); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(locally_relevant_solution, "solution"); 
    data_out.add_data_vector(fe_degrees, "fe_degree"); 
    data_out.add_data_vector(subdomain, "subdomain"); 
    data_out.add_data_vector(estimated_error_per_cell, "error"); 
    data_out.add_data_vector(hp_decision_indicators, "hp_indicator"); 
    data_out.build_patches(mapping_collection); 

    data_out.write_vtu_with_pvtu_record( 
      "./", "solution", cycle, mpi_communicator, 2, 1); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::run() 
  { 
    pcout << "Running with Trilinos on " 
          << Utilities::MPI::n_mpi_processes(mpi_communicator) 
          << " MPI rank(s)..." << std::endl; 

    { 
      pcout << "Calculating transformation matrices..." << std::endl; 
      TimerOutput::Scope t(computing_timer, "calculate transformation"); 
      legendre->precalculate_all_transformation_matrices(); 
    } 

    for (unsigned int cycle = 0; cycle < prm.n_cycles; ++cycle) 
      { 
        pcout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          initialize_grid(); 
        else 
          adapt_resolution(); 

        setup_system(); 

        print_diagnostics(); 

        solve_system(); 

        compute_indicators(); 

        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32) 
          output_results(cycle); 

        computing_timer.print_summary(); 
        computing_timer.reset(); 

        pcout << std::endl; 
      } 
  } 
} // namespace Step75 



int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step75; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

      Parameters        prm; 
      LaplaceProblem<2> laplace_problem(prm); 
      laplace_problem.run(); 
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

