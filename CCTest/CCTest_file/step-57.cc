CCTest_file/step-57.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2008 - 2020 by the deal.II authors 
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
 * Author: Liang Zhao and Timo Heister, Clemson University, 2016 
 */ 


// @sect3{Include files}  

// 像往常一样，我们从包括一些著名的文件开始。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/tensor.h> 

#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_tools.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

// 为了在网格之间传输解决方案，包括这个文件。

#include <deal.II/numerics/solution_transfer.h> 

// 这个文件包括UMFPACK：直接求解器。

#include <deal.II/lac/sparse_direct.h> 

// 还有一个ILU预处理程序。

#include <deal.II/lac/sparse_ilu.h> 

#include <fstream> 
#include <iostream> 

namespace Step57 
{ 
  using namespace dealii; 
// @sect3{The <code>NavierStokesProblem</code> class template}  

// 该类管理介绍中描述的矩阵和向量：特别是，我们为当前的解决方案、当前的牛顿更新和直线搜索更新存储了一个BlockVector。 我们还存储了两个AffineConstraints对象：一个是强制执行Dirichlet边界条件的对象，另一个是将所有边界值设为0的对象。第一个约束解向量，第二个约束更新（也就是说，我们从不更新边界值，所以我们强制相关的更新向量值为零）。

  template <int dim> 
  class StationaryNavierStokes 
  { 
  public: 
    StationaryNavierStokes(const unsigned int degree); 
    void run(const unsigned int refinement); 

  private: 
    void setup_dofs(); 

    void initialize_system(); 

    void assemble(const bool initial_step, const bool assemble_matrix); 

    void assemble_system(const bool initial_step); 

    void assemble_rhs(const bool initial_step); 

    void solve(const bool initial_step); 

    void refine_mesh(); 

    void process_solution(unsigned int refinement); 

    void output_results(const unsigned int refinement_cycle) const; 

    void newton_iteration(const double       tolerance, 
                          const unsigned int max_n_line_searches, 
                          const unsigned int max_n_refinements, 
                          const bool         is_initial_step, 
                          const bool         output_result); 

    void compute_initial_guess(double step_size); 

    double                               viscosity; 
    double                               gamma; 
    const unsigned int                   degree; 
    std::vector<types::global_dof_index> dofs_per_block; 

    Triangulation<dim> triangulation; 
    FESystem<dim>      fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> zero_constraints; 
    AffineConstraints<double> nonzero_constraints; 

    BlockSparsityPattern      sparsity_pattern; 
    BlockSparseMatrix<double> system_matrix; 
    SparseMatrix<double>      pressure_mass_matrix; 

    BlockVector<double> present_solution; 
    BlockVector<double> newton_update; 
    BlockVector<double> system_rhs; 
    BlockVector<double> evaluation_point; 
  }; 
// @sect3{Boundary values and right hand side}  

// 在这个问题中，我们设定沿空腔上表面的速度为1，其他三面墙的速度为0。右边的函数为零，所以我们在本教程中不需要设置右边的函数。边界函数的分量数为  <code>dim+1</code>  。我们最终将使用 VectorTools::interpolate_boundary_values 来设置边界值，这就要求边界值函数的分量数与解相同，即使没有全部使用。换个说法：为了让这个函数高兴，我们为压力定义了边界值，尽管我们实际上永远不会用到它们。

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    BoundaryValues() 
      : Function<dim>(dim + 1) 
    {} 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> & p, 
                                    const unsigned int component) const 
  { 
    Assert(component < this->n_components, 
           ExcIndexRange(component, 0, this->n_components)); 
    if (component == 0 && std::abs(p[dim - 1] - 1.0) < 1e-10) 
      return 1.0; 

    return 0; 
  } 
// @sect3{BlockSchurPreconditioner for Navier Stokes equations}  

// 正如介绍中所讨论的，Krylov迭代方法中的预处理器是作为一个矩阵-向量乘积算子实现的。在实践中，舒尔补码预处理器被分解为三个矩阵的乘积（如第一节所述）。第一个因素中的 $\tilde{A}^{-1}$ 涉及到对线性系统 $\tilde{A}x=b$ 的求解。在这里，为了简单起见，我们通过一个直接求解器来解决这个系统。第二个因素中涉及的计算是一个简单的矩阵-向量乘法。舒尔补码 $\tilde{S}$ 可以被压力质量矩阵很好地近似，其逆值可以通过不精确求解器得到。因为压力质量矩阵是对称和正定的，我们可以用CG来解决相应的线性系统。

  template <class PreconditionerMp> 
  class BlockSchurPreconditioner : public Subscriptor 
  { 
  public: 
    BlockSchurPreconditioner(double                           gamma, 
                             double                           viscosity, 
                             const BlockSparseMatrix<double> &S, 
                             const SparseMatrix<double> &     P, 
                             const PreconditionerMp &         Mppreconditioner); 

    void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const; 

  private: 
    const double                     gamma; 
    const double                     viscosity; 
    const BlockSparseMatrix<double> &stokes_matrix; 
    const SparseMatrix<double> &     pressure_mass_matrix; 
    const PreconditionerMp &         mp_preconditioner; 
    SparseDirectUMFPACK              A_inverse; 
  }; 

// 我们可以注意到，左上角的矩阵逆的初始化是在构造函数中完成的。如果是这样，那么预处理程序的每一次应用就不再需要计算矩阵因子了。

  template <class PreconditionerMp> 
  BlockSchurPreconditioner<PreconditionerMp>::BlockSchurPreconditioner( 
    double                           gamma, 
    double                           viscosity, 
    const BlockSparseMatrix<double> &S, 
    const SparseMatrix<double> &     P, 
    const PreconditionerMp &         Mppreconditioner) 
    : gamma(gamma) 
    , viscosity(viscosity) 
    , stokes_matrix(S) 
    , pressure_mass_matrix(P) 
    , mp_preconditioner(Mppreconditioner) 
  { 
    A_inverse.initialize(stokes_matrix.block(0, 0)); 
  } 

  template <class PreconditionerMp> 
  void BlockSchurPreconditioner<PreconditionerMp>::vmult( 
    BlockVector<double> &      dst, 
    const BlockVector<double> &src) const 
  { 
    Vector<double> utmp(src.block(0)); 

    { 
      SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 

      dst.block(1) = 0.0; 
      cg.solve(pressure_mass_matrix, 
               dst.block(1), 
               src.block(1), 
               mp_preconditioner); 
      dst.block(1) *= -(viscosity + gamma); 
    } 

    { 
      stokes_matrix.block(0, 1).vmult(utmp, dst.block(1)); 
      utmp *= -1.0; 
      utmp += src.block(0); 
    } 

    A_inverse.vmult(dst.block(0), utmp); 
  } 
// @sect3{StationaryNavierStokes class implementation}  
// @sect4{StationaryNavierStokes::StationaryNavierStokes}  

// 该类的构造函数看起来与  step-22  中的构造函数非常相似。唯一的区别是粘度和增强的拉格朗日系数  <code>gamma</code>  。

  template <int dim> 
  StationaryNavierStokes<dim>::StationaryNavierStokes(const unsigned int degree) 
    : viscosity(1.0 / 7500.0) 
    , gamma(1.0) 
    , degree(degree) 
    , triangulation(Triangulation<dim>::maximum_smoothing) 
    , fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1) 
    , dof_handler(triangulation) 
  {} 
// @sect4{StationaryNavierStokes::setup_dofs}  

// 这个函数初始化DoFHandler，列举当前网格上的自由度和约束。

  template <int dim> 
  void StationaryNavierStokes<dim>::setup_dofs() 
  { 
    system_matrix.clear(); 
    pressure_mass_matrix.clear(); 

// 第一步是将DoFs与给定的网格联系起来。

    dof_handler.distribute_dofs(fe); 

// 我们对组件重新编号，使所有的速度DoF在压力DoF之前，以便能够将解向量分成两个块，在块预处理程序中分别访问。

    std::vector<unsigned int> block_component(dim + 1, 0); 
    block_component[dim] = 1; 
    DoFRenumbering::component_wise(dof_handler, block_component); 

    dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
    unsigned int dof_u = dofs_per_block[0]; 
    unsigned int dof_p = dofs_per_block[1]; 

// 在牛顿方案中，我们首先将边界条件应用于从初始步骤得到的解。为了确保边界条件在牛顿迭代过程中保持满足，在更新时使用零边界条件  $\delta u^k$  。因此我们设置了两个不同的约束对象。

    FEValuesExtractors::Vector velocities(0); 
    { 
      nonzero_constraints.clear(); 

      DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               0, 
                                               BoundaryValues<dim>(), 
                                               nonzero_constraints, 
                                               fe.component_mask(velocities)); 
    } 
    nonzero_constraints.close(); 

    { 
      zero_constraints.clear(); 

      DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               0, 
                                               Functions::ZeroFunction<dim>( 
                                                 dim + 1), 
                                               zero_constraints, 
                                               fe.component_mask(velocities)); 
    } 
    zero_constraints.close(); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (" << dof_u << " + " << dof_p << ')' << std::endl; 
  } 
// @sect4{StationaryNavierStokes::initialize_system}  

// 在每个网格上，SparsityPattern和线性系统的大小是不同的。这个函数在网格细化后初始化它们。

  template <int dim> 
  void StationaryNavierStokes<dim>::initialize_system() 
  { 
    { 
      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 
      DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints); 
      sparsity_pattern.copy_from(dsp); 
    } 

    system_matrix.reinit(sparsity_pattern); 

    present_solution.reinit(dofs_per_block); 
    newton_update.reinit(dofs_per_block); 
    system_rhs.reinit(dofs_per_block); 
  } 
// @sect4{StationaryNavierStokes::assemble}  

// 这个函数建立了我们目前工作的系统矩阵和右手边。 @p initial_step 参数用于确定我们应用哪一组约束（初始步骤为非零，其他为零）。 @p assemble_matrix 参数分别决定了是组装整个系统还是只组装右手边的向量。

  template <int dim> 
  void StationaryNavierStokes<dim>::assemble(const bool initial_step, 
                                             const bool assemble_matrix) 
  { 
    if (assemble_matrix) 
      system_matrix = 0; 

    system_rhs = 0; 

    QGauss<dim> quadrature_formula(degree + 2); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values | update_gradients); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 对于线性化系统，我们为当前速度和梯度以及当前压力创建临时存储。在实践中，它们都是通过正交点的形状函数获得的。

    std::vector<Tensor<1, dim>> present_velocity_values(n_q_points); 
    std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points); 
    std::vector<double>         present_pressure_values(n_q_points); 

    std::vector<double>         div_phi_u(dofs_per_cell); 
    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell); 
    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell); 
    std::vector<double>         phi_p(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 

        local_matrix = 0; 
        local_rhs    = 0; 

        fe_values[velocities].get_function_values(evaluation_point, 
                                                  present_velocity_values); 

        fe_values[velocities].get_function_gradients( 
          evaluation_point, present_velocity_gradients); 

        fe_values[pressure].get_function_values(evaluation_point, 
                                                present_pressure_values); 

//装配类似于  step-22  。一个以gamma为系数的附加项是增强拉格朗日（AL），它是通过grad-div稳定化组装的。 正如我们在介绍中所讨论的，系统矩阵的右下块应该为零。由于压力质量矩阵是在创建预处理程序时使用的，所以我们在这里组装它，然后在最后把它移到一个单独的SparseMatrix中（与 step-22 相同）。

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                div_phi_u[k]  = fe_values[velocities].divergence(k, q); 
                grad_phi_u[k] = fe_values[velocities].gradient(k, q); 
                phi_u[k]      = fe_values[velocities].value(k, q); 
                phi_p[k]      = fe_values[pressure].value(k, q); 
              } 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              { 
                if (assemble_matrix) 
                  { 
                    for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                      { 
                        local_matrix(i, j) += 
                          (viscosity * 
                             scalar_product(grad_phi_u[j], grad_phi_u[i]) + 
                           present_velocity_gradients[q] * phi_u[j] * phi_u[i] + 
                           grad_phi_u[j] * present_velocity_values[q] * 
                             phi_u[i] - 
                           div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] + 
                           gamma * div_phi_u[j] * div_phi_u[i] + 
                           phi_p[i] * phi_p[j]) * 
                          fe_values.JxW(q); 
                      } 
                  } 

                double present_velocity_divergence = 
                  trace(present_velocity_gradients[q]); 
                local_rhs(i) += 
                  (-viscosity * scalar_product(present_velocity_gradients[q], 
                                               grad_phi_u[i]) - 
                   present_velocity_gradients[q] * present_velocity_values[q] * 
                     phi_u[i] + 
                   present_pressure_values[q] * div_phi_u[i] + 
                   present_velocity_divergence * phi_p[i] - 
                   gamma * present_velocity_divergence * div_phi_u[i]) * 
                  fe_values.JxW(q); 
              } 
          } 

        cell->get_dof_indices(local_dof_indices); 

        const AffineConstraints<double> &constraints_used = 
          initial_step ? nonzero_constraints : zero_constraints; 

        if (assemble_matrix) 
          { 
            constraints_used.distribute_local_to_global(local_matrix, 
                                                        local_rhs, 
                                                        local_dof_indices, 
                                                        system_matrix, 
                                                        system_rhs); 
          } 
        else 
          { 
            constraints_used.distribute_local_to_global(local_rhs, 
                                                        local_dof_indices, 
                                                        system_rhs); 
          } 
      } 

    if (assemble_matrix) 
      { 

// 最后我们把压力质量矩阵移到一个单独的矩阵中。

        pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1)); 
        pressure_mass_matrix.copy_from(system_matrix.block(1, 1)); 

// 注意，将这个压力块设置为零并不等同于不在这个块中装配任何东西，因为这里的操作将（错误地）删除从压力作用力的悬挂节点约束中进来的对角线条目。这意味着，我们的整个系统矩阵将有完全为零的行。幸运的是，FGMRES处理这些行没有任何问题。

        system_matrix.block(1, 1) = 0; 
      } 
  } 

  template <int dim> 
  void StationaryNavierStokes<dim>::assemble_system(const bool initial_step) 
  { 
    assemble(initial_step, true); 
  } 

  template <int dim> 
  void StationaryNavierStokes<dim>::assemble_rhs(const bool initial_step) 
  { 
    assemble(initial_step, false); 
  } 
// @sect4{StationaryNavierStokes::solve}  

// 在这个函数中，我们使用FGMRES和程序开始时定义的块状预处理程序来解决线性系统。我们在这一步得到的是解向量。如果这是初始步骤，解向量为我们提供了纳维尔-斯托克斯方程的初始猜测。对于初始步骤，非零约束被应用，以确保边界条件得到满足。在下面的步骤中，我们将求解牛顿更新，所以使用零约束。

  template <int dim> 
  void StationaryNavierStokes<dim>::solve(const bool initial_step) 
  { 
    const AffineConstraints<double> &constraints_used = 
      initial_step ? nonzero_constraints : zero_constraints; 

    SolverControl solver_control(system_matrix.m(), 
                                 1e-4 * system_rhs.l2_norm(), 
                                 true); 

    SolverFGMRES<BlockVector<double>> gmres(solver_control); 
    SparseILU<double>                 pmass_preconditioner; 
    pmass_preconditioner.initialize(pressure_mass_matrix, 
                                    SparseILU<double>::AdditionalData()); 

    const BlockSchurPreconditioner<SparseILU<double>> preconditioner( 
      gamma, 
      viscosity, 
      system_matrix, 
      pressure_mass_matrix, 
      pmass_preconditioner); 

    gmres.solve(system_matrix, newton_update, system_rhs, preconditioner); 
    std::cout << "FGMRES steps: " << solver_control.last_step() << std::endl; 

    constraints_used.distribute(newton_update); 
  } 
// @sect4{StationaryNavierStokes::refine_mesh}  

// 在粗略的网格上找到一个好的初始猜测后，我们希望通过细化网格来减少误差。这里我们做了类似于 step-15 的自适应细化，只是我们只使用了速度上的Kelly估计器。我们还需要使用SolutionTransfer类将当前的解转移到下一个网格。

  template <int dim> 
  void StationaryNavierStokes<dim>::refine_mesh() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
    FEValuesExtractors::Vector velocity(0); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      present_solution, 
      estimated_error_per_cell, 
      fe.component_mask(velocity)); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.0); 

    triangulation.prepare_coarsening_and_refinement(); 
    SolutionTransfer<dim, BlockVector<double>> solution_transfer(dof_handler); 
    solution_transfer.prepare_for_coarsening_and_refinement(present_solution); 
    triangulation.execute_coarsening_and_refinement(); 

// 首先，DoFHandler被设置，约束被生成。然后我们创建一个临时的BlockVector  <code>tmp</code>  ，其大小与新网格上的解决方案一致。

    setup_dofs(); 

    BlockVector<double> tmp(dofs_per_block); 

// 将解决方案从粗网格转移到细网格，并对新转移的解决方案应用边界值约束。注意，present_solution仍然是对应于旧网格的一个向量。

    solution_transfer.interpolate(present_solution, tmp); 
    nonzero_constraints.distribute(tmp); 

// 最后设置矩阵和向量，并将present_solution设置为插值后的数据。

    initialize_system(); 
    present_solution = tmp; 
  } 
// @sect4{StationaryNavierStokes<dim>::newton_iteration}  

// 这个函数实现了牛顿迭代，给定了公差、最大迭代次数和要做的网格细化次数。

// 参数 <code>is_initial_step</code> 告诉我们是否需要 <code>setup_system</code> ，以及应该装配哪一部分，系统矩阵或右手边的矢量。如果我们做直线搜索，在最后一次迭代中检查残差准则时，右手边已经被组装起来了。因此，我们只需要在当前迭代中装配系统矩阵。最后一个参数 <code>output_result</code> 决定了是否应该产生图形输出。

  template <int dim> 
  void StationaryNavierStokes<dim>::newton_iteration( 
    const double       tolerance, 
    const unsigned int max_n_line_searches, 
    const unsigned int max_n_refinements, 
    const bool         is_initial_step, 
    const bool         output_result) 
  { 
    bool first_step = is_initial_step; 

    for (unsigned int refinement_n = 0; refinement_n < max_n_refinements + 1; 
         ++refinement_n) 
      { 
        unsigned int line_search_n = 0; 
        double       last_res      = 1.0; 
        double       current_res   = 1.0; 
        std::cout << "grid refinements: " << refinement_n << std::endl 
                  << "viscosity: " << viscosity << std::endl; 

        while ((first_step || (current_res > tolerance)) && 
               line_search_n < max_n_line_searches) 
          { 
            if (first_step) 
              { 
                setup_dofs(); 
                initialize_system(); 
                evaluation_point = present_solution; 
                assemble_system(first_step); 
                solve(first_step); 
                present_solution = newton_update; 
                nonzero_constraints.distribute(present_solution); 
                first_step       = false; 
                evaluation_point = present_solution; 
                assemble_rhs(first_step); 
                current_res = system_rhs.l2_norm(); 
 
 
 
 
 
 
                evaluation_point = present_solution; 
                assemble_system(first_step); 
                solve(first_step); 

// 为了确保我们的解决方案越来越接近精确的解决方案，我们让解决方案用权重 <code>alpha</code> 更新，使新的残差小于上一步的残差，这是在下面的循环中完成。这与  step-15  中使用的线搜索算法相同。

                for (double alpha = 1.0; alpha > 1e-5; alpha *= 0.5) 
                  { 
                    evaluation_point = present_solution; 
                    evaluation_point.add(alpha, newton_update); 
                    nonzero_constraints.distribute(evaluation_point); 
                    assemble_rhs(first_step); 
                    current_res = system_rhs.l2_norm(); 
                    std::cout << "  alpha: " << std::setw(10) << alpha 
                              << std::setw(0) << "  residual: " << current_res 
                              << std::endl; 
                    if (current_res < last_res) 
                      break; 
                  } 
                { 
                  present_solution = evaluation_point; 
                  std::cout << "  number of line searches: " << line_search_n 
                            << "  residual: " << current_res << std::endl; 
                  last_res = current_res; 
                } 
                ++line_search_n; 
              } 

            if (output_result) 
              { 
                output_results(max_n_line_searches * refinement_n + 
                               line_search_n); 

                if (current_res <= tolerance) 
                  process_solution(refinement_n); 
              } 
          } 

        if (refinement_n < max_n_refinements) 
          { 
            refine_mesh(); 
          } 
      } 
  } 
// @sect4{StationaryNavierStokes::compute_initial_guess}  

// 这个函数将通过使用延续法为我们提供一个初始猜测，正如我们在介绍中讨论的那样。雷诺数被逐级增加 step- ，直到我们达到目标值。通过实验，斯托克斯的解足以成为雷诺数为1000的NSE的初始猜测，所以我们从这里开始。 为了确保前一个问题的解决方案与下一个问题足够接近，步长必须足够小。

  template <int dim> 
  void StationaryNavierStokes<dim>::compute_initial_guess(double step_size) 
  { 
    const double target_Re = 1.0 / viscosity; 

    bool is_initial_step = true; 

    for (double Re = 1000.0; Re < target_Re; 
         Re        = std::min(Re + step_size, target_Re)) 
      { 
        viscosity = 1.0 / Re; 
        std::cout << "Searching for initial guess with Re = " << Re 
                  << std::endl; 
        newton_iteration(1e-12, 50, 0, is_initial_step, false); 
        is_initial_step = false; 
      } 
  } 
// @sect4{StationaryNavierStokes::output_results}  

// 这个函数与 step-22 中的函数相同，只是我们为输出文件选择了一个同时包含雷诺数（即当前环境下的粘度的倒数）的名称。

  template <int dim> 
  void StationaryNavierStokes<dim>::output_results( 
    const unsigned int output_index) const 
  { 
    std::vector<std::string> solution_names(dim, "velocity"); 
    solution_names.emplace_back("pressure"); 

    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(present_solution, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.build_patches(); 

    std::ofstream output(std::to_string(1.0 / viscosity) + "-solution-" + 
                         Utilities::int_to_string(output_index, 4) + ".vtk"); 
    data_out.write_vtk(output); 
  } 
// @sect4{StationaryNavierStokes::process_solution}  

// 在我们的测试案例中，我们不知道分析解。该函数输出沿 $x=0.5$ 和 $0 \leq y \leq 1$ 的速度分量，以便与文献中的数据进行比较。

  template <int dim> 
  void StationaryNavierStokes<dim>::process_solution(unsigned int refinement) 
  { 
    std::ofstream f(std::to_string(1.0 / viscosity) + "-line-" + 
                    std::to_string(refinement) + ".txt"); 
    f << "# y u_x u_y" << std::endl; 

    Point<dim> p; 
    p(0) = 0.5; 
    p(1) = 0.5; 

    f << std::scientific; 

    for (unsigned int i = 0; i <= 100; ++i) 
      { 
        p(dim - 1) = i / 100.0; 

        Vector<double> tmp_vector(dim + 1); 
        VectorTools::point_value(dof_handler, present_solution, p, tmp_vector); 
        f << p(dim - 1); 

        for (int j = 0; j < dim; j++) 
          f << " " << tmp_vector(j); 
        f << std::endl; 
      } 
  } 
// @sect4{StationaryNavierStokes::run}  

// 这是本程序的最后一步。在这一部分，我们分别生成网格和运行其他函数。最大细化度可以通过参数来设置。

  template <int dim> 
  void StationaryNavierStokes<dim>::run(const unsigned int refinement) 
  { 
    GridGenerator::hyper_cube(triangulation); 
    triangulation.refine_global(5); 

    const double Re = 1.0 / viscosity; 

// 如果粘度小于 $1/1000$ ，我们必须首先通过延续法搜索初始猜测。我们应该注意的是，搜索总是在初始网格上进行的，也就是这个程序中的 $8 \times 8$ 网格。之后，我们只需做与粘度大于 $1/1000$ 时相同的工作：运行牛顿迭代，细化网格，转移解决方案，并重复。

    if (Re > 1000.0) 
      { 
        std::cout << "Searching for initial guess ..." << std::endl; 
        const double step_size = 2000.0; 
        compute_initial_guess(step_size); 
        std::cout << "Found initial guess." << std::endl; 
        std::cout << "Computing solution with target Re = " << Re << std::endl; 
        viscosity = 1.0 / Re; 
        newton_iteration(1e-12, 50, refinement, false, true); 
      } 
    else 
      { 

// 当粘度大于1/1000时，斯托克斯方程的解作为初始猜测已经足够好。如果是这样，我们就不需要用延续法来搜索初始猜测了。牛顿迭代可以直接开始。

        newton_iteration(1e-12, 50, refinement, true, true); 
      } 
  } 
} // namespace Step57 

int main() 
{ 
  try 
    { 
      using namespace Step57; 

      StationaryNavierStokes<2> flow(/* degree = */ 

1); 
      flow.run(4); 
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


