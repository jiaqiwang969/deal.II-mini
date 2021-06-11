CCTest_file/step-15.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2012 - 2021 by the deal.II authors 
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
 * Author: Sven Wetterauer, University of Heidelberg, 2012 
 */ 


// @sect3{Include files}  

// 前面几个文件已经在前面的例子中讲过了，因此不再做进一步的评论。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <fstream> 
#include <iostream> 

// 我们将在牛顿迭代之间使用自适应网格细化技术。要做到这一点，我们需要能够在新的网格上使用解决方案，尽管它是在旧的网格上计算出来的。SolutionTransfer类将解决方案从旧网格转移到新网格。

#include <deal.II/numerics/solution_transfer.h> 

// 然后，我们为这个程序打开一个命名空间，像以前的程序一样，将dealii命名空间中的所有东西导入其中。

namespace Step15 
{ 
  using namespace dealii; 
// @sect3{The <code>MinimalSurfaceProblem</code> class template}  

// 类模板与  step-6  中的基本相同。 增加了三个内容。

// - 有两个解决方案向量，一个用于牛顿更新  $\delta u^n$  ，另一个用于当前迭代  $u^n$  。

// -  <code>setup_system</code> 函数需要一个参数，表示这是否是第一次被调用。不同的是，第一次我们需要分配自由度，并将 $u^n$ 的解向量设置为正确的大小。接下来的几次，该函数是在我们已经完成了这些步骤，作为细化 <code>refine_mesh</code> 中网格的一部分之后被调用的。

// - 然后我们还需要新的函数。  <code>set_boundary_values()</code> 负责正确设置解向量的边界值，这在介绍的最后已经讨论过了。  <code>compute_residual()</code> 是一个计算非线性（离散）残差规范的函数。我们用这个函数来监测牛顿迭代的收敛性。该函数以步长 $\alpha^n$ 为参数来计算 $u^n + \alpha^n \; \delta u^n$ 的残差。这是人们通常需要的步长控制，尽管我们在这里不会使用这个功能。最后， <code>determine_step_length()</code> 计算每个牛顿迭代中的步长 $\alpha^n$ 。正如介绍中所讨论的，我们在这里使用一个固定的步长，并把实现一个更好的策略作为一个练习。(  step-77 的做法不同。它只是在整个求解过程中使用了一个外部包，而一个好的直线搜索策略是该包所提供的一部分）。)

  template <int dim> 
  class MinimalSurfaceProblem 
  { 
  public: 
    MinimalSurfaceProblem(); 
    void run(); 

  private: 
    void   setup_system(const bool initial_step); 
    void   assemble_system(); 
    void   solve(); 
    void   refine_mesh(); 
    void   set_boundary_values(); 
    double compute_residual(const double alpha) const; 
    double determine_step_length() const; 
    void   output_results(const unsigned int refinement_cycle) const; 

    Triangulation<dim> triangulation; 

    DoFHandler<dim> dof_handler; 
    FE_Q<dim>       fe; 

    AffineConstraints<double> hanging_node_constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> current_solution; 
    Vector<double> newton_update; 
    Vector<double> system_rhs; 
  }; 
// @sect3{Boundary condition}  

// 边界条件的实现就像在  step-4  中一样。 它被选为  $g(x,y)=\sin(2 \pi (x+y))$  。

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> &p, 
                                    const unsigned int /*component*/) const 
  { 
    return std::sin(2 * numbers::PI * (p[0] + p[1])); 
  } 
// @sect3{The <code>MinimalSurfaceProblem</code> class implementation}  
// @sect4{MinimalSurfaceProblem::MinimalSurfaceProblem}  

// 该类的构造函数和析构函数与前几篇教程中的相同。

  template <int dim> 
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem() 
    : dof_handler(triangulation) 
    , fe(2) 
  {} 
// @sect4{MinimalSurfaceProblem::setup_system}  

// 在setup-system函数中，我们总是设置有限元方法的变量。与 step-6 有相同的区别，因为在那里我们在每个细化周期中都要从头开始求解PDE，而在这里我们需要把以前的网格的解放到当前的网格上。因此，我们不能只是重置解向量。因此，传递给这个函数的参数表明我们是否可以分布自由度（加上计算约束）并将解向量设置为零，或者这在其他地方已经发生过了（特别是在 <code>refine_mesh()</code> ）。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step) 
  { 
    if (initial_step) 
      { 
        dof_handler.distribute_dofs(fe); 
        current_solution.reinit(dof_handler.n_dofs()); 

        hanging_node_constraints.clear(); 
        DoFTools::make_hanging_node_constraints(dof_handler, 
                                                hanging_node_constraints); 
        hanging_node_constraints.close(); 
      } 

// 该函数的其余部分与  step-6  中的相同。

    newton_update.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

    hanging_node_constraints.condense(dsp); 

    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 
  } 
// @sect4{MinimalSurfaceProblem::assemble_system}  

// 这个函数的作用与前面的教程相同，当然，现在矩阵和右手边的函数取决于上一次迭代的解。正如在介绍中所讨论的，我们需要使用牛顿更新的零边界值；我们在这个函数的最后计算它们。

// 该函数的顶部包含了通常的模板代码，设置了允许我们在正交点评估形状函数的对象，以及本地矩阵和向量的临时存储位置，以及正交点上先前解的梯度。然后我们开始在所有单元格上进行循环。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::assemble_system() 
  { 
    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    system_matrix = 0; 
    system_rhs    = 0; 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_gradients | update_quadrature_points | 
                              update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        fe_values.reinit(cell); 

// 为了组装线性系统，我们必须在正交点上获得前一个解的梯度值。有一个标准的方法： FEValues::get_function_gradients 函数接收一个代表定义在DoFHandler上的有限元场的向量，并评估这个场在FEValues对象最后被重新初始化的单元的正交点的梯度。然后将所有正交点的梯度值写入第二个参数中。

        fe_values.get_function_gradients(current_solution, 
                                         old_solution_gradients); 

// 有了这个，我们就可以对所有的正交点和形状函数进行积分循环。 在刚刚计算了正交点中旧解的梯度后，我们就可以计算这些点中的系数 $a_{n}$ 。 然后，系统本身的组装看起来与我们一贯的做法相似，除了非线性项之外，将结果从局部对象复制到全局对象中也是如此。

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double coeff = 
              1.0 / std::sqrt(1 + old_solution_gradients[q] * 
                                    old_solution_gradients[q]); 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              { 
                for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                  cell_matrix(i, j) += 
                    (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i 
                       * coeff                         //   * a_n 
                       * fe_values.shape_grad(j, q))   //   * \nabla \phi_j) 
                      -                                //  - 
                      (fe_values.shape_grad(i, q)      //  (\nabla \phi_i 
                       * coeff * coeff * coeff         //   * a_n^3 
                       * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j 
                          * old_solution_gradients[q]) //      * \nabla u_n) 
                       * old_solution_gradients[q]))   //   * \nabla u_n))) 
                     * fe_values.JxW(q));              // * dx 

                cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i 
                                * coeff                     // * a_n 
                                * old_solution_gradients[q] // * u_n 
                                * fe_values.JxW(q));        // * dx 
              } 
          } 

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          { 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              system_matrix.add(local_dof_indices[i], 
                                local_dof_indices[j], 
                                cell_matrix(i, j)); 

            system_rhs(local_dof_indices[i]) += cell_rhs(i); 
          } 
      } 

// 最后，我们从系统中移除悬挂的节点，并将零边界值应用到定义牛顿更新的线性系统中  $\delta u^n$  。

    hanging_node_constraints.condense(system_matrix); 
    hanging_node_constraints.condense(system_rhs); 

    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(), 
                                             boundary_values); 
    MatrixTools::apply_boundary_values(boundary_values, 
                                       system_matrix, 
                                       newton_update, 
                                       system_rhs); 
  } 

//  @sect4{MinimalSurfaceProblem::solve}  

// 解算函数和以往一样。在求解过程的最后，我们通过设置 $u^{n+1}=u^n+\alpha^n\;\delta u^n$ 来更新当前的解决方案。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::solve() 
  { 
    SolverControl            solver_control(system_rhs.size(), 
                                 system_rhs.l2_norm() * 1e-6); 
    SolverCG<Vector<double>> solver(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    solver.solve(system_matrix, newton_update, system_rhs, preconditioner); 

    hanging_node_constraints.distribute(newton_update); 

    const double alpha = determine_step_length(); 
    current_solution.add(alpha, newton_update); 
  } 
// @sect4{MinimalSurfaceProblem::refine_mesh}  

// 这个函数的第一部分与 step-6 中的内容相同 ... 然而，在细化网格后，我们必须将旧的解决方案转移到新的解决方案中，我们在SolutionTransfer类的帮助下完成。这个过程稍微有点复杂，所以让我们详细描述一下。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::refine_mesh() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      current_solution, 
      estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.03); 

// 然后我们需要一个额外的步骤：例如，如果你标记了一个比它的邻居更精炼一次的单元，而这个邻居没有被标记为精炼，我们最终会在一个单元界面上跳过两个精炼级别。 为了避免这些情况，库将默默地也要对邻居单元进行一次细化。它通过在实际进行细化和粗化之前调用 Triangulation::prepare_coarsening_and_refinement 函数来实现。 这个函数标志着一组额外的单元格进行细化或粗化，以执行像单悬节点规则这样的规则。 调用此函数后，被标记为细化和粗化的单元格正是那些将被实际细化或粗化的单元格。通常情况下，你不需要手工操作 (Triangulation::execute_coarsening_and_refinement 为你做这个）。) 然而，我们需要初始化SolutionTransfer类，它需要知道最终将被粗化或细化的单元集，以便存储旧网格的数据并转移到新网格。因此，我们手动调用这个函数。

    triangulation.prepare_coarsening_and_refinement(); 

// 有了这个方法，我们用现在的DoFHandler初始化一个SolutionTransfer对象，并将解决方案向量附加到它上面，然后在新网格上进行实际的细化和自由度分配

    SolutionTransfer<dim> solution_transfer(dof_handler); 
    solution_transfer.prepare_for_coarsening_and_refinement(current_solution); 

    triangulation.execute_coarsening_and_refinement(); 

    dof_handler.distribute_dofs(fe); 

// 最后，我们找回插值到新网格的旧解。由于SolutionTransfer函数实际上并不存储旧的解决方案的值，而是索引，我们需要保留旧的解决方案向量，直到我们得到新的内插值。因此，我们将新的数值写入一个临时的向量中，之后才将其写入解决方案向量对象中。

    Vector<double> tmp(dof_handler.n_dofs()); 
    solution_transfer.interpolate(current_solution, tmp); 
    current_solution = tmp; 

// 在新的网格上，有不同的悬挂节点，对于这些节点，我们必须在扔掉之前的对象内容后，重新计算约束。为了安全起见，我们还应该确保当前解决方案的向量条目满足悬空节点的约束条件（参见SolutionTransfer类文档中的讨论，了解为什么必须这样做）。我们可以通过明确调用`hanging_node_constraints.distribution(current_solution)`来做到这一点；我们省略这一步，因为这将在下面调用`set_boundary_values()`的最后发生，而且没有必要做两次。

    hanging_node_constraints.clear(); 

    DoFTools::make_hanging_node_constraints(dof_handler, 
                                            hanging_node_constraints); 
    hanging_node_constraints.close(); 

// 一旦我们有了内插的解决方案和所有关于悬挂节点的信息，我们必须确保我们现在的 $u^n$ 实际上有正确的边界值。正如在介绍的最后所解释的，即使细化前的解决方案有正确的边界值，也不会自动出现这种情况，因此我们必须明确地确保它现在有。

    set_boundary_values(); 

// 我们通过更新所有剩余的数据结构来结束这个函数，向 <code>setup_dofs()</code> 表明这不是第一次了，它需要保留解向量的内容。

    setup_system(false); 
  } 

//  @sect4{MinimalSurfaceProblem::set_boundary_values}  

// 下一个函数确保解向量的条目尊重我们问题的边界值。 在细化了网格之后（或者刚刚开始计算），边界上可能会出现新的节点。这些节点的数值是在`refine_mesh()`中从之前的网格中简单插值出来的，而不是正确的边界值。这个问题可以通过将当前解决方案向量的所有边界节点明确设置为正确的值来解决。

// 但是有一个问题我们必须注意：如果我们有一个挂起的节点紧挨着一个新的边界节点，那么它的值也必须被调整以确保有限元场保持连续。这就是这个函数最后一行的调用所做的。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::set_boundary_values() 
  { 
    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             BoundaryValues<dim>(), 
                                             boundary_values); 
    for (auto &boundary_value : boundary_values) 
      current_solution(boundary_value.first) = boundary_value.second; 

    hanging_node_constraints.distribute(current_solution); 
  } 
// @sect4{MinimalSurfaceProblem::compute_residual}  

// 为了监测收敛性，我们需要一种方法来计算（离散）残差的规范，即在介绍中讨论的向量 $\left<F(u^n),\varphi_i\right>$ 与 $F(u)=-\nabla \cdot \left(\frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right)$ 的规范。事实证明，（尽管我们在当前版本的程序中没有使用这个功能）在确定最佳步长时需要计算残差 $\left<F(u^n+\alpha^n\;\delta u^n),\varphi_i\right>$ ，因此这就是我们在这里实现的：该函数将步长 $\alpha^n$ 作为参数。原有的功能当然是通过传递一个零作为参数得到的。

// 在下面的函数中，我们首先为残差设置一个向量，然后为评估点设置一个向量  $u^n+\alpha^n\;\delta u^n$  。接下来是我们在所有的积分操作中使用的相同的模板代码。

  template <int dim> 
  double MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const 
  { 
    Vector<double> residual(dof_handler.n_dofs()); 

    Vector<double> evaluation_point(dof_handler.n_dofs()); 
    evaluation_point = current_solution; 
    evaluation_point.add(alpha, newton_update); 

    const QGauss<dim> quadrature_formula(fe.degree + 1); 
    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_gradients | update_quadrature_points | 
                              update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double>              cell_residual(dofs_per_cell); 
    std::vector<Tensor<1, dim>> gradients(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_residual = 0; 
        fe_values.reinit(cell); 

// 实际的计算与  <code>assemble_system()</code>  中的计算差不多。我们首先评估 $u^n+\alpha^n\,\delta u^n$ 在正交点的梯度，然后计算系数 $a_n$ ，然后将其全部插入残差公式中。

        fe_values.get_function_gradients(evaluation_point, gradients); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double coeff = 
              1. / std::sqrt(1 + gradients[q] * gradients[q]); 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i 
                                   * coeff                    // * a_n 
                                   * gradients[q]             // * u_n 
                                   * fe_values.JxW(q));       // * dx 
          } 

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          residual(local_dof_indices[i]) += cell_residual(i); 
      } 

// 在这个函数的最后，我们还必须处理悬挂节点的约束和边界值的问题。关于后者，我们必须将所有对应于位于边界的自由度的条目的残差向量元素设置为零。原因是，由于那里的解的值是固定的，它们当然不是 "真正的 "自由度，因此，严格来说，我们不应该在残差向量中为它们集合条目。然而，正如我们一直所做的那样，我们想在每个单元上做完全相同的事情，因此我们并不想在上面的积分中处理某个自由度是否位于边界的问题。相反，我们将简单地在事后将这些条目设置为零。为此，我们需要确定哪些自由度实际上属于边界，然后在所有这些自由度上进行循环，并将剩余条目设置为零。这发生在以下几行中，我们已经在 step-11 中看到了使用DoFTools命名空间的适当函数。

    hanging_node_constraints.condense(residual); 

    for (types::global_dof_index i : 
         DoFTools::extract_boundary_dofs(dof_handler)) 
      residual(i) = 0; 

// 在函数的最后，我们返回残差的常数。

    return residual.l2_norm(); 
  } 

//  @sect4{MinimalSurfaceProblem::determine_step_length}  

// 正如介绍中所讨论的，如果我们总是采取全步，即计算 $u^{n+1}=u^n+\delta u^n$ ，牛顿方法经常不收敛。相反，我们需要一个阻尼参数（步长）  $\alpha^n$  并设置  $u^{n+1}=u^n+\alpha^n\delta u^n$  。这个函数是用来计算 $\alpha^n$  的。

// 在这里，我们简单地总是返回0.1。这当然是一个次优的选择：理想情况下，人们希望的是，当我们越来越接近解的时候，步长变成1，这样我们就可以享受牛顿方法的快速二次收敛。我们将在下面的结果部分讨论更好的策略， step-77 也涉及这方面的内容。

  template <int dim> 
  double MinimalSurfaceProblem<dim>::determine_step_length() const 
  { 
    return 0.1; 
  } 

//  @sect4{MinimalSurfaceProblem::output_results}  

// 从`run()`调用的最后一个函数以图形形式输出当前的解决方案（和牛顿更新），作为VTU文件。它与之前教程中使用的完全相同。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::output_results( 
    const unsigned int refinement_cycle) const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(current_solution, "solution"); 
    data_out.add_data_vector(newton_update, "update"); 
    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu"; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 
// @sect4{MinimalSurfaceProblem::run}  

// 在运行函数中，我们建立第一个网格，然后有牛顿迭代的顶层逻辑。

// 正如在介绍中所描述的，领域是围绕原点的单位圆盘，创建方式与 step-6 中所示相同。网格经过两次全局细化，然后再进行若干次适应性循环。

// 在开始牛顿循环之前，我们还需要做一些设置工作。我们需要创建基本的数据结构，并确保第一个牛顿迭代已经有了正确的边界值，这在介绍中已经讨论过了。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::run() 
  { 
    GridGenerator::hyper_ball(triangulation); 
    triangulation.refine_global(2); 

    setup_system(/*first time=*/true); 
    set_boundary_values(); 

// 接下来开始牛顿迭代。我们一直迭代到上一次迭代结束时计算的残差（规范）小于 $10^{-3}$ ，正如在 "do{ ... } while "循环结束时的检查。因为我们没有一个合理的值来初始化这个变量，所以我们只是使用可以表示为`双数'的最大值。

    double       last_residual_norm = std::numeric_limits<double>::max(); 
    unsigned int refinement_cycle   = 0; 
    do 
      { 
        std::cout << "Mesh refinement step " << refinement_cycle << std::endl; 

        if (refinement_cycle != 0) 
          refine_mesh(); 

// 在每个网格上，我们正好做五个牛顿步骤。我们在这里打印初始残差，然后在这个网格上开始迭代。

// 在每一个牛顿步骤中，首先要计算系统矩阵和右手边，然后我们存储右手边的规范作为残差，以便在决定是否停止迭代时进行检查。然后我们求解线性系统（该函数也会更新 $u^{n+1}=u^n+\alpha^n\;\delta u^n$ ），并在这个牛顿步骤结束时输出残差的准则。

// 在这个循环结束后，我们还将以图形形式输出当前网格上的解，并增加网格细化循环的计数器。

        std::cout << "  Initial residual: " << compute_residual(0) << std::endl; 

        for (unsigned int inner_iteration = 0; inner_iteration < 5; 
             ++inner_iteration) 
          { 
            assemble_system(); 
            last_residual_norm = system_rhs.l2_norm(); 

            solve(); 

            std::cout << "  Residual: " << compute_residual(0) << std::endl; 
          } 

        output_results(refinement_cycle); 

        ++refinement_cycle; 
        std::cout << std::endl; 
      } 
    while (last_residual_norm > 1e-3); 
  } 
} // namespace Step15 
// @sect4{The main function}  

// 最后是主函数。这遵循了所有其他主函数的方案。

int main() 
{ 
  try 
    { 
      using namespace Step15; 

      MinimalSurfaceProblem<2> laplace_problem_2d; 
      laplace_problem_2d.run(); 
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



