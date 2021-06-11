

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2013 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2013 
 */ 



// 程序以通常的包含文件开始，所有这些文件你现在应该都见过了。

#include <deal.II/base/utilities.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/solution_transfer.h> 
#include <deal.II/numerics/matrix_tools.h> 

#include <fstream> 
#include <iostream> 

// 然后照例将这个程序的所有内容放入一个命名空间，并将deal.II命名空间导入到我们将要工作的命名空间中。

namespace Step26 
{ 
  using namespace dealii; 
// @sect3{The <code>HeatEquation</code> class}  

// 下一个部分是这个程序的主类的声明。它沿用了以前的例子中公认的路径。如果你看过 step-6 ，例如，这里唯一值得注意的是，我们需要建立两个矩阵（质量和拉普拉斯矩阵），并保存当前和前一个时间步骤的解。然后，我们还需要存储当前时间、时间步长和当前时间步长的编号。最后一个成员变量表示介绍中讨论的theta参数，它允许我们在一个程序中处理显式和隐式欧拉方法，以及Crank-Nicolson方法和其他通用方法。

// 就成员函数而言，唯一可能的惊喜是 <code>refine_mesh</code> 函数需要最小和最大的网格细化级别的参数。这样做的目的在介绍中已经讨论过了。

  template <int dim> 
  class HeatEquation 
  { 
  public: 
    HeatEquation(); 
    void run(); 

  private: 
    void setup_system(); 
    void solve_time_step(); 
    void output_results() const; 
    void refine_mesh(const unsigned int min_grid_level, 
                     const unsigned int max_grid_level); 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> mass_matrix; 
    SparseMatrix<double> laplace_matrix; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> old_solution; 
    Vector<double> system_rhs; 

    double       time; 
    double       time_step; 
    unsigned int timestep_number; 

    const double theta; 
  }; 

//  @sect3{Equation data}  

// 在下面的类和函数中，我们实现了定义这个问题的各种数据（右手边和边界值），这些数据在这个程序中使用，我们需要函数对象。右手边的选择是在介绍的最后讨论的。对于边界值，我们选择零值，但这很容易在下面改变。

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    RightHandSide() 
      : Function<dim>() 
      , period(0.2) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

  private: 
    const double period; 
  }; 

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> & p, 
                                   const unsigned int component) const 
  { 
    (void)component; 
    AssertIndexRange(component, 1); 
    Assert(dim == 2, ExcNotImplemented()); 

    const double time = this->get_time(); 
    const double point_within_period = 
      (time / period - std::floor(time / period)); 

    if ((point_within_period >= 0.0) && (point_within_period <= 0.2)) 
      { 
        if ((p[0] > 0.5) && (p[1] > -0.5)) 
          return 1; 
        else 
          return 0; 
      } 
    else if ((point_within_period >= 0.5) && (point_within_period <= 0.7)) 
      { 
        if ((p[0] > -0.5) && (p[1] > 0.5)) 
          return 1; 
        else 
          return 0; 
      } 
    else 
      return 0; 
  } 

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> & /*p*/, 
                                    const unsigned int component) const 
  { 
    (void)component; 
    Assert(component == 0, ExcIndexRange(component, 0, 1)); 
    return 0; 
  } 

//  @sect3{The <code>HeatEquation</code> implementation}  

// 现在是实现主类的时候了。让我们从构造函数开始，它选择了一个线性元素，一个时间步长为1/500的常数（记得上面把右边的源的一个周期设置为0.2，所以我们用100个时间步长来解决每个周期），并通过设置  $\theta=1/2$  选择了Crank Nicolson方法.

  template <int dim> 
  HeatEquation<dim>::HeatEquation() 
    : fe(1) 
    , dof_handler(triangulation) 
    , time_step(1. / 500) 
    , theta(0.5) 
  {} 

//  @sect4{<code>HeatEquation::setup_system</code>}  

// 下一个函数是设置DoFHandler对象，计算约束，并将线性代数对象设置为正确的大小。我们还在这里通过简单地调用库中的两个函数来计算质量和拉普拉斯矩阵。

// 注意我们在组装矩阵时不考虑悬挂节点的约束（两个函数都有一个AffineConstraints参数，默认为一个空对象）。这是因为我们要在结合当前时间步长的矩阵后，在run()中浓缩约束。

  template <int dim> 
  void HeatEquation<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 

    std::cout << std::endl 
              << "===========================================" << std::endl 
              << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl 
              << std::endl; 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    dsp, 
                                    constraints, 
                                    /*keep_constrained_dofs =  */ true);

    sparsity_pattern.copy_from(dsp); 

    mass_matrix.reinit(sparsity_pattern); 
    laplace_matrix.reinit(sparsity_pattern); 
    system_matrix.reinit(sparsity_pattern); 

    MatrixCreator::create_mass_matrix(dof_handler, 
                                      QGauss<dim>(fe.degree + 1), 
                                      mass_matrix); 
    MatrixCreator::create_laplace_matrix(dof_handler, 
                                         QGauss<dim>(fe.degree + 1), 
                                         laplace_matrix); 

    solution.reinit(dof_handler.n_dofs()); 
    old_solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 
// @sect4{<code>HeatEquation::solve_time_step</code>}  

// 下一个函数是解决单个时间步骤的实际线性系统的函数。这里没有什么值得惊讶的。

  template <int dim> 
  void HeatEquation<dim>::solve_time_step() 
  { 
    SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.0); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    constraints.distribute(solution); 

    std::cout << "     " << solver_control.last_step() << " CG iterations." 
              << std::endl; 
  } 

//  @sect4{<code>HeatEquation::output_results</code>}  

// 在生成图形输出方面也没有什么新东西，只是我们告诉DataOut对象当前的时间和时间步长是多少，以便将其写入输出文件中。

  template <int dim> 
  void HeatEquation<dim>::output_results() const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "U"); 

    data_out.build_patches(); 

    data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number)); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtk"; 
    std::ofstream output(filename); 
    data_out.write_vtk(output); 
  } 
// @sect4{<code>HeatEquation::refine_mesh</code>}  

// 这个函数是程序中最有趣的部分。它负责自适应网格细化的工作。这个函数执行的三个任务是：首先找出需要细化/粗化的单元，然后实际进行细化，最后在两个不同的网格之间传输解向量。第一个任务是通过使用成熟的凯利误差估计器来实现的。第二项任务是实际进行再细化。这也只涉及到基本的函数，例如 <code>refine_and_coarsen_fixed_fraction</code> ，它可以细化那些具有最大估计误差的单元，这些误差加起来占60%，并粗化那些具有最小误差的单元，这些单元加起来占40%的误差。请注意，对于像当前这样的问题，即有事发生的区域正在四处移动，我们希望积极地进行粗化，以便我们能够将单元格移动到有必要的地方。

// 正如在介绍中已经讨论过的，太小的网格会导致太小的时间步长，而太大的网格会导致太小的分辨率。因此，在前两个步骤之后，我们有两个循环，将细化和粗化限制在一个允许的单元范围内。

  template <int dim> 
  void HeatEquation<dim>::refine_mesh(const unsigned int min_grid_level, 
                                      const unsigned int max_grid_level) 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, 
                                                      estimated_error_per_cell, 
                                                      0.6, 
                                                      0.4); 

    if (triangulation.n_levels() > max_grid_level) 
      for (const auto &cell : 
           triangulation.active_cell_iterators_on_level(max_grid_level)) 
        cell->clear_refine_flag(); 
    for (const auto &cell : 
         triangulation.active_cell_iterators_on_level(min_grid_level)) 
      cell->clear_coarsen_flag(); 

// 上面这两个循环略有不同，但这很容易解释。在第一个循环中，我们没有调用  <code>triangulation.end()</code>  ，而是调用  <code>triangulation.end_active(max_grid_level)</code>  。这两个调用应该产生相同的迭代器，因为迭代器是按级别排序的，不应该有任何级别高于 <code>max_grid_level</code> 的单元格。事实上，这段代码确保了这种情况的发生。

// 作为网格细化的一部分，我们需要将旧的网格中的解向量转移到新的网格中。为此，我们使用了SolutionTransfer类，我们必须准备好需要转移到新网格的解向量（一旦完成细化，我们将失去旧的网格，所以转移必须与细化同时发生）。在我们调用这个函数的时候，我们将刚刚计算出解决方案，所以我们不再需要old_solution变量（它将在网格被细化后被解决方案覆盖，也就是在时间步长结束时；见下文）。换句话说，我们只需要一个求解向量，并将其复制到一个临时对象中，当我们进一步向下调用 <code>setup_system()</code> 时，它就不会被重置。

// 因此，我们将一个SolutionTransfer对象附加到旧的DoF处理程序中，以初始化它。然后，我们准备好三角形和数据向量，以便进行细化（按照这个顺序）。

    SolutionTransfer<dim> solution_trans(dof_handler); 

    Vector<double> previous_solution; 
    previous_solution = solution; 
    triangulation.prepare_coarsening_and_refinement(); 
    solution_trans.prepare_for_coarsening_and_refinement(previous_solution); 

// 现在一切都准备好了，所以进行细化并在新网格上重新创建DoF结构，最后在 <code>setup_system</code> 函数中初始化矩阵结构和新的向量。接下来，我们实际执行从旧网格到新网格的插值解。最后一步是对解向量应用悬空节点约束，即确保位于悬空节点上的自由度值，使解是连续的。这是必要的，因为SolutionTransfer只对单元格进行局部操作，不考虑邻域。

    triangulation.execute_coarsening_and_refinement(); 
    setup_system(); 

    solution_trans.interpolate(previous_solution, solution); 
    constraints.distribute(solution); 
  } 

//  @sect4{<code>HeatEquation::run</code>}  

// 这是程序的主要驱动，我们在这里循环所有的时间步骤。在函数的顶部，我们通过重复第一个时间步长，设置初始全局网格细化的数量和自适应网格细化的初始周期数量。然后，我们创建一个网格，初始化我们要处理的各种对象，设置一个标签，说明我们在重新运行第一个时间步长时应该从哪里开始，并将初始解插值到网格上（我们在这里选择了零函数，当然，我们可以用更简单的方法，直接将解向量设置为零）。我们还输出一次初始时间步长。

//  @note  如果你是一个有经验的程序员，你可能会对我们在这段代码中使用 <code>goto</code> 语句感到吃惊   <code>goto</code> 语句现在已经不是特别受人欢迎了，因为计算机科学界的大师之一Edsgar Dijkstra在1968年写了一封信，叫做 "Go To Statement considered harmful"（见<a href="http:en.wikipedia.org/wiki/Considered_harmful">here</a>）。这段代码的作者全心全意地赞同这一观念。  <code>goto</code> 是难以理解的。事实上，deal.II几乎不包含任何出现的情况：不包括基本上是从书本上转录的代码，也不计算重复的代码片断，在写这篇笔记时，大约60万行代码中有3个位置；我们还在4个教程程序中使用它，其背景与这里完全相同。与其在这里试图证明这种情况的出现，不如先看看代码，我们在函数的最后再来讨论这个问题。

  template <int dim> 
  void HeatEquation<dim>::run() 
  { 
    const unsigned int initial_global_refinement       = 2; 
    const unsigned int n_adaptive_pre_refinement_steps = 4; 

    GridGenerator::hyper_L(triangulation); 
    triangulation.refine_global(initial_global_refinement); 

    setup_system(); 

    unsigned int pre_refinement_step = 0; 

    Vector<double> tmp; 
    Vector<double> forcing_terms; 

  start_time_iteration: 

    time            = 0.0; 
    timestep_number = 0; 

    tmp.reinit(solution.size()); 
    forcing_terms.reinit(solution.size()); 

    VectorTools::interpolate(dof_handler, 
                             Functions::ZeroFunction<dim>(), 
                             old_solution); 
    solution = old_solution; 

    output_results(); 

// 然后我们开始主循环，直到计算的时间超过我们的结束时间0.5。第一个任务是建立我们需要在每个时间步骤中解决的线性系统的右手边。回顾一下，它包含项 $MU^{n-1}-(1-\theta)k_n AU^{n-1}$  。我们把这些项放到变量system_rhs中，借助于一个临时矢量。

    while (time <= 0.5) 
      { 
        time += time_step; 
        ++timestep_number; 

        std::cout << "Time step " << timestep_number << " at t=" << time 
                  << std::endl; 

        mass_matrix.vmult(system_rhs, old_solution); 

        laplace_matrix.vmult(tmp, old_solution); 
        system_rhs.add(-(1 - theta) * time_step, tmp); 

// 第二块是计算源项的贡献。这与术语  $k_n \left[ (1-\theta)F^{n-1} + \theta F^n \right]$  相对应。下面的代码调用  VectorTools::create_right_hand_side  来计算向量  $F$  ，在这里我们在评估之前设置了右侧（源）函数的时间。这一切的结果最终都在forcing_terms变量中。

        RightHandSide<dim> rhs_function; 
        rhs_function.set_time(time); 
        VectorTools::create_right_hand_side(dof_handler, 
                                            QGauss<dim>(fe.degree + 1), 
                                            rhs_function, 
                                            tmp); 
        forcing_terms = tmp; 
        forcing_terms *= time_step * theta; 

        rhs_function.set_time(time - time_step); 
        VectorTools::create_right_hand_side(dof_handler, 
                                            QGauss<dim>(fe.degree + 1), 
                                            rhs_function, 
                                            tmp); 

        forcing_terms.add(time_step * (1 - theta), tmp); 

// 接下来，我们将强迫项加入到来自时间步长的强迫项中，同时建立矩阵 $M+k_n\theta A$ ，我们必须在每个时间步长中进行反转。这些操作的最后一块是消除线性系统中悬挂的节点约束自由度。

        system_rhs += forcing_terms; 

        system_matrix.copy_from(mass_matrix); 
        system_matrix.add(theta * time_step, laplace_matrix); 

        constraints.condense(system_matrix, system_rhs); 

// 在解决这个问题之前，我们还需要做一个操作：边界值。为此，我们创建一个边界值对象，将适当的时间设置为当前时间步长的时间，并像以前多次那样对其进行评估。其结果也被用来在线性系统中设置正确的边界值。

        { 
          BoundaryValues<dim> boundary_values_function; 
          boundary_values_function.set_time(time); 

          std::map<types::global_dof_index, double> boundary_values; 
          VectorTools::interpolate_boundary_values(dof_handler, 
                                                   0, 
                                                   boundary_values_function, 
                                                   boundary_values); 

          MatrixTools::apply_boundary_values(boundary_values, 
                                             system_matrix, 
                                             solution, 
                                             system_rhs); 
        } 

// 有了这些，我们要做的就是解决这个系统，生成图形数据，以及......

        solve_time_step(); 

        output_results(); 

// ...负责网格的细化。在这里，我们要做的是：(i)在求解过程的最开始，细化所要求的次数，之后我们跳到顶部重新开始时间迭代，(ii)之后每隔五步细化一次。

// 时间循环和程序的主要部分以开始进入下一个时间步骤结束，将old_solution设置为我们刚刚计算出的解决方案。

        if ((timestep_number == 1) && 
            (pre_refinement_step < n_adaptive_pre_refinement_steps)) 
          { 
            refine_mesh(initial_global_refinement, 
                        initial_global_refinement + 
                          n_adaptive_pre_refinement_steps); 
            ++pre_refinement_step; 

            tmp.reinit(solution.size()); 
            forcing_terms.reinit(solution.size()); 

            std::cout << std::endl; 

            goto start_time_iteration; 
          } 
        else if ((timestep_number > 0) && (timestep_number % 5 == 0)) 
          { 
            refine_mesh(initial_global_refinement, 
                        initial_global_refinement + 
                          n_adaptive_pre_refinement_steps); 
            tmp.reinit(solution.size()); 
            forcing_terms.reinit(solution.size()); 
          } 

        old_solution = solution; 
      } 
  } 
} // namespace Step26 

// 现在你已经看到了这个函数的作用，让我们再来看看  <code>goto</code>  的问题。从本质上讲，代码所做的事情是这样的。
// @code
//    void run ()
//    {
//      initialize;
//    start_time_iteration:
//      for (timestep=1...)
//      {
//         solve timestep;
//         if (timestep==1 && not happy with the result)
//         {
//           adjust some data structures;
//           goto start_time_iteration; simply try again
//         }
//         postprocess;
//      }
//    }
//  @endcode 
//  这里，"对结果满意 "的条件是我们想保留当前的网格，还是宁愿细化网格并在新网格上重新开始。我们当然可以用下面的方法来取代  <code>goto</code>  的使用。
//  @code
//    void run ()
//    {
//      initialize;
//      while (true)
//      {
//         solve timestep;
//         if (not happy with the result)
//            adjust some data structures;
//         else
//            break;
//      }
//      postprocess;


//      for (timestep=2...)
//      {
//         solve timestep;
//         postprocess;
//      }
//    }
//  @endcode 
//  这样做的好处是摆脱了 <code>goto</code> ，但缺点是必须在两个不同的地方重复实现 "解算时间步长 "和 "后处理 "操作的代码。这可以通过将这些部分的代码（在上面的实际实现中是相当大的块）放到自己的函数中来解决，但是一个带有 <code>break</code> 语句的 <code>while(true)</code> 循环并不真的比 <code>goto</code> 容易阅读或理解。

// 最后，人们可能会简单地同意，<i>in general</i> 。
// <code>goto</code>  语句是个坏主意，但要务实地指出，在某些情况下，它们可以帮助避免代码重复和尴尬的控制流。这可能就是其中之一，它与Steve McConnell在他关于良好编程实践的优秀书籍 "Code Complete"  @cite CodeComplete 中采取的立场一致（见 step-1 的介绍中提到的这本书），该书花了惊人的10页来讨论一般的 <code>goto</code> 问题。

//  @sect3{The <code>main</code> function}  

// 走到这一步，这个程序的主函数又没有什么好讨论的了：它看起来就像自 step-6 以来的所有此类函数一样。

int main() 
{ 
  try 
    { 
      using namespace Step26; 

      HeatEquation<2> heat_equation_solver; 
      heat_equation_solver.run(); 
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

