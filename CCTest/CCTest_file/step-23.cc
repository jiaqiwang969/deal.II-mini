

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2006 - 2020 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2006 
 */ 


// @sect3{Include files}  

// 我们从通常的各种各样的包含文件开始，我们在以前的许多测试中都看到过。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/data_out.h> 

#include <fstream> 
#include <iostream> 

// 这里是仅有的三个有一些新兴趣的包含文件。第一个文件已经被使用了，例如，用于 VectorTools::interpolate_boundary_values 和 MatrixTools::apply_boundary_values 函数。然而，我们在这里使用该类中的另一个函数， VectorTools::project 来计算我们的初始值，作为连续初始值的 $L^2$ 投影。此外，我们使用  VectorTools::create_right_hand_side  来生成积分  $(f^n,\phi^n_i)$  。这些以前总是由 <code>assemble_system</code> 或应用程序代码中的类似函数手工生成。然而，我们太懒了，不能在这里这么做，所以干脆使用库函数。

#include <deal.II/numerics/vector_tools.h> 

// 与此非常相似，我们也懒得写代码来组装质量矩阵和拉普拉斯矩阵，尽管这只需要从以前的任何一个教程程序中复制相关代码。相反，我们想把重点放在这个程序中真正新的东西上，因此使用了 MatrixCreator::create_mass_matrix 和 MatrixCreator::create_laplace_matrix 函数。它们被声明在这里。

#include <deal.II/numerics/matrix_tools.h> 

// 最后，这里有一个include文件，它包含了人们有时需要的各种工具函数。特别是，我们需要 Utilities::int_to_string 类，该类在给定一个整数参数后，返回它的字符串表示。它特别有用，因为它允许第二个参数，表明我们希望结果用前导零填充的数字数。我们将用它来写输出文件，其形式为 <code>solution-XXX.vtu</code> where <code>XXX</code> 表示时间步数，并且总是由三位数组成，即使我们仍然处于个位或两位数的时间步数中。

#include <deal.II/base/utilities.h> 

// 最后一步和以前所有的程序一样。

namespace Step23 
{ 
  using namespace dealii; 
// @sect3{The <code>WaveEquation</code> class}  

// 接下来是主类的声明。它的公共函数接口与其他大多数教程程序一样。值得一提的是，我们现在必须存储四个矩阵，而不是一个：质量矩阵  $M$  ，拉普拉斯矩阵  $A$  ，用于求解  $U^n$  的矩阵  $M+k^2\theta^2A$  ，以及用于求解  $V^n$  的带有边界条件的质量矩阵副本。请注意，在周围有一个额外的质量矩阵副本是有点浪费的。我们将在可能的改进部分讨论如何避免这种情况的策略。

// 同样，我们需要 $U^n,V^n$ 的解向量，以及前一个时间步骤 $U^{n-1},V^{n-1}$ 的相应向量。 <code>system_rhs</code> 将用于我们在每个时间步骤中求解两个线性系统之一时的任何右手向量。这些将在两个函数  <code>solve_u</code>  和  <code>solve_v</code>  中解决。

// 最后，变量 <code>theta</code> 用来表示参数 $\theta$ ，该参数用于定义使用哪种时间步进方案，这在介绍中已经说明。剩下的就不言而喻了。

  template <int dim> 
  class WaveEquation 
  { 
  public: 
    WaveEquation(); 
    void run(); 

  private: 
    void setup_system(); 
    void solve_u(); 
    void solve_v(); 
    void output_results() const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> mass_matrix; 
    SparseMatrix<double> laplace_matrix; 
    SparseMatrix<double> matrix_u; 
    SparseMatrix<double> matrix_v; 

    Vector<double> solution_u, solution_v; 
    Vector<double> old_solution_u, old_solution_v; 
    Vector<double> system_rhs; 

    double       time_step; 
    double       time; 
    unsigned int timestep_number; 
    const double theta; 
  }; 

//  @sect3{Equation data}  

// 在我们继续填写主类的细节之前，让我们定义与问题相对应的方程数据，即解 $u$ 及其时间导数 $v$ 的初始值和边界值，以及一个右手类。我们使用从Function类模板派生出来的类来做这件事，这个模板以前已经用过很多次了，所以下面的内容不应该是一个惊喜。

// 我们从初始值开始，对数值 $u$ 以及它的时间导数，即速度 $v$ 都选择零。

  template <int dim> 
  class InitialValuesU : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & /*p*/, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      Assert(component == 0, ExcIndexRange(component, 0, 1)); 
      return 0; 
    } 
  }; 

  template <int dim> 
  class InitialValuesV : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & /*p*/, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      Assert(component == 0, ExcIndexRange(component, 0, 1)); 
      return 0; 
    } 
  }; 

// 其次，我们有右手边的强制项。无聊的是，我们在这里也选择零。

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & /*p*/, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      Assert(component == 0, ExcIndexRange(component, 0, 1)); 
      return 0; 
    } 
  }; 

// 最后，我们有  $u$  和  $v$  的边界值。它们与介绍中描述的一样，一个是另一个的时间导数。

  template <int dim> 
  class BoundaryValuesU : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      Assert(component == 0, ExcIndexRange(component, 0, 1)); 

      if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) && 
          (p[1] > -1. / 3)) 
        return std::sin(this->get_time() * 4 * numbers::PI); 
      else 
        return 0; 
    } 
  }; 

  template <int dim> 
  class BoundaryValuesV : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      Assert(component == 0, ExcIndexRange(component, 0, 1)); 

      if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) && 
          (p[1] > -1. / 3)) 
        return (std::cos(this->get_time() * 4 * numbers::PI) * 4 * numbers::PI); 
      else 
        return 0; 
    } 
  }; 

//  @sect3{Implementation of the <code>WaveEquation</code> class}  

// 实际逻辑的实现实际上是相当短的，因为我们把组装矩阵和右手边的向量等事情交给了库。其余的实际代码不超过130行，其中相当一部分是可以从以前的例子程序中获取的模板代码（例如，解决线性系统的函数，或生成输出的函数）。

// 我们从构造函数开始（关于时间步长的选择的解释，请参见介绍中关于Courant, Friedrichs, and Lewy的部分）。

  template <int dim> 
  WaveEquation<dim>::WaveEquation() 
    : fe(1) 
    , dof_handler(triangulation) 
    , time_step(1. / 64) 
    , time(time_step) 
    , timestep_number(1) 
    , theta(0.5) 
  {} 
// @sect4{WaveEquation::setup_system}  

// 下一个函数是在程序开始时，也就是在第一个时间步骤之前，设置网格、DoFHandler以及矩阵和向量。如果你已经阅读了至少到 step-6 为止的教程程序，那么前几行是相当标准的。

  template <int dim> 
  void WaveEquation<dim>::setup_system() 
  { 
    GridGenerator::hyper_cube(triangulation, -1, 1); 
    triangulation.refine_global(7); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl; 

    dof_handler.distribute_dofs(fe); 

    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl 
              << std::endl; 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

// 然后，我们必须初始化程序过程中需要的3个矩阵：质量矩阵、拉普拉斯矩阵和在每个时间步长中求解 $M+k^2\theta^2A$ 时使用的矩阵 $U^n$ 。

// 在设置这些矩阵时，请注意它们都是利用了相同的稀疏模式对象。最后，在deal.II中矩阵和稀疏模式是独立对象的原因（与其他许多有限元或线性代数类不同）变得很清楚：在相当一部分应用中，我们必须持有几个恰好具有相同稀疏模式的矩阵，它们没有理由不共享这一信息，而不是重新建立并多次浪费内存。

// 在初始化所有这些矩阵后，我们调用库函数来建立拉普拉斯和质量矩阵。它们所需要的只是一个DoFHandler对象和一个将用于数值积分的正交公式对象。请注意，在许多方面，这些函数比我们通常在应用程序中做的要好，例如，如果一台机器有多个处理器，它们会自动并行构建矩阵：更多信息见WorkStream的文档或 @ref threads "多处理器并行计算 "模块。解决线性系统的矩阵将在run()方法中被填充，因为我们需要在每个时间步长中重新应用边界条件。

    mass_matrix.reinit(sparsity_pattern); 
    laplace_matrix.reinit(sparsity_pattern); 
    matrix_u.reinit(sparsity_pattern); 
    matrix_v.reinit(sparsity_pattern); 

    MatrixCreator::create_mass_matrix(dof_handler, 
                                      QGauss<dim>(fe.degree + 1), 
                                      mass_matrix); 
    MatrixCreator::create_laplace_matrix(dof_handler, 
                                         QGauss<dim>(fe.degree + 1), 
                                         laplace_matrix); 

// 该函数的其余部分用于将矢量大小设置为正确的值。最后一行关闭了悬挂的节点约束对象。由于我们在一个均匀细化的网格上工作，所以不存在或没有计算过约束条件（即没有必要像其他程序那样调用 DoFTools::make_hanging_node_constraints ），但无论如何，我们需要在下面的一个地方进一步设置一个约束对象。

    solution_u.reinit(dof_handler.n_dofs()); 
    solution_v.reinit(dof_handler.n_dofs()); 
    old_solution_u.reinit(dof_handler.n_dofs()); 
    old_solution_v.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.close(); 
  } 

//  @sect4{WaveEquation::solve_u and WaveEquation::solve_v}  

// 接下来的两个函数是解决与  $U^n$  和  $V^n$  的方程有关的线性系统。这两个函数并不特别有趣，因为它们基本沿用了前面所有教程程序中的方案。

// 我们可以对我们要反转的两个矩阵的预处理程序做一些小实验。然而，事实证明，对于这里的矩阵，使用雅可比或SSOR预处理器可以稍微减少解决线性系统所需的迭代次数，但由于应用预处理器的成本，在运行时间方面并不占优势。这也不是什么损失，但让我们保持简单，只做不做。

  template <int dim> 
  void WaveEquation<dim>::solve_u() 
  { 
    SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    cg.solve(matrix_u, solution_u, system_rhs, PreconditionIdentity()); 

    std::cout << "   u-equation: " << solver_control.last_step() 
              << " CG iterations." << std::endl; 
  } 

  template <int dim> 
  void WaveEquation<dim>::solve_v() 
  { 
    SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    cg.solve(matrix_v, solution_v, system_rhs, PreconditionIdentity()); 

    std::cout << "   v-equation: " << solver_control.last_step() 
              << " CG iterations." << std::endl; 
  } 

//  @sect4{WaveEquation::output_results}  

// 同样地，下面的函数也和我们之前做的差不多。唯一值得一提的是，这里我们使用 Utilities::int_to_string 函数的第二个参数，生成了一个用前导零填充的时间步长的字符串表示，长度为3个字符。

  template <int dim> 
  void WaveEquation<dim>::output_results() const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution_u, "U"); 
    data_out.add_data_vector(solution_v, "V"); 

    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 

// 像  step-15  一样，由于我们在每个时间步长写输出（而且我们要解决的系统相对简单），我们指示DataOut使用zlib压缩算法，该算法针对速度而不是磁盘使用进行了优化，因为否则绘制输出会成为一个瓶颈。

    DataOutBase::VtkFlags vtk_flags; 
    vtk_flags.compression_level = 
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
    data_out.set_flags(vtk_flags); 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 

//  @sect4{WaveEquation::run}  

// 下面是程序中唯一有趣的功能。它包含了所有时间步骤的循环，但在这之前我们必须设置网格、DoFHandler和矩阵。此外，我们必须以某种方式从初始值开始。为此，我们使用 VectorTools::project 函数，该函数接收一个描述连续函数的对象，并计算该函数在DoFHandler对象所描述的有限元空间的 $L^2$ 投影。没有比这更简单的了。

  template <int dim> 
  void WaveEquation<dim>::run() 
  { 
    setup_system(); 

    VectorTools::project(dof_handler, 
                         constraints, 
                         QGauss<dim>(fe.degree + 1), 
                         InitialValuesU<dim>(), 
                         old_solution_u); 
    VectorTools::project(dof_handler, 
                         constraints, 
                         QGauss<dim>(fe.degree + 1), 
                         InitialValuesV<dim>(), 
                         old_solution_v); 

// 接下来是循环所有的时间步骤，直到我们到达结束时间（本例中为 $T=5$ ）。在每个时间步骤中，我们首先要解决 $U^n$ ，使用方程  $(M^n + k^2\theta^2 A^n)U^n =$  。
// $(M^{n,n-1} - k^2\theta(1-\theta) A^{n,n-1})U^{n-1} + kM^{n,n-1}V^{n-1} +$  
// $k\theta \left[k \theta F^n + k(1-\theta) F^{n-1} \right]$  . 请注意，我们在所有的时间步骤中使用相同的网格，因此， $M^n=M^{n,n-1}=M$  和  $A^n=A^{n,n-1}=A$  。因此，我们首先要做的是将 $MU^{n-1} - k^2\theta(1-\theta) AU^{n-1} + kMV^{n-1}$ 和强制项相加，并将结果放入 <code>system_rhs</code> 向量中。(对于这些加法，我们需要在循环之前声明一个临时向量，以避免在每个时间步骤中重复分配内存)。

// 这里需要意识到的是我们如何将时间变量传达给描述右手边的对象：每个从函数类派生出来的对象都有一个时间字段，可以用 Function::set_time 来设置，用 Function::get_time. 来读取。 实质上，使用这种机制，所有空间和时间的函数因此被认为是在某个特定时间评估的空间的函数。这与我们在有限元程序中的典型需求非常吻合，在有限元程序中，我们几乎总是在一个时间步长上工作，而且从来没有发生过，例如，人们想在任何给定的空间位置上为所有时间评估一个时空函数。

    Vector<double> tmp(solution_u.size()); 
    Vector<double> forcing_terms(solution_u.size()); 

    for (; time <= 5; time += time_step, ++timestep_number) 
      { 
        std::cout << "Time step " << timestep_number << " at t=" << time 
                  << std::endl; 

        mass_matrix.vmult(system_rhs, old_solution_u); 

        mass_matrix.vmult(tmp, old_solution_v); 
        system_rhs.add(time_step, tmp); 

        laplace_matrix.vmult(tmp, old_solution_u); 
        system_rhs.add(-theta * (1 - theta) * time_step * time_step, tmp); 

        RightHandSide<dim> rhs_function; 
        rhs_function.set_time(time); 
        VectorTools::create_right_hand_side(dof_handler, 
                                            QGauss<dim>(fe.degree + 1), 
                                            rhs_function, 
                                            tmp); 
        forcing_terms = tmp; 
        forcing_terms *= theta * time_step; 

        rhs_function.set_time(time - time_step); 
        VectorTools::create_right_hand_side(dof_handler, 
                                            QGauss<dim>(fe.degree + 1), 
                                            rhs_function, 
                                            tmp); 

        forcing_terms.add((1 - theta) * time_step, tmp); 

        system_rhs.add(theta * time_step, forcing_terms); 

// 如此构建了第一个方程的右手向量后，我们要做的就是应用正确的边界值。至于右手边，这是一个在特定时间评估的时空函数，我们在边界节点插值，然后像通常那样用结果来应用边界值。然后将结果交给solve_u()函数。

        { 
          BoundaryValuesU<dim> boundary_values_u_function; 
          boundary_values_u_function.set_time(time); 

          std::map<types::global_dof_index, double> boundary_values; 
          VectorTools::interpolate_boundary_values(dof_handler, 
                                                   0, 
                                                   boundary_values_u_function, 
                                                   boundary_values); 

// solve_u()的矩阵在每个时间步骤中都是相同的，所以人们可以认为只在模拟开始时做一次就足够了。然而，由于我们需要对线性系统应用边界值（消除了一些矩阵的行和列，并对右手边做出了贡献），在实际应用边界数据之前，我们必须在每个时间步骤中重新填充该矩阵。实际内容非常简单：它是质量矩阵和加权拉普拉斯矩阵的总和。

          matrix_u.copy_from(mass_matrix); 
          matrix_u.add(theta * theta * time_step * time_step, laplace_matrix); 
          MatrixTools::apply_boundary_values(boundary_values, 
                                             matrix_u, 
                                             solution_u, 
                                             system_rhs); 
        } 
        solve_u(); 

// 第二步，即求解 $V^n$ ，工作原理类似，只是这次左边的矩阵是质量矩阵（我们再次复制，以便能够应用边界条件，而右边是 $MV^{n-1} - k\left[ \theta A U^n + (1-\theta) AU^{n-1}\right]$ 加上强制项。边界值的应用方式与之前相同，只是现在我们必须使用BoundaryValuesV类。

        laplace_matrix.vmult(system_rhs, solution_u); 
        system_rhs *= -theta * time_step; 

        mass_matrix.vmult(tmp, old_solution_v); 
        system_rhs += tmp; 

        laplace_matrix.vmult(tmp, old_solution_u); 
        system_rhs.add(-time_step * (1 - theta), tmp); 

        system_rhs += forcing_terms; 

        { 
          BoundaryValuesV<dim> boundary_values_v_function; 
          boundary_values_v_function.set_time(time); 

          std::map<types::global_dof_index, double> boundary_values; 
          VectorTools::interpolate_boundary_values(dof_handler, 
                                                   0, 
                                                   boundary_values_v_function, 
                                                   boundary_values); 
          matrix_v.copy_from(mass_matrix); 
          MatrixTools::apply_boundary_values(boundary_values, 
                                             matrix_v, 
                                             solution_v, 
                                             system_rhs); 
        } 
        solve_v(); 

// 最后，在计算完两个解的组成部分后，我们输出结果，计算解中的能量，并在将现在的解移入持有上一个时间步长的解的向量后，继续下一个时间步长。注意函数 SparseMatrix::matrix_norm_square 可以在一个步骤中计算 $\left<V^n,MV^n\right>$ 和 $\left<U^n,AU^n\right>$ ，为我们节省了一个临时向量和几行代码的费用。

        output_results(); 

        std::cout << "   Total energy: " 
                  << (mass_matrix.matrix_norm_square(solution_v) + 
                      laplace_matrix.matrix_norm_square(solution_u)) / 
                       2 
                  << std::endl; 

        old_solution_u = solution_u; 
        old_solution_v = solution_v; 
      } 
  } 
} // namespace Step23 
// @sect3{The <code>main</code> function}  

//剩下的就是程序的主要功能了。这里没有什么是在前面几个程序中没有展示过的。

int main() 
{ 
  try 
    { 
      using namespace Step23; 

      WaveEquation<2> wave_equation_solver; 
      wave_equation_solver.run(); 
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


