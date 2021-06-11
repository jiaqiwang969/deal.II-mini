CCTest_file/step-25.cc

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
 * Author: Ivan Christov, Wolfgang Bangerth, Texas A&M University, 2006 
 */ 


// @sect3{Include files and global variables}  

// 关于包含文件的解释，读者应该参考示例程序  step-1  到  step-4  。它们的标准顺序是  <code>base</code> -- <code>lac</code> -- <code>grid</code>  --  <code>dofs</code> -- <code>fe</code> -- <code>numerics</code>  （因为每一类大致都是建立在前面的基础上），然后是一些用于文件输入/输出和字符串流的C++头文件。

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

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <fstream> 
#include <iostream> 

// 最后一步和以前所有的程序一样。

namespace Step25 
{ 
  using namespace dealii; 
// @sect3{The <code>SineGordonProblem</code> class template}  

// 解决问题的整个算法被封装在这个类中。和以前的例子程序一样，这个类在声明时有一个模板参数，就是空间维度，这样我们就可以在一个、两个或三个空间维度上解决正弦-戈登方程。关于这个问题的独立于维度的类封装的更多信息，读者应该参考  step-3  和  step-4  。

// 与 step-23 和 step-24 相比，在程序的总体结构中没有任何值得注意的地方（当然，在各种函数的内部运作中也有！）。最明显的区别是出现了两个新的函数 <code>compute_nl_term</code> 和 <code>compute_nl_matrix</code> ，计算系统矩阵的非线性贡献和第一个方程的右手边，正如在介绍中讨论的那样。此外，我们还必须有一个向量 <code>solution_update</code> ，它包含在每个牛顿步骤中对解向量的非线性更新。

// 正如介绍中也提到的，我们在这个程序中不存储速度变量，而是质量矩阵乘以速度。这是在 <code>M_x_velocity</code> 变量中完成的（"x "是代表 "次数"）。

// 最后， <code>output_timestep_skip</code> 变量存储了在生成图形输出前每次所需的时间步数。这一点在使用精细网格（因此时间步数较小）时非常重要，在这种情况下，我们会运行大量的时间步数，并创建大量的输出文件，这些文件中的解决方案在后续文件中看起来几乎是一样的。这只会堵塞我们的可视化程序，我们应该避免创建比我们真正感兴趣的更多的输出。因此，如果这个变量被设置为大于1的值 $n$ ，那么只有在每一个 $n$ 的时间步长时才会产生输出。

  template <int dim> 
  class SineGordonProblem 
  { 
  public: 
    SineGordonProblem(); 
    void run(); 

  private: 
    void         make_grid_and_dofs(); 
    void         assemble_system(); 
    void         compute_nl_term(const Vector<double> &old_data, 
                                 const Vector<double> &new_data, 
                                 Vector<double> &      nl_term) const; 
    void         compute_nl_matrix(const Vector<double> &old_data, 
                                   const Vector<double> &new_data, 
                                   SparseMatrix<double> &nl_matrix) const; 
    unsigned int solve(); 
    void         output_results(const unsigned int timestep_number) const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 
    SparseMatrix<double> mass_matrix; 
    SparseMatrix<double> laplace_matrix; 

    const unsigned int n_global_refinements; 

    double       time; 
    const double final_time, time_step; 
    const double theta; 

    Vector<double> solution, solution_update, old_solution; 
    Vector<double> M_x_velocity; 
    Vector<double> system_rhs; 

    const unsigned int output_timestep_skip; 
  }; 
// @sect3{Initial conditions}  

// 在下面两类中，我们首先实现了本程序介绍中提到的一维、二维和三维的精确解。如果想通过比较数值解和分析解来测试程序的准确性，那么这个时空解可能会有独立的意义（但是请注意，程序使用的是有限域，而这些是无界域的分析解）。例如，这可以用 VectorTools::integrate_difference 函数来完成。再次注意（正如在 step-23 中已经讨论过的），我们如何将时空函数描述为依赖于时间变量的空间函数，该变量可以使用FunctionTime基类的 FunctionTime::set_time() 和 FunctionTime::get_time() 成员函数进行设置和查询。

  template <int dim> 
  class ExactSolution : public Function<dim> 
  { 
  public: 
    ExactSolution(const unsigned int n_components = 1, const double time = 0.) 
      : Function<dim>(n_components, time) 
    {} 

    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      const double t = this->get_time(); 

      switch (dim) 
        { 
          case 1: 
            { 
              const double m  = 0.5; 
              const double c1 = 0.; 
              const double c2 = 0.; 
              return -4. * std::atan(m / std::sqrt(1. - m * m) * 
                                     std::sin(std::sqrt(1. - m * m) * t + c2) / 
                                     std::cosh(m * p[0] + c1)); 
            } 

          case 2: 
            { 
              const double theta  = numbers::PI / 4.; 
              const double lambda = 1.; 
              const double a0     = 1.; 
              const double s      = 1.; 
              const double arg    = p[0] * std::cos(theta) + 
                                 std::sin(theta) * (p[1] * std::cosh(lambda) + 
                                                    t * std::sinh(lambda)); 
              return 4. * std::atan(a0 * std::exp(s * arg)); 
            } 

          case 3: 
            { 
              const double theta = numbers::PI / 4; 
              const double phi   = numbers::PI / 4; 
              const double tau   = 1.; 
              const double c0    = 1.; 
              const double s     = 1.; 
              const double arg   = p[0] * std::cos(theta) + 
                                 p[1] * std::sin(theta) * std::cos(phi) + 
                                 std::sin(theta) * std::sin(phi) * 
                                   (p[2] * std::cosh(tau) + t * std::sinh(tau)); 
              return 4. * std::atan(c0 * std::exp(s * arg)); 
            } 

          default: 
            Assert(false, ExcNotImplemented()); 
            return -1e8; 
        } 
    } 
  }; 

// 在本节的第二部分，我们提供初始条件。我们很懒惰（也很谨慎），不想第二次实现与上面相同的函数。相反，如果我们被查询到初始条件，我们会创建一个对象 <code>ExactSolution</code> ，将其设置为正确的时间，并让它计算当时的精确解的任何值。

  template <int dim> 
  class InitialValues : public Function<dim> 
  { 
  public: 
    InitialValues(const unsigned int n_components = 1, const double time = 0.) 
      : Function<dim>(n_components, time) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override 
    { 
      return ExactSolution<dim>(1, this->get_time()).value(p, component); 
    } 
  }; 
// @sect3{Implementation of the <code>SineGordonProblem</code> class}  

// 让我们继续讨论主类的实现，因为它实现了介绍中概述的算法。

//  @sect4{SineGordonProblem::SineGordonProblem}  

// 这是 <code>SineGordonProblem</code> 类的构造函数。它指定了所需的有限元的多项式程度，关联了一个 <code>DoFHandler</code> to the <code>triangulation</code> 对象（就像在示例程序 step-3 和 step-4 中一样），初始化了当前或初始时间、最终时间、时间步长，以及时间步长方案的 $\theta$ 值。由于我们在这里计算的解是时间周期性的，所以开始时间的实际值并不重要，我们选择它是为了让我们在一个有趣的时间开始。

// 请注意，如果我们选择显式欧拉时间步进方案（ $\theta = 0$ ），那么我们必须选择一个时间步进 $k \le h$ ，否则该方案不稳定，解中可能出现振荡。Crank-Nicolson方案（ $\theta = \frac{1}{2}$ ）和隐式Euler方案（ $\theta=1$ ）不存在这个缺陷，因为它们是无条件稳定的。然而，即使如此，时间步长也应选择在 $h$ 的数量级上，以获得一个好的解决方案。由于我们知道我们的网格是由矩形的均匀细分而来，我们可以很容易地计算出这个时间步长；如果我们有一个不同的域， step-24 中的技术使用 GridTools::minimal_cell_diameter 也是可以的。

  template <int dim> 
  SineGordonProblem<dim>::SineGordonProblem() 
    : fe(1) 
    , dof_handler(triangulation) 
    , n_global_refinements(6) 
    , time(-5.4414) 
    , final_time(2.7207) 
    , time_step(10 * 1. / std::pow(2., 1. * n_global_refinements)) 
    , theta(0.5) 
    , output_timestep_skip(1) 
  {} 
// @sect4{SineGordonProblem::make_grid_and_dofs}  

// 这个函数创建了一个 <code>dim</code> 维度的矩形网格，并对其进行了多次细化。同时，一旦自由度被集合起来， <code>SineGordonProblem</code> 类的所有矩阵和向量成员都被初始化为相应的大小。像 step-24 一样，我们使用 <code>MatrixCreator</code> 函数来生成质量矩阵 $M$ 和拉普拉斯矩阵 $A$ ，并在程序的剩余时间里将它们存储在适当的变量中。

  template <int dim> 
  void SineGordonProblem<dim>::make_grid_and_dofs() 
  { 
    GridGenerator::hyper_cube(triangulation, -10, 10); 
    triangulation.refine_global(n_global_refinements); 

    std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "   Total number of cells: " << triangulation.n_cells() 
              << std::endl; 

    dof_handler.distribute_dofs(fe); 

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
    mass_matrix.reinit(sparsity_pattern); 
    laplace_matrix.reinit(sparsity_pattern); 

    MatrixCreator::create_mass_matrix(dof_handler, 
                                      QGauss<dim>(fe.degree + 1), 
                                      mass_matrix); 
    MatrixCreator::create_laplace_matrix(dof_handler, 
                                         QGauss<dim>(fe.degree + 1), 
                                         laplace_matrix); 

    solution.reinit(dof_handler.n_dofs()); 
    solution_update.reinit(dof_handler.n_dofs()); 
    old_solution.reinit(dof_handler.n_dofs()); 
    M_x_velocity.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 
// @sect4{SineGordonProblem::assemble_system}  

// 这个函数为牛顿方法的每次迭代组装系统矩阵和右手向量。关于系统矩阵和右手边的明确公式，读者应该参考导论。

// 请注意，在每个时间步长中，我们必须把对矩阵和右手边的各种贡献加起来。与 step-23 和 step-24 相比，这需要集合更多的项，因为它们取决于前一个时间步骤或前一个非线性步骤的解。我们使用函数 <code>compute_nl_matrix</code> 和 <code>compute_nl_term</code> 来做到这一点，而本函数提供了顶层逻辑。

  template <int dim> 
  void SineGordonProblem<dim>::assemble_system() 
  { 

// 首先我们组装雅各布矩阵 $F'_h(U^{n,l})$  ，其中 $U^{n,l}$ 为方便起见被储存在向量 <code>solution</code> 中。

    system_matrix.copy_from(mass_matrix); 
    system_matrix.add(std::pow(time_step * theta, 2), laplace_matrix); 

    SparseMatrix<double> tmp_matrix(sparsity_pattern); 
    compute_nl_matrix(old_solution, solution, tmp_matrix); 
    system_matrix.add(std::pow(time_step * theta, 2), tmp_matrix); 

// 接下来我们计算右手边的向量。这只是介绍中对 $-F_h(U^{n,l})$ 的描述所暗示的矩阵-向量的组合。

    system_rhs = 0.; 

    Vector<double> tmp_vector(solution.size()); 

    mass_matrix.vmult(system_rhs, solution); 
    laplace_matrix.vmult(tmp_vector, solution); 
    system_rhs.add(std::pow(time_step * theta, 2), tmp_vector); 

    mass_matrix.vmult(tmp_vector, old_solution); 
    system_rhs.add(-1.0, tmp_vector); 
    laplace_matrix.vmult(tmp_vector, old_solution); 
    system_rhs.add(std::pow(time_step, 2) * theta * (1 - theta), tmp_vector); 

    system_rhs.add(-time_step, M_x_velocity); 

    compute_nl_term(old_solution, solution, tmp_vector); 
    system_rhs.add(std::pow(time_step, 2) * theta, tmp_vector); 

    system_rhs *= -1.; 
  } 
// @sect4{SineGordonProblem::compute_nl_term}  

// 这个函数计算向量 $S(\cdot,\cdot)$ ，它出现在分裂公式的两个方程的非线性项中。这个函数不仅简化了这个项的重复计算，而且也是我们在时间步长为隐式时使用的非线性迭代求解器的基本组成部分（即 $\theta\ne 0$  ）。此外，我们必须允许该函数接收一个 "旧 "和一个 "新 "的解决方案作为输入。这些可能不是存储在 <code>old_solution</code> and <code>solution</code> 中的问题的实际解决方案，而只是我们线性化的两个函数。为了这个函数的目的，让我们在下面这个类的文档中分别调用前两个参数  $w_{\mathrm{old}}$  和  $w_{\mathrm{new}}$  。

// 作为一个旁注，也许值得研究一下什么阶次的正交公式最适合这种类型的积分。由于 $\sin(\cdot)$ 不是一个多项式，可能没有正交公式可以准确地积分这些项。通常只需确保右手边的积分达到与离散化方案相同的精度即可，但通过选择更精确的正交公式，也许可以改善渐近收敛声明中的常数。

  template <int dim> 
  void SineGordonProblem<dim>::compute_nl_term(const Vector<double> &old_data, 
                                               const Vector<double> &new_data, 
                                               Vector<double> &nl_term) const 
  { 
    nl_term = 0; 
    const QGauss<dim> quadrature_formula(fe.degree + 1); 
    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_JxW_values | 
                              update_quadrature_points); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double>                       local_nl_term(dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
    std::vector<double>                  old_data_values(n_q_points); 
    std::vector<double>                  new_data_values(n_q_points); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        local_nl_term = 0; 

// 一旦我们将 <code>FEValues</code> 实例化重新初始化到当前单元格，我们就利用 <code>get_function_values</code> 例程来获取 "旧 "数据（大概在 $t=t_{n-1}$ ）和 "新 "数据（大概在 $t=t_n$ ）在所选正交公式节点的值。

        fe_values.reinit(cell); 
        fe_values.get_function_values(old_data, old_data_values); 
        fe_values.get_function_values(new_data, new_data_values); 

// 现在，我们可以用所需的正交公式来评估  $\int_K \sin\left[\theta w_{\mathrm{new}} + (1-\theta) w_{\mathrm{old}}\right] \,\varphi_j\,\mathrm{d}x$  。

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            local_nl_term(i) += 
              (std::sin(theta * new_data_values[q_point] + 
                        (1 - theta) * old_data_values[q_point]) * 
               fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)); 

// 我们通过将各单元的积分对全局积分的贡献相加来得出结论。

        cell->get_dof_indices(local_dof_indices); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          nl_term(local_dof_indices[i]) += local_nl_term(i); 
      } 
  } 
// @sect4{SineGordonProblem::compute_nl_matrix}  

// 这是处理非线性方案的第二个函数。它计算矩阵  $N(\cdot,\cdot)$  ，它出现在  $F(\cdot)$  的雅各布项的非线性项中。正如 <code>compute_nl_term</code> 一样，我们必须让这个函数接收一个 "旧 "和一个 "新 "的解决方案作为输入，我们再次将其分别称为 $w_{\mathrm{old}}$ 和 $w_{\mathrm{new}}$ ，如下。

  template <int dim> 
  void SineGordonProblem<dim>::compute_nl_matrix( 
    const Vector<double> &old_data, 
    const Vector<double> &new_data, 
    SparseMatrix<double> &nl_matrix) const 
  { 
    QGauss<dim>   quadrature_formula(fe.degree + 1); 
    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_JxW_values | 
                              update_quadrature_points); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> local_nl_matrix(dofs_per_cell, dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
    std::vector<double>                  old_data_values(n_q_points); 
    std::vector<double>                  new_data_values(n_q_points); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        local_nl_matrix = 0; 

// 同样，首先我们将我们的 <code>FEValues</code> 实例化重新初始化为当前单元。

        fe_values.reinit(cell); 
        fe_values.get_function_values(old_data, old_data_values); 
        fe_values.get_function_values(new_data, new_data_values); 

// 然后，我们用所需的正交公式评估 $\int_K \cos\left[\theta w_{\mathrm{new}} + (1-\theta) w_{\mathrm{old}}\right]\, \varphi_i\, \varphi_j\,\mathrm{d}x$ 。

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              local_nl_matrix(i, j) += 
                (std::cos(theta * new_data_values[q_point] + 
                          (1 - theta) * old_data_values[q_point]) * 
                 fe_values.shape_value(i, q_point) * 
                 fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 

// 最后，我们将各单元的积分对全局积分的贡献相加。

        cell->get_dof_indices(local_dof_indices); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            nl_matrix.add(local_dof_indices[i], 
                          local_dof_indices[j], 
                          local_nl_matrix(i, j)); 
      } 
  } 

//  @sect4{SineGordonProblem::solve}  

// 正如在介绍中所讨论的，这个函数在线性方程组上使用CG迭代求解器，该方程组是由牛顿方法的每个迭代的有限元空间离散化产生的，用于分割公式中的（非线性）第一个方程。该系统的解实际上是 $\delta U^{n,l}$ ，所以它被存储在 <code>solution_update</code> and used to update <code>solution</code> 的 <code>run</code> 函数中。

// 注意，我们在求解前将解的更新值重新设置为零。这是没有必要的：迭代求解器可以从任何一点开始并收敛到正确的解。如果对线性系统的解有一个很好的估计，那么从这个向量开始可能是值得的，但是作为一个一般的观察，起点并不是很重要：它必须是一个非常非常好的猜测，以减少超过几个迭代的次数。事实证明，对于这个问题，使用之前的非线性更新作为起点实际上会损害收敛性并增加所需的迭代次数，所以我们简单地将其设置为零。

// 该函数返回收敛到一个解决方案所需的迭代次数。这个数字以后将被用来在屏幕上生成输出，显示每次非线性迭代需要多少次迭代。

  template <int dim> 
  unsigned int SineGordonProblem<dim>::solve() 
  { 
    SolverControl            solver_control(1000, 1e-12 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution_update, system_rhs, preconditioner); 

    return solver_control.last_step(); 
  } 
// @sect4{SineGordonProblem::output_results}  

// 这个函数将结果输出到一个文件。它与  step-23  和  step-24  中的相应函数基本相同。

  template <int dim> 
  void SineGordonProblem<dim>::output_results( 
    const unsigned int timestep_number) const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "u"); 
    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 
    DataOutBase::VtkFlags vtk_flags; 
    vtk_flags.compression_level = 
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
    data_out.set_flags(vtk_flags); 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 
// @sect4{SineGordonProblem::run}  

// 这个函数对一切都有最高级别的控制：它运行（外部）时间步长循环，（内部）非线性求解器循环，并在每个时间步长后输出解。

  template <int dim> 
  void SineGordonProblem<dim>::run() 
  { 
    make_grid_and_dofs(); 

// 为了确认初始条件，我们必须使用函数  $u_0(x)$  来计算  $U^0$  。为此，下面我们将创建一个 <code>InitialValues</code> 类型的对象；注意，当我们创建这个对象（它来自 <code>Function</code> 类）时，我们将其内部的时间变量设置为 $t_0$ ，以表明初始条件是在 $t=t_0$ 处评估的空间和时间的函数。

// 然后我们通过使用 <code>VectorTools::project</code> 将 $u_0(x)$ 投影到网格上，产生 $U^0$ 。我们必须使用与 step-21 相同的悬挂节点约束结构： VectorTools::project 函数需要一个悬挂节点约束对象，但为了使用它，我们首先需要关闭它。

    { 
      AffineConstraints<double> constraints; 
      constraints.close(); 
      VectorTools::project(dof_handler, 
                           constraints, 
                           QGauss<dim>(fe.degree + 1), 
                           InitialValues<dim>(1, time), 
                           solution); 
    } 

// 为了完整起见，我们像其他时间步长一样，将第2个时间步长输出到一个文件。

    output_results(0); 

// 现在我们进行时间步进：在每个时间步进中，我们解决与问题的有限元离散化相对应的矩阵方程，然后根据我们在介绍中讨论的时间步进公式推进我们的解决方案。

    unsigned int timestep_number = 1; 
    for (time += time_step; time <= final_time; 
         time += time_step, ++timestep_number) 
      { 
        old_solution = solution; 

        std::cout << std::endl 
                  << "Time step #" << timestep_number << "; " 
                  << "advancing to t = " << time << "." << std::endl; 

// 在每个时间步长的开始，我们必须通过牛顿方法求解拆分公式中的非线性方程---即先求解 $\delta U^{n,l}$ ，然后再计算 $U^{n,l+1}$ ，如此反复。这种非线性迭代的停止标准是： $\|F_h(U^{n,l})\|_2 \le 10^{-6} \|F_h(U^{n,0})\|_2$  。因此，我们需要记录第一次迭代中残差的规范。

// 在每次迭代结束时，我们向控制台输出我们花了多少次线性求解器的迭代。当下面的循环完成后，我们有（一个近似的） $U^n$  。

        double initial_rhs_norm = 0.; 
        bool   first_iteration  = true; 
        do 
          { 
            assemble_system(); 

            if (first_iteration == true) 
              initial_rhs_norm = system_rhs.l2_norm(); 

            const unsigned int n_iterations = solve(); 

            solution += solution_update; 

            if (first_iteration == true) 
              std::cout << "    " << n_iterations; 
            else 
              std::cout << '+' << n_iterations; 
            first_iteration = false; 
          } 
        while (system_rhs.l2_norm() > 1e-6 * initial_rhs_norm); 

        std::cout << " CG iterations per nonlinear step." << std::endl; 

// 在得到问题的第一个方程 $t=t_n$ 的解后，我们必须更新辅助速度变量  $V^n$  。然而，我们不计算和存储 $V^n$ ，因为它不是我们在问题中直接使用的数量。因此，为了简单起见，我们直接更新 $MV^n$ 。

        Vector<double> tmp_vector(solution.size()); 
        laplace_matrix.vmult(tmp_vector, solution); 
        M_x_velocity.add(-time_step * theta, tmp_vector); 

        laplace_matrix.vmult(tmp_vector, old_solution); 
        M_x_velocity.add(-time_step * (1 - theta), tmp_vector); 

        compute_nl_term(old_solution, solution, tmp_vector); 
        M_x_velocity.add(-time_step, tmp_vector); 

// 很多时候，特别是对于细网格，我们必须选择相当小的时间步长，以使方案稳定。因此，有很多时间步长，在解的过程中 "没有什么有趣的事情发生"。为了提高整体效率--特别是加快程序速度和节省磁盘空间--我们每隔 <code>output_timestep_skip</code> 个时间步数才输出解。

        if (timestep_number % output_timestep_skip == 0) 
          output_results(timestep_number); 
      } 
  } 
} // namespace Step25 
// @sect3{The <code>main</code> function}  

// 这是该程序的主函数。它创建一个顶层类的对象并调用其主函数。如果在执行 <code>SineGordonProblem</code> 类的运行方法时抛出了异常，我们会在这里捕获并报告它们。关于异常的更多信息，读者应该参考  step-6  。

int main() 
{ 
  try 
    { 
      using namespace Step25; 

      SineGordonProblem<1> sg_problem; 
      sg_problem.run(); 
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

