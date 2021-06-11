

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
 * Author: Tyler Anderson, Colorado State University, 2021
 */


// @sect3{Include files}

// 程序以通常的包含文件开始，所有这些文件你现在应该都见过了。

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_stack.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

// 然后照例将这个程序的所有内容放入一个命名空间，并将deal.II命名空间导入到我们将要工作的命名空间中。我们还定义了一个标识符，以便在
// <code>MMS</code> 被定义时可以运行MMS代码。否则，该程序就会解决原来的问题。

namespace BlackScholesSolver
{
  using namespace dealii;

#define MMS
  // @sect3{Solution Class}

  // 在使用MMS进行测试时，这部分为已知的解决方案创建一个类。这里我们使用
  // $v(\tau,S) = -\tau^2 -S^2 + 6$
  // 作为解决方案。我们需要包括求解方程和梯度，以便进行H1半规范计算。

  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution(const double maturity_time);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

  private:
    const double maturity_time;
  };

  template <int dim>
  Solution<dim>::Solution(const double maturity_time)
    : maturity_time(maturity_time)
  {
    Assert(dim == 1, ExcNotImplemented());
  }

  template <int dim>
  double
  Solution<dim>::value(const Point<dim> &p, const unsigned int component) const
  {
    return -Utilities::fixed_power<2, double>(p(component)) -
           Utilities::fixed_power<2, double>(this->get_time()) + 6;
  }

  template <int dim>
  Tensor<1, dim>
  Solution<dim>::gradient(const Point<dim> & p,
                          const unsigned int component) const
  {
    return Point<dim>(-2 * p(component));
  }

  //  @sect3{Equation Data}

  // 在下面的类和函数中，我们实现了定义这个问题的右手边和边界值，为此我们需要函数对象。右手边的选择是在介绍的最后讨论的。

  // 首先，我们处理初始条件。

  template <int dim>
  class InitialConditions : public Function<dim>
  {
  public:
    InitialConditions(const double strike_price);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    const double strike_price;
  };

  template <int dim>
  InitialConditions<dim>::InitialConditions(const double strike_price)
    : strike_price(strike_price)
  {}

  template <int dim>
  double
  InitialConditions<dim>::value(const Point<dim> & p,
                                const unsigned int component) const
  {
#ifdef MMS
    return -Utilities::fixed_power<2, double>(p(component)) + 6;
#else
    return std::max(p(component) - strike_price, 0.);
#endif
  }

  // 接下来，我们处理左边的边界条件。

  template <int dim>
  class LeftBoundaryValues : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };

  template <int dim>
  double
  LeftBoundaryValues<dim>::value(const Point<dim> &,
                                 const unsigned int /*component*/) const
  {
#ifdef MMS
    return -Utilities::fixed_power<2, double>(this->get_time()) + 6;
#else
    return 0.;
#endif
  }

  // 然后，我们处理右边的边界条件。

  template <int dim>
  class RightBoundaryValues : public Function<dim>
  {
  public:
    RightBoundaryValues(const double strike_price, const double interest_rate);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    const double strike_price;
    const double interest_rate;
  };

  template <int dim>
  RightBoundaryValues<dim>::RightBoundaryValues(const double strike_price,
                                                const double interest_rate)
    : strike_price(strike_price)
    , interest_rate(interest_rate)
  {}

  template <int dim>
  double
  RightBoundaryValues<dim>::value(const Point<dim> & p,
                                  const unsigned int component) const
  {
#ifdef MMS
    return -Utilities::fixed_power<2, double>(p(component)) -
           Utilities::fixed_power<2, double>(this->get_time()) + 6;
#else
    return (p(component) - strike_price) *
           exp((-interest_rate) * (this->get_time()));
#endif
  }

  // 最后，我们处理右边的问题。

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide(const double asset_volatility, const double interest_rate);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    const double asset_volatility;
    const double interest_rate;
  };

  template <int dim>
  RightHandSide<dim>::RightHandSide(const double asset_volatility,
                                    const double interest_rate)
    : asset_volatility(asset_volatility)
    , interest_rate(interest_rate)
  {}

  template <int dim>
  double
  RightHandSide<dim>::value(const Point<dim> & p,
                            const unsigned int component) const
  {
#ifdef MMS
    return 2 * (this->get_time()) -
           Utilities::fixed_power<2, double>(asset_volatility * p(component)) -
           2 * interest_rate * Utilities::fixed_power<2, double>(p(component)) -
           interest_rate *
             (-Utilities::fixed_power<2, double>(p(component)) -
              Utilities::fixed_power<2, double>(this->get_time()) + 6);
#else
    (void)p;
    (void)component;
    return 0.0;
#endif
  }

  //  @sect3{The <code>BlackScholes</code> Class}

  // 下一块是这个程序的主类的声明。这与 Step-26
  // 的教程非常相似，只是做了一些修改。必须添加新的矩阵来计算A和B矩阵，以及介绍中提到的
  // $V_{diff}$ 向量。我们还定义了问题中使用的参数。



  // -  <code>maximum_stock_price</code>
  // ：空间域的强加上限。这是允许的最大股票价格。

  // -  <code>maturity_time</code>  ：时间域的上限。这是期权到期的时间。

  // -  <code>asset_volatility</code>  ：股票价格的波动率。

  // -  <code>interest_rate</code>  : 无风险利率。

  // -  <code>strike_price</code>  ：买方在到期时可以选择购买股票的约定价格。

  // 本程序与 step-26 之间的一些细微差别是创建了 <code>a_matrix</code> and the
  // <code>b_matrix</code>
  // ，这在介绍中已经说明。然后，我们还需要存储当前时间、时间步长和当前时间步长的数字。接下来，我们将把输出存储到一个
  // <code>DataOutStack</code>
  // 的变量中，因为我们将把每个时间的解分层在上面，以创建解流形。然后，我们有一个变量来存储当前的周期和我们在计算解决方案时将运行的周期数。循环是给定一个网格的一个完整的解决方案计算。我们在每个周期之间细化一次网格，以展示我们程序的收敛特性。最后，我们将收敛数据存储到一个收敛表中。

  // 就成员函数而言，我们有一个函数可以计算每个周期的收敛信息，称为
  // <code>process_solution</code>  。这就像在  step-7  中所做的那样。

  template <int dim>
  class BlackScholes
  {
  public:
    BlackScholes();

    void
    run();

  private:
    void
    setup_system();
    void
    solve_time_step();
    void
    refine_grid();
    void
    process_solution();
    void
    add_results_for_output();
    void
    write_convergence_table();

    const double maximum_stock_price;
    const double maturity_time;
    const double asset_volatility;
    const double interest_rate;
    const double strike_price;

    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> a_matrix;
    SparseMatrix<double> b_matrix;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    double time;
    double time_step;



    const unsigned int n_time_steps;

    DataOutStack<dim>        data_out_stack;
    std::vector<std::string> solution_names;

    ConvergenceTable convergence_table;
  };
  // @sect3{The <code>BlackScholes</code> Implementation}

  // 现在，我们进入主类的实现阶段。我们将为问题中使用的各种参数设置数值。选择这些是因为它们是这些参数的相当正常的值。尽管股票价格在现实中没有上限（事实上是无限的），但我们规定了一个上限，即行权价格的两倍。两倍于行权价的选择有些武断，但它足够大，可以看到解决方案的有趣部分。

  template <int dim>
  BlackScholes<dim>::BlackScholes()
    : maximum_stock_price(1.)
    , maturity_time(1.)
    , asset_volatility(.2)
    , interest_rate(0.05)
    , strike_price(0.5)
    , fe(1)
    , dof_handler(triangulation)
    , time(0.0)
    , theta(0.5)
    , n_cycles(4)
    , n_time_steps(5000)
  {
    Assert(dim == 1, ExcNotImplemented());
  }
  // @sect4{<code>BlackScholes::setup_system</code>}

  // 下一个函数设置了DoFHandler对象，计算了约束条件，并将线性代数对象设置为正确的大小。我们还在这里通过调用库中的一个函数来计算质量矩阵。接下来我们将计算其他三个矩阵，因为这些矩阵需要
  // "手工 "计算。

  // 注意，时间步长在这里被初始化，因为计算时间步长需要成熟的时间。

  template <int dim>
  void
  BlackScholes<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    time_step = maturity_time / n_time_steps;

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,

                                    // keep_constrained_dofs =  */

                                    true）。)

      sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);
    a_matrix.reinit(sparsity_pattern);
    b_matrix.reinit(sparsity_pattern);
    system_matrix.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(dof_handler,
                                      QGauss<dim>(fe.degree + 1),
                                      mass_matrix);

    // 下面是创建非恒定系数的拉普拉斯矩阵的代码。这与介绍中的矩阵D相对应。这个非恒定系数在
    // <code>current_coefficient</code> 变量中表示。

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    QGauss<dim>        quadrature_formula(fe.degree + 1);
    FEValues<dim>      fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0.;
        fe_values.reinit(cell);
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const double current_coefficient =
              fe_values.quadrature_point(q_index).square();
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  cell_matrix(i, j) +=
                    (current_coefficient *              // (x_q)^2
                     fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                     fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                     fe_values.JxW(q_index));           // dx
              }
          }
        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              laplace_matrix.add(local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i, j));
          }
      }

    // 现在我们将创建A矩阵。下面是创建矩阵A的代码，在介绍中已经讨论过。非恒定系数再次用
    // <code>current_coefficient</code> 这个变量表示。

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0.;
        fe_values.reinit(cell);
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const Tensor<1, dim> current_coefficient =
              fe_values.quadrature_point(q_index);
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  {
                    cell_matrix(i, j) +=
                      (current_coefficient *               // x_q
                       fe_values.shape_grad(i, q_index) *  // grad phi_i(x_q)
                       fe_values.shape_value(j, q_index) * // phi_j(x_q)
                       fe_values.JxW(q_index));            // dx
                  }
              }
          }
        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              a_matrix.add(local_dof_indices[i],
                           local_dof_indices[j],
                           cell_matrix(i, j));
          }
      }

    // 最后我们将创建矩阵B。下面是创建矩阵B的代码，在介绍中已经讨论过。非恒定系数再次用
    // <code>current_coefficient</code> 这个变量表示。

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0.;
        fe_values.reinit(cell);
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const Tensor<1, dim> current_coefficient =
              fe_values.quadrature_point(q_index);
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  cell_matrix(i, j) +=
                    (current_coefficient *               // x_q
                     fe_values.shape_value(i, q_index) * // phi_i(x_q)
                     fe_values.shape_grad(j, q_index) *  // grad phi_j(x_q)
                     fe_values.JxW(q_index));            // dx
              }
          }
        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              b_matrix.add(local_dof_indices[i],
                           local_dof_indices[j],
                           cell_matrix(i, j));
          }
      }

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }
  // @sect4{<code>BlackScholes::solve_time_step</code>}

  // 下一个函数是解决单个时间步长的实际线性系统的函数。这里唯一有趣的是，我们建立的矩阵是对称正定的，所以我们可以使用共轭梯度法。

  template <int dim>
  void
  BlackScholes<dim>::solve_time_step()
  {
    SolverControl                          solver_control(1000, 1e-12);
    SolverCG<Vector<double>>               cg(solver_control);
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.0);
    cg.solve(system_matrix, solution, system_rhs, preconditioner);
    constraints.distribute(solution);
  }
  // @sect4{<code>BlackScholes::add_results_for_output</code>}

  // 这是简单地将解决方案的碎片拼接起来的功能。为此，我们在每个时间段创建一个新的层，然后添加该时间段的解决方案向量。然后，该函数使用'build_patches'将其与旧的解决方案缝合在一起。

  template <int dim>
  void
  BlackScholes<dim>::add_results_for_output()
  {
    data_out_stack.new_parameter_value(time, time_step);
    data_out_stack.attach_dof_handler(dof_handler);
    data_out_stack.add_data_vector(solution, solution_names);
    data_out_stack.build_patches(2);
    data_out_stack.finish_parameter_value();
  }
  // @sect4{<code>BlackScholes::refine_grid</code>}

  // 对于我们所做的全局细化来说，有一个函数是有些不必要的。之所以有这个函数，是为了允许以后有可能进行适应性细化。

  template <int dim>
  void
  BlackScholes<dim>::refine_grid()
  {
    triangulation.refine_global(1);
  }
  // @sect4{<code>BlackScholes::process_solution</code>}

  // 这就是我们计算收敛和误差数据的地方，以评估程序的有效性。在这里，我们计算
  // $L^2$  、 $H^1$  和 $L^{\infty}$ 的准则。

  template <int dim>
  void
  BlackScholes<dim>::process_solution()
  {
    Solution<dim> sol(maturity_time);
    sol.set_time(time);
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      sol,
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 1),
                                      VectorTools::L2_norm);
    const double L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      sol,
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 1),
                                      VectorTools::H1_seminorm);
    const double H1_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::H1_seminorm);
    const QTrapezoid<1>  q_trapezoid;
    const QIterated<dim> q_iterated(q_trapezoid, fe.degree * 2 + 1);
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      sol,
                                      difference_per_cell,
                                      q_iterated,
                                      VectorTools::Linfty_norm);
    const double Linfty_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::Linfty_norm);
    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("L2", L2_error);
    convergence_table.add_value("H1", H1_error);
    convergence_table.add_value("Linfty", Linfty_error);
  }
  // @sect4{<code>BlackScholes::write_convergence_table</code> }

  // 接下来的部分是建立收敛和误差表。通过这个，我们需要设置如何输出在
  // <code>BlackScholes::process_solution</code>
  // 期间计算的数据。首先，我们将创建标题并正确设置单元格。在这期间，我们还将规定结果的精度。然后，我们将根据
  // $L^2$  、  $H^1$  和  $L^{\infty}$
  // 规范把计算出来的误差写到控制台和错误的LaTeX文件中。

  template <int dim>
  void
  BlackScholes<dim>::write_convergence_table()
  {
    convergence_table.set_precision("L2", 3);
    convergence_table.set_precision("H1", 3);
    convergence_table.set_precision("Linfty", 3);
    convergence_table.set_scientific("L2", true);
    convergence_table.set_scientific("H1", true);
    convergence_table.set_scientific("Linfty", true);
    convergence_table.set_tex_caption("cells", "\\# cells");
    convergence_table.set_tex_caption("dofs", "\\# dofs");
    convergence_table.set_tex_caption("L2", "@f$L^2@f$-error");
    convergence_table.set_tex_caption("H1", "@f$H^1@f$-error");
    convergence_table.set_tex_caption("Linfty", "@f$L^\\infty@f$-error");
    convergence_table.set_tex_format("cells", "r");
    convergence_table.set_tex_format("dofs", "r");
    std::cout << std::endl;
    convergence_table.write_text(std::cout);
    std::string error_filename = "error";
    error_filename += "-global";
    error_filename += ".tex";
    std::ofstream error_table_file(error_filename);
    convergence_table.write_tex(error_table_file);

    // 接下来，我们将制作收敛表。我们将再次把它写到控制台和收敛LaTeX文件中。

    convergence_table.add_column_to_supercolumn("cells", "n cells");
    std::vector<std::string> new_order;
    new_order.emplace_back("n cells");
    new_order.emplace_back("H1");
    new_order.emplace_back("L2");
    convergence_table.set_column_order(new_order);
    convergence_table.evaluate_convergence_rates(
      "L2", ConvergenceTable::reduction_rate);
    convergence_table.evaluate_convergence_rates(
      "L2", ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates(
      "H1", ConvergenceTable::reduction_rate);
    convergence_table.evaluate_convergence_rates(
      "H1", ConvergenceTable::reduction_rate_log2);
    std::cout << std::endl;
    convergence_table.write_text(std::cout);
    std::string conv_filename = "convergence";
    conv_filename += "-global";
    switch (fe.degree)
      {
        case 1:
          conv_filename += "-q1";
          break;
        case 2:
          conv_filename += "-q2";
          break;
        default:
          Assert(false, ExcNotImplemented());
      }
    conv_filename += ".tex";
    std::ofstream table_file(conv_filename);
    convergence_table.write_tex(table_file);
  }
  // @sect4{<code>BlackScholes::run</code>}

  // 现在我们进入了程序的主要驱动部分。在这里我们要做的是在时间步数中循环往复，并在每次计算解向量的工作。在这里的顶部，我们设置初始细化值，然后创建一个网格。然后我们对这个网格进行一次细化。接下来，我们设置了data_out_stack对象来存储我们的解决方案。最后，我们启动一个for循环来循环处理这些循环。这让我们为每一个连续的网格细化重新计算出一个解决方案。在每次迭代开始时，我们需要重新设置时间和时间步长。我们引入一个if语句来完成这个任务，因为我们不想在第一次迭代时就这样做。

  template <int dim>
  void
  BlackScholes<dim>::run()
  {
    GridGenerator::hyper_cube(triangulation, 0.0, maximum_stock_price, true);
    triangulation.refine_global(0);

    solution_names.emplace_back("u");
    data_out_stack.declare_data_vector(solution_names,
                                       DataOutStack<dim>::dof_vector);

    Vector<double> vmult_result;
    Vector<double> forcing_terms;

    for (unsigned int cycle = 0; cycle < n_cycles; cycle++)
      {
        if (cycle != 0)
          {
            refine_grid();
            time = 0.0;
          }

        setup_system();

        std::cout << std::endl
                  << "===========================================" << std::endl
                  << "Cycle " << cycle << ':' << std::endl
                  << "Number of active cells: "
                  << triangulation.n_active_cells() << std::endl
                  << "Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl
                  << std::endl;

        VectorTools::interpolate(dof_handler,
                                 InitialConditions<dim>(strike_price),
                                 solution);

        if (cycle == (n_cycles - 1))
          {
            add_results_for_output();
          }

        // 接下来，我们运行主循环，该循环一直运行到超过成熟时间。我们首先计算方程的右侧，这在介绍中有所描述。回顾一下，它包含术语
        // $\left[-\frac{1}{4}k_n\sigma^2\mathbf{D}-k_nr\mathbf{M}+k_n\sigma^2
        // \mathbf{B}-k_nr\mathbf{A}+\mathbf{M}\right]V^{n-1}$
        // 。我们把这些项放到变量system_rhs中，借助于一个临时向量。

        vmult_result.reinit(dof_handler.n_dofs());
        forcing_terms.reinit(dof_handler.n_dofs());
        for (unsigned int timestep_number = 0; timestep_number < n_time_steps;
             ++timestep_number)
          {
            time += time_step;

            if (timestep_number % 1000 == 0)
              std::cout << "Time step " << timestep_number << " at t=" << time
                        << std::endl;

            mass_matrix.vmult(system_rhs, solution);

            laplace_matrix.vmult(vmult_result, solution);
            system_rhs.add(
              (-1) * (1 - theta) * time_step *
                Utilities::fixed_power<2, double>(asset_volatility) * 0.5,
              vmult_result);
            mass_matrix.vmult(vmult_result, solution);

            system_rhs.add((-1) * (1 - theta) * time_step * interest_rate * 2,
                           vmult_result);

            a_matrix.vmult(vmult_result, solution);
            system_rhs.add((-1) * time_step * interest_rate, vmult_result);

            b_matrix.vmult(vmult_result, solution);
            system_rhs.add(
              (-1) * Utilities::fixed_power<2, double>(asset_volatility) *
                time_step * 1,
              vmult_result);

            // 第二块是计算源项的贡献。这与术语  $-k_n\left[\frac{1}{2}F^{n-1}
            // +\frac{1}{2}F^n\right]$  相对应。下面的代码调用
            // VectorTools::create_right_hand_side  来计算向量  $F$
            // ，在这里我们在评估之前设置了右侧（源）函数的时间。这一切的结果最终都在forcing_terms变量中。

            RightHandSide<dim> rhs_function(asset_volatility, interest_rate);
            rhs_function.set_time(time);
            VectorTools::create_right_hand_side(dof_handler,
                                                QGauss<dim>(fe.degree + 1),
                                                rhs_function,
                                                forcing_terms);
            forcing_terms *= time_step * theta;
            system_rhs -= forcing_terms;

            rhs_function.set_time(time - time_step);
            VectorTools::create_right_hand_side(dof_handler,
                                                QGauss<dim>(fe.degree + 1),
                                                rhs_function,
                                                forcing_terms);
            forcing_terms *= time_step * (1 - theta);
            system_rhs -= forcing_terms;

            // 接下来，我们将强迫项添加到来自时间步长的强迫项中，同时建立矩阵
            // $\left[\mathbf{M}+
            // \frac{1}{4}k_n\sigma^2\mathbf{D}+k_nr\mathbf{M}\right]$
            // ，我们必须在每个时间步长中进行反转。这些操作的最后一块是消除线性系统中悬挂的节点约束自由度。

            system_matrix.copy_from(mass_matrix);
            system_matrix.add(
              (theta)*time_step *
                Utilities::fixed_power<2, double>(asset_volatility) * 0.5,
              laplace_matrix);
            system_matrix.add((time_step)*interest_rate * theta * (1 + 1),
                              mass_matrix);

            constraints.condense(system_matrix, system_rhs);

            // 在解决这个问题之前，我们还需要做一个操作：边界值。为此，我们创建一个边界值对象，将适当的时间设置为当前时间步长的时间，并像以前多次那样对其进行评估。其结果也被用来在线性系统中设置正确的边界值。

            {
              RightBoundaryValues<dim> right_boundary_function(strike_price,
                                                               interest_rate);
              LeftBoundaryValues<dim>  left_boundary_function;
              right_boundary_function.set_time(time);
              left_boundary_function.set_time(time);
              std::map<types::global_dof_index, double> boundary_values;
              VectorTools::interpolate_boundary_values(dof_handler,
                                                       0,
                                                       left_boundary_function,
                                                       boundary_values);
              VectorTools::interpolate_boundary_values(dof_handler,
                                                       1,
                                                       right_boundary_function,
                                                       boundary_values);
              MatrixTools::apply_boundary_values(boundary_values,
                                                 system_matrix,
                                                 solution,
                                                 system_rhs);
            }

            // 解决了这个问题，我们要做的就是求解系统，生成最后一个周期的图形数据，并创建收敛表数据。

            solve_time_step();

            if (cycle == (n_cycles - 1))
              {
                add_results_for_output();
              }
          }
#ifdef MMS
        process_solution();
#endif
      }

    const std::string filename = "solution.vtk";
    std::ofstream     output(filename);
    data_out_stack.write_vtk(output);

#ifdef MMS
    write_convergence_table();
#endif
  }

} // namespace BlackScholesSolver
// @sect3{The <code>main</code> Function}

// 走到这一步，这个程序的主函数又没有什么好讨论的了：看起来自 step-6
// 以来的所有此类函数。

int
main()
{
  try
    {
      using namespace BlackScholesSolver;

      BlackScholes<1> black_scholes_solver;
      black_scholes_solver.run();
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
