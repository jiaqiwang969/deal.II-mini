

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
 * Authors: Jean-Paul Pelteret,
 *          Wolfgang Bangerth, Colorado State University, 2021.
 * Based on step-15, authored by Sven Wetterauer, University of Heidelberg, 2012
 */



// 本教程的大部分内容是对  step-15
// 的完全复制。因此，为了简洁起见，并保持对这里所实现的变化的关注，我们将只记录新的内容，并简单地指出哪些部分的代码是对以前内容的重复。

//  @sect3{Include files}

// 本教程中包含了几个新的头文件。第一个是提供ParameterAcceptor类的声明的文件。

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// 这是第二个，这是一个包罗万象的头，它将使我们能够在这段代码中纳入自动区分（AD）功能。

#include <deal.II/differentiation/ad.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

// 而接下来的三个提供了一些使用通用 MeshWorker::mesh_loop() 框架的多线程能力。

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

// 然后，我们为这个程序打开一个命名空间，像以前的程序一样，将dealii命名空间中的所有东西导入其中。

namespace Step72
{
  using namespace dealii;
  // @sect3{The <code>MinimalSurfaceProblemParameters</code> class}

  // 在本教程中，我们将实现三种不同的方法来组装线性系统。其中一种反映了最初在
  // step-15
  // 中提供的手工实现，而另外两种则使用作为Trilinos框架的一部分提供的Sacado自动微分库。

  // 为了方便在三种实现之间进行切换，我们有这个非常基本的参数类，它只有两个可配置的选项。

  class MinimalSurfaceProblemParameters : public ParameterAcceptor
  {
  public:
    MinimalSurfaceProblemParameters();

    // 选择要使用的配方和相应的AD框架。

    // - formulation = 0 : 无辅助执行（全手工线性化）。

    // - 配方 = 1 : 有限元残差的自动线性化。

    // - formulation = 2 : 使用变量公式自动计算有限元残差和线性化。

    unsigned int formulation = 0;

    // 线性系统残差的最大可接受公差。我们将看到，一旦我们使用AD框架，装配时间就会变得很明显，所以我们将
    // step-15
    // 中选择的公差提高了一个数量级。这样，计算就不会花费太长时间来完成。

    double tolerance = 1e-2;
  };

  MinimalSurfaceProblemParameters::MinimalSurfaceProblemParameters()
    : ParameterAcceptor("Minimal Surface Problem/")
  {
    add_parameter(
      "Formulation", formulation, "", this->prm, Patterns::Integer(0, 2));
    add_parameter("Tolerance", tolerance, "", this->prm, Patterns::Double(0.0));
  }

  //  @sect3{The <code>MinimalSurfaceProblem</code> class template}

  // 该类模板与  step-15  中的内容基本相同。该类的唯一功能变化是：。

  // -
  // run()函数现在接收两个参数：一个是选择采用哪种装配方式，一个是允许的最终残差的公差，以及

  // -
  // 现在有三个不同的装配函数来实现线性系统的三种装配方法。我们将在后面提供关于这些的细节。

  template <int dim>
  class MinimalSurfaceProblem
  {
  public:
    MinimalSurfaceProblem();

    void
    run(const int formulation, const double tolerance);

  private:
    void
    setup_system(const bool initial_step);
    void
    assemble_system_unassisted();
    void
    assemble_system_with_residual_linearization();
    void
    assemble_system_using_energy_functional();
    void
    solve();
    void
    refine_mesh();
    void
    set_boundary_values();
    double
    compute_residual(const double alpha) const;
    double
    determine_step_length() const;
    void
    output_results(const unsigned int refinement_cycle) const;

    Triangulation<dim> triangulation;

    DoFHandler<dim> dof_handler;
    FE_Q<dim>       fe;
    QGauss<dim>     quadrature_formula;

    AffineConstraints<double> hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> current_solution;
    Vector<double> newton_update;
    Vector<double> system_rhs;
  };
  // @sect3{Boundary condition}

  //应用于该问题的边界条件没有变化。

  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };

  template <int dim>
  double
  BoundaryValues<dim>::value(const Point<dim> &p,
                             const unsigned int /*component*/) const
  {
    return std::sin(2 * numbers::PI * (p[0] + p[1]));
  }
  // @sect3{The <code>MinimalSurfaceProblem</code> class implementation}
  // @sect4{MinimalSurfaceProblem::MinimalSurfaceProblem}

  // 对类的构造函数没有做任何修改。

  template <int dim>
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
    : dof_handler(triangulation)
    , fe(2)
    , quadrature_formula(fe.degree + 1)
  {}
  // @sect4{MinimalSurfaceProblem::setup_system}

  // 设置类数据结构的函数没有任何变化，即DoFHandler、应用于问题的悬挂节点约束以及线性系统。

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
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

    newton_update.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    hanging_node_constraints.condense(dsp);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }
  // @sect4{Assembling the linear system}
  // @sect5{Manual assembly}

  // 汇编函数是本教程的有趣贡献。assemble_system_unassisted()方法实现了与
  // step-15 中详述的完全相同的装配函数，但在这个例子中，我们使用
  // MeshWorker::mesh_loop()
  // 函数来多线程装配过程。这样做的原因很简单。当使用自动分化时，我们知道会有一些额外的计算开销产生。为了减轻这种性能损失，我们希望尽可能多地利用（容易获得的）计算资源。
  // MeshWorker::mesh_loop()
  // 的概念使这成为一个相对简单的任务。同时，为了公平比较，我们需要对在计算残差或其线性化时不使用任何援助的实现做同样的事情。(
  // MeshWorker::mesh_loop() 函数首先在 step-12 和 step-16
  // 中讨论，如果你想阅读它的话。)

  // 实现多线程所需的步骤在这三个函数中是相同的，所以我们将利用assemble_system_unassisted()函数的机会，重点讨论多线程本身。

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::assemble_system_unassisted()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    //  MeshWorker::mesh_loop()
    //  希望我们提供两个示范性的数据结构。第一个，`ScratchData`，是用来存储所有要在线程间重复使用的大数据。`CopyData`将保存来自每个单元的对线性系统的贡献。这些独立的矩阵-向量对必须按顺序累积到全局线性系统中。由于我们不需要
    //  MeshWorker::ScratchData 和 MeshWorker::CopyData
    //  类已经提供的东西，所以我们使用这些确切的类定义来解决我们的问题。请注意，我们只需要一个局部矩阵、局部右手向量和单元自由度索引向量的单个实例--因此
    //  MeshWorker::CopyData 的三个模板参数都是`1`。

    using ScratchData = MeshWorker::ScratchData<dim>;
    using CopyData    = MeshWorker::CopyData<1, 1, 1>;

    // 我们还需要知道我们在装配过程中要处理的迭代器的类型。为了简单起见，我们只要求编译器使用decltype()指定器为我们解决这个问题，知道我们将在由  @p dof_handler.  拥有的活动单元上迭代。
    using CellIteratorType = decltype(dof_handler.begin_active());

    // 在这里我们初始化示例的数据结构。因为我们知道我们需要计算形状函数梯度、加权雅各布和四分位点在实空间的位置，所以我们把这些标志传给类的构造函数。

    const ScratchData sample_scratch_data(fe,
                                          quadrature_formula,
                                          update_gradients |
                                            update_quadrature_points |
                                            update_JxW_values);
    const CopyData    sample_copy_data(dofs_per_cell);

    // 现在我们定义一个lambda函数，它将在一个单元格上执行装配。三个参数是由于我们将传递给该最终调用的参数，将被 MeshWorker::mesh_loop(), 所期望的参数。我们还捕获了 @p this 指针，这意味着我们将可以访问 "this"（即当前的`MinimalSurfaceProblem<dim>`）类实例，以及它的私有成员数据（因为lambda函数被定义在MinimalSurfaceProblem<dim>方法中）。

    // 在函数的顶部，我们初始化了依赖于正在执行工作的单元的数据结构。请注意，重新初始化的调用实际上返回了一个FEValues对象的实例，该对象被初始化并存储在`scratch_data`对象中（因此，被重复使用）。

    // 同样地，我们从 MeshWorker::mesh_loop()
    // 提供的`copy_data`实例中获得本地矩阵、本地RHS向量和本地单元格DoF指数的别名。然后我们初始化单元格的DoF指数，因为我们知道本地矩阵和向量的大小已经正确。

    const auto cell_worker = [this](const CellIteratorType &cell,
                                    ScratchData &           scratch_data,
                                    CopyData &              copy_data) {
      const auto &fe_values = scratch_data.reinit(cell);

      FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
      Vector<double> &                      cell_rhs    = copy_data.vectors[0];
      std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];
      cell->get_dof_indices(local_dof_indices);

      // 对于牛顿方法，我们需要问题被线性化的那一点的解的梯度。

      // 一旦我们有了这个梯度，我们就可以用通常的方法对这个单元进行装配。 与
      // step-15
      // 的一个小区别是，我们使用了（相当方便的）基于范围的循环来迭代所有的正交点和自由度。

      std::vector<Tensor<1, dim>> old_solution_gradients(
        fe_values.n_quadrature_points);
      fe_values.get_function_gradients(current_solution,
                                       old_solution_gradients);

      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const double coeff =
            1.0 / std::sqrt(1.0 + old_solution_gradients[q] *
                                    old_solution_gradients[q]);

          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
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
    };

    //  MeshWorker::mesh_loop()
    //  要求的第二个lambda函数是一个执行累积全局线性系统中的局部贡献的任务。这正是这个函数所做的，实现的细节在前面已经看到过。需要认识的主要一点是，局部贡献被存储在传入该函数的`copy_data`实例中。这个`copy_data`在
    //  @a 对`cell_worker`的一些调用中已经被填满了数据。

    const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
      const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
      const Vector<double> &    cell_rhs    = copy_data.vectors[0];
      const std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    };

    // 我们已经有了所有需要的函数定义，所以现在我们调用 MeshWorker::mesh_loop()
    // 来执行实际的装配。
    // 我们传递一个标志作为最后的参数，说明我们只想对单元格进行装配。在内部，
    // MeshWorker::mesh_loop()
    // 然后将可用的工作分配给不同的线程，有效地利用当今几乎所有的处理器所提供的多核。

    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          cell_worker,
                          copier,
                          sample_scratch_data,
                          sample_copy_data,
                          MeshWorker::assemble_own_cells);

    // 最后，正如在  step-15
    // 中所做的那样，我们从系统中移除悬空的节点，并对定义牛顿更新的线性系统应用零边界值
    // $\delta u^n$  。

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
  // @sect5{Assembly via differentiation of the residual vector}

  // 正如介绍中所述，我们需要为第二种方法做的是实现 $F(U)^K$
  // 单元对残差向量的局部贡献，然后让AD机器处理如何计算它的导数
  // $J(U)_{ij}^K=\frac{\partial F(U)^K_i}{\partial U_j}$ 。

  // 对于下面的内容，请记住，
  // @f[
  //    F(U)_i^K \dealcoloneq
  //    \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  //    u|^{2}}} \nabla u \right] \, dV ,
  //  @f]
  //  其中 $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$  。

  // 我们来看看这在实践中是如何实现的。

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::assemble_system_with_residual_linearization()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    using ScratchData      = MeshWorker::ScratchData<dim>;
    using CopyData         = MeshWorker::CopyData<1, 1, 1>;
    using CellIteratorType = decltype(dof_handler.begin_active());

    const ScratchData sample_scratch_data(fe,
                                          quadrature_formula,
                                          update_gradients |
                                            update_quadrature_points |
                                            update_JxW_values);
    const CopyData    sample_copy_data(dofs_per_cell);

    // 我们将利用  step-71
    // 中所示的技术，预先定义我们要使用的AD数据结构。在这种情况下，我们选择辅助类，它将使用Sacado向前自动微分类型自动计算有限元残差的线性化。这些数字类型可以只用来计算一阶导数。这正是我们想要的，因为我们知道我们将只对残差进行线性化，这意味着我们只需要计算一阶导数。计算的返回值将是`double`类型。

    // 我们还需要一个提取器来检索一些与问题的现场解决方案有关的数据。

    using ADHelper = Differentiation::AD::ResidualLinearization<
      Differentiation::AD::NumberTypes::sacado_dfad,
      double>;
    using ADNumberType = typename ADHelper::ad_type;

    const FEValuesExtractors::Scalar u_fe(0);

    // 有了这个，让我们定义lambda函数，它将被用来计算单元格对雅各布矩阵和右手边的贡献。

    const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
                                           ScratchData &           scratch_data,
                                           CopyData &              copy_data) {
      const auto &       fe_values     = scratch_data.reinit(cell);
      const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell();

      FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
      Vector<double> &                      cell_rhs    = copy_data.vectors[0];
      std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];
      cell->get_dof_indices(local_dof_indices);

      // 我们现在要创建并初始化一个AD辅助类的实例。要做到这一点，我们需要指定有多少个自变量和因变量。自变量将是我们的解向量所具有的局部自由度的数量，即离散化解向量
      // $u (\mathbf{x})|_K = \sum\limits_{j} U^K_i \varphi_j(\mathbf{x})$
      // 的每元素表示中的数字 $j$
      // ，它表示每个有限元素有多少个解系数。在deal.II中，这等于
      // FiniteElement::dofs_per_cell.
      // ，自变量的数量将是我们要形成的局部残差向量的条目数。在这个特定的问题中（就像许多其他采用[标准Galerkin方法](https:en.wikipedia.org/wiki/Galerkin_method)的问题一样），局部求解系数的数量与局部残差方程的数量相符。

      const unsigned int n_independent_variables = local_dof_indices.size();
      const unsigned int n_dependent_variables   = dofs_per_cell;
      ADHelper ad_helper(n_independent_variables, n_dependent_variables);

      // 接下来，我们将解决方案的值告知帮助器，即我们希望线性化的 $U_j$
      // 的实际值。由于这是在每个元素上单独进行的，我们必须从全局解决方案向量中提取解决方案的系数。换句话说，我们将所有这些系数
      // $U_j$ （其中 $j$ 是一个局部自由度）定义为进入向量 $F(U)^{K}$
      // （因果函数）计算的自变量。
      //然后，
      //我们就得到了由可自动微分的数字表示的自由度值的完整集合。对这些变量进行的操作从这一点开始被AD库跟踪，直到对象超出范围。所以正是这些变量
      //<em>  </em> ，我们将对其计算残差项的导数。

      ad_helper.register_dof_values(current_solution, local_dof_indices);

      const std::vector<ADNumberType> &dof_values_ad =
        ad_helper.get_sensitive_dof_values();

      // 然后我们做一些特定问题的任务，首先是根据 "敏感 "的AD自由度值计算所有数值、（空间）梯度等。在这个例子中，我们要检索每个正交点的解梯度。请注意，现在解梯度对自由度值很敏感，因为它们使用 @p ADNumberType 作为标量类型， @p dof_values_ad 矢量提供局部自由度值。

      std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
        fe_values.n_quadrature_points);
      fe_values[u_fe].get_function_gradients_from_local_dof_values(
        dof_values_ad, old_solution_gradients);

      // 我们声明的下一个变量将存储单元格残余向量贡献。这是相当不言自明的，除了一个<b>very important</b>的细节。请注意，向量中的每个条目都是手工初始化的，数值为0。这是一个 <em> 强烈推荐的 </em> 做法，因为一些AD库似乎没有安全地初始化这些数字类型的内部数据结构。不这样做可能会导致一些非常难以理解或检测的错误（感谢这个程序的作者出于一般的坏经验而提到这一点）。因此，出于谨慎考虑，值得明确地将初始值归零。在这之后，除了符号的改变，残差集看起来和我们之前看到的单元格RHS向量差不多。我们在所有正交点上循环，确保系数现在通过使用正确的`ADNumberType'来编码它对（敏感的）有限元DoF值的依赖性，最后我们组装残差向量的组件。为了完全清楚，有限元形状函数（及其梯度等）以及 "JxW "值仍然是标量值，但每个正交点的 @p coeff 和 @p old_solution_gradients 是以独立变量计算的。

      std::vector<ADNumberType> residual_ad(n_dependent_variables,
                                            ADNumberType(0.0));
      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const ADNumberType coeff =
            1.0 / std::sqrt(1.0 + old_solution_gradients[q] *
                                    old_solution_gradients[q]);

          for (const unsigned int i : fe_values.dof_indices())
            {
              residual_ad[i] += (fe_values.shape_grad(i, q)   // \nabla \phi_i
                                 * coeff                      // * a_n
                                 * old_solution_gradients[q]) // * u_n
                                * fe_values.JxW(q);           // * dx
            }
        }

      // 一旦我们计算出完整的单元格残差向量，我们就可以将其注册到辅助类。

      // 此后，我们在评估点计算残差值（基本上是从我们已经计算出来的东西中提取出真实的值）和它们的Jacobian（每个残差分量相对于所有单元DoF的线性化）。为了组装成全局线性系统，我们必须尊重残差和RHS贡献之间的符号差异。对于牛顿方法，右手边的向量需要等于*负的残差向量。

      ad_helper.register_residual_vector(residual_ad);

      ad_helper.compute_residual(cell_rhs);
      cell_rhs *= -1.0;

      ad_helper.compute_linearization(cell_matrix);
    };

    // 该函数的剩余部分等于我们之前的内容。

    const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
      const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
      const Vector<double> &    cell_rhs    = copy_data.vectors[0];
      const std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    };

    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          cell_worker,
                          copier,
                          sample_scratch_data,
                          sample_copy_data,
                          MeshWorker::assemble_own_cells);

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
  // @sect5{Assembly via differentiation of the energy functional}

  // 在这第三种方法中，我们将残差和雅各布作为局部能量函数
  // @f[
  //     E\left( U \right)^K
  //      \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
  //      \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
  //      \mathbf{X}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times
  //      W_{q}}_{\text{JxW(q)}}
  //  @f]
  //  的第一和第二导数来计算，能量密度由
  //  @f[
  //    \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}} .
  //  @f]给出。

  // 我们再来看看这是如何做到的。

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::assemble_system_using_energy_functional()
  {
    system_matrix = 0;
    system_rhs    = 0;

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    using ScratchData      = MeshWorker::ScratchData<dim>;
    using CopyData         = MeshWorker::CopyData<1, 1, 1>;
    using CellIteratorType = decltype(dof_handler.begin_active());

    const ScratchData sample_scratch_data(fe,
                                          quadrature_formula,
                                          update_gradients |
                                            update_quadrature_points |
                                            update_JxW_values);
    const CopyData    sample_copy_data(dofs_per_cell);

    // 在这个装配过程的实现中，我们选择了辅助类，它将使用嵌套的Sacado前向自动微分类型自动计算残差及其从单元贡献到能量函数的线性化。所选的数字类型可以用来计算第一和第二导数。我们需要这样做，因为残差定义为势能对DoF值的敏感性（即其梯度）。然后我们需要将残差线性化，这意味着必须计算势能的二阶导数。你可能想把这与之前函数中使用的
    // "ADHelper "的定义进行比较，在那里我们使用
    // `Differentiation::AD::ResidualLinearization<Differentiation::AD::NumberTypes::sacado_dfad,double>`.
    // 。
    using ADHelper = Differentiation::AD::EnergyFunctional<
      Differentiation::AD::NumberTypes::sacado_dfad_dfad,
      double>;
    using ADNumberType = typename ADHelper::ad_type;

    const FEValuesExtractors::Scalar u_fe(0);

    // 然后让我们再次定义lambda函数，对一个单元进行积分。

    // 为了初始化辅助类的实例，我们现在只需要预先知道自变量的数量（即与元素解向量相关的自由度数量）。这是因为由能量函数产生的二阶导数矩阵必然是平方的（顺便说一下，也是对称的）。

    const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
                                           ScratchData &           scratch_data,
                                           CopyData &              copy_data) {
      const auto &fe_values = scratch_data.reinit(cell);

      FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
      Vector<double> &                      cell_rhs    = copy_data.vectors[0];
      std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];
      cell->get_dof_indices(local_dof_indices);

      const unsigned int n_independent_variables = local_dof_indices.size();
      ADHelper           ad_helper(n_independent_variables);

      // 再一次，我们将所有的单元格DoFs值注册到帮助器中，然后提取这些值的 "敏感
      // "变体，用于后续必须区分的操作--其中之一是计算解决方案的梯度。

      ad_helper.register_dof_values(current_solution, local_dof_indices);

      const std::vector<ADNumberType> &dof_values_ad =
        ad_helper.get_sensitive_dof_values();

      std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
        fe_values.n_quadrature_points);
      fe_values[u_fe].get_function_gradients_from_local_dof_values(
        dof_values_ad, old_solution_gradients);

      // 我们接下来创建一个变量来存储电池的总能量。我们再一次强调，我们明确地对这个值进行零初始化，从而确保这个起始值的数据的完整性。

      // 我们的目的是计算细胞总能量，它是内部能量（由于右手函数，通常是 $U$
      // 的线性）和外部能量的总和。在这种特殊情况下，我们没有外部能量（例如，来自源项或诺伊曼边界条件），所以我们将关注内部能量部分。

      // 事实上，计算 $E(U)^K$ 几乎是微不足道的，只需要以下几行。

      ADNumberType energy_ad = ADNumberType(0.0);
      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const ADNumberType psi = std::sqrt(1.0 + old_solution_gradients[q] *
                                                     old_solution_gradients[q]);

          energy_ad += psi * fe_values.JxW(q);
        }

      // 在我们计算出这个单元的总能量后，我们将把它注册到帮助器上。
      // 在此基础上，我们现在可以计算出所需的数量，即残差值和它们在评估点的雅各布系数。和以前一样，牛顿的右手边需要是残差的负数。

      ad_helper.register_energy_functional(energy_ad);

      ad_helper.compute_residual(cell_rhs);
      cell_rhs *= -1.0;
    };

    // 与前两个函数一样，函数的剩余部分与之前一样。

    const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
      const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
      const Vector<double> &    cell_rhs    = copy_data.vectors[0];
      const std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices[0];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    };

    MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                          cell_worker,
                          copier,
                          sample_scratch_data,
                          sample_copy_data,
                          MeshWorker::assemble_own_cells);

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
  // @sect4{MinimalSurfaceProblem::solve}

  // 解算函数与  step-15  中使用的相同。

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::solve()
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

  //自 step-15
  //以来，在网格细化程序和适应性网格之间的解决方案的转移方面没有任何变化。

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::refine_mesh()
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

    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim> solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(current_solution);
    triangulation.execute_coarsening_and_refinement();

    dof_handler.distribute_dofs(fe);

    Vector<double> tmp(dof_handler.n_dofs());
    solution_transfer.interpolate(current_solution, tmp);
    current_solution = tmp;

    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    set_boundary_values();
  }

  //  @sect4{MinimalSurfaceProblem::set_boundary_values}

  // 边界条件的选择仍然与 step-15 相同 ...

  template <int dim>
  void
  MinimalSurfaceProblem<dim>::set_boundary_values()
  {
    std::map<types::global_dof_index, double> boundary_values;
  };
  template <int dim> 
                                             BoundaryValues<dim>(), 
                                             boundary_values);
  for (auto &boundary_value : boundary_values)
    current_solution(boundary_value.first) = boundary_value.second;

  hanging_node_constraints.distribute(current_solution);
}
// @sect4{MinimalSurfaceProblem::compute_residual}

// ...就像在求解迭代过程中用来计算残差的函数一样。如果真的需要，我们可以用能量函数的微分来代替它，但是为了简单起见，我们在这里只是简单地复制我们在
// step-15 中已经有的东西。

template <int dim>
double
MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const
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

      fe_values.get_function_gradients(evaluation_point, gradients);

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          const double coeff =
            1.0 / std::sqrt(1.0 + gradients[q] * gradients[q]);

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

  hanging_node_constraints.condense(residual);

  for (types::global_dof_index i : DoFTools::extract_boundary_dofs(dof_handler))
    residual(i) = 0;

  return residual.l2_norm();
}

//  @sect4{MinimalSurfaceProblem::determine_step_length}

// 非线性迭代程序的步长（或欠松系数）的选择仍然固定在  step-15
// 中选择和讨论的值。

template <int dim>
double
MinimalSurfaceProblem<dim>::determine_step_length() const
{
  return 0.1;
}

//  @sect4{MinimalSurfaceProblem::output_results}

// 从`run()`调用的最后一个函数以图形形式输出当前的解决方案（和牛顿更新），作为VTU文件。它与之前教程中使用的完全相同。

template <int dim>
void
MinimalSurfaceProblem<dim>::output_results(
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

// 在运行函数中，大部分内容与最初在  step-15
// 中实现的相同。唯一可以观察到的变化是，我们现在可以（通过参数文件）选择系统残差的最终可接受的公差是什么，并且我们可以选择我们希望利用的装配方法。为了使第二个选择明确，我们向控制台输出一些信息，表明选择。由于我们对比较三种方法中每一种的装配时间感兴趣，我们还添加了一个计时器，跟踪装配过程中所花费的时间。我们还跟踪了解决线性系统所需的时间，这样我们就可以将这些数字与通常需要最长时间执行的那部分代码进行对比。

template <int dim>
void
MinimalSurfaceProblem<dim>::run(const int formulation, const double tolerance)
{
  std::cout << "******** Assembly approach ********" << std::endl;
  const std::array<std::string, 3> method_descriptions = {
    {"Unassisted implementation (full hand linearization).",
     "Automated linearization of the finite element residual.",
     "Automated computation of finite element residual and linearization using a variational formulation."}};
  AssertIndexRange(formulation, method_descriptions.size());
  std::cout << method_descriptions[formulation] << std::endl << std::endl;

  TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);

  GridGenerator::hyper_ball(triangulation);
  triangulation.refine_global(2);

  setup_system(/*first time=*/true);
  set_boundary_values();

  double       last_residual_norm = std::numeric_limits<double>::max();
  unsigned int refinement_cycle   = 0;
  do
    {
      std::cout << "Mesh refinement step " << refinement_cycle << std::endl;

      if (refinement_cycle != 0)
        refine_mesh();

      std::cout << "  Initial residual: " << compute_residual(0) << std::endl;

      for (unsigned int inner_iteration = 0; inner_iteration < 5;
           ++inner_iteration)
        {
          {
            TimerOutput::Scope t(timer, "Assemble");

            if (formulation == 0)
              assemble_system_unassisted();
            else if (formulation == 1)
              assemble_system_with_residual_linearization();
            else if (formulation == 2)
              assemble_system_using_energy_functional();
            else
              AssertThrow(false, ExcNotImplemented());
          }

          last_residual_norm = system_rhs.l2_norm();

          {
            TimerOutput::Scope t(timer, "Solve");
            solve();
          }

          std::cout << "  Residual: " << compute_residual(0) << std::endl;
        }

      output_results(refinement_cycle);

      ++refinement_cycle;
      std::cout << std::endl;
  } while (last_residual_norm > tolerance);
}
} // namespace Step72
// @sect4{The main function}

// 最后是主函数。它遵循大多数其他主函数的方案，但有两个明显的例外。

// - 我们调用 Utilities::MPI::MPI_InitFinalize
// ，以便（通过一个隐藏的默认参数）设置使用多线程任务执行的线程数。

// -
// 我们还有几行专门用于读取或初始化用户定义的参数，这些参数将在程序执行过程中被考虑。

int
main(int argc, char *argv[])
{
  try
    {
      using namespace Step72;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

      std::string prm_file;
      if (argc > 1)
        prm_file = argv[1];
      else
        prm_file = "parameters.prm";

      const MinimalSurfaceProblemParameters parameters;
      ParameterAcceptor::initialize(prm_file);

      MinimalSurfaceProblem<2> minimal_surface_problem_2d;
      minimal_surface_problem_2d.run(parameters.formulation,
                                     parameters.tolerance);
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
