

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2005 - 2021 by the deal.II authors 
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
 */ 


// @sect3{Include files}  

// 由于这个程序只是对 step-4 的改编，所以在头文件方面没有太多的新东西。在deal.II中，我们通常按照base-lac-grid-dofs-fe-numerics的顺序列出包含文件，然后是C++标准包含文件。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 

#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 

// 唯一值得关注的两个新头文件是LinearOperator和PackagedOperation类的文件。

#include <deal.II/lac/linear_operator.h> 
#include <deal.II/lac/packaged_operation.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <fstream> 
#include <iostream> 

// 这是唯一重要的新标题，即声明Raviart-Thomas有限元的标题。

#include <deal.II/fe/fe_raviart_thomas.h> 

// 最后，作为本程序中的一项奖励，我们将使用一个张量系数。由于它可能具有空间依赖性，我们认为它是一个张量值的函数。下面的include文件提供了 <code>TensorFunction</code> 类，提供了这样的功能。

#include <deal.II/base/tensor_function.h> 

// 最后一步和以前所有的程序一样。我们把所有与这个程序相关的代码放到一个命名空间中。(这个想法在  step-7  中首次提出) 。

namespace Step20 
{ 
  using namespace dealii; 
// @sect3{The <code>MixedLaplaceProblem</code> class template}  

// 同样，由于这是对 step-6 的改编，主类与该教程程序中的主类几乎相同。就成员函数而言，主要区别在于构造函数将Raviart-Thomas元素的度数作为参数（并且有一个相应的成员变量来存储这个值），并且增加了 <code>compute_error</code> 函数，在这个函数中，不出意外，我们将计算精确解和数值解之间的差异，以确定我们计算的收敛性。

  template <int dim> 
  class MixedLaplaceProblem 
  { 
  public: 
    MixedLaplaceProblem(const unsigned int degree); 
    void run(); 

  private: 
    void make_grid_and_dofs(); 
    void assemble_system(); 
    void solve(); 
    void compute_errors() const; 
    void output_results() const; 

    const unsigned int degree; 

    Triangulation<dim> triangulation; 
    FESystem<dim>      fe; 
    DoFHandler<dim>    dof_handler; 

// 第二个区别是疏散模式、系统矩阵、解和右手向量现在被封锁了。这意味着什么，人们可以用这些对象做什么，在本程序的介绍中已经解释过了，下面我们在解释这个问题的线性求解器和预处理器时也会进一步解释。

    BlockSparsityPattern      sparsity_pattern; 
    BlockSparseMatrix<double> system_matrix; 

    BlockVector<double> solution; 
    BlockVector<double> system_rhs; 
  }; 
// @sect3{Right hand side, boundary values, and exact solution}  

// 我们的下一个任务是定义我们问题的右手边（即原始拉普拉斯方程中压力的标量右手边），压力的边界值，以及一个描述压力和精确解的速度的函数，以便以后计算误差。请注意，这些函数分别有一个、一个和 <code>dim+1</code> 个分量，我们将分量的数量传递给 <code>Function@<dim@></code> 基类。对于精确解，我们只声明实际一次性返回整个解向量（即其中的所有成分）的函数。下面是各自的声明。

  namespace PrescribedSolution 
  { 
    constexpr double alpha = 0.3; 
    constexpr double beta  = 1; 

    template <int dim> 
    class RightHandSide : public Function<dim> 
    { 
    public: 
      RightHandSide() 
        : Function<dim>(1) 
      {} 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 
    }; 

    template <int dim> 
    class PressureBoundaryValues : public Function<dim> 
    { 
    public: 
      PressureBoundaryValues() 
        : Function<dim>(1) 
      {} 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 
    }; 

    template <int dim> 
    class ExactSolution : public Function<dim> 
    { 
    public: 
      ExactSolution() 
        : Function<dim>(dim + 1) 
      {} 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  value) const override; 
    }; 

// 然后我们还必须定义这些各自的函数，当然了。鉴于我们在介绍中讨论了解决方案应该是怎样的，下面的计算应该是很简单的。

    template <int dim> 
    double RightHandSide<dim>::value(const Point<dim> & /*p*/, 
                                     const unsigned int /*component*/) const 
    { 
      return 0; 
    } 

    template <int dim> 
    double 
    PressureBoundaryValues<dim>::value(const Point<dim> &p, 
                                       const unsigned int /*component*/) const 
    { 
      return -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] - 
               alpha * p[0] * p[0] * p[0] / 6); 
    } 

    template <int dim> 
    void ExactSolution<dim>::vector_value(const Point<dim> &p, 
                                          Vector<double> &  values) const 
    { 
      Assert(values.size() == dim + 1, 
             ExcDimensionMismatch(values.size(), dim + 1)); 

      values(0) = alpha * p[1] * p[1] / 2 + beta - alpha * p[0] * p[0] / 2; 
      values(1) = alpha * p[0] * p[1]; 
      values(2) = -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] - 
                    alpha * p[0] * p[0] * p[0] / 6); 
    } 

//  @sect3{The inverse permeability tensor}  

// 除了其他方程数据外，我们还想使用渗透性张量，或者更好的是--因为这是在弱形式中出现的全部内容--渗透性张量的逆，  <code>KInverse</code>  。对于验证解的精确性和确定收敛顺序的目的来说，这个张量的作用大于帮助。因此，我们将简单地把它设置为同一矩阵。

// 然而，在现实生活中的多孔介质流动模拟中，空间变化的渗透率张量是不可缺少的，我们想利用这个机会来展示使用张量值函数的技术。

// 可能不足为奇，deal.II也有一个基类，不仅适用于标量和一般的矢量值函数（ <code>Function</code> 基类），也适用于返回固定维度和等级的张量的函数， <code>TensorFunction</code> 模板。在这里，所考虑的函数返回一个dim-by-dim矩阵，即一个等级为2、维度为 <code>dim</code> 的张量。然后我们适当地选择基类的模板参数。

//  <code>TensorFunction</code> 类提供的接口本质上等同于 <code>Function</code> 类。特别是，存在一个 <code>value_list</code> 函数，它接收一个评估函数的点的列表，并在第二个参数中返回函数的值，一个张量的列表。

    template <int dim> 
    class KInverse : public TensorFunction<2, dim> 
    { 
    public: 
      KInverse() 
        : TensorFunction<2, dim>() 
      {} 

      virtual void 
      value_list(const std::vector<Point<dim>> &points, 
                 std::vector<Tensor<2, dim>> &  values) const override; 
    }; 

// 实现起来就不那么有趣了。和以前的例子一样，我们在类的开头添加一个检查，以确保输入和输出参数的大小是相同的（关于这个技术的讨论见 step-5 ）。然后我们在所有的评估点上循环，对于每一个评估点，将输出张量设置为身份矩阵。

// 在函数的顶部有一个奇怪的地方（`(void)point;`语句），值得讨论。我们放到输出`values`数组中的值实际上并不取决于函数被评估的坐标`points`数组。换句话说，`points'参数实际上是不用的，如果我们想的话，可以不给它起名字。但是我们想用`points`对象来检查`values`对象是否有正确的大小。问题是，在发布模式下，`AssertDimension`被定义为一个宏，扩展为空；然后编译器会抱怨`points`对象没有使用。消除这个警告的习惯方法是有一个评估（读取）变量的语句，但实际上不做任何事情：这就是`(void)points;`所做的：它从`points`中读取，然后将读取的结果转换为`void`，也就是什么都没有。换句话说，这句话是完全没有意义的，除了向编译器解释是的，这个变量事实上是被使用的，即使是在发布模式下。(在调试模式下，`AssertDimension`宏会扩展为从变量中读出的东西，所以在调试模式下，这个有趣的语句是没有必要的)。

    template <int dim> 
    void KInverse<dim>::value_list(const std::vector<Point<dim>> &points, 
                                   std::vector<Tensor<2, dim>> &  values) const 
    { 
      (void)points; 
      AssertDimension(points.size(), values.size()); 

      for (auto &value : values) 
        value = unit_symmetric_tensor<dim>(); 
    } 
  } // namespace PrescribedSolution 

//  @sect3{MixedLaplaceProblem class implementation}  
// @sect4{MixedLaplaceProblem::MixedLaplaceProblem}  

// 在这个类的构造函数中，我们首先存储传入的关于我们将使用的有限元的度数的值（例如，度数为0，意味着使用RT(0)和DG(0)），然后构造属于介绍中描述的空间 $X_h$ 的向量值的元素。构造函数的其余部分与早期的教程程序一样。

// 这里唯一值得描述的是，这个变量所属的 <code>fe</code> variable. The <code>FESystem</code> 类的构造函数调用有很多不同的构造函数，它们都是指将较简单的元素绑定在一起，成为一个较大的元素。在目前的情况下，我们想把一个RT(度)元素与一个DQ(度)元素结合起来。这样做的 <code>FESystem</code> 构造函数要求我们首先指定第一个基本元素（给定程度的 <code>FE_RaviartThomas</code> 对象），然后指定这个基本元素的副本数量，然后类似地指定 <code>FE_DGQ</code> 元素的种类和数量。注意Raviart-Thomas元素已经有 <code>dim</code> 个矢量分量，所以耦合元素将有 <code>dim+1</code> 个矢量分量，其中第一个 <code>dim</code> 个对应于速度变量，最后一个对应于压力。

// 我们从基本元素中构建这个元素的方式与我们在 step-8 中的方式也值得比较：在那里，我们将其构建为 <code>fe (FE_Q@<dim@>(1), dim)</code> ，即我们简单地使用 <code>dim</code> copies of the <code>FE_Q(1)</code> 元素，每个坐标方向上的位移都有一份。

  template <int dim> 
  MixedLaplaceProblem<dim>::MixedLaplaceProblem(const unsigned int degree) 
    : degree(degree) 
    , fe(FE_RaviartThomas<dim>(degree), 1, FE_DGQ<dim>(degree), 1) 
    , dof_handler(triangulation) 
  {} 

//  @sect4{MixedLaplaceProblem::make_grid_and_dofs}  

// 接下来的函数开始于众所周知的函数调用，创建和细化一个网格，然后将自由度与之关联。

  template <int dim> 
  void MixedLaplaceProblem<dim>::make_grid_and_dofs() 
  { 
    GridGenerator::hyper_cube(triangulation, -1, 1); 
    triangulation.refine_global(5); 

    dof_handler.distribute_dofs(fe); 

// 然而，接下来事情就变得不同了。正如介绍中提到的，我们要将矩阵细分为对应于速度和压力这两种不同的变量的块。为此，我们首先要确保与速度和压力相对应的指数不会混在一起。首先是所有速度自由度，然后是所有压力自由度。这样一来，全局矩阵就很好地分离成一个 $2 \times 2$ 系统。为了达到这个目的，我们必须根据自由度的矢量分量对其重新编号，这个操作已经很方便地实现了。

    DoFRenumbering::component_wise(dof_handler); 

// 接下来，我们要弄清楚这些块的大小，以便我们可以分配适当的空间量。为此，我们调用了 DoFTools::count_dofs_per_fe_component() 函数，该函数计算了某个向量分量的形状函数非零的数量。我们有 <code>dim+1</code> 个向量分量， DoFTools::count_dofs_per_fe_component() 将计算有多少个形状函数属于这些分量中的每个。

// 这里有一个问题。正如该函数的文档所描述的，它 <i>wants</i> 将  $x$  -速度形状函数的数量放入  <code>dofs_per_component[0]</code>  中，将  $y$  -速度形状函数的数量放入  <code>dofs_per_component[1]</code>  中（以及类似的3d），并将压力形状函数的数量放入  <code>dofs_per_component[dim]</code>  中 。但是，Raviart-Thomas元素的特殊性在于它是非 @ref GlossPrimitive "原始 "的，也就是说，对于Raviart-Thomas元素，所有的速度形状函数在所有分量中都是非零。换句话说，该函数不能区分 $x$ 和 $y$ 速度函数，因为<i>is</i>没有这种区分。因此，它将速度的总体数量放入 <code>dofs_per_component[c]</code>  ,  $0\le c\le \text{dim}$ 中的每一个。另一方面，压力变量的数量等于在dim-th分量中不为零的形状函数的数量。

// 利用这些知识，我们可以从 <code>dofs_per_component</code> 的第一个 <code>dim</code> 元素中的任何一个得到速度形状函数的数量，然后用下面这个来初始化向量和矩阵块的大小，以及创建输出。

//  @note  如果你觉得这个概念难以理解，你可以考虑用函数  DoFTools::count_dofs_per_fe_block()  来代替，就像我们在  step-22  的相应代码中做的那样。你可能还想阅读一下术语表中 @ref GlossBlock "块 "和 @ref GlossComponent "组件 "的区别。

    const std::vector<types::global_dof_index> dofs_per_component = 
      DoFTools::count_dofs_per_fe_component(dof_handler); 
    const unsigned int n_u = dofs_per_component[0], 
                       n_p = dofs_per_component[dim]; 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "Total number of cells: " << triangulation.n_cells() 
              << std::endl 
              << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (" << n_u << '+' << n_p << ')' << std::endl; 

// 下一个任务是为我们将要创建的矩阵分配一个稀疏模式。我们使用与前面步骤一样的压缩稀疏模式，但是由于 <code>system_matrix</code> 是一个块状矩阵，我们使用 <code>BlockDynamicSparsityPattern</code> 类，而不仅仅是 <code>DynamicSparsityPattern</code>  。这种块状稀疏模式在 $2 \times 2$ 模式下有四个块。块的大小取决于 <code>n_u</code> and <code>n_p</code> ，它持有速度和压力变量的数量。在第二步中，我们必须指示块系统更新它所管理的块的大小的知识；这发生在 <code>dsp.collect_sizes ()</code> 的调用中。

    BlockDynamicSparsityPattern dsp(2, 2); 
    dsp.block(0, 0).reinit(n_u, n_u); 
    dsp.block(1, 0).reinit(n_p, n_u); 
    dsp.block(0, 1).reinit(n_u, n_p); 
    dsp.block(1, 1).reinit(n_p, n_p); 
    dsp.collect_sizes(); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

// 我们以与非区块版本相同的方式使用压缩的区块稀疏模式，以创建稀疏模式，然后创建系统矩阵。

    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 

// 然后，我们必须以与块压缩稀疏度模式完全相同的方式调整解决方案和右侧向量的大小。

    solution.reinit(2); 
    solution.block(0).reinit(n_u); 
    solution.block(1).reinit(n_p); 
    solution.collect_sizes(); 

    system_rhs.reinit(2); 
    system_rhs.block(0).reinit(n_u); 
    system_rhs.block(1).reinit(n_p); 
    system_rhs.collect_sizes(); 
  } 
// @sect4{MixedLaplaceProblem::assemble_system}  

// 同样地，组装线性系统的函数在这个例子的介绍中已经讨论过很多了。在它的顶部，发生的是所有常见的步骤，此外，我们不仅为单元项分配正交和 <code>FEValues</code> 对象，而且还为面项分配。之后，我们为变量定义通常的缩写，并为本地矩阵和右手贡献分配空间，以及保存当前单元的全局自由度数的数组。

  template <int dim> 
  void MixedLaplaceProblem<dim>::assemble_system() 
  { 
    QGauss<dim>     quadrature_formula(degree + 2); 
    QGauss<dim - 1> face_quadrature_formula(degree + 2); 

    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 下一步是声明代表方程中源项、压力边界值和系数的对象。除了这些代表连续函数的对象外，我们还需要数组来保存它们在各个单元格（或面，对于边界值）的正交点的值。请注意，在系数的情况下，数组必须是矩阵的一种。

    const PrescribedSolution::RightHandSide<dim> right_hand_side; 
    const PrescribedSolution::PressureBoundaryValues<dim> 
                                            pressure_boundary_values; 
    const PrescribedSolution::KInverse<dim> k_inverse; 

    std::vector<double>         rhs_values(n_q_points); 
    std::vector<double>         boundary_values(n_face_q_points); 
    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 

// 最后，我们需要几个提取器，用来获取矢量值形状函数的速度和压力成分。它们的功能和使用在 @ref vector_valued报告中有详细描述。基本上，我们将把它们作为下面FEValues对象的下标：FEValues对象描述了形状函数的所有矢量分量，而在订阅后，它将只指速度（一组从零分量开始的 <code>dim</code> 分量）或压力（位于 <code>dim</code> 位置的标量分量）。

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

// 有了这些，我们就可以继续对所有单元进行循环。这个循环的主体已经在介绍中讨论过了，这里就不再做任何评论了。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        local_matrix = 0; 
        local_rhs    = 0; 

        right_hand_side.value_list(fe_values.get_quadrature_points(), 
                                   rhs_values); 
        k_inverse.value_list(fe_values.get_quadrature_points(), 
                             k_inverse_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q); 
              const double div_phi_i_u = fe_values[velocities].divergence(i, q); 
              const double phi_i_p     = fe_values[pressure].value(i, q); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const Tensor<1, dim> phi_j_u = 
                    fe_values[velocities].value(j, q); 
                  const double div_phi_j_u = 
                    fe_values[velocities].divergence(j, q); 
                  const double phi_j_p = fe_values[pressure].value(j, q); 

                  local_matrix(i, j) += 
                    (phi_i_u * k_inverse_values[q] * phi_j_u // 
                     - phi_i_p * div_phi_j_u                 // 
                     - div_phi_i_u * phi_j_p)                // 
                    * fe_values.JxW(q); 
                } 

              local_rhs(i) += -phi_i_p * rhs_values[q] * fe_values.JxW(q); 
            } 

        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary()) 
            { 
              fe_face_values.reinit(cell, face); 

              pressure_boundary_values.value_list( 
                fe_face_values.get_quadrature_points(), boundary_values); 

              for (unsigned int q = 0; q < n_face_q_points; ++q) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  local_rhs(i) += -(fe_face_values[velocities].value(i, q) * // 
                                    fe_face_values.normal_vector(q) *        // 
                                    boundary_values[q] *                     // 
                                    fe_face_values.JxW(q)); 
            } 

// 循环所有单元的最后一步是将局部贡献转移到全局矩阵和右手向量中。请注意，我们使用的接口与之前的例子完全相同，尽管我们现在使用的是块状矩阵和向量，而不是常规的。换句话说，对于外界来说，块对象具有与矩阵和向量相同的接口，但它们还允许访问单个块。

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(local_dof_indices[i], 
                              local_dof_indices[j], 
                              local_matrix(i, j)); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          system_rhs(local_dof_indices[i]) += local_rhs(i); 
      } 
  } 
// @sect3{Implementation of linear solvers and preconditioners}  

// 我们在这个例子中使用的线性求解器和预处理器已经在介绍中进行了详细的讨论。因此，我们在这里不再讨论我们的方法的原理，而只是对剩下的一些实现方面进行评论。

//  @sect4{MixedLaplace::solve}  

// 正如在介绍中所概述的那样，求解函数基本上由两个步骤组成。首先，我们必须形成涉及舒尔补数的第一个方程，并求解压力（解决方案的第一部分）。然后，我们可以从第二个方程（解的第0部分）中重构速度。

  template <int dim> 
  void MixedLaplaceProblem<dim>::solve() 
  { 

// 作为第一步，我们声明对矩阵的所有块状成分、右手边和我们将需要的解向量的引用。

    const auto &M = system_matrix.block(0, 0); 
    const auto &B = system_matrix.block(0, 1); 

    const auto &F = system_rhs.block(0); 
    const auto &G = system_rhs.block(1); 

    auto &U = solution.block(0); 
    auto &P = solution.block(1); 

// 然后，我们将创建相应的LinearOperator对象并创建 <code>op_M_inv</code> 运算器。

    const auto op_M = linear_operator(M); 
    const auto op_B = linear_operator(B); 

    ReductionControl         reduction_control_M(2000, 1.0e-18, 1.0e-10); 
    SolverCG<Vector<double>> solver_M(reduction_control_M); 
    PreconditionJacobi<SparseMatrix<double>> preconditioner_M; 

    preconditioner_M.initialize(M); 

    const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M); 

// 这样我们就可以声明舒尔补数  <code>op_S</code>  和近似舒尔补数  <code>op_aS</code>  。

    const auto op_S = transpose_operator(op_B) * op_M_inv * op_B; 
    const auto op_aS = 
      transpose_operator(op_B) * linear_operator(preconditioner_M) * op_B; 

// 我们现在从 <code>op_aS</code> 中创建一个预处理程序，应用固定数量的30次（便宜的）CG迭代。

    IterationNumberControl   iteration_number_control_aS(30, 1.e-18); 
    SolverCG<Vector<double>> solver_aS(iteration_number_control_aS); 

    const auto preconditioner_S = 
      inverse_operator(op_aS, solver_aS, PreconditionIdentity()); 

// 现在来看看第一个方程。它的右边是 $B^TM^{-1}F-G$  ，这就是我们在前几行计算的结果。然后我们用CG求解器和我们刚刚声明的预处理程序来解决第一个方程。

    const auto schur_rhs = transpose_operator(op_B) * op_M_inv * F - G; 

    SolverControl            solver_control_S(2000, 1.e-12); 
    SolverCG<Vector<double>> solver_S(solver_control_S); 

    const auto op_S_inv = inverse_operator(op_S, solver_S, preconditioner_S); 

    P = op_S_inv * schur_rhs; 

    std::cout << solver_control_S.last_step() 
              << " CG Schur complement iterations to obtain convergence." 
              << std::endl; 

// 得到压力后，我们可以计算速度。方程为 $MU=-BP+F$  ，我们通过首先计算右手边，然后与代表质量矩阵逆的对象相乘来解决这个问题。

    U = op_M_inv * (F - op_B * P); 
  } 
// @sect3{MixedLaplaceProblem class implementation (continued)}  
// @sect4{MixedLaplace::compute_errors}  

// 在我们处理完线性求解器和预处理器之后，我们继续实现我们的主类。特别是，下一个任务是计算我们数值解的误差，包括压力和速度。

// 为了计算解的误差，我们已经在  step-7  和  step-11  中介绍了  <code>VectorTools::integrate_difference</code>  函数。然而，在那里我们只处理了标量解，而在这里我们有一个矢量值的解，其组成部分甚至表示不同的量，并且可能有不同的收敛阶数（由于所使用的有限元的选择，这里不是这种情况，但在混合有限元应用中经常出现这种情况）。因此，我们要做的是 "掩盖 "我们感兴趣的成分。这很容易做到： <code>VectorTools::integrate_difference</code> 函数将一个指向权重函数的指针作为其参数之一（该参数默认为空指针，意味着单位权重）。我们要做的是传递一个函数对象，在我们感兴趣的成分中等于1，而在其他成分中等于0。例如，为了计算压力误差，我们应该传入一个函数，该函数在分量 <code>dim</code> 中代表单位值的常数向量，而对于速度，常数向量在第一个 <code>dim</code> 分量中应该是1，而在压力的位置是0。

// 在deal.II中， <code>ComponentSelectFunction</code> 正是这样做的：它想知道它要表示的函数应该有多少个向量分量（在我们的例子中，这将是 <code>dim+1</code> ，用于联合速度-压力空间），哪个个体或范围的分量应该等于1。因此，我们在函数的开头定义了两个这样的掩码，接下来是一个代表精确解的对象和一个向量，我们将在其中存储由 <code>integrate_difference</code> 计算的单元误差。

  template <int dim> 
  void MixedLaplaceProblem<dim>::compute_errors() const 
  { 
    const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1); 
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim), 
                                                     dim + 1); 

    PrescribedSolution::ExactSolution<dim> exact_solution; 
    Vector<double> cellwise_errors(triangulation.n_active_cells()); 

// 正如在 step-7 中已经讨论过的那样，我们必须认识到，不可能精确地整合误差。我们所能做的就是用正交法对这个积分进行近似。这实际上在这里提出了一个小小的转折：如果我们像人们可能倾向于做的那样天真地选择一个 <code>QGauss@<dim@>(degree+1)</code> 类型的对象（这就是我们用于积分线性系统的对象），就会发现误差非常小，根本不遵循预期的收敛曲线。现在的情况是，对于这里使用的混合有限元，高斯点恰好是超收敛点，其中的点误差要比其他地方小得多（而且收敛的阶数更高）。因此，这些点不是特别好的积分点。为了避免这个问题，我们只需使用梯形法则，并在每个坐标方向上迭代 <code>degree+2</code> 次（同样如 step-7 中的解释）。

    QTrapezoid<1>  q_trapez; 
    QIterated<dim> quadrature(q_trapez, degree + 2); 

// 有了这个，我们就可以让库计算出误差并将其输出到屏幕上。

    VectorTools::integrate_difference(dof_handler, 
                                      solution, 
                                      exact_solution, 
                                      cellwise_errors, 
                                      quadrature, 
                                      VectorTools::L2_norm, 
                                      &pressure_mask); 
    const double p_l2_error = 
      VectorTools::compute_global_error(triangulation, 
                                        cellwise_errors, 
                                        VectorTools::L2_norm); 

    VectorTools::integrate_difference(dof_handler, 
                                      solution, 
                                      exact_solution, 
                                      cellwise_errors, 
                                      quadrature, 
                                      VectorTools::L2_norm, 
                                      &velocity_mask); 
    const double u_l2_error = 
      VectorTools::compute_global_error(triangulation, 
                                        cellwise_errors, 
                                        VectorTools::L2_norm); 

    std::cout << "Errors: ||e_p||_L2 = " << p_l2_error 
              << ",   ||e_u||_L2 = " << u_l2_error << std::endl; 
  } 
// @sect4{MixedLaplace::output_results}  

// 最后一个有趣的函数是我们生成图形输出的函数。请注意，所有的速度分量都得到相同的解名 "u"。再加上使用 DataComponentInterpretation::component_is_part_of_vector ，这将导致 DataOut<dim>::write_vtu() 生成各个速度分量的矢量表示，更多信息请参见 step-22 或 @ref VVOutput 模块中的 "生成图形输出 "部分。最后，对于高阶元素来说，在图形输出中每个单元只显示一个双线性四边形似乎不合适。因此，我们生成大小为(度数+1)x(度数+1)的斑块来捕捉解决方案的全部信息内容。有关这方面的更多信息，请参见 step-7 的教程程序。

  template <int dim> 
  void MixedLaplaceProblem<dim>::output_results() const 
  { 
    std::vector<std::string> solution_names(dim, "u"); 
    solution_names.emplace_back("p"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      interpretation(dim, 
                     DataComponentInterpretation::component_is_part_of_vector); 
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 

    DataOut<dim> data_out; 
    data_out.add_data_vector(dof_handler, 
                             solution, 
                             solution_names, 
                             interpretation); 

    data_out.build_patches(degree + 1); 

    std::ofstream output("solution.vtu"); 
    data_out.write_vtu(output); 
  } 

//  @sect4{MixedLaplace::run}  

// 这是我们主类的最后一个函数。它唯一的工作是按照自然顺序调用其他函数。

  template <int dim> 
  void MixedLaplaceProblem<dim>::run() 
  { 
    make_grid_and_dofs(); 
    assemble_system(); 
    solve(); 
    compute_errors(); 
    output_results(); 
  } 
} // namespace Step20 
// @sect3{The <code>main</code> function}  

// 我们从  step-6  而不是  step-4  那里偷来的主函数。它几乎等同于 step-6 中的函数（当然，除了改变的类名），唯一的例外是我们将有限元空间的度数传递给混合拉普拉斯问题的构造函数（这里，我们使用零阶元素）。

int main() 
{ 
  try 
    { 
      using namespace Step20; 

      const unsigned int     fe_degree = 0; 
      MixedLaplaceProblem<2> mixed_laplace_problem(fe_degree); 
      mixed_laplace_problem.run(); 
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



