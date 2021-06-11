CCTest_file/step-65.cc

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
 * This tutorial program was contributed by Martin Kronbichler 
 */ 


// @sect3{Include files}  

// 本教程的包含文件与  step-6  中的基本相同。重要的是，我们将使用的TransfiniteInterpolationManifold类是由`deal.II/grid/manifold_lib.h`提供。

#include <deal.II/base/timer.h> 

#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q_generic.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/vector_tools.h> 

#include <fstream> 

// 唯一的新include文件是MappingQCache类的文件。

#include <deal.II/fe/mapping_q_cache.h> 

namespace Step65 
{ 
  using namespace dealii; 
// @sect3{Analytical solution and coefficient}  

// 在这个教程程序中，我们要解决泊松方程，其系数沿半径为0.5的球体跳跃，并使用一个恒定的右手边值 $f(\mathbf{x}) = -3$  。（这个设置与 step-5 和 step-6 相似，但系数和右手边的具体数值不同）。由于系数的跳跃，分析解必须有一个结点，即系数从一个值切换到另一个值。为了保持简单，我们选择了一个在所有分量中都是二次的分析解，即在半径为0.5的球中为 $u(x,y,z) = x^2 + y^2 + z^2$ ，在域的外部为 $u(x,y,z) = 0.1(x^2 + y^2 + z^2) + 0.25-0.025$ 。这个分析解在内球的系数为0.5，外球的系数为5的情况下与右手边兼容。它也是沿着半径为0.5的圆连续的。

  template <int dim> 
  class ExactSolution : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      if (p.norm_square() < 0.25) 
        return p.norm_square(); 
      else 
        return 0.1 * p.norm_square() + (0.25 - 0.025); 
    } 

    virtual Tensor<1, dim> 
    gradient(const Point<dim> &p, 
             const unsigned int /*component*/ = 0) const override 
    { 
      if (p.norm_square() < 0.25) 
        return 2. * p; 
      else 
        return 0.2 * p; 
    } 
  }; 

  template <int dim> 
  double coefficient(const Point<dim> &p) 
  { 
    if (p.norm_square() < 0.25) 
      return 0.5; 
    else 
      return 5.0; 
  } 

//  @sect3{The PoissonProblem class}  

// 泊松问题的实现与我们在  step-5  教程中使用的非常相似。两个主要的区别是，我们向程序中的各个步骤传递了一个映射对象，以便在两种映射表示法之间进行切换，正如介绍中所解释的那样，还有一个`计时器'对象（TimerOutput类型），将用于测量各种情况下的运行时间。(映射对象的概念在 step-10 和 step-11 中首次提出，如果你想查一下这些类的用途的话)。

  template <int dim> 
  class PoissonProblem 
  { 
  public: 
    PoissonProblem(); 
    void run(); 

  private: 
    void create_grid(); 
    void setup_system(const Mapping<dim> &mapping); 
    void assemble_system(const Mapping<dim> &mapping); 
    void solve(); 
    void postprocess(const Mapping<dim> &mapping); 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> constraints; 
    SparsityPattern           sparsity_pattern; 
    SparseMatrix<double>      system_matrix; 
    Vector<double>            solution; 
    Vector<double>            system_rhs; 

    TimerOutput timer; 
  }; 

// 在构造函数中，我们设置了定时器对象来记录墙的时间，但在正常执行过程中是安静的。我们将在 `PoissonProblem::run()` 函数中查询它的计时细节。此外，我们为正在使用的有限元选择了一个相对较高的多项式三度。

  template <int dim> 
  PoissonProblem<dim>::PoissonProblem() 
    : fe(3) 
    , dof_handler(triangulation) 
    , timer(std::cout, TimerOutput::never, TimerOutput::wall_times) 
  {} 

//  @sect3{Grid creation and initialization of the manifolds}  

// 接下来的函数介绍了TransfiniteInterpolationManifold的典型用法。第一步是创建所需的网格，这可以通过GridGenerator的两个网格的组合来完成。内球网格是很简单的。我们以原点为中心运行 GridGenerator::hyper_cube() ，半径为0.5（第三个函数参数）。第二个网格更有趣，构建方法如下。我们希望有一个在内部是球形的，但在外表面是平的网格。此外，内球的网格拓扑结构应该与外球的网格兼容，即它们的顶点重合，这样才能使两个网格合并起来。从 GridGenerator::hyper_shell 出来的网格满足了内侧的要求，如果它是用 $2d$ 的粗大单元创建的（在3D中我们将使用6个粗大单元）&ndash；这与球的边界面的单元数量相同。对于外表面，我们利用这样一个事实：没有流形附着的壳表面的6个面将退化为立方体的表面。我们仍然缺少的是外壳边界的半径。由于我们想要一个范围为 $[-1, 1]$ 的立方体，而6单元壳将其8个外顶点放在8条对角线上，我们必须将点 $(\pm 1, \pm 1, \pm 1)$ 转化为半径。显然，在 $d$ 维度上，半径必须是 $\sqrt{d}$ ，也就是说，对于我们要考虑的三维情况，半径是 $\sqrt{3}$ 。

// 这样，我们就有了一个计划。在创建了球的内部三角形和外壳的三角形之后，我们将这两个网格合并，但是将GridGenerator中的函数可能从产生的三角形中设置的所有流形移除，以确保我们对流形有充分的控制。特别是，我们希望在细化过程中在边界上添加的额外点能够遵循平坦的流形描述。为了开始添加更合适的流形ID的过程，我们给所有的网格实体（单元、面、线）分配流形ID 0，这些实体以后将与TransfiniteInterpolationManifold相关联。然后，我们必须识别沿着半径为0.5的球体的面和线，并给它们标记一个不同的流形ID，以便随后给这些面和线分配一个SphericalManifold。由于我们在调用 GridGenerator::hyper_ball(), 后丢弃了所有预先存在的流形，我们手动检查了网格的单元格和所有的面。如果四个顶点的半径都是0.5，我们就在球体上找到了一个面，或者像我们在程序中写的那样，有  $r^2-0.25 \approx 0$  。注意，我们调用`cell->face(f)->set_all_manifold_ids(1)`来设置面和周围线上的流形id。此外，我们希望通过一个材料ID来区分球内和球外的单元，以便于可视化，对应于介绍中的图片。

  template <int dim> 
  void PoissonProblem<dim>::create_grid() 
  { 
    Triangulation<dim> tria_inner; 
    GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5); 

    Triangulation<dim> tria_outer; 
    GridGenerator::hyper_shell( 
      tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim); 

    GridGenerator::merge_triangulations(tria_inner, tria_outer, triangulation); 

    triangulation.reset_all_manifolds(); 
    triangulation.set_all_manifold_ids(0); 

    for (const auto &cell : triangulation.cell_iterators()) 
      { 
        for (const auto &face : cell->face_iterators()) 
          { 
            bool face_at_sphere_boundary = true; 
            for (const auto v : face->vertex_indices()) 
              { 
                if (std::abs(face->vertex(v).norm_square() - 0.25) > 1e-12) 
                  face_at_sphere_boundary = false; 
              } 
            if (face_at_sphere_boundary) 
              face->set_all_manifold_ids(1); 
          } 
        if (cell->center().norm_square() < 0.25) 
          cell->set_material_id(1); 
        else 
          cell->set_material_id(0); 
      } 

// 有了所有单元格、面和线的适当标记，我们可以将流形对象附加到这些数字上。流形ID为1的实体将得到一个球形流形，而流形ID为0的其他实体将被分配到TransfiniteInterpolationManifold。正如介绍中提到的，我们必须通过调用 TransfiniteInterpolationManifold::initialize() 显式初始化当前网格的流形，以获取粗略的网格单元和连接到这些单元边界的流形。我们还注意到，我们在这个函数中本地创建的流形对象是允许超出范围的（就像它们在函数范围结束时那样），因为Triangulation对象在内部复制它们。

// 在连接了所有的流形之后，我们最后将去细化网格几次，以创建一个足够大的测试案例。

    triangulation.set_manifold(1, SphericalManifold<dim>()); 

    TransfiniteInterpolationManifold<dim> transfinite_manifold; 
    transfinite_manifold.initialize(triangulation); 
    triangulation.set_manifold(0, transfinite_manifold); 

    triangulation.refine_global(9 - 2 * dim); 
  } 

//  @sect3{Setup of data structures}  

// 下面的函数在其他教程中是众所周知的，它枚举了自由度，创建了一个约束对象并为线性系统设置了一个稀疏矩阵。唯一值得一提的是，该函数接收了一个映射对象的引用，然后我们将其传递给 VectorTools::interpolate_boundary_values() 函数，以确保我们的边界值在用于装配的高阶网格上被评估。在本例中，这并不重要，因为外表面是平的，但对于弯曲的外单元，这将导致边界值的更精确的近似。

  template <int dim> 
  void PoissonProblem<dim>::setup_system(const Mapping<dim> &mapping) 
  { 
    dof_handler.distribute_dofs(fe); 
    std::cout << "   Number of active cells:       " 
              << triangulation.n_global_active_cells() << std::endl; 
    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 

    { 
      TimerOutput::Scope scope(timer, "Compute constraints"); 

      constraints.clear(); 

      DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
      VectorTools::interpolate_boundary_values( 
        mapping, dof_handler, 0, ExactSolution<dim>(), constraints); 

      constraints.close(); 
    } 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 

    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 
// @sect3{Assembly of the system matrix and right hand side}  

// 组装线性系统的函数在前面的教程程序中也是众所周知的。有一点需要注意的是，我们将正交点的数量设置为多项式的度数加2，而不是像其他大多数教程中的度数加1。这是因为我们期望有一些额外的精度，因为映射也涉及到比解的多项式多一度的程度。

// 汇编中唯一有点不寻常的代码是我们计算单元格矩阵的方式。我们没有使用正交点索引、行和矩阵列的三个嵌套循环，而是首先收集形状函数的导数，乘以系数和积分因子`JxW`的乘积的平方根，放在一个单独的矩阵`partial_matrix`中。为了计算单元矩阵，我们在 "partial_matrix.mTmult(cell_matrix, partial_matrix); "一行中执行 "cell_matrix = partial_matrix * transpose(partial_matrix)"。为了理解这一点，我们要知道矩阵与矩阵的乘法是对`partial_matrix`的各列进行求和。如果我们用 
// $a(\mathbf{x}_q)$ 表示系数，临时矩阵的条目是 $\sqrt{\text{det}(J) w_q a(x)} \frac{\partial \varphi_i(\boldsymbol
//  \xi_q)}{\partial x_k}$ 。如果我们将该矩阵的第<i>i</i>行与第<i>j</i>列相乘，我们计算出一个涉及 $\sum_q \sum_{k=1}^d \sqrt{\text{det}(J) w_q a(x)} \frac{\partial
//  \varphi_i(\boldsymbol \xi_q)}{\partial x_k} \sqrt{\text{det}(J) w_q a(x)}
//  \frac{\partial \varphi_j(\boldsymbol \xi_q)}{\partial x_k} = \sum_q
//  \sum_{k=1}^d\text{det}(J) w_q a(x)\frac{\partial \varphi_i(\boldsymbol
//  \xi_q)}{\partial x_k} \frac{\partial \varphi_j(\boldsymbol
//  \xi_q)}{\partial x_k}$ 的嵌套和，这正是拉普拉斯方程的双线性形式所需的条款。

// 选择这种有点不寻常的方案的原因是由于计算三维中相对较高的多项式程度的单元矩阵所涉及的繁重工作。由于我们想在这个教程程序中强调映射的成本，我们最好以优化的方式进行装配，以便不追逐已经被社区解决的瓶颈。矩阵-矩阵乘法是HPC背景下最好的优化内核之一， FullMatrix::mTmult() 函数将调用到那些优化的BLAS函数。如果用户在配置deal.II时提供了一个好的BLAS库（如OpenBLAS或英特尔的MKL），那么单元矩阵的计算将执行到接近处理器的峰值算术性能。顺便提一下，尽管有优化的矩阵-矩阵乘法，但目前的策略在复杂性方面是次优的，因为要做的工作与 $(p+1)^9$ 度 $p$ 的运算成正比（这也适用于用FEValues的通常评估）。我们可以通过利用形状函数的张量乘积结构，用 $\mathcal O((p+1)^7)$ 的操作来计算单元格矩阵，就像交易二中的无矩阵框架那样。我们参考 step-37 和张量积感知评估器FEEvaluation的文档，以了解如何实现更有效的单元矩阵计算的细节。

  template <int dim> 
  void PoissonProblem<dim>::assemble_system(const Mapping<dim> &mapping) 
  { 
    TimerOutput::Scope scope(timer, "Assemble linear system"); 

    const QGauss<dim> quadrature_formula(fe.degree + 2); 
    FEValues<dim>     fe_values(mapping, 
                            fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 
    FullMatrix<double> partial_matrix(dofs_per_cell, dim * n_q_points); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_rhs = 0.; 
        fe_values.reinit(cell); 

        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) 
          { 
            const double current_coefficient = 
              coefficient(fe_values.quadrature_point(q_index)); 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              { 
                for (unsigned int d = 0; d < dim; ++d) 
                  partial_matrix(i, q_index * dim + d) = 
                    std::sqrt(fe_values.JxW(q_index) * current_coefficient) * 
                    fe_values.shape_grad(i, q_index)[d]; 
                cell_rhs(i) += 
                  (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
                   (-dim) *                            // f(x_q) 
                   fe_values.JxW(q_index));            // dx 
              } 
          } 

        partial_matrix.mTmult(cell_matrix, partial_matrix); 

        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global( 
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 

//  @sect3{Solution of the linear system}  

// 对于线性系统的求解，我们选择一个简单的雅可比条件共轭梯度求解器，类似于早期教程中的设置。

  template <int dim> 
  void PoissonProblem<dim>::solve() 
  { 
    TimerOutput::Scope scope(timer, "Solve linear system"); 

    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> solver(solver_control); 

    PreconditionJacobi<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix); 

    solver.solve(system_matrix, solution, system_rhs, preconditioner); 
    constraints.distribute(solution); 

    std::cout << "   Number of solver iterations:  " 
              << solver_control.last_step() << std::endl; 
  } 

//  @sect3{Output of the solution and computation of errors}  

// 在下一个函数中，我们对解决方案做了各种后处理步骤，所有这些步骤都以这种或那种方式涉及映射。

// 我们做的第一个操作是把解决方案以及材料ID写到VTU文件中。这与其他许多教程程序中的做法类似。这个教程程序中提出的新内容是，我们要确保写到文件中用于可视化的数据实际上是deal.II内部使用的数据的忠实代表。这是因为大多数可视化数据格式只用顶点坐标表示单元，但没有办法表示deal.II中使用高阶映射时的曲线边界--换句话说，你在可视化工具中看到的东西实际上不是你正在计算的东西。顺带一提，在使用高阶形状函数时也是如此。大多数可视化工具只呈现双线性/三线性的表示。这在 DataOut::build_patches().) 中有详细的讨论。

// 所以我们需要确保高阶表示被写入文件中。我们需要考虑两个特别的话题。首先，我们通过 DataOutBase::VtkFlags 告诉DataOut对象，我们打算把元素的细分解释为高阶拉格朗日多项式，而不是双线性补丁的集合。最近的可视化程序，如ParaView 5.5版或更新的程序，然后可以呈现高阶解决方案（更多细节见<a
//  href="https:github.com/dealii/dealii/wiki/Notes-on-visualizing-high-order-output">wiki
//  page</a>）。其次，我们需要确保映射被传递给 DataOut::build_patches() 方法。最后，DataOut类默认只打印<i>boundary</i>单元的曲面，所以我们需要确保通过映射将内部单元也打印成曲面。

  template <int dim> 
  void PoissonProblem<dim>::postprocess(const Mapping<dim> &mapping) 
  { 
    { 
      TimerOutput::Scope scope(timer, "Write output"); 

      DataOut<dim> data_out; 

      DataOutBase::VtkFlags flags; 
      flags.write_higher_order_cells = true; 
      data_out.set_flags(flags); 

      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution, "solution"); 

      Vector<double> material_ids(triangulation.n_active_cells()); 
      for (const auto &cell : triangulation.active_cell_iterators()) 
        material_ids[cell->active_cell_index()] = cell->material_id(); 
      data_out.add_data_vector(material_ids, "material_ids"); 

      data_out.build_patches(mapping, 
                             fe.degree, 
                             DataOut<dim>::curved_inner_cells); 

      std::ofstream file( 
        ("solution-" + 
         std::to_string(triangulation.n_global_levels() - 10 + 2 * dim) + 
         ".vtu") 
          .c_str()); 

      data_out.write_vtu(file); 
    } 

// 后处理函数的下一个操作是对照分析解计算 $L_2$ 和 $H^1$ 误差。由于分析解是一个二次多项式，我们期望在这一点上得到一个非常准确的结果。如果我们是在一个具有平面面的简单网格上求解，并且系数的跳动与单元间的面对齐，那么我们会期望数值结果与分析解相吻合，直至舍去精度。然而，由于我们使用的是跟随球体的变形单元，这些单元只能由4度的多项式跟踪（比有限元的度数多一个），我们会发现在 $10^{-7}$ 附近有一个误差。我们可以通过增加多项式的度数或细化网格来获得更多的精度。

    { 
      TimerOutput::Scope scope(timer, "Compute error norms"); 

      Vector<double> norm_per_cell_p(triangulation.n_active_cells()); 

      VectorTools::integrate_difference(mapping, 
                                        dof_handler, 
                                        solution, 
                                        ExactSolution<dim>(), 
                                        norm_per_cell_p, 
                                        QGauss<dim>(fe.degree + 2), 
                                        VectorTools::L2_norm); 
      std::cout << "   L2 error vs exact solution:   " 
                << norm_per_cell_p.l2_norm() << std::endl; 

      VectorTools::integrate_difference(mapping, 
                                        dof_handler, 
                                        solution, 
                                        ExactSolution<dim>(), 
                                        norm_per_cell_p, 
                                        QGauss<dim>(fe.degree + 2), 
                                        VectorTools::H1_norm); 
      std::cout << "   H1 error vs exact solution:   " 
                << norm_per_cell_p.l2_norm() << std::endl; 
    } 

// 我们在这里做的最后一个后处理操作是用KellyErrorEstimator计算出一个误差估计。我们使用了与 step-6 教程程序中完全相同的设置，只是我们还交出了映射，以确保误差是沿着曲线元素评估的，与程序的其余部分一致。然而，我们并没有真正使用这里的结果来驱动网格适应步骤（会沿着球体细化材料界面周围的网格），因为这里的重点是这个操作的成本。

    { 
      TimerOutput::Scope scope(timer, "Compute error estimator"); 

      Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
      KellyErrorEstimator<dim>::estimate( 
        mapping, 
        dof_handler, 
        QGauss<dim - 1>(fe.degree + 1), 
        std::map<types::boundary_id, const Function<dim> *>(), 
        solution, 
        estimated_error_per_cell); 
      std::cout << "   Max cell-wise error estimate: " 
                << estimated_error_per_cell.linfty_norm() << std::endl; 
    } 
  } 

//  @sect3{The PoissonProblem::run() function}  

// 最后，我们定义了`run()`函数，控制我们如何执行这个程序（由main()函数以常规方式调用）。我们首先调用`create_grid()`函数，用适当的流形设置我们的几何体。然后我们运行两个求解器链的实例，从方程的设置开始，组装线性系统，用一个简单的迭代求解器求解，以及上面讨论的后处理。这两个实例在使用映射的方式上有所不同。第一个使用传统的MappingQGeneric映射对象，我们将其初始化为比有限元多一级的程度；毕竟，我们期望几何表示是瓶颈，因为分析解只是二次多项式。实际上，事情在相当程度上是相互关联的，因为实坐标中多项式的评估涉及到高阶多项式的映射，而高阶多项式代表一些光滑的有理函数。因此，高阶多项式还是有回报的，所以进一步增加映射的度数是没有意义的)。一旦第一遍完成，我们就让定时器打印出各个阶段的计算时间的摘要。

  template <int dim> 
  void PoissonProblem<dim>::run() 
  { 
    create_grid(); 

    { 
      std::cout << std::endl 
                << "====== Running with the basic MappingQGeneric class ====== " 
                << std::endl 
                << std::endl; 

      MappingQGeneric<dim> mapping(fe.degree + 1); 
      setup_system(mapping); 
      assemble_system(mapping); 
      solve(); 
      postprocess(mapping); 

      timer.print_summary(); 
      timer.reset(); 
    } 

// 对于第二个实例，我们转而设置了MappingQCache类。它的使用非常简单。在构建好它之后（考虑到我们希望它在其他情况下显示正确的度数功能，所以用度数），我们通过 MappingQCache::initialize() 函数填充缓存。在这个阶段，我们为缓存指定我们想要使用的映射（很明显，与之前的MappingQGeneric相同，以便重复相同的计算），然后再次运行相同的函数，现在交出修改后的映射。最后，我们再次打印重置后的累计壁挂时间，看看这些时间与原来的设置相比如何。

    { 
      std::cout 
        << "====== Running with the optimized MappingQCache class ====== " 
        << std::endl 
        << std::endl; 

      MappingQCache<dim> mapping(fe.degree + 1); 
      { 
        TimerOutput::Scope scope(timer, "Initialize mapping cache"); 
        mapping.initialize(MappingQGeneric<dim>(fe.degree + 1), triangulation); 
      } 
      std::cout << "   Memory consumption cache:     " 
                << 1e-6 * mapping.memory_consumption() << " MB" << std::endl; 

      setup_system(mapping); 
      assemble_system(mapping); 
      solve(); 
      postprocess(mapping); 

      timer.print_summary(); 
    } 
  } 
} // namespace Step65 

int main() 
{ 
  Step65::PoissonProblem<3> test_program; 
  test_program.run(); 
  return 0; 
} 


