

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2010 - 2021 by the deal.II authors 
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
 * Authors: Andrea Bonito, Sebastian Pauletti. 
 */ 


// @sect3{Include files}  

// 如果你读过 step-4 和 step-7 ，你会认识到我们已经在那里使用了以下所有的包含文件。因此，我们不会在这里再次解释它们的含义。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 

#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/solver_control.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 

#include <fstream> 
#include <iostream> 

namespace Step38 
{ 
  using namespace dealii; 
// @sect3{The <code>LaplaceBeltramiProblem</code> class template}  

//这个类几乎与 step-4 中的 <code>LaplaceProblem</code> 类完全相似。

//本质上的区别是这样的。



// - 模板参数现在表示嵌入空间的维度，它不再与域和我们计算的三角形的维度相同。我们通过调用参数 @p spacedim, 并引入一个等于域的维度的常数 @p dim 来表明这一点--这里等于 <code>spacedim-1</code>  。

// - 所有具有几何特征的成员变量现在都需要知道它们自己的维度以及嵌入空间的维度。因此，我们需要指定它们的模板参数，一个是网格的维度 @p dim, ，另一个是嵌入空间的维度， @p spacedim.  这正是我们在 step-34 中所做的，请看那里有更深的解释。

// - 我们需要一个对象来描述从参考单元到三角形组成的单元所使用的哪种映射。从Mapping基类派生出来的类正是这样做的。在deal.II的大部分时间里，如果你不做任何事情，图书馆会假定你想要一个使用（双，三）线性映射的MappingQ1对象。在许多情况下，这就足够了，这就是为什么这些对象的使用大多是可选的：例如，如果你有一个二维空间中的多边形二维域，参考单元到三角形单元的双线性映射会产生该域的精确表示。如果你有一个弯曲的域，你可能想对那些位于域的边界的单元使用一个高阶映射--例如，这就是我们在 step-11 中所做的。然而，在这里我们有一个弯曲的域，而不仅仅是一个弯曲的边界，虽然我们可以用双线性映射的单元来近似它，但对所有单元使用高阶映射才是真正谨慎的。因此，这个类有一个MappingQ类型的成员变量；我们将选择映射的多项式程度等于计算中使用的有限元的多项式程度，以确保最佳近似，尽管这种等参数性不是必须的。

  template <int spacedim> 
  class LaplaceBeltramiProblem 
  { 
  public: 
    LaplaceBeltramiProblem(const unsigned degree = 2); 
    void run(); 

  private: 
    static constexpr unsigned int dim = spacedim - 1; 

    void make_grid_and_dofs(); 
    void assemble_system(); 
    void solve(); 
    void output_results() const; 
    void compute_error() const; 

    Triangulation<dim, spacedim> triangulation; 
    FE_Q<dim, spacedim>          fe; 
    DoFHandler<dim, spacedim>    dof_handler; 
    MappingQ<dim, spacedim>      mapping; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 
  }; 
// @sect3{Equation data}  

// 接下来，让我们定义描述问题的精确解和右手边的类。这与 step-4 和 step-7 相类似，在那里我们也定义了此类对象。鉴于介绍中的讨论，实际的公式应该是不言自明的。值得关注的一点是，我们是如何使用一般模板的明确特化，分别定义2D和3D情况下的值和梯度函数的。另一种方法是定义通用模板，并为空间维度的每个可能的值设置一个 <code>switch</code> 语句（或一串 <code>if</code> s）。

  template <int dim> 
  class Solution : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual Tensor<1, dim> 
    gradient(const Point<dim> & p, 
             const unsigned int component = 0) const override; 
  }; 

  template <> 
  double Solution<2>::value(const Point<2> &p, const unsigned int) const 
  { 
    return (-2. * p(0) * p(1)); 
  } 

  template <> 
  Tensor<1, 2> Solution<2>::gradient(const Point<2> &p, 
                                     const unsigned int) const 
  { 
    Tensor<1, 2> return_value; 
    return_value[0] = -2. * p(1) * (1 - 2. * p(0) * p(0)); 
    return_value[1] = -2. * p(0) * (1 - 2. * p(1) * p(1)); 

    return return_value; 
  } 

  template <> 
  double Solution<3>::value(const Point<3> &p, const unsigned int) const 
  { 
    return (std::sin(numbers::PI * p(0)) * std::cos(numbers::PI * p(1)) * 
            exp(p(2))); 
  } 

  template <> 
  Tensor<1, 3> Solution<3>::gradient(const Point<3> &p, 
                                     const unsigned int) const 
  { 
    using numbers::PI; 

    Tensor<1, 3> return_value; 

    return_value[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
    return_value[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
    return_value[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 

    return return_value; 
  } 

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <> 
  double RightHandSide<2>::value(const Point<2> &p, 
                                 const unsigned int /*component*/) const 
  { 
    return (-8. * p(0) * p(1)); 
  } 

  template <> 
  double RightHandSide<3>::value(const Point<3> &p, 
                                 const unsigned int /*component*/) const 
  { 
    using numbers::PI; 

    Tensor<2, 3> hessian; 

    hessian[0][0] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
    hessian[1][1] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
    hessian[2][2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 

    hessian[0][1] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
    hessian[1][0] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 

    hessian[0][2] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
    hessian[2][0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 

    hessian[1][2] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
    hessian[2][1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 

    Tensor<1, 3> gradient; 
    gradient[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
    gradient[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
    gradient[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 

    Point<3> normal = p; 
    normal /= p.norm(); 

    return (-trace(hessian) + 2 * (gradient * normal) + 
            (hessian * normal) * normal); 
  } 
// @sect3{Implementation of the <code>LaplaceBeltramiProblem</code> class}  

// 如果你知道  step-4  ，程序的其余部分实际上是很不引人注目的。我们的第一步是定义构造函数，设置有限元和映射的多项式程度，并将DoF处理程序与三角形关联。

  template <int spacedim> 
  LaplaceBeltramiProblem<spacedim>::LaplaceBeltramiProblem( 
    const unsigned degree) 
    : fe(degree) 
    , dof_handler(triangulation) 
    , mapping(degree) 
  {} 
// @sect4{LaplaceBeltramiProblem::make_grid_and_dofs}  

// 下一步是创建网格，分配自由度，并设置描述线性系统的各种变量。所有这些步骤都是标准的，只有如何创建一个描述曲面的网格除外。我们可以为我们感兴趣的领域生成一个网格，用一个网格生成器生成一个三角形，然后用GridIn类将其读入。或者，就像我们在这里做的那样，我们使用GridGenerator命名空间的设施来生成网格。

// 具体来说，我们要做的是这样的（在下面的大括号中）：我们使用 <code>spacedim</code> 函数为半圆盘（2D）或半球（3D）生成一个 GridGenerator::half_hyper_ball 维度的网格。这个函数将位于圆盘/球周边的所有面的边界指标设置为零，而在将整个圆盘/球分成两半的直线部分设置为零。下一步是主要的一点。 GridGenerator::extract_boundary_mesh 函数创建的网格是由那些作为前一个网格的面的单元组成的，也就是说，它描述了原始（体积）网格的<i>surface</i>单元。然而，我们不需要所有的面：只需要那些在圆盘或球的周边，边界指示器为零的面；我们可以使用一组边界指示器来选择这些单元，并传递给 GridGenerator::extract_boundary_mesh. 。

// 有一点需要提及。为了在流形是弯曲的情况下适当地细化表面网格（类似于细化与弯曲边界相邻的单元面），三角形必须要有一个对象附加在上面，描述新顶点应该位于何处。如果你不附加这样的边界对象，它们将位于现有顶点之间的中间位置；如果你有一个具有直线边界的域（例如多边形），这是很合适的，但如果像这里一样，流形具有曲率，则不合适。因此，为了让事情正常进行，我们需要将流形对象附加到我们的（表面）三角形上，其方式与我们在1d中为边界所做的大致相同。我们创建这样一个对象，并将其附加到三角剖面上。

// 创建网格的最后一步是对其进行多次细化。该函数的其余部分与之前的教程程序中相同。

  template <int spacedim> 
  void LaplaceBeltramiProblem<spacedim>::make_grid_and_dofs() 
  { 
    { 
      Triangulation<spacedim> volume_mesh; 
      GridGenerator::half_hyper_ball(volume_mesh); 

      std::set<types::boundary_id> boundary_ids; 
      boundary_ids.insert(0); 

      GridGenerator::extract_boundary_mesh(volume_mesh, 
                                           triangulation, 
                                           boundary_ids); 
    } 
    triangulation.set_all_manifold_ids(0); 
    triangulation.set_manifold(0, SphericalManifold<dim, spacedim>()); 

    triangulation.refine_global(4); 

    std::cout << "Surface mesh has " << triangulation.n_active_cells() 
              << " cells." << std::endl; 

    dof_handler.distribute_dofs(fe); 

    std::cout << "Surface mesh has " << dof_handler.n_dofs() 
              << " degrees of freedom." << std::endl; 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 
// @sect4{LaplaceBeltramiProblem::assemble_system}  

// 下面是这个程序的中心函数，即组装与表面拉普拉斯（Laplace-Beltrami算子）相对应的矩阵。也许令人惊讶的是，它实际上与例如在  step-4  中讨论的普通拉普拉斯算子看起来完全一样。关键是 FEValues::shape_grad() 函数发挥了魔力：它返回 $i$ 第1个形状函数在 $q$ 第1个正交点的表面梯度 $\nabla_K \phi_i(x_q)$ 。其余的也不需要任何改变。

  template <int spacedim> 
  void LaplaceBeltramiProblem<spacedim>::assemble_system() 
  { 
    system_matrix = 0; 
    system_rhs    = 0; 

    const QGauss<dim>       quadrature_formula(2 * fe.degree); 
    FEValues<dim, spacedim> fe_values(mapping, 
                                      fe, 
                                      quadrature_formula, 
                                      update_values | update_gradients | 
                                        update_quadrature_points | 
                                        update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<double>                  rhs_values(n_q_points); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    RightHandSide<spacedim> rhs; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        fe_values.reinit(cell); 

        rhs.value_list(fe_values.get_quadrature_points(), rhs_values); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
              cell_matrix(i, j) += fe_values.shape_grad(i, q_point) * 
                                   fe_values.shape_grad(j, q_point) * 
                                   fe_values.JxW(q_point); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            cell_rhs(i) += fe_values.shape_value(i, q_point) * 
                           rhs_values[q_point] * fe_values.JxW(q_point); 

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

    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values( 
      mapping, dof_handler, 0, Solution<spacedim>(), boundary_values); 

    MatrixTools::apply_boundary_values( 
      boundary_values, system_matrix, solution, system_rhs, false); 
  } 

//  @sect4{LaplaceBeltramiProblem::solve}  

// 下一个函数是解决线性系统的函数。在这里，也不需要做任何改变。

  template <int spacedim> 
  void LaplaceBeltramiProblem<spacedim>::solve() 
  { 
    SolverControl solver_control(solution.size(), 1e-7 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 
  } 

//  @sect4{LaplaceBeltramiProblem::output_result}  

// 这是一个从解决方案中生成图形输出的函数。它的大部分都是模板代码，但有两点值得指出。



// -  DataOut::add_data_vector() 函数可以接受两种向量。  一种是之前通过 DataOut::attach_dof_handler(); 连接的DoFHandler对象定义的每个自由度有一个值的向量，另一种是三角测量的每个单元有一个值的向量，例如，输出每个单元的估计误差。通常，DataOut类知道如何区分这两种向量：自由度几乎总是比单元格多，所以我们可以通过两种向量的长度来区分。我们在这里也可以这样做，但只是因为我们很幸运：我们使用了一个半球体。如果我们用整个球体作为域和 $Q_1$ 元素，我们将有相同数量的单元格作为顶点，因此这两种向量将有相同数量的元素。为了避免由此产生的混乱，我们必须告诉 DataOut::add_data_vector() 函数我们有哪种矢量。DoF数据。这就是该函数的第三个参数的作用。

// -  DataOut::build_patches() 函数可以生成细分每个单元的输出，这样可视化程序可以更好地解决弯曲流形或更高的多项式程度的形状函数。在这里，我们在每个坐标方向上对每个元素进行细分，细分的次数与使用的有限元的多项式程度相同。

  template <int spacedim> 
  void LaplaceBeltramiProblem<spacedim>::output_results() const 
  { 
    DataOut<dim, spacedim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, 
                             "solution", 
                             DataOut<dim, spacedim>::type_dof_data); 
    data_out.build_patches(mapping, mapping.get_degree()); 

    const std::string filename = 
      "solution-" + std::to_string(spacedim) + "d.vtk"; 
    std::ofstream output(filename); 
    data_out.write_vtk(output); 
  } 

//  @sect4{LaplaceBeltramiProblem::compute_error}  

// 这是最后一块功能：我们要计算数值解的误差。它是之前在  step-7  中展示和讨论的代码的逐字复制。正如介绍中提到的， <code>Solution</code> 类提供了解决方案的（切向）梯度。为了避免只评估超收敛点的误差，我们选择一个足够高阶的正交规则。

  template <int spacedim> 
  void LaplaceBeltramiProblem<spacedim>::compute_error() const 
  { 
    Vector<float> difference_per_cell(triangulation.n_active_cells()); 
    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      solution, 
                                      Solution<spacedim>(), 
                                      difference_per_cell, 
                                      QGauss<dim>(2 * fe.degree + 1), 
                                      VectorTools::H1_norm); 

    double h1_error = VectorTools::compute_global_error(triangulation, 
                                                        difference_per_cell, 
                                                        VectorTools::H1_norm); 
    std::cout << "H1 error = " << h1_error << std::endl; 
  } 

//  @sect4{LaplaceBeltramiProblem::run}  

// 最后一个函数提供了顶层的逻辑。它的内容是不言自明的。

  template <int spacedim> 
  void LaplaceBeltramiProblem<spacedim>::run() 
  { 
    make_grid_and_dofs(); 
    assemble_system(); 
    solve(); 
    output_results(); 
    compute_error(); 
  } 
} // namespace Step38 
// @sect3{The main() function}  

// 该程序的其余部分由 <code>main()</code> 函数占据。它完全遵循首次在 step-6 中介绍的一般布局，并在随后的所有教程程序中使用。

int main() 
{ 
  try 
    { 
      using namespace Step38; 

      LaplaceBeltramiProblem<3> laplace_beltrami; 
      laplace_beltrami.run(); 
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


