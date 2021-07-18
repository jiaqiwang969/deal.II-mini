

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2018 - 2021 by the deal.II authors 
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
 *      Author: Zhuoran Wang, Colorado State University, 2018 
 */ 


// @sect3{Include files}  这个程序是基于 step-7  、 step-20  和  step-51  ，所以下面的头文件大部分是熟悉的。我们需要以下文件，其中只有导入FE_DGRaviartThomas类的文件（即`deal.II/fe/fe_dg_vector.h`）是真正的新文件；FE_DGRaviartThomas实现了介绍中讨论的 "破碎 "Raviart-Thomas空间。

#include <deal.II/base/quadrature.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/tensor_function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/point.h> 
#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_raviart_thomas.h> 
#include <deal.II/fe/fe_dg_vector.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_face.h> 
#include <deal.II/fe/component_mask.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/data_out_faces.h> 

#include <fstream> 
#include <iostream> 

// 我们的第一步，像往常一样，是把所有与本教程程序有关的东西放到自己的命名空间中。

namespace Step61 
{ 
  using namespace dealii; 
// @sect3{The WGDarcyEquation class template}  

// 这是本程序的主类。我们将使用弱加勒金（WG）方法求解内部和面上的数值压力，并计算出压力的 $L_2$ 误差。在后处理步骤中，我们还将计算速度和通量的 $L_2$  误差。

// 该类的结构与以前的教程程序没有根本的不同，所以除了一个例外，没有必要对细节进行评论。该类有一个成员变量`fe_dgrt`，对应于介绍中提到的 "破碎 "的Raviart-Thomas空间。还有一个与之匹配的`dof_handler_dgrt`，表示从这个元素创建的有限元场的全局枚举，还有一个向量`darcy_velocity`，用于保持这个场的节点值。在求解压力后，我们将使用这三个变量来计算一个后处理的速度场，然后我们可以对其进行误差评估，并将其输出用于可视化。

  template <int dim> 
  class WGDarcyEquation 
  { 
  public: 
    WGDarcyEquation(const unsigned int degree); 
    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void compute_postprocessed_velocity(); 
    void compute_velocity_errors(); 
    void compute_pressure_error(); 
    void output_results() const; 

    Triangulation<dim> triangulation; 

    FESystem<dim>   fe; 
    DoFHandler<dim> dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    FE_DGRaviartThomas<dim> fe_dgrt; 
    DoFHandler<dim>         dof_handler_dgrt; 
    Vector<double>          darcy_velocity; 
  }; 

//  @sect3{Right hand side, boundary values, and exact solution}  

// 接下来，我们定义系数矩阵 $\mathbf{K}$ （这里是身份矩阵），迪里希特边界条件，右手边 $f = 2\pi^2 \sin(\pi x) \sin(\pi y)$  ，以及与这些选择相对应的 $K$ 和 $f$ 的精确解，即 $p = \sin(\pi x) \sin(\pi y)$  。

  template <int dim> 
  class Coefficient : public TensorFunction<2, dim> 
  { 
  public: 
    Coefficient() 
      : TensorFunction<2, dim>() 
    {} 

    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<Tensor<2, dim>> &values) const override; 
  }; 

  template <int dim> 
  void Coefficient<dim>::value_list(const std::vector<Point<dim>> &points, 
                                    std::vector<Tensor<2, dim>> &  values) const 
  { 
    Assert(points.size() == values.size(), 
           ExcDimensionMismatch(points.size(), values.size())); 
    for (unsigned int p = 0; p < points.size(); ++p) 
      values[p] = unit_symmetric_tensor<dim>(); 
  } 

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    BoundaryValues() 
      : Function<dim>(2) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> & /*p*/, 
                                    const unsigned int /*component*/) const 
  { 
    return 0; 
  } 

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  };

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> &p, 
                                   const unsigned int /*component*/) const 
  { 
    return (2 * numbers::PI * numbers::PI * std::sin(numbers::PI * p[0]) * 
            std::sin(numbers::PI * p[1])); 
  } 

// 实现精确压力解决方案的类有一个奇怪的地方，我们把它作为一个有两个分量的向量值来实现。(我们在构造函数中说它有两个分量，在这里我们调用基函数类的构造函数)。在`value()`函数中，我们不测试`component`参数的值，这意味着我们为向量值函数的两个分量返回相同的值。我们这样做是因为我们将本程序中使用的有限元描述为一个包含内部和界面压力的矢量值系统，当我们计算误差时，我们希望使用相同的压力解来测试这两个分量。

  template <int dim> 
  class ExactPressure : public Function<dim> 
  { 
  public: 
    ExactPressure() 
      : Function<dim>(2) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component) const override; 
  }; 

  template <int dim> 
  double ExactPressure<dim>::value(const Point<dim> &p, 
                                   const unsigned int /*component*/) const 
  { 
    return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]); 
  } 

  template <int dim> 
  class ExactVelocity : public TensorFunction<1, dim> 
  { 
  public: 
    ExactVelocity() 
      : TensorFunction<1, dim>() 
    {} 

    virtual Tensor<1, dim> value(const Point<dim> &p) const override; 
  }; 

  template <int dim> 
  Tensor<1, dim> ExactVelocity<dim>::value(const Point<dim> &p) const 
  { 
    Tensor<1, dim> return_value; 
    return_value[0] = -numbers::PI * std::cos(numbers::PI * p[0]) * 
                      std::sin(numbers::PI * p[1]); 
    return_value[1] = -numbers::PI * std::sin(numbers::PI * p[0]) * 
                      std::cos(numbers::PI * p[1]); 
    return return_value; 
  } 

//  @sect3{WGDarcyEquation class implementation}  
// @sect4{WGDarcyEquation::WGDarcyEquation}  

// 在这个构造函数中，我们创建了一个矢量值函数的有限元空间，这里将包括用于内部和界面压力的函数， $p^\circ$  和  $p^\partial$  。

  template <int dim> 
  WGDarcyEquation<dim>::WGDarcyEquation(const unsigned int degree) 
    : fe(FE_DGQ<dim>(degree), 1, FE_FaceQ<dim>(degree), 1) 
    , dof_handler(triangulation) 
    , fe_dgrt(degree) 
    , dof_handler_dgrt(triangulation) 
  {} 

//  @sect4{WGDarcyEquation::make_grid}  

// 我们在单位平方域上生成一个网格并对其进行细化。

  template <int dim> 
  void WGDarcyEquation<dim>::make_grid() 
  { 
    GridGenerator::hyper_cube(triangulation, 0, 1); 
    triangulation.refine_global(5); 

    std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "   Total number of cells: " << triangulation.n_cells() 
              << std::endl; 
  } 

//  @sect4{WGDarcyEquation::setup_system}  

// 在我们创建了上面的网格后，我们分配自由度并调整矩阵和向量的大小。这个函数中唯一值得关注的部分是我们如何插值压力的边界值。由于压力由内部和界面分量组成，我们需要确保我们只插值到矢量值解空间中与界面压力相对应的分量上（因为这些分量是唯一定义在域的边界上的）。我们通过一个只针对界面压力的分量屏蔽对象来做到这一点。

  template <int dim> 
  void WGDarcyEquation<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    dof_handler_dgrt.distribute_dofs(fe_dgrt); 

    std::cout << "   Number of pressure degrees of freedom: " 
              << dof_handler.n_dofs() << std::endl; 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    { 
      constraints.clear(); 
      const FEValuesExtractors::Scalar interface_pressure(1); 
      const ComponentMask              interface_pressure_mask = 
        fe.component_mask(interface_pressure); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               0, 
                                               BoundaryValues<dim>(), 
                                               constraints, 
                                               interface_pressure_mask); 
      constraints.close(); 
    } 

// 在双线性形式中，在两个相邻单元之间的面上没有积分项，所以我们可以直接使用 <code>DoFTools::make_sparsity_pattern</code> 来计算稀疏矩阵。

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
  } 

//  @sect4{WGDarcyEquation::assemble_system}  

// 这个函数比较有趣。正如介绍中所详述的，线性系统的装配要求我们评估形状函数的弱梯度，这是Raviart-Thomas空间的一个元素。因此，我们需要定义一个Raviart-Thomas有限元对象，并有FEValues对象在正交点评估它。然后我们需要计算每个单元 $K$ 上的矩阵 $C^K$ ，为此我们需要介绍中提到的矩阵 $M^K$ 和 $G^K$ 。

// 有一点可能不是很明显，在之前所有的教程程序中，我们总是用DoFHandler的单元格迭代器来调用 FEValues::reinit() 。这样就可以调用诸如 FEValuesBase::get_function_values() 这样的函数，在单元格的正交点上提取有限元函数的值（用DoF值的矢量表示）。为了使这种操作发挥作用，人们需要知道哪些向量元素对应于给定单元上的自由度--也就是说，正是DoFHandler类所提供的那种信息和操作。

// 我们可以为 "破碎的 "Raviart-Thomas空间创建一个DoFHandler对象（使用FE_DGRT类），但是我们在这里真的不想这样做。至少在当前函数中，我们不需要任何与这个破碎空间相关的全局定义的自由度，而只需要引用当前单元上的这种空间的形状函数。因此，我们利用这样一个事实，即人们也可以用单元格迭代器来调用 FEValues::reinit() 的Triangulation对象（而不是DoFHandler对象）。在这种情况下，FEValues当然只能为我们提供只引用单元格的信息，而不是这些单元格上列举的自由度。所以我们不能使用 FEValuesBase::get_function_values(), ，但我们可以使用 FEValues::shape_value() 来获取当前单元上正交点的形状函数值。下面我们要利用的就是这种功能。下面给我们提供Raviart-Thomas函数信息的变量是`fe_values_rt`（和相应的`fe_face_values_rt`）对象。

// 鉴于上述介绍，下面的声明应该是非常明显的。

  template <int dim> 
  void WGDarcyEquation<dim>::assemble_system() 
  { 
    const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 

    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values); 
    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 

    FEValues<dim>     fe_values_dgrt(fe_dgrt, 
                                 quadrature_formula, 
                                 update_values | update_gradients | 
                                   update_quadrature_points | 
                                   update_JxW_values); 
    FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
                                          face_quadrature_formula, 
                                          update_values | 
                                            update_normal_vectors | 
                                            update_quadrature_points | 
                                            update_JxW_values); 

    const unsigned int dofs_per_cell      = fe.n_dofs_per_cell(); 
    const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell(); 

    const unsigned int n_q_points      = fe_values.get_quadrature().size(); 
    const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 

    const unsigned int n_face_q_points = fe_face_values.get_quadrature().size(); 

    RightHandSide<dim>  right_hand_side; 
    std::vector<double> right_hand_side_values(n_q_points); 

    const Coefficient<dim>      coefficient; 
    std::vector<Tensor<2, dim>> coefficient_values(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 接下来，让我们声明介绍中讨论的各种单元格矩阵。

    FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell); 
    FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt); 
    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 
    Vector<double>     cell_solution(dofs_per_cell); 

// 我们需要  <code>FEValuesExtractors</code>  来访问形状函数的  @p interior  和  @p face  部分。

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure_interior(0); 
    const FEValuesExtractors::Scalar pressure_face(1); 

// 这最终让我们在所有单元格上进行循环。在每个单元中，我们将首先计算用于构建局部矩阵的各种单元矩阵--因为它们取决于相关的单元，所以它们需要在每个单元中重新计算。我们还需要Raviart-Thomas空间的形状函数，为此我们需要首先创建一个通往三角化单元的迭代器，我们可以通过从指向DoFHandler的单元中的赋值来获得。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 

        const typename Triangulation<dim>::active_cell_iterator cell_dgrt = 
          cell; 
        fe_values_dgrt.reinit(cell_dgrt); 

        right_hand_side.value_list(fe_values.get_quadrature_points(), 
                                   right_hand_side_values); 
        coefficient.value_list(fe_values.get_quadrature_points(), 
                               coefficient_values); 

// 我们要计算的第一个单元矩阵是拉维-托马斯空间的质量矩阵。 因此，我们需要循环计算速度FEValues对象的所有正交点。

        cell_matrix_M = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q); 
              for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
                { 
                  const Tensor<1, dim> v_k = 
                    fe_values_dgrt[velocities].value(k, q); 
                  cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q)); 
                } 
            } 

// 接下来我们通过使用 FullMatrix::gauss_jordan(). 对这个矩阵进行求逆 它将被用来计算后面的系数矩阵 $C^K$ 。值得一提的是，后面的 "cell_matrix_M "实际上包含了*的逆*。
//在这个调用之后的 $M^K$ 的*逆。

        cell_matrix_M.gauss_jordan(); 

// 从介绍中，我们知道定义 $C^K$ 的方程的右边 $G^K$ 是面积分和单元积分的区别。在这里，我们对内部的贡献的负值进行了近似。这个矩阵的每个分量都是多项式空间的一个基函数与拉维-托马斯空间的一个基函数的发散之间的乘积的积分。这些基函数是在内部定义的。

        cell_matrix_G = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const double div_v_i = 
                fe_values_dgrt[velocities].divergence(i, q); 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const double phi_j_interior = 
                    fe_values[pressure_interior].value(j, q); 

                  cell_matrix_G(i, j) -= 
                    (div_v_i * phi_j_interior * fe_values.JxW(q)); 
                } 
            } 

// 接下来，我们用正交法对面的积分进行近似。每个分量都是多项式空间的基函数与Raviart-Thomas空间的基函数与法向量的点积的积分。所以我们在元素的所有面上循环，得到法向量。

        for (const auto &face : cell->face_iterators()) 
          { 
            fe_face_values.reinit(cell, face); 
            fe_face_values_dgrt.reinit(cell_dgrt, face); 

            for (unsigned int q = 0; q < n_face_q_points; ++q) 
              { 
                const Tensor<1, dim> &normal = fe_face_values.normal_vector(q); 

                for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
                  { 
                    const Tensor<1, dim> v_i = 
                      fe_face_values_dgrt[velocities].value(i, q); 
                    for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                      { 
                        const double phi_j_face = 
                          fe_face_values[pressure_face].value(j, q); 

                        cell_matrix_G(i, j) += 
                          ((v_i * normal) * phi_j_face * fe_face_values.JxW(q)); 
                      } 
                  } 
              } 
          } 
// @p cell_matrix_C 是 $G^K$ 的转置与质量矩阵的逆之间的矩阵乘积（该逆存储在 @p cell_matrix_M): 中）。
        cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M); 

// 最后我们可以计算出本地矩阵  $A^K$  。 元素  $A^K_{ij}$  由  $\int_{E} \sum_{k,l} C_{ik} C_{jl} (\mathbf{K} \mathbf{v}_k) \cdot \mathbf{v}_l \mathrm{d}x$  得到。我们在上一步已经计算了系数 $C$ ，因此在适当地重新排列循环后得到以下结果。

        local_matrix = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
              { 
                const Tensor<1, dim> v_k = 
                  fe_values_dgrt[velocities].value(k, q); 
                for (unsigned int l = 0; l < dofs_per_cell_dgrt; ++l) 
                  { 
                    const Tensor<1, dim> v_l = 
                      fe_values_dgrt[velocities].value(l, q); 

                    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                      for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                        local_matrix(i, j) += 
                          (coefficient_values[q] * cell_matrix_C[i][k] * v_k) * 
                          cell_matrix_C[j][l] * v_l * fe_values_dgrt.JxW(q); 
                  } 
              } 
          } 

// 接下来，我们计算右手边， $\int_{K} f q \mathrm{d}x$  。

        cell_rhs = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              cell_rhs(i) += (fe_values[pressure_interior].value(i, q) * 
                              right_hand_side_values[q] * fe_values.JxW(q)); 
            } 

// 最后一步是将本地矩阵的组件分配到系统矩阵中，并将单元格右侧的组件转移到系统右侧。

        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global( 
          local_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 

//  @sect4{WGDarcyEquation<dim>::solve}  

// 这一步相当琐碎，与之前的许多教程程序相同。

  template <int dim> 
  void WGDarcyEquation<dim>::solve() 
  { 
    SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> solver(solver_control); 
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 
    constraints.distribute(solution); 
  } 
// @sect4{WGDarcyEquation<dim>::compute_postprocessed_velocity}  

// 在这个函数中，根据之前计算的压力解计算出速度场。速度被定义为 $\mathbf{u}_h = \mathbf{Q}_h \left(-\mathbf{K}\nabla_{w,d}p_h \right)$ ，这需要我们计算许多与系统矩阵组装相同的项。还有一些矩阵 $E^K,D^K$ 我们也需要组装（见介绍），但它们实际上只是遵循相同的模式。

// 在这里计算与我们在`assemble_system()`函数中已经完成的相同的矩阵，当然是浪费CPU时间的。同样地，我们把那里的一些代码复制到这个函数中，这通常也是一个糟糕的主意。一个更好的实现可能会提供一个函数来封装这些重复的代码。我们也可以考虑使用计算效率和内存效率之间的经典权衡，在装配过程中每个单元只计算一次 $C^K$ 矩阵，把它们存储在边上的某个地方，然后在这里重新使用它们。例如， step-51 就是这样做的，`assemble_system()`函数需要一个参数来决定是否重新计算本地矩阵，类似的方法--也许是将本地矩阵存储在其他地方--可以适用于当前的程序）。

  template <int dim> 
  void WGDarcyEquation<dim>::compute_postprocessed_velocity() 
  { 
    darcy_velocity.reinit(dof_handler_dgrt.n_dofs()); 

    const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values); 

    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 

    FEValues<dim> fe_values_dgrt(fe_dgrt, 
                                 quadrature_formula, 
                                 update_values | update_gradients | 
                                   update_quadrature_points | 
                                   update_JxW_values); 

    FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
                                          face_quadrature_formula, 
                                          update_values | 
                                            update_normal_vectors | 
                                            update_quadrature_points | 
                                            update_JxW_values); 

    const unsigned int dofs_per_cell      = fe.n_dofs_per_cell(); 
    const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell(); 

    const unsigned int n_q_points      = fe_values.get_quadrature().size(); 
    const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 

    const unsigned int n_face_q_points = fe_face_values.get_quadrature().size(); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices_dgrt( 
      dofs_per_cell_dgrt); 

    FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell); 
    FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_D(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_E(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 

    Vector<double> cell_solution(dofs_per_cell); 
    Vector<double> cell_velocity(dofs_per_cell_dgrt); 

    const Coefficient<dim>      coefficient; 
    std::vector<Tensor<2, dim>> coefficient_values(n_q_points_dgrt); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure_interior(0); 
    const FEValuesExtractors::Scalar pressure_face(1); 

// 在介绍中，我们解释了如何计算单元上的数值速度。我们需要每个单元上的压力解值、格拉姆矩阵的系数和 $L_2$ 投影的系数。我们已经计算了全局解，所以我们将从全局解中提取单元解。格拉姆矩阵的系数在我们计算压力的系统矩阵时已经计算过了。我们在这里也要这样做。对于投影的系数，我们做矩阵乘法，即用格拉姆矩阵的倒数乘以 $(\mathbf{K} \mathbf{w}, \mathbf{w})$ 的矩阵作为组成部分。然后，我们将所有这些系数相乘，称之为β。数值速度是贝塔和拉维尔特-托马斯空间的基础函数的乘积。

    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(), 
      endc = dof_handler.end(), cell_dgrt = dof_handler_dgrt.begin_active(); 
    for (; cell != endc; ++cell, ++cell_dgrt) 
      { 
        fe_values.reinit(cell); 
        fe_values_dgrt.reinit(cell_dgrt); 

        coefficient.value_list(fe_values_dgrt.get_quadrature_points(), 
                               coefficient_values); 

// 这个 <code>cell_matrix_E</code> 的分量是 $(\mathbf{K} \mathbf{w}, \mathbf{w})$ 的积分。  <code>cell_matrix_M</code> 是格拉姆矩阵。

        cell_matrix_M = 0; 
        cell_matrix_E = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q); 
              for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
                { 
                  const Tensor<1, dim> v_k = 
                    fe_values_dgrt[velocities].value(k, q); 

                  cell_matrix_E(i, k) += 
                    (coefficient_values[q] * v_i * v_k * fe_values_dgrt.JxW(q)); 

                  cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q)); 
                } 
            } 

// 为了计算介绍中提到的矩阵 $D$ ，我们就需要按照介绍中的解释来评估 $D=M^{-1}E$ 。

        cell_matrix_M.gauss_jordan(); 
        cell_matrix_M.mmult(cell_matrix_D, cell_matrix_E); 

// 然后，我们还需要再次计算矩阵 $C$ ，用于评估弱离散梯度。这与组装系统矩阵时使用的代码完全相同，所以我们只需从那里复制它。

        cell_matrix_G = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const double div_v_i = 
                fe_values_dgrt[velocities].divergence(i, q); 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const double phi_j_interior = 
                    fe_values[pressure_interior].value(j, q); 

                  cell_matrix_G(i, j) -= 
                    (div_v_i * phi_j_interior * fe_values.JxW(q)); 
                } 
            } 

        for (const auto &face : cell->face_iterators()) 
          { 
            fe_face_values.reinit(cell, face); 
            fe_face_values_dgrt.reinit(cell_dgrt, face); 

            for (unsigned int q = 0; q < n_face_q_points; ++q) 
              { 
                const Tensor<1, dim> &normal = fe_face_values.normal_vector(q); 

                for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
                  { 
                    const Tensor<1, dim> v_i = 
                      fe_face_values_dgrt[velocities].value(i, q); 
                    for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                      { 
                        const double phi_j_face = 
                          fe_face_values[pressure_face].value(j, q); 

                        cell_matrix_G(i, j) += 
                          ((v_i * normal) * phi_j_face * fe_face_values.JxW(q)); 
                      } 
                  } 
              } 
          } 
        cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M); 

// 最后，我们需要提取对应于当前单元的压力未知数。

        cell->get_dof_values(solution, cell_solution); 

// 我们现在可以计算当地的速度未知数（相对于我们将 $-\mathbf K \nabla_{w,d} p_h$ 项投影到的Raviart-Thomas空间而言）。

        cell_velocity = 0; 
        for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
          for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j) 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              cell_velocity(k) += 
                -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j)); 

// 我们计算达西速度。这与cell_velocity相同，但用于绘制Darcy速度图。

        cell_dgrt->get_dof_indices(local_dof_indices_dgrt); 
        for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
          for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j) 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              darcy_velocity(local_dof_indices_dgrt[k]) += 
                -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j)); 
      } 
  } 

//  @sect4{WGDarcyEquation<dim>::compute_pressure_error}  

// 这一部分是为了计算压力的 $L_2$ 误差。 我们定义一个向量，用来保存每个单元上的误差规范。接下来，我们使用 VectorTool::integrate_difference() 来计算每个单元上的 $L_2$ 准则的误差。然而，我们实际上只关心解向量的内部分量的误差（我们甚至不能评估正交点的界面压力，因为这些都位于单元格的内部），因此必须使用一个权重函数，确保解变量的界面分量被忽略。这是通过使用ComponentSelectFunction来实现的，其参数表明我们要选择哪个分量（零分量，即内部压力）以及总共有多少分量（两个）。

  template <int dim> 
  void WGDarcyEquation<dim>::compute_pressure_error() 
  { 
    Vector<float> difference_per_cell(triangulation.n_active_cells()); 
    const ComponentSelectFunction<dim> select_interior_pressure(0, 2); 
    VectorTools::integrate_difference(dof_handler, 
                                      solution, 
                                      ExactPressure<dim>(), 
                                      difference_per_cell, 
                                      QGauss<dim>(fe.degree + 2), 
                                      VectorTools::L2_norm, 
                                      &select_interior_pressure); 

    const double L2_error = difference_per_cell.l2_norm(); 
    std::cout << "L2_error_pressure " << L2_error << std::endl; 
  } 

//  @sect4{WGDarcyEquation<dim>::compute_velocity_error}  

// 在这个函数中，我们评估每个单元的速度的 $L_2$ 误差，以及面的流量的 $L_2$ 误差。该函数依赖于之前计算过的`compute_postprocessed_velocity()`函数，该函数根据之前计算过的压力解来计算速度场。

// 我们将评估每个单元的速度，并计算数值速度和精确速度之间的差异。

  template <int dim> 
  void WGDarcyEquation<dim>::compute_velocity_errors() 
  { 
    const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 

    FEValues<dim> fe_values_dgrt(fe_dgrt, 
                                 quadrature_formula, 
                                 update_values | update_gradients | 
                                   update_quadrature_points | 
                                   update_JxW_values); 

    FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
                                          face_quadrature_formula, 
                                          update_values | 
                                            update_normal_vectors | 
                                            update_quadrature_points | 
                                            update_JxW_values); 

    const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 
    const unsigned int n_face_q_points_dgrt = 
      fe_face_values_dgrt.get_quadrature().size(); 

    std::vector<Tensor<1, dim>> velocity_values(n_q_points_dgrt); 
    std::vector<Tensor<1, dim>> velocity_face_values(n_face_q_points_dgrt); 

    const FEValuesExtractors::Vector velocities(0); 

    const ExactVelocity<dim> exact_velocity; 

    double L2_err_velocity_cell_sqr_global = 0; 
    double L2_err_flux_sqr                 = 0; 

// 在之前计算了后处理的速度之后，我们在这里只需要提取每个单元和面的相应数值，并与精确的数值进行比较。

    for (const auto &cell_dgrt : dof_handler_dgrt.active_cell_iterators()) 
      { 
        fe_values_dgrt.reinit(cell_dgrt); 

// 首先计算后处理的速度场与精确速度场之间的 $L_2$ 误差。

        fe_values_dgrt[velocities].get_function_values(darcy_velocity, 
                                                       velocity_values); 
        double L2_err_velocity_cell_sqr_local = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          { 
            const Tensor<1, dim> velocity = velocity_values[q]; 
            const Tensor<1, dim> true_velocity = 
              exact_velocity.value(fe_values_dgrt.quadrature_point(q)); 

            L2_err_velocity_cell_sqr_local += 
              ((velocity - true_velocity) * (velocity - true_velocity) * 
               fe_values_dgrt.JxW(q)); 
          } 
        L2_err_velocity_cell_sqr_global += L2_err_velocity_cell_sqr_local; 

// 为了重建通量，我们需要单元格和面的大小。由于通量是按面计算的，我们必须在每个单元的所有四个面上进行循环。为了计算面的速度，我们从之前计算的`darcy_velocity`中提取正交点的值。然后，我们计算法线方向的速度平方误差。最后，我们通过对面和单元面积的适当缩放来计算单元上的 $L_2$ 通量误差，并将其加入全局误差。

        const double cell_area = cell_dgrt->measure(); 
        for (const auto &face_dgrt : cell_dgrt->face_iterators()) 
          { 
            const double face_length = face_dgrt->measure(); 
            fe_face_values_dgrt.reinit(cell_dgrt, face_dgrt); 
            fe_face_values_dgrt[velocities].get_function_values( 
              darcy_velocity, velocity_face_values); 

            double L2_err_flux_face_sqr_local = 0; 
            for (unsigned int q = 0; q < n_face_q_points_dgrt; ++q) 
              { 
                const Tensor<1, dim> velocity = velocity_face_values[q]; 
                const Tensor<1, dim> true_velocity = 
                  exact_velocity.value(fe_face_values_dgrt.quadrature_point(q)); 

                const Tensor<1, dim> &normal = 
                  fe_face_values_dgrt.normal_vector(q); 

 
                  ((velocity * normal - true_velocity * normal) * 
                   (velocity * normal - true_velocity * normal) * 
                   fe_face_values_dgrt.JxW(q)); 
              } 
            const double err_flux_each_face = 
              L2_err_flux_face_sqr_local / face_length * cell_area; 
            L2_err_flux_sqr += err_flux_each_face; 
          } 
      } 

// 将所有单元和面的误差相加后，我们进行平方根计算，得到速度和流量的 $L_2$ 误差。我们将这些数据输出到屏幕上。

    const double L2_err_velocity_cell = 
      std::sqrt(L2_err_velocity_cell_sqr_global); 
    const double L2_err_flux_face = std::sqrt(L2_err_flux_sqr); 

    std::cout << "L2_error_vel:  " << L2_err_velocity_cell << std::endl 
              << "L2_error_flux: " << L2_err_flux_face << std::endl; 
  } 
// @sect4{WGDarcyEquation::output_results}  

// 我们有两组结果要输出：内部解和骨架解。我们使用 <code>DataOut</code> 来显示内部结果。骨架结果的图形输出是通过使用DataOutFaces类完成的。

// 在这两个输出文件中，内部和面的变量都被存储。对于界面输出，输出文件只是包含了内部压力对面的插值，但是因为没有确定从两个相邻的单元中得到的是哪一个内部压力变量，所以在界面输出文件中最好是忽略内部压力。相反，对于单元格内部输出文件，当然不可能显示任何界面压力 $p^\partial$ ，因为这些压力只适用于界面，而不是单元格内部。因此，你会看到它们被显示为一个无效的值（比如一个无穷大）。

// 对于单元内部的输出，我们还想输出速度变量。这有点棘手，因为它生活在同一个网格上，但使用不同的DoFHandler对象（压力变量生活在`dof_handler`对象上，达西速度生活在`dof_handler_dgrt`对象上）。幸运的是， DataOut::add_data_vector() 函数有一些变化，允许指定一个矢量对应的DoFHandler，因此我们可以在同一个文件中对两个DoFHandler对象的数据进行可视化。

  template <int dim> 
  void WGDarcyEquation<dim>::output_results() const 
  { 
    { 
      DataOut<dim> data_out; 

// 首先将压力解决方案附加到DataOut对象上。

      const std::vector<std::string> solution_names = {"interior_pressure", 
                                                       "interface_pressure"}; 
      data_out.add_data_vector(dof_handler, solution, solution_names); 

// 然后对达西速度场做同样的处理，并继续将所有内容写进文件。

      const std::vector<std::string> velocity_names(dim, "velocity"); 
      const std::vector< 
        DataComponentInterpretation::DataComponentInterpretation> 
        velocity_component_interpretation( 
          dim, DataComponentInterpretation::component_is_part_of_vector); 
      data_out.add_data_vector(dof_handler_dgrt, 
                               darcy_velocity, 
                               velocity_names, 
                               velocity_component_interpretation); 

      data_out.build_patches(fe.degree); 
      std::ofstream output("solution_interior.vtu"); 
      data_out.write_vtu(output); 
    } 

    { 
      DataOutFaces<dim> data_out_faces(false); 
      data_out_faces.attach_dof_handler(dof_handler); 
      data_out_faces.add_data_vector(solution, "Pressure_Face"); 
      data_out_faces.build_patches(fe.degree); 
      std::ofstream face_output("solution_interface.vtu"); 
      data_out_faces.write_vtu(face_output); 
    } 
  } 
// @sect4{WGDarcyEquation::run}  

// 这是主类的最后一个函数。它调用我们类的其他函数。

  template <int dim> 
  void WGDarcyEquation<dim>::run() 
  { 
    std::cout << "Solving problem in " << dim << " space dimensions." 
              << std::endl; 
    make_grid(); 
    setup_system(); 
    assemble_system(); 
    solve(); 
    compute_postprocessed_velocity(); 
    compute_pressure_error(); 
    compute_velocity_errors(); 
    output_results(); 
  } 

} // namespace Step61 
// @sect3{The <code>main</code> function}  

// 这是主函数。我们可以在这里改变维度以在3D中运行。

int main() 
{ 
  try 
    { 
      Step61::WGDarcyEquation<2> wg_darcy(0); 
      wg_darcy.run(); 
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

