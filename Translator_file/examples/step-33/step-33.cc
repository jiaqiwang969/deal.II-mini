CCTest_file/step-33.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2007 - 2021 by the deal.II authors 
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
 * Author: David Neckels, Boulder, Colorado, 2007, 2008 
 */ 


// @sect3{Include files}  

// 首先是一套标准的deal.II包括。这里没有什么特别需要评论的。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/function_parser.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/conditional_ostream.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_in.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/mapping_q1.h> 
#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/solution_transfer.h> 

// 然后，正如介绍中提到的，我们使用各种Trilinos软件包作为线性求解器以及自动微分。这些都在以下的包含文件中。

// 由于deal.II提供了基本的Trilinos矩阵、预处理程序和求解器的接口，我们把它们作为deal.II线性代数结构类似地包括在内。

#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_precondition.h> 
#include <deal.II/lac/trilinos_solver.h> 

// Sacado是Trilinos中的自动微分包，用于寻找全隐式牛顿迭代的雅各布系数。

#include <Sacado.hpp> 

// 这又是C++语言。

#include <iostream> 
#include <fstream> 
#include <vector> 
#include <memory> 
#include <array> 

// 在本节结束时，将dealii库中的所有内容引入本程序内容将进入的命名空间。

namespace Step33 
{ 
  using namespace dealii; 
// @sect3{Euler equation specifics}  

// 这里我们定义了这个特定的守恒定律系统的通量函数，以及几乎所有其他的气体动力学欧拉方程所特有的东西，原因在介绍中讨论过。我们将所有这些归入一个结构，该结构定义了所有与通量有关的东西。这个结构的所有成员都是静态的，也就是说，这个结构没有由实例成员变量指定的实际状态。更好的方法是使用命名空间，而不是一个拥有所有静态成员的结构--但是命名空间不能被模板化，而且我们希望结构中的一些成员变量取决于空间维度，我们以通常的方式用模板参数来引入。

  template <int dim> 
  struct EulerEquations 
  { 
// @sect4{Component description}  

// 首先是几个变量，它们以一种通用的方式描述了我们的解向量的各个组成部分。这包括系统中分量的数量（欧拉方程中每个空间方向的动量都有一个条目，加上能量和密度分量，总共有 <code>dim+2</code> 个分量），以及描述第一个动量分量、密度分量和能量密度分量在解向量中的索引的函数。请注意，所有这些%数都取决于空间维度；以通用的方式定义它们（而不是以隐含的惯例）使我们的代码更加灵活，并使以后的扩展更加容易，例如，在方程中加入更多的分量。

    static const unsigned int n_components             = dim + 2; 
    static const unsigned int first_momentum_component = 0; 
    static const unsigned int density_component        = dim; 
    static const unsigned int energy_component         = dim + 1; 

// 在这个程序中一路生成图形输出时，我们需要指定解变量的名称，以及各种成分如何分组为矢量和标量场。我们可以在这里进行描述，但是为了使与欧拉方程有关的事情在这里得到解决，并使程序的其他部分尽可能地通用，我们在以下两个函数中提供了这类信息。

    static std::vector<std::string> component_names() 
    { 
      std::vector<std::string> names(dim, "momentum"); 
      names.emplace_back("density"); 
      names.emplace_back("energy_density"); 

      return names; 
    } 

    static std::vector<DataComponentInterpretation::DataComponentInterpretation> 
    component_interpretation() 
    { 
      std::vector<DataComponentInterpretation::DataComponentInterpretation> 
        data_component_interpretation( 
          dim, DataComponentInterpretation::component_is_part_of_vector); 
      data_component_interpretation.push_back( 
        DataComponentInterpretation::component_is_scalar); 
      data_component_interpretation.push_back( 
        DataComponentInterpretation::component_is_scalar); 

      return data_component_interpretation; 
    } 
// @sect4{Transformations between variables}  

// 接下来，我们定义气体常数。我们将在紧接着这个类的声明之后的定义中把它设置为1.4（与整数变量不同，比如上面的变量，静态常量浮点成员变量不能在C++的类声明中被初始化）。这个1.4的值代表了由两个原子组成的分子的气体，比如空气，它几乎完全由 $N_2$ 和 $O_2$ 组成，痕迹很小。

    static const double gas_gamma; 

// 在下文中，我们将需要从保守变量的矢量中计算动能和压力。我们可以根据能量密度和动能 $\frac 12 \rho |\mathbf v|^2= \frac{|\rho \mathbf v|^2}{2\rho}$ 来做这件事（注意，独立变量包含动量分量 $\rho v_i$ ，而不是速度 $v_i$ ）。

    template <typename InputVector> 
    static typename InputVector::value_type 
    compute_kinetic_energy(const InputVector &W) 
    { 
      typename InputVector::value_type kinetic_energy = 0; 
      for (unsigned int d = 0; d < dim; ++d) 
        kinetic_energy += 
          W[first_momentum_component + d] * W[first_momentum_component + d]; 
      kinetic_energy *= 1. / (2 * W[density_component]); 

      return kinetic_energy; 
    } 

    template <typename InputVector> 
    static typename InputVector::value_type 
    compute_pressure(const InputVector &W) 
    { 
      return ((gas_gamma - 1.0) * 
              (W[energy_component] - compute_kinetic_energy(W))); 
    } 
// @sect4{EulerEquations::compute_flux_matrix}  

// 我们将通量函数 $F(W)$ 定义为一个大矩阵。 这个矩阵的每一行都代表了该行成分的标量守恒定律。 这个矩阵的确切形式在介绍中给出。请注意，我们知道这个矩阵的大小：它的行数与系统的分量一样多， <code>dim</code> 列数一样多；我们没有为这样的矩阵使用FullMatrix对象（它的行数和列数是可变的，因此每次创建这样的矩阵时必须在堆上分配内存），而是马上使用一个矩形的数字阵列。

// 我们将通量函数的数值类型模板化，这样我们就可以在这里使用自动微分类型。 同样地，我们将用不同的输入矢量数据类型来调用该函数，所以我们也对其进行模板化。

    template <typename InputVector> 
    static void compute_flux_matrix(const InputVector &W, 
                                    ndarray<typename InputVector::value_type, 
                                            EulerEquations<dim>::n_components, 
                                            dim> &     flux) 
    { 

// 首先计算出现在通量矩阵中的压力，然后计算矩阵中对应于动量项的前 <code>dim</code> 列。

      const typename InputVector::value_type pressure = compute_pressure(W); 

      for (unsigned int d = 0; d < dim; ++d) 
        { 
          for (unsigned int e = 0; e < dim; ++e) 
            flux[first_momentum_component + d][e] = 
              W[first_momentum_component + d] * 
              W[first_momentum_component + e] / W[density_component]; 

          flux[first_momentum_component + d][d] += pressure; 
        } 

// 然后是密度（即质量守恒）的条款，最后是能量守恒。

      for (unsigned int d = 0; d < dim; ++d) 
        flux[density_component][d] = W[first_momentum_component + d]; 

      for (unsigned int d = 0; d < dim; ++d) 
        flux[energy_component][d] = W[first_momentum_component + d] / 
                                    W[density_component] * 
                                    (W[energy_component] + pressure); 
    } 
// @sect4{EulerEquations::compute_normal_flux}  

// 在域的边界和跨挂节点上，我们使用一个数值通量函数来强制执行边界条件。 这个程序是基本的Lax-Friedrich的通量，有一个稳定的参数  $\alpha$  。它的形式也已经在介绍中给出。

    template <typename InputVector> 
    static void numerical_normal_flux( 
      const Tensor<1, dim> &                                      normal, 
      const InputVector &                                         Wplus, 
      const InputVector &                                         Wminus, 
      const double                                                alpha, 
      std::array<typename InputVector::value_type, n_components> &normal_flux) 
    { 
      ndarray<typename InputVector::value_type, 
              EulerEquations<dim>::n_components, 
              dim> 
        iflux, oflux; 

      compute_flux_matrix(Wplus, iflux); 
      compute_flux_matrix(Wminus, oflux); 

      for (unsigned int di = 0; di < n_components; ++di) 
        { 
          normal_flux[di] = 0; 
          for (unsigned int d = 0; d < dim; ++d) 
            normal_flux[di] += 0.5 * (iflux[di][d] + oflux[di][d]) * normal[d]; 

          normal_flux[di] += 0.5 * alpha * (Wplus[di] - Wminus[di]); 
        } 
    } 
// @sect4{EulerEquations::compute_forcing_vector}  

// 与描述通量函数 $\mathbf F(\mathbf w)$ 的方式相同，我们也需要有一种方法来描述右侧的强迫项。正如介绍中提到的，我们在这里只考虑重力，这导致了具体的形式 $\mathbf G(\mathbf w) = \left( g_1\rho, g_2\rho, g_3\rho, 0, \rho \mathbf g \cdot \mathbf v \right)^T$ ，这里显示的是三维情况。更具体地说，我们将只考虑三维的 $\mathbf g=(0,0,-1)^T$ ，或二维的 $\mathbf g=(0,-1)^T$ 。这自然导致了以下函数。

    template <typename InputVector> 
    static void compute_forcing_vector( 
      const InputVector &                                         W, 
      std::array<typename InputVector::value_type, n_components> &forcing) 
    { 
      const double gravity = -1.0; 

      for (unsigned int c = 0; c < n_components; ++c) 
        switch (c) 
          { 
            case first_momentum_component + dim - 1: 
              forcing[c] = gravity * W[density_component]; 
              break; 
            case energy_component: 
              forcing[c] = gravity * W[first_momentum_component + dim - 1]; 
              break; 
            default: 
              forcing[c] = 0; 
          } 
    } 
// @sect4{Dealing with boundary conditions}  

// 我们必须处理的另一件事是边界条件。为此，让我们首先定义一下我们目前知道如何处理的各种边界条件。

    enum BoundaryKind 
    { 
      inflow_boundary, 
      outflow_boundary, 
      no_penetration_boundary, 
      pressure_boundary 
    }; 

// 接下来的部分是实际决定在每一种边界上做什么。为此，请记住，从介绍中可以看出，边界条件是通过在给定的不均匀性 $\mathbf j$ 的边界外侧选择一个值 $\mathbf w^-$ ，以及可能在内部选择解的值 $\mathbf w^+$ 来指定的。然后，两者都被传递给数值通量 $\mathbf H(\mathbf{w}^+, \mathbf{w}^-, \mathbf{n})$ ，以定义边界对双线性形式的贡献。

// 边界条件在某些情况下可以为解矢量的每个分量独立指定。例如，如果分量 $c$ 被标记为流入，那么 $w^-_c = j_c$  。如果是流出，那么 $w^-_c = w^+_c$  。这两种简单的情况在下面的函数中首先得到处理。

// 有一个小插曲，从C++语言的角度来看，这个函数是不愉快的。输出向量  <code>Wminus</code>  当然会被修改，所以它不应该是  <code>const</code>  的参数。然而，在下面的实现中，它却成为了参数，而且为了使代码能够编译，它必须成为参数。原因是我们在 <code>Wminus</code> 类型为 <code>Table@<2,Sacado::Fad::DFad@<double@> @></code> 的地方调用这个函数，这是一个2d表，其指数分别代表正交点和向量分量。我们用 <code>Wminus[q]</code> 作为最后一个参数来调用这个函数；对2d表进行下标会产生一个代表1d向量的临时访问器对象，这正是我们在这里想要的。问题是，根据C++ 1998和2003标准，临时访问器对象不能被绑定到一个函数的非静态引用参数上，就像我们在这里希望的那样（这个问题将在下一个标准中以rvalue引用的形式得到解决）。 我们在这里把输出参数变成常量，是因为<i>accessor</i>对象是常量，而不是它所指向的表：那个表仍然可以被写到。然而，这个黑客是不愉快的，因为它限制了可以作为这个函数的模板参数的数据类型：一个普通的向量是不行的，因为当标记为  <code>const</code>  时，不能被写入。由于目前没有好的解决方案，我们将采用这里显示的务实的，甚至是不漂亮的解决方案。

    template <typename DataVector> 
    static void 
    compute_Wminus(const std::array<BoundaryKind, n_components> &boundary_kind, 
                   const Tensor<1, dim> &                        normal_vector, 
                   const DataVector &                            Wplus, 
                   const Vector<double> &boundary_values, 
                   const DataVector &    Wminus) 
    { 
      for (unsigned int c = 0; c < n_components; c++) 
        switch (boundary_kind[c]) 
          { 
            case inflow_boundary: 
              { 
                Wminus[c] = boundary_values(c); 
                break; 
              } 

            case outflow_boundary: 
              { 
                Wminus[c] = Wplus[c]; 
                break; 
              } 

// 规定的压力边界条件有点复杂，因为即使压力是规定的，我们在这里真正设定的是能量分量，它将取决于速度和压力。因此，尽管这似乎是一个Dirichlet类型的边界条件，但我们得到了能量对速度和密度的敏感性（除非这些也被规定了）。

            case pressure_boundary: 
              { 
                const typename DataVector::value_type density = 
                  (boundary_kind[density_component] == inflow_boundary ? 
                     boundary_values(density_component) : 
                     Wplus[density_component]); 

                typename DataVector::value_type kinetic_energy = 0; 
                for (unsigned int d = 0; d < dim; ++d) 
                  if (boundary_kind[d] == inflow_boundary) 
                    kinetic_energy += boundary_values(d) * boundary_values(d); 
                  else 
                    kinetic_energy += Wplus[d] * Wplus[d]; 
                kinetic_energy *= 1. / 2. / density; 

                Wminus[c] = 
                  boundary_values(c) / (gas_gamma - 1.0) + kinetic_energy; 

                break; 
              } 

            case no_penetration_boundary: 
              { 

// 我们规定了速度（我们在这里处理的是一个特定的分量，所以速度的平均值是与表面法线正交的。 这就形成了整个速度分量的敏感度。

                typename DataVector::value_type vdotn = 0; 
                for (unsigned int d = 0; d < dim; d++) 
                  { 
                    vdotn += Wplus[d] * normal_vector[d]; 
                  } 

                Wminus[c] = Wplus[c] - 2.0 * vdotn * normal_vector[c]; 
                break; 
              } 

            default: 
              Assert(false, ExcNotImplemented()); 
          } 
    } 
// @sect4{EulerEquations::compute_refinement_indicators}  

// 在这个类中，我们也要指定如何细化网格。这个类 <code>ConservationLaw</code> 将使用我们在 <code>EulerEquation</code> 类中提供的所有信息，对于它所求解的特定守恒定律是不可知的：它甚至不关心一个求解向量有多少个分量。因此，它不可能知道合理的细化指标是什么。另一方面，在这里我们知道，或者至少我们可以想出一个合理的选择：我们简单地看一下密度的梯度，然后计算  $\eta_K=\log\left(1+|\nabla\rho(x_K)|\right)$  ，其中  $x_K$  是单元格  $K$  的中心。

// 当然也有很多同样合理的细化指标，但这个指标确实如此，而且很容易计算。

    static void 
    compute_refinement_indicators(const DoFHandler<dim> &dof_handler, 
                                  const Mapping<dim> &   mapping, 
                                  const Vector<double> & solution, 
                                  Vector<double> &       refinement_indicators) 
    { 
      const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 
      std::vector<unsigned int> dofs(dofs_per_cell); 

      const QMidpoint<dim> quadrature_formula; 
      const UpdateFlags    update_flags = update_gradients; 
      FEValues<dim>        fe_v(mapping, 
                         dof_handler.get_fe(), 
                         quadrature_formula, 
                         update_flags); 

      std::vector<std::vector<Tensor<1, dim>>> dU( 
        1, std::vector<Tensor<1, dim>>(n_components)); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        { 
          const unsigned int cell_no = cell->active_cell_index(); 
          fe_v.reinit(cell); 
          fe_v.get_function_gradients(solution, dU); 

          refinement_indicators(cell_no) = std::log( 
            1 + std::sqrt(dU[0][density_component] * dU[0][density_component])); 
        } 
    } 

//  @sect4{EulerEquations::Postprocessor}  

// 最后，我们声明一个实现数据组件后处理的类。这个类解决的问题是，我们使用的欧拉方程的表述中的变量是保守的而不是物理形式的：它们是动量密度  $\mathbf m=\rho\mathbf v$  ，密度  $\rho$  ，和能量密度  $E$  。我们还想把速度  $\mathbf v=\frac{\mathbf m}{\rho}$  和压力  $p=(\gamma-1)(E-\frac{1}{2} \rho |\mathbf v|^2)$  放入我们的输出文件中。

// 此外，我们还想增加生成Schlieren图的可能性。Schlieren图是一种将冲击和其他尖锐界面可视化的方法。Schlieren "这个词是一个德语单词，可以翻译成 "条纹"--不过，用一个例子来解释可能更简单：比如说，当你把高浓度的酒精或透明的盐水倒入水中时，你会看到schlieren；这两种物质的颜色相同，但它们的折射率不同，所以在它们完全混合之前，光线会沿着弯曲的光线穿过混合物，如果你看它，会导致亮度变化。这就是 "分光"。类似的效果发生在可压缩流中，因为折射率取决于气体的压力（以及因此的密度）。

// 这个词的起源是指三维体积的二维投影（我们看到的是三维流体的二维图片）。在计算流体力学中，我们可以通过考虑其原因来了解这种效应：密度变化。因此，Schlieren图是通过绘制 $s=|\nabla \rho|^2$ 产生的；显然， $s$ 在冲击和其他高度动态的地方很大。如果用户需要（通过在输入文件中指定），我们希望除了上面列出的其他派生量之外，还能生成这些希里伦图。

// 从解决我们问题的数量中计算出派生数量，并将其输出到数据文件中的算法的实现依赖于DataPostprocessor类。它有大量的文档，该类的其他用途也可以在  step-29  中找到。因此，我们避免了大量的评论。

    class Postprocessor : public DataPostprocessor<dim> 
    { 
    public: 
      Postprocessor(const bool do_schlieren_plot); 

      virtual void evaluate_vector_field( 
        const DataPostprocessorInputs::Vector<dim> &inputs, 
        std::vector<Vector<double>> &computed_quantities) const override; 

      virtual std::vector<std::string> get_names() const override; 

      virtual std::vector< 
        DataComponentInterpretation::DataComponentInterpretation> 
      get_data_component_interpretation() const override; 

      virtual UpdateFlags get_needed_update_flags() const override; 

    private: 
      const bool do_schlieren_plot; 
    }; 
  }; 

  template <int dim> 
  const double EulerEquations<dim>::gas_gamma = 1.4; 

  template <int dim> 
  EulerEquations<dim>::Postprocessor::Postprocessor( 
    const bool do_schlieren_plot) 
    : do_schlieren_plot(do_schlieren_plot) 
  {} 

// 这是唯一值得评论的函数。在生成图形输出时，DataOut和相关的类将在每个单元格上调用这个函数，以获取每个正交点的值、梯度、Hessians和法向量（如果我们在处理面）。请注意，每个正交点的数据本身就是矢量值，即保守变量。我们在这里要做的是计算每个正交点上我们感兴趣的量。注意，为此我们可以忽略Hessians（"inputs.solution_hessians"）和法向量（"inputs.normals"）。

  template <int dim> 
  void EulerEquations<dim>::Postprocessor::evaluate_vector_field( 
    const DataPostprocessorInputs::Vector<dim> &inputs, 
    std::vector<Vector<double>> &               computed_quantities) const 
  { 

// 在函数的开始，让我们确保所有的变量都有正确的大小，这样我们就可以访问各个向量元素，而不必怀疑我们是否可能读或写无效的元素；我们还检查 <code>solution_gradients</code> 向量只包含我们真正需要的数据（系统知道这个，因为我们在下面的 <code>get_needed_update_flags()</code> 函数中这样说）。对于内向量，我们检查至少外向量的第一个元素有正确的内部大小。

    const unsigned int n_quadrature_points = inputs.solution_values.size(); 

    if (do_schlieren_plot == true) 
      Assert(inputs.solution_gradients.size() == n_quadrature_points, 
             ExcInternalError()); 

    Assert(computed_quantities.size() == n_quadrature_points, 
           ExcInternalError()); 

    Assert(inputs.solution_values[0].size() == n_components, 
           ExcInternalError()); 

    if (do_schlieren_plot == true) 
      { 
        Assert(computed_quantities[0].size() == dim + 2, ExcInternalError()); 
      } 
    else 
      { 
        Assert(computed_quantities[0].size() == dim + 1, ExcInternalError()); 
      } 

// 然后在所有的正交点上循环，在那里做我们的工作。这段代码应该是不言自明的。输出变量的顺序首先是 <code>dim</code> 速度，然后是压力，如果需要的话，还可以是SCHLIEREN图。请注意，我们尝试使用 <code>first_momentum_component</code> 和 <code>density_component</code> 的信息，对输入向量中的变量顺序进行通用处理。

    for (unsigned int q = 0; q < n_quadrature_points; ++q) 
      { 
        const double density = inputs.solution_values[q](density_component); 

        for (unsigned int d = 0; d < dim; ++d) 
          computed_quantities[q](d) = 
            inputs.solution_values[q](first_momentum_component + d) / density; 

        computed_quantities[q](dim) = 
          compute_pressure(inputs.solution_values[q]); 

        if (do_schlieren_plot == true) 
          computed_quantities[q](dim + 1) = 
            inputs.solution_gradients[q][density_component] * 
            inputs.solution_gradients[q][density_component]; 
      } 
  } 

  template <int dim> 
  std::vector<std::string> EulerEquations<dim>::Postprocessor::get_names() const 
  { 
    std::vector<std::string> names; 
    for (unsigned int d = 0; d < dim; ++d) 
      names.emplace_back("velocity"); 
    names.emplace_back("pressure"); 

    if (do_schlieren_plot == true) 
      names.emplace_back("schlieren_plot"); 

    return names; 
  } 

  template <int dim> 
  std::vector<DataComponentInterpretation::DataComponentInterpretation> 
  EulerEquations<dim>::Postprocessor::get_data_component_interpretation() const 
  { 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      interpretation(dim, 
                     DataComponentInterpretation::component_is_part_of_vector); 

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 

    if (do_schlieren_plot == true) 
      interpretation.push_back( 
        DataComponentInterpretation::component_is_scalar); 

    return interpretation; 
  } 

  template <int dim> 
  UpdateFlags 
  EulerEquations<dim>::Postprocessor::get_needed_update_flags() const 
  { 
    if (do_schlieren_plot == true) 
      return update_values | update_gradients; 
    else 
      return update_values; 
  } 
// @sect3{Run time parameter handling}  

// 我们接下来的工作是定义一些包含运行时参数的类（例如求解器的公差、迭代次数、稳定参数等等）。我们可以在主类中完成这项工作，但我们将其与主类分开，以使程序更加模块化和易于阅读。所有与运行时参数有关的东西都在以下命名空间中，而程序逻辑则在主类中。

// 我们将把运行时参数分成几个独立的结构，我们将把这些结构全部放在一个命名空间  <code>Parameters</code>  中。在这些类中，有几个类将参数分组，用于单独的组，比如用于求解器、网格细化或输出。这些类中的每一个都有函数  <code>declare_parameters()</code>  和  <code>parse_parameters()</code>  ，分别在ParameterHandler对象中声明参数子段和条目，并从这样的对象中检索实际参数值。这些类在ParameterHandler的子段中声明它们的所有参数。

// 以下命名空间的最后一个类结合了前面所有的类，从它们派生出来，并负责处理输入文件顶层的一些条目，以及其他一些奇特的条目，这些条目在子段中太短了，不值得本身有一个结构。

// 这里值得指出的是一件事。下面这些类中没有一个构造函数可以初始化各种成员变量。不过这不是问题，因为我们将从输入文件中读取这些类中声明的所有变量（或者间接地：一个ParameterHandler对象将从那里读取，我们将从这个对象中获取数值），它们将以这种方式被初始化。如果输入文件中根本没有指定某个变量，这也不是问题。在这种情况下，ParameterHandler类将简单地采取默认值，这个默认值是在声明下面这些类的 <code>declare_parameters()</code> 函数中的一个条目时指定的。

  namespace Parameters 
  { 
// @sect4{Parameters::Solver}  

// 这些类中的第一个是关于线性内部求解器的参数。它提供的参数表明使用哪种求解器（GMRES作为一般非对称不定式系统的求解器，或稀疏直接求解器），要产生的输出量，以及各种调整阈值不完全LU分解（ILUT）的参数，我们使用它作为GMRES的预处理器。

// 特别是，ILUT需要以下参数。

// - ilut_fill：形成ILU分解时要增加的额外条目数

// - ilut_atol, ilut_rtol: 在形成预处理程序时，对于某些问题，不好的条件（或者只是运气不好）会导致预处理程序的条件很差。 因此，将对角线扰动添加到原始矩阵中，并为这个稍好的矩阵形成预处理程序会有帮助。 ATOL是一个绝对扰动，在形成预处理之前加到对角线上，RTOL是一个比例因子  $rtol \geq 1$  。

// - ilut_drop: ILUT将放弃任何幅度小于此值的数值。 这是一种管理该预处理程序所使用的内存量的方法。

// 每个参数的含义在  ParameterHandler::declare_entry  调用的第三个参数中也有简要说明  <code>declare_parameters()</code>  。

    struct Solver 
    { 
      enum SolverType 
      { 
        gmres, 
        direct 
      }; 
      SolverType solver; 

      enum OutputType 
      { 
        quiet, 
        verbose 
      }; 
      OutputType output; 

      double linear_residual; 
      int    max_iterations; 

      double ilut_fill; 
      double ilut_atol; 
      double ilut_rtol; 
      double ilut_drop; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Solver::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("linear solver"); 
      { 
        prm.declare_entry( 
          "output", 
          "quiet", 
          Patterns::Selection("quiet|verbose"), 
          "State whether output from solver runs should be printed. " 
          "Choices are <quiet|verbose>."); 
        prm.declare_entry("method", 
                          "gmres", 
                          Patterns::Selection("gmres|direct"), 
                          "The kind of solver for the linear system. " 
                          "Choices are <gmres|direct>."); 
        prm.declare_entry("residual", 
                          "1e-10", 
                          Patterns::Double(), 
                          "Linear solver residual"); 
        prm.declare_entry("max iters", 
                          "300", 
                          Patterns::Integer(), 
                          "Maximum solver iterations"); 
        prm.declare_entry("ilut fill", 
                          "2", 
                          Patterns::Double(), 
                          "Ilut preconditioner fill"); 
        prm.declare_entry("ilut absolute tolerance", 
                          "1e-9", 
                          Patterns::Double(), 
                          "Ilut preconditioner tolerance"); 
        prm.declare_entry("ilut relative tolerance", 
                          "1.1", 
                          Patterns::Double(), 
                          "Ilut relative tolerance"); 
        prm.declare_entry("ilut drop tolerance", 
                          "1e-10", 
                          Patterns::Double(), 
                          "Ilut drop tolerance"); 
      } 
      prm.leave_subsection(); 
    } 

    void Solver::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("linear solver"); 
      { 
        const std::string op = prm.get("output"); 
        if (op == "verbose") 
          output = verbose; 
        if (op == "quiet") 
          output = quiet; 

        const std::string sv = prm.get("method"); 
        if (sv == "direct") 
          solver = direct; 
        else if (sv == "gmres") 
          solver = gmres; 

        linear_residual = prm.get_double("residual"); 
        max_iterations  = prm.get_integer("max iters"); 
        ilut_fill       = prm.get_double("ilut fill"); 
        ilut_atol       = prm.get_double("ilut absolute tolerance"); 
        ilut_rtol       = prm.get_double("ilut relative tolerance"); 
        ilut_drop       = prm.get_double("ilut drop tolerance"); 
      } 
      prm.leave_subsection(); 
    } 

//  @sect4{Parameters::Refinement}  

// 同样的，这里有几个参数决定了网格如何被细化（以及是否要被细化）。关于冲击参数的具体作用，请看下面的网格细化函数。

    struct Refinement 
    { 
      bool   do_refine; 
      double shock_val; 
      double shock_levels; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Refinement::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("refinement"); 
      { 
        prm.declare_entry("refinement", 
                          "true", 
                          Patterns::Bool(), 
                          "Whether to perform mesh refinement or not"); 
        prm.declare_entry("refinement fraction", 
                          "0.1", 
                          Patterns::Double(), 
                          "Fraction of high refinement"); 
        prm.declare_entry("unrefinement fraction", 
                          "0.1", 
                          Patterns::Double(), 
                          "Fraction of low unrefinement"); 
        prm.declare_entry("max elements", 
                          "1000000", 
                          Patterns::Double(), 
                          "maximum number of elements"); 
        prm.declare_entry("shock value", 
                          "4.0", 
                          Patterns::Double(), 
                          "value for shock indicator"); 
        prm.declare_entry("shock levels", 
                          "3.0", 
                          Patterns::Double(), 
                          "number of shock refinement levels"); 
      } 
      prm.leave_subsection(); 
    } 

    void Refinement::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("refinement"); 
      { 
        do_refine    = prm.get_bool("refinement"); 
        shock_val    = prm.get_double("shock value"); 
        shock_levels = prm.get_double("shock levels"); 
      } 
      prm.leave_subsection(); 
    } 

//  @sect4{Parameters::Flux}  

// 接下来是关于通量修改的部分，使其更加稳定。特别是提供了两个选项来稳定Lax-Friedrichs通量：要么选择 $\mathbf{H}(\mathbf{a},\mathbf{b},\mathbf{n}) = \frac{1}{2}(\mathbf{F}(\mathbf{a})\cdot \mathbf{n} + \mathbf{F}(\mathbf{b})\cdot \mathbf{n} + \alpha (\mathbf{a} - \mathbf{b}))$ ，其中 $\alpha$ 是在输入文件中指定的一个固定数字，要么 $\alpha$ 是一个与网格有关的值。在后一种情况下，它被选择为 $\frac{h}{2\delta T}$ ，其中 $h$ 是施加流量的面的直径， $\delta T$ 是当前的时间步长。

    struct Flux 
    { 
      enum StabilizationKind 
      { 
        constant, 
        mesh_dependent 
      }; 
      StabilizationKind stabilization_kind; 

      double stabilization_value; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Flux::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("flux"); 
      { 
        prm.declare_entry( 
          "stab", 
          "mesh", 
          Patterns::Selection("constant|mesh"), 
          "Whether to use a constant stabilization parameter or " 
          "a mesh-dependent one"); 
        prm.declare_entry("stab value", 
                          "1", 
                          Patterns::Double(), 
                          "alpha stabilization"); 
      } 
      prm.leave_subsection(); 
    } 

    void Flux::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("flux"); 
      { 
        const std::string stab = prm.get("stab"); 
        if (stab == "constant") 
          stabilization_kind = constant; 
        else if (stab == "mesh") 
          stabilization_kind = mesh_dependent; 
        else 
          AssertThrow(false, ExcNotImplemented()); 

        stabilization_value = prm.get_double("stab value"); 
      } 
      prm.leave_subsection(); 
    } 

//  @sect4{Parameters::Output}  

// 然后是关于输出参数的部分。我们提供产生Schlieren图（密度的平方梯度，一种可视化冲击前沿的工具），以及图形输出的时间间隔，以防我们不希望每个时间步骤都有输出文件。

    struct Output 
    { 
      bool   schlieren_plot; 
      double output_step; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Output::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("output"); 
      { 
        prm.declare_entry("schlieren plot", 
                          "true", 
                          Patterns::Bool(), 
                          "Whether or not to produce schlieren plots"); 
        prm.declare_entry("step", 
                          "-1", 
                          Patterns::Double(), 
                          "Output once per this period"); 
      } 
      prm.leave_subsection(); 
    } 

    void Output::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("output"); 
      { 
        schlieren_plot = prm.get_bool("schlieren plot"); 
        output_step    = prm.get_double("step"); 
      } 
      prm.leave_subsection(); 
    } 

//  @sect4{Parameters::AllParameters}  

// 最后是将这一切结合起来的类。它自己声明了一些参数，主要是参数文件顶层的参数，以及一些太小的部分，以至于不值得有自己的类。它还包含了所有实际上与空间维度有关的东西，比如初始或边界条件。

// 因为这个类是由上面所有的类派生出来的，所以 <code>declare_parameters()</code> and <code>parse_parameters()</code> 函数也会调用基类的相应函数。

// 注意这个类也处理输入文件中指定的初始和边界条件的声明。为此，在这两种情况下，都有像 "w_0值 "这样的条目，它代表了 $x,y,z$ 方面的表达式，将初始或边界条件描述为一个公式，随后将由FunctionParser类来解析。类似的表达方式还有 "w_1"、"w_2 "等，表示欧拉系统的 <code>dim+2</code> 守恒变量。同样，我们允许在输入文件中最多使用 <code>max_n_boundaries</code> 个边界指标，这些边界指标中的每一个都可以与流入、流出或压力边界条件相关联，同质的边界条件要分别为每个组件和每个边界指标指定。

// 用来存储边界指标的数据结构有点复杂。它是一个 <code>max_n_boundaries</code> 元素的数组，表示将被接受的边界指标的范围。对于这个数组中的每个条目，我们在 <code>BoundaryCondition</code> 结构中存储一对数据：首先是一个大小为 <code>n_components</code> 的数组，对于解向量的每个分量，它表明它是流入、流出还是其他类型的边界，其次是一个FunctionParser对象，它一次描述了这个边界ID的解向量的所有分量。

//  <code>BoundaryCondition</code> 结构需要一个构造器，因为我们需要在构造时告诉函数解析器对象它要描述多少个向量分量。因此，这个初始化不能等到我们在后面的 <code>AllParameters::parse_parameters()</code> 中实际设置FunctionParser对象所代表的公式。

// 由于必须在构造时告诉Function对象其向量大小的同样原因，我们必须有一个 <code>AllParameters</code> 类的构造函数，至少要初始化另一个FunctionParser对象，即描述初始条件的对象。

    template <int dim> 
    struct AllParameters : public Solver, 
                           public Refinement, 
                           public Flux, 
                           public Output 
    { 
      static const unsigned int max_n_boundaries = 10; 

      struct BoundaryConditions 
      { 
        std::array<typename EulerEquations<dim>::BoundaryKind, 
                   EulerEquations<dim>::n_components> 
          kind; 

        FunctionParser<dim> values; 

        BoundaryConditions(); 
      }; 

      AllParameters(); 

      double diffusion_power; 

      double time_step, final_time; 
      double theta; 
      bool   is_stationary; 

      std::string mesh_filename; 

 
      BoundaryConditions  boundary_conditions[max_n_boundaries]; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    template <int dim> 
    AllParameters<dim>::BoundaryConditions::BoundaryConditions() 
      : values(EulerEquations<dim>::n_components) 
    { 
      std::fill(kind.begin(), 
                kind.end(), 
                EulerEquations<dim>::no_penetration_boundary); 
    } 

    template <int dim> 
    AllParameters<dim>::AllParameters() 
      : diffusion_power(0.) 
      , time_step(1.) 
      , final_time(1.) 
      , theta(.5) 
      , is_stationary(true) 
      , initial_conditions(EulerEquations<dim>::n_components) 
    {} 

    template <int dim> 
    void AllParameters<dim>::declare_parameters(ParameterHandler &prm) 
    { 
      prm.declare_entry("mesh", 
                        "grid.inp", 
                        Patterns::Anything(), 
                        "input file name"); 

      prm.declare_entry("diffusion power", 
                        "2.0", 
                        Patterns::Double(), 
                        "power of mesh size for diffusion"); 

      prm.enter_subsection("time stepping"); 
      { 
        prm.declare_entry("time step", 
                          "0.1", 
                          Patterns::Double(0), 
                          "simulation time step"); 
        prm.declare_entry("final time", 
                          "10.0", 
                          Patterns::Double(0), 
                          "simulation end time"); 
        prm.declare_entry("theta scheme value", 
                          "0.5", 
                          Patterns::Double(0, 1), 
                          "value for theta that interpolated between explicit " 
                          "Euler (theta=0), Crank-Nicolson (theta=0.5), and " 
                          "implicit Euler (theta=1)."); 
      } 
      prm.leave_subsection(); 

      for (unsigned int b = 0; b < max_n_boundaries; ++b) 
        { 
          prm.enter_subsection("boundary_" + Utilities::int_to_string(b)); 
          { 
            prm.declare_entry("no penetration", 
                              "false", 
                              Patterns::Bool(), 
                              "whether the named boundary allows gas to " 
                              "penetrate or is a rigid wall"); 

            for (unsigned int di = 0; di < EulerEquations<dim>::n_components; 
                 ++di) 
              { 
                prm.declare_entry("w_" + Utilities::int_to_string(di), 
                                  "outflow", 
                                  Patterns::Selection( 
                                    "inflow|outflow|pressure"), 
                                  "<inflow|outflow|pressure>"); 

                prm.declare_entry("w_" + Utilities::int_to_string(di) + 
                                    " value", 
                                  "0.0", 
                                  Patterns::Anything(), 
                                  "expression in x,y,z"); 
              } 
          } 
          prm.leave_subsection(); 
        } 

      prm.enter_subsection("initial condition"); 
      { 
        for (unsigned int di = 0; di < EulerEquations<dim>::n_components; ++di) 
          prm.declare_entry("w_" + Utilities::int_to_string(di) + " value", 
                            "0.0", 
                            Patterns::Anything(), 
                            "expression in x,y,z"); 
      } 
      prm.leave_subsection(); 

      Parameters::Solver::declare_parameters(prm); 
      Parameters::Refinement::declare_parameters(prm); 
      Parameters::Flux::declare_parameters(prm); 
      Parameters::Output::declare_parameters(prm); 
    } 

    template <int dim> 
    void AllParameters<dim>::parse_parameters(ParameterHandler &prm) 
    { 
      mesh_filename   = prm.get("mesh"); 
      diffusion_power = prm.get_double("diffusion power"); 

      prm.enter_subsection("time stepping"); 
      { 
        time_step = prm.get_double("time step"); 
        if (time_step == 0) 
          { 
            is_stationary = true; 
            time_step     = 1.0; 
            final_time    = 1.0; 
          } 
        else 
          is_stationary = false; 

        final_time = prm.get_double("final time"); 
        theta      = prm.get_double("theta scheme value"); 
      } 
      prm.leave_subsection(); 

      for (unsigned int boundary_id = 0; boundary_id < max_n_boundaries; 
           ++boundary_id) 
        { 
          prm.enter_subsection("boundary_" + 
                               Utilities::int_to_string(boundary_id)); 
          { 
            std::vector<std::string> expressions( 
              EulerEquations<dim>::n_components, "0.0"); 

            const bool no_penetration = prm.get_bool("no penetration"); 

            for (unsigned int di = 0; di < EulerEquations<dim>::n_components; 
                 ++di) 
              { 
                const std::string boundary_type = 
                  prm.get("w_" + Utilities::int_to_string(di)); 

                if ((di < dim) && (no_penetration == true)) 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::no_penetration_boundary; 
                else if (boundary_type == "inflow") 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::inflow_boundary; 
                else if (boundary_type == "pressure") 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::pressure_boundary; 
                else if (boundary_type == "outflow") 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::outflow_boundary; 
                else 
                  AssertThrow(false, ExcNotImplemented()); 

                expressions[di] = 
                  prm.get("w_" + Utilities::int_to_string(di) + " value"); 
              } 

            boundary_conditions[boundary_id].values.initialize( 
              FunctionParser<dim>::default_variable_names(), 
              expressions, 
              std::map<std::string, double>()); 
          } 
          prm.leave_subsection(); 
        } 

      prm.enter_subsection("initial condition"); 
      { 
        std::vector<std::string> expressions(EulerEquations<dim>::n_components, 
                                             "0.0"); 
        for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) 
          expressions[di] = 
            prm.get("w_" + Utilities::int_to_string(di) + " value"); 
        initial_conditions.initialize( 
          FunctionParser<dim>::default_variable_names(), 
          expressions, 
          std::map<std::string, double>()); 
      } 
      prm.leave_subsection(); 

      Parameters::Solver::parse_parameters(prm); 
      Parameters::Refinement::parse_parameters(prm); 
      Parameters::Flux::parse_parameters(prm); 
      Parameters::Output::parse_parameters(prm); 
    } 
  } // namespace Parameters 

//  @sect3{Conservation law class}  

// 这里终于出现了一个类，它实际上是对我们上面定义的所有欧拉方程和参数的具体内容做了一些事情。公共接口与以往基本相同（构造函数现在需要一个文件名来读取参数，这个文件名在命令行中传递）。私有函数接口也与通常的安排非常相似， <code>assemble_system</code> 函数被分成三个部分：一个包含所有单元的主循环，然后分别调用另外两个单元和面的积分。

  template <int dim> 
  class ConservationLaw 
  { 
  public: 
    ConservationLaw(const char *input_filename); 
    void run(); 

  private: 
    void setup_system(); 

    void assemble_system(); 
    void assemble_cell_term(const FEValues<dim> &                       fe_v, 
                            const std::vector<types::global_dof_index> &dofs); 
    void assemble_face_term( 
      const unsigned int                          face_no, 
      const FEFaceValuesBase<dim> &               fe_v, 
      const FEFaceValuesBase<dim> &               fe_v_neighbor, 
      const std::vector<types::global_dof_index> &dofs, 
      const std::vector<types::global_dof_index> &dofs_neighbor, 
      const bool                                  external_face, 
      const unsigned int                          boundary_id, 
      const double                                face_diameter); 

    std::pair<unsigned int, double> solve(Vector<double> &solution); 

    void compute_refinement_indicators(Vector<double> &indicator) const; 
    void refine_grid(const Vector<double> &indicator); 

    void output_results() const; 

// 前面的几个成员变量也相当标准。请注意，我们定义了一个映射对象，在整个程序中组装术语时使用（我们将把它交给每个FEValues和FEFaceValues对象）；我们使用的映射只是标准的 $Q_1$ 映射--换句话说，没有什么花哨的东西--但是在这里声明一个映射并在整个程序中使用它将使以后在有必要时改变它更加简单。事实上，这一点相当重要：众所周知，对于欧拉方程的跨音速模拟，如果边界近似没有足够高的阶数，计算就不会收敛，即使像 $h\rightarrow 0$ 那样。

    Triangulation<dim>   triangulation; 
    const MappingQ1<dim> mapping; 

    const FESystem<dim> fe; 
    DoFHandler<dim>     dof_handler; 

    const QGauss<dim>     quadrature; 
    const QGauss<dim - 1> face_quadrature; 

// 接下来是一些数据向量，对应于前一个时间步骤的解决方案（ <code>old_solution</code> ），当前解决方案的最佳猜测（ <code>current_solution</code> ；我们说<i>guess</i>是因为计算它的牛顿迭代可能还没有收敛，而 <code>old_solution</code> 是指前一个时间步骤的完全收敛的最终结果），以及下一个时间步骤的解决方案的预测器，通过将当前和之前的解决方案推算到未来一个时间步骤计算。

    Vector<double> old_solution; 
    Vector<double> current_solution; 
    Vector<double> predictor; 

    Vector<double> right_hand_side; 

// 这一组最后的成员变量（除了最下面的保存所有运行时参数的对象和一个屏幕输出流，它只在要求verbose输出的情况下打印一些东西）处理我们在这个程序中与Trilinos库的接口，该库为我们提供了线性求解器。与在 step-17 和 step-18 中包括PETSc矩阵类似，我们需要做的是创建一个Trilinos稀疏矩阵而不是标准的deal.II类。该系统矩阵在每个牛顿步骤中被用于雅各布系数。由于我们不打算并行运行这个程序（不过用Trilinos数据结构也不难），所以我们不必考虑其他的事情，比如分配自由度。

    TrilinosWrappers::SparseMatrix system_matrix; 

    Parameters::AllParameters<dim> parameters; 
    ConditionalOStream             verbose_cout; 
  }; 
// @sect4{ConservationLaw::ConservationLaw}  

// 关于构造函数没有什么可说的。基本上，它读取输入文件并将解析后的值填充到参数对象中。

  template <int dim> 
  ConservationLaw<dim>::ConservationLaw(const char *input_filename) 
    : mapping() 
    , fe(FE_Q<dim>(1), EulerEquations<dim>::n_components) 
    , dof_handler(triangulation) 
    , quadrature(fe.degree + 1) 
    , face_quadrature(fe.degree + 1) 
    , verbose_cout(std::cout, false) 
  { 
    ParameterHandler prm; 
    Parameters::AllParameters<dim>::declare_parameters(prm); 

    prm.parse_input(input_filename); 
    parameters.parse_parameters(prm); 

    verbose_cout.set_condition(parameters.output == 
                               Parameters::Solver::verbose); 
  } 

//  @sect4{ConservationLaw::setup_system}  

// 每次改变网格时都会调用下面这个（简单的）函数。它所做的就是根据我们在之前所有的教程程序中生成的稀疏模式来调整特里诺斯矩阵的大小。

  template <int dim> 
  void ConservationLaw<dim>::setup_system() 
  { 
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

    system_matrix.reinit(dsp); 
  } 
// @sect4{ConservationLaw::assemble_system}  

// 这个和下面两个函数是这个程序的核心。它们将牛顿方法应用于非线性守恒方程组所产生的线性系统组合起来。

// 第一个函数将所有的装配部件放在一个例行程序中，为每个单元格/面分配正确的部件。 对这些对象的装配的实际实现是在以下函数中完成的。

// 在函数的顶部，我们做了常规的内务处理：分配FEValues、FEFaceValues和FESubfaceValues对象，这些对象对单元、面和子面（在不同细化级别的相邻单元的情况下）进行积分。请注意，我们并不需要所有这些对象的所有信息（如值、梯度或正交点的实际位置），所以我们只让FEValues类通过指定最小的UpdateFlags集来获得实际需要的信息。例如，当使用邻接单元的FEFaceValues对象时，我们只需要形状值。给定一个特定的面，正交点和 <code>JxW</code> 值与当前单元格相同，法向量已知为当前单元格的法向量的负值。

  template <int dim> 
  void ConservationLaw<dim>::assemble_system() 
  { 
    const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell); 
    std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell); 

    const UpdateFlags update_flags = update_values | update_gradients | 
                                     update_quadrature_points | 
                                     update_JxW_values, 
                      face_update_flags = 
                        update_values | update_quadrature_points | 
                        update_JxW_values | update_normal_vectors, 
                      neighbor_face_update_flags = update_values; 

    FEValues<dim>        fe_v(mapping, fe, quadrature, update_flags); 
    FEFaceValues<dim>    fe_v_face(mapping, 
                                fe, 
                                face_quadrature, 
                                face_update_flags); 
    FESubfaceValues<dim> fe_v_subface(mapping, 
                                      fe, 
                                      face_quadrature, 
                                      face_update_flags); 
    FEFaceValues<dim>    fe_v_face_neighbor(mapping, 
                                         fe, 
                                         face_quadrature, 
                                         neighbor_face_update_flags); 
    FESubfaceValues<dim> fe_v_subface_neighbor(mapping, 
                                               fe, 
                                               face_quadrature, 
                                               neighbor_face_update_flags); 

// 然后循环所有单元，初始化当前单元的FEValues对象，并调用在此单元上组装问题的函数。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_v.reinit(cell); 
        cell->get_dof_indices(dof_indices); 

        assemble_cell_term(fe_v, dof_indices); 

// 然后在这个单元的所有面上循环。 如果一个面是外部边界的一部分，那么就在那里集合边界条件（ <code>assemble_face_terms</code> 的第五个参数表示我们是在外部面还是内部面工作；如果是外部面，表示邻居自由度指数的第四个参数被忽略，所以我们传递一个空向量）。

        for (const auto face_no : cell->face_indices()) 
          if (cell->at_boundary(face_no)) 
            { 
              fe_v_face.reinit(cell, face_no); 
              assemble_face_term(face_no, 
                                 fe_v_face, 
                                 fe_v_face, 
                                 dof_indices, 
                                 std::vector<types::global_dof_index>(), 
                                 true, 
                                 cell->face(face_no)->boundary_id(), 
                                 cell->face(face_no)->diameter()); 
            } 

// 另一种情况是，我们正在处理一个内部面。我们需要区分两种情况：这是在同一细化水平的两个单元之间的正常面，和在不同细化水平的两个单元之间的面。

// 在第一种情况下，我们不需要做什么：我们使用的是连续有限元，在这种情况下，面条款不会出现在双线性表格中。第二种情况通常也不会导致面条款，如果我们强烈地执行悬挂节点约束的话（就像到目前为止，只要我们使用连续有限元的所有教程程序一样--这种执行是由AffineConstraints类和 DoFTools::make_hanging_node_constraints). 一起完成的）。 然而，在当前程序中，我们选择在不同细化水平的单元之间的面弱地执行连续性，原因有二。(i)因为我们可以，更重要的是(ii)因为我们必须通过AffineConstraints类的操作，将我们用来计算牛顿矩阵元素的自动微分穿起来。这是有可能的，但不是微不足道的，所以我们选择了这种替代方法。

// 需要决定的是我们坐在两个不同细化水平的单元之间的接口的哪一边。

// 让我们先来看看邻居更精细的情况。然后，我们必须在当前单元格的面的子代上循环，并在每个子代上进行整合。我们在代码中加入了几个断言，以确保我们试图找出邻居的哪个子面与当前单元格的某个子面重合的推理是正确的--有点防御性的编程永远不会有坏处。

// 然后我们调用对面进行整合的函数；由于这是一个内部面，第五个参数是假的，第六个参数被忽略了，所以我们再次传递一个无效的值。

          else 
            { 
              if (cell->neighbor(face_no)->has_children()) 
                { 
                  const unsigned int neighbor2 = 
                    cell->neighbor_of_neighbor(face_no); 

                  for (unsigned int subface_no = 0; 
                       subface_no < cell->face(face_no)->n_children(); 
                       ++subface_no) 
                    { 
                      const typename DoFHandler<dim>::active_cell_iterator 
                        neighbor_child = 
                          cell->neighbor_child_on_subface(face_no, subface_no); 

                      Assert(neighbor_child->face(neighbor2) == 
                               cell->face(face_no)->child(subface_no), 
                             ExcInternalError()); 
                      Assert(neighbor_child->is_active(), ExcInternalError()); 

                      fe_v_subface.reinit(cell, face_no, subface_no); 
                      fe_v_face_neighbor.reinit(neighbor_child, neighbor2); 

                      neighbor_child->get_dof_indices(dof_indices_neighbor); 

                      assemble_face_term( 
                        face_no, 
                        fe_v_subface, 
                        fe_v_face_neighbor, 
                        dof_indices, 
                        dof_indices_neighbor, 
                        false, 
                        numbers::invalid_unsigned_int, 
                        neighbor_child->face(neighbor2)->diameter()); 
                    } 
                } 

// 我们必须关注的另一种可能性是邻居是否比当前单元更粗（特别是，由于每个面只有一个悬挂节点的通常限制，邻居必须正好比当前单元更粗一级，这是我们用断言检查的）。同样，我们在这个接口上进行整合。

              else if (cell->neighbor(face_no)->level() != cell->level()) 
                { 
                  const typename DoFHandler<dim>::cell_iterator neighbor = 
                    cell->neighbor(face_no); 
                  Assert(neighbor->level() == cell->level() - 1, 
                         ExcInternalError()); 

                  neighbor->get_dof_indices(dof_indices_neighbor); 

                  const std::pair<unsigned int, unsigned int> faceno_subfaceno = 
                    cell->neighbor_of_coarser_neighbor(face_no); 
                  const unsigned int neighbor_face_no = faceno_subfaceno.first, 
                                     neighbor_subface_no = 
                                       faceno_subfaceno.second; 

                  Assert(neighbor->neighbor_child_on_subface( 
                           neighbor_face_no, neighbor_subface_no) == cell, 
                         ExcInternalError()); 

                  fe_v_face.reinit(cell, face_no); 
                  fe_v_subface_neighbor.reinit(neighbor, 
                                               neighbor_face_no, 
                                               neighbor_subface_no); 

                  assemble_face_term(face_no, 
                                     fe_v_face, 
                                     fe_v_subface_neighbor, 
                                     dof_indices, 
                                     dof_indices_neighbor, 
                                     false, 
                                     numbers::invalid_unsigned_int, 
                                     cell->face(face_no)->diameter()); 
                } 
            } 
      } 
  } 
// @sect4{ConservationLaw::assemble_cell_term}  

// 这个函数通过计算残差的单元部分来组装单元项，将其负数加到右手边的向量上，并将其相对于局部变量的导数加到雅各布系数（即牛顿矩阵）上。回顾一下，单元格对残差的贡献为 $R_i = \left(\frac{\mathbf{w}^{k}_{n+1} - \mathbf{w}_n}{\delta t} , \mathbf{z}_i \right)_K $ 。
// $ + \theta \mathbf{B}(\mathbf{w}^{k}_{n+1})(\mathbf{z}_i)_K $  
// $ + (1-\theta) \mathbf{B}(\mathbf{w}_{n}) (\mathbf{z}_i)_K $ ，其中 $\mathbf{B}(\mathbf{w})(\mathbf{z}_i)_K = - \left(\mathbf{F}(\mathbf{w}),\nabla\mathbf{z}_i\right)_K $  。
// $ + h^{\eta}(\nabla \mathbf{w} , \nabla \mathbf{z}_i)_K $  
// $ - (\mathbf{G}(\mathbf {w}), \mathbf{z}_i)_K $ 为 $\mathbf{w} = \mathbf{w}^k_{n+1}$ 和 $\mathbf{w} = \mathbf{w}_{n}$  ， $\mathbf{z}_i$ 为 $i$ 的第1个向量值测试函数。  此外，标量积 $\left(\mathbf{F}(\mathbf{w}), \nabla\mathbf{z}_i\right)_K$ 可以理解为 $\int_K \sum_{c=1}^{\text{n\_components}}  \sum_{d=1}^{\text{dim}} \mathbf{F}(\mathbf{w})_{cd} \frac{\partial z^c_i}{x_d}$ ，其中 $z^c_i$ 是 $i$ 第1个测试函数的 $c$ 分量。

// 在这个函数的顶部，我们做了一些常规的内务工作，即分配一些我们以后需要的局部变量。特别是，我们将分配一些变量来保存 $k$ 次牛顿迭代后的当前解 $W_{n+1}^k$ （变量 <code>W</code> ）和上一时间步长的解 $W_{n}$ （变量 <code>W_old</code> ）的值。

// 除此以外，我们还需要当前变量的梯度。 我们必须计算这些是有点遗憾的，我们几乎不需要。 一个简单的守恒定律的好处是，通量一般不涉及任何梯度。 然而，我们确实需要这些梯度，用于扩散稳定化。

// 我们存储这些变量的实际格式需要一些解释。首先，我们需要解向量的 <code>EulerEquations::n_components</code> 分量在每个正交点的数值。这就构成了一个二维表，我们使用deal.II的表类（这比 <code>std::vector@<std::vector@<T@> @></code> 更有效，因为它只需要分配一次内存，而不是为外向量的每个元素分配一次）。同样地，梯度是一个三维表，Table类也支持。

// 其次，我们想使用自动微分。为此，我们使用 Sacado::Fad::DFad 模板来计算所有我们想计算导数的变量。这包括当前解和正交点的梯度（是自由度的线性组合），以及由它们计算出来的所有东西，如残差，但不包括前一个时间步长的解。这些变量都可以在函数的第一部分找到，同时还有一个变量，我们将用它来存储残差的一个分量的导数。

  template <int dim> 
  void ConservationLaw<dim>::assemble_cell_term( 
    const FEValues<dim> &                       fe_v, 
    const std::vector<types::global_dof_index> &dof_indices) 
  { 
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell; 
    const unsigned int n_q_points    = fe_v.n_quadrature_points; 

    Table<2, Sacado::Fad::DFad<double>> W(n_q_points, 
                                          EulerEquations<dim>::n_components); 

    Table<2, double> W_old(n_q_points, EulerEquations<dim>::n_components); 

    Table<3, Sacado::Fad::DFad<double>> grad_W( 
      n_q_points, EulerEquations<dim>::n_components, dim); 

    Table<3, double> grad_W_old(n_q_points, 
                                EulerEquations<dim>::n_components, 
                                dim); 

    std::vector<double> residual_derivatives(dofs_per_cell); 

// 接下来，我们必须定义自变量，我们将尝试通过解决一个牛顿步骤来确定自变量。这些自变量是局部自由度的值，我们在这里提取。

    std::vector<Sacado::Fad::DFad<double>> independent_local_dof_values( 
      dofs_per_cell); 
    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
      independent_local_dof_values[i] = current_solution(dof_indices[i]); 

// 下一步包含了所有的魔法：我们宣布自分变量的一个子集为独立自由度，而所有其他的变量仍然是依存函数。这些正是刚刚提取的局部自由度。所有引用它们的计算（无论是直接还是间接）都将积累与这些变量有关的敏感度。

// 为了将这些变量标记为独立变量，下面的方法可以起到作用，将 <code>independent_local_dof_values[i]</code> 标记为总共 <code>dofs_per_cell</code> 中的 $i$ 个独立变量。

    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
      independent_local_dof_values[i].diff(i, dofs_per_cell); 

// 在所有这些声明之后，让我们实际计算一些东西。首先， <code>W</code>, <code>W_old</code>, <code>grad_W</code> 和 <code>grad_W_old</code> 的值，我们可以通过使用公式 $W(x_q)=\sum_i \mathbf W_i \Phi_i(x_q)$ 从局部DoF值计算出来，其中 $\mathbf W_i$ 是解向量（局部部分）的第 $i$ 项，而 $\Phi_i(x_q)$ 是在正交点 $x_q$ 评估的第 $i$ 个矢量值的形状函数的值。梯度可以用类似的方法来计算。

// 理想情况下，我们可以通过调用类似 FEValues::get_function_values 和 FEValues::get_function_gradients, 的东西来计算这些信息，但是由于（i）我们必须为此扩展FEValues类，以及（ii）我们不想让整个 <code>old_solution</code> 矢量fad类型，只有局部单元变量，我们明确编码上面的循环。在这之前，我们增加一个循环，将所有的fad变量初始化为零。

    for (unsigned int q = 0; q < n_q_points; ++q) 
      for (unsigned int c = 0; c < EulerEquations<dim>::n_components; ++c) 
        { 
          W[q][c]     = 0; 
          W_old[q][c] = 0; 
          for (unsigned int d = 0; d < dim; ++d) 
            { 
              grad_W[q][c][d]     = 0; 
              grad_W_old[q][c][d] = 0; 
            } 
        } 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        { 
          const unsigned int c = 
            fe_v.get_fe().system_to_component_index(i).first; 

          W[q][c] += independent_local_dof_values[i] * 
                     fe_v.shape_value_component(i, q, c); 
          W_old[q][c] += 
            old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, c); 

          for (unsigned int d = 0; d < dim; d++) 
            { 
              grad_W[q][c][d] += independent_local_dof_values[i] * 
                                 fe_v.shape_grad_component(i, q, c)[d]; 
              grad_W_old[q][c][d] += old_solution(dof_indices[i]) * 
                                     fe_v.shape_grad_component(i, q, c)[d]; 
            } 
        } 

// 接下来，为了计算单元贡献，我们需要在所有正交点评估 $\mathbf{F}({\mathbf w}^k_{n+1})$  ,  $\mathbf{G}({\mathbf w}^k_{n+1})$  和  $\mathbf{F}({\mathbf w}_n)$  ,  $\mathbf{G}({\mathbf w}_n)$  。为了存储这些，我们还需要分配一点内存。请注意，我们以自分变量的方式计算通量矩阵和右手边，这样以后就可以很容易地从中计算出雅各布贡献。

    std::vector<ndarray<Sacado::Fad::DFad<double>, 
                        EulerEquations<dim>::n_components, 
                        dim>> 
      flux(n_q_points); 

    std::vector<ndarray<double, EulerEquations<dim>::n_components, dim>> 
      flux_old(n_q_points); 

 
      std::array<Sacado::Fad::DFad<double>, EulerEquations<dim>::n_components>> 
      forcing(n_q_points); 

    std::vector<std::array<double, EulerEquations<dim>::n_components>> 
      forcing_old(n_q_points); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        EulerEquations<dim>::compute_flux_matrix(W_old[q], flux_old[q]); 
        EulerEquations<dim>::compute_forcing_vector(W_old[q], forcing_old[q]); 
        EulerEquations<dim>::compute_flux_matrix(W[q], flux[q]); 
        EulerEquations<dim>::compute_forcing_vector(W[q], forcing[q]); 
      } 

// 我们现在已经有了所有的部件，所以进行组装。 我们有一个通过系统组件的外循环，和一个通过正交点的内循环，在那里我们积累了对 $i$ 的残差 $R_i$ 的贡献。这个残差的一般公式在引言和本函数的顶部给出。然而，考虑到  $i$  第三个（矢量值）测试函数  $\mathbf{z}_i$  实际上只有一个非零分量（关于这个主题的更多信息可以在  @ref  矢量值模块中找到），我们可以把它简化一下。它将由下面的变量 <code>component_i</code> 表示。有了这个，残差项可以重新写成
// @f{eqnarray*}
//  R_i &=&
//  \left(\frac{(\mathbf{w}_{n+1} -
//  \mathbf{w}_n)_{\text{component\_i}}}{\delta
//  t},(\mathbf{z}_i)_{\text{component\_i}}\right)_K
//  \\ &-& \sum_{d=1}^{\text{dim}} \left(  \theta \mathbf{F}
//  ({\mathbf{w}^k_{n+1}})_{\text{component\_i},d} + (1-\theta)
//  \mathbf{F} ({\mathbf{w}_{n}})_{\text{component\_i},d}  ,
//  \frac{\partial(\mathbf{z}_i)_{\text{component\_i}}} {\partial
//  x_d}\right)_K
//  \\ &+& \sum_{d=1}^{\text{dim}} h^{\eta} \left( \theta \frac{\partial
//  (\mathbf{w}^k_{n+1})_{\text{component\_i}}}{\partial x_d} + (1-\theta)
//  \frac{\partial (\mathbf{w}_n)_{\text{component\_i}}}{\partial x_d} ,
//  \frac{\partial (\mathbf{z}_i)_{\text{component\_i}}}{\partial x_d}
//  \right)_K
//  \\ &-& \left( \theta\mathbf{G}({\mathbf{w}^k_n+1} )_{\text{component\_i}}
//  + (1-\theta)\mathbf{G}({\mathbf{w}_n})_{\text{component\_i}} ,
//  (\mathbf{z}_i)_{\text{component\_i}} \right)_K ,
//  @f}
//  ，其中积分可以理解为通过对正交点求和来评估。

// 我们最初对残差的所有贡献进行正向求和，这样我们就不需要对雅各布项进行负数。 然后，当我们对 <code>right_hand_side</code> 矢量进行求和时，我们就否定了这个残差。

    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
      { 
        Sacado::Fad::DFad<double> R_i = 0; 

        const unsigned int component_i = 
          fe_v.get_fe().system_to_component_index(i).first; 

// 每一行（i）的残差将被累积到这个fad变量中。 在这一行的装配结束时，我们将查询这个变量的敏感度，并将其加入到雅各布系数中。

        for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
          { 
            if (parameters.is_stationary == false) 
              R_i += 1.0 / parameters.time_step * 
                     (W[point][component_i] - W_old[point][component_i]) * 
                     fe_v.shape_value_component(i, point, component_i) * 
                     fe_v.JxW(point); 

            for (unsigned int d = 0; d < dim; d++) 
              R_i -= 
                (parameters.theta * flux[point][component_i][d] + 
                 (1.0 - parameters.theta) * flux_old[point][component_i][d]) * 
                fe_v.shape_grad_component(i, point, component_i)[d] * 
                fe_v.JxW(point); 

            for (unsigned int d = 0; d < dim; d++) 
              R_i += 
                1.0 * 
                std::pow(fe_v.get_cell()->diameter(), 
                         parameters.diffusion_power) * 
                (parameters.theta * grad_W[point][component_i][d] + 
                 (1.0 - parameters.theta) * grad_W_old[point][component_i][d]) * 
                fe_v.shape_grad_component(i, point, component_i)[d] * 
                fe_v.JxW(point); 

            R_i -= 
              (parameters.theta * forcing[point][component_i] + 
               (1.0 - parameters.theta) * forcing_old[point][component_i]) * 
              fe_v.shape_value_component(i, point, component_i) * 
              fe_v.JxW(point); 
          } 

// 在循环结束时，我们必须将敏感度加到矩阵上，并从右手边减去残差。Trilinos FAD数据类型让我们可以使用  <code>R_i.fastAccessDx(k)</code>  访问导数，所以我们将数据存储在一个临时数组中。然后，这些关于整行本地道夫的信息被一次性添加到特里诺斯矩阵中（支持我们选择的数据类型）。

        for (unsigned int k = 0; k < dofs_per_cell; ++k) 
          residual_derivatives[k] = R_i.fastAccessDx(k); 
        system_matrix.add(dof_indices[i], dof_indices, residual_derivatives); 
        right_hand_side(dof_indices[i]) -= R_i.val(); 
      } 
  } 
// @sect4{ConservationLaw::assemble_face_term}  

// 在这里，我们做的事情与前面的函数基本相同。在顶部，我们引入自变量。因为如果我们在两个单元格之间的内部面上工作，也会使用当前的函数，所以自变量不仅是当前单元格上的自由度，而且在内部面上的情况下，也是邻近单元格上的自由度。

  template <int dim> 
  void ConservationLaw<dim>::assemble_face_term( 
    const unsigned int                          face_no, 
    const FEFaceValuesBase<dim> &               fe_v, 
    const FEFaceValuesBase<dim> &               fe_v_neighbor, 
    const std::vector<types::global_dof_index> &dof_indices, 
    const std::vector<types::global_dof_index> &dof_indices_neighbor, 
    const bool                                  external_face, 
    const unsigned int                          boundary_id, 
    const double                                face_diameter) 
  { 
    const unsigned int n_q_points    = fe_v.n_quadrature_points; 
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell; 

    std::vector<Sacado::Fad::DFad<double>> independent_local_dof_values( 
      dofs_per_cell), 
      independent_neighbor_dof_values(external_face == false ? dofs_per_cell : 
                                                               0); 

    const unsigned int n_independent_variables = 
      (external_face == false ? 2 * dofs_per_cell : dofs_per_cell); 

    for (unsigned int i = 0; i < dofs_per_cell; i++) 
      { 
        independent_local_dof_values[i] = current_solution(dof_indices[i]); 
        independent_local_dof_values[i].diff(i, n_independent_variables); 
      } 

    if (external_face == false) 
      for (unsigned int i = 0; i < dofs_per_cell; i++) 
        { 
          independent_neighbor_dof_values[i] = 
            current_solution(dof_indices_neighbor[i]); 
          independent_neighbor_dof_values[i].diff(i + dofs_per_cell, 
                                                  n_independent_variables); 
        } 

// 接下来，我们需要定义保守变量  ${\mathbf W}$  在面的这一侧（  $ {\mathbf W}^+$  ）和另一侧（  ${\mathbf W}^-$  ）的值，对于  ${\mathbf W} = {\mathbf W}^k_{n+1}$  和  ${\mathbf W} = {\mathbf W}_n$  。"这一边 "的值可以用与前一个函数完全相同的方式计算，但注意 <code>fe_v</code> 变量现在是FEFaceValues或FESubfaceValues的类型。

    Table<2, Sacado::Fad::DFad<double>> Wplus( 
      n_q_points, EulerEquations<dim>::n_components), 
      Wminus(n_q_points, EulerEquations<dim>::n_components); 
    Table<2, double> Wplus_old(n_q_points, EulerEquations<dim>::n_components), 
      Wminus_old(n_q_points, EulerEquations<dim>::n_components); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        { 
          const unsigned int component_i = 
            fe_v.get_fe().system_to_component_index(i).first; 
          Wplus[q][component_i] += 
            independent_local_dof_values[i] * 
            fe_v.shape_value_component(i, q, component_i); 
          Wplus_old[q][component_i] += 
            old_solution(dof_indices[i]) * 
            fe_v.shape_value_component(i, q, component_i); 
        } 

// 计算 "对立面 "就比较复杂了。如果这是一个内部面，我们可以像上面那样，简单地使用邻居的独立变量来计算它。

    if (external_face == false) 
      { 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const unsigned int component_i = 
                fe_v_neighbor.get_fe().system_to_component_index(i).first; 
              Wminus[q][component_i] += 
                independent_neighbor_dof_values[i] * 
                fe_v_neighbor.shape_value_component(i, q, component_i); 
              Wminus_old[q][component_i] += 
                old_solution(dof_indices_neighbor[i]) * 
                fe_v_neighbor.shape_value_component(i, q, component_i); 
            } 
      } 

// 另一方面，如果这是一个外部边界面，那么 $\mathbf{W}^-$ 的值将是 $\mathbf{W}^+$ 的函数，或者它们将是规定的，这取决于这里施加的边界条件的种类。

// 为了开始评估，让我们确保为这个边界指定的边界ID是我们在参数对象中实际有数据的一个。接下来，我们对不均匀性的函数对象进行评估。 这有点棘手：一个给定的边界可能同时有规定的和隐含的值。 如果一个特定的成分没有被规定，那么这些值就会被评估为零，并在下面被忽略。

// 剩下的部分由一个实际了解欧拉方程边界条件具体内容的函数完成。请注意，由于我们在这里使用的是fad变量，敏感度将被适当地更新，否则这个过程将是非常复杂的。

    else 
      { 
        Assert(boundary_id < Parameters::AllParameters<dim>::max_n_boundaries, 
               ExcIndexRange(boundary_id, 
                             0, 
                             Parameters::AllParameters<dim>::max_n_boundaries)); 

        std::vector<Vector<double>> boundary_values( 
          n_q_points, Vector<double>(EulerEquations<dim>::n_components)); 
        parameters.boundary_conditions[boundary_id].values.vector_value_list( 
          fe_v.get_quadrature_points(), boundary_values); 

        for (unsigned int q = 0; q < n_q_points; q++) 
          { 
            EulerEquations<dim>::compute_Wminus( 
              parameters.boundary_conditions[boundary_id].kind, 
              fe_v.normal_vector(q), 
              Wplus[q], 
              boundary_values[q], 
              Wminus[q]); 

// 这里我们假设边界类型、边界法向量和边界数据值在时间推进中保持不变。

            EulerEquations<dim>::compute_Wminus( 
              parameters.boundary_conditions[boundary_id].kind, 
              fe_v.normal_vector(q), 
              Wplus_old[q], 
              boundary_values[q], 
              Wminus_old[q]); 
          } 
      } 

// 现在我们有了 $\mathbf w^+$ 和 $\mathbf w^-$ ，我们可以去计算每个正交点的数值通量函数 $\mathbf H(\mathbf w^+,\mathbf w^-, \mathbf n)$ 。在调用这个函数之前，我们还需要确定Lax-Friedrich的稳定性参数。

    std::vector< 
      std::array<Sacado::Fad::DFad<double>, EulerEquations<dim>::n_components>> 
      normal_fluxes(n_q_points); 
    std::vector<std::array<double, EulerEquations<dim>::n_components>> 
      normal_fluxes_old(n_q_points); 

    double alpha; 

 
      { 
        case Parameters::Flux::constant: 
          alpha = parameters.stabilization_value; 
          break; 
        case Parameters::Flux::mesh_dependent: 
          alpha = face_diameter / (2.0 * parameters.time_step); 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
          alpha = 1; 
      } 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        EulerEquations<dim>::numerical_normal_flux( 
          fe_v.normal_vector(q), Wplus[q], Wminus[q], alpha, normal_fluxes[q]); 
        EulerEquations<dim>::numerical_normal_flux(fe_v.normal_vector(q), 
                                                   Wplus_old[q], 
                                                   Wminus_old[q], 
                                                   alpha, 
                                                   normal_fluxes_old[q]); 
      } 

// 现在以与前面函数中的单元格贡献完全相同的方式组装面项。唯一不同的是，如果这是一个内部面，我们还必须考虑到剩余贡献对相邻单元自由度的敏感性。

    std::vector<double> residual_derivatives(dofs_per_cell); 
    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true) 
        { 
          Sacado::Fad::DFad<double> R_i = 0; 

          for (unsigned int point = 0; point < n_q_points; ++point) 
            { 
              const unsigned int component_i = 
                fe_v.get_fe().system_to_component_index(i).first; 

              R_i += (parameters.theta * normal_fluxes[point][component_i] + 
                      (1.0 - parameters.theta) * 
                        normal_fluxes_old[point][component_i]) * 
                     fe_v.shape_value_component(i, point, component_i) * 
                     fe_v.JxW(point); 
            } 

          for (unsigned int k = 0; k < dofs_per_cell; ++k) 
            residual_derivatives[k] = R_i.fastAccessDx(k); 
          system_matrix.add(dof_indices[i], dof_indices, residual_derivatives); 

          if (external_face == false) 
            { 
              for (unsigned int k = 0; k < dofs_per_cell; ++k) 
                residual_derivatives[k] = R_i.fastAccessDx(dofs_per_cell + k); 
              system_matrix.add(dof_indices[i], 
                                dof_indices_neighbor, 
                                residual_derivatives); 
            } 

          right_hand_side(dof_indices[i]) -= R_i.val(); 
        } 
  } 
// @sect4{ConservationLaw::solve}  

// 在这里，我们实际解决线性系统，使用Trilinos的Aztec或Amesos线性求解器。计算的结果将被写入传递给这个函数的参数向量中。其结果是一对迭代次数和最终的线性残差。

  template <int dim> 
  std::pair<unsigned int, double> 
  ConservationLaw<dim>::solve(Vector<double> &newton_update) 
  { 
    switch (parameters.solver) 
      { 

// 如果参数文件指定要使用直接求解器，那么我们就到这里。这个过程很简单，因为deal.II在Trilinos中为Amesos直接求解器提供了一个封装类。我们所要做的就是创建一个求解器控制对象（这里只是一个虚拟对象，因为我们不会进行任何迭代），然后创建直接求解器对象。在实际进行求解时，注意我们没有传递一个预处理程序。无论如何，这对直接求解器来说没有什么意义。 最后我们返回求解器的控制统计信息&mdash;它将告诉我们没有进行任何迭代，并且最终的线性残差为零，这里没有任何可能提供的更好的信息。

        case Parameters::Solver::direct: 
          { 
            SolverControl                                  solver_control(1, 0); 
            TrilinosWrappers::SolverDirect::AdditionalData data( 
              parameters.output == Parameters::Solver::verbose); 
            TrilinosWrappers::SolverDirect direct(solver_control, data); 

            direct.solve(system_matrix, newton_update, right_hand_side); 

            return {solver_control.last_step(), solver_control.last_value()}; 
          } 

// 同样地，如果我们要使用一个迭代求解器，我们使用Aztec的GMRES求解器。我们也可以在这里使用Trilinos的迭代求解器和预处理类，但是我们选择直接使用Aztec的求解器。对于给定的问题，Aztec的内部预处理实现优于deal.II的包装类，所以我们在AztecOO求解器中使用ILU-T预处理，并设置了一堆可以从参数文件中修改的选项。

// 还有两个实际问题。由于我们将右手边和求解向量建立为deal.II向量对象（而不是矩阵，它是一个Trilinos对象），我们必须将Trilinos Epetra向量交给求解器。 幸运的是，他们支持 "视图 "的概念，所以我们只需发送一个指向deal.II向量的指针。我们必须为设置平行分布的向量提供一个Epetra_Map，这只是一个串行的假对象。最简单的方法是要求矩阵提供它的地图，我们要用它为矩阵-向量乘积做好准备。

// 其次，Aztec求解器希望我们传入一个Trilinos Epetra_CrsMatrix，而不是 deal.II包装类本身。所以我们通过trilinos_matrix()命令来访问Trilinos包装类中的实际Trilinos矩阵。Trilinos希望矩阵是非常量的，所以我们必须使用const_cast手动删除常量。

        case Parameters::Solver::gmres: 
          { 
            Epetra_Vector x(View, 
                            system_matrix.trilinos_matrix().DomainMap(), 
                            newton_update.begin()); 
            Epetra_Vector b(View, 
                            system_matrix.trilinos_matrix().RangeMap(), 
                            right_hand_side.begin()); 

            AztecOO solver; 
            solver.SetAztecOption( 
              AZ_output, 
              (parameters.output == Parameters::Solver::quiet ? AZ_none : 
                                                                AZ_all)); 
            solver.SetAztecOption(AZ_solver, AZ_gmres); 
            solver.SetRHS(&b); 
            solver.SetLHS(&x); 

            solver.SetAztecOption(AZ_precond, AZ_dom_decomp); 
            solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut); 
            solver.SetAztecOption(AZ_overlap, 0); 
            solver.SetAztecOption(AZ_reorder, 0); 

            solver.SetAztecParam(AZ_drop, parameters.ilut_drop); 
            solver.SetAztecParam(AZ_ilut_fill, parameters.ilut_fill); 
            solver.SetAztecParam(AZ_athresh, parameters.ilut_atol); 
            solver.SetAztecParam(AZ_rthresh, parameters.ilut_rtol); 

            solver.SetUserMatrix( 
              const_cast<Epetra_CrsMatrix *>(&system_matrix.trilinos_matrix())); 

            solver.Iterate(parameters.max_iterations, 
                           parameters.linear_residual); 

            return {solver.NumIters(), solver.TrueResidual()}; 
          } 
      } 

    Assert(false, ExcNotImplemented()); 
    return {0, 0}; 
  } 
// @sect4{ConservationLaw::compute_refinement_indicators}  

// 这个函数是真正的简单。我们在这里并不假装知道一个好的细化指标会是什么。相反，我们认为 <code>EulerEquation</code> 类会知道这个问题，所以我们只是简单地服从于我们在那里实现的相应函数。

  template <int dim> 
  void ConservationLaw<dim>::compute_refinement_indicators( 
    Vector<double> &refinement_indicators) const 
  { 
    EulerEquations<dim>::compute_refinement_indicators(dof_handler, 
                                                       mapping, 
                                                       predictor, 
                                                       refinement_indicators); 
  } 

//  @sect4{ConservationLaw::refine_grid}  

// 在这里，我们使用之前计算的细化指标来细化网格。在开始的时候，我们在所有的单元格上循环，并标记那些我们认为应该被细化的单元格。

  template <int dim> 
  void 
  ConservationLaw<dim>::refine_grid(const Vector<double> &refinement_indicators) 
  { 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int cell_no = cell->active_cell_index(); 
        cell->clear_coarsen_flag(); 
        cell->clear_refine_flag(); 

        if ((cell->level() < parameters.shock_levels) && 
            (std::fabs(refinement_indicators(cell_no)) > parameters.shock_val)) 
          cell->set_refine_flag(); 
        else if ((cell->level() > 0) && 
                 (std::fabs(refinement_indicators(cell_no)) < 
                  0.75 * parameters.shock_val)) 
          cell->set_coarsen_flag(); 
      } 

// 然后，我们需要在进行细化的同时，将各种解决方案向量从旧网格转移到新网格。SolutionTransfer类是我们的朋友；它有相当丰富的文档，包括例子，所以我们不会对下面的代码做太多评论。最后三行只是把其他一些向量的大小重新设置为现在的正确大小。

    std::vector<Vector<double>> transfer_in; 
    std::vector<Vector<double>> transfer_out; 

    transfer_in.push_back(old_solution); 
    transfer_in.push_back(predictor); 

    triangulation.prepare_coarsening_and_refinement(); 

    SolutionTransfer<dim> soltrans(dof_handler); 
    soltrans.prepare_for_coarsening_and_refinement(transfer_in); 

    triangulation.execute_coarsening_and_refinement(); 

    dof_handler.clear(); 
    dof_handler.distribute_dofs(fe); 

    { 
      Vector<double> new_old_solution(1); 
      Vector<double> new_predictor(1); 

      transfer_out.push_back(new_old_solution); 
      transfer_out.push_back(new_predictor); 
      transfer_out[0].reinit(dof_handler.n_dofs()); 
      transfer_out[1].reinit(dof_handler.n_dofs()); 
    } 

    soltrans.interpolate(transfer_in, transfer_out); 

 
    old_solution = transfer_out[0]; 

    predictor.reinit(transfer_out[1].size()); 
    predictor = transfer_out[1]; 

    current_solution.reinit(dof_handler.n_dofs()); 
    current_solution = old_solution; 
    right_hand_side.reinit(dof_handler.n_dofs()); 
  } 
// @sect4{ConservationLaw::output_results}  

// 现在的这个函数是相当直接的。所有的魔法，包括将数据从保守变量转化为物理变量，都已经被抽象化，并被移到EulerEquations类中，这样在我们想要解决其他双曲守恒定律时就可以被替换。

// 请注意，输出文件的数量是通过保持一个静态变量形式的计数器来确定的，这个计数器在我们第一次来到这个函数时被设置为零，并在每次调用结束时被增加一。

  template <int dim> 
  void ConservationLaw<dim>::output_results() const 
  { 
    typename EulerEquations<dim>::Postprocessor postprocessor( 
      parameters.schlieren_plot); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 

    data_out.add_data_vector(current_solution, 
                             EulerEquations<dim>::component_names(), 
                             DataOut<dim>::type_dof_data, 
                             EulerEquations<dim>::component_interpretation()); 

    data_out.add_data_vector(current_solution, postprocessor); 

    data_out.build_patches(); 

    static unsigned int output_file_number = 0; 
    std::string         filename = 
      "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk"; 
    std::ofstream output(filename); 
    data_out.write_vtk(output); 

    ++output_file_number; 
  } 

//  @sect4{ConservationLaw::run}  

// 这个函数包含了这个程序的顶层逻辑：初始化，时间循环，以及牛顿内部迭代。

// 在开始时，我们读取参数文件指定的网格文件，设置DoFHandler和各种向量，然后在这个网格上插值给定的初始条件。然后我们在初始条件的基础上进行一系列的网格细化，以获得一个已经很适应起始解的网格。在这个过程结束时，我们输出初始解。

  template <int dim> 
  void ConservationLaw<dim>::run() 
  { 
    { 
      GridIn<dim> grid_in; 
      grid_in.attach_triangulation(triangulation); 

      std::ifstream input_file(parameters.mesh_filename); 
      Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str())); 

      grid_in.read_ucd(input_file); 
    } 

    dof_handler.clear(); 
    dof_handler.distribute_dofs(fe); 

// 所有字段的大小。

    old_solution.reinit(dof_handler.n_dofs()); 
    current_solution.reinit(dof_handler.n_dofs()); 
    predictor.reinit(dof_handler.n_dofs()); 
    right_hand_side.reinit(dof_handler.n_dofs()); 

    setup_system(); 

    VectorTools::interpolate(dof_handler, 
                             parameters.initial_conditions, 
                             old_solution); 
    current_solution = old_solution; 
    predictor        = old_solution; 

    if (parameters.do_refine == true) 
      for (unsigned int i = 0; i < parameters.shock_levels; ++i) 
        { 
          Vector<double> refinement_indicators(triangulation.n_active_cells()); 

          compute_refinement_indicators(refinement_indicators); 
          refine_grid(refinement_indicators); 

          setup_system(); 

          VectorTools::interpolate(dof_handler, 
                                   parameters.initial_conditions, 
                                   old_solution); 
          current_solution = old_solution; 
          predictor        = old_solution; 
        } 

    output_results(); 

// 然后我们进入主时间步进循环。在顶部，我们简单地输出一些状态信息，这样就可以跟踪计算的位置，以及显示非线性内部迭代进展的表格的标题。

    Vector<double> newton_update(dof_handler.n_dofs()); 

    double time        = 0; 
    double next_output = time + parameters.output_step; 

    predictor = old_solution; 
    while (time < parameters.final_time) 
      { 
        std::cout << "T=" << time << std::endl 
                  << "   Number of active cells:       " 
                  << triangulation.n_active_cells() << std::endl 
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl 
                  << std::endl; 

        std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl 
                  << "   _____________________________________" << std::endl; 

// 然后是内牛顿迭代，解决每个时间步长的非线性问题。它的工作方式是将矩阵和右手边重置为零，然后组装线性系统。如果右手边的规范足够小，那么我们就宣布牛顿迭代已经收敛了。否则，我们求解线性系统，用牛顿增量更新当前解，并输出收敛信息。最后，我们检查牛顿迭代的次数是否超过了10次的限制--如果超过了，就说明迭代有可能出现了发散，继续迭代也没有什么好处。如果发生这种情况，我们就抛出一个异常，这个异常将在 <code>main()</code> 中被捕获，并在程序终止前显示状态信息。

// 注意，我们写AssertThrow宏的方式基本上等同于写<code>if (!(nonlin_iter  @<=  10)) throw ExcMessage ("No convergence in nonlinear solver");</code>这样的话。唯一显著的区别是，AssertThrow还确保被抛出的异常带有它产生的位置（文件名和行号）的信息。这在这里不是太关键，因为只有一个地方可能发生这种异常；然而，当人们想找出错误发生的地方时，它通常是一个非常有用的工具。

        unsigned int nonlin_iter = 0; 
        current_solution         = predictor; 
        while (true) 
          { 
            system_matrix = 0; 

            right_hand_side = 0; 
            assemble_system(); 

            const double res_norm = right_hand_side.l2_norm(); 
            if (std::fabs(res_norm) < 1e-10) 
              { 
                std::printf("   %-16.3e (converged)\n\n", res_norm); 
                break; 
              } 
            else 
              { 
                newton_update = 0; 

                std::pair<unsigned int, double> convergence = 
                  solve(newton_update); 

                current_solution += newton_update; 

                std::printf("   %-16.3e %04d        %-5.2e\n", 
                            res_norm, 
                            convergence.first, 
                            convergence.second); 
              } 

            ++nonlin_iter; 
            AssertThrow(nonlin_iter <= 10, 
                        ExcMessage("No convergence in nonlinear solver")); 
          } 

// 只有在牛顿迭代已经收敛的情况下，我们才会到达这一点，所以在这里做各种收敛后的任务。

// 首先，我们更新时间，如果需要的话，产生图形输出。然后，我们通过近似 $\mathbf w^{n+1}\approx \mathbf w^n + \delta t \frac{\partial \mathbf w}{\partial t} \approx \mathbf w^n + \delta t \; \frac{\mathbf w^n-\mathbf w^{n-1}}{\delta t} = 2 \mathbf w^n - \mathbf w^{n-1}$ 来更新下一个时间步长的解决方案的预测器，以尝试使适应性更好地工作。 我们的想法是尝试在前面进行细化，而不是步入一个粗略的元素集并抹去旧的解决方案。 这个简单的时间推断器可以完成这个工作。有了这个，如果用户需要的话，我们就可以对网格进行细化，最后继续进行下一个时间步骤。

        time += parameters.time_step; 

        if (parameters.output_step < 0) 
          output_results(); 
        else if (time >= next_output) 
          { 
            output_results(); 
            next_output += parameters.output_step; 
          } 

        predictor = current_solution; 
        predictor.sadd(2.0, -1.0, old_solution); 

        old_solution = current_solution; 

        if (parameters.do_refine == true) 
          { 
            Vector<double> refinement_indicators( 
              triangulation.n_active_cells()); 
            compute_refinement_indicators(refinement_indicators); 

            refine_grid(refinement_indicators); 
            setup_system(); 

            newton_update.reinit(dof_handler.n_dofs()); 
          } 
      } 
  } 
} // namespace Step33 
// @sect3{main()}  

// 下面的``main''函数与前面的例子类似，不需要进行注释。请注意，如果在命令行上没有给出输入文件名，程序就会中止。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step33; 

      if (argc != 2) 
        { 
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl; 
          std::exit(1); 
        } 

      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, dealii::numbers::invalid_unsigned_int); 

      ConservationLaw<2> cons(argv[1]); 
      cons.run(); 
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
    }; 

  return 0; 
} 


