//include/deal.II-translator/numerics/vector_tools_constraints_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_vector_tools_constraints_h
#define dealii_vector_tools_constraints_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <map>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
struct StaticMappingQ1;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;


namespace VectorTools
{
  /**
   * @name  内插和投影
   *
   */
  //@{

  /**
   * 这个函数计算对应于 $\vec u \cdot \vec n=\vec u_\Gamma \cdot \vec
   * n$ 形式的边界条件的约束，即法向通量约束，其中 $\vec
   * u$ 是一个矢量值的解变量， $\vec u_\Gamma$
   * 是一个规定的矢量场，我们希望其法向分量等于解的法向分量。
   * 这些条件完全具有AffineConstraints类所处理的形式，即它们将边界自由度的<i>linear
   * combination</i>与相应的值（约束的不均匀性）联系起来。因此，当前函数创建了一个写进AffineConstraints容器的约束列表。这个对象可能已经有了一些内容，例如来自悬挂节点的约束，但仍未触及。这些约束条件必须像其他此类约束条件一样应用于线性系统，也就是说，你必须在求解前用约束条件浓缩线性系统，而且你必须在求解后分配解向量。
   * 这个函数处理的情况比
   * VectorTools::compute_no_normal_flux_constraints()
   * 更普遍（它只能处理 $\vec u_\Gamma \cdot \vec n = 0$
   * 的情况，并在 step-31 和 step-32
   * 中使用）。然而，由于适用于该函数的所有内容也作为特例适用于当前的函数，所以下面的讨论与两者都有关。
   * @note  这个函数在1d中没有什么意义，所以如果 @p dim
   * 等于1，它会抛出一个异常。      <h4>Arguments to this
   * function</h4>
   * 这个函数的第二个参数表示有限元中的第一个矢量分量，对应于你要约束的矢量函数。例如，如果我们正在求解2d中的斯托克斯方程，有限元有分量
   * $(u,v,p)$  ，那么如果你打算约束矢量 $(u,v)^T \cdot \vec n =
   * \vec u_\Gamma \cdot \vec n$ ， @p  first_vector_component需要为0。
   * 另一方面，如果我们在三维中求解麦克斯韦方程，有限元有分量
   * $(E_x,E_y,E_z,B_x,B_y,B_z)$ ，我们想要边界条件 $\vec B\cdot \vec
   * n=\vec B_\Gamma\cdot \vec n$ ，那么 @p first_vector_component
   * 将是3。隐含地假设矢量正好有 <code>dim</code>
   * 个分量，这些分量的排序方式与我们通常对坐标方向的排序方式相同，即
   * $x$  -， $y$  -，最后是 $z$
   * -分量。该函数假设，但不能检查，
   * <code>[first_vector_component,first_vector_component+dim)</code>
   * 范围内的矢量分量来自同一个基础有限元。例如，在上面的Stokes例子中，使用
   * <code>FESystem@<dim@>(FE_Q@<dim@>(2), 1, FE_Q@<dim@>(1), dim)</code>
   * 是没有意义的（注意第一个速度矢量分量是 $Q_2$
   * 元素，而其他的都是 $Q_1$
   * 元素），因为在边界上会有定义了 $x$
   * -速度的点，但没有相应的 $y$ -或 $z$ -velocities。
   * 第三个参数表示要强制执行边界条件的边界指标集。请注意，正如下面所解释的，这是少数几个函数之一，在这些函数中，我们只用一个边界指标多次调用该函数，还是用整个边界指标一次调用该函数，是有区别的。
   * 参数四（  @p function_map)  描述了每个边界id的边界函数
   * $\vec u_\Gamma$ 。函数 <code>function_map[id]</code> 用于id为 @p id
   * 的边界，取自集合 @p boundary_ids. 中的每个函数预计都有
   * @p dim 成分，这些成分的使用与 @p first_vector_component.
   * 无关
   * 映射参数用于计算边界点，在边界点上，函数需要从边界描述中请求法向量
   * $\vec n$ 。
   * @note  当在一个AffineConstraints对象中结合具有悬挂节点约束和边界条件的自适应细化网格时，悬挂节点约束应该总是首先被设置，然后是边界条件，因为边界条件不会在已经被约束的自由度的第二次操作中被设置。这可以确保离散化保持所需的一致性。参见  @ref constraints  模块中关于冲突约束的讨论。      <h4>Computing constraints in 2d</h4> 计算这些约束需要一些智慧。主要问题是围绕着法向量是什么这个问题。考虑一下下面的情况。      <p ALIGN="center"
   * >
   @image html no_normal_flux_1.png
   * </p>
   * 这里，我们有两个使用双线性映射（即MappingQGeneric(1)）的单元格。因此，对于每一个单元格来说，法向量都是垂直于直边的。如果顶部和右侧的两条边是为了近似弯曲的边界（如虚线所示），那么两个计算出来的法向量都不等于准确的法向量（尽管随着网格的进一步细化，它们近似于法向量）。更糟糕的是，如果我们用两个单元的法线向量来约束公共顶点的
   * $\vec u \cdot \vec n= \vec u_\Gamma \cdot \vec n$
   * ，那么我们就用两个线性独立的向量来约束向量 $\vec u$
   * ；因此，在这一点上的约束将是 $\vec u=\vec u_\Gamma$
   * （即向量的<i>all</i>分量），这并不是我们想要的结果。
   * 为了处理这种情况，算法的工作方式如下：在我们想要约束
   * $\vec u$
   * 的每一点，我们首先收集相邻单元在这一点可能计算的所有法向量。然后我们不对这些法向量中的<i>each</i>进行约束
   * $\vec u \cdot \vec n=\vec u_\Gamma \cdot \vec n$
   * ，而只对法向量中的<i>average</i>进行约束。在上面的例子中，我们因此只记录了一个单一的约束
   * $\vec u \cdot \vec {\bar n}=\vec u_\Gamma \cdot \vec {\bar n}$  ，其中
   * $\vec {\bar n}$ 是两个指定法向量的平均值。
   * 不幸的是，这还远远不够。考虑一下这里的情况。
   * <p ALIGN="center">
   @image html no_normal_flux_2.png
   * </p>
   * 如果顶边和右边又近似于一个弯曲的边界，而左边的边界是一个单独的边界（例如直的），所以确切的边界在左顶点确实有一个角，那么上面的构造就不起作用了：在这里，我们确实希望在这一点上的约束是
   * $\vec u$
   * （因为相对于左边的法线以及上面的法线矢量的法线速度应该是零），而不是说平均法线矢量方向的速度是零。
   * 因此，我们使用以下启发式方法来确定在某一点计算的所有法向量是否要被平均化：如果同一点的两个法向量是在<i>different</i>单元上计算的，那么它们就要被平均化。这包括上面的第一个例子。如果它们是从同一个单元计算出来的，那么它们不同的事实被认为是表明它们来自边界的不同部分，可能被一个真正的角所连接，因此必须不被平均化。
   * 这个方案有一个问题。例如，如果我们在上面考虑的同一个域，用下面的网格进行离散，那么我们就会陷入麻烦。
   * <p ALIGN="center">
   @image html no_normal_flux_3.png
   * </p> 这里，算法假设边界在面 $F1$ 和 $F2$
   * 连接处没有角，因为在这一点上有两个不同的法向量从不同的单元计算出来。如果你想在这一点上有一个确切的边界角，处理这个问题的唯一方法是给边界的两个部分分配不同的边界指标，并调用这个函数两次，一次针对每个边界指标；这样做每次调用只能得到一个法向量（因为我们一次只考虑一个边界部分），结果是法向量不会被平均化。在笛卡尔网格上的重心角周围使用该函数时，也需要考虑到这种情况。如果法向流边界条件要在非直角坐标系网格上的重心角处执行，我们甚至可能在约束中得到循环，因为一般来说，我们会约束来自两边的不同分量。在这种情况下，首先在重心顶点上设置一个无滑移约束。
   * <h4>Computing constraints in 3d</h4>
   * 在三维中，情况更为复杂。考虑以下情况，我们要在被标记的顶点计算约束。
   * <p ALIGN="center">
   @image html no_normal_flux_4.png
   * </p>
   * 在这里，我们得到四个不同的法向量，其中一个来自在顶点相遇的四个面。尽管它们可能形成一个完整的向量集，但我们的目的不是在这一点上约束向量场的所有成分。相反，我们希望仍然允许切向流动，其中
   * "切向 "一词必须得到适当的定义。
   * 在这种情况下，算法进行如下：对于在这一点上计算了两个切向矢量的每个单元，我们计算不受约束的方向，作为两个切向矢量的外积（如果有必要，乘以减一）。然后我们对这些切向矢量进行平均。最后，我们计算垂直于这个平均切线方向的两个方向的约束。
   * 有些情况下，一个单元贡献了两个切向，而另一个单元只贡献了一个切向；例如，如果左边单元的顶面和正面都属于选定的边界，而只有右边单元的顶面属于边界，就会发生这种情况，也许这表明域的整个正面部分是一个平滑流形，而顶部确实形成了两个独立的流形，在一个山脊上相遇，而且只希望在正面流形和顶部的右边流形上设置法向流动边界条件。在这样的情况下，很难定义应该发生什么。目前的实现只是忽略了只贡献了一个法向的单元的贡献。在所示的例子中，这是可以接受的，因为左边单元的正面的法向量与右边单元的正面提供的法向量是一样的（表面是平面的），但是如果前面的流形是弯曲的，这将是一个问题。无论如何，不清楚在这种情况下该如何进行，忽略单细胞可能是最好的办法了。
   * <h4>Results</h4>因为它能产生很好的图片，这里有两张圆和球体上的矢量场的图片，这个函数计算的约束被应用于此（为了说明问题，我们强制执行零法线通量，这可以更容易地使用
   * VectorTools::compute_no_normal_flux_constraints(),
   * 计算，因为这必须导致<i>tangential</i>矢量场）。      <p
   * ALIGN="center">
   @image html no_normal_flux_5.png
   @image html no_normal_flux_6.png
   * </p>
   * 矢量场在物理上是不合理的，但切向性约束显然是被执行的。事实上，向量场在边界上的某些点是零的，这是创建方式的一个伪命题，它没有被约束在这些点上为零。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim>
  void
  compute_nonzero_normal_flux_constraints(
    const DoFHandler<dim, spacedim> &   dof_handler,
    const unsigned int                  first_vector_component,
    const std::set<types::boundary_id> &boundary_ids,
    const std::map<types::boundary_id, const Function<spacedim, double> *>
      &                           function_map,
    AffineConstraints<double> &   constraints,
    const Mapping<dim, spacedim> &mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 这个函数与compute_nonzero_normal_flux_constraints()函数的作用相同（更多信息见那里），但用于同质法向流约束的更简单情况，即用于施加条件
   * $\vec u \cdot \vec n= 0$  。这个函数在  step-31  和  step-32
   * 中使用。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim>
  void
  compute_no_normal_flux_constraints(
    const DoFHandler<dim, spacedim> &   dof_handler,
    const unsigned int                  first_vector_component,
    const std::set<types::boundary_id> &boundary_ids,
    AffineConstraints<double> &         constraints,
    const Mapping<dim, spacedim> &      mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 计算对应于 $\vec u \times \vec n=\vec u_\Gamma \times \vec n$
   * 形式的边界条件的约束，即切向流约束，其中 $\vec u$
   * 是一个矢量值的解变量， $\vec u_\Gamma$
   * 是规定的矢量场，我们希望其切向分量与解的切向分量相等。这个函数正好约束那些不受
   * VectorTools::compute_no_normal_flux_constraints(),
   * 约束的dim-1矢量值分量，并留下一个不受约束的分量，这个分量受到该函数的约束。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim>
  void
  compute_nonzero_tangential_flux_constraints(
    const DoFHandler<dim, spacedim> &   dof_handler,
    const unsigned int                  first_vector_component,
    const std::set<types::boundary_id> &boundary_ids,
    const std::map<types::boundary_id, const Function<spacedim, double> *>
      &                           function_map,
    AffineConstraints<double> &   constraints,
    const Mapping<dim, spacedim> &mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 与上述同质切向流约束相同。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim>
  void
  compute_normal_flux_constraints(
    const DoFHandler<dim, spacedim> &   dof_handler,
    const unsigned int                  first_vector_component,
    const std::set<types::boundary_id> &boundary_ids,
    AffineConstraints<double> &         constraints,
    const Mapping<dim, spacedim> &      mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_constraints_h


