//include/deal.II-translator/numerics/vector_tools_integrate_difference_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_vector_tools_integrate_difference_h
#define dealii_vector_tools_integrate_difference_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/vector_tools_common.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim>
class Quadrature;
template <int dim, int spacedim>
class Triangulation;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
  template <int dim>
  class QCollection;
} // namespace hp


namespace VectorTools
{
  /**
   * @name 函数的评估和错误
   *
   */
  //@{

  /**
   * 计算有限元解的单元误差。
   * 对作为连续函数对象的参考函数和有限元函数之间的差异进行积分。这个函数的结果是向量
   * @p difference ，包含三角形的每个活动单元 $K$
   * 的一个值。这个向量 $d$ 的每个值都等于
   * @f{align*}{
   * d_K = \| u-u_h \|_X
   * @f}
   * 其中 $X$ 表示选择的规范， $u$ 表示精确解。
   * 假设函数 @p exact_solution的分量数与 @p dof.
   * 使用的有限元的分量数一致。
   * 为了计算有限元解的全局误差准则，使用
   * VectorTools::compute_global_error()
   * 与用该函数计算的输出向量。      @param[in]  映射
   * 在积分差分时使用的映射  $u-u_h$  。    @param[in]  dof
   * 描述解向量所在的有限元空间的DoFHandler对象。
   * @param[in]  fe_function 代表数值逼近的节点值的向量  $u_h$
   * 。这个向量需要对应于  @p dof.  所代表的有限元空间
   * @param[in]  exact_solution 用于计算误差的精确解。
   * @param[out]  difference 如上所述计算的值 $d_K$ 的向量。
   * @param[in]  q
   * 用来近似上述积分的正交公式。注意有些正交公式在积分中比其他公式更有用
   * $u-u_h$  。例如，已知 $Q_1$ 对拉普拉斯方程精确解 $u$
   * 的近似 $u_h$
   * 在对应于QGauss(2)对象的2D单元的4个高斯点（或3D的8个点）上特别精确（实际上是超融合，即精确到高阶）。因此，由于QGauss(2)公式只对这些特定点上的两个解进行评估，选择这个正交公式可能表明误差远远小于实际情况。
   * @param[in]  法线 上面显示的应该计算的法线 $X$
   * 。如果规范是 NormType::Hdiv_seminorm,
   * ，那么调用此函数的有限元需要至少有dim向量分量，分歧将在第一个div分量上计算。例如，这对用于混合拉普拉斯方程（
   * step-20 ）和斯托克斯方程（ step-22 ）的有限元有效。
   * @param[in]  权重 附加参数 @p weight 允许评估加权规范。
   * 权重函数可以是标量的，在域中为所有分量平等地建立一个空间上的可变权重。例如，这可用于只对域的部分进行积分。
   * 权重函数也可以是矢量值的，有和有限元一样多的成分。然后，不同的分量得到不同的权重。一个典型的应用是当只计算与一个或一个解变量子集有关的误差时，在这种情况下，其他组件的权重值等于零。ComponentSelectFunction类对这一目的特别有用，因为它提供了这样一个
   * "屏蔽
   * "权重。该权重函数被期望为正值，但负值不会被过滤。这个函数的默认值是一个空指针，被解释为
   * "没有加权函数"，即在整个域中，所有向量成分的权重统一为1。
   * @param[in]  指数 这个值表示在计算 $L^p$  -norms和 $W^{1,p}$
   * -norms时使用的 $p$  。如果选择了 NormType::Lp_norm,
   * NormType::W1p_norm, 或 NormType::W1p_seminorm 以外的 @p norm
   * ，该值会被忽略。
   * 更多信息请参见该命名空间的一般文档。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const Mapping<dim, spacedim> &                           mapping,
    const DoFHandler<dim, spacedim> &                        dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const Quadrature<dim> &                                  q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * 调用 integrate_difference()
   * 函数，见上文，<tt>mapping=MappingQGeneric  @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const DoFHandler<dim, spacedim> &                        dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const Quadrature<dim> &                                  q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * 与上面的hp相同。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const hp::MappingCollection<dim, spacedim> &             mapping,
    const DoFHandler<dim, spacedim> &                        dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const hp::QCollection<dim> &                             q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * 调用 integrate_difference()
   * 函数，见上文，<tt>mapping=MappingQGeneric  @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference(
    const DoFHandler<dim, spacedim> &                        dof,
    const InVector &                                         fe_function,
    const Function<spacedim, typename InVector::value_type> &exact_solution,
    OutVector &                                              difference,
    const hp::QCollection<dim> &                             q,
    const NormType &                                         norm,
    const Function<spacedim, double> *                       weight   = nullptr,
    const double                                             exponent = 2.);

  /**
   * 在每个有<tt>tria.n_active_cells()</tt>项的单元格上取一个错误向量
   * @p cellwise_error ，并返回由 @p norm. 给出的全局错误。  @p
   * cellwise_error 向量通常是由 VectorTools::integrate_difference()
   * 产生的输出，你通常希望为 @p norm  ]的值与你在
   * VectorTools::integrate_difference().
   * 中使用的值相同。如果给定的三角形是
   * parallel::TriangulationBase,
   * 中不对应于本地拥有的单元格的条目被假定为0.0，并使用MPI进行并行还原以计算全局误差。
   * @param  tria 与 @p cellwise_error.
   * 中的条目相对应的活动单元的三角形  @param  cellwise_error
   * 每个活动单元的误差矢量。    @param  norm
   * 要计算的规范类型。    @param  exponent 用于 $L^p$  -norms和
   * $W^{1,p}$  -norms的指数 $p$ 。如果选择 @p norm 以外的
   * NormType::Lp_norm,  NormType::W1p_norm, 或 NormType::W1p_seminorm
   * ，该值会被忽略。
   * @note  为Vector<double>和Vector<float>类型实例化。
   *
   */
  template <int dim, int spacedim, class InVector>
  double
  compute_global_error(const Triangulation<dim, spacedim> &tria,
                       const InVector &                    cellwise_error,
                       const NormType &                    norm,
                       const double                        exponent = 2.);

  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_integrate_difference_h


