//include/deal.II-translator/fe/fe_rannacher_turek_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


#ifndef dealii_fe_rannacher_turek_h
#define dealii_fe_rannacher_turek_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_rannacher_turek.h>

#include <deal.II/fe/fe_base.h>
#include <deal.II/fe/fe_poly.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * Rannacher-Turek元素的实现。该元用于为斯托克斯方程生成一对稳定的函数空间，而不必像稳定的Taylor-Hood元那样增加速度空间的多项式度数，后者使用
 * $Q_2^d\times Q_1$
 * 对速度和压力进行处理。也就是说，像许多其他不符合要求的元素一样，它也可以用于拉普拉斯方程的离散化。该元素首次描述于R.Rannacher和S.Turek:
 * "Simple non-conforming quadrilateral Stokes element", Numerical Methods for
 * Partial Differential Equations, vol. 8, pp.97-112, 1992。
 * 该元素产生的形状函数一般是不连续的，因此该元素不是
 * $H^1$ 符合的（即，它是一个 "不符合的
 * "元素）。然而，形状函数是以这样一种方式构造的，即沿面的跳动的平均值为零，因此，该元素有<i>some</i>种符合性：符合性元素会有一个点状的零跳动，像FE_DGQ元素这样完全不连续的元素可以有完全任意的跨面跳动值，而当前的元素处于中间位置，因为它的跳动是非零的，但至少是平均值为零。
 * 该元素目前只在维度2中实现，为最低的多项式阶数，并且没有悬挂节点和限制/延长。
 *
 *  <h3>Interpolation</h3> <h4>Node values</h4>
 * @ref GlossNodes  "节点值 "
 * 是面的时刻。 <h4>Generalized support points</h4>
 * 为了计算节点值，我们在每个面上使用QGauss规则。默认情况下，我们使用两点规则来精确地整合Rannacher-Turek函数。但是为了能够以足够的精度插值其他函数，可以在构造函数中调整一个面上使用的正交点的数量。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int dim>
class FE_RannacherTurek : public FE_Poly<dim>
{
public:
  /**
   * 给定 @p order, 的Rannacher-Turek元素的构造函数，使用 @p
   * 每个面上的n_face_support_points正交点进行插值。
   * 注意0阶的元素包含2度的多项式。
   * 该元素目前只在二维中实现了0阶。
   *
   */
  FE_RannacherTurek(const unsigned int order                 = 0,
                    const unsigned int n_face_support_points = 2);

  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

private:
  /**
   * 这个元素的阶数。
   *
   */
  const unsigned int order;

  /**
   * 在插值过程中，每个面上用于评估节点函数的正交点的数量。
   *
   */
  const unsigned int n_face_support_points;

  /**
   * 用于评估节点函数的面的权重。
   *
   */
  std::vector<double> weights;

  /**
   * 计算广义支持点和它们的权重。
   *
   */
  void
  initialize_support_points();
  /**
   * 在建造过程中根据需要返回每个物体的自由度信息。
   *
   */
  std::vector<unsigned int>
  get_dpo_vector();
};


DEAL_II_NAMESPACE_CLOSE

#endif


