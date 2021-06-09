//include/deal.II-translator/grid/tensor_product_manifold_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_tensor_product_manifold_h
#define dealii_tensor_product_manifold_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/manifold.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN



/**
 * @brief Tensor product manifold of two ChartManifolds.
 * 这个流形将结合构造函数中给出的ChartManifold  @p A  和  @p B
 * ，通过建立张量积  $A\otimes B$  ，形成一个新的ChartManifold
 * 。实空间的第一个 @p spacedim_A 维度和图表的第一个 @p
 * chartdim_A 维度将由流形 @p A, 给出，而其余坐标由 @p B.
 * 给出。流形将由<tt>三角法 @<dim,  space_dim_A+space_dim_B
 * @></tt>.  使用。
 * 一个例子是将空间维度为2的SphericalManifold和空间维度为1的FlatManifold结合起来，形成一个圆柱形流形。
 * pull_back()、push_forward()和push_forward_gradient()是通过将输入参数按照给定的维度分割成
 * @p A 和 @p B
 * 的输入并在连接结果之前应用相应的操作来实现的。
 *
 *
 * @note  尺寸参数 @p dim_A 和 @p dim_B 不被使用。
 * @tparam  dim
 * 单元的尺寸（需要与要连接的三角结构的第一个模板参数相匹配。
 * @tparam  dim_A ChartManifold A的尺寸。  @tparam  spacedim_A
 * ChartManifold A的空间尺寸。  @tparam  chartdim_A ChartManifold
 * A的图表尺寸。
 *
 *
 */
template <int dim,
          int dim_A,
          int spacedim_A,
          int chartdim_A,
          int dim_B,
          int spacedim_B,
          int chartdim_B>
class TensorProductManifold
  : public ChartManifold<dim, spacedim_A + spacedim_B, chartdim_A + chartdim_B>
{
public:
  /**
   * 图表尺寸是流形 @p A 和 @p B. 的图表尺寸之和。
   *
   */
  static const unsigned int chartdim = chartdim_A + chartdim_B;
  /**
   * 空间维度是流形 @p A 和 @p B. 的空间维度之和。
   *
   */
  static const unsigned int spacedim = spacedim_A + spacedim_B;

  /**
   * 构建器。
   *
   */
  TensorProductManifold(
    const ChartManifold<dim_A, spacedim_A, chartdim_A> &manifold_A,
    const ChartManifold<dim_B, spacedim_B, chartdim_B> &manifold_B);

  /**
   * 克隆此流形。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim_A + spacedim_B>>
  clone() const override;

  /**
   * 拉回操作。
   *
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const override;

  /**
   * 前推操作。
   *
   */
  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const override;

  /**
   * 梯度。
   *
   */
  virtual DerivativeForm<1, chartdim, spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const override;

private:
  std::unique_ptr<const ChartManifold<dim_A, spacedim_A, chartdim_A>>
    manifold_A;

  std::unique_ptr<const ChartManifold<dim_B, spacedim_B, chartdim_B>>
    manifold_B;
};



 /*------------------Template Implementations------------------------*/ 



namespace internal
{
  namespace TensorProductManifoldImplementation
  {
    template <int dim1, int dim2>
    Tensor<1, dim1 + dim2>
    concat(const Tensor<1, dim1> &p1, const Tensor<1, dim2> &p2)
    {
      Tensor<1, dim1 + dim2> r;
      for (unsigned int d = 0; d < dim1; ++d)
        r[d] = p1[d];
      for (unsigned int d = 0; d < dim2; ++d)
        r[dim1 + d] = p2[d];
      return r;
    }

    template <int dim1, int dim2>
    Point<dim1 + dim2>
    concat(const Point<dim1> &p1, const Point<dim2> &p2)
    {
      Point<dim1 + dim2> r;
      for (unsigned int d = 0; d < dim1; ++d)
        r[d] = p1[d];
      for (unsigned int d = 0; d < dim2; ++d)
        r[dim1 + d] = p2[d];
      return r;
    }

    template <int dim1, int dim2>
    void
    split_point(const Point<dim1 + dim2> &source,
                Point<dim1> &             p1,
                Point<dim2> &             p2)
    {
      for (unsigned int d = 0; d < dim1; ++d)
        p1[d] = source[d];
      for (unsigned int d = 0; d < dim2; ++d)
        p2[d] = source[dim1 + d];
    }

  } // namespace TensorProductManifoldImplementation
} // namespace internal

template <int dim,
          int dim_A,
          int spacedim_A,
          int chartdim_A,
          int dim_B,
          int spacedim_B,
          int chartdim_B>
TensorProductManifold<dim,
                      dim_A,
                      spacedim_A,
                      chartdim_A,
                      dim_B,
                      spacedim_B,
                      chartdim_B>::
  TensorProductManifold(
    const ChartManifold<dim_A, spacedim_A, chartdim_A> &manifold_A,
    const ChartManifold<dim_B, spacedim_B, chartdim_B> &manifold_B)
  : ChartManifold<dim, spacedim_A + spacedim_B, chartdim_A + chartdim_B>(
      internal::TensorProductManifoldImplementation::concat(
        manifold_A.get_periodicity(),
        manifold_B.get_periodicity()))
  , manifold_A(Utilities::dynamic_unique_cast<
               ChartManifold<dim_A, spacedim_A, chartdim_A>,
               Manifold<dim_A, spacedim_A>>(manifold_A.clone()))
  , manifold_B(Utilities::dynamic_unique_cast<
               ChartManifold<dim_B, spacedim_B, chartdim_B>,
               Manifold<dim_B, spacedim_B>>(manifold_B.clone()))
{}

template <int dim,
          int dim_A,
          int spacedim_A,
          int chartdim_A,
          int dim_B,
          int spacedim_B,
          int chartdim_B>
std::unique_ptr<Manifold<dim, spacedim_A + spacedim_B>>
TensorProductManifold<dim,
                      dim_A,
                      spacedim_A,
                      chartdim_A,
                      dim_B,
                      spacedim_B,
                      chartdim_B>::clone() const
{
  return std::make_unique<TensorProductManifold<dim,
                                                dim_A,
                                                spacedim_A,
                                                chartdim_A,
                                                dim_B,
                                                spacedim_B,
                                                chartdim_B>>(*manifold_A,
                                                             *manifold_B);
}

template <int dim,
          int dim_A,
          int spacedim_A,
          int chartdim_A,
          int dim_B,
          int spacedim_B,
          int chartdim_B>
Point<TensorProductManifold<dim,
                            dim_A,
                            spacedim_A,
                            chartdim_A,
                            dim_B,
                            spacedim_B,
                            chartdim_B>::chartdim>
TensorProductManifold<dim,
                      dim_A,
                      spacedim_A,
                      chartdim_A,
                      dim_B,
                      spacedim_B,
                      chartdim_B>::
  pull_back(
    const Point<TensorProductManifold<dim,
                                      dim_A,
                                      spacedim_A,
                                      chartdim_A,
                                      dim_B,
                                      spacedim_B,
                                      chartdim_B>::spacedim> &space_point) const
{
  Point<spacedim_A> space_point_A;
  Point<spacedim_B> space_point_B;
  internal::TensorProductManifoldImplementation::split_point(space_point,
                                                             space_point_A,
                                                             space_point_B);

  Point<chartdim_A> result_A = manifold_A->pull_back(space_point_A);
  Point<chartdim_B> result_B = manifold_B->pull_back(space_point_B);

  return internal::TensorProductManifoldImplementation::concat(result_A,
                                                               result_B);
}

template <int dim,
          int dim_A,
          int spacedim_A,
          int chartdim_A,
          int dim_B,
          int spacedim_B,
          int chartdim_B>
Point<TensorProductManifold<dim,
                            dim_A,
                            spacedim_A,
                            chartdim_A,
                            dim_B,
                            spacedim_B,
                            chartdim_B>::spacedim>
TensorProductManifold<dim,
                      dim_A,
                      spacedim_A,
                      chartdim_A,
                      dim_B,
                      spacedim_B,
                      chartdim_B>::
  push_forward(
    const Point<TensorProductManifold<dim,
                                      dim_A,
                                      spacedim_A,
                                      chartdim_A,
                                      dim_B,
                                      spacedim_B,
                                      chartdim_B>::chartdim> &chart_point) const
{
  Point<chartdim_A> chart_point_A;
  Point<chartdim_B> chart_point_B;
  internal::TensorProductManifoldImplementation::split_point(chart_point,
                                                             chart_point_A,
                                                             chart_point_B);

  Point<spacedim_A> result_A = manifold_A->push_forward(chart_point_A);
  Point<spacedim_B> result_B = manifold_B->push_forward(chart_point_B);

  return internal::TensorProductManifoldImplementation::concat(result_A,
                                                               result_B);
}

template <int dim,
          int dim_A,
          int spacedim_A,
          int chartdim_A,
          int dim_B,
          int spacedim_B,
          int chartdim_B>
DerivativeForm<1,
               TensorProductManifold<dim,
                                     dim_A,
                                     spacedim_A,
                                     chartdim_A,
                                     dim_B,
                                     spacedim_B,
                                     chartdim_B>::chartdim,
               TensorProductManifold<dim,
                                     dim_A,
                                     spacedim_A,
                                     chartdim_A,
                                     dim_B,
                                     spacedim_B,
                                     chartdim_B>::spacedim>

TensorProductManifold<dim,
                      dim_A,
                      spacedim_A,
                      chartdim_A,
                      dim_B,
                      spacedim_B,
                      chartdim_B>::
  push_forward_gradient(
    const Point<TensorProductManifold<dim,
                                      dim_A,
                                      spacedim_A,
                                      chartdim_A,
                                      dim_B,
                                      spacedim_B,
                                      chartdim_B>::chartdim> &chart_point) const
{
  Point<chartdim_A> chart_point_A;
  Point<chartdim_B> chart_point_B;
  internal::TensorProductManifoldImplementation::split_point(chart_point,
                                                             chart_point_A,
                                                             chart_point_B);

  DerivativeForm<1, chartdim_A, spacedim_A> result_A =
    manifold_A->push_forward_gradient(chart_point_A);
  DerivativeForm<1, chartdim_B, spacedim_B> result_B =
    manifold_B->push_forward_gradient(chart_point_B);


  DerivativeForm<1, chartdim, spacedim> result;
  for (unsigned int i = 0; i < chartdim_A; ++i)
    for (unsigned int j = 0; j < spacedim_A; ++j)
      result[j][i] = result_A[j][i];
  for (unsigned int i = 0; i < chartdim_B; ++i)
    for (unsigned int j = 0; j < spacedim_B; ++j)
      result[j + spacedim_A][i + chartdim_A] = result_B[j][i];

  return result;
}



DEAL_II_NAMESPACE_CLOSE

#endif


