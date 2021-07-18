//include/deal.II-translator/base/function_spherical_0.txt
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

#ifndef dealii_function_spherical_h
#define dealii_function_spherical_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * 一个抽象的基类，用于定义球面坐标的标量值函数
   * $f=f(r,\theta,\phi)$
   * 。这个类包裹了从球面坐标到Function基类所使用的笛卡尔坐标系统的数值、梯度和
   * hessians的转换。因此，派生类只需要在球面坐标中实现这些函数（特别是svalue(),
   * sgradient() 和 shessian() ）。角度的约定与
   * GeometricUtilities::Coordinates. 中相同。
   * @note  这个函数目前只对dim==3实现。
   * @ingroup functions
   *
   */
  template <int dim>
  class Spherical : public Function<dim>
  {
  public:
    /**
     * 构造函数，应提供 @p center ，定义坐标系的原点。
     * 请注意，这个函数的组成部分被视为完全独立的数量。
     *
     * - 而不是作为将在不同坐标系中重新解释的向量的组成部分。
     *
     */
    Spherical(const Point<dim> & center       = Point<dim>(),
              const unsigned int n_components = 1);

    /**
     * 返回该函数在给定点的值。
     * 这个函数将给定的点转换为球面坐标，用它调用svalue()，并返回结果。
     *
     */
    virtual double
    value(const Point<dim> & point,
          const unsigned int component = 0) const override;

    /**
     * 返回相对于点 @p p. 的笛卡尔坐标的梯度
     * 这个函数将给定的点转换为球面坐标，用它调用sgradient()，并将结果转换为笛卡尔坐标。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 返回关于点 @p p. 的笛卡尔坐标的Hessian
     * 这个函数将给定的点转换为球面坐标，用它调用sgradient和Shessian()，并将结果转换为笛卡尔坐标。
     *
     */
    virtual SymmetricTensor<2, dim>
    hessian(const Point<dim> & p,
            const unsigned int component = 0) const override;

    /**
     * 返回这个对象的内存消耗估计值，以字节为单位。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

  private:
    /**
     * 返回 @p sp. 点的值 这里， @p sp 是以球面坐标提供的。
     *
     */
    virtual double
    svalue(const std::array<double, dim> &sp,
           const unsigned int             component) const;

    /**
     * 返回球面坐标中的梯度。
     * 返回的对象应该按照以下顺序包含导数。      $\{
     * f_{,r},\, f_{,\theta},\, f_{,\phi}\}$  .
     *
     */
    virtual std::array<double, dim>
    sgradient(const std::array<double, dim> &sp,
              const unsigned int             component) const;

    /**
     * 返回球面坐标中的Hessian。
     * 返回的对象应包含按以下顺序排列的导数。      $\{
     * f_{,rr},\, f_{,\theta\theta},\, f_{,\phi\phi},\, f_{,r\theta},\,
     * f_{,r\phi},\, f_{,\theta\phi}\}$  .
     *
     */
    virtual std::array<double, 6>
    shessian(const std::array<double, dim> &sp,
             const unsigned int             component) const;

    /**
     * 从原点到球面坐标系中心的一个矢量。
     *
     */
    const Tensor<1, dim> coordinate_system_offset;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif


