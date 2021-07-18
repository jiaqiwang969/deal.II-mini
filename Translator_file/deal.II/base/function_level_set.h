//include/deal.II-translator/base/function_level_set_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_function_level_set_h
#define dealii_function_level_set_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  namespace LevelSet
  {
    /**
     * 球体的有符号距离水平集函数。      $\psi(x) = \| x
     *
     * - x^c \|
     *
     * - R$  .     这里， $x^c$ 是球体的中心， $R$
     * 是其半径。因此，这个函数在球体上为零，在以球体为边界的球体内为负数，在
     * $\mathbb{R}^{dim}$ 的其余部分为正数。
     * 这个函数的梯度和Hessian等于  $\partial_i \psi(x) = (x
     *
     * - x^c)/\| x
     *
     * - x^c \|$  ,  $\partial_i \partial_j \psi = \delta_{ij}/\| x
     *
     * - x^c \|
     *
     * - (x_i
     *
     * - x_i^c)(x_j
     *
     * - x_j^c)/\| x
     *
     * - x^c \|^3$  , 其中  $\delta_{ij}$  是Kronecker delta函数。
     * @ingroup functions
     *
     */
    template <int dim>
    class Sphere : public Function<dim>
    {
    public:
      /**
       * 构造函数，获取球体的中心和半径。
       *
       */
      Sphere(const Point<dim> &center = Point<dim>(), const double radius = 1);

      double
      value(const Point<dim> & point,
            const unsigned int component = 0) const override;

      /**
       * @copydoc   Function::gradient()  *  Function::gradient() .
       * @note  梯度在球体的中心是单数。如果输入的 @p point
       * 离中心太近，可能会出现浮点异常，或者梯度中的条目可能是+inf/-inf或+nan/-nan，这取决于该点相对于奇异点的位置。
       *
       */
      Tensor<1, dim>
      gradient(const Point<dim> & point,
               const unsigned int component = 0) const override;

      /**
       * @copydoc   Function::hessian() .
       * @note  Hessian在球体中心是奇点。如果输入的 @p point
       * 离中心太近，可能会出现浮点异常，或者Hessian中的条目可能是+inf/inf或+nan/nan，这取决于该点相对于奇异点的位置。
       *
       */
      SymmetricTensor<2, dim>
      hessian(const Point<dim> & point,
              const unsigned int component = 0) const override;

    private:
      const Point<dim> center;
      const double     radius;
    };


    /**
     * 在 $\mathbb{R}^{dim}$ 中的平面的有符号水平集函数：
     * $\psi(x) = n \cdot (x
     *
     * - x_p)$  。    这里， $n$ 是平面法线， $x_p$
     * 是平面上的一个点。
     * 因此，关于法线的方向，这个函数在平面上是正的，在平面上是零，在平面下是负的。如果法线被归一化，
     * $\psi$ 将是到平面内最近点的有符号距离。
     * @ingroup functions
     *
     */
    template <int dim>
    class Plane : public Function<dim>
    {
    public:
      /**
       * 构造函数，接收平面内的一个点和平面法线。
       *
       */
      Plane(const Point<dim> &point, const Tensor<1, dim> &normal);

      double
      value(const Point<dim> & point,
            const unsigned int component = 0) const override;

      Tensor<1, dim>
      gradient(const Point<dim> &,
               const unsigned int component = 0) const override;

      SymmetricTensor<2, dim>
      hessian(const Point<dim> &,
              const unsigned int component = 0) const override;

    private:
      const Point<dim>     point_in_plane;
      const Tensor<1, dim> normal;
    };

  } // namespace LevelSet
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif  /* CE8BBB3E_B726_40A7_B963_561AE7B84973 */ 


