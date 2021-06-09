//include/deal.II-translator/base/function_tools_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2021 by the deal.II authors
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

#ifndef dealii_function_tools_h
#define dealii_function_tools_h


#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/function.h>

DEAL_II_NAMESPACE_OPEN

namespace FunctionTools
{
  /**
   * 估算一个函数 $f$
   * 的值和每个梯度分量的边界，在一个BoundingBox上，通过从盒子中心开始的二阶泰勒多项式对其进行逼近。
   * 每个下限和上限都以 <code>std::pair<double, double></code>
   * 的形式返回，第一项是下限， $L$ ，第二项是上限， $U$
   * ，即 $f(x) \in [L, U]$  。
   * 函数值、梯度和Hessian是在盒子中心计算的。
   * 然后，函数值的边界被估计为  $f(x) \in [f(x_c)
   *
   * - F, f(x_c) + F]$  ，其中  $F = \sum_i |\partial_i f(x_c)| h_i + 1/2
   * \sum_i \sum_j |\partial_i \partial_j f(x_c)| h_i h_j$  。    这里，
   * $h_i$ 是盒子在 $i$
   * 个坐标方向上边长的一半，这是我们外推的距离。梯度分量的界限估计类似于
   * $\partial_i f \in [\partial_i f(x_c)
   *
   * - G_i, \partial_i f(x_c) + G_i]$  ，其中  $G_i = \sum_j |\partial_i
   * \partial_j f(x_c)| h_j$  。    如果函数有1个以上的分量， @p
   * component 参数可以用来指定计算哪个函数分量的边界。
   *
   */
  template <int dim>
  void
  taylor_estimate_function_bounds(
    const Function<dim> &                       function,
    const BoundingBox<dim> &                    box,
    std::pair<double, double> &                 value_bounds,
    std::array<std::pair<double, double>, dim> &gradient_bounds,
    const unsigned int                          component = 0);

} // namespace FunctionTools
DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_function_tools_h */ 


