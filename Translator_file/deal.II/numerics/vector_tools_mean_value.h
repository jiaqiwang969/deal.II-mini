//include/deal.II-translator/numerics/vector_tools_mean_value_0.txt
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

#ifndef dealii_vector_tools_mean_value_h
#define dealii_vector_tools_mean_value_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q1.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class DoFHandler;

namespace VectorTools
{
  /**
   * 均值运算
   *
   */
  //@{

  /**
   * 从一个向量中减去（代数）均值。
   * 这个函数最经常被用作斯托克斯的均值过滤器。
   * 在斯托克斯方程中，速度只有Dirichlet边界的压力只确定到一个常数。这个函数允许减去压力的平均值。它通常在预处理程序中调用，并产生平均值为零的更新。平均值是按照输入矢量给出的自由度值的平均值计算的；它们不以单元面积为权重，也就是说，平均值被计算为
   * $\sum_i v_i$ ，而不是 $\int_\Omega v(x) =
   * \int_\Omega \sum_i
   * v_i \phi_i(x)$ 。然而，后者可以从
   * VectorTools::compute_mean_function, 中得到。    除了对向量 @p v
   * 进行操作外，该函数还需要一个布尔掩码 @p p_select
   * ，该掩码对向量的每个元素都有一个真条目，对这些元素将计算出平均值并随后减去。该参数用于表示解向量中与压力相对应的分量，而避免触及向量的所有其他分量，如速度分量。(但是请注意，该掩码不是对解向量
   * @p v 可能与之相关的有限元的向量分量进行操作的 @ref
   * GlossComponentMask
   * ；相反，它是对整个向量的掩码，而不参考向量元素的含义)。
   * 布尔掩码 @p p_select
   * 有一个空矢量作为默认值，这将被解释为选择所有矢量元素，因此，减去整个矢量上的代数平均值。如果要处理整个向量，这允许在没有布尔掩码的情况下调用这个函数。
   * @note
   * 在使用这个函数过滤掉一个算子的内核的情况下（例如由常数压力组成的斯托克斯算子的空空间），这个函数只对空空间确实由矢量组成的有限元素有意义
   * $(1,1,\ldots,1)^T$
   * 。例如，对于通常的拉格朗日元素就是这种情况，所有形状函数的总和等于恒一的函数。然而，对于其他一些函数则不然：例如，对于FE_DGP元素（Stokes离散中压力的另一个有效选择），每个单元上的第一个形状函数是常数，而更多的元素与它正交（在参考单元上）；因此，所有形状函数之和不等于一，与常数模式相关的向量也不等于
   * $(1,1,\ldots,1)^T$
   * 。对于这样的元素，在减去平均值时，必须使用不同的程序。
   * @warning
   * 这个函数只能用于分布式向量类，前提是布尔掩码为空，即选择整个向量。
   *
   */
  template <typename VectorType>
  void
  subtract_mean_value(VectorType &v, const std::vector<bool> &p_select = {});


  /**
   * 计算解决方案中一个分量的平均值。
   * 这个函数将所选分量在整个域上进行积分并返回结果，即计算
   * $\frac{1}{|\Omega|}\int_\Omega
   * @note
   * 该函数最常用于解决其解只定义到一个常数的问题，例如纯诺伊曼问题或斯托克斯或纳维尔斯托克斯问题中的压力。在这两种情况下，从节点向量中减去当前函数计算出的均值，一般不会产生均值为零的有限元函数的理想结果。事实上，它只对拉格朗日元素有效。对于所有其他元素，你需要计算均值并在评估程序中减去它。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  compute_mean_value(const Mapping<dim, spacedim> &   mapping,
                     const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim> &          quadrature,
                     const VectorType &               v,
                     const unsigned int               component);

  /**
   * 调用另一个compute_mean_value()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  compute_mean_value(const DoFHandler<dim, spacedim> &dof,
                     const Quadrature<dim> &          quadrature,
                     const VectorType &               v,
                     const unsigned int               component);
  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_mean_value_h


