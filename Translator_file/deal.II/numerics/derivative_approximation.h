//include/deal.II-translator/numerics/derivative_approximation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_derivative_approximation_h
#define dealii_derivative_approximation_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/synchronous_iterator.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/lac/vector.h>
#ifdef _MSC_VER
#  include <deal.II/dofs/dof_accessor.h>
#endif
#include <utility>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个命名空间提供了一些函数，这些函数通过取相邻单元之间的差值商来计算有限元场的导数的规范的逐一逼近。这是一种相当简单而有效的错误指示形式，因为它可以用相对较少的数值努力来计算，但却能给出一个合理的近似值。
 * 在单元格 $K$
 * 上计算差商的方式如下（这里描述的是有限元场梯度的近似，但对于更高的导数见下文）：让
 * $K'$ 是一个相邻的单元格，让 $y_{K'}=x_{K'}-x_K$
 * ]是两个单元中心之间的距离向量，那么 $ \frac{u_h(x_{K'})
 *
 * - u_h(x_K)}{ \|y_{K'}\| }$ 是方向性导数的近似值 $ \nabla u(x_K)
 * \cdot \frac{y_{K'}}{ \|y_{K'}\| }.$  通过从左边开始乘以
 * $\frac{y_{K'}}{ \|y_{K'}\| }$ 并对所有邻居 $K'$ 求和，我们得到
 * $ \sum_{K'} \left( \frac{y_{K'}}{ \|y_{K'}\|} \frac{y_{K'}^T}{ \|y_{K'}\| }
 * \right) \nabla u(x_K) \approx \sum_{K'} \left( \frac{y_{K'}}{ \|y_{K'}\|}
 * \frac{u_h(x_{K'})
 *
 * - u_h(x_K)}{ \|y_{K'}\| }  \right).$ 因此，如果矩阵 $ Y =  \sum_{K'}
 * \left( \frac{y_{K'}}{\|y_{K'}\|} \frac{y_{K'}^T}{ \|y_{K'}\| } \right)$
 * 是规则的（当向量 $y_{K'}$
 * 到所有邻居跨越整个空间时就是这种情况），我们可以通过
 * $ \nabla u(x_K) \approx Y^{-1} \sum_{K'} \left( \frac{y_{K'}}{\|y_{K'}\|}
 * \frac{u_h(x_{K'})
 *
 * - u_h(x_K)}{ \|y_{K'}\| } \right).$
 * 得到真实梯度的近似值，这是一个容易计算的量。当调用该类的
 * @p approximate_gradient
 * 函数时，每个单元格的返回值是这个梯度近似值的 $l_2$
 * 准则。为了使这个量成为一个有用的量，你可能想把每个元素按各自的单元格大小的正确幂进行缩放。
 * 如果一个单元只有方向向量 $y_K$
 * 不跨越整个空间的邻居，那么这个量的计算一定会失败，因为此时矩阵
 * $Y$
 * 不再是可逆的。如果发生这种情况，你会得到一个类似于这个的错误。
 *
 * @code
 *
 *
 *
 *
 *
 * --------------------------------------------------------
 * An error occurred in line <749>
 * of file <source/numerics/derivative_approximation.cc> in function
 *   void DerivativeApproximation::approximate(...)
 * [with DerivativeDescription = DerivativeApproximation::Gradient<3>, int
 * dim = 3, InputVector = Vector<double>, spacedim = 3]
 * The violated condition was:
 *   determinant(Y) != 0
 * The name and call sequence of the exception was:
 *   ExcInsufficientDirections()
 * Additional Information:
 * (none)
 *
 *
 *
 *
 *
 * --------------------------------------------------------
 * @endcode
 * 正如很容易验证的那样，这种情况只可能发生在非常粗糙的网格上，当一些单元和它们的所有邻居甚至没有被精炼过一次。因此，你应该只在所有单元至少被精炼过一次的情况下才调用这个类别的函数。在实践中，这并不是一个很大的限制。
 *
 *  <h3>Approximation of higher derivatives</h3>
 * 与上面的推理类似，高阶导数的近似值也可以用类似的方式计算。例如，二阶导数的张量由公式
 * $ \nabla^2 u(x_K) \approx Y^{-1} \sum_{K'} \left( \frac{y_{K'}}{\|y_{K'}\|}
 * \otimes \frac{\nabla u_h(x_{K'})
 *
 *
 *
 * - \nabla u_h(x_K)}{ \|y_{K'}\| } \right), $ 近似，其中 $\otimes$
 * 表示两个向量的外积。请注意，与真正的二阶导数张量不同，其近似值不一定是对称的。这是由于在推导过程中，不清楚我们是将
 * $\nabla^2 u y_{KK'}$ 还是 $y_{KK'}^T \nabla^2 u$
 * 项作为投射二阶导数。根据我们的选择，我们会得到二阶导数张量的一个近似值或其转置。为了避免这种模糊性，作为结果，我们采取对称的形式，即近似值及其转置的平均值。
 * 每个单元格的返回值是二阶导数的近似张量的谱准则，即绝对值最大的特征值。这等于每个单元的有限元场的最大曲率，谱规范是与
 * $l_2$ 向量规范相关的矩阵规范。
 * 甚至比二阶导数更高的导数也可以沿着上面暴露的相同思路得到。
 *
 *  <h3>Refinement indicators based on the derivatives</h3>
 * 如果你想在这些导数的近似值的基础上建立一个细化标准，你必须将这一类的结果按网格宽度的一个适当的幂进行缩放。例如，由于
 * $\|u-u_h\|^2_{L_2} \le C h^2 \|\nabla u\|^2_{L_2}$
 * ，可能正确的做法是将指标缩放为 $\eta_K = h \|\nabla u\|_K$
 * ，即 $\eta_K = h^{1+d/2} \|\nabla u\|_{\infty;K}$ ，即正确的幂是
 * $1+d/2$  。
 * 同样的，对于二阶导数，应该选择比梯度高一个的网格大小
 * $h$ 的幂。
 *
 *  <h3>Implementation</h3>
 * 上面显示的梯度和二阶导数张量的近似计算公式是非常相似的。基本区别在于，在一种情况下，有限差分商是一个标量，而在另一种情况下，它是一个矢量。对于更高的导数，这将是一个更高等级的张量。然后我们必须用距离向量
 * $y_{KK'}$ 形成这个差分商的外积，对其进行对称，用矩阵
 * $Y^{-1}$
 * 收缩，并计算其规范。为了使实现更简单并允许代码重用，所有这些取决于要逼近的导数的实际顺序的操作，以及进入差商的量的计算，都被分离到辅助嵌套类中（名称为
 * @p Gradient 和 @p SecondDerivative)
 * ），主算法只是传递一个或其他数据类型并要求它们执行取决于顺序的操作。独立于此的主要框架，如寻找所有活动邻居，或设置矩阵
 * $Y$ 是在主函数 @p approximate. 中完成的。
 * 由于这种操作方式，该类可以很容易地扩展到比目前实现的更高阶导数。基本上，只需要按照导数描述符类
 * @p 梯度和 @p SecondDerivative
 * 的思路实现一个额外的类，用要近似的导数的适当类似物替换各自的别名和函数。
 *
 *
 * @ingroup numerics
 *
 *
 */
namespace DerivativeApproximation
{
  /**
   * 该函数用于获得梯度的近似值。将描述有限元场的DoF处理程序对象和一个节点值向量传递给它，并接收近似梯度的单元格欧氏规范。
   * 最后一个参数是指要计算梯度的解分量。它默认为第一个分量。对于标量元素，这是唯一有效的选择；对于矢量值元素，可以在这里给出零和矢量分量数量之间的任何分量。
   * 在并行计算中， @p solution
   * 向量需要包含本地相关的未知数。
   *
   */
  template <int dim, class InputVector, int spacedim>
  void
  approximate_gradient(const Mapping<dim, spacedim> &   mapping,
                       const DoFHandler<dim, spacedim> &dof,
                       const InputVector &              solution,
                       Vector<float> &                  derivative_norm,
                       const unsigned int               component = 0);

  /**
   * 用<tt>mapping=MappingQGeneric  @<dim@>(1)</tt>. 调用上述函数。
   *
   */
  template <int dim, class InputVector, int spacedim>
  void
  approximate_gradient(const DoFHandler<dim, spacedim> &dof,
                       const InputVector &              solution,
                       Vector<float> &                  derivative_norm,
                       const unsigned int               component = 0);

  /**
   * 这个函数是上面那个函数的类似物，计算二阶导数张量的有限差分近似值。将描述有限元场的DoF处理程序对象、节点值向量传递给它，并接收第二导数张量近似值的单元格谱准则。谱准则是与
   * $l_2$ 向量准则相关的矩阵准则。
   * 最后一个参数是指要计算梯度的解分量。它默认为第一个分量。对于标量元素，这是唯一有效的选择；对于矢量值元素，这里可以给出零和矢量分量数量之间的任何分量。
   * 在并行计算中， @p solution
   * 向量需要包含本地相关的未知数。
   *
   */
  template <int dim, class InputVector, int spacedim>
  void
  approximate_second_derivative(const Mapping<dim, spacedim> &   mapping,
                                const DoFHandler<dim, spacedim> &dof,
                                const InputVector &              solution,
                                Vector<float> &    derivative_norm,
                                const unsigned int component = 0);

  /**
   * 用<tt>mapping=MappingQGeneric  @<dim@>(1)</tt>. 调用上述函数。
   *
   */
  template <int dim, class InputVector, int spacedim>
  void
  approximate_second_derivative(const DoFHandler<dim, spacedim> &dof,
                                const InputVector &              solution,
                                Vector<float> &    derivative_norm,
                                const unsigned int component = 0);

  /**
   * 这个函数计算<tt>阶</tt>-阶近似导数，并返回单个单元的完整张量。
   * 最后一个参数是指要计算梯度的解分量。它默认为第一个分量。对于标量元素，这是唯一有效的选择；对于矢量值元素，可以在这里给出零和矢量分量数量之间的任何分量。
   * 在并行计算中， @p solution
   * 向量需要包含本地相关的未知数。
   *
   */
  template <int dim, int spacedim, class InputVector, int order>
  void
  approximate_derivative_tensor(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const InputVector &              solution,
#ifndef _MSC_VER
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
#else
    const TriaActiveIterator<dealii::DoFCellAccessor<dim, spacedim, false>>
      &cell,
#endif
    Tensor<order, dim> &derivative,
    const unsigned int  component = 0);

  /**
   * 同上，<tt>mapping=MappingQGeneric  @<dim@>(1)</tt>. 。
   *
   */
  template <int dim, int spacedim, class InputVector, int order>
  void
  approximate_derivative_tensor(
    const DoFHandler<dim, spacedim> &dof,
    const InputVector &              solution,
#ifndef _MSC_VER
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
#else
    const TriaActiveIterator<dealii::DoFCellAccessor<dim, spacedim, false>>
      &cell,
#endif
    Tensor<order, dim> &derivative,
    const unsigned int  component = 0);

  /**
   * 返回导数的规范。
   *
   */
  template <int dim, int order>
  double
  derivative_norm(const Tensor<order, dim> &derivative);

  /**
   * 异常情况
   *
   */
  DeclException2(ExcVectorLengthVsNActiveCells,
                 int,
                 int,
                 << "The output vector needs to have a size equal "
                    "to the number of active cells of your triangulation "
                    "but has length "
                 << arg1 << "There are " << arg2
                 << " active cells in your triangulation.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcInsufficientDirections,
                   "While computing a finite difference approximation to "
                   "derivatives, the algorithm encountered a cell on which "
                   "the number of linearly "
                   "independent directions that span the matrix Y (discussed "
                   "in the documentation of the DerivativeApproximation "
                   "class) is not equal to dim. The matrix Y then is "
                   "rank deficient and can not be inverted. A common reason "
                   "why this might be happening is if a cell has neither "
                   "left/right (or up/down, or front/back) neighbors, for "
                   "example because the mesh is too coarse.");
} // namespace DerivativeApproximation



DEAL_II_NAMESPACE_CLOSE

#endif


