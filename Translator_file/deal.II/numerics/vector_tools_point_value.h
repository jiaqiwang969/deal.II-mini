//include/deal.II-translator/numerics/vector_tools_point_value_0.txt
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

#ifndef dealii_vector_tools_point_value_h
#define dealii_vector_tools_point_value_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim, typename Number>
class Point;
template <typename Number>
class Vector;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
} // namespace hp

namespace VectorTools
{
  /**
   * @name  右手边的组装
   *
   */
  //@{

  /**
   * 在点 @p p. 处为点源创建一个右手边矢量
   * 换句话说，它创建一个矢量 $F$ ，以便 $F_i = \int_\Omega
   * \delta(x-p) \varphi_i(x) dx$ 其中 $\varphi_i$ 是 @p dof_handler
   * 描述的形状函数， @p p 是delta函数所在的点。给定的 @p
   * rhs_vector
   * 矢量的事先内容被删除。这个函数是针对标量有限元的情况。
   * 这个函数通常在这两种情况下使用。
   *
   *
   *
   * - 假设你想多次解决同一类问题，用不同的右手边或系数的值，然后每次都在同一点上评估解。你可以通过在每次解题后调用 VectorTools::point_value() 来实现，或者你可以意识到在某一点 $p$ 评估解 $u_h$ ，你可以像这样重新安排操作。
   * @f{align*}{
   *   u_h(p) &= \sum_j U_j \varphi_j(p) = \sum_j U_j F_j
   *     \\   &= U \cdot F
   * @f}
   * 与上面定义的向量。换句话说，只需一个向量-向量的乘积就可以实现点的评估，而向量
   * $F$
   * 可以一次性计算出来，并在每次求解时重复使用，而不必每次都要翻阅网格，找出点
   * $p$ 所在的单元（以及单元中的哪个位置）。
   *
   *
   *
   *
   *
   *
   * - 如果你想计算你要解决的问题的格林函数，这个函数也很有用。这是因为格林函数 $G(x,p)$ 的定义为
   * @f{align*}{
   *   L G(x,p) &= \delta(x-p)
   * @f}
   * 其中 $L$
   * 是你问题的微分算子。离散版本需要计算右手边的向量
   * $F_i = \int_\Omega \varphi_i(x) \delta(x-p)$
   * ，这正是由当前函数计算的向量。
   * 虽然可能与记录这个函数所做的<i>what</i>无关，但值得注意的是，delta函数在现实中并不存在，因此，使用这个函数并不能模拟任何真实情况。这是因为，没有一个真实的物体能够将无限大的力密度集中在领域的一个无限小的部分（相反，所有真实的设备都会将力分散在一个有限的区域内）；也不可能在单个点上测量值（但所有的测量值都会以某种方式在小区域内取平均值）。只有当这个区域非常小，以至于不能被任何网格所解决时，用一个具有相同的整体力或灵敏度的delta函数的方式来建模才有意义。另一方面，用delta函数模拟的情况可能更有成效，那就是点源的电动势；在这种情况下，已知解有一个对数奇点（在2D中）或一个
   * $\frac{1}{r}$ 奇点（在3D中），这两个奇点都不受约束。
   * 在数学上，使用delta函数通常会导致精确的解，而在数值上得到的近似解并不收敛。这是因为，以拉普拉斯方程为例，精确解和数值解之间的误差可以用以下表达式来限定
   * @f{align*}{
   * \| u-u_h \|_{L_2} \le C h \| \nabla u \|_{L_2}
   * @f}
   * 但当在右侧使用delta函数时，项 $\| \nabla u \|_{L_2} =
   * |u|_{H^1}$ 不是有限的。这可以通过使用拉普拉斯方程
   * $-\Delta u = f$ 的解的先验约束看出，该约束指出 $|u|_{H^1}
   * \le \|f\|_{H^{-1}}$  。  当使用delta函数作为右手时，
   * $f(x)=\delta(x-p)$  ，我们需要取delta函数的 $H^{-1}$
   * 规范，然而这不是有限的，因为 $\delta(\cdot-p) \not\in
   * H^{-1}$  。
   * 所有这些的结果是，拉普拉斯方程的精确解在右手边有一个delta函数
   *
   * --即<i>Green's function</i>。
   *
   * 在 $p$
   * 处有一个奇点，这个奇点非常强，不能用有限元解来解决，因此有限元近似值不能以任何通常的准则收敛于精确解。
   * 所有这些对于所有其他二阶偏微分方程在二维或更高维度的情况也是如此。(因为在二维或更高维度上，
   * $H^1$
   * 函数不一定是连续的，因此，delta函数不在对偶空间
   * $H^{-1}$ 中) 。
   *
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const Mapping<dim, spacedim> &   mapping,
                             const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             Vector<double> &                 rhs_vector);

  /**
   * 像前面的函数一样，但对hp-对象而言。
   *
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof_handler,
    const Point<spacedim, double> &             p,
    Vector<double> &                            rhs_vector);

  /**
   * 调用create_point_source_vector()函数，见上文，有一个隐含的默认
   * $Q_1$ 映射对象。
   * 注意，如果你的DoFHandler使用了0以外的任何活动FE索引，那么你需要调用上面的函数，为每个活动FE索引提供一个映射对象。
   *
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             Vector<double> &                 rhs_vector);

  /**
   * 为点 @p p. 处的点源创建右手向量
   * 这个函数的变化是为了解决恰好有二维分量的向量值问题（它也适用于有多于二维分量的问题，在这种情况下只需考虑形状函数的第一个二维分量）。它计算出一个对应于强迫函数的右手边，该强迫函数等于德尔塔函数乘以一个给定的方向。换句话说，它创建了一个向量
   * $F$ ，以便 $F_i = \int_\Omega [\mathbf d \delta(x-p)] \cdot
   * \varphi_i(x) dx$  。这里要注意， $\varphi_i$
   * 是一个矢量值的函数。  $\mathbf d$ 是源项 $\mathbf d
   * \delta(x-p)$ 的给定方向，对应于要传递给该函数的 @p
   * 方向参数。    给定的 @p rhs_vector
   * 矢量的先前内容被删除。
   * 关于delta函数的使用，请参见第一个create_point_source_vector()变量的讨论。
   *
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const Mapping<dim, spacedim> &   mapping,
                             const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             const Point<dim, double> &       direction,
                             Vector<double> &                 rhs_vector);

  /**
   * 和前面的函数一样，但是针对hp-objects。
   *
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof_handler,
    const Point<spacedim, double> &             p,
    const Point<dim, double> &                  direction,
    Vector<double> &                            rhs_vector);

  /**
   * 为矢量值的有限元调用create_point_source_vector()函数，见上文，有一个隐含的默认
   * $Q_1$ 映射对象。
   * 注意，如果你的DoFHandler使用了零以外的任何活动FE索引，那么你需要调用上面的函数，为每个活动FE索引提供一个映射对象。
   *
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             const Point<dim, double> &       direction,
                             Vector<double> &                 rhs_vector);
  // @}

  /**
   * @name  函数的评估和错误
   *
   */
  //@{

  /**
   * 点误差评估。找到包含给定点的第一个单元，并计算一个（可能是矢量值的）有限元函数和一个连续函数（具有与有限元一样多的矢量分量）在该点的差值。
   * 这是一个使用Q1映射的单元格边界的封装函数，用于调用其他point_difference()函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_difference(
    const DoFHandler<dim, spacedim> &                          dof,
    const VectorType &                                         fe_function,
    const Function<spacedim, typename VectorType::value_type> &exact_solution,
    Vector<typename VectorType::value_type> &                  difference,
    const Point<spacedim, double> &                            point);

  /**
   * 点的错误评估。找到包含给定点的第一个单元，并计算一个（可能是矢量值的）有限元函数和一个连续函数（具有与有限元一样多的矢量分量）在该点的差值。
   * 与另一个同名的函数相比，这个函数使用一个任意的映射来评估差异。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_difference(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const VectorType &                                         fe_function,
    const Function<spacedim, typename VectorType::value_type> &exact_solution,
    Vector<typename VectorType::value_type> &                  difference,
    const Point<spacedim, double> &                            point);

  /**
   * 在给定的点 @p 点评估一个由给定的DoFHandler和节点矢量
   * @p fe_function
   * 定义的可能是矢量值的有限元函数，并通过最后一个参数返回这个函数的（矢量）值。
   * 这个函数对点被评估的单元格使用 $Q_1$
   * -映射。如果你需要使用不同的映射来评估（例如，当使用弯曲的边界时），请使用接受映射的point_difference()函数。
   * 这个函数不是特别便宜。这是因为它首先需要找到给定点所在的单元格，然后在参考单元格上找到与给定评估点相匹配的点，然后评估那里的形状函数。你可能不想用这个函数来评估<i>many</i>点的解。对于这种应用，FEFieldFunction类至少提供了一些优化。另一方面，如果你想在同一个点上评估<i>many
   * solutions</i>，你可能想看一下
   * VectorTools::create_point_source_vector() 函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const DoFHandler<dim, spacedim> &        dof,
              const VectorType &                       fe_function,
              const Point<spacedim, double> &          point,
              Vector<typename VectorType::value_type> &value);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const DoFHandler<dim, spacedim> &        dof,
              const VectorType &                       fe_function,
              const Point<spacedim, double> &          point,
              Vector<typename VectorType::value_type> &value);

  /**
   * @note
   * 如果找到的点所在的单元格不是本地拥有的，会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然这只能在一定的数字公差内完成。
   * 因此，对于处于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const DoFHandler<dim, spacedim> &dof,
              const VectorType &               fe_function,
              const Point<spacedim, double> &  point);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于处于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const DoFHandler<dim, spacedim> &dof,
              const VectorType &               fe_function,
              const Point<spacedim, double> &  point);

  /**
   * 评估由给定的DoFHandler和节点向量 @p fe_function 在给定点
   * @p
   * 点定义的可能是矢量值的有限元函数，并通过最后一个参数返回该函数的（矢量）值。
   * 与另一个同名的函数相比，这个函数使用一个任意的映射来评估点值。
   * 这个函数不是特别便宜。这是因为它首先需要找到给定的点在哪个单元格中，然后在参考单元格上找到与给定的评估点相匹配的点，然后评估那里的形状函数。你可能不想用这个函数来评估<i>many</i>点的解。对于这种应用，FEFieldFunction类至少提供了一些优化。另一方面，如果你想在同一个点上评估<i>many
   * solutions</i>，你可能想看一下
   * VectorTools::create_point_source_vector() 函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const Mapping<dim, spacedim> &           mapping,
              const DoFHandler<dim, spacedim> &        dof,
              const VectorType &                       fe_function,
              const Point<spacedim, double> &          point,
              Vector<typename VectorType::value_type> &value);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const hp::MappingCollection<dim, spacedim> &mapping,
              const DoFHandler<dim, spacedim> &           dof,
              const VectorType &                          fe_function,
              const Point<spacedim, double> &             point,
              Vector<typename VectorType::value_type> &   value);

  /**
   * 评估由给定的DoFHandler和节点向量 @p fe_function 在给定点
   * @p point, 定义的标量有限元函数，并返回这个函数的值。
   * 与另一个同名的函数相比，这个函数使用一个任意的映射来评估差异。
   * 这个函数不是特别便宜。这是因为它首先需要找到给定点所在的单元格，然后在参考单元格上找到与给定评估点相匹配的点，然后评估那里的形状函数。你可能不想用这个函数来评估<i>many</i>点的解。对于这种应用，FEFieldFunction类至少提供了一些优化。另一方面，如果你想在同一个点上评估<i>many
   * solutions</i>，你可能想看一下
   * VectorTools::create_point_source_vector() 函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然这只能在一定的数字公差内完成。
   * 因此，对于处于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const Mapping<dim, spacedim> &   mapping,
              const DoFHandler<dim, spacedim> &dof,
              const VectorType &               fe_function,
              const Point<spacedim, double> &  point);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于处于或接近单元边界的点，你可能会在这里或那里得到有限元场的值，这取决于该点是在哪个单元中找到的。如果有限元场是连续的，这并不重要（在相同的公差内）。另一方面，如果使用的有限元是<i>not</i>连续的，那么你会在单元格的边界上或接近边界的地方得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const hp::MappingCollection<dim, spacedim> &mapping,
              const DoFHandler<dim, spacedim> &           dof,
              const VectorType &                          fe_function,
              const Point<spacedim, double> &             point);
  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_point_value_h


