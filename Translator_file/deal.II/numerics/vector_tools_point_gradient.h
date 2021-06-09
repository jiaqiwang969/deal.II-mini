//include/deal.II-translator/numerics/vector_tools_point_gradient_0.txt
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

#ifndef dealii_vector_tools_point_gradient_h
#define dealii_vector_tools_point_gradient_h


#include <deal.II/base/config.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim, typename Number>
class Point;
template <int rank_, int dim, typename Number>
class Tensor;
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
   * @name  函数的评估和错误
   *
   */
  //@{

  /**
   * 评估由给定DoFHandler和节点向量在给定点定义的可能是矢量值的有限元函数，并通过最后一个参数返回该函数的（矢量）梯度。
   * 这是一个使用单元格边界的Q1映射来调用其他point_gradient()函数的包装函数。
   * 这个函数不是特别便宜。这是因为它首先需要找到给定点所在的单元格，然后在参考单元格上找到与给定评估点相匹配的点，然后在那里评估形状函数。你可能不想用这个函数来评估<i>many</i>点的解。对于这种应用，FEFieldFunction类至少提供了一些优化。另一方面，如果你想在同一个点上评估<i>many
   * solutions</i>，你可能想看看
   * VectorTools::create_point_source_vector() 函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const DoFHandler<dim, spacedim> &dof,
    const VectorType &               fe_function,
    const Point<spacedim, double> &  point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果找到的点所在的单元格不是本地拥有的，会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于处于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const DoFHandler<dim, spacedim> &dof,
    const VectorType &               fe_function,
    const Point<spacedim, double> &  point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * 评估由给定DoFHandler和节点向量在给定点定义的标量有限元函数，并返回这个函数的梯度。
   * 与另一个同名的函数相比，这是一个使用Q1映射的单元格的包装函数。
   * 这个函数不是特别便宜。这是因为它首先需要找到一个给定的点在哪个单元格中，然后在参考单元格上找到与给定的评估点相匹配的点，然后评估那里的形状函数。你可能不想用这个函数来评估<i>many</i>点的解。对于这种应用，FEFieldFunction类至少提供了一些优化。另一方面，如果你想在同一个点上评估<i>many
   * solutions</i>，你可能想看一下
   * VectorTools::create_point_source_vector() 函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const DoFHandler<dim, spacedim> &dof,
                 const VectorType &               fe_function,
                 const Point<spacedim, double> &  point);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果找到的点所在的单元格不是本地拥有的，会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const DoFHandler<dim, spacedim> &dof,
                 const VectorType &               fe_function,
                 const Point<spacedim, double> &  point);

  /**
   * 评估由给定DoFHandler和节点向量在给定点定义的可能是矢量值的有限元函数，并通过最后一个参数返回这个函数的梯度。
   * 与另一个同名的函数相比，这个函数使用一个任意的映射进行评估。
   * 这个函数不是特别便宜。这是因为它首先需要找到给定点所在的单元格，然后在参考单元格上找到与给定评价点相匹配的点，然后在那里评价形状函数。你可能不想用这个函数来评估<i>many</i>点的解。对于这种应用，FEFieldFunction类至少提供了一些优化。另一方面，如果你想在同一个点上评估<i>many
   * solutions</i>，你可能想看一下
   * VectorTools::create_point_source_vector() 函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const VectorType &               fe_function,
    const Point<spacedim, double> &  point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果找到的点所在的单元格不是本地拥有的，会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const VectorType &                          fe_function,
    const Point<spacedim, double> &             point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * 评估由给定DoFHandler和节点向量在给定点定义的标量有限元函数，并返回这个函数的梯度。
   * 与另一个同名的函数相比，这个函数使用一个任意的映射进行评估。
   * 这个函数不是特别便宜。这是因为它首先需要找到给定点所在的单元格，然后在参考单元格上找到与给定评价点相匹配的点，然后在那里评价形状函数。你可能不想用这个函数来评估<i>many</i>点的解。对于这种应用，FEFieldFunction类至少提供了一些优化。另一方面，如果你想在同一个点上评估<i>many
   * solutions</i>，你可能想看一下
   * VectorTools::create_point_source_vector() 函数。
   * @note
   * 如果发现点所在的单元格不是本地拥有的，就会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const Mapping<dim, spacedim> &   mapping,
                 const DoFHandler<dim, spacedim> &dof,
                 const VectorType &               fe_function,
                 const Point<spacedim, double> &  point);

  /**
   * 与上述hp的情况相同。
   * @note
   * 如果找到的点所在的单元格不是本地拥有的，会抛出一个
   * VectorTools::ExcPointNotAvailableHere 类型的异常。
   * @note
   * 这个函数需要找到一个点所在的单元格，当然，这只能在一定的数字公差内完成。
   * 因此，对于位于或接近单元边界的点，你可能会在这里或那里得到有限元场的梯度，这取决于该点是在哪个单元中找到的。因为对于大多数元素来说，梯度从一个单元到另一个单元都是不连续的，所以对于单元边界上的点或接近单元边界的点，你会得到不可预测的值，正如人们在试图评估不连续函数的点值时所期望的那样。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const hp::MappingCollection<dim, spacedim> &mapping,
                 const DoFHandler<dim, spacedim> &           dof,
                 const VectorType &                          fe_function,
                 const Point<spacedim, double> &             point);

  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_point_gradient_h


