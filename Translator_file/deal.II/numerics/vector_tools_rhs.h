//include/deal.II-translator/numerics/vector_tools_rhs_0.txt
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

#ifndef dealii_vector_tools_rhs_h
#define dealii_vector_tools_rhs_h

#include <deal.II/base/config.h>

#include <set>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim>
class Quadrature;
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
   * @name  右手边的组装
   *
   */
  //@{

  /**
   * 创建一个右手边的向量。先前给定的 @p rhs_vector
   * 向量的内容被删除。
   * 更多信息请参见该命名空间的一般文档。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_right_hand_side(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const Quadrature<dim> &                                    q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const AffineConstraints<typename VectorType::value_type> & constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * 调用create_right_hand_side()函数，见上文，<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_right_hand_side(
    const DoFHandler<dim, spacedim> &                          dof,
    const Quadrature<dim> &                                    q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const AffineConstraints<typename VectorType::value_type> & constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * 和前面的一组函数一样，但是针对hp-objects。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_right_hand_side(
    const hp::MappingCollection<dim, spacedim> &               mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const hp::QCollection<dim> &                               q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const AffineConstraints<typename VectorType::value_type> & constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * 像以前的一组函数，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_right_hand_side(
    const DoFHandler<dim, spacedim> &                          dof,
    const hp::QCollection<dim> &                               q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const AffineConstraints<typename VectorType::value_type> & constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * 从边界力创建一个右手边的向量。之前给定的 @p rhs_vector 矢量的内容被删除。    更多信息请参见该命名空间的一般文档。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_boundary_right_hand_side(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const Quadrature<dim - 1> &                                q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const std::set<types::boundary_id> &                       boundary_ids =
      std::set<types::boundary_id>());

  /**
   * 调用create_boundary_right_hand_side()函数，见上文，<tt>mapping=MappingQGeneric  @<dim@>(1)</tt>.   @see   @ref GlossBoundaryIndicator  "关于边界指标的术语条目"
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_boundary_right_hand_side(
    const DoFHandler<dim, spacedim> &                          dof,
    const Quadrature<dim - 1> &                                q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const std::set<types::boundary_id> &                       boundary_ids =
      std::set<types::boundary_id>());

  /**
   * 与上面的函数集相同，但用于hp-objects。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_boundary_right_hand_side(
    const hp::MappingCollection<dim, spacedim> &               mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const hp::QCollection<dim - 1> &                           q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const std::set<types::boundary_id> &                       boundary_ids =
      std::set<types::boundary_id>());

  /**
   * 调用create_boundary_right_hand_side()函数，见上文，用一个Q1映射作为集合。因此，这个函数只有在使用中的唯一活跃的FE指标为0时才会起作用。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  create_boundary_right_hand_side(
    const DoFHandler<dim, spacedim> &                          dof,
    const hp::QCollection<dim - 1> &                           q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType &                                               rhs_vector,
    const std::set<types::boundary_id> &                       boundary_ids =
      std::set<types::boundary_id>());
  // @}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_rhs_h


