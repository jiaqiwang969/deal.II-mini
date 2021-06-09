//include/deal.II-translator/base/types_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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

#ifndef dealii_types_h
#define dealii_types_h


#include <deal.II/base/config.h>

#include <cstdint>


DEAL_II_NAMESPACE_OPEN

/**
 * 一个命名空间，我们在其中为deal.II中使用的类型声明别名，以及这些类型的特殊值。
 *
 *
 */
namespace types
{
  /**
   * 用来表示单元格的子域_ID的类型。    参见 @ref GlossSubdomainId "词汇表 "
   * 以了解更多信息。    有一个特殊的值，
   * numbers::invalid_subdomain_id ，用来表示这种类型的无效值。
   *
   */
  using subdomain_id = unsigned int;

  /**
   * 用于顶点的全局索引的类型。
   *
   */
  using global_vertex_index = uint64_t;

  /**
   * 表示与 types::global_vertex_index.
   * 相关的MPI类型的一个标识符。
   *
   */
#define DEAL_II_VERTEX_INDEX_MPI_TYPE MPI_UINT64_T

  /**
   * 用来表示自由度的全局索引的类型。这个类型也用于查询全局自由度的数量*，因为这个数量只是最大的索引加1。
   * 虽然在顺序计算中，32位无符号整数的40亿索引是足够的，但使用（例如）
   * parallel::distributed::Triangulation
   * 类的并行计算会溢出这个数字，因此，当配置为使用64位索引时，deal.II会选择一个更大的整数类型。
   * 数据类型总是与无符号整数类型相对应。    参见 @ref
   * GlobalDoFIndex
   * 页，了解何时应该或不应该使用这种类型的指导。
   *
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
  using global_dof_index = uint64_t;
#else
  using global_dof_index  = unsigned int;
#endif

  /**
   * 表示与 types::global_dof_index.
   * 相关的MPI类型的一个标识符。
   *
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
#  define DEAL_II_DOF_INDEX_MPI_TYPE MPI_UINT64_T
#else
#  define DEAL_II_DOF_INDEX_MPI_TYPE MPI_UNSIGNED
#endif

  /**
   * 用于表示一个单元格的全局索引的类型。这个类型也用于查询三角形中单元格的全局数量*，因为这个数量只是最大的索引加1。
   * 虽然在顺序计算中，32位无符号整数的40亿索引是足够的，但使用（例如）
   * parallel::distributed::Triangulation
   * 类的并行计算可能会溢出这个数字，因此，当配置为使用64位索引时，deal.II选择了一个更大的整数类型。
   * 数据类型总是与无符号整数类型相对应。
   *
   */
#ifdef DEAL_II_WITH_64BIT_INDICES
  using global_cell_index = uint64_t;
#else
  using global_cell_index = unsigned int;
#endif

  /**
   * 用于粗略的单元ID的类型。更多信息请参见术语表中关于 @ref GlossCoarseCellId "粗略单元ID "
   * 的条目。
   *
   */
  using coarse_cell_id = global_cell_index;

  /**
   * 用于表示与边界的每一块相关的边界指标的类型，在描述高维流形的网格中，与每一个单元相关。    有一个特殊的值， numbers::internal_face_boundary_id ，用来表示这种类型的无效值，它被用作位于域内部的面的边界指示器，因此不是任何可寻址的边界组件的一部分。      @see   @ref GlossBoundaryIndicator  "关于边界指示器的词汇条目"
   *
   */
  using boundary_id = unsigned int;

  /**
   * 用于表示与网格的每个对象相关的流形指标的类型。    有一个特殊的值， numbers::flat_manifold_id ，用于表示标准的笛卡尔流形。      @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  using manifold_id = unsigned int;

  /**
   * 用于表示与每个单元相关的材料指标的类型。    有一个特殊的值， numbers::invalid_material_id ，用来表示该类型的无效值。      @see   @ref GlossMaterialId  "关于材料指标的词汇条目"
   *
   */
  using material_id = unsigned int;

  /**
   * 用于表示几何实体类型的类型。
   *
   */
  using geometric_entity_type = std::uint8_t;
} // namespace types

/**
 * 在Epetra中使用的声明类型。
 *
 *
 */
using TrilinosScalar = double;


namespace TrilinosWrappers
{
  namespace types
  {
#ifdef DEAL_II_WITH_64BIT_INDICES
    /**
     * 声明在Trilinos的Epetra包中使用的整数类型。
     *
     */
    using int_type = long long int;
#else
    /**
     * 声明在Trilinos的Epetra包中使用的整数类型。
     *
     */
    using int_type = int;
#endif
  } // namespace types
} // namespace TrilinosWrappers


// this part of the namespace numbers got moved to the bottom types.h file,
// because otherwise we get a circular inclusion of config.h, types.h, and
// numbers.h
namespace numbers
{
  /**
   * 表示可以放入无符号整数的最大数字。这个值在整个库中被广泛使用，作为无效的无符号整数值的标记，例如无效的数组索引，无效的数组大小，等等。
   *
   */
  static const unsigned int invalid_unsigned_int =
    static_cast<unsigned int>(-1);

  /**
   * 表示可以放入一个size_type的最大数字。
   * 这个值在整个库中被用作无效的size_type值的标记，例如无效的数组索引，无效的数组大小，以及类似的。Invalid_size_type等同于invalid_dof_index。
   *
   */
  const types::global_dof_index invalid_size_type =
    static_cast<types::global_dof_index>(-1);

  /**
   * 一个无效的自由度索引值。
   *
   */
  const types::global_dof_index invalid_dof_index =
    static_cast<types::global_dof_index>(-1);

  /**
   * 一个无效的粗略单元ID的值。更多信息请参见术语表中的 @ref GlossCoarseCellId "粗略单元ID "
   * 条目。
   *
   */
  const types::coarse_cell_id invalid_coarse_cell_id =
    static_cast<types::coarse_cell_id>(-1);

  /**
   * 无效的 material_id，我们在几个地方需要它作为默认值。
   * 我们假设所有的 material_id 都在 [0, invalid_material_id)
   * 范围内。
   *
   */
  const types::material_id invalid_material_id =
    static_cast<types::material_id>(-1);

  /**
   * 无效的边界_id，我们需要在几个地方作为默认值。  我们假设所有有效的边界ID都在[0, invalid_boundary_id]的范围内。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  const types::boundary_id invalid_boundary_id =
    static_cast<types::boundary_id>(-1);

  /**
   * 我们为内部面保留的一个边界指示器号码。 我们假设所有有效的边界指标都在[0, internal_face_boundary_id]的范围内。    这是一个内部使用的指标（由库使用），用于区分位于域的边界的面和位于域的内部的面。你不应该在用户代码中尝试将这个边界指标分配给任何东西。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  const types::boundary_id internal_face_boundary_id =
    static_cast<types::boundary_id>(-1);

  /**
   * 我们为默认的平直角坐标流形保留的流形_ID。      @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  const types::manifold_id flat_manifold_id =
    static_cast<types::manifold_id>(-1);

  /**
   * 一个无效的子域id的特殊id。这个值不能作为有效的id使用，但用于例如默认参数，表示不使用的子域id。    参见 @ref GlossSubdomainId "词汇表 "
   * 以了解更多信息。
   *
   */
  const types::subdomain_id invalid_subdomain_id =
    static_cast<types::subdomain_id>(-1);

  /**
   *
   */
  const types::subdomain_id artificial_subdomain_id =
    static_cast<types::subdomain_id>(-2);
} // namespace numbers

DEAL_II_NAMESPACE_CLOSE

#endif


