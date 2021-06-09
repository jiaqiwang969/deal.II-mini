//include/deal.II-translator/grid/cell_id_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_cell_id_h
#define dealii_cell_id_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#ifdef DEAL_II_WITH_P4EST
#  include <deal.II/distributed/p4est_wrappers.h>
#endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int, int>
class Triangulation;
#endif

/**
 * 一个表示三角结构中单元格的唯一ID的类。它由`cell->id()`返回（即，
 * CellAccessor::id()) ，其中`cell`被认为是一个单元格迭代器。
 * 这个类存储了一个单元的后裔的粗略单元的索引（或者更具体地说， @ref GlossCoarseCellId "粗略单元ID "
 * 的条目），以及如何从该粗略单元到达该单元的信息（即，当从一个单元移动到其子单元时，在三角结构的每一层采取哪个子单元索引）。关于这个类的重要的一点是，当前类的对象在三角剖分中唯一地识别一个单元，它甚至在类型为
 * parallel::distributed::Triangulation
 * 的对象的背景下也是如此，因为网格的局部部分可能不存储所有单元。例如，在一个处理器上为一个幽灵单元计算的CellId与在实际拥有该单元的处理器上为同一单元计算的CellId完全相同，尽管指向该单元<i>within
 * the triangulation stored on each of the
 * processors</i>的迭代器的级别和索引可能（而且一般会）不同。换句话说，CellId提供了一个工具，用它可以全局地、唯一地识别平行三角形中的单元，从而使处理器之间交换与单个单元相关的数据成为可能。
 *
 *
 * @note
 * 这些数据在内部如何表示并不重要（也没有故意暴露）。
 *
 *
 */
class CellId
{
public:
  /**
   * 一个用于以紧凑和快速的方式编码CellId数据的类型（例如，用于MPI传输到其他进程）。请注意，它限制了可以传输的子节点数量，在三维中为20个，在二维中为30个（使用2倍的32位进行存储），这个限制与p4est使用的限制相同。
   *
   */
  using binary_type = std::array<unsigned int, 4>;

  /**
   * 用给定的 @p coarse_cell_id
   * 和子指数的向量构造一个CellId对象。  @p child_indices
   * 的解释与同名的成员变量相同，即每个条目表示从一个细化级别到下一个细化级别挑选哪个子单元，从粗略的单元开始，直到我们到达当前对象代表的单元。因此，每个条目应该是一个介于0和当前空间维度中单元格的子代数之间的数字（即，
   * GeometryInfo<dim>::max_children_per_cell). 。
   *
   */
  CellId(const types::coarse_cell_id      coarse_cell_id,
         const std::vector<std::uint8_t> &child_indices);

  /**
   * 用给定的 @p coarse_cell_id 和 @p child_indices.
   * 中提供的子节点数组构造一个CellId对象 @p child_indices
   * ，其解释与同名的成员变量相同，即每个条目表示从一个细化层到下一个细化层要挑选哪个子节点，从粗大的单元开始，直到我们到达当前对象所代表的单元。因此，每个条目应该是一个介于0和当前空间维度中单元格的子数之间的数字（即
   * GeometryInfo<dim>::max_children_per_cell).  数组 @p child_indices
   * 必须至少有 @p n_child_indices 个有效条目。
   *
   */
  CellId(const types::coarse_cell_id coarse_cell_id,
         const unsigned int          n_child_indices,
         const std::uint8_t *        child_indices);

  /**
   * 用给定的二进制表示法构建一个CellId对象，该对象之前是由
   * CellId::to_binary. 构建的。
   *
   */
  CellId(const binary_type &binary_representation);

  /**
   * 从一个与to_string()产生的格式相同的字符串创建一个CellId。
   *
   */
  explicit CellId(const std::string &string_representation);

  /**
   * 构建一个无效的CellId。
   *
   */
  CellId();

  /**
   * 返回这个CellId的人可读的字符串表示。
   * 该函数返回的字符串仅由ASCII字符组成，例如，看起来像这样。`"0_3:006"`.
   * 它可以*被人类解释为："这个单元来自于第四个粗略网格单元，生活在细化水平3，从粗略网格单元到其子代和孙代的路径由006给出"。
   * 但这并不意味着*可以用任何有意义的方式来解释。它只是一种表示当前对象的内部状态的方法，只使用可打印范围内的ASCII字符。
   *
   */
  std::string
  to_string() const;

  /**
   * 返回这个CellId的一个紧凑而快速的二进制表示。
   *
   */
  template <int dim>
  binary_type
  to_binary() const;

  /**
   * 返回一个到此CellId所代表的单元格的cell_iterator。
   * @deprecated  使用 Triangulation::create_cell_iterator() 代替。
   *
   */
  template <int dim, int spacedim>
  DEAL_II_DEPRECATED typename Triangulation<dim, spacedim>::cell_iterator
  to_cell(const Triangulation<dim, spacedim> &tria) const;

  /**
   * 比较两个CellId对象是否相等。
   *
   */
  bool
  operator==(const CellId &other) const;

  /**
   * 比较两个CellIds的不等式。
   *
   */
  bool
  operator!=(const CellId &other) const;

  /**
   * 比较两个CellIds的排序。这个排序的细节是未指定的，只是该操作提供了所有单元格之间的总排序。
   *
   */
  bool
  operator<(const CellId &other) const;

  /**
   * 确定这个单元格ID是否是输入单元格ID的直接父级。
   *
   */
  bool
  is_parent_of(const CellId &other) const;

  /**
   * 确定此单元格ID是否为输入单元格ID的祖先。
   *
   */
  bool
  is_ancestor_of(const CellId &other) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * 返回粗体单元的id。
   *
   */
  types::coarse_cell_id
  get_coarse_cell_id() const;

  /**
   * 返回一个整数的只读容器，表示从一个细化级别到下一个细化级别，从粗略单元开始，直到我们到达当前对象所代表的单元，要挑选哪个孩子。
   * 这个容器中的元素数对应于当前单元的（第1级）。
   *
   */
  ArrayView<const std::uint8_t>
  get_child_indices() const;

private:
  /**
   * 当前对象所代表的单元格位于其树中的粗单元格的编号。
   *
   */
  types::coarse_cell_id coarse_cell_id;

  /**
   * 存储在child_indices数组中的子索引的数量。这相当于当前单元格的（第1级）。
   *
   */
  unsigned int n_child_indices;

  /**
   * 一个整数数组，表示从一个细化级别到下一个细化级别，从粗略的单元开始，直到我们到达当前对象所代表的单元，要挑选哪个孩子。
   * 只有最初的n_child_indices条目被使用，但我们使用静态分配的数组，而不是n_child_indices大小的向量，以加快该对象的创建速度。如果给定的尺寸成为一种限制，那么数组可以被扩展。
   *
   */
#ifdef DEAL_II_WITH_P4EST
  std::array<std::uint8_t, internal::p4est::functions<2>::max_level>
    child_indices;
#else
  std::array<std::uint8_t, 30> child_indices;
#endif

  friend std::istream &
  operator>>(std::istream &is, CellId &cid);
  friend std::ostream &
  operator<<(std::ostream &os, const CellId &cid);
};



/**
 * 将一个CellId对象写入一个流中。
 *
 *
 */
inline std::ostream &
operator<<(std::ostream &os, const CellId &cid)
{
  os << cid.coarse_cell_id << '_' << cid.n_child_indices << ':';
  for (unsigned int i = 0; i < cid.n_child_indices; ++i)
    // write the child indices. because they are between 0 and 2^dim-1, they all
    // just have one digit, so we could write them as one character
    // objects. it's probably clearer to write them as one-digit characters
    // starting at '0'
    os << static_cast<unsigned char>('0' + cid.child_indices[i]);
  return os;
}



/**
 * 序列化功能
 *
 *
 */
template <class Archive>
void
CellId::serialize(Archive &ar, const unsigned int  /*version*/ )
{
  ar &coarse_cell_id;
  ar &n_child_indices;
  ar &child_indices;
}

/**
 * 从一个流中读取一个CellId对象。
 *
 *
 */
inline std::istream &
operator>>(std::istream &is, CellId &cid)
{
  unsigned int cellid;
  is >> cellid;
  if (is.eof())
    return is;

  cid.coarse_cell_id = cellid;
  char dummy;
  is >> dummy;
  Assert(dummy == '_', ExcMessage("invalid CellId"));
  is >> cid.n_child_indices;
  is >> dummy;
  Assert(dummy == ':', ExcMessage("invalid CellId"));

  unsigned char value;
  for (unsigned int i = 0; i < cid.n_child_indices; ++i)
    {
      // read the one-digit child index (as an integer number) and
      // convert it back into unsigned integer type
      is >> value;
      cid.child_indices[i] = value - '0';
    }
  return is;
}



inline bool
CellId::operator==(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;
  if (n_child_indices != other.n_child_indices)
    return false;

  for (unsigned int i = 0; i < n_child_indices; ++i)
    if (child_indices[i] != other.child_indices[i])
      return false;

  return true;
}



inline bool
CellId::operator!=(const CellId &other) const
{
  return !(*this == other);
}



inline bool
CellId::operator<(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return this->coarse_cell_id < other.coarse_cell_id;

  unsigned int idx = 0;
  while (idx < n_child_indices)
    {
      if (idx >= other.n_child_indices)
        return false;

      if (child_indices[idx] != other.child_indices[idx])
        return child_indices[idx] < other.child_indices[idx];

      ++idx;
    }

  if (n_child_indices == other.n_child_indices)
    return false;
  return true; // other.id is longer
}



inline bool
CellId::is_parent_of(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;

  if (n_child_indices + 1 != other.n_child_indices)
    return false;

  for (unsigned int idx = 0; idx < n_child_indices; ++idx)
    if (child_indices[idx] != other.child_indices[idx])
      return false;

  return true; // other.id is longer
}



inline bool
CellId::is_ancestor_of(const CellId &other) const
{
  if (this->coarse_cell_id != other.coarse_cell_id)
    return false;

  if (n_child_indices >= other.n_child_indices)
    return false;

  for (unsigned int idx = 0; idx < n_child_indices; ++idx)
    if (child_indices[idx] != other.child_indices[idx])
      return false;

  return true; // other.id is longer
}



inline types::coarse_cell_id
CellId::get_coarse_cell_id() const
{
  return coarse_cell_id;
}



inline ArrayView<const std::uint8_t>
CellId::get_child_indices() const
{
  return {child_indices.data(), n_child_indices};
}


DEAL_II_NAMESPACE_CLOSE

#endif


