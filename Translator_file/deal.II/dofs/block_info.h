//include/deal.II-translator/dofs/block_info_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
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

#ifndef dealii_block_info_h
#define dealii_block_info_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/block_indices.h>

#include <iomanip>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class DoFHandler;
#endif


/**
 *
 * @brief A small class collecting the different BlockIndices involved in
 * 全局、多层次和局部计算。
 * 一旦DoFHandler被初始化为FESystem，一个BlockInfo类型的数据对象（通过
 * DoFHandler::block_info()
 * 访问）就被填充，它反映了自由度的块结构。
 * BlockInfo由几个BlockIndices对象组成。成员global()反映了活动单元层面上的系统块结构，通常被称为全局系统。一旦
 * DoFHandler::distribute_dofs() 被调用，global()中的函数
 * BlockIndices::block_size() 将返回每个块的正确尺寸。在
 * DoFRenumbering::block_wise(), 之后， BlockIndices::block_start()
 * 将返回每个块的起始索引。
 * 当使用带有级别的DoFHandler时，每个级别都会自动生成相同的结构。级别块可以通过level()访问。
 * 最后，还有local()
 * BlockIndices，它描述了单个单元上的块结构。例如，这被
 * MeshWorker::Assembler::MatrixLocalBlocksToGlobalBlocks.
 * 所使用。本地索引不是自动填充的，因为它们改变了依赖BlockInfo的
 * MeshWorker::Assembler
 * 类的行为。它们必须通过initialize_local()手工初始化。
 * <h3>Usage</h3>
 * 这个对象最常见的用法是初始化向量，如下面的代码。
 *
 *
 * @code
 * DoFHandler<dim> dof_handler(triangulation);
 * dof_handler.distribute_dofs(fe_system);
 * dof_handler.distribute_mg_dofs(fe_system);
 * DoFRenumbering::block_wise(dof_handler);
 *
 * BlockVector<double> solution(dof_handler.block_info().global());
 *
 * MGLevelObject<BlockVector<double> > mg_vector(0, triangulation.n_levels()-1);
 * for (unsigned int i = 0; i < triangulation.n_levels(); ++i)
 * {
 *   mg_vector[i].reinit(dof_handler.block_info().level(i));
 * }
 * @endcode
 *
 * 在这个例子中，<tt>solution</tt>获得了在DoFHandler上表示一个有限元函数所需的块结构。同样地，<tt>mg_vector</tt>的所有层次都会有该层次所需的块结构。
 * @todo  扩展函数local()和renumber()以允许hp-capablilites。
 *
 *
 * @ingroup dofs
 *
 *
 */
class BlockInfo : public Subscriptor
{
public:
  /**
   * @brief Fill the object with values describing block structure of the
   * DoFHandler。
   * 默认情况下，这个函数将尝试初始化任何可能的东西。如果在DoFHandler参数中已经分配了活动的Dofs，它们的BlockIndices将被生成。对水平方向也是如此。
   * 这个默认行为可以被两个参数覆盖，这两个参数可以关闭活动道夫或水平道夫。
   * 这个函数也将清除local()指数。
   *
   */
  template <int dim, int spacedim>
  void
  initialize(const DoFHandler<dim, spacedim> &,
             bool levels_only = false,
             bool active_only = false);

  /**
   * @brief Initialize block structure on cells and compute renumbering
   * 在单元格道夫和块单元格道夫之间。
   *
   */
  template <int dim, int spacedim>
  void
  initialize_local(const DoFHandler<dim, spacedim> &);

  /**
   * 访问全局系统的BlockIndices结构。
   *
   */
  const BlockIndices &
  global() const;

  /**
   * 访问一个单元上的本地系统的BlockIndices。
   *
   */
  const BlockIndices &
  local() const;

  /**
   * 访问多级层次结构中某一级的BlockIndices结构。
   *
   */
  const BlockIndices &
  level(unsigned int level) const;

  /**
   * 返回局部重新编号后的索引。
   * 这个函数的输入是一个介于零和每个单元的道夫数之间的索引，按本地块的顺序编号，即首先是第一个系统块的所有索引，然后是第二个块的所有索引，依此类推。然后该函数以DoFAccessor的标准本地编号输出该索引。
   *
   */
  types::global_dof_index
  renumber(const unsigned int i) const;

  /**
   * 基础元素的数量。
   *
   */
  unsigned int
  n_base_elements() const;

  /**
   * 返回这个索引的基数元素。
   *
   */
  unsigned int
  base_element(const unsigned int i) const;

  /**
   * 将块结构的摘要写到流中。
   *
   */
  template <class OS>
  void
  print(OS &stream) const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int  /*version*/ );

private:
  /**
   * @brief The block structure of the global system.
   *
   */
  BlockIndices bi_global;
  /**
   * @brief The multilevel block structure.
   *
   */
  std::vector<BlockIndices> levels;

  /**
   * @brief The block structure of the cell systems.
   *
   */
  BlockIndices bi_local;

  /**
   * 与每个区块相关的基本元素。
   *
   */
  std::vector<unsigned int> base_elements;

  /**
   * 一个向量，包含从单元上的标准自由度顺序到组件明智顺序的重新编号。由initialize()填充。
   *
   */
  std::vector<types::global_dof_index> local_renumbering;
};



//----------------------------------------------------------------------//

inline const BlockIndices &
BlockInfo::global() const
{
  return bi_global;
}


inline const BlockIndices &
BlockInfo::local() const
{
  return bi_local;
}


inline const BlockIndices &
BlockInfo::level(const unsigned int l) const
{
  AssertIndexRange(l, levels.size());
  return levels[l];
}


inline types::global_dof_index
BlockInfo::renumber(const unsigned int i) const
{
  AssertIndexRange(i, static_cast<unsigned int>(local_renumbering.size()));
  return local_renumbering[i];
}


inline unsigned int
BlockInfo::base_element(const unsigned int i) const
{
  AssertIndexRange(i, base_elements.size());

  return base_elements[i];
}


inline unsigned int
BlockInfo::n_base_elements() const
{
  return base_elements.size();
}



template <class OS>
inline void
BlockInfo::print(OS &os) const
{
  os << "global   dofs " << std::setw(5) << global().total_size() << " blocks";
  for (unsigned int i = 0; i < global().size(); ++i)
    os << ' ' << std::setw(5) << global().block_size(i);
  os << std::endl;

  if (local().size() == 0)
    {
      os << "local dofs not initialized" << std::endl;
    }
  else
    {
      os << "local    dofs " << std::setw(5) << local().total_size()
         << " blocks";
      for (unsigned int i = 0; i < local().size(); ++i)
        os << ' ' << std::setw(5) << local().block_size(i);
      os << std::endl;
    }

  for (unsigned int l = 0; l < levels.size(); ++l)
    {
      os << "level " << std::setw(2) << l << " dofs " << std::setw(5)
         << level(l).total_size() << " blocks";
      for (unsigned int i = 0; i < level(l).size(); ++i)
        os << ' ' << std::setw(5) << level(l).block_size(i);
      os << std::endl;
    }
}


inline std::size_t
BlockInfo::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(bi_global) +
          MemoryConsumption::memory_consumption(levels) +
          MemoryConsumption::memory_consumption(bi_local) +
          MemoryConsumption::memory_consumption(base_elements));
}


template <class Archive>
void
BlockInfo::serialize(Archive &ar, const unsigned int  /*version*/ )
{
  ar &bi_global;
  ar &levels;
  ar &bi_local;
  ar &base_elements;
  ar &local_renumbering;
}


DEAL_II_NAMESPACE_CLOSE

#endif


