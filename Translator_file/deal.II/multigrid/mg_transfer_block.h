//include/deal.II-translator/multigrid/mg_transfer_block_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#ifndef dealii_mg_transfer_block_h
#define dealii_mg_transfer_block_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
class DoFHandler;
#endif

/* MGTransferBase在mg_base.h中定义。

* 
*
*/

 /*!@addtogroup mg */ 
 /*@{*/ 

/**
 * 实现MGTransferBlock的矩阵生成。
 * 这是MGTransfer对象的基类，用于多网格只应用于一个或一些块的方程系统，其中一个 @ref
 * GlossBlock 包括由一个基元产生的所有自由度。
 *
 *
 */
class MGTransferBlockBase
{
public:
  /**
   * 没有约束矩阵的构造函数。只在不连续的有限元或没有局部细化的情况下使用此构造函数。
   *
   */
  MGTransferBlockBase();

  /**
   * 带有约束矩阵以及mg_constrained_dofs的构造函数。
   *
   */
  MGTransferBlockBase(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 这个对象所使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 实际建立每个层次的延长矩阵。
   * 这个函数只被派生类调用。这些也可以设置成员变量#selected和其他来限制传送矩阵到某些区块。
   *
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 所选块的标志。
   * 转移运算符只作用于这里有<tt>true</tt>项的块。
   *
   */
  // TODO: rename this to block_mask, in the same way as has already been done
  // in MGTransferComponent, and give it type BlockMask
  std::vector<bool> selected;

  /**
   * 多网格矢量的块数。
   *
   */
  unsigned int n_mg_blocks;

  /**
   * 对于整个块向量的每个块，列出它被映射到多网格向量的哪个块。因为根据#selected，多层次块可能比原始块少，所以有些条目可能是非法的无符号整数。
   *
   */
  // TODO: rename this to mg_block_mask, in the same way as has already been
  // done in MGTransferComponent, and give it type BlockMask
  std::vector<unsigned int> mg_block;

  /**
   * 多级向量的大小。
   *
   */
  mutable std::vector<std::vector<types::global_dof_index>> sizes;

  /**
   * 每个区块的起始索引。
   *
   */
  std::vector<types::global_dof_index> block_start;

  /**
   * 所有级别上每个块的起始索引。
   *
   */
  std::vector<std::vector<types::global_dof_index>> mg_block_start;

  /**
   * 先调用build()函数。
   *
   */
  DeclException0(ExcMatricesNotBuilt);

private:
  std::vector<std::shared_ptr<BlockSparsityPattern>> prolongation_sparsities;

protected:
  /**
   * 实际的延长矩阵。列指数属于母单元的道夫指数，即粗级别。而行指数属于子单元，即细级别。
   *
   */
  std::vector<std::shared_ptr<BlockSparseMatrix<double>>> prolongation_matrices;

  /**
   * 用于<tt>copy_to/from_mg</tt>-函数的映射。进入这个向量的索引是（按这个顺序）：全局块号，级别号。数据首先是块内的全局索引，然后是块内的级别索引。
   *
   */
  std::vector<std::vector<std::vector<std::pair<unsigned int, unsigned int>>>>
    copy_indices;

  /**
   * 级系统的mg_constrained_dofs。
   *
   */

  SmartPointer<const MGConstrainedDoFs, MGTransferBlockBase>
    mg_constrained_dofs;
};

/**
 * 实现块矩阵和块向量的MGTransferBase接口。
 * @warning
 * 该类处于未测试状态。如果你使用它时遇到了问题，请联系Guido
 * Kanschat。
 * 除了MGTransferPrebuilt的功能外，操作可以被限制在向量的某些块上。
 * 如果选择了限制模式，在转移例程中使用的块向量可能只有选择字段中的
 * @p trues 那么多块。
 * 请参阅MGTransferBase，了解哪种转移类最适合你的需要。
 *
 *
 */
template <typename number>
class MGTransferBlock : public MGTransferBase<BlockVector<number>>,
                        private MGTransferBlockBase
{
public:
  /**
   * 默认构造函数。
   *
   */
  MGTransferBlock();

  /**
   * 解构器。
   *
   */
  virtual ~MGTransferBlock() override;

  /**
   * 如果要对块的限制进行不同的加权，则初始化额外的#因素和#内存。
   *
   */
  void
  initialize(const std::vector<number> &   factors,
             VectorMemory<Vector<number>> &memory);

  /**
   * 为每一层建立延长矩阵。
   * 这个函数是MGTransferBlockBase中同一函数的前端。
   *
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof_handler,
        const std::vector<bool> &        selected);

  virtual void
  prolongate(const unsigned int         to_level,
             BlockVector<number> &      dst,
             const BlockVector<number> &src) const override;

  virtual void
  restrict_and_add(const unsigned int         from_level,
                   BlockVector<number> &      dst,
                   const BlockVector<number> &src) const override;

  /**
   * 从全局网格上的向量转移到活动自由度的多级向量。特别是，对于全局细化的网格，只有
   * @p dst 中最细的层次被填充为 @p src.
   * 的普通拷贝，其他的层次对象都没有被触动。
   * 对不连续元素的操作如下：在一个活动的网格单元上，全局向量条目被简单地复制到水平向量的相应条目上。然后，这些值被限制到最粗的级别。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &   dof_handler,
             MGLevelObject<BlockVector<number>> &dst,
             const BlockVector<number2> &        src) const;

  /**
   * 从多级向量转移到普通向量。
   * 将数据从多级向量的活动部分复制到全局向量的相应位置。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &         dof_handler,
               BlockVector<number2> &                    dst,
               const MGLevelObject<BlockVector<number>> &src) const;

  /**
   * 将一个多级向量添加到一个法向量中。
   * 和前面的函数一样工作，但可能不适合连续元素。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim> &         dof_handler,
                   BlockVector<number2> &                    dst,
                   const MGLevelObject<BlockVector<number>> &src) const;

  using MGTransferBlockBase::memory_consumption;

private:
  /**
   * 每个区块的可选乘法因子。需要初始化#memory。
   *
   */
  std::vector<number> factors;

  /**
   * 如果需要使用#因子的额外乘法，则需要内存池。
   *
   */
  SmartPointer<VectorMemory<Vector<number>>, MGTransferBlock<number>> memory;
};


// TODO:[GK] Update documentation for copy_* functions

/**
 * MGTransferBase接口的实现，用于块状矩阵和简单向量。该类使用MGTransferBlockBase选择单一块。对Vector对象实现了网格间转移操作，对Vector和BlockVector实现了常规和多网格向量间的复制功能。
 * 请参阅MGTransferBase来了解哪一个转移类最适合你的需要。
 *
 *
 */
template <typename number>
class MGTransferBlockSelect : public MGTransferBase<Vector<number>>,
                              private MGTransferBlockBase
{
public:
  /**
   * 没有约束矩阵的构造函数。只在不连续的有限元或没有局部细化的情况下使用这个构造函数。
   *
   */
  MGTransferBlockSelect();

  /**
   * 带有约束矩阵以及mg_constrained_dofs的构造函数。
   *
   */
  MGTransferBlockSelect(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 解构器。
   *
   */
  virtual ~MGTransferBlockSelect() override = default;

  /**
   * 实际上是为分组的块建立延长矩阵。
   * 这个函数是MGTransferBlockBase中相同函数的前端。      @param
   * dof_handler 要使用的DoFHandler。    @param  selected
   * 应建立转移矩阵的块的编号。
   *
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof_handler, unsigned int selected);

  /**
   * 改变选定的区块。小心处理!
   *
   */
  void
  select(const unsigned int block);

  virtual void
  prolongate(const unsigned int    to_level,
             Vector<number> &      dst,
             const Vector<number> &src) const override;

  virtual void
  restrict_and_add(const unsigned int    from_level,
                   Vector<number> &      dst,
                   const Vector<number> &src) const override;

  /**
   * 将单一块从全局网格上的矢量转移到活动自由度的多级矢量上。特别是，对于一个全局细化的网格，只有
   * @p dst 中最细的层次被填充为 @p src.
   * 的普通拷贝，其他的层次对象都没有被触动。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<Vector<number>> &  dst,
             const Vector<number2> &          src) const;

  /**
   * 从多级向量转移到普通向量。
   * 将一个多级向量的活动部分的数据复制到向量的相应位置。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &    dof_handler,
               Vector<number2> &                    dst,
               const MGLevelObject<Vector<number>> &src) const;

  /**
   * 将一个多级向量添加到一个正常向量。
   * 和前面的函数一样工作，但可能不适合连续元素。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim> &    dof_handler,
                   Vector<number2> &                    dst,
                   const MGLevelObject<Vector<number>> &src) const;

  /**
   * 从全局网格上的一个向量转移一个块到多级向量。
   * 只有所选块的活动自由度的值被转移。特别是，对于一个全局细化的网格，只有
   * @p dst 中最细的层次被填充为 @p src.
   * 的普通拷贝，所有其他的层次对象都不被触动。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<Vector<number>> &  dst,
             const BlockVector<number2> &     src) const;

  /**
   * 从多级向量转移到普通向量。
   * 将多级向量的活动部分的数据复制到全局BlockVector的相应位置。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &    dof_handler,
               BlockVector<number2> &               dst,
               const MGLevelObject<Vector<number>> &src) const;

  /**
   * 将一个多级向量添加到一个正常向量。
   * 和前面的函数一样工作，但可能不适合连续元素。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim> &    dof_handler,
                   BlockVector<number2> &               dst,
                   const MGLevelObject<Vector<number>> &src) const;

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 公共函数的实现。
   *
   */
  template <int dim, class OutVector, int spacedim>
  void
  do_copy_from_mg(const DoFHandler<dim, spacedim> &    dof_handler,
                  OutVector &                          dst,
                  const MGLevelObject<Vector<number>> &src,
                  const unsigned int                   offset) const;

  /**
   * 公共函数的实现。
   *
   */
  template <int dim, class OutVector, int spacedim>
  void
  do_copy_from_mg_add(const DoFHandler<dim, spacedim> &    dof_handler,
                      OutVector &                          dst,
                      const MGLevelObject<Vector<number>> &src,
                      const unsigned int                   offset) const;

  /**
   * copy_to_mg()的实际实现。
   *
   */
  template <int dim, class InVector, int spacedim>
  void
  do_copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
                MGLevelObject<Vector<number>> &  dst,
                const InVector &                 src,
                const unsigned int               offset) const;
  /**
   * 选定的块。
   *
   */
  unsigned int selected_block;
};

 /*@}*/ 

//------------------------- inline function definition ------------------------
template <typename number>
inline void
MGTransferBlockSelect<number>::select(const unsigned int block)
{
  selected_block = block;
}

DEAL_II_NAMESPACE_CLOSE

#endif


