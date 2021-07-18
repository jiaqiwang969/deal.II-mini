//include/deal.II-translator/multigrid/mg_transfer_component_0.txt
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

#ifndef dealii_mg_transfer_component_h
#define dealii_mg_transfer_component_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_base.h>

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
 * 为组件明智的多栅格传输实现矩阵生成。
 *
 *
 * @note
 * MGTransferBlockBase可能是一个更合理的类。但最终还是应该开发一个允许选择多个组件的类。
 *
 *
 */
class MGTransferComponentBase
{
public:
  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;


protected:
  /**
   * 实际上是为每个层次建立延长矩阵。
   * 这个函数只被派生类调用。这些也可以设置成员变量
   * <code>selected_component</code> 和 <code>mg_selected_component</code>
   * 成员变量来限制传递矩阵的某些成分。此外，它们使用
   * <code>target_component</code> and <code>mg_target_component</code>
   * 对组件进行重新排序和分组。
   *
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 所选组件的标志。
   * 转移运算符只作用于这里有<tt>true</tt>条目的组件。如果使用#target_component的重新编号，这指的是<b>renumbered</b>组件。
   *
   */
  ComponentMask component_mask;

  /**
   * 所选组件的标志。
   * 转移操作只作用于这里有<tt>true</tt>项的组件。如果使用#mg_target_component的重新编号，这指的是<b>renumbered</b>组件。
   *
   */
  ComponentMask mg_component_mask;

  /**
   * 如果需要重新编号，细级向量的目标分量。
   *
   */
  std::vector<unsigned int> target_component;

  /**
   * 如果需要对水平向量进行重新编号，则为目标分量。
   *
   */
  std::vector<unsigned int> mg_target_component;

  /**
   * 多级向量的大小。
   *
   */
  mutable std::vector<std::vector<types::global_dof_index>> sizes;

  /**
   * 每个组件的起始索引。
   *
   */
  std::vector<types::global_dof_index> component_start;

  /**
   * 所有层次上的每个组件的起始索引。
   *
   */
  std::vector<std::vector<types::global_dof_index>> mg_component_start;

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
   * 这个变量持有<tt>copy_to/from_mg</tt>-函数的映射。
   * 该数据首先是全局索引，然后是级别索引。
   *
   */
  std::vector<std::vector<std::pair<types::global_dof_index, unsigned int>>>
    copy_to_and_from_indices;

  /**
   * 存储边界指数。这些是限制矩阵中的边界值所需要的。
   *
   */
  std::vector<std::set<types::global_dof_index>> boundary_indices;
};

// TODO:[GK] Update documentation for copy_* functions

// TODO: Use same kind of template argument as MGTransferSelect

/**
 * MGTransferBase接口的实现，用于块状矩阵和简单向量。这个类使用MGTransferComponentBase选择一个单独的组件或将几个组件分组到一个块中。转移操作符本身是为Vector和BlockVector对象实现的。
 * 请参阅MGTransferBase来了解哪一个转移类最适合你的需要。
 *
 *
 */
template <typename number>
class MGTransferSelect : public MGTransferBase<Vector<number>>,
                         private MGTransferComponentBase
{
public:
  /**
   * 没有约束矩阵的构造函数。只在不连续的有限元或没有局部细化的情况下使用这个构造函数。
   *
   */
  MGTransferSelect();

  /**
   * 带有约束矩阵的构造函数。
   *
   */
  MGTransferSelect(const AffineConstraints<double> &constraints);

  /**
   * 解构器。
   *
   */
  virtual ~MGTransferSelect() override = default;

  // TODO: rewrite docs; make sure defaulted args are actually allowed
  /**
   * 实际上是为分组的组件建立延长矩阵。
   * 这个函数是MGTransferComponentBase中相同函数的前端。
   * @arg 选定
   * 要从多级向量复制到全局向量的块的编号。这个数字是指通过<tt>target_component</tt>重新编号。
   * @arg  mg_selected 应建立转移矩阵的块的编号。
   * 如果<tt>mg_target_component</tt>存在，这指的是重新编号的组件。
   * @arg
   * target_component这个参数允许对细级向量中的组件进行分组和重新编号（见
   * DoFRenumbering::component_wise).   @arg
   * mg_target_component这个参数允许对级向量中的组件进行分组和重新编号（见
   * DoFRenumbering::component_wise).
   * 它也会影响<tt>selected</tt>参数的行为  @arg
   * boundary_indices持有每一级的边界指数。
   *
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof,
        unsigned int                     selected,
        unsigned int                     mg_selected,
        const std::vector<unsigned int> &target_component =
          std::vector<unsigned int>(),
        const std::vector<unsigned int> &mg_target_component =
          std::vector<unsigned int>(),
        const std::vector<std::set<types::global_dof_index>> &boundary_indices =
          std::vector<std::set<types::global_dof_index>>());

  /**
   * 改变选定的组件。小心处理!
   *
   */
  void
  select(const unsigned int component,
         const unsigned int mg_component = numbers::invalid_unsigned_int);

  virtual void
  prolongate(const unsigned int    to_level,
             Vector<number> &      dst,
             const Vector<number> &src) const override;

  virtual void
  restrict_and_add(const unsigned int    from_level,
                   Vector<number> &      dst,
                   const Vector<number> &src) const override;

  /**
   * 从全局网格上的向量转移到活动自由度的多级向量。特别是，对于一个全局细化的网格，只有
   * @p dst 中最细的层次被填充为 @p src.
   * 的普通拷贝，所有其他的层次对象都不被触动。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof,
             MGLevelObject<Vector<number>> &  dst,
             const Vector<number2> &          src) const;

  /**
   * 从多级向量转移到普通向量。
   * 将一个多级向量的活动部分的数据复制到向量的相应位置。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &    mg_dof,
               Vector<number2> &                    dst,
               const MGLevelObject<Vector<number>> &src) const;

  /**
   * 将一个多级向量添加到一个正常向量中。
   * 和前面的函数一样工作，但可能不适合连续元素。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim> &    mg_dof,
                   Vector<number2> &                    dst,
                   const MGLevelObject<Vector<number>> &src) const;

  /**
   * 从全局网格上的一个向量转移到活动自由度的多级向量。特别是，对于全局细化的网格，只有
   * @p dst 中最细的层次被填充为 @p src.
   * 的普通拷贝，其他的层次对象都没有被触动。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof,
             MGLevelObject<Vector<number>> &  dst,
             const BlockVector<number2> &     src) const;

  /**
   * 从多级向量转移到普通向量。
   * 将多级向量的活动部分的数据复制到全局BlockVector的相应位置。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &    mg_dof,
               BlockVector<number2> &               dst,
               const MGLevelObject<Vector<number>> &src) const;

  /**
   * 将一个多级向量添加到一个正常向量。
   * 和前面的函数一样工作，但可能不适合连续元素。
   *
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim> &    mg_dof,
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
  do_copy_from_mg(const DoFHandler<dim, spacedim> &    mg_dof,
                  OutVector &                          dst,
                  const MGLevelObject<Vector<number>> &src) const;

  /**
   * 公共函数的实现。
   *
   */
  template <int dim, class OutVector, int spacedim>
  void
  do_copy_from_mg_add(const DoFHandler<dim, spacedim> &    mg_dof,
                      OutVector &                          dst,
                      const MGLevelObject<Vector<number>> &src) const;

  /**
   * copy_to_mg()的实际实现。
   *
   */
  template <int dim, class InVector, int spacedim>
  void
  do_copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof,
                MGLevelObject<Vector<number>> &  dst,
                const InVector &                 src) const;
  /**
   * 全局向量的选定组件。
   *
   */
  unsigned int selected_component;
  /**
   * 多网格内的选定分量。
   *
   */
  unsigned int mg_selected_component;

  /**
   * 细化边上的自由度。对于每一层，索引集表示哪一层的自由度在朝向下一层的细化边上，不包括边界自由度。
   *
   */
  std::vector<IndexSet> interface_dofs;

  /**
   * 全局系统的约束。
   *
   */
public:
  SmartPointer<const AffineConstraints<double>> constraints;
};

 /*@}*/ 

//---------------------------------------------------------------------------
template <typename number>
inline void
MGTransferSelect<number>::select(const unsigned int component,
                                 const unsigned int mg_component)
{
  selected_component = component;
  mg_selected_component =
    (mg_component == numbers::invalid_unsigned_int) ? component : mg_component;
}

DEAL_II_NAMESPACE_CLOSE

#endif


