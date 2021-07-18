//include/deal.II-translator/multigrid/mg_transfer_matrix_free_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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

#ifndef dealii_mg_transfer_matrix_free_h
#define dealii_mg_transfer_matrix_free_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_internal.h>


DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup mg */ 
 /*@{*/ 

/**
 * MGTransferBase接口的实现，其转移操作是基于底层有限元的插值矩阵，以无矩阵的方式实现。这比MGTransferPrebuilt需要的内存少得多，也比该变体快得多。
 * 该类目前仅适用于基于FE_Q和FE_DGQ元素的张量积有限元，包括涉及这些元素之一的多个组件的系统。带有不同元素或其他元素的系统目前还没有实现。
 *
 *
 */
template <int dim, typename Number>
class MGTransferMatrixFree
  : public MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * 没有约束矩阵的构造函数。只在不连续的有限元或没有局部细化的情况下使用这个构造函数。
   *
   */
  MGTransferMatrixFree();

  /**
   * 带约束的构造器。相当于默认的构造函数，后面加上initialize_constraints()。
   *
   */
  MGTransferMatrixFree(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 解构器。
   *
   */
  virtual ~MGTransferMatrixFree() override = default;

  /**
   * 初始化将在build()中使用的约束。
   *
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 将对象重置为默认构造函数之后的状态。
   *
   */
  void
  clear();

  /**
   * 实际构建每一层的延长线的信息。
   * 外部分区器的可选第二个参数允许用户建议对各层进行矢量分区。如果发现分区器包含所有通过转移访问的幽灵未知数，则选择给定的分区器。这就保证了在延长和限制过程中矢量与用户给出的外部分区器的兼容性，这反过来又节省了一些复制操作。然而，在有未知数丢失的情况下
   *
   * 在h-coarsening过程中通常会出现这种情况，因为处理器需要退出，因此某个处理器上的子单元的未知数需要作为另一个处理器上的父单元的幽灵。
   *
   * - 提供的外部分区器被忽略，而使用内部的变体。
   *
   */
  void
  build(const DoFHandler<dim, dim> &dof_handler,
        const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          &external_partitioners =
            std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>());

  /**
   * 使用底层有限元的嵌入矩阵将一个向量从<tt>to_level-1</tt>延长到<tt>to_level</tt>。<tt>dst</tt>的先前内容被覆盖。
   * @param  to_level 要延长到的层次的索引，即 @p dst.   @param
   * src是一个向量，其元素数与涉及的较粗层次上的自由度相同。
   * @param
   * dst有多少个元素，就有多少个更细层次上的自由度。
   *
   */
  virtual void
  prolongate(
    const unsigned int                                to_level,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  virtual void
  prolongate_and_add(
    const unsigned int                                to_level,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * 使用prolongate()方法的转置操作，将一个向量从<tt>from_level</tt>级限制到<tt>from_level-1</tt>级。如果<tt>from_level</tt>层的单元所覆盖的区域小于<tt>from_level-1</tt>层的区域（局部细化），那么<tt>dst</tt>中的一些自由度是有效的，将不会被改变。对于其他自由度，限制的结果被添加。
   * @param  from_level 要限制的层次的索引，即 @p src.   @param
   * src是一个向量，其元素数量与所涉及的更细层次上的自由度相同。
   * @param
   * dst有多少个元素，就有多少个在较粗层次上的自由度。
   *
   */
  virtual void
  restrict_and_add(
    const unsigned int                                from_level,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * 将精细网格场 @p src 插值到 @p dof_handler
   * 中的每个多网格层次，并将结果存储在 @p dst.
   * 中。这个函数与限制不同，限制是将加权残差转移到较粗的层次（延长矩阵的转置）。
   * 参数 @p dst
   * 必须根据三角结构的层数以正确的大小进行初始化。
   * 如果 @p dst
   * 的内向量是空的或者有不正确的局部拥有的大小，它将被调整为每个层次上的局部相关自由度。
   * 这个函数的使用在 step-66 中得到了证明。
   *
   */
  template <typename Number2, int spacedim>
  void
  interpolate_to_mg(
    const DoFHandler<dim, spacedim> &                          dof_handler,
    MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
    const LinearAlgebra::distributed::Vector<Number2> &        src) const;

  /**
   * 有限元不提供延长矩阵。
   *
   */
  DeclException0(ExcNoProlongation);

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 一个变量，存储传递给build()的DoFHandler中包含的有限元的程度。计算内核的选择是基于这个数字的。
   *
   */
  unsigned int fe_degree;

  /**
   * 一个变量，存储该元素是否是连续的，并且在一维线的中心有一个联合自由度。
   *
   */
  bool element_is_continuous;

  /**
   * 一个变量，用于存储传递给build()的DoFHandler中包含的有限元的分量数量。
   *
   */
  unsigned int n_components;

  /**
   * 一个存储所有子单元的自由度数量的变量。对于DG元素来说，它是<tt>2<sup>dim</sup>*fe.n_dofs_per_cell()</tt>，对于连续元素来说，它的数量要少一些。
   *
   */
  unsigned int n_child_cell_dofs;

  /**
   * 这个变量保存了给定层次上的单元格的索引，从DoFHandler中提取出来以便快速访问。一个给定层次上的所有DoF指数被存储为一个普通数组（因为这个类假设每个单元的DoF是恒定的）。要索引到这个数组，使用单元格编号乘以dofs_per_cell。
   * 这个数组首先被安排成所有本地拥有的级别的单元首先出现（在变量n_owned_level_cells中找到），然后是转移到下一级别所需的其他单元。
   *
   */
  std::vector<std::vector<unsigned int>> level_dof_indices;

  /**
   * 一个变量存储了每个级别的父级到子级单元格编号的连接。
   *
   */
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
    parent_child_connect;

  /**
   * 一个存储每一级在给定进程上拥有的单元格数量的变量（为工作者循环设定界限）。
   *
   */
  std::vector<unsigned int> n_owned_level_cells;

  /**
   * 这个变量持有从母元素到所有子元素的一维嵌入（延长）矩阵。
   *
   */
  AlignedVector<VectorizedArray<Number>> prolongation_matrix_1d;

  /**
   * 这个变量保存了张量评估的临时值。
   *
   */
  mutable AlignedVector<VectorizedArray<Number>> evaluation_data;

  /**
   * 对于连续元素，限制不是加法的，我们需要在延长结束时（和限制开始时）根据自由度的价位对结果进行加权，也就是说，根据它们出现的元素数量进行加权。我们以矢量的形式存储数据，以允许廉价的访问。此外，我们利用了我们只需要存储<tt>3<sup>dim</sup></tt>指数的事实。
   * 数据以每一层（外向量）和每一层的单元格（内向量）为单位组织。
   *
   */
  std::vector<AlignedVector<VectorizedArray<Number>>> weights_on_refined;

  /**
   * 一个存储所有层次（外指数）、层次内的单元（第二指数）和单元上的指数（内指数）的Dirichlet边界条件的局部指数的变量。
   *
   */
  std::vector<std::vector<std::vector<unsigned short>>> dirichlet_indices;

  /**
   * 一个向量，持有转移的分区器的共享指针。这些分区器可能与通过build()从外部传入的东西共享，或者与从MGLevelGlobalTransfer继承的级别向量共享。
   *
   */
  MGLevelObject<std::shared_ptr<const Utilities::MPI::Partitioner>>
    vector_partitioners;

  /**
   * 执行延长操作。
   *
   */
  template <int degree>
  void
  do_prolongate_add(
    const unsigned int                                to_level,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * 执行限制性操作。
   *
   */
  template <int degree>
  void
  do_restrict_add(const unsigned int                                from_level,
                  LinearAlgebra::distributed::Vector<Number> &      dst,
                  const LinearAlgebra::distributed::Vector<Number> &src) const;
};


/**
 * MGTransferBase接口的实现，其转移操作是基于底层有限元的插值矩阵，以无矩阵的方式实现。这比MGTransferPrebuilt需要的内存少得多，也比该变体快得多。
 * 该类与 LinearAlgebra::distributed::BlockVector
 * 一起工作，对每个区块执行与MGTransferMatrixFree完全相同的转移操作。支持所有块使用相同的DoFHandler和每个块使用自己的DoFHandler两种情况。
 *
 *
 */
template <int dim, typename Number>
class MGTransferBlockMatrixFree
  : public MGTransferBase<LinearAlgebra::distributed::BlockVector<Number>>
{
public:
  /**
   * 没有约束矩阵的构造函数。只在不连续的有限元或没有局部细化的情况下使用这个构造函数。
   *
   */
  MGTransferBlockMatrixFree() = default;

  /**
   * 带约束的构造器。相当于默认的构造函数，后面加上initialize_constraints()。
   *
   */
  MGTransferBlockMatrixFree(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 与上面的情况相同，每个块都有自己的DoFHandler。
   *
   */
  MGTransferBlockMatrixFree(
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * 解构器。
   *
   */
  virtual ~MGTransferBlockMatrixFree() override = default;

  /**
   * 初始化将在build()中使用的约束。
   *
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 对于每个块都有自己的DoFHandler的情况，与上述相同。
   *
   */
  void
  initialize_constraints(
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * 将对象重置为默认构造函数后的状态。
   *
   */
  void
  clear();

  /**
   * 实际建立每一层的延长线的信息。
   *
   */
  void
  build(const DoFHandler<dim, dim> &dof_handler);

  /**
   * 在每个区块都有自己的DoFHandler的情况下，与上面一样。
   *
   */
  void
  build(const std::vector<const DoFHandler<dim, dim> *> &dof_handler);

  /**
   * 使用底层有限元的嵌入矩阵，将一个向量从<tt>to_level-1</tt>延长到<tt>to_level</tt>。<tt>dst</tt>的先前内容被覆盖。
   * @param  to_level 要延长到的层次的索引，也就是 @p dst.
   * @param
   * src是一个向量，其元素数量与所涉及的较粗层次上的自由度相同。
   * @param
   * dst有多少个元素，就有多少个更细层次上的自由度。
   *
   */
  virtual void
  prolongate(
    const unsigned int                                     to_level,
    LinearAlgebra::distributed::BlockVector<Number> &      dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  virtual void
  prolongate_and_add(
    const unsigned int                                     to_level,
    LinearAlgebra::distributed::BlockVector<Number> &      dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  /**
   * 使用prolongate()方法的转置操作，将一个向量从<tt>from_level</tt>级限制到<tt>from_level-1</tt>级。如果<tt>from_level</tt>层的单元所覆盖的区域小于<tt>from_level-1</tt>层的区域（局部细化），那么<tt>dst</tt>中的一些自由度是有效的，将不会被改变。对于其他自由度，限制的结果被添加。
   * @param  from_level 要限制的层次的索引，即 @p src.   @param
   * src是一个向量，其元素数量与所涉及的更细层次上的自由度相同。
   * @param
   * dst有多少个元素，就有多少个在较粗层次上的自由度。
   *
   */
  virtual void
  restrict_and_add(
    const unsigned int                                     from_level,
    LinearAlgebra::distributed::BlockVector<Number> &      dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  /**
   * 从全局网格上的块向量转移到每个层次上的块向量，分别定义为活动自由度。
   * 特别是，对于一个全局细化的网格，只有 @p dst
   * 中最细的层次被填充为 @p src.
   * 的普通拷贝，其他的层次对象都没有被触动。
   * 如果需要的话，这个函数会根据Multigrid类的要求，相应地初始化
   * @p dst 。
   *
   */
  template <typename Number2, int spacedim>
  void
  copy_to_mg(
    const DoFHandler<dim, spacedim> &                               dof_handler,
    MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
    const LinearAlgebra::distributed::BlockVector<Number2> &        src) const;

  /**
   * 对于每个块都有自己的DoFHandler的情况，与上述相同。
   *
   */
  template <typename Number2, int spacedim>
  void
  copy_to_mg(
    const std::vector<const DoFHandler<dim, spacedim> *> &          dof_handler,
    MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
    const LinearAlgebra::distributed::BlockVector<Number2> &        src) const;

  /**
   * 从多级块-向量转移到法向量。
   *
   */
  template <typename Number2, int spacedim>
  void
  copy_from_mg(
    const DoFHandler<dim, spacedim> &                 dof_handler,
    LinearAlgebra::distributed::BlockVector<Number2> &dst,
    const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
    const;

  /**
   * 与上述每个块都有自己的DoFHandler的情况相同。
   *
   */
  template <typename Number2, int spacedim>
  void
  copy_from_mg(
    const std::vector<const DoFHandler<dim, spacedim> *> &dof_handler,
    LinearAlgebra::distributed::BlockVector<Number2> &    dst,
    const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
    const;

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 这个类既可以用一个DoFHandler，也可以为每个区块用一个单独的DoFHandler。
   *
   */
  static const bool supports_dof_handler_vector = true;

private:
  /**
   * 非块的无矩阵版本的转移操作。
   *
   */
  std::vector<MGTransferMatrixFree<dim, Number>> matrix_free_transfer_vector;

  /**
   * 一个标志，表示所有组件使用同一个DoFHandler，还是每个块都有自己的DoFHandler。
   *
   */
  const bool same_for_all;
};


 /*@}*/ 


//------------------------ templated functions -------------------------
#ifndef DOXYGEN


template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferMatrixFree<dim, Number>::interpolate_to_mg(
  const DoFHandler<dim, spacedim> &                          dof_handler,
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
  const LinearAlgebra::distributed::Vector<Number2> &        src) const
{
  const unsigned int min_level = dst.min_level();
  const unsigned int max_level = dst.max_level();

  Assert(max_level == dof_handler.get_triangulation().n_global_levels() - 1,
         ExcDimensionMismatch(
           max_level, dof_handler.get_triangulation().n_global_levels() - 1));

  const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();

  for (unsigned int level = min_level; level <= max_level; ++level)
    if (dst[level].size() != dof_handler.n_dofs(level) ||
        dst[level].locally_owned_size() !=
          dof_handler.locally_owned_mg_dofs(level).n_elements())
      dst[level].reinit(this->vector_partitioners[level]);

  // copy fine level vector to active cells in MG hierarchy
  this->copy_to_mg(dof_handler, dst, src, true);

  // FIXME: maybe need to store hanging nodes constraints per level?
  // MGConstrainedDoFs does NOT keep this info right now, only periodicity
  // constraints...

  // do the transfer from level to level-1:
  dst[max_level].update_ghost_values();
  for (unsigned int level = max_level; level > min_level; --level)
    {
      // auxiliary vector which always has ghost elements
      const LinearAlgebra::distributed::Vector<Number> *input = nullptr;
      LinearAlgebra::distributed::Vector<Number>        ghosted_fine;
      if (dst[level].get_partitioner().get() ==
          this->vector_partitioners[level].get())
        input = &dst[level];
      else
        {
          ghosted_fine.reinit(this->vector_partitioners[level]);
          ghosted_fine.copy_locally_owned_data_from(dst[level]);
          ghosted_fine.update_ghost_values();
          input = &ghosted_fine;
        }

      std::vector<Number> dof_values_coarse(fe.n_dofs_per_cell());
      Vector<Number>      dof_values_fine(fe.n_dofs_per_cell());
      Vector<Number>      tmp(fe.n_dofs_per_cell());
      std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());
      for (const auto &cell : dof_handler.cell_iterators_on_level(level - 1))
        if (cell->is_locally_owned_on_level())
          {
            // if we get to a cell without children (== active), we can
            // skip it as there values should be already set by the
            // equivalent of copy_to_mg()
            if (cell->is_active())
              continue;

            std::fill(dof_values_coarse.begin(), dof_values_coarse.end(), 0.);
            for (unsigned int child = 0; child < cell->n_children(); ++child)
              {
                cell->child(child)->get_mg_dof_indices(dof_indices);
                for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
                  dof_values_fine(i) = (*input)(dof_indices[i]);
                fe.get_restriction_matrix(child, cell->refinement_case())
                  .vmult(tmp, dof_values_fine);
                for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
                  if (fe.restriction_is_additive(i))
                    dof_values_coarse[i] += tmp[i];
                  else if (tmp(i) != 0.)
                    dof_values_coarse[i] = tmp[i];
              }
            cell->get_mg_dof_indices(dof_indices);
            for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
              if (dof_handler.locally_owned_mg_dofs(level - 1).is_element(
                    dof_indices[i]))
                dst[level - 1](dof_indices[i]) = dof_values_coarse[i];
          }

      dst[level - 1].update_ghost_values();
    }
}



template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_to_mg(
  const DoFHandler<dim, spacedim> &                               dof_handler,
  MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
  const LinearAlgebra::distributed::BlockVector<Number2> &        src) const
{
  AssertDimension(matrix_free_transfer_vector.size(), 1);
  Assert(same_for_all,
         ExcMessage(
           "This object was initialized with support for usage with one "
           "DoFHandler for each block, but this method assumes that "
           "the same DoFHandler is used for all the blocks!"));
  const std::vector<const DoFHandler<dim, spacedim> *> mg_dofs(src.n_blocks(),
                                                               &dof_handler);

  copy_to_mg(mg_dofs, dst, src);
}



template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_to_mg(
  const std::vector<const DoFHandler<dim, spacedim> *> &          dof_handler,
  MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
  const LinearAlgebra::distributed::BlockVector<Number2> &        src) const
{
  const unsigned int n_blocks = src.n_blocks();
  AssertDimension(dof_handler.size(), n_blocks);

  if (n_blocks == 0)
    return;

  const unsigned int min_level = dst.min_level();
  const unsigned int max_level = dst.max_level();

  // this function is normally called within the Multigrid class with
  // dst == defect level block vector. At first run this vector is not
  // initialized. Do this below:
  {
    const parallel::TriangulationBase<dim, spacedim> *tria =
      (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &(dof_handler[0]->get_triangulation())));
    for (unsigned int i = 1; i < n_blocks; ++i)
      AssertThrow(
        (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
           &(dof_handler[0]->get_triangulation())) == tria),
        ExcMessage("The DoFHandler use different Triangulations!"));

    MGLevelObject<bool> do_reinit;
    do_reinit.resize(min_level, max_level);
    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        do_reinit[level] = false;
        if (dst[level].n_blocks() != n_blocks)
          {
            do_reinit[level] = true;
            continue; // level
          }
        for (unsigned int b = 0; b < n_blocks; ++b)
          {
            LinearAlgebra::distributed::Vector<Number> &v = dst[level].block(b);
            if (v.size() !=
                  dof_handler[b]->locally_owned_mg_dofs(level).size() ||
                v.locally_owned_size() !=
                  dof_handler[b]->locally_owned_mg_dofs(level).n_elements())
              {
                do_reinit[level] = true;
                break; // b
              }
          }
      }

    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        if (do_reinit[level])
          {
            dst[level].reinit(n_blocks);
            for (unsigned int b = 0; b < n_blocks; ++b)
              {
                LinearAlgebra::distributed::Vector<Number> &v =
                  dst[level].block(b);
                v.reinit(dof_handler[b]->locally_owned_mg_dofs(level),
                         dof_handler[b]->get_communicator());
              }
            dst[level].collect_sizes();
          }
        else
          dst[level] = 0;
      }
  }

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> dst_non_block(
    min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (unsigned int l = min_level; l <= max_level; ++l)
        dst_non_block[l].reinit(dst[l].block(b));
      const unsigned int data_block = same_for_all ? 0 : b;
      matrix_free_transfer_vector[data_block].copy_to_mg(*dof_handler[b],
                                                         dst_non_block,
                                                         src.block(b));

      for (unsigned int l = min_level; l <= max_level; ++l)
        dst[l].block(b) = dst_non_block[l];
    }
}

template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_from_mg(
  const DoFHandler<dim, spacedim> &                 dof_handler,
  LinearAlgebra::distributed::BlockVector<Number2> &dst,
  const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
  const
{
  AssertDimension(matrix_free_transfer_vector.size(), 1);
  const std::vector<const DoFHandler<dim, spacedim> *> mg_dofs(dst.n_blocks(),
                                                               &dof_handler);

  copy_from_mg(mg_dofs, dst, src);
}

template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_from_mg(
  const std::vector<const DoFHandler<dim, spacedim> *> &dof_handler,
  LinearAlgebra::distributed::BlockVector<Number2> &    dst,
  const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
  const
{
  const unsigned int n_blocks = dst.n_blocks();
  AssertDimension(dof_handler.size(), n_blocks);

  if (n_blocks == 0)
    return;

  const unsigned int min_level = src.min_level();
  const unsigned int max_level = src.max_level();

  for (unsigned int l = min_level; l <= max_level; ++l)
    AssertDimension(src[l].n_blocks(), dst.n_blocks());

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> src_non_block(
    min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          src_non_block[l].reinit(src[l].block(b));
          src_non_block[l] = src[l].block(b);
        }
      const unsigned int data_block = same_for_all ? 0 : b;
      matrix_free_transfer_vector[data_block].copy_from_mg(*dof_handler[b],
                                                           dst.block(b),
                                                           src_non_block);
    }
}



#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif


