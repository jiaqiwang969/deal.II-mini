//include/deal.II-translator/multigrid/mg_transfer_internal_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


#ifndef dealii_mg_transfer_internal_h
#define dealii_mg_transfer_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MGTransfer
  {
    /**
     * 内部函数，用于填充从全局索引到层次索引的复制索引
     * 如果 @p skip_interface_dofs
     * 为假，映射也将包含层次间接口的DoF。这在传输解向量而不是残差时是可取的。
     *
     */
    template <int dim, int spacedim>
    void
    fill_copy_indices(
      const DoFHandler<dim, spacedim> &dof_handler,
      const MGConstrainedDoFs *        mg_constrained_dofs,
      std::vector<std::vector<
        std::pair<types::global_dof_index, types::global_dof_index>>>
        &copy_indices,
      std::vector<std::vector<
        std::pair<types::global_dof_index, types::global_dof_index>>>
        &copy_indices_global_mine,
      std::vector<std::vector<
        std::pair<types::global_dof_index, types::global_dof_index>>>
        &        copy_indices_level_mine,
      const bool skip_interface_dofs = true);



    /**
     * 给出从父体看到的子体单元的集合，这个函数计算出给定子体的第一个索引
     *
     */
    template <int dim>
    unsigned int
    compute_shift_within_children(const unsigned int child,
                                  const unsigned int fe_shift_1d,
                                  const unsigned int fe_degree);

    /**
     * 一个存储与DoFHandler中包含的有限元相关的数据的结构。仅用于使用<tt>setup_transfer</tt>的初始化。
     *
     */
    template <typename Number>
    struct ElementInfo
    {
      /**
       * 一个存储有限元程度的变量。计算内核的选择是基于这个数字的。
       *
       */
      unsigned int fe_degree;

      /**
       * 一个变量，存储元素是否是连续的，在一维线的中心有一个联合自由度。
       *
       */
      bool element_is_continuous;

      /**
       * 一个存储有限元中分量数量的变量。
       *
       */
      unsigned int n_components;

      /**
       * 一个存储所有子单元上自由度数量的变量。
       * 对于DG元素来说，它是<tt>2<sup>dim</sup>*fe.n_dofs_per_cell()</tt>，对于连续元素来说，它略少一些。
       *
       */
      unsigned int n_child_cell_dofs;

      /**
       * 一个数组，用于保存有限元中自由度的编号和张量乘积应用所需的词法编号之间的编号。
       *
       */
      std::vector<unsigned int> lexicographic_numbering;

      /**
       * 这个变量持有从母元素到所有子元素的一维嵌入（延长）矩阵。
       *
       */
      std::vector<Number> prolongation_matrix_1d;
    };

    /**
     * 设置MGTransferMatrixFree的大部分内部数据结构
     *
     */
    template <int dim, typename Number>
    void
    setup_transfer(
      const DoFHandler<dim> &  dof_handler,
      const MGConstrainedDoFs *mg_constrained_dofs,
      const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        &                                     external_partitioners,
      ElementInfo<Number> &                   elem_info,
      std::vector<std::vector<unsigned int>> &level_dof_indices,
      std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
        &                        parent_child_connect,
      std::vector<unsigned int> &n_owned_level_cells,
      std::vector<std::vector<std::vector<unsigned short>>> &dirichlet_indices,
      std::vector<std::vector<Number>> &                     weights_on_refined,
      std::vector<Table<2, unsigned int>> &copy_indices_global_mine,
      MGLevelObject<std::shared_ptr<const Utilities::MPI::Partitioner>>
        &vector_partitioners);

  } // namespace MGTransfer
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


