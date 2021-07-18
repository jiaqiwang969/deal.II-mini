//include/deal.II-translator/numerics/cell_data_transfer_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_cell_data_transfer_h
#define dealii_cell_data_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/grid/tria.h>

#include <deal.II/numerics/adaptation_strategies.h>

#include <algorithm>
#include <functional>
#include <set>


DEAL_II_NAMESPACE_OPEN

/**
 * 在细化和/或粗化三角测量时，传输与每个活动单元相关的数据（如误差指标）。
 * 因此，这个类对于单元相关的信息，就像SolutionTransfer对于定义在三角剖分上的自由度值一样。
 * 必须提供一个非分布式的容器（如Vector或 `std::vector`)
 * ），它以活动单元被遍历的相同顺序持有单元的数据。换句话说，每个条目都对应于具有相同索引
 * CellAccessor::active_cell_index(), 的单元，容器的大小必须是
 * Triangulation::n_active_cells().  。 <h3>Transferring cell-wise data</h3>
 * 下面的代码片断演示了如何在细化/粗化注册三角的过程中传输单元的相关数据。
 *
 *
 * @code
 * // prepare the triangulation,
 * triangulation.prepare_coarsening_and_refinement();
 *
 * // prepare the CellDataTransfer object for coarsening and refinement
 * // and give the cell data vector that we intend to unpack later,
 * Vector<double> data_to_transfer(triangulation.n_active_cells());
 * //[fill data_to_transfer with cell-wise values...]
 *
 * CellDataTransfer<dim, spacedim, Vector<double>>
 * cell_data_trans(triangulation);
 * cell_data_trans.prepare_for_coarsening_and_refinement();
 *
 * // actually execute the refinement,
 * triangulation.execute_coarsening_and_refinement();
 *
 * // unpack transferred data,
 * Vector<double> transferred_data(triangulation.n_active_cells());
 * cell_data_trans.unpack(data_to_transfer, transferred_data);
 * @endcode
 *
 * 当使用 parallel::shared::Triangulation,
 * 时，我们需要确保在细化发生之前，我们的本地向量中有全局数据。我们可以通过以下方式实现这一点。
 *
 *
 * @code
 * Vector<double> data_to_transfer(triangulation.n_active_cells());
 * //[fill data_to_transfer with cell-wise values...]
 *
 * PETScWrappers::MPI::Vector
 * distributed_data_to_transfer(mpi_communicator,
 *                            triangulation.n_active_cells(),
 *                            triangulation.n_locally_owned_active_cells());
 * for (const auto &cell : triangulation.active_cell_iterators())
 * if (cell->is_locally_owned())
 *   {
 *     const unsigned int index = cell->active_cell_index();
 *     distributed_data_to_transfer(index) = data_to_transfer(index);
 *   }
 * distributed_data_to_transfer.compress(VectorOperation::insert);
 *
 * data_to_transfer = distributed_data_to_transfer;
 * @endcode
 *
 * 对于并行分布的情况，有一个指定的类
 * parallel::distributed::CellDataTransfer 。在使用
 * parallel::distributed::Triangulation. 时，请参考这个特定的类。
 *
 *
 * @note
 * 参见SolutionTransfer的文档，以了解传输的匹配代码片段。
 *
 *
 * @ingroup numerics
 *
 *
 */
template <int dim, int spacedim = dim, typename VectorType = Vector<double>>
class CellDataTransfer
{
private:
  /**
   * 一个别名，定义了所提供的容器模板的数据类型。
   *
   */
  using value_type = typename VectorType::value_type;

public:
  /**
   * 构造函数。      @param[in]  triangulation
   * 所有操作都将发生在这个三角形上。当这个构造函数被调用时，有关的细化还没有发生。
   * @param[in]  refinement_strategy
   * %函数，决定如何将数据从其父单元存储到细化单元上。
   * @param[in]  coarsening_strategy
   * %函数决定在子单元上存储哪些数据，其子单元将被粗化为。
   *
   */
  CellDataTransfer(
    const Triangulation<dim, spacedim> &triangulation,
    const std::function<std::vector<value_type>(
      const typename Triangulation<dim, spacedim>::cell_iterator &parent,
      const value_type parent_value)>   refinement_strategy =
      &AdaptationStrategies::Refinement::preserve<dim, spacedim, value_type>,
    const std::function<value_type(
      const typename Triangulation<dim, spacedim>::cell_iterator &parent,
      const std::vector<value_type> &children_values)> coarsening_strategy =
      &AdaptationStrategies::Coarsening::
        check_equality<dim, spacedim, value_type>);

  /**
   * 为粗化和细化当前对象做准备。
   * 存储相关三角形上所有活动单元的active_cell_indices，并将它们归属于持久化、细化或粗化的单元。
   *
   */
  void
  prepare_for_coarsening_and_refinement();

  /**
   * 将上一个网格的信息转移到更新的网格中。    由 @p in
   * 提供的前一个网格的数据将被转移到更新的网格，并存储在
   * @p out.  @p out
   * 必须提供足够的空间来容纳转移的数据，即必须是`triangulation.n_active_cells()`的大小。
   *
   */
  void
  unpack(const VectorType &in, VectorType &out);

private:
  /**
   * 指向要处理的三角结构的指针。
   *
   */
  SmartPointer<const Triangulation<dim, spacedim>,
               CellDataTransfer<dim, spacedim, VectorType>>
    triangulation;

  /**
   * 决定数据如何从其父单元存储到细化单元的函数。
   *
   */
  const std::function<std::vector<value_type>(
    const typename Triangulation<dim, spacedim>::cell_iterator &parent,
    const value_type                                            parent_value)>
    refinement_strategy;

  /**
   * 决定如何处理来自子单元的数据以存储在父单元上的函数。
   *
   */
  const std::function<value_type(
    const typename Triangulation<dim, spacedim>::cell_iterator &parent,
    const std::vector<value_type> &children_indices)>
    coarsening_strategy;

  /**
   * 容器，用于临时存储持续存在的单元格的迭代器和活动单元格索引。
   *
   */
  std::map<const typename Triangulation<dim, spacedim>::cell_iterator,
           const unsigned int>
    persisting_cells_active_index;

  /**
   * 用于临时存储将被精炼的单元格的迭代器和活动单元格索引的容器。
   *
   */
  std::map<const typename Triangulation<dim, spacedim>::cell_iterator,
           const unsigned int>
    refined_cells_active_index;

  /**
   * 容器用于临时存储粗化后将保留的父单元格的迭代器以及相应子单元格的活动单元格索引。
   *
   */
  std::map<const typename Triangulation<dim, spacedim>::cell_iterator,
           const std::set<unsigned int>>
    coarsened_cells_active_index;

  /**
   * 初始三角形上尚未被细化的活动单元的数量。
   * 它将在prepare_for_coarsening_and_refinement()中设置，并在细化发生后用于验证用户输入（仅在调试模式）。
   *
   */
  unsigned int n_active_cells_pre;
};


DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_cell_data_transfer_h */ 


