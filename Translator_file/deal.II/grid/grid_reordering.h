//include/deal.II-translator/grid/grid_reordering_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_grid_reordering_h
#define dealii_grid_reordering_h


#include <deal.II/base/config.h>

#include <deal.II/grid/tria.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个实现各种网格重排算法的类。更多信息见 @ref reordering "重新排序模块"
 * 。
 * @deprecated 使用 GridTools::invert_all_negative_measure_cells() 或
 * GridTools::consistently_order_cells()
 * 代替这个类所提供的函数。旧式编号的使用已被弃用。
 *
 *
 */
template <int dim, int spacedim = dim>
class DEAL_II_DEPRECATED GridReordering
{
public:
  /**
   * 这是主函数，在dim=2和3的情况下做这个类的一般文档中所宣布的事情，在dim=1的情况下不做任何事情。
   * 如果在dim=3的情况下不可能进行一致的重新排序，则恢复原始连接数据。
   * @param  original_cells 一个包含描述网格的数据的对象。
   * @param  use_new_style_ordering
   * 如果是true，则使用单元格内顶点的标准排序。如果是false（默认），则使用deal.II在5.2版本之前使用的单元格内顶点的
   * "旧式 "排序，在这个类的文档中也有解释。
   * @deprecated  使用 GridTools::consistently_order_cells() 代替。
   *
   */
  DEAL_II_DEPRECATED
  static void
  reorder_cells(std::vector<CellData<dim>> &original_cells,
                const bool                  use_new_style_ordering = false);

  /**
   * 由网格生成器生成的网格可能具有与deal.II所要求的方向相反的单元格方向。
   * 在2D和3D中，这个函数检查所有单元是否有负的或正的度量/体积。在前一种情况下，所有的单元格都是倒置的。在1d中，它没有任何作用。
   * 当所有单元格中只有一个子集的体积为负时，单元格的反转也可能起作用。然而，由负向和正向单元混合组成的网格很可能被打破。因此，如果单元格的方向不一致，就会抛出一个异常。
   * 注意，这个函数应该在reorder_cells()之前调用。      @param
   * all_vertices 网格的顶点。    @param  original_cells
   * 一个包含描述网格的数据的对象。    @param
   * use_new_style_ordering
   * 如果为true，则使用单元格内顶点的标准排序。如果是false（默认），则使用deal.II在5.2版本之前使用的单元格内顶点的
   * "旧式 "排序，在这个类的文档中也有解释。
   * @deprecated  使用 GridTools::invert_all_negative_measure_cells()
   * 代替。
   *
   */
  DEAL_II_DEPRECATED
  static void
  invert_all_cells_of_negative_grid(
    const std::vector<Point<spacedim>> &all_vertices,
    std::vector<CellData<dim>> &        original_cells,
    const bool                          use_new_style_ordering = false);
};


// declaration of explicit specializations
template <>
void
GridReordering<2>::invert_all_cells_of_negative_grid(
  const std::vector<Point<2>> &all_vertices,
  std::vector<CellData<2>> &   cells,
  const bool                   use_new_style_ordering);

template <>
void
GridReordering<2, 3>::invert_all_cells_of_negative_grid(
  const std::vector<Point<3>> &all_vertices,
  std::vector<CellData<2>> &   cells,
  const bool                   use_new_style_ordering);

template <>
void
GridReordering<3>::invert_all_cells_of_negative_grid(
  const std::vector<Point<3>> &all_vertices,
  std::vector<CellData<3>> &   cells,
  const bool                   use_new_style_ordering);

DEAL_II_NAMESPACE_CLOSE

#endif


