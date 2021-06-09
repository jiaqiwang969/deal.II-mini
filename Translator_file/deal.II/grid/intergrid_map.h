//include/deal.II-translator/grid/intergrid_map_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_intergrid_map_h
#define dealii_intergrid_map_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

DEAL_II_NAMESPACE_OPEN


/**
 * 这个类提供了两个网格之间的映射，这两个网格来自于同一个粗略的网格。对于源地图的每个单元格迭代器，它通过其<tt>operator
 * []</tt>，提供了目标地图上相应的单元格迭代器。
 * 通常情况下，两个网格的细化程度是不同的。那么，源网格上的迭代器返回的值将是：。  <ul>   <li>  目标网格上的同一个单元，如果它存在的话；  <li>  目标网格上最精炼的单元，通过精炼可以得到源单元的垂线。该单元始终处于活动状态，其细化程度小于源单元的细化程度。  </ul>  这个地图的键是源网格上的所有单元，无论是否活动。
 * 例如，考虑这两个一维网格。
 *
 * @verbatim
 * Grid 1:
 * x--x--x-----x-----------x
 *  1  2    3        4
 *
 * Grid 2:
 * x-----x-----x-----x-----x
 *    1     2     3     4
 * @endverbatim
 * （单元格编号只是作为一个例子给出，不会与真实的单元格迭代器的索引相对应）。那么从网格1到网格2的映射将如下。
 *
 * @verbatim
 *  Cell on grid 1         Cell on grid 2
 *        1
 *
 *
 *
 *
 *
 * ------------------>  1
 *        2
 *
 *
 *
 *
 *
 * ------------------>  1
 *        3
 *
 *
 *
 *
 *
 * ------------------>  2
 *        4
 *
 *
 *
 *
 *
 * ------------------>  mother cell of cells 3 and 4
 *                                (a non-active cell, not shown here)
 * @endverbatim
 * 除了这里显示的映射，网格1上的非活动单元也是有效的键。例如，第一个网格上的单元格1和2的母单元格的映射将指向第二个网格上的单元格1。
 * @tparam  MeshType 这个类可以和任何满足 @ref ConceptMeshType  "MeshType概念 "
 * 的类一起使用。对其他提供迭代器函数的类的扩展和一些小的附加要求很简单。
 * 请注意，这个类原则上可以基于C++  <tt>std::map<Key,Value></tt>
 * 数据类型。相反，它使用了另一种数据格式，在访问的计算时间和内存消耗方面都更有效。
 *
 *  <h3>Usage</h3> 在实践中，该类的使用情况如下。
 *
 * @code
 * // have two grids, which are derived from the same coarse grid
 * Triangulation<dim> tria1, tria2;
 * DoFHandler<dim> dof_handler_1 (tria1), dof_handler_2 (tria2);
 * ...
 * // do something with these objects, e.g. refine the triangulations
 * // differently, distribute degrees of freedom, etc
 * ...
 * // create the mapping
 * InterGridMap<DoFHandler<dim> > grid_1_to_2_map;
 * grid_1_to_2_map.make_mapping (dof_handler_1,
 *                               dof_handler_2);
 * ...
 * typename DoFHandler<dim>::cell_iterator cell = dof_handler_1.begin(),
 *                                         endc = dof_handler_1.end();
 * for (; cell!=endc; ++cell)
 *   // now do something with the cell of dof_handler_2 corresponding to
 *   // cell (which is one of dof_handler_1's cells)
 *   f (grid_1_to_2_map[cell]);
 * @endcode
 *
 * 注意这个类的模板参数必须以<tt>InterGridMap<DoFHandler<2>></tt>的形式给出，这里是DoFHandler（也同样可以是Triangulation或PersistentTriangulation）。
 *
 *
 * @ingroup grid
 *
 */
template <class MeshType>
class InterGridMap : public Subscriptor
{
public:
  /**
   * 对所考虑的网格类的迭代器类型的类型化定义。
   *
   */
  using cell_iterator = typename MeshType::cell_iterator;

  /**
   * 构造函数设置SmartPointer成员中的类名参数。
   *
   */
  InterGridMap();

  /**
   * 创建两个网格之间的映射。
   *
   */
  void
  make_mapping(const MeshType &source_grid, const MeshType &destination_grid);

  /**
   * 访问操作符：给出源网格上的一个单元，并接收另一个网格上相应的单元，如果不存在，则是源单元进一步细化后将创建的最细化单元。
   *
   */
  cell_iterator operator[](const cell_iterator &source_cell) const;

  /**
   * 删除该类的所有数据。
   *
   */
  void
  clear();

  /**
   * 返回一个对源网格的引用。
   *
   */
  const MeshType &
  get_source_grid() const;

  /**
   * 返回一个对目标网格的引用。
   *
   */
  const MeshType &
  get_destination_grid() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidKey,
                 cell_iterator,
                 << "The iterator " << arg1 << " is not valid as key for "
                 << "this map.");
  /**
   * 异常情况
   *
   */
  DeclException0(ExcIncompatibleGrids);

private:
  /**
   * 实际的数据。为每一层的每个单元格保留一个迭代器。
   *
   */
  std::vector<std::vector<cell_iterator>> mapping;

  /**
   * 存储一个指向源网格的指针。
   *
   */
  SmartPointer<const MeshType, InterGridMap<MeshType>> source_grid;

  /**
   * 同样，对于目标网格也是如此。
   *
   */
  SmartPointer<const MeshType, InterGridMap<MeshType>> destination_grid;

  /**
   * 为给定的一对单元格设置映射。这些应该在细化程度和所有其他属性上匹配。
   *
   */
  void
  set_mapping(const cell_iterator &src_cell, const cell_iterator &dst_cell);

  /**
   * 将键 @p src_cell 的值设置为 @p dst_cell.  对 @p src_cell.
   * 的所有子单元和它们的子单元也是如此 这个函数用于在
   * @p src_grid 上比在 @p dst_grid;
   * 上更精细的单元，那么所有单元的层次及其子单元的值都指向
   * @p dst_grid.  上的一个单元
   *
   */
  void
  set_entries_to_cell(const cell_iterator &src_cell,
                      const cell_iterator &dst_cell);
};


DEAL_II_NAMESPACE_CLOSE

#endif


