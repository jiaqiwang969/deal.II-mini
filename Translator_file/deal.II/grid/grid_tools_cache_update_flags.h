//include/deal.II-translator/grid/grid_tools_cache_update_flags_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_grid_tria_info_cache_update_flags_h
#define dealii_grid_tria_info_cache_update_flags_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  /**
   * 给予Cache类的枚举类型，用于选择要更新的信息。
   * 你可以用位法或
   * <code>operator|(CacheUpdateFlags,CacheUpdateFlags)</code>
   * 串联选择一个以上的标志。
   *
   */
  enum CacheUpdateFlags
  {
    /**
     * 更新无。
     *
     */
    update_nothing = 0x000,

    /**
     * 更新Vertex_to_cell_map，如 GridTools::vertex_to_cell_map().
     * 所返回的。
     *
     */
    update_vertex_to_cell_map = 0x001,

    /**
     * 更新由 GridTools::vertex_to_cell_centers_directions()
     * 返回的Vertex_to_cell_centers_directions。
     *
     */
    update_vertex_to_cell_centers_directions =
      update_vertex_to_cell_map | 0x002,

    /**
     * 更新已用顶点的映射。
     *
     */
    update_used_vertices = 0x008,

    /**
     * 更新一个已使用顶点的RTree。
     *
     */
    update_used_vertices_rtree = 0x010,

    /**
     * 更新单元格包围盒的RT树。
     *
     */
    update_cell_bounding_boxes_rtree = 0x020,

    /**
     * 更新覆盖的rtree对象，初始化为界线盒和无符号int对。边界框用于描述网格的哪一部分包含本地拥有的单元，通过对的第二个元素的排序过程。
     *
     */
    update_covering_rtree = 0x040,

    /**
     * 更新本地拥有的单元格边界盒的RTree。
     *
     */
    update_locally_owned_cell_bounding_boxes_rtree = 0x080,

    /**
     * 更新顶点到邻居子域
     *
     */
    update_vertex_to_neighbor_subdomain = 0x100,

    /**
     * 更新所有对象。
     *
     */
    update_all = 0xFFF,
  };


  /**
   * 输出操作符，将集合标志作为一组或的文本值输出。
   * @ref CacheUpdateFlags
   *
   */
  template <class StreamType>
  inline StreamType &
  operator<<(StreamType &s, const CacheUpdateFlags u)
  {
    s << " CacheUpdateFlags";
    if (u & update_vertex_to_cell_map)
      s << "|vertex_to_cell_map";
    if (u & update_vertex_to_cell_centers_directions)
      s << "|vertex_to_cells_centers_directions";
    if (u & update_covering_rtree)
      s << "|covering_rtree";
    return s;
  }


  /**
   * 全局运算符，它返回一个对象，其中所有的位都被设置为第一或第二个参数中的设置。这个操作符的存在是因为如果它不存在，那么bit-or <tt>操作符|</tt>的结果将是一个整数，当我们试图将其分配给CacheUpdateFlags类型的对象时，又会引发编译器警告。
   * @ref CacheUpdateFlags
   *
   */
  inline CacheUpdateFlags
  operator|(const CacheUpdateFlags f1, const CacheUpdateFlags f2)
  {
    return static_cast<CacheUpdateFlags>(static_cast<unsigned int>(f1) |
                                         static_cast<unsigned int>(f2));
  }

  /**
   * 全局操作符，它返回一个对象，其中所有的位都被设置了，而在参数中没有设置。这个操作符的存在是因为如果它不存在，那么位负的<tt>操作符~</tt>的结果将是一个整数，当我们试图将其分配给CacheUpdateFlags类型的对象时，会引发编译器警告。
   * @ref CacheUpdateFlags
   *
   */
  inline CacheUpdateFlags
  operator~(const CacheUpdateFlags f1)
  {
    return static_cast<CacheUpdateFlags>(static_cast<unsigned int>(f1) ^
                                         static_cast<unsigned int>(update_all));
  }



  /**
   * 全局操作符，将第二个参数的位也设置在第一个参数中。
   * @ref CacheUpdateFlags
   *
   */
  inline CacheUpdateFlags &
  operator|=(CacheUpdateFlags &f1, const CacheUpdateFlags f2)
  {
    f1 = f1 | f2;
    return f1;
  }


  /**
   * 全局操作符，它返回一个对象，其中所有位都被设置在第一个和第二个参数中。这个操作符的存在是因为如果它不存在，那么位和<tt>操作符&</tt>的结果将是一个整数，当我们试图将其分配给CacheUpdateFlags类型的对象时，会引发编译器警告。
   * @ref CacheUpdateFlags
   *
   */
  inline CacheUpdateFlags operator&(const CacheUpdateFlags f1,
                                    const CacheUpdateFlags f2)
  {
    return static_cast<CacheUpdateFlags>(static_cast<unsigned int>(f1) &
                                         static_cast<unsigned int>(f2));
  }


  /**
   * 全局操作符，如果第一个参数中的所有位没有在第二个参数中设置，则将其清除。
   * @ref CacheUpdateFlags
   *
   */
  inline CacheUpdateFlags &
  operator&=(CacheUpdateFlags &f1, const CacheUpdateFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }

} // namespace GridTools
DEAL_II_NAMESPACE_CLOSE

#endif


