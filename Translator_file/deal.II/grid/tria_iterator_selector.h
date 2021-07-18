//include/deal.II-translator/grid/tria_iterator_selector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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

#ifndef dealii_tria_iterator_selector_h
#define dealii_tria_iterator_selector_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class CellAccessor;
template <int, int, int>
class InvalidAccessor;
template <int, int, int>
class TriaAccessor;
template <int dim, int spacedim>
class TriaAccessor<0, dim, spacedim>;
template <typename Accessor>
class TriaRawIterator;
template <typename Accessor>
class TriaIterator;
template <typename Accessor>
class TriaActiveIterator;
#endif

namespace internal
{
  namespace TriangulationImplementation
  {
    template <int dim, int spacedim>
    struct Iterators;

    /**
     * 这个类实现了一些不同维度的类型。
     * 这些是只针对一维情况的声明。参见 @ref Iterators
     * 模块以了解更多信息。        一个 @p line_iterator
     * 别名为在<tt>Triangulation<1></tt>对象的 @p
     * 行成员变量上操作的一个迭代器。一个 @p
     * active_line_iterator只对活动线进行操作。  @p
     * raw_line_iterator对象对所有的线进行操作，无论是否使用。
     * 由于我们是在一维中，所以声明了以下身份。
     * @code
     *  using raw_cell_iterator = raw_line_iterator;
     *  using cell_iterator = line_iterator;
     *  using active_cell_iterator = active_line_iterator;
     * @endcode
     * 为了能够在<tt>Triangulation<1></tt>中声明 @p begin_quad 等，
     * @p quad_iterators
     * 被声明为InvalidAccessor上的迭代器。因此，这些类型存在，但没有用，而且肯定会使任何非自愿的使用变得明显。这同样适用于六面体迭代器。
     * 这同样适用于 @p face_iterator
     * 类型，因为线条除了顶点之外没有任何子结构，不过这是以不同的方式处理的。
     *
     */
    template <int spacedim>
    struct Iterators<1, spacedim>
    {
      using raw_line_iterator =
        TriaRawIterator<dealii::CellAccessor<1, spacedim>>;
      using line_iterator = TriaIterator<dealii::CellAccessor<1, spacedim>>;
      using active_line_iterator =
        TriaActiveIterator<dealii::CellAccessor<1, spacedim>>;

      using raw_quad_iterator =
        TriaRawIterator<dealii::InvalidAccessor<2, 1, spacedim>>;
      using quad_iterator =
        TriaIterator<dealii::InvalidAccessor<2, 1, spacedim>>;
      using active_quad_iterator =
        TriaActiveIterator<dealii::InvalidAccessor<2, 1, spacedim>>;

      using raw_hex_iterator =
        TriaRawIterator<dealii::InvalidAccessor<3, 1, spacedim>>;
      using hex_iterator =
        TriaIterator<dealii::InvalidAccessor<3, 1, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<dealii::InvalidAccessor<3, 1, spacedim>>;
    };



    /**
     * 这个类实现了一些不同维度的类型。
     * 这些是只针对二维情况的声明。参见 @ref Iterators
     * 模块以了解更多信息。        一个 @p line_iterator
     * 别名为在<tt>Triangulation<2></tt>对象的 @p
     * 行成员变量上操作的一个迭代器。一个 @p
     * active_line_iterator只对活动线进行操作。  @p
     * raw_line_iterator对象对所有的线进行操作，无论是否使用。使用
     * @p
     * active_line_iterator可能在2D中不是特别有用，因为它只对未精炼的线进行操作。然而，如果相邻的单元格比现在的单元格多精炼了一次，那么精炼过的线也可能与未精炼过的单元格绑定。
     * 与线条迭代器类似， @p quad_iterator,  @p raw_quad_iterator 和
     * @p active_quad_iterator 也被声明。
     * 为了能够在<tt>Triangulation<[12]></tt>中声明 @p begin_hex
     * 等， @p hex_iterators
     * 被声明为InvalidAccessor的迭代器。因此，这些类型存在，但没有用，而且肯定会使任何非自愿的使用变得明显。
     * 由于我们是在二维空间，所以声明了以下身份。
     * @code
     *  using raw_cell_iterator = raw_quad_iterator;
     *  using cell_iterator = quad_iterator;
     *  using active_cell_iterator = active_quad_iterator;
     *
     *  using raw_face_iterator = raw_line_iterator;
     *  using face_iterator = line_iterator;
     *  using active_face_iterator = active_line_iterator;
     * @endcode
     *
     *
     */
    template <int spacedim>
    struct Iterators<2, spacedim>
    {
      using raw_line_iterator =
        TriaRawIterator<dealii::TriaAccessor<1, 2, spacedim>>;
      using line_iterator = TriaIterator<dealii::TriaAccessor<1, 2, spacedim>>;
      using active_line_iterator =
        TriaActiveIterator<dealii::TriaAccessor<1, 2, spacedim>>;

      using raw_quad_iterator =
        TriaRawIterator<dealii::CellAccessor<2, spacedim>>;
      using quad_iterator = TriaIterator<dealii::CellAccessor<2, spacedim>>;
      using active_quad_iterator =
        TriaActiveIterator<dealii::CellAccessor<2, spacedim>>;

      using raw_hex_iterator =
        TriaRawIterator<dealii::InvalidAccessor<3, 2, spacedim>>;
      using hex_iterator =
        TriaIterator<dealii::InvalidAccessor<3, 2, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<dealii::InvalidAccessor<3, 2, spacedim>>;
    };


    /**
     * 这个类实现了一些不同维度的类型。
     * 这些是只针对三维情况的声明。更多信息见 @ref
     * Iterators 模块。
     * 对于数据类型的声明，或多或少与低维度的声明相同（见<tt>Iterators<[12]></tt>）。特定维度的数据类型在这里，因为我们是在三个维度。
     * @code
     *  using raw_cell_iterator = raw_hex_iterator;
     *  using cell_iterator = hex_iterator;
     *  using active_cell_iterator = active_hex_iterator;
     *
     *  using raw_face_iterator = raw_quad_iterator;
     *  using face_iterator = quad_iterator;
     *  using active_face_iterator = active_quad_iterator;
     * @endcode
     *
     *
     */
    template <int spacedim>
    struct Iterators<3, spacedim>
    {
      using raw_line_iterator =
        TriaRawIterator<dealii::TriaAccessor<1, 3, spacedim>>;
      using line_iterator = TriaIterator<dealii::TriaAccessor<1, 3, spacedim>>;
      using active_line_iterator =
        TriaActiveIterator<dealii::TriaAccessor<1, 3, spacedim>>;

      using raw_quad_iterator =
        TriaRawIterator<dealii::TriaAccessor<2, 3, spacedim>>;
      using quad_iterator = TriaIterator<dealii::TriaAccessor<2, 3, spacedim>>;
      using active_quad_iterator =
        TriaActiveIterator<dealii::TriaAccessor<2, 3, spacedim>>;

      using raw_hex_iterator =
        TriaRawIterator<dealii::CellAccessor<3, spacedim>>;
      using hex_iterator = TriaIterator<dealii::CellAccessor<3, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<dealii::CellAccessor<3, spacedim>>;
    };

  } // namespace TriangulationImplementation

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_tria_iterator_selector_h


