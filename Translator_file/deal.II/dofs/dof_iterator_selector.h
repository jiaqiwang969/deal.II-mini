//include/deal.II-translator/dofs/dof_iterator_selector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_dof_iterators_h
#define dealii_dof_iterators_h

#include <deal.II/base/config.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int, int, int>
class DoFInvalidAccessor;

template <int structdim, int dim, int spacedim, bool lda>
class DoFAccessor;
template <int dim, int spacedim, bool lda>
class DoFCellAccessor;

template <int dim, int spacedim>
class DoFHandler;

template <typename Accessor>
class TriaRawIterator;
template <typename Accessor>
class TriaIterator;
template <typename Accessor>
class TriaActiveIterator;
#endif

namespace internal
{
  namespace DoFHandlerImplementation
  {
    template <int dim, int spacedim, bool lda = false>
    struct Iterators;


    /**
     * 在一个维度上定义一些用于DoF处理的类型。
     * 这些类型与
     * internal::TriangulationImplementation::Iterators<1,spacedim>,
     * 中声明的类型含义相同，只是对模板的处理要复杂一些。更多信息见
     * @ref Iterators 模块。
     *
     */
    template <int spacedim, bool lda>
    struct Iterators<1, spacedim, lda>
    {
      using CellAccessor = dealii::DoFCellAccessor<1, spacedim, lda>;
      using FaceAccessor = dealii::DoFAccessor<0, 1, spacedim, lda>;

      using raw_line_iterator    = TriaRawIterator<CellAccessor>;
      using line_iterator        = TriaIterator<CellAccessor>;
      using active_line_iterator = TriaActiveIterator<CellAccessor>;

      using raw_quad_iterator =
        TriaRawIterator<DoFInvalidAccessor<2, 1, spacedim>>;
      using quad_iterator = TriaIterator<DoFInvalidAccessor<2, 1, spacedim>>;
      using active_quad_iterator =
        TriaActiveIterator<DoFInvalidAccessor<2, 1, spacedim>>;

      using raw_hex_iterator =
        TriaRawIterator<DoFInvalidAccessor<3, 1, spacedim>>;
      using hex_iterator = TriaIterator<DoFInvalidAccessor<3, 1, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<DoFInvalidAccessor<3, 1, spacedim>>;

      using raw_cell_iterator    = raw_line_iterator;
      using cell_iterator        = line_iterator;
      using active_cell_iterator = active_line_iterator;

      using raw_face_iterator    = TriaRawIterator<FaceAccessor>;
      using face_iterator        = TriaIterator<FaceAccessor>;
      using active_face_iterator = TriaActiveIterator<FaceAccessor>;
    };



    /**
     * 为二维的DoF处理定义一些类型。
     * 这些类型与
     * internal::TriangulationImplementation::Iterators<2,spacedim>,
     * 中声明的类型具有相同的含义，只是对模板的处理要复杂一些。参见
     * @ref Iterators 模块以获得更多信息。
     *
     */
    template <int spacedim, bool lda>
    struct Iterators<2, spacedim, lda>
    {
      using CellAccessor = dealii::DoFCellAccessor<2, spacedim, lda>;
      using FaceAccessor = dealii::DoFAccessor<1, 2, spacedim, lda>;

      using raw_line_iterator    = TriaRawIterator<FaceAccessor>;
      using line_iterator        = TriaIterator<FaceAccessor>;
      using active_line_iterator = TriaActiveIterator<FaceAccessor>;

      using raw_quad_iterator    = TriaRawIterator<CellAccessor>;
      using quad_iterator        = TriaIterator<CellAccessor>;
      using active_quad_iterator = TriaActiveIterator<CellAccessor>;

      using raw_hex_iterator =
        TriaRawIterator<DoFInvalidAccessor<3, 2, spacedim>>;
      using hex_iterator = TriaIterator<DoFInvalidAccessor<3, 2, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<DoFInvalidAccessor<3, 2, spacedim>>;

      using raw_cell_iterator    = raw_quad_iterator;
      using cell_iterator        = quad_iterator;
      using active_cell_iterator = active_quad_iterator;

      using raw_face_iterator    = raw_line_iterator;
      using face_iterator        = line_iterator;
      using active_face_iterator = active_line_iterator;
    };



    /**
     * 为三维空间的DoF处理定义了一些类型。
     * 这些类型与
     * internal::TriangulationImplementation::Iterators<3,spacedim>,
     * 中声明的类型具有相同的含义，只是对模板的处理要复杂一些。更多信息见
     * @ref Iterators 模块。
     *
     */
    template <int spacedim, bool lda>
    struct Iterators<3, spacedim, lda>
    {
      using CellAccessor = dealii::DoFCellAccessor<3, spacedim, lda>;
      using FaceAccessor = dealii::DoFAccessor<2, 3, spacedim, lda>;

      using raw_line_iterator =
        TriaRawIterator<dealii::DoFAccessor<1, 3, spacedim, lda>>;
      using line_iterator =
        TriaIterator<dealii::DoFAccessor<1, 3, spacedim, lda>>;
      using active_line_iterator =
        TriaActiveIterator<dealii::DoFAccessor<1, 3, spacedim, lda>>;

      using raw_quad_iterator    = TriaRawIterator<FaceAccessor>;
      using quad_iterator        = TriaIterator<FaceAccessor>;
      using active_quad_iterator = TriaActiveIterator<FaceAccessor>;

      using raw_hex_iterator    = TriaRawIterator<CellAccessor>;
      using hex_iterator        = TriaIterator<CellAccessor>;
      using active_hex_iterator = TriaActiveIterator<CellAccessor>;

      using raw_cell_iterator    = raw_hex_iterator;
      using cell_iterator        = hex_iterator;
      using active_cell_iterator = active_hex_iterator;

      using raw_face_iterator    = raw_quad_iterator;
      using face_iterator        = quad_iterator;
      using active_face_iterator = active_quad_iterator;
    };
  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_dof_iterator_selector_h


