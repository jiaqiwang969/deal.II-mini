//include/deal.II-translator/grid/tria_faces_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2021 by the deal.II authors
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

#ifndef dealii_tria_faces_h
#define dealii_tria_faces_h

#include <deal.II/base/config.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_objects.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * 这个类包含了属于三角形的面的信息。这些类类似于TriaLevel类。由于细胞是以层次结构组织起来的，每个三角形由几个这样的TriaLevels组成。然而，三角形的面，低维物体，如二维的线或三维的线和四边形，不一定要基于这样的层次结构。事实上，如果我们想实现各向异性的细化，我们必须将它们组织在一个对象中。因此，TriaFaces类将属于三角形的面的信息与TriaLevel类分开存储。
     *
     */
    class TriaFaces
    {
    public:
      /**
       * 构造函数。
       *
       */
      TriaFaces(const unsigned int dim);

      /**
       * Boost::serialization. 的默认构造函数。
       *
       */
      TriaFaces() = default;

      /**
       * 底层三角结构的尺寸。
       *
       */
      unsigned int dim;

      /**
       * 包含四边形数据的TriaObject。
       * @note 只用于dim=3。
       *
       */
      TriaObjects quads;

      /**
       * 每个四边形的每一行的方向。
       * @note 仅用于dim=3。
       *
       */
      std::vector<unsigned char> quads_line_orientations;

      /**
       * 每个四边形的参考单元类型。
       * @note  仅用于dim=3。
       *
       */
      std::vector<dealii::ReferenceCell> quad_reference_cell;

      /**
       * 包含线的数据的TriaObject。
       * @note  仅用于dim>1。
       *
       */
      TriaObjects lines;

      /**
       * 确定此对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，将此对象的数据读入或写入一个流中。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };



    template <class Archive>
    void
    TriaFaces::serialize(Archive &ar, const unsigned int)
    {
      ar &dim;

      if (dim == 2)
        ar &lines;

      if (dim == 3)
        ar &quads &lines &quads_line_orientations &quad_reference_cell;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


