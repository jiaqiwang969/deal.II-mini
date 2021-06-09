//include/deal.II-translator/dofs/dof_faces_0.txt
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

#ifndef dealii_dof_faces_h
#define dealii_dof_faces_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/dofs/dof_objects.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * 用于DoFHandler类组的内部数据结构的命名空间。
   * @ingroup dofs
   *
   */
  namespace DoFHandlerImplementation
  {
    /**
     * <h4>DoFFaces</h4>
     * 这些类与DoFLevel类类似。我们在这里存储的是与面孔相关的信息，而不是细胞，因为这些信息与细胞的层次结构无关，而细胞是以层次来组织的。在二维中我们存储位于线上的自由度信息，而在三维中我们存储位于四边形和线上的自由度信息。在一维中我们什么都不做，因为线的面是顶点，被单独处理。
     * 除了DoFObjects对象包含要存储的数据（自由度指数）外，我们不存储任何数据或提供任何功能。然而，我们确实实现了一个函数来确定所含DoFObjects对象的内存消耗估计值。
     * 所包含的数据通常不被直接访问。相反，除了来自DoFHandler类的一些访问，通常是通过
     * DoFAccessor::set_dof_index() 和 DoFAccessor::dof_index()
     * 函数或派生类的类似函数来访问，而派生类又使用
     * DoFHandler::get_dof_index()
     * 和相应的setter函数访问成员变量。因此，实际数据格式的知识被封装到本层次的类以及
     * dealii::DoFHandler 类中。
     *
     */
    template <int dim>
    class DoFFaces
    {
    public:
      /**
       * 构造函数。这个构造函数被删除，以防止使用这个模板，因为只应该使用特殊化。
       *
       */
      DoFFaces() = delete;
    };

    /**
     * 存储1D中面的自由度指数。因为这些将是顶点，被单独处理，所以不要做任何事情。
     *
     */
    template <>
    class DoFFaces<1>
    {
    public:
      /**
       * 确定这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };

    /**
     * 存储2D中面的自由度的指数，这些自由度是线。
     *
     */
    template <>
    class DoFFaces<2>
    {
    public:
      /**
       * 含有线上自由度数据的对象。
       *
       */
      internal::DoFHandlerImplementation::DoFObjects<1> lines;

      /**
       * 确定这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };

    /**
     * 存储三维中面的自由度指数，这些面是四边形的，另外也存储在线上。
     *
     */
    template <>
    class DoFFaces<3>
    {
    public:
      /**
       * 包含线上自由度数据的对象。
       *
       */
      internal::DoFHandlerImplementation::DoFObjects<1> lines;

      /**
       * 包含四边形上的DoFs数据的对象。
       *
       */
      internal::DoFHandlerImplementation::DoFObjects<2> quads;

      /**
       * 确定此对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };



    template <class Archive>
    void
    DoFFaces<1>::serialize(Archive &, const unsigned int)
    {}


    template <class Archive>
    void
    DoFFaces<2>::serialize(Archive &ar, const unsigned int)
    {
      ar &lines;
    }


    template <class Archive>
    void
    DoFFaces<3>::serialize(Archive &ar, const unsigned int)
    {
      ar &lines &quads;
    }

  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


