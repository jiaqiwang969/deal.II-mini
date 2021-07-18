//include/deal.II-translator/dofs/dof_levels_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_dof_levels_h
#define dealii_dof_levels_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <deal.II/dofs/dof_objects.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DoFHandlerImplementation
  {
    /**
     * 用于存储细胞自由度信息的结构，按级别组织。
     * 我们在#cell_dof_indices_cache中存储了每个单元的自由度指数的缓存值，因为这是一个经常被要求的操作。这些值由
     * DoFCellAccessor::update_cell_dof_indices_cache 设置，并由
     * DoFCellAccessor::get_dof_indices.
     * 使用。请注意，顶点是独立的，事实上与单元无关。因此，位于顶点上的自由度指数不存储在这里，而是存储在
     * dealii::DoFHandler 类的成员变量中。
     * 位于低维物体上的自由度指数，即二维的线上和三维的四边形和线上的自由度指数的处理与单元上的类似。然而，这些几何对象，作为一种概括，被称为面，并不是以层次结构的方式组织的。因此，位于这些对象上的自由度被存储在单独的类中，即<tt>DoFFaces</tt>类。
     * 对该对象的访问通常是通过 DoFAccessor::set_dof_index() 和
     * DoFAccessor::dof_index()
     * 函数或派生类的类似函数，而派生类又使用
     * DoFHandler::get_dof_index()
     * 和相应的setter函数访问成员变量。因此，实际数据格式的知识被封装到目前的类的层次结构以及
     * dealii::DoFHandler 类中。
     *
     */
    template <int dim>
    class DoFLevel
    {
    public:
      /**
       * 用于单元格上的DoF指数的缓存。这个数组的大小等于某一层的单元格数量乘以selected_fe.n_dofs_per_cell()。
       *
       */
      std::vector<types::global_dof_index> cell_dof_indices_cache;

      /**
       * 包含dof-indices和相关访问函数的对象。
       *
       */
      DoFObjects<dim> dof_object;

      /**
       * 返回一个指针，指向给定单元的DoF指数缓存的开头。
       * @param  obj_index 我们正在查看的单元的编号。
       * @param  dofs_per_cell 该单元的每个DoFs的数量。
       * @return
       * 指向当前单元的第一个DoF索引的指针。接下来的
       * dofs_per_cell 指数是针对当前单元的。
       *
       */
      const types::global_dof_index *
      get_cell_cache_start(const unsigned int obj_index,
                           const unsigned int dofs_per_cell) const;

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



    template <int dim>
    inline const types::global_dof_index *
    DoFLevel<dim>::get_cell_cache_start(const unsigned int obj_index,
                                        const unsigned int dofs_per_cell) const
    {
      Assert(obj_index * dofs_per_cell + dofs_per_cell <=
               cell_dof_indices_cache.size(),
             ExcInternalError());

      return cell_dof_indices_cache.data() + (obj_index * dofs_per_cell);
    }



    template <int dim>
    inline std::size_t
    DoFLevel<dim>::memory_consumption() const
    {
      return (MemoryConsumption::memory_consumption(cell_dof_indices_cache) +
              MemoryConsumption::memory_consumption(dof_object));
    }


    template <int dim>
    template <class Archive>
    inline void
    DoFLevel<dim>::serialize(Archive &ar, const unsigned int)
    {
      ar &cell_dof_indices_cache;
      ar &dof_object;
    }
  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


