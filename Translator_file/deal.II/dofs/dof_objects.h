//include/deal.II-translator/dofs/dof_objects_0.txt
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

#ifndef dealii_dof_objects_h
#define dealii_dof_objects_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int, int>
class DoFHandler;
#endif

namespace internal
{
  namespace DoFHandlerImplementation
  {
#ifndef DOXYGEN
    template <int>
    class DoFLevel;
    template <int>
    class DoFFaces;
#endif

    /**
     * 存储位于维度为 @p dim. <h3>Information for all DoFObjects
     * classes</h3>的对象上的自由度的索引。DoFObjects类存储了某一层面上每个单元的自由度的全局索引。一个自由度的全局索引或数字是解向量中相应值的零基索引，以及该层的全局矩阵或多网格矩阵的行和列索引。这些指数指的是无约束的向量和矩阵，其中我们没有考虑到悬挂节点引入的约束。
     * 由于顶点不与某一特定层次相关联，与顶点相关的指数不存储在DoFObjects类中，而是存储在
     * DoFHandler::vertex_dofs 数组中。
     * DoFObjects类不被直接使用，但这些类的对象被包含在DoFLevel和DoFFaces类中。
     * @ingroup dofs
     *
     */
    template <int dim>
    class DoFObjects
    {
    public:
      /**
       * 存储自由度的全局指数。
       *
       */
      std::vector<types::global_dof_index> dofs;

    public:
      /**
       * 将位于编号为 @p obj_index 的对象上的 @p local_index-th
       * 自由度的全局索引设置为最后一个参数所给的值。参数
       * @p dof_handler
       * 用于访问要用来计算存储该数据的位置的有限元。
       * 第三个参数， @p fe_index,
       * 必须等于零。否则它是未使用的，但我们保留这个参数，以便我们可以为非hp-和hp-有限元方法使用相同的接口，实际上使hp-和非hp-类之间共享DoFAccessor类的层次结构成为可能。
       *
       */
      template <int dh_dim, int spacedim>
      void
      set_dof_index(const dealii::DoFHandler<dh_dim, spacedim> &dof_handler,
                    const unsigned int                          obj_index,
                    const unsigned int                          fe_index,
                    const unsigned int                          local_index,
                    const types::global_dof_index               global_index);

      /**
       * 返回位于编号为 @p obj_index. 的对象上的 @p local_index-th
       * 自由度的全局索引。  @p dof_handler
       * 参数用于访问要用来计算该数据存储位置的有限元。
       * 第三个参数， @p fe_index,
       * 必须等于零。否则它将不被使用，但我们保留这个参数，以便我们可以为非hp-和hp-有限元方法使用相同的接口，实际上使hp-和非hp-类之间共享DoFAccessor类的层次结构成为可能。
       *
       */
      template <int dh_dim, int spacedim>
      types::global_dof_index
      get_dof_index(const dealii::DoFHandler<dh_dim, spacedim> &dof_handler,
                    const unsigned int                          obj_index,
                    const unsigned int                          fe_index,
                    const unsigned int local_index) const;

      /**
       * 返回值为1。这个函数的含义通过查看类
       * internal::hp::DoFObjects
       * 中的相应函数的内容而变得清晰。
       *
       */
      template <int dh_dim, int spacedim>
      unsigned int
      n_active_fe_indices(
        const dealii::DoFHandler<dh_dim, spacedim> &dof_handler,
        const types::global_dof_index               index) const;

      /**
       * 与上面的函数类似。断言给定的索引是零，然后返回真。
       *
       */
      template <int dh_dim, int spacedim>
      bool
      fe_index_is_active(
        const dealii::DoFHandler<dh_dim, spacedim> &dof_handler,
        const types::global_dof_index               index,
        const unsigned int                          fe_index) const;

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

      // Declare the classes that store levels and faces of DoFs friends so
      // that they can resize arrays.
      template <int>
      friend class DoFLevel;
      template <int>
      friend class DoFFaces;
    };


    // --------------------- template and inline functions ------------------

    template <int dim>
    template <int dh_dim, int spacedim>
    inline unsigned int
    DoFObjects<dim>::n_active_fe_indices(
      const dealii::DoFHandler<dh_dim, spacedim> &,
      const types::global_dof_index) const
    {
      return 1;
    }



    template <int dim>
    template <int dh_dim, int spacedim>
    inline bool
    DoFObjects<dim>::fe_index_is_active(
      const dealii::DoFHandler<dh_dim, spacedim> &,
      const types::global_dof_index,
      const unsigned int fe_index) const
    {
      (void)fe_index;
      Assert((fe_index ==
              dealii::DoFHandler<dh_dim, spacedim>::default_fe_index),
             ExcMessage("Only zero fe_index values are allowed for "
                        "non-hp-DoFHandlers."));
      return true;
    }



    template <int dim>
    template <int dh_dim, int spacedim>
    inline types::global_dof_index
    DoFObjects<dim>::get_dof_index(
      const dealii::DoFHandler<dh_dim, spacedim> &dof_handler,
      const unsigned int                          obj_index,
      const unsigned int                          fe_index,
      const unsigned int                          local_index) const
    {
      (void)fe_index;
      Assert(
        (fe_index == dealii::DoFHandler<dh_dim, spacedim>::default_fe_index),
        ExcMessage(
          "Only the default FE index is allowed for non-hp-DoFHandler objects"));
      Assert(
        local_index < dof_handler.get_fe().template n_dofs_per_object<dim>(),
        ExcIndexRange(local_index,
                      0,
                      dof_handler.get_fe().template n_dofs_per_object<dim>()));
      Assert(obj_index *
                   dof_handler.get_fe().template n_dofs_per_object<dim>() +
                 local_index <
               dofs.size(),
             ExcInternalError());

      return dofs[obj_index *
                    dof_handler.get_fe().template n_dofs_per_object<dim>() +
                  local_index];
    }


    template <int dim>
    template <class Archive>
    void
    DoFObjects<dim>::serialize(Archive &ar, const unsigned int)
    {
      ar &dofs;
    }

  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


