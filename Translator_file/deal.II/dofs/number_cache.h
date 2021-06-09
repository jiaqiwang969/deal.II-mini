//include/deal.II-translator/dofs/number_cache_0.txt
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

#ifndef dealii_number_cache_h
#define dealii_number_cache_h

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandlerImplementation
  {
    /**
     * 一个由DoFHandler类使用的结构，用于存储它们所处理的自由度信息。
     *
     */
    struct NumberCache
    {
      /**
       * 默认构造函数。
       *
       */
      NumberCache();

      /**
       * 复制构造函数。简单地将引用对象的所有成员复制到当前对象中。
       *
       */
      NumberCache(const NumberCache &) = default;

      /**
       * 移动构造函数。简单地将被引用对象的所有成员移到当前对象中。
       *
       */
      NumberCache(NumberCache &&) = default;

      /**
       * 创建一个NumberCache对象，对应于一个顺序DoFHandler对象，其中一个处理器存储所有自由度。(这里，"顺序
       * "意味着要么整个程序不使用MPI，要么使用MPI但只使用一个MPI进程，要么有多个MPI进程，但这个DoFHandler建立的Triangulation只在一个MPI进程上工作。)
       *
       */
      NumberCache(const types::global_dof_index n_global_dofs);


      /**
       * 创建一个NumberCache对象，该对象对应于一个并行的DoFHandler对象，其处理器数量与给定参数的大小相同，其中每个处理器存储作为第一个参数传递的向量的相应元素所指示的自由度。第二个参数表示当前处理器在所有参与处理器中的等级，这样我们就可以设置
       * @p locally_owned_dofs 和 @p n_locally_owned_dofs 字段。
       * 所有其他由当前对象存储的字段都可以并且是由参数计算出来的。
       *
       */
      NumberCache(const std::vector<IndexSet> &locally_owned_dofs_per_processor,
                  const unsigned int           my_rank);

      /**
       * 复制操作符。简单地将被引用对象的所有成员复制到当前对象中。
       *
       */
      NumberCache &
      operator=(const NumberCache &) = default;

      /**
       * 移动赋值运算符。简单地将被引用对象的所有成员移到当前对象上。
       *
       */
      NumberCache &
      operator=(NumberCache &&) = default;

      /**
       * 确定这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 这个函数重置了所有的存储信息。
       *
       */
      void
      clear();

      /**
       * 返回 @p n_locally_owned_dofs_per_processor
       * 的表示，既是在它被设置的情况下（直接返回数组），也是在我们需要在所有处理器上积累一些信息的情况下。后一种情况涉及到全局通信，通常设置起来比较昂贵，因为它调用了MPI_Allgather。
       *
       */
      std::vector<types::global_dof_index>
      get_n_locally_owned_dofs_per_processor(
        const MPI_Comm &mpi_communicator) const;

      /**
       * 返回 @p locally_owned_dofs_per_processor
       * 的表示，既是在它被设置的情况下（直接返回IndexSet字段的数组），也是在我们需要在所有处理器上积累一些信息的情况下。后一种情况涉及到全局通信，通常设置起来比较昂贵，因为它调用了MPI_Allgather。
       *
       */
      std::vector<IndexSet>
      get_locally_owned_dofs_per_processor(
        const MPI_Comm &mpi_communicator) const;

      /**
       * 道夫总数，在所有可能参与这个网格的处理器上累积。
       *
       */
      types::global_dof_index n_global_dofs;

      /**
       * 这个MPI进程所拥有的道夫数量。如果这是一个连续的计算，那么这等于n_global_dofs。这里，"顺序
       * "意味着要么整个程序不使用MPI，要么使用MPI但只使用一个MPI进程，要么有多个MPI进程，但这个DoFHandler建立的三角计算只在一个MPI进程上工作）。
       *
       */
      types::global_dof_index n_locally_owned_dofs;

      /**
       * 一个索引集，表示本地拥有的DoF的集合。如果这是一个顺序计算，那么它包含整个[0,n_global_dofs]范围。这里，"顺序
       * "意味着要么整个程序不使用MPI，要么使用MPI但只使用一个MPI进程，要么有多个MPI进程，但这个DoFHandler建立的Triangulation只在一个MPI进程上工作）。
       *
       */
      IndexSet locally_owned_dofs;

      /**
       * 各个MPI进程所拥有的DoF的数量。如果这是一个连续的计算，那么该向量包含一个等于n_global_dofs的单元素。这里，"顺序
       * "意味着要么整个程序不使用MPI，要么使用MPI但只使用一个MPI进程，要么有多个MPI进程，但这个DoFHandler建立的三角计算只在一个MPI进程上工作）。
       *
       */
      std::vector<types::global_dof_index> n_locally_owned_dofs_per_processor;

      /**
       * 各个MPI进程所拥有的DoF。如果这是一个连续的DoFHandler，那么该向量有一个等于local_owned_dofs的单元素。(这里，"顺序
       * "意味着要么整个程序不使用MPI，要么使用MPI但只使用一个MPI进程，要么有多个MPI进程，但这个DoFHandler建立的Triangulation只在一个MPI进程上工作。)
       *
       */
      std::vector<IndexSet> locally_owned_dofs_per_processor;

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
    NumberCache::serialize(Archive &ar, const unsigned int  /*version*/ )
    {
      ar &n_global_dofs &n_locally_owned_dofs;
      ar &               locally_owned_dofs;
      ar &               n_locally_owned_dofs_per_processor;
      ar &               locally_owned_dofs_per_processor;
    }

  } // namespace DoFHandlerImplementation
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_dof_iterator_selector_h


