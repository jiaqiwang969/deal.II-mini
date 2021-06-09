//include/deal.II-translator/distributed/cell_data_transfer_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#ifndef dealii_distributed_cell_data_transfer_h
#define dealii_distributed_cell_data_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/adaptation_strategies.h>

#include <boost/range/iterator_range.hpp>

#include <functional>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    /**
     * 在细化和/或粗化分布式三角测量时，传输与每个活动单元相关的数据（如错误指示器）并处理必要的通信。
     * 因此，该类对单元相关信息的作用就像
     * parallel::distributed::SolutionTransfer 对定义在
     * parallel::distributed::Triangulation.
     * 上的自由度值的作用一样。该类被设计用来操作任何一种可序列化的数据类型。必须提供一个非分布式的容器（如Vector或
     * `std::vector`)
     * ），它以活动单元被遍历的相同顺序来保存单元的数据。换句话说，每个条目都对应于具有相同索引
     * CellAccessor::active_cell_index(), 的单元，容器的大小必须是
     * Triangulation::n_active_cells().  <h3>Transferring cell-wise data</h3>
     * 下面的代码片段演示了如何在细化/粗化注册三角的过程中传输单元相关数据。
     * @code
     * // prepare the triangulation,
     * triangulation.prepare_coarsening_and_refinement();
     *
     * // prepare the CellDataTransfer object for coarsening and refinement
     * // and give the cell data vector that we intend to unpack later,
     * Vector<double> data_to_transfer(triangulation.n_active_cells());
     * //[fill data_to_transfer with cell-wise values...]
     *
     * parallel::distributed::CellDataTransfer<dim, spacedim, Vector<double>>
     * cell_data_trans(triangulation);
     * cell_data_trans.prepare_for_coarsening_and_refinement(data_to_transfer);
     *
     * // actually execute the refinement,
     * triangulation.execute_coarsening_and_refinement();
     *
     * // unpack transferred data,
     * Vector<double> transferred_data(triangulation.n_active_cells());
     * cell_data_trans.unpack(transferred_data);
     *
     * @endcode
     * <h3>Use for serialization</h3>
     * 这个类可以用来序列化和随后反序列化一个带有附加数据的分布式网格到不同的文件。
     * 对于序列化，下面的代码片段不仅保存了三角图本身，而且还保存了附加的单元数据。
     * @code
     * Vector<double> data_to_transfer(triangulation.n_active_cells());
     * //[fill data_to_transfer with cell-wise values...]
     *
     * parallel::distributed::CellDataTransfer<dim, spacedim, Vector<double>>
     * cell_data_trans(triangulation);
     * cell_data_trans.prepare_for_serialization(data_to_transfer);
     *
     * triangulation.save(filename);
     * @endcode
     * 之后，在反序列化过程中，三角图和数据都可以按以下方式恢复。
     * @code
     * //[create coarse mesh...]
     * triangulation.load(filename);
     *
     * parallel::distributed::CellDataTransfer<dim, spacedim, Vector<double>>
     * cell_data_trans(triangulation);
     * Vector<double> transferred_data(triangulation.n_active_cells());
     * cell_data_trans.deserialize(transferred_data);
     * @endcode
     *
     * @note  如果你通过
     * parallel::distributed::Triangulation::register_data_attach() 和
     * parallel::distributed::Triangulation::notify_ready_for_unpack()
     * 接口使用一个以上的对象来传输数据，目的是为了序列化，对相应的prepare_for_serialization()和deserialize()函数的调用需要分别以相同的顺序发生。依靠这个接口的类，例如
     * parallel::distributed::CellDataTransfer,
     * parallel::distributed::SolutionTransfer, 和
     * Particles::ParticleHandler.  。
     * @note  参见 parallel::distributed::SolutionTransfer
     * 的文档，了解传输和序列化的匹配代码片段。
     * @ingroup distributed
     *
     */
    template <int dim, int spacedim = dim, typename VectorType = Vector<double>>
    class CellDataTransfer
    {
    private:
      /**
       * 一个别名，定义了所提供的容器模板的数据类型。
       *
       */
      using value_type = typename VectorType::value_type;

    public:
      /**
       * 构造函数。              @param[in]  triangulation
       * 所有的操作都将发生在这个三角形上。当这个构造函数被调用时，有关的细化还没有发生。
       * @param[in]  transfer_variable_size_data
       * 指定你的VectorType容器是否存储有不同大小的值。不同数量的数据可能会被打包到每个单元格中，例如，如果VectorType容器的底层ValueType是一个容器本身。
       * @param[in]  refinement_strategy
       * %函数决定数据将如何从其父单元格存储在精炼单元格上。
       * @param[in]  coarsening_strategy
       * %函数决定在一个单元格上存储哪些数据，其子单元格将被粗化为。
       *
       */
      CellDataTransfer(
        const parallel::distributed::Triangulation<dim, spacedim>
          &                               triangulation,
        const bool                        transfer_variable_size_data = false,
        const std::function<std::vector<value_type>(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator
            &              parent,
          const value_type parent_value)> refinement_strategy =
          &dealii::AdaptationStrategies::Refinement::
            preserve<dim, spacedim, value_type>,
        const std::function<value_type(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator
            &                            parent,
          const std::vector<value_type> &children_values)> coarsening_strategy =
          &dealii::AdaptationStrategies::Coarsening::
            check_equality<dim, spacedim, value_type>);

      /**
       * 为粗化和细化的当前对象做准备。
       * 它在底层三角形上注册了 @p in 的数据传输。        @p
       * in
       * 包括要插值到新的（细化和/或粗化的）网格上的数据。关于如何使用这个功能的更多信息，请参见这个类的文档。
       * 这个函数对指定的容器只能调用一次，直到数据传输完成。如果要通过这个类传输多个向量，请使用下面这个函数。
       *
       */
      void
      prepare_for_coarsening_and_refinement(const VectorType &in);

      /**
       * 与上面的函数相同，只是针对一个向量的列表。
       *
       */
      void
      prepare_for_coarsening_and_refinement(
        const std::vector<const VectorType *> &all_in);

      /**
       * 准备对给定的向量进行序列化。            序列化是由
       * Triangulation::save().
       * 完成的，关于如何使用这一功能的更多信息，请参见该类的文档。
       * 这个函数对指定的容器只能调用一次，直到数据传输完成。如果多个向量应通过这个类来传输，请使用下面的函数。
       *
       */
      void
      prepare_for_serialization(const VectorType &in);

      /**
       * 与上面的函数相同，只是针对一个向量的列表。
       *
       */
      void
      prepare_for_serialization(const std::vector<const VectorType *> &all_in);

      /**
       * 在网格被细化或粗化到当前的单元格集合之前，解压之前存储在这个对象中的信息。
       *
       */
      void
      unpack(VectorType &out);

      /**
       * 与上面的函数相同，只是针对一个向量列表。
       *
       */
      void
      unpack(std::vector<VectorType *> &all_out);

      /**
       * 执行存储信息的反序列化。      这需要在调用
       * Triangulation::load(). 后进行。
       *
       */
      void
      deserialize(VectorType &out);

      /**
       * 和上面的函数一样，只是针对一个向量的列表。
       *
       */
      void
      deserialize(std::vector<VectorType *> &all_out);

    private:
      /**
       * 指向要处理的三角结构的指针。
       *
       */
      SmartPointer<const parallel::distributed::Triangulation<dim, spacedim>,
                   CellDataTransfer<dim, spacedim, VectorType>>
        triangulation;

      /**
       * 指定要传输的数据的大小是否因单元格而异。
       *
       */
      const bool transfer_variable_size_data;

      /**
       * 决定数据如何从其父单元存储到细化单元的函数。
       *
       */
      const std::function<std::vector<value_type>(
        const typename Triangulation<dim, spacedim>::cell_iterator &parent,
        const value_type parent_value)>
        refinement_strategy;

      /**
       * 决定如何处理来自子单元的数据以存储在父单元上的函数。
       *
       */
      const std::function<value_type(
        const typename Triangulation<dim, spacedim>::cell_iterator &parent,
        const std::vector<value_type> &children_values)>
        coarsening_strategy;

      /**
       * 一个向量，存储所有我们应该从旧网格复制到新网格的向量的指针。
       *
       */
      std::vector<const VectorType *> input_vectors;

      /**
       * triangulation分配给这个对象的句柄，我们可以用它来访问我们的内存偏移和我们的打包函数。
       *
       */
      unsigned int handle;

      /**
       * 将pack_callback()函数注册给triangulation并存储返回的句柄。
       *
       */
      void
      register_data_attach();

      /**
       * 一个回调函数，用于将当前网格上的数据打包成对象，以后可以在细化、粗化和重新划分后取回。
       *
       */
      std::vector<char>
      pack_callback(const typename parallel::distributed::
                      Triangulation<dim, spacedim>::cell_iterator &cell,
                    const typename parallel::distributed::
                      Triangulation<dim, spacedim>::CellStatus status);

      /**
       * 一个回调函数，用来解开当前网格上的数据，这些数据在细化、粗化和重新划分之前，已经被打包在网格上了。
       *
       */
      void
      unpack_callback(
        const typename parallel::distributed::Triangulation<dim, spacedim>::
          cell_iterator &cell,
        const typename parallel::distributed::Triangulation<dim, spacedim>::
          CellStatus status,
        const boost::iterator_range<std::vector<char>::const_iterator>
          &                        data_range,
        std::vector<VectorType *> &all_out);
    };
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_distributed_cell_data_transfer_h */ 


