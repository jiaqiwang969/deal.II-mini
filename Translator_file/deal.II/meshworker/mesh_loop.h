//include/deal.II-translator/meshworker/mesh_loop_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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


#ifndef dealii_mesh_worker_mesh_loop_h
#define dealii_mesh_worker_mesh_loop_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/loop.h>

#include <functional>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename>
class TriaActiveIterator;
#endif

namespace MeshWorker
{
  namespace internal
  {
    /**
     * 一个帮助类，为底层单元格迭代器类型提供类型定义。
     *
     */
    template <class CellIteratorType>
    struct CellIteratorBaseType
    {
      /**
       * 用于单元格迭代器类型的类型定义。
       *
       */
      using type = CellIteratorType;
    };

    /**
     * 一个辅助类，为底层单元格迭代器类型提供类型定义。
     * 这个特殊化是针对IteratorRange的，它可以有TriaActiveIterator或FilteredIterator作为其基础类型。
     *
     */
    template <class CellIteratorType>
    struct CellIteratorBaseType<IteratorOverIterators<CellIteratorType>>
    {
      /**
       * 单元迭代器类型的类型定义。
       *
       */
      // Since we can have filtered iterators and the like as template
      // arguments, we recursivelyremove the template layers to retrieve the
      // underlying iterator type.
      using type = typename CellIteratorBaseType<CellIteratorType>::type;
    };

    /**
     * 一个帮助类，为底层单元格迭代器类型提供类型定义。
     * 这个特殊化是针对FilteredIterator的，它可以有一个TriaActiveIterator作为它的基础类型，或者可以与另一个FilteredIterator嵌套作为迭代的类型。
     *
     */
    template <class CellIteratorType>
    struct CellIteratorBaseType<FilteredIterator<CellIteratorType>>
    {
      /**
       * 单元格迭代器类型的定义。
       *
       */
      // Since we can have nested filtered iterators, we recursively
      // remove the template layers to retrieve the underlying iterator type.
      using type = typename CellIteratorBaseType<CellIteratorType>::type;
    };
  } // namespace internal

#ifdef DOXYGEN
  /**
   * 这个别名为Mesh_loop()中使用的单元格工作者的函数类型引入了一个友好而简短的名称。
   *
   */
  using CellWorkerFunctionType = std::function<
    void(const CellIteratorBaseType &, ScratchData &, CopyData &)>;

  /**
   * 这个别名为Mesh_loop()中使用的单元格工作者的函数类型引入了一个友好而简短的名字。
   *
   */
  using CopierFunctionType = std::function<void(const CopyData &)>;

  /**
   * 这个别名为Mesh_loop()中使用的边界工作者的函数类型引入了一个友好而简短的名字。
   *
   */
  using BoundaryWorkerFunctionType =
    std::function<void(const CellIteratorBaseType &,
                       const unsigned int,
                       ScratchData &,
                       CopyData &)>;

  /**
   * 这个别名为Mesh_loop()中使用的面片工作者的函数类型引入了一个友好而简短的名称。
   *
   */
  using FaceWorkerFunctionType =
    std::function<void(const CellIteratorBaseType &,
                       const unsigned int,
                       const unsigned int,
                       const CellIteratorBaseType &,
                       const unsigned int,
                       const unsigned int,
                       ScratchData &,
                       CopyData &)>;
#endif

  /**
   * 这个函数扩展了WorkStream的概念，可以在网格（单元和/或面）上工作，并处理自适应细化面的工作和并行计算的复杂逻辑（例如，对面的工作到幽灵邻居）。
   * @p mesh_loop
   * 可用于简化对单元（例如装配）、边界（诺伊曼型边界条件）或内部面（例如在非连续加尔金方法中）的操作。该函数在许多教程中使用，包括
   * step-12 ,  step-16 , 和 step-47 ，仅举几例。
   * 对于均匀细化的网格，使用 WorkStream::run() 和 @p cell_worker
   * 会比较容易， @p cell_worker
   * 也会在面上循环，并根据当前和相邻单元的情况来处理面的组合条件。所有进行这些循环的用户代码都需要手动插入逻辑，为当前单元的每个面确定相邻的单元，以及相邻单元上与当前面相对应的面的索引。
   * 如果启用了局部细化，并且当前或相邻单元有悬挂的节点，这就更复杂了。在这种情况下，还需要确定当前或邻近面的相应子面。
   * 这个方法将该逻辑外部化（独立于用户代码），并将面的装配条款（内部面、边界面或平行计算中不同子域ID之间的面）与单元上的装配分开，允许用户指定两个额外的工作者（一个
   * @p cell_worker, 一个 @p boundary_worker, 和一个 @p face_worker)
   * ，根据传递的特定AssembleFlags @p flags 在每个 @p cell,
   * 中自动调用。 @p cell_worker
   * 被传递给单元格标识符、ScratchData对象和CopyData对象，遵循与
   * WorkStream::run(). 相同的原则。
   * 在内部，该函数除了传递给 @p boundary_worker,
   * 之外，还传递了一个 @p face_no
   * 参数，用于识别应该被执行集成的面。而 @p face_worker
   * 则需要在单元格和相邻单元格上明确识别当前面，因此它被调用时有六个参数（每个单元格有三个参数：实际单元格、面的索引和子面_索引。如果不需要子面集成，那么除了通常的ScratchData和CopyData对象外，子面_索引为
   * numbers::invalid_unsigned_int) 。    如果传递了标志
   * AssembleFlags::assemble_own_cells
   * ，那么默认的行为是先在面上循环并在那里做功，然后再计算单元上的实际功。如果指定了标志
   * AssembleFlags::assemble_own_interior_faces_once
   * ，那么每个内部面只被访问一次， @p face_worker
   * 被假定为一次集成所有面的条款（并在不连续Galerkin设置中增加面的两边贡献）。
   * 当AssembleFlags只包含 @p assemble_own_cells, 时，该方法等同于
   * WorkStream::run() 方法，可以作为该方法的替代品使用。
   * 两个数据类型ScratchData和CopyData需要有一个工作拷贝构造函数。ScratchData只在worker函数中使用，而CopyData是由worker传递给copyer的对象。每次这个函数访问一个新的单元时，CopyData对象就会被重置为提供给这个函数的值（在这里它就会调用单元和面对工作者）。换句话说，在一个单元格上调用`copier`和在下一个单元格上调用`cell_worker`/`face_worker`/`boundary_worker`函数之间，没有任何状态的延续，用户代码不需要在单元格整合的开始或拷贝操作的结束时重置拷贝对象。在
   * "face_worker "或 "boundary_worker "内部重置 "copier
   * "的状态构成一个bug，并可能导致一些意外的结果。下面的例子显示了什么是不允许的，因为复印机有可能在一个单元上的许多面之间共享。
   * @code
   *
   * using ScratchData      = MeshWorker::ScratchData<dim, spacedim>;
   * using CopyData         = MeshWorker::CopyData<1, 1, 1>;
   * using CellIteratorType = decltype(dof_handler.begin_active());
   *
   * ScratchData            scratch(...);
   * CopyData               copy(...);
   *
   * std::function<void(const CellIteratorType &, ScratchData &, CopyData &)>
   * empty_cell_worker;
   *
   * auto boundary_worker = [...] (
   * const CellIteratorType &cell,
   * const unsigned int      face,
   * ScratchData            &scratch_data,
   * CopyData               &copy_data)
   * {
   * const auto &fe_face_values = scratch_data.reinit(cell, face);
   * copy_data = CopyData(...); // This is an error, as we lose the
   *                           // accumulation that has been performed on
   *                           // other boundary faces of the same cell.
   *
   * for (unsigned int q_point = 0;
   *     q_point < fe_face_values.n_quadrature_points;
   *     ++q_point)
   *  {
   *    copy_data.vectors[0][0] += 1.0 fe_face_values.JxW(q_point);
   *  }
   * };
   *
   * double value = 0;
   * auto copier = [...](const CopyData &copy_data)
   * {
   * value += copy_data.vectors[0][0]; // Contributions from some faces may
   *                                   // be missing.
   * };
   *
   * MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
   *                     empty_cell_worker, copier,
   *                     scratch, copy,
   *                     MeshWorker::assemble_boundary_faces,
   *                     boundary_worker);
   * @endcode
   * queue_length参数表示在任何给定时间内可以活用的项目的数量。每个项目由输入流的chunk_size元素组成，这些元素将被worker和copyer函数在同一个线程上一个接一个地处理。
   * 如果你的数据对象很大，或者它们的构造函数很昂贵，记住ScratchData对象的queue_length副本和CopyData对象的`queue_length*chunk_size`副本是有帮助的。
   * @note
   * 这里的Doxygen文档中显示的函数参数类型和默认值（空工作者函数）与实际类型相比略有简化。
   * @note  关于模板类型的要求以及 @p queue_length 和 @p
   * chunk_size
   * 的含义的更多信息可以在WorkStream命名空间及其成员的文档中找到。
   * @ingroup MeshWorker
   *
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class CellIteratorBaseType =
              typename internal::CellIteratorBaseType<CellIteratorType>::type>
  void
  mesh_loop(
#ifdef DOXYGEN
    const CellIteratorType &begin,
    const CellIteratorType &end,

    const CellWorkerFunctionType &cell_worker,
    const CopierType &            copier,

    const ScratchData &sample_scratch_data,
    const CopyData &   sample_copy_data,

    const AssembleFlags flags = assemble_own_cells,

    const BoundaryWorkerFunctionType &boundary_worker =
      BoundaryWorkerFunctionType(),

    const FaceWorkerFunctionType &face_worker = FaceWorkerFunctionType(),
    const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
    const unsigned int chunk_size   = 8
#else
    const CellIteratorType &                         begin,
    const typename identity<CellIteratorType>::type &end,

    const typename identity<std::function<
      void(const CellIteratorBaseType &, ScratchData &, CopyData &)>>::type
      &cell_worker,
    const typename identity<std::function<void(const CopyData &)>>::type
      &copier,

    const ScratchData &sample_scratch_data,
    const CopyData &   sample_copy_data,

    const AssembleFlags flags = assemble_own_cells,

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &boundary_worker = std::function<void(const CellIteratorBaseType &,
                                            const unsigned int,
                                            ScratchData &,
                                            CopyData &)>(),

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &face_worker = std::function<void(const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        ScratchData &,
                                        CopyData &)>(),

    const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
    const unsigned int chunk_size   = 8
#endif
  )
  {
    Assert(
      (!cell_worker) == !(flags & work_on_cells),
      ExcMessage(
        "If you provide a cell worker function, you also need to request "
        "that work should be done on cells by setting the 'work_on_cells' flag. "
        "Conversely, if you don't provide a cell worker function, you "
        "cannot set the 'work_on_cells' flag. One of these two "
        "conditions is not satisfied."));

    Assert((flags & (assemble_own_interior_faces_once |
                     assemble_own_interior_faces_both)) !=
             (assemble_own_interior_faces_once |
              assemble_own_interior_faces_both),
           ExcMessage(
             "If you provide a face worker function, you also need to request "
             "that work should be done on interior faces by setting either the "
             "'assemble_own_interior_faces_once' flag or the "
             "'assemble_own_interior_faces_both' flag. "
             "Conversely, if you don't provide a face worker function, you "
             "cannot set either of these two flags. One of these two "
             "conditions is not satisfied."));

    Assert((flags & (assemble_ghost_faces_once | assemble_ghost_faces_both)) !=
             (assemble_ghost_faces_once | assemble_ghost_faces_both),
           ExcMessage(
             "You can only 'specify assemble_ghost_faces_once' "
             "OR 'assemble_ghost_faces_both', but not both of these flags."));

    Assert(
      !(flags & cells_after_faces) ||
        (flags & (assemble_own_cells | assemble_ghost_cells)),
      ExcMessage(
        "The option 'cells_after_faces' only makes sense if you assemble on cells."));

    Assert(
      (!face_worker) == !(flags & work_on_faces),
      ExcMessage(
        "If you provide a face worker function, you also need to request "
        "that work should be done on faces by setting the 'work_on_faces' flag. "
        "Conversely, if you don't provide a face worker function, you "
        "cannot set the 'work_on_faces' flag. One of these two "
        "conditions is not satisfied."));

    Assert(
      (!boundary_worker) == !(flags & assemble_boundary_faces),
      ExcMessage(
        "If you provide a boundary face worker function, you also need to request "
        "that work should be done on boundary faces by setting the 'assemble_boundary_faces' flag. "
        "Conversely, if you don't provide a boundary face worker function, you "
        "cannot set the 'assemble_boundary_faces' flag. One of these two "
        "conditions is not satisfied."));

    auto cell_action = [&](const CellIteratorBaseType &cell,
                           ScratchData &               scratch,
                           CopyData &                  copy) {
      // First reset the CopyData class to the empty copy_data given by the
      // user.
      copy = sample_copy_data;

      // Store the dimension in which we are working for later use
      const auto dim = cell->get_triangulation().dimension;

      const bool ignore_subdomain =
        (cell->get_triangulation().locally_owned_subdomain() ==
         numbers::invalid_subdomain_id);

      types::subdomain_id current_subdomain_id =
        (cell->is_level_cell() ? cell->level_subdomain_id() :
                                 cell->subdomain_id());

      const bool own_cell =
        ignore_subdomain ||
        (current_subdomain_id ==
         cell->get_triangulation().locally_owned_subdomain());

      if ((!ignore_subdomain) &&
          (current_subdomain_id == numbers::artificial_subdomain_id))
        return;

      if (!(flags & (cells_after_faces)) &&
          (((flags & (assemble_own_cells)) && own_cell) ||
           ((flags & assemble_ghost_cells) && !own_cell)))
        cell_worker(cell, scratch, copy);

      if (flags & (work_on_faces | work_on_boundary))
        for (const unsigned int face_no : cell->face_indices())
          {
            if (cell->at_boundary(face_no) &&
                !cell->has_periodic_neighbor(face_no))
              {
                // only integrate boundary faces of own cells
                if ((flags & assemble_boundary_faces) && own_cell)
                  boundary_worker(cell, face_no, scratch, copy);
              }
            else
              {
                // interior face, potentially assemble
                TriaIterator<typename CellIteratorBaseType::AccessorType>
                  neighbor = cell->neighbor_or_periodic_neighbor(face_no);

                types::subdomain_id neighbor_subdomain_id =
                  numbers::artificial_subdomain_id;
                if (neighbor->is_level_cell())
                  neighbor_subdomain_id = neighbor->level_subdomain_id();
                // subdomain id is only valid for active cells
                else if (neighbor->is_active())
                  neighbor_subdomain_id = neighbor->subdomain_id();

                const bool own_neighbor =
                  ignore_subdomain ||
                  (neighbor_subdomain_id ==
                   cell->get_triangulation().locally_owned_subdomain());

                // skip all faces between two ghost cells
                if (!own_cell && !own_neighbor)
                  continue;

                // skip if the user doesn't want faces between own cells
                if (own_cell && own_neighbor &&
                    !(flags & (assemble_own_interior_faces_both |
                               assemble_own_interior_faces_once)))
                  continue;

                // skip face to ghost
                if (own_cell != own_neighbor &&
                    !(flags &
                      (assemble_ghost_faces_both | assemble_ghost_faces_once)))
                  continue;

                // Deal with refinement edges from the refined side. Assuming
                // one-irregular meshes, this situation should only occur if
                // both cells are active.
                const bool periodic_neighbor =
                  cell->has_periodic_neighbor(face_no);

                if (dim > 1 && ((!periodic_neighbor &&
                                 cell->neighbor_is_coarser(face_no) &&
                                 neighbor->is_active()) ||
                                (periodic_neighbor &&
                                 cell->periodic_neighbor_is_coarser(face_no) &&
                                 neighbor->is_active())))
                  {
                    Assert(cell->is_active(), ExcInternalError());

                    // skip if only one processor needs to assemble the face
                    // to a ghost cell and the fine cell is not ours.
                    if (!own_cell && (flags & assemble_ghost_faces_once))
                      continue;

                    const std::pair<unsigned int, unsigned int>
                      neighbor_face_no =
                        periodic_neighbor ?
                          cell->periodic_neighbor_of_coarser_periodic_neighbor(
                            face_no) :
                          cell->neighbor_of_coarser_neighbor(face_no);

                    face_worker(cell,
                                face_no,
                                numbers::invalid_unsigned_int,
                                neighbor,
                                neighbor_face_no.first,
                                neighbor_face_no.second,
                                scratch,
                                copy);

                    if (flags & assemble_own_interior_faces_both)
                      {
                        // If own faces are to be assembled from both sides,
                        // call the faceworker again with swapped arguments.
                        // This is because we won't be looking at an adaptively
                        // refined edge coming from the other side.
                        face_worker(neighbor,
                                    neighbor_face_no.first,
                                    neighbor_face_no.second,
                                    cell,
                                    face_no,
                                    numbers::invalid_unsigned_int,
                                    scratch,
                                    copy);
                      }
                  }
                else if (dim == 1 && cell->level() > neighbor->level())
                  {
                    // In one dimension, there is no other check to do
                    const unsigned int neighbor_face_no =
                      periodic_neighbor ?
                        cell->periodic_neighbor_face_no(face_no) :
                        cell->neighbor_face_no(face_no);
                    Assert(periodic_neighbor ||
                             neighbor->face(neighbor_face_no) ==
                               cell->face(face_no),
                           ExcInternalError());

                    face_worker(cell,
                                face_no,
                                numbers::invalid_unsigned_int,
                                neighbor,
                                neighbor_face_no,
                                numbers::invalid_unsigned_int,
                                scratch,
                                copy);

                    if (flags & assemble_own_interior_faces_both)
                      {
                        // If own faces are to be assembled from both sides,
                        // call the faceworker again with swapped arguments.
                        face_worker(neighbor,
                                    neighbor_face_no,
                                    numbers::invalid_unsigned_int,
                                    cell,
                                    face_no,
                                    numbers::invalid_unsigned_int,
                                    scratch,
                                    copy);
                      }
                  }
                else
                  {
                    // If iterator is active and neighbor is refined, skip
                    // internal face.
                    if (dealii::internal::is_active_iterator(cell) &&
                        neighbor->has_children())
                      continue;

                    // Now neighbor is on the same refinement level.
                    // Double check.
                    Assert(!cell->neighbor_is_coarser(face_no),
                           ExcInternalError());

                    // If we own both cells only do faces from one side (unless
                    // AssembleFlags says otherwise). Here, we rely on cell
                    // comparison that will look at cell->index().
                    if (own_cell && own_neighbor &&
                        (flags & assemble_own_interior_faces_once) &&
                        (neighbor < cell))
                      continue;

                    // We only look at faces to ghost on the same level once
                    // (only where own_cell=true and own_neighbor=false)
                    if (!own_cell)
                      continue;

                    // now only one processor assembles faces_to_ghost. We let
                    // the processor with the smaller (level-)subdomain id
                    // assemble the face.
                    if (own_cell && !own_neighbor &&
                        (flags & assemble_ghost_faces_once) &&
                        (neighbor_subdomain_id < current_subdomain_id))
                      continue;

                    const unsigned int neighbor_face_no =
                      periodic_neighbor ?
                        cell->periodic_neighbor_face_no(face_no) :
                        cell->neighbor_face_no(face_no);
                    Assert(periodic_neighbor ||
                             neighbor->face(neighbor_face_no) ==
                               cell->face(face_no),
                           ExcInternalError());

                    face_worker(cell,
                                face_no,
                                numbers::invalid_unsigned_int,
                                neighbor,
                                neighbor_face_no,
                                numbers::invalid_unsigned_int,
                                scratch,
                                copy);
                  }
              }
          } // faces

      // Execute the cell_worker if faces are handled before cells
      if ((flags & cells_after_faces) &&
          (((flags & assemble_own_cells) && own_cell) ||
           ((flags & assemble_ghost_cells) && !own_cell)))
        cell_worker(cell, scratch, copy);
    };

    // Submit to workstream
    WorkStream::run(begin,
                    end,
                    cell_action,
                    copier,
                    sample_scratch_data,
                    sample_copy_data,
                    queue_length,
                    chunk_size);
  }

  /**
   * 与上面的函数相同，但适用于迭代器范围（因此也适用于过滤的迭代器）。
   * 在串行情况下，该函数的一个使用实例是这样给出的
   * @code
   *
   * using ScratchData      = MeshWorker::ScratchData<dim, spacedim>;
   * using CopyData         = MeshWorker::CopyData<1, 1, 1>;
   * using CellIteratorType = decltype(dof_handler.begin_active());
   *
   * ScratchData            scratch(...);
   * CopyData               copy(...);
   *
   * auto cell_worker = [...] (
   * const CellIteratorType &cell,
   * ScratchData            &scratch_data,
   * CopyData               &copy_data)
   * {
   * ...
   * };
   *
   * auto copier = [...](const CopyData &copy_data)
   * {
   * ...
   * };
   *
   * MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
   *                     cell_worker, copier,
   *                     scratch, copy,
   *                     MeshWorker::assemble_own_cells);
   * @endcode
   * 该函数在并行分布式情况下的使用示例，其中复制器只在本地拥有的单元上调用，由以下内容给出
   * @code
   *
   * using ScratchData      = MeshWorker::ScratchData<dim, spacedim>;
   * using CopyData         = MeshWorker::CopyData<1, 1, 1>;
   * using CellIteratorType = decltype(dof_handler.begin_active());
   *
   * ScratchData            scratch(...);
   * CopyData               copy(...);
   *
   * auto cell_worker = [...] (
   * const CellIteratorType &cell,
   * ScratchData            &scratch_data,
   * CopyData               &copy_data)
   * {
   * ...
   * };
   *
   * auto copier = [...](const CopyData &copy_data)
   * {
   * ...
   * };
   *
   * const auto filtered_iterator_range =
   * filter_iterators(dof_handler.active_cell_iterators(),
   *                  IteratorFilters::LocallyOwnedCell());
   *
   * MeshWorker::mesh_loop(filtered_iterator_range,
   *                     cell_worker, copier,
   *                     scratch, copy,
   *                     MeshWorker::assemble_own_cells);
   * @endcode
   *
   * @ingroup MeshWorker
   *
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class CellIteratorBaseType =
              typename internal::CellIteratorBaseType<CellIteratorType>::type>
  void
  mesh_loop(
    IteratorRange<CellIteratorType> iterator_range,
    const typename identity<std::function<
      void(const CellIteratorBaseType &, ScratchData &, CopyData &)>>::type
      &cell_worker,
    const typename identity<std::function<void(const CopyData &)>>::type
      &copier,

    const ScratchData &sample_scratch_data,
    const CopyData &   sample_copy_data,

    const AssembleFlags flags = assemble_own_cells,

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &boundary_worker = std::function<void(const CellIteratorBaseType &,
                                            const unsigned int,
                                            ScratchData &,
                                            CopyData &)>(),

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &face_worker = std::function<void(const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        ScratchData &,
                                        CopyData &)>(),

    const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
    const unsigned int chunk_size   = 8)
  {
    // Call the function above
    mesh_loop<typename IteratorRange<CellIteratorType>::IteratorOverIterators,
              ScratchData,
              CopyData,
              CellIteratorBaseType>(iterator_range.begin(),
                                    iterator_range.end(),
                                    cell_worker,
                                    copier,
                                    sample_scratch_data,
                                    sample_copy_data,
                                    flags,
                                    boundary_worker,
                                    face_worker,
                                    queue_length,
                                    chunk_size);
  }

  /**
   * 这是Mesh_loop()函数的一个变体，可用于作为类成员函数的worker和copier函数。
   * 作为 @p end 传递的参数必须可以转换为与 @p
   * 开始相同的类型，但本身不一定是同一类型。这允许编写类似<code>mesh_loop(dof_handler.begin_active(),
   * dof_handler.end(), ...)</code>的代码，其中第一个是
   * DoFHandler::active_cell_iterator 类型，而第二个是
   * DoFHandler::raw_cell_iterator. 类型。  @p queue_length
   * 参数表示在任何特定时间可以活着的项目数量。每个项目由输入流的
   * @p chunk_size
   * 个元素组成，这些元素将由工作者和复制者函数在同一线程上一个接一个地处理。
   * @note
   * 如果你的数据对象很大，或者它们的构造函数很昂贵，记住<tt>queue_length</tt>拷贝的<tt>ScratchData</tt>对象和<tt>queue_length*chunk_size</tt>拷贝的<tt>CopyData</tt>对象是很有用的。
   * 该函数的一个使用例子是这样给出的
   * @code
   *
   * struct ScratchData;
   * struct CopyData;
   *
   * template <int dim, int spacedim>
   * class MyClass
   * {
   * public:
   * void
   * cell_worker(const CellIteratorType &cell, ScratchData &, CopyData &);
   *
   * void
   * copier(const CopyData &);
   *
   * ...
   * };
   *
   * ...
   *
   * MyClass<dim, spacedim> my_class;
   * ScratchData            scratch;
   * CopyData               copy;
   *
   * mesh_loop(tria.begin_active(),
   *         tria.end(),
   *         my_class,
   *         &MyClass<dim, spacedim>::cell_worker,
   *         &MyClass<dim, spacedim>::copier,
   *         scratch,
   *         copy,
   *         assemble_own_cells);
   * @endcode
   * @ingroup MeshWorker
   *
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class MainClass>
  void
  mesh_loop(const CellIteratorType &                         begin,
            const typename identity<CellIteratorType>::type &end,
            MainClass &                                      main_class,
            void (MainClass::*cell_worker)(const CellIteratorType &,
                                           ScratchData &,
                                           CopyData &),
            void (MainClass::*copier)(const CopyData &),
            const ScratchData & sample_scratch_data,
            const CopyData &    sample_copy_data,
            const AssembleFlags flags                      = assemble_own_cells,
            void (MainClass::*boundary_worker)(const CellIteratorType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &) = nullptr,
            void (MainClass::*face_worker)(const CellIteratorType &,
                                           const unsigned int,
                                           const unsigned int,
                                           const CellIteratorType &,
                                           const unsigned int,
                                           const unsigned int,
                                           ScratchData &,
                                           CopyData &)     = nullptr,
            const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
            const unsigned int chunk_size   = 8)
  {
    std::function<void(const CellIteratorType &, ScratchData &, CopyData &)>
      f_cell_worker;

    std::function<void(
      const CellIteratorType &, const unsigned int, ScratchData &, CopyData &)>
      f_boundary_worker;

    std::function<void(const CellIteratorType &,
                       const unsigned int,
                       const unsigned int,
                       const CellIteratorType &,
                       const unsigned int,
                       const unsigned int,
                       ScratchData &,
                       CopyData &)>
      f_face_worker;

    if (cell_worker != nullptr)
      f_cell_worker = [&main_class,
                       cell_worker](const CellIteratorType &cell_iterator,
                                    ScratchData &           scratch_data,
                                    CopyData &              copy_data) {
        (main_class.*cell_worker)(cell_iterator, scratch_data, copy_data);
      };

    if (boundary_worker != nullptr)
      f_boundary_worker =
        [&main_class, boundary_worker](const CellIteratorType &cell_iterator,
                                       const unsigned int      face_no,
                                       ScratchData &           scratch_data,
                                       CopyData &              copy_data) {
          (main_class.*
           boundary_worker)(cell_iterator, face_no, scratch_data, copy_data);
        };

    if (face_worker != nullptr)
      f_face_worker = [&main_class,
                       face_worker](const CellIteratorType &cell_iterator_1,
                                    const unsigned int      face_index_1,
                                    const unsigned int      subface_index_1,
                                    const CellIteratorType &cell_iterator_2,
                                    const unsigned int      face_index_2,
                                    const unsigned int      subface_index_2,
                                    ScratchData &           scratch_data,
                                    CopyData &              copy_data) {
        (main_class.*face_worker)(cell_iterator_1,
                                  face_index_1,
                                  subface_index_1,
                                  cell_iterator_2,
                                  face_index_2,
                                  subface_index_2,
                                  scratch_data,
                                  copy_data);
      };

    mesh_loop(begin,
              end,
              f_cell_worker,
              [&main_class, copier](const CopyData &copy_data) {
                (main_class.*copier)(copy_data);
              },
              sample_scratch_data,
              sample_copy_data,
              flags,
              f_boundary_worker,
              f_face_worker,
              queue_length,
              chunk_size);
  }

  /**
   * 和上面的函数一样，但是对于迭代器的范围（以及，因此，过滤的迭代器）。
   * 在串行情况下，该函数的一个使用实例是这样给出的
   * @code
   *
   * struct ScratchData;
   * struct CopyData;
   *
   * template <int dim, int spacedim>
   * class MyClass
   * {
   * public:
   * void
   * cell_worker(const CellIteratorType &cell, ScratchData &, CopyData &);
   *
   * void
   * copier(const CopyData &);
   *
   * ...
   * };
   *
   * ...
   *
   * MyClass<dim, spacedim> my_class;
   * ScratchData            scratch;
   * CopyData               copy;
   *
   * mesh_loop(tria.active_cell_iterators(),
   *         my_class,
   *         &MyClass<dim, spacedim>::cell_worker,
   *         &MyClass<dim, spacedim>::copier,
   *         scratch,
   *         copy,
   *         assemble_own_cells);
   * @endcode
   * 该函数在并行分布式情况下的使用示例，其中复制器只在本地拥有的单元上调用，由以下内容给出
   * @code
   *
   * struct ScratchData;
   * struct CopyData;
   *
   * template <int dim, int spacedim>
   * class MyClass
   * {
   * public:
   * void
   * cell_worker(const CellIteratorType &cell, ScratchData &, CopyData &);
   *
   * void
   * copier(const CopyData &);
   *
   * ...
   * };
   *
   * ...
   *
   * MyClass<dim, spacedim> my_class;
   * ScratchData            scratch;
   * CopyData               copy;
   *
   * const auto filtered_iterator_range =
   * filter_iterators(distributed_tria.active_cell_iterators(),
   *                  IteratorFilters::LocallyOwnedCell());
   *
   * mesh_loop(filtered_iterator_range,
   *         my_class,
   *         &MyClass<dim, spacedim>::cell_worker,
   *         &MyClass<dim, spacedim>::copier,
   *         scratch,
   *         copy,
   *         assemble_own_cells);
   * @endcode
   *
   * @ingroup MeshWorker
   *
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class MainClass,
            class CellIteratorBaseType =
              typename internal::CellIteratorBaseType<CellIteratorType>::type>
  void
  mesh_loop(IteratorRange<CellIteratorType> iterator_range,
            MainClass &                     main_class,
            void (MainClass::*cell_worker)(const CellIteratorBaseType &,
                                           ScratchData &,
                                           CopyData &),
            void (MainClass::*copier)(const CopyData &),
            const ScratchData & sample_scratch_data,
            const CopyData &    sample_copy_data,
            const AssembleFlags flags                      = assemble_own_cells,
            void (MainClass::*boundary_worker)(const CellIteratorBaseType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &) = nullptr,
            void (MainClass::*face_worker)(const CellIteratorBaseType &,
                                           const unsigned int,
                                           const unsigned int,
                                           const CellIteratorBaseType &,
                                           const unsigned int,
                                           const unsigned int,
                                           ScratchData &,
                                           CopyData &)     = nullptr,
            const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
            const unsigned int chunk_size   = 8)
  {
    // Call the function above
    mesh_loop<typename IteratorRange<CellIteratorType>::IteratorOverIterators,
              ScratchData,
              CopyData,
              MainClass,
              CellIteratorBaseType>(iterator_range.begin(),
                                    iterator_range.end(),
                                    main_class,
                                    cell_worker,
                                    copier,
                                    sample_scratch_data,
                                    sample_copy_data,
                                    flags,
                                    boundary_worker,
                                    face_worker,
                                    queue_length,
                                    chunk_size);
  }
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


