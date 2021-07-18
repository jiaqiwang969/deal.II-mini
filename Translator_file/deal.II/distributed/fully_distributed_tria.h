//include/deal.II-translator/distributed/fully_distributed_tria_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_fully_distributed_tria_h
#define dealii_fully_distributed_tria_h


#include <deal.II/base/config.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/grid_tools.h>

#include <vector>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
// forward declaration of the data type for periodic face pairs
namespace GridTools
{
  template <typename CellIterator>
  struct PeriodicFacePair;
}
#endif

namespace parallel
{
  /**
   * 一个用于全分布式三角测量的命名空间。
   * @ingroup parallel
   *
   */
  namespace fullydistributed
  {
    /**
     * 一个具有分布式粗网格的分布式三角测量。
     * parallel::fullydistributed::Triangulation
     * 的动机源于以下关于复杂几何体和/或由外部网格生成器创建的给定网格的观察。我们认为复杂的几何体是指只能用不可忽略的粗单元数量（>10,000）来进行网格划分的几何体。
     *
     *
     *
     *
     *
     *
     * - 从内存的角度来看，在每个进程上存储粗网格信息过于昂贵（如 parallel::distributed::Triangulation). 所做的那样 通常，一个进程只需要全局三角图的一小部分，即粗网格的一小部分，这样，粗网格的分区确实是必要的。每个进程上存储的单元由 @ref GlossLocallyOwnedCell "本地拥有的单元 "和 @ref GlossGhostCell "幽灵单元 "组成。
     *
     *
     *
     *
     *
     * - 活动细胞的分布
     *
     * - 在最精细的层面上
     *
     * - 通过简单地划分空间填充曲线在所有进程中的分布，对于源自大的粗网格的三角形来说，可能不会导致最佳结果：例如，属于同一进程的分区可能是不连续的，导致通信量增加（在一个节点内部和外部）。基于图的分区算法可能是 parallel::distributed::Triangulation. 所使用的空间填充曲线的一个合理的替代方案。为了能够构建一个完全分区的三角形，分配粗网格并给出关于分区的灵活性，需要以下成分。
     *
     *
     *
     *
     *
     * - 一个本地相关的粗网格三角图（顶点、单元格定义；包括一层幽灵单元格
     *
     * - 本地相关的粗网格三角图到全球粗网格三角图的映射
     *
     *
     *
     *
     * - 关于哪个单元应该被细化的信息，以及关于每个单元的子域_id、水平_子域_id、流形_id和边界_id的信息。        上面列出的成分被捆绑在结构 TriangulationDescription::Description. 中，用户必须填写这个数据结构
     *
     * - 在一个预处理步骤中
     *
     * - 在实际创建三角测量之前。创建 TriangulationDescription::Description 的预定义函数可以在命名空间 TriangulationDescription::Utilities. 中找到。一旦 TriangulationDescription::Description  `construction_data`被构建，就可以通过调用`tria.create_triangulation(construction_data);`来创建三角形`。
     * @note
     * 这个三角图支持。1D/2D/3D，悬挂节点，几何多网格，以及周期性。
     * @note
     * 你可以用create_triangulation()创建一个具有悬空节点和多网格层次的三角形。然而，一旦它被创建，就不能再被改变，也就是说，你不能在之后进行粗化或细化。
     * @note
     * 目前只有简单的周期性条件（即没有偏移量和旋转矩阵
     *
     * - 也可参见 GridTools::collect_periodic_faces()) 的文档。
     *
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public parallel::DistributedTriangulationBase<dim, spacedim>
    {
    public:
      using cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::cell_iterator;

      using active_cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::active_cell_iterator;

      using CellStatus =
        typename dealii::Triangulation<dim, spacedim>::CellStatus;

      /**
       * 构造函数。              @param  mpi_communicator
       * 用于三角测量的MPI通信器。
       *
       */
      explicit Triangulation(const MPI_Comm &mpi_communicator);

      /**
       * 解构器。
       *
       */
      virtual ~Triangulation() = default;

      /**
       * @copydoc   dealii::Triangulation::create_triangulation()
       * dealii::Triangulation::create_triangulation() .
       * @note  这是用于代替 Triangulation::create_triangulation()
       * 的函数，用于处理.II的一些其他三角计算。
       *
       */
      void
      create_triangulation(
        const TriangulationDescription::Description<dim, spacedim>
          &construction_data) override;

      /**
       * @note
       * 这个函数没有为这个类实现，并抛出一个断言。相反，使用其他create_triangulation()函数来创建三角形。
       *
       */
      virtual void
      create_triangulation(const std::vector<Point<spacedim>> &      vertices,
                           const std::vector<dealii::CellData<dim>> &cells,
                           const SubCellData &subcelldata) override;

      /**
       * 实现与基类中相同的函数。              @param  other_tria
       * 要复制的三角图。它可以是一个序列三角形，也可以是一个
       * parallel::distributed::Triangulation.
       * 两者都可以已经被完善。
       * @note  这个函数使用用set_partitioner()注册的分区器。
       *
       */
      void
      copy_triangulation(
        const dealii::Triangulation<dim, spacedim> &other_tria) override;

      /**
       * 注册一个分区器，在copy_triangulation方法中使用。
       * @param  partitioner
       * 一个分区函数，它的输入参数是对要分区的三角形的引用和要创建的分区的数量。
       * 该函数需要为给定三角形的每个活动单元设置子域ID，其值在零（包括）和函数的第二个参数（排他）之间。
       * @param  设置 参见设置枚举器的描述。
       * @note 默认情况下， GridTools::partition_triangulation_zorder()
       * 被用作分区器，不设置多网格层次上的数据结构。
       *
       */
      void
      set_partitioner(
        const std::function<void(dealii::Triangulation<dim, spacedim> &,
                                 const unsigned int)> &partitioner,
        const TriangulationDescription::Settings &     settings);

      /**
       * 根据设置的细化和粗化标志，粗化和细化网格。
       * @note 还没有实现。
       *
       */
      virtual void
      execute_coarsening_and_refinement() override;

      /**
       * 覆盖基类中prepare_coarsening_and_refinement的实现。
       * @note  还没有实现。
       *
       */
      virtual bool
      prepare_coarsening_and_refinement() override;

      /**
       * 如果三角形有悬空节点，返回true。
       * @note  还没有实现。
       *
       */
      virtual bool
      has_hanging_nodes() const override;

      /**
       * 返回本地内存消耗的字节数。
       *
       */
      virtual std::size_t
      memory_consumption() const override;

      virtual bool
      is_multilevel_hierarchy_constructed() const override;

      /**
       * 将三角图保存到给定的文件中。这个文件需要在一个共享的网络文件系统上，从计算中的所有节点都可以到达。参见SolutionTransfer类，了解如何将解决方案向量存储到这个文件中。其他基于单元的数据可以用register_data_attach()来保存。
       *
       */
      virtual void
      save(const std::string &filename) const override;

      /**
       * 将用save()保存的三角剖面加载回来。在调用这个函数之前，网格必须是空的。
       * 你需要用与保存时相同数量的MPI进程来加载，因此自动分区功能被禁用。
       * 用register_data_attach()保存的基于单元的数据可以在调用load()后用notify_ready_to_unpack()读入。
       *
       */
      virtual void
      load(const std::string &filename,
           const bool         autopartition = false) override;

    private:
      virtual unsigned int
      coarse_cell_id_to_coarse_cell_index(
        const types::coarse_cell_id coarse_cell_id) const override;

      virtual types::coarse_cell_id
      coarse_cell_index_to_coarse_cell_id(
        const unsigned int coarse_cell_index) const override;

      /**
       * 遍历所有属于本地的活动单元，并在私有成员向量local_cell_relations中记录它们将如何变化。
       * 由于该类目前不支持自适应网格细化，所有的单元将被标记为CellStatus
       * CELL_PERSIST。      这些关系目前只用于序列化。
       * 存储的向量的大小将等于本地拥有的活动单元的数量，并将按照这些单元的出现率排序。
       *
       */
      virtual void
      update_cell_relations() override;

      /**
       * 存储设置。
       *
       */
      TriangulationDescription::Settings settings;

      /**
       * copy_triangulation()中使用的分区器。
       *
       */
      std::function<void(dealii::Triangulation<dim, spacedim> &,
                         const unsigned int)>
        partitioner;

      /**
       * 成对的粗略单元ID及其索引的排序列表。
       *
       */
      std::vector<std::pair<types::coarse_cell_id, unsigned int>>
        coarse_cell_id_to_coarse_cell_index_vector;

      /**
       * 每个粗放单元的粗放单元ID的列表（存储在cell->index()）。
       *
       */
      std::vector<types::coarse_cell_id>
        coarse_cell_index_to_coarse_cell_id_vector;

      /**
       * 表示函数create_triangulation()被调用为内部使用的布尔值。
       *
       */
      bool currently_processing_create_triangulation_for_internal_usage;

      /**
       * 表示函数prepare_coarsening_and_refinement()被调用用于内部使用的布尔值。
       *
       */
      bool
        currently_processing_prepare_coarsening_and_refinement_for_internal_usage;
    };

  } // namespace fullydistributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif


