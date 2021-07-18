//include/deal.II-translator/distributed/shared_tria_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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

#ifndef dealii_distributed_shared_tria_h
#define dealii_distributed_shared_tria_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/tria.h>

#include <functional>
#include <list>
#include <set>
#include <utility>
#include <vector>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
#ifdef DEAL_II_WITH_MPI


  namespace shared
  {
    /**
     * 这个类提供了一个并行的三角形，每个处理器都知道全局网格的每个单元（与
     * parallel::distributed::Triangulation
     * 类不同），但在用MPI运行时，单元被自动分割，这样每个处理器就
     * "拥有 "一个单元子集。这个类的使用在  step-18
     * 中得到了证明。        与 parallel::distributed::Triangulation
     * 和 parallel::fullydistributed::Triangulation
     * 类不同，这意味着整个网格被存储在每个处理器上。虽然这显然是一个内存瓶颈，将这个类的使用限制在几十个或几百个MPI进程中，但网格的划分可以用来划分参与处理器之间的工作，如装配或后处理，它也可以用来划分哪个处理器存储矩阵和向量的哪一部分。因此，使用这个类通常是比更多的
     * parallel::distributed::Triangulation
     * 类更温和的代码并行化介绍，在后者中，处理器只知道他们自己的那部分网格，但对其他处理器拥有的单元却一无所知，除了他们自己那部分域周围的单层幽灵单元之外。
     * 作为在每个处理器上存储整个网格的结果，活动单元需要在所有处理器上被标记为细化或粗化，如果你想调整它们，无论它们被归类为本地所有、幽灵或人工。
     * 在计算时间和内存的考虑决定了程序需要并行运行，但在算法上需要每个处理器都知道整个网格的情况下，该类也很有用。一个例子是，一个应用程序必须同时拥有体积和表面网格，然后都可以独立分区，但很难确保本地拥有的表面网格单元集与本地拥有的体积网格单元集相邻，反之亦然。在这种情况下，知道两个网格的<i>entirety</i>可以确保耦合项的装配可以实现，而不需要实施过于复杂的方案来在处理器之间传递相邻单元的信息。
     * 单元在处理器之间的划分是根据一些不同的可能性在内部完成的。通过向该类的构造函数传递适当的标志（见
     * parallel::shared::Triangulation::Settings
     * 枚举），可以选择不同的网格划分方式，包括由应用决定的方式，而不是由希望最小化处理器所拥有的子域之间的接口长度（如METIS和Zoltan包所做的，两者都是划分的选项）。DoFHandler类知道如何以适合于分区网格的方式列举自由度。
     * @ingroup distributed
     *
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public dealii::parallel::TriangulationBase<dim, spacedim>
    {
    public:
      using active_cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::active_cell_iterator;
      using cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::cell_iterator;

      /**
       * 分布式三角计算的配置标志要在构造函数中设置。设置可以用位法OR来组合。
       * 构造函数要求 <code>partition_auto</code>,
       * <code>partition_metis</code> 、 <code>partition_zorder</code>,
       * <code>partition_zoltan</code> 和
       * <code>partition_custom_signal</code>
       * 中正好有一个被设置。如果 <code>partition_auto</code>
       * 被选中，它将使用 <code>partition_zoltan</code>
       * （如果可用），然后是 <code>partition_metis</code>
       * （如果可用），最后是 <code>partition_zorder</code>  。
       *
       */
      enum Settings
      {
        /**
         * 根据配置deal.II时发现的启用的依赖关系，选择分区器。
         * 特别是，如果发现了Trilinos软件包Zoltan，那么就使用
         * @p partition_zoltan
         * 策略。如果没有找到Zoltan，但找到了METIS包，那么就使用partition_metis策略。如果这两个都没有找到，那么就使用partition_zorder分区策略。
         *
         */
        partition_auto = 0x0,

        /**
         * 使用METIS分区器对活动单元进行分区。
         *
         */
        partition_metis = 0x1,

        /**
         *
         */
        partition_zorder = 0x2,

        /**
         * 使用Zoltan来划分活动单元。
         *
         */
        partition_zoltan = 0x3,

        /**
         * 使用一个自定义的、用户定义的函数来划分单元格。这可以通过在三角图首次创建时将post_refinement信号连接到三角图，并通过信号传递用户定义的函数来实现
         * <code>std::bind</code>  。        下面是一个例子。
         * @code
         * template <int dim>
         * void mypartition(parallel::shared::Triangulation<dim> &tria)
         * {
         * // user defined partitioning scheme: assign subdomain_ids
         * // round-robin in a mostly random way:
         * std::vector<unsigned int> assignment =
         *   {0,0,1,2,0,0,2,1,0,2,2,1,2,2,0,0};
         * unsigned int index = 0;
         * for (const auto &cell : tria.active_cell_iterators())
         *   cell->set_subdomain_id(assignment[(index++)%16]);
         * }
         *
         * int main ()
         * {
         * parallel::shared::Triangulation<dim> tria(
         *   ...,
         *   parallel::shared::Triangulation<dim>::partition_custom_signal);
         * tria.signals.post_refinement.connect(std::bind(&mypartition<dim>,
         *                                      std::ref(tria)));
         * }
         * @endcode
         * 一个使用lambda函数的等效代码看起来是这样的。
         * @code
         * int main ()
         * {
         * parallel::shared::Triangulation<dim> tria(
         *   ...,
         *   parallel::shared::Triangulation<dim>::partition_custom_signal);
         * tria.signals.post_refinement.connect (
         *   [&tria]()
         *   {
         *     // user defined partitioning scheme as above
         *     ...
         *   });
         * }
         * @endcode
         *
         * @note
         * 如果你打算使用几何多网格的自定义分区，除了活动单元外，你必须手动分区水平单元。
         *
         */
        partition_custom_signal = 0x4,

        /**
         * 这个标志需要被设置以使用几何多网格功能。这个选项需要额外的计算和通信。
         * 注意：这个标志应该总是和活动单元划分方法的标志一起设置。
         *
         */
        construct_multigrid_hierarchy = 0x8,
      };


      /**
       * 构造器。            标志 @p allow_artificial_cells
       * 可以用来启用人工细胞。如果启用，这个类的行为与
       * parallel::distributed::Triangulation 和
       * parallel::fullydistributed::Triangulation
       * 类似，即会有本地拥有的细胞、单层的幽灵细胞和人工细胞。然而，我们不应该忘记，与那些平行三角形相比，所有的细胞都是在所有的过程中重复的，在大多数情况下，导致人工细胞明显增多。
       * 如果人工细胞被禁用，所有非本地拥有的细胞都被视为幽灵细胞。这可能会导致非常昂贵的幽灵值更新步骤。虽然在人工单元的情况下，幽灵值更新只导致与直接进程邻居的点对点通信，但如果没有人工单元，这些就会退化为每个进程与其他每个进程的通信（"全对全
       * "通信）。如果这样的幽灵值更新是你代码中的瓶颈，你可能想考虑启用人工单元。
       *
       */
      Triangulation(
        const MPI_Comm &mpi_communicator,
        const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing =
          (dealii::Triangulation<dim, spacedim>::none),
        const bool     allow_artificial_cells = false,
        const Settings settings               = partition_auto);

      /**
       * 销毁器。
       *
       */
      virtual ~Triangulation() override = default;

      /**
       * 如果支持多级层次结构并且已经构建完成，则返回。
       *
       */
      virtual bool
      is_multilevel_hierarchy_constructed() const override;

      /**
       * 根据设置的细化和粗化标志，对网格进行粗化和细化。
       * 这一步等同于 dealii::Triangulation 类，在最后增加调用
       * dealii::GridTools::partition_triangulation() 。
       *
       */
      virtual void
      execute_coarsening_and_refinement() override;

      /**
       * 创建一个三角剖面。
       * 这个函数也是根据提供给构造函数的MPI通信器来划分三角形的。
       *
       */
      virtual void
      create_triangulation(const std::vector<Point<spacedim>> &vertices,
                           const std::vector<CellData<dim>> &  cells,
                           const SubCellData &subcelldata) override;

      /**
       * @copydoc   Triangulation::create_triangulation() .
       * @note  还没有实现。
       *
       */
      virtual void
      create_triangulation(
        const TriangulationDescription::Description<dim, spacedim>
          &construction_data) override;

      /**
       * 将 @p other_tria 复制到这个三角形中。
       * 这个函数也是根据提供给构造器的MPI通信器来划分三角结构。
       * @note 这个函数不能与 parallel::distributed::Triangulation,
       * 一起使用，因为它只存储它拥有的那些单元，在它本地拥有的单元周围的一层幽灵单元，以及一些人造单元。
       *
       */
      virtual void
      copy_triangulation(
        const dealii::Triangulation<dim, spacedim> &other_tria) override;

      /**
       * 为了序列化的目的，从一个流中读取此对象的数据。扔掉之前的内容。
       * 这个函数首先做与 dealii::Triangulation::load,
       * 相同的工作，然后根据提供给构造函数的MPI通信器来划分三角形。
       *
       */
      template <class Archive>
      void
      load(Archive &ar, const unsigned int version);

      /**
       * 返回一个长度为 Triangulation::n_active_cells()
       * 的向量，其中每个元素都存储了这个单元的所有者的子域ID。向量的元素显然与本地拥有的和幽灵单元的子域id相同，但对于那些不在其子域_id字段中存储单元所有者是谁的人工单元也是正确的。
       *
       */
      const std::vector<types::subdomain_id> &
      get_true_subdomain_ids_of_cells() const;

      /**
       * 返回一个长度为 Triangulation::n_cells(level)
       * 的向量，其中每个元素都存储了这个单元的所有者的水平子域id。该向量的元素显然与本地拥有的和幽灵单元的水平子域ID相同，但对于那些不在其水平子域_id字段中存储单元的所有者的人工单元也是正确的。
       *
       */
      const std::vector<types::subdomain_id> &
      get_true_level_subdomain_ids_of_cells(const unsigned int level) const;

      /**
       * 返回allow_artificial_cells
       * ，即如果允许人工细胞，则为true。
       *
       */
      bool
      with_artificial_cells() const;

    private:
      /**
       * 设置
       *
       */
      const Settings settings;

      /**
       * 一个决定是否允许人造细胞的标志。
       *
       */
      const bool allow_artificial_cells;

      /**
       * 这个函数调用 GridTools::partition_triangulation
       * ()，如果在类的构造函数中要求，就会标记人工细胞。
       *
       */
      void
      partition();

      /**
       * 一个包含子域ID的向量，这些子域是通过使用zorder、METIS或用户定义的分区方案获得的细胞。
       * 在allow_artificial_cells为false的情况下，这个向量与存储在三角形类的cell->subdomain_id()的ID一致。当allow_artificial_cells为true时，属于人工的单元将有cell->subdomain_id()
       * ==  numbers::artificial;
       * 原始分区信息被存储，以允许使用半人工的单元的顺序DoF分布和分区函数。
       *
       */
      std::vector<types::subdomain_id> true_subdomain_ids_of_cells;

      /**
       * 一个包含通过划分每个级别获得的单元的级别子域ID的向量。
       * 原始的分区信息被存储起来，以允许使用半人工单元的顺序DoF分布和分区函数。
       *
       */
      std::vector<std::vector<types::subdomain_id>>
        true_level_subdomain_ids_of_cells;
    };

    template <int dim, int spacedim>
    template <class Archive>
    void
    Triangulation<dim, spacedim>::load(Archive &ar, const unsigned int version)
    {
      dealii::Triangulation<dim, spacedim>::load(ar, version);
      partition();
      this->update_number_cache();
    }
  } // namespace shared
#else

  namespace shared
  {
    /**
     * 如果我们没有用MPI库实际配置deal.II，编译器为并行共享三角计算选择的假类。这个类的存在使得我们可以在整个库中引用
     * parallel::shared::Triangulation 对象，即使它被禁用。
     * 由于这个类的构造函数被删除，实际上不能创建这样的对象，因为考虑到MPI不可用，这将是毫无意义的。
     *
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public dealii::parallel::TriangulationBase<dim, spacedim>
    {
    public:
      /**
       * 构造器。已删除，以确保不能构造这种类型的对象（也可参见类的文档）。
       *
       */
      Triangulation() = delete;

      /**
       * 如果支持多级层次结构并已构建，则返回。
       *
       */
      virtual bool
      is_multilevel_hierarchy_constructed() const override;

      /**
       * 一个假函数，用于返回空向量。
       *
       */
      const std::vector<types::subdomain_id> &
      get_true_subdomain_ids_of_cells() const;

      /**
       * 一个假函数，用于返回空向量。
       *
       */
      const std::vector<types::subdomain_id> &
      get_true_level_subdomain_ids_of_cells(const unsigned int level) const;

      /**
       * 一个假函数，它总是返回真。
       *
       */
      bool
      with_artificial_cells() const;

    private:
      /**
       * 一个假的向量。
       *
       */
      std::vector<types::subdomain_id> true_subdomain_ids_of_cells;

      /**
       * 一个假的向量。
       *
       */
      std::vector<types::subdomain_id> true_level_subdomain_ids_of_cells;
    };
  } // namespace shared


#endif
} // namespace parallel


namespace internal
{
  namespace parallel
  {
    namespace shared
    {
      /**
       *
       */
      template <int dim, int spacedim = dim>
      class TemporarilyRestoreSubdomainIds : public Subscriptor
      {
      public:
        /**
         * 构造函数。                如果提供的Triangulation是
         * parallel::shared::Triangulation.
         * 类型，则存储所有活动单元的子域ID，用它们的真实子域ID等价物替换。
         *
         */
        TemporarilyRestoreSubdomainIds(
          const Triangulation<dim, spacedim> &tria);

        /**
         * 解构器。                将 parallel::shared::Triangulation
         * 上所有活动单元的子域ID返回到它们之前的状态。
         *
         */
        ~TemporarilyRestoreSubdomainIds();

      private:
        /**
         * 修改后的 parallel::shared::Triangulation. 。
         *
         */
        const SmartPointer<
          const dealii::parallel::shared::Triangulation<dim, spacedim>>
          shared_tria;

        /**
         * 一个向量，在 parallel::shared::Triangulation.
         * 上的所有活动单元被修改之前，临时存储它们的子域ID。
         *
         */
        std::vector<unsigned int> saved_subdomain_ids;
      };
    } // namespace shared
  }   // namespace parallel
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


