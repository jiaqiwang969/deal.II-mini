//include/deal.II-translator/distributed/repartitioning_policy_tools_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_distributed_repartitioning_policy_tools_h
#define dealii_distributed_repartitioning_policy_tools_h

#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个具有重新分区策略的命名空间。这些类会返回三角形对象中活跃的本地拥有单元和幽灵单元的新主人的向量。返回的向量可用于，例如，在
 * TriangulationDescription::Utilities::create_description_from_triangulation()
 * 中，基于给定的Triangulation和预先描述的分区创建一个
 * TriangulationDescription::Description ，这可用于建立一个
 * parallel::fullydistributed::Triangulation 对象。
 * 这些策略也可以在
 * MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence()
 * 的背景下使用，以规定全局粗化多网格划分的多网格层次中的任意分区。
 *
 *
 */
namespace RepartitioningPolicyTools
{
  /**
   * 一个重新分区策略的基类。
   *
   */
  template <int dim, int spacedim = dim>
  class Base
  {
  public:
    /**
     * 返回一个活跃的本地拥有单元和幽灵单元的新主人的向量。
     *
     */
    virtual LinearAlgebra::distributed::Vector<double>
    partition(const Triangulation<dim, spacedim> &tria_coarse_in) const = 0;
  };

  /**
   * 一个简单地返回空向量的假策略，在
   * MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence()
   * 中被解释为三角结构没有被重新分区的方式。
   *
   */
  template <int dim, int spacedim = dim>
  class DefaultPolicy : public Base<dim, spacedim>
  {
  public:
    virtual LinearAlgebra::distributed::Vector<double>
    partition(
      const Triangulation<dim, spacedim> &tria_coarse_in) const override;
  };

  /**
   * 一个策略，根据第一子策略，基于基础三角形划分粗大的网格。被分割的三角形应该能够通过一连串的（全局）粗化步骤获得。
   *
   */
  template <int dim, int spacedim = dim>
  class FirstChildPolicy : public Base<dim, spacedim>
  {
  public:
    /**
     * 构造函数获取基础（精细）三角结构。
     *
     */
    FirstChildPolicy(const Triangulation<dim, spacedim> &tria_fine);

    virtual LinearAlgebra::distributed::Vector<double>
    partition(
      const Triangulation<dim, spacedim> &tria_coarse_in) const override;

  private:
    /**
     * 粗略单元的数量。
     *
     */
    const unsigned int n_coarse_cells;

    /**
     * 全局水平的数量。
     *
     */
    const unsigned int n_global_levels;

    /**
     * 从传递给构造函数的三角结构中构建的索引集。
     * 它包含了所有的单元，如果层级按照第一子策略进行划分的话，这些单元将被当前进程所拥有。
     *
     */
    IndexSet is_level_partitions;
  };

  /**
   * 一个允许指定每个进程的最小单元数的策略。如果达到一个阈值，进程就可能没有单元。
   *
   */
  template <int dim, int spacedim = dim>
  class MinimalGranularityPolicy : public Base<dim, spacedim>
  {
  public:
    /**
     * 构造器获取每个进程的最小单元数。
     *
     */
    MinimalGranularityPolicy(const unsigned int n_min_cells);

    virtual LinearAlgebra::distributed::Vector<double>
    partition(const Triangulation<dim, spacedim> &tria_in) const override;

  private:
    /**
     * 每个进程的最小单元数。
     *
     */
    const unsigned int n_min_cells;
  };

} // namespace RepartitioningPolicyTools

DEAL_II_NAMESPACE_CLOSE

#endif


