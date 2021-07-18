//include/deal.II-translator/distributed/grid_refinement_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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

#ifndef dealii_distributed_grid_refinement_h
#define dealii_distributed_grid_refinement_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/vector_tools_common.h>

#include <limits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace parallel
  {
    namespace distributed
    {
      namespace GridRefinement
      {
        /**
         * 计算标准向量的全局最大值和最小值。这些数据只在等级为0的处理器上返回，所有其他处理器都得到一对零。
         *
         */
        template <typename number>
        std::pair<number, number>
        compute_global_min_and_max_at_root(
          const dealii::Vector<number> &criteria,
          const MPI_Comm &              mpi_communicator);

        namespace RefineAndCoarsenFixedNumber
        {
          /**
           * 计算一个阈值，以便正好n_target_cells有一个较大的值。
           *
           */
          template <typename number>
          number
          compute_threshold(const dealii::Vector<number> &   criteria,
                            const std::pair<double, double> &global_min_and_max,
                            const types::global_cell_index   n_target_cells,
                            const MPI_Comm &                 mpi_communicator);
        } // namespace RefineAndCoarsenFixedNumber

        namespace RefineAndCoarsenFixedFraction
        {
          /**
           * 计算一个阈值，使所有标准[i]上累积的误差，使标准[i]>阈值大于target_error。
           *
           */
          template <typename number>
          number
          compute_threshold(const dealii::Vector<number> &   criteria,
                            const std::pair<double, double> &global_min_and_max,
                            const double                     target_error,
                            const MPI_Comm &                 mpi_communicator);
        } // namespace RefineAndCoarsenFixedFraction
      }   // namespace GridRefinement
    }     // namespace distributed
  }       // namespace parallel
} // namespace internal



namespace parallel
{
  namespace distributed
  {
    /**
     * 这个命名空间提供了一个函数集合，帮助细化和粗化三角形。尽管命名空间的名字，这些函数实际上并不<i>refine</i>三角化，而只是<i>mark
     * cells for refinement or
     * coarsening</i>。换句话说，它们执行自适应有限元循环中典型的
     * "求解-估计-标记-细化 "循环中的 "标记 "部分。
     * 与命名空间 dealii::GridRefinement,
     * 中的函数相反，当前命名空间中的函数旨在用于分布式网格，即
     * parallel::distributed::Triangulation. 类型的对象。
     * @ingroup grid
     *
     */
    namespace GridRefinement
    {
      /**
       * 就像 dealii::GridRefinement::refine_and_coarsen_fixed_number,
       * 一样，但用于并行分布式三角计算。
       * 标准向量需要是当前三角形上所有活动单元的细化标准向量，即需要长度为
       * <code>tria.n_active_cells()</code>  （而不是
       * <code>tria.n_locally_owned_active_cells()</code>
       * ）。换句话说，这个向量需要包括鬼魂和人工单元的条目。然而，目前的函数将只查看与那些实际属于本地的单元对应的指标，而忽略所有其他单元的指标。然后，该函数将在所有存储部分三角形的处理器之间进行协调，以便在最后对所有
       * Triangulation::n_global_active_cells() 活动单元的一部分 @p
       * top_fraction_of_cells
       * 进行精炼，而不是每个处理器上的一部分
       * Triangulation::n_locally_active_cells 。
       * 换句话说，可能在某些处理器上，根本就没有单元被精炼。
       * 对于被粗化的单元的部分也是如此。
       * @param[in,out]  tria 三角形
       * 这个函数应该对其单元进行粗化和细化的标记。
       * @param[in]  criteria
       * 当前三角形上的每个网格单元的细化标准。条目不可以是负数。
       * @param[in]  top_fraction_of_cells 要精简的单元的比例。
       * 如果这个数字是零，没有单元会被精简。如果它等于1，结果将被标记为全局细化。
       * @param[in]  bottom_fraction_of_cells
       * 要被粗化的单元格的比例。如果这个数字为0，则没有单元会被粗化。
       * @param[in]  max_n_cells
       * 这个参数可以用来指定一个最大的细胞数。如果细化时超过这个数字，那么细化和粗化的比例将被调整，以达到最大的细胞数。但是要注意，由于
       * Triangulation::MeshSmoothing,
       * 的细化扩散，这个数字只是一个指标。这个参数的默认值是对单元格的数量没有限制。
       *
       */
      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_number(
        parallel::distributed::Triangulation<dim, spacedim> &tria,
        const dealii::Vector<Number> &                       criteria,
        const double                   top_fraction_of_cells,
        const double                   bottom_fraction_of_cells,
        const types::global_cell_index max_n_cells =
          std::numeric_limits<types::global_cell_index>::max());

      /**
       * 与 dealii::GridRefinement::refine_and_coarsen_fixed_fraction,
       * 类似，但用于平行分布式三角计算。
       * 标准向量需要是当前三角形上所有活动单元的细化标准向量，即需要长度为
       * <code>tria.n_active_cells()</code>  （而不是
       * <code>tria.n_locally_owned_active_cells()</code>
       * ）。换句话说，这个向量需要包括鬼魂和人工单元的条目。然而，目前的函数将只查看与那些实际属于本地的单元对应的指标，而忽略所有其他单元的指标。然后，该函数将在所有存储部分三角图的处理器之间进行协调，以便在最后提炼出
       * Triangulation::n_global_active_cells
       * （而不是每个处理器上的
       * Triangulation::n_locally_owned_active_cells()
       * ）的最小部分，这些最小部分共同构成了 @p
       * top_fraction_of_error
       * 的总误差。换句话说，可能在某些处理器上，根本就没有单元被精炼。
       * 对于被粗化的单元的部分也是如此。
       * @param[in,out]  tria 三角形
       * 这个函数应该对其单元进行粗化和细化的标记。
       * @param[in]  criteria
       * 对当前三角形上的每个网格单元计算的细化准则。条目不可以是负数。
       * @param[in]  top_fraction_of_error
       * 应该被细化的总估算值的分数。如果这个数字为零，则没有单元被精炼。如果它等于1，结果将被标记为全局细化。
       * @param[in]  bottom_fraction_of_error
       * 粗化的估计值的分数。如果这个数字为零，则没有单元会被粗化。
       * @param[in]  norm_type
       * 为了确定阈值，单元格子集上的综合误差被计算为这些单元格上的准则的规范。不同类型的准则可用于此目的，目前支持
       * VectorTools::NormType::L1_norm 和 VectorTools::NormType::L2_norm 。
       *
       */
      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_fraction(
        parallel::distributed::Triangulation<dim, spacedim> &tria,
        const dealii::Vector<Number> &                       criteria,
        const double                top_fraction_of_error,
        const double                bottom_fraction_of_error,
        const VectorTools::NormType norm_type = VectorTools::NormType::L1_norm);
    } // namespace GridRefinement
  }   // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_distributed_grid_refinement_h


