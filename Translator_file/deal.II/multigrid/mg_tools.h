//include/deal.II-translator/multigrid/mg_tools_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_mg_tools_h
#define dealii_mg_tools_h

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <set>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class DoFHandler;
class MGConstrainedDoFs;
#endif

 /* !@addtogroup mg */ 
 /* 
     * @{ */ 

/**
 * 这是在多级三角中操作和处理自由度数量的函数集合。它在目的和功能上与
 * @p DoFTools
 * 命名空间相似，但对DoFHandler对象的层次进行操作。更多信息见那里和成员函数的文档。
 *
 *
 */
namespace MGTools
{
  /**
   * 计算多层次方法的行长向量。
   *
   */
  template <int dim, int spacedim>
  void
  compute_row_length_vector(
    const DoFHandler<dim, spacedim> &dofs,
    const unsigned int               level,
    std::vector<unsigned int> &      row_lengths,
    const DoFTools::Coupling         flux_couplings = DoFTools::none);

  /**
   * 计算多级方法的行长向量，对块状耦合进行优化。
   *
   */
  template <int dim, int spacedim>
  void
  compute_row_length_vector(const DoFHandler<dim, spacedim> &   dofs,
                            const unsigned int                  level,
                            std::vector<unsigned int> &         row_lengths,
                            const Table<2, DoFTools::Coupling> &couplings,
                            const Table<2, DoFTools::Coupling> &flux_couplings);

  /**
   * 写出属于指定 @p
   * 级的矩阵的稀疏性结构。稀疏性模式没有被压缩，所以在创建实际的矩阵之前，你必须自己压缩矩阵，使用
   * <tt>SparsityPatternType::compress()</tt>.
   * 可选的AffineConstraints参数允许定义水平矩阵的约束，如Dirichlet边界条件。请注意，需要考虑典型层次矩阵上的悬空节点，因为只考虑一个层次。关于参数的更多细节见
   * DoFTools::make_sparsity_pattern() 。
   *
   */
  template <int dim,
            int spacedim,
            typename SparsityPatternType,
            typename number = double>
  void
  make_sparsity_pattern(
    const DoFHandler<dim, spacedim> &dof_handler,
    SparsityPatternType &            sparsity,
    const unsigned int               level,
    const AffineConstraints<number> &constraints = AffineConstraints<number>(),
    const bool                       keep_constrained_dofs = true);

  /**
   * 做一个包括通量的不连续Galerkin方法的稀疏模式。
   * @see   @ref make_sparsity_pattern  和  @ref DoFTools 。
   *
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_handler,
                             SparsityPatternType &            sparsity,
                             const unsigned int               level);

  /**
   * 为细化边缘的通量创建稀疏模式。矩阵将精细层次空间
   * @p level 的函数映射到更粗的空间。
   * make_flux_sparsity_pattern()
   *
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern_edge(const DoFHandler<dim, spacedim> &dof_handler,
                                  SparsityPatternType &            sparsity,
                                  const unsigned int               level);
  /**
   * 这个函数与另一个同名函数的作用相同，但它得到两个额外的系数矩阵。如果在系数矩阵中有一个非零的条目连接它们的相关组件，则只会为两个基函数生成一个矩阵条目。
   * 一个单元中的耦合有一个矩阵，通量中出现的耦合有一个矩阵。
   *
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &   dof,
                             SparsityPatternType &               sparsity,
                             const unsigned int                  level,
                             const Table<2, DoFTools::Coupling> &int_mask,
                             const Table<2, DoFTools::Coupling> &flux_mask);

  /**
   * 为细化边缘的通量创建稀疏模式。该矩阵将精细层次空间
   * @p level 的一个函数映射到较粗的空间。
   * 这是一个将模式限制在实际需要的元素上的版本。
   * make_flux_sparsity_pattern()
   *
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern_edge(
    const DoFHandler<dim, spacedim> &   dof_handler,
    SparsityPatternType &               sparsity,
    const unsigned int                  level,
    const Table<2, DoFTools::Coupling> &flux_mask);


  /**
   * 为多网格计算中使用的interface_in/out矩阵创建稀疏性模式。这些矩阵包含一个条目，代表细化边上的自由度与不在某一级细化边上的自由度的耦合。
   *
   */
  template <int dim, int spacedim, typename SparsityPatternType>
  void
  make_interface_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_handler,
                                  const MGConstrainedDoFs &mg_constrained_dofs,
                                  SparsityPatternType &    sparsity,
                                  const unsigned int       level);


  /**
   * 计算每个层次上的自由度块。
   * 结果是一个向量，其中包含每个级别的向量，包含每个块的自由度数量（访问是<tt>result[level][block]</tt>）。
   *
   */
  template <int dim, int spacedim>
  void
  count_dofs_per_block(
    const DoFHandler<dim, spacedim> &                  dof_handler,
    std::vector<std::vector<types::global_dof_index>> &dofs_per_block,
    std::vector<unsigned int>                          target_block = {});

  /**
   * 在每一级上逐一计算道夫。
   * 结果是一个向量，其中包含每个级别的向量，包含每个组件的道夫数（访问是<tt>result[level][component]</tt>）。
   *
   */
  template <int dim, int spacedim>
  void
  count_dofs_per_component(
    const DoFHandler<dim, spacedim> &                  mg_dof,
    std::vector<std::vector<types::global_dof_index>> &result,
    const bool                                         only_once        = false,
    std::vector<unsigned int>                          target_component = {});

  /**
   * 生成一个在域的边界上的自由度的列表，这些自由度应该从矩阵中消除，因为它们将受到Dirichlet边界条件的限制。
   * 这相当于 VectorTools::interpolate_boundary_values,
   * 的多级方法，但由于多级方法没有自己的右手边，作为function_map参数一部分的函数对象所返回的函数值被忽略了。
   * @arg
   * <tt>boundary_indices</tt>是一个向量，返回时包含每个层次的自由度的所有索引，这些自由度在function_map参数确定的边界部分。它的长度必须与自由度处理程序对象中的级别数相匹配。
   * @p boundary_indices
   * 中以前的内容不会被覆盖，而是被添加到。
   *
   */
  template <int dim, int spacedim>
  void
  make_boundary_list(
    const DoFHandler<dim, spacedim> &mg_dof,
    const std::map<types::boundary_id, const Function<spacedim> *>
      &                                             function_map,
    std::vector<std::set<types::global_dof_index>> &boundary_indices,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 与上面的函数相同，但在每个级别上返回一个IndexSet而不是一个
   * std::set<unsigned  int>。     @p boundary_indices
   * 中以前的内容不会被覆盖，而是被添加到。
   *
   */
  template <int dim, int spacedim>
  void
  make_boundary_list(const DoFHandler<dim, spacedim> &           mg_dof,
                     const std::map<types::boundary_id,
                                    const Function<spacedim> *> &function_map,
                     std::vector<IndexSet> &boundary_indices,
                     const ComponentMask &  component_mask = ComponentMask());

  /**
   * 与上面的函数相同，但在每一层上返回一个IndexSet而不是一个
   * std::set<unsigned  int>，并使用一个 std::set
   * 的boundary_ids作为输入。    以前在 @p boundary_indices
   * 中的内容不会被覆盖，而是添加到。
   *
   */
  template <int dim, int spacedim>
  void
  make_boundary_list(const DoFHandler<dim, spacedim> &   mg_dof,
                     const std::set<types::boundary_id> &boundary_ids,
                     std::vector<IndexSet> &             boundary_indices,
                     const ComponentMask &component_mask = ComponentMask());

  /**
   * 对于多网格层次结构中的每一层，产生一个IndexSet，表明哪些自由度是沿着这一层的界面到只存在于更粗层的单元。
   *
   */
  template <int dim, int spacedim>
  void
  extract_inner_interface_dofs(const DoFHandler<dim, spacedim> &mg_dof_handler,
                               std::vector<IndexSet> &          interface_dofs);

  /**
   * 返回可以作为多重网格计算中最粗层次的最高层次，即层次结构中网格覆盖整个领域的最高层次。这相当于活动网格上一个单元的最小级别。由于每个处理器只有一个网格的局部视图，每个处理器都必须调用这个函数。请注意，这是一个覆盖整个网格的全局最小值，因此每个处理器将返回相同的值。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  max_level_for_coarse_mesh(const Triangulation<dim, spacedim> &tria);

  /**
   * 返回多网格层次的平行分布的不平衡度。理想情况下，这个值等于1（每个处理器在每个层次上拥有相同数量的单元，对于大多数全局细化的网格来说，大约是真的）。大于1的值估计在一个几何多网格v-cycle中，与在一个完全分布的网格层次上的相同计算相比，应该看到的速度下降。
   * 这个函数是Triangulation所有等级之间的MPI集体调用，因此需要从所有等级调用。
   * @note  这个函数要求
   * parallel::TriangulationBase::is_multilevel_hierarchy_constructed()
   * 为真，这可以通过在构造三角结构时设置construct_multigrid_hierarchy标志来控制。
   *
   */
  template <int dim, int spacedim>
  double
  workload_imbalance(const Triangulation<dim, spacedim> &tria);

} // namespace MGTools

 /* @} */ 

DEAL_II_NAMESPACE_CLOSE

#endif


