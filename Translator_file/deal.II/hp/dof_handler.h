//include/deal.II-translator/hp/dof_handler_0.txt
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

#ifndef dealii_hp_dof_handler_h
#define dealii_hp_dof_handler_h

#include <deal.II/dofs/dof_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   *
   *
   *
   *
   *
   *
   * - 网格细化时，子单元会继承父单元的活动FE索引。
   *
   * - 当粗化单元时，（现在活动的）父单元将被分配一个活动的FE索引，该索引由其（不再活动的）子单元决定，遵循FiniteElementDomination逻辑。在以前分配给前子女的元素集合中，我们选择一个由所有子女支配的元素作为父单元。如果没有找到，我们就在整个集合中挑选一个被所有以前的孩子支配的最主要的元素。关于这个主题的进一步信息，请参见 hp::FECollection::find_dominated_fe_extended() 。
   * @note  有限元素需要先通过调用set_fe()或distribution_dofs()来分配给每个单元格，以使这个功能可用。      <h3>Active FE indices and parallel meshes</h3> 当这个类与 parallel::shared::Triangulation 或 parallel::distributed::Triangulation, 一起使用时，你只能在本地拥有的单元上设置活动FE指数，使用诸如 <code>cell-@>set_active_fe_index(...)</code> 的调用。  另一方面，不允许在幽灵或人工单元上设置活动FE指数。    然而，幽灵单元确实获得了什么元素在其上处于活动状态的信息：每当你调用 hp::DoFHandler::distribute_dofs(), 时，所有参与并行网格的处理器都会以这样的方式交换信息，幽灵单元上的活动FE指数等于在拥有该特定幽灵单元的处理器上设置的活动FE指数。  因此，人们可以<i>query</i>幽灵单元上的 @p active_fe_index ，只是不能用手去设置。    在人工单元上，没有关于那里使用的 @p active_fe_index 的信息可用。这是因为我们甚至不知道这些细胞是否存在，即使存在，当前的处理器也不知道关于它们的任何具体信息。  更多信息见 @ref GlossArtificialCell "人工细胞的词汇表条目"
   * 。    在细化和粗化过程中，关于每个单元的 @p
   * active_fe_index 的信息将被自动转移。    然而，使用
   * parallel::distributed::Triangulation 和 hp::DoFHandler
   * 需要在序列化过程中额外注意，因为活动FE指数的信息不会被自动转移。这必须使用prepare_for_serialization_of_active_fe_indices()和deserialize_active_fe_indices()函数手动完成。前者必须在调用
   * parallel::distributed::Triangulation::save()
   * 之前调用，后者需要在
   * parallel::distributed::Triangulation::load(). 之后运行。
   * 如果进一步的数据将通过
   * parallel::distributed::CellDataTransfer,
   * parallel::distributed::SolutionTransfer, 或 Particles::ParticleHandler
   * 类附加到三角形上，所有相应的准备和反序列化函数调用需要以相同的顺序发生。更多信息请参考
   * parallel::distributed::SolutionTransfer 的文档。
   * @ingroup dofs   @deprecated  基本的 dealii::DoFHandler
   * 现在能够进行hp-adaptation。
   *
   */
  template <int dim, int spacedim = dim>
  using DoFHandler DEAL_II_DEPRECATED = dealii::DoFHandler<dim, spacedim>;
} // namespace hp

DEAL_II_NAMESPACE_CLOSE

#endif


