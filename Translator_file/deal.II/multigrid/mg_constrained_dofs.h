//include/deal.II-translator/multigrid/mg_constrained_dofs_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_mg_constrained_dofs_h
#define dealii_mg_constrained_dofs_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/multigrid/mg_tools.h>

#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
class DoFHandler;
#endif


/**
 * 水平向量的边界约束和细化边缘约束的集合。
 *
 *
 * @ingroup mg
 *
 *
 */
class MGConstrainedDoFs : public Subscriptor
{
public:
  using size_dof = std::vector<std::set<types::global_dof_index>>::size_type;
  /**
   * 用从dof
   * handler对象中提取的悬挂节点约束来填充内部数据结构。只在自然边界条件下工作。
   * 也存在一个设置边界约束的姐妹函数。
   * 这个函数确保在每个层次上，细化层次的内部边缘的自由度被正确处理，但假设不存在迪里希特边界条件，则不触及域的边界自由度。
   * 此外，这个调用在每个层次上设置了一个AffineConstraints对象，该对象包含可能的周期性约束，以防这些约束被添加到底层三角形中。该AffineConstraints对象可以通过get_level_constraints(level)查询。请注意，目前这个类中周期性约束的实现不支持周期性定义中的旋转矩阵，也就是说，
   * GridTools::collect_periodic_faces()
   * 中的相应参数可能与身份矩阵不同。
   * 如果没有传递level_relevant_dofs作为第二个参数，该函数将使用由
   * DoFTools::extract_locally_relevant_level_dofs().
   * 提取的本地相关级别DoFs。否则，用户提供的IndexSets（应定义本地相关DoFs的超集）将在每个级别上使用，以允许用户向受限DoFs集添加额外的指数。
   *
   */
  template <int dim, int spacedim>
  void
  initialize(const DoFHandler<dim, spacedim> &dof,
             const MGLevelObject<IndexSet> &  level_relevant_dofs =
               MGLevelObject<IndexSet>());

  /**
   * 用有关Dirichlet边界DoF的信息填充内部数据结构。
   * 在设置悬挂节点约束之前，必须调用initialize()函数。
   * 这个函数可以被多次调用，以允许考虑不同组件的不同边界_ID集。
   *
   */
  template <int dim, int spacedim>
  void
  make_zero_boundary_constraints(
    const DoFHandler<dim, spacedim> &   dof,
    const std::set<types::boundary_id> &boundary_ids,
    const ComponentMask &               component_mask = ComponentMask());

  /**
   * 添加用户定义的约束条件，用于水平 @p level.
   * 用户可以多次调用这个函数，任何新的、冲突的约束条件将覆盖该DoF的先前约束。
   * 在转移之前，用户定义的约束将被分配到源矢量，然后任何使用make_zero_boundary_constraints()设置的DoF索引将被覆盖为零值。
   * @note 目前这只在MGTransferMatrixFree中实现。
   *
   */
  void
  add_user_constraints(const unsigned int               level,
                       const AffineConstraints<double> &constraints_on_level);

  /**
   * 用无法线通量边界道夫的信息填充内部数据结构。
   * 这个函数只限于那些无法线通量边界的面与x轴、y轴或z轴为法线的网格。同时，对于一个特定的边界ID，所有的面必须朝向同一个方向，也就是说，对X轴法线的边界必须与对Y轴或Z轴法线的边界具有不同的边界ID，以此类推。如果网格是用
   * <tt>GridGenerator::hyper_cube()</tt>
   * 函数生成的，在网格生成时设置<tt>colorize=true</tt>，并为每个无法线通量边界调用make_no_normal_flux_constraints()就可以了。
   *
   */
  template <int dim, int spacedim>
  void
  make_no_normal_flux_constraints(const DoFHandler<dim, spacedim> &dof,
                                  const types::boundary_id         bid,
                                  const unsigned int first_vector_component);

  /**
   * 重置数据结构。
   *
   */
  void
  clear();

  /**
   * 判断一个道夫指数是否受到边界约束。
   *
   */
  bool
  is_boundary_index(const unsigned int            level,
                    const types::global_dof_index index) const;

  /**
   * 判断一个道夫指数是否处于细化边缘。
   *
   */
  bool
  at_refinement_edge(const unsigned int            level,
                     const types::global_dof_index index) const;


  /**
   * 确定是否应该设置给定层次上的界面矩阵的（i,j）项。这是就dof
   * i而言的，也就是说，如果i在细化边缘，j不在细化边缘，并且两者都不在外部边界上，则返回true。
   *
   */
  bool
  is_interface_matrix_entry(const unsigned int            level,
                            const types::global_dof_index i,
                            const types::global_dof_index j) const;

  /**
   * 返回给定层次上受Dirichlet边界条件约束的层次道夫的索引（如initialize()中的
   * @p function_map 参数所设定）。
   * 指数被限制在本地相关的层面道夫的集合中。
   *
   */
  const IndexSet &
  get_boundary_indices(const unsigned int level) const;


  /**
   * 返回位于细化边上的给定层次上的道夫的指数（道夫在面与邻居的关系更粗）。
   *
   */
  const IndexSet &
  get_refinement_edge_indices(unsigned int level) const;


  /**
   * 返回是否在initialize()中设置了Dirichlet边界指数。
   *
   */
  bool
  have_boundary_indices() const;

  /**
   * 返回给定层次的AffineConstraints对象，包含周期性约束（如果在三角剖面上启用）。
   *
   */
  const AffineConstraints<double> &
  get_level_constraints(const unsigned int level) const;

  /**
   * 返回给定水平的用户定义的约束矩阵。这些约束是用函数add_user_constraints()设置的，不应该包含在make_zero_boundary_constraints()中设置的DoF指数的约束，因为它们将在传输过程中被覆盖。
   *
   */
  const AffineConstraints<double> &
  get_user_constraint_matrix(const unsigned int level) const;

private:
  /**
   * 每一层的边界DoF的指数。
   *
   */
  std::vector<IndexSet> boundary_indices;

  /**
   * 某一层次上的自由度，它们生活在该层次和更粗层次上的单元之间的细化边缘上。
   *
   */
  std::vector<IndexSet> refinement_edge_indices;

  /**
   * 包含每个层的潜在周期性边界条件信息的约束矩阵。
   *
   */
  std::vector<AffineConstraints<double>> level_constraints;

  /**
   * 由用户定义的约束矩阵。
   *
   */
  std::vector<AffineConstraints<double>> user_constraints;
};


template <int dim, int spacedim>
inline void
MGConstrainedDoFs::initialize(
  const DoFHandler<dim, spacedim> &dof,
  const MGLevelObject<IndexSet> &  level_relevant_dofs)
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
  level_constraints.clear();
  user_constraints.clear();

  const unsigned int nlevels   = dof.get_triangulation().n_global_levels();
  const unsigned int min_level = level_relevant_dofs.min_level();
  const unsigned int max_level = (level_relevant_dofs.max_level() == 0) ?
                                   nlevels - 1 :
                                   level_relevant_dofs.max_level();
  const bool user_level_dofs =
    (level_relevant_dofs.max_level() == 0) ? false : true;

  // At this point level_constraint and refinement_edge_indices are empty.
  refinement_edge_indices.resize(nlevels);
  level_constraints.resize(nlevels);
  user_constraints.resize(nlevels);
  for (unsigned int l = min_level; l <= max_level; ++l)
    {
      if (user_level_dofs)
        {
          level_constraints[l].reinit(level_relevant_dofs[l]);
        }
      else
        {
          IndexSet relevant_dofs;
          DoFTools::extract_locally_relevant_level_dofs(dof, l, relevant_dofs);
          level_constraints[l].reinit(relevant_dofs);
        }

      // Loop through relevant cells and faces finding those which are periodic
      // neighbors.
      typename DoFHandler<dim, spacedim>::cell_iterator cell = dof.begin(l),
                                                        endc = dof.end(l);
      for (; cell != endc; ++cell)
        if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
          {
            for (auto f : cell->face_indices())
              if (cell->has_periodic_neighbor(f) &&
                  cell->periodic_neighbor(f)->level() == cell->level())
                {
                  if (cell->is_locally_owned_on_level())
                    {
                      Assert(
                        cell->periodic_neighbor(f)->level_subdomain_id() !=
                          numbers::artificial_subdomain_id,
                        ExcMessage(
                          "Periodic neighbor of a locally owned cell must either be owned or ghost."));
                    }
                  // Cell is a level-ghost and its neighbor is a
                  // level-artificial cell nothing to do here
                  else if (cell->periodic_neighbor(f)->level_subdomain_id() ==
                           numbers::artificial_subdomain_id)
                    {
                      Assert(cell->is_locally_owned_on_level() == false,
                             ExcInternalError());
                      continue;
                    }

                  const unsigned int dofs_per_face =
                    dof.get_fe(0).n_dofs_per_face(f);
                  std::vector<types::global_dof_index> dofs_1(dofs_per_face);
                  std::vector<types::global_dof_index> dofs_2(dofs_per_face);

                  cell->periodic_neighbor(f)
                    ->face(cell->periodic_neighbor_face_no(f))
                    ->get_mg_dof_indices(l, dofs_1, 0);
                  cell->face(f)->get_mg_dof_indices(l, dofs_2, 0);
                  // Store periodicity information in the level
                  // AffineConstraints object. Skip DoFs for which we've
                  // previously entered periodicity constraints already; this
                  // can happen, for example, for a vertex dof at a periodic
                  // boundary that we visit from more than one cell
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    if (level_constraints[l].can_store_line(dofs_2[i]) &&
                        level_constraints[l].can_store_line(dofs_1[i]) &&
                        !level_constraints[l].is_constrained(dofs_2[i]) &&
                        !level_constraints[l].is_constrained(dofs_1[i]))
                      {
                        level_constraints[l].add_line(dofs_2[i]);
                        level_constraints[l].add_entry(dofs_2[i],
                                                       dofs_1[i],
                                                       1.);
                      }
                }
          }
      level_constraints[l].close();

      // Initialize with empty IndexSet of correct size
      refinement_edge_indices[l] = IndexSet(dof.n_dofs(l));
    }

  MGTools::extract_inner_interface_dofs(dof, refinement_edge_indices);
}


template <int dim, int spacedim>
inline void
MGConstrainedDoFs::make_zero_boundary_constraints(
  const DoFHandler<dim, spacedim> &   dof,
  const std::set<types::boundary_id> &boundary_ids,
  const ComponentMask &               component_mask)
{
  // allocate an IndexSet for each global level. Contents will be
  // overwritten inside make_boundary_list.
  const unsigned int n_levels = dof.get_triangulation().n_global_levels();
  Assert(boundary_indices.size() == 0 || boundary_indices.size() == n_levels,
         ExcInternalError());
  boundary_indices.resize(n_levels);

  MGTools::make_boundary_list(dof,
                              boundary_ids,
                              boundary_indices,
                              component_mask);
}



template <int dim, int spacedim>
inline void
MGConstrainedDoFs::make_no_normal_flux_constraints(
  const DoFHandler<dim, spacedim> &dof,
  const types::boundary_id         bid,
  const unsigned int               first_vector_component)
{
  // For a given boundary id, find which vector component is on the boundary
  // and set a zero boundary constraint for those degrees of freedom.
  const unsigned int n_components = dof.get_fe_collection().n_components();
  AssertIndexRange(first_vector_component + dim - 1, n_components);

  ComponentMask comp_mask(n_components, false);


  typename Triangulation<dim>::face_iterator
    face = dof.get_triangulation().begin_face(),
    endf = dof.get_triangulation().end_face();
  for (; face != endf; ++face)
    if (face->at_boundary() && face->boundary_id() == bid)
      for (unsigned int d = 0; d < dim; ++d)
        {
          Tensor<1, dim, double> unit_vec;
          unit_vec[d] = 1.0;

          const Tensor<1, dim> normal_vec =
            face->get_manifold().normal_vector(face, face->center());

          if (std::abs(std::abs(unit_vec * normal_vec) - 1.0) < 1e-10)
            comp_mask.set(d + first_vector_component, true);
          else
            Assert(
              std::abs(unit_vec * normal_vec) < 1e-10,
              ExcMessage(
                "We can currently only support no normal flux conditions "
                "for a specific boundary id if all faces are normal to the "
                "x, y, or z axis."));
        }

  Assert(comp_mask.n_selected_components() == 1,
         ExcMessage(
           "We can currently only support no normal flux conditions "
           "for a specific boundary id if all faces are facing in the "
           "same direction, i.e., a boundary normal to the x-axis must "
           "have a different boundary id than a boundary normal to the "
           "y- or z-axis and so on. If the mesh here was produced using "
           "GridGenerator::..., setting colorize=true during mesh generation "
           "and calling make_no_normal_flux_constraints() for each no normal "
           "flux boundary will fulfill the condition."));

  this->make_zero_boundary_constraints(dof, {bid}, comp_mask);
}


inline void
MGConstrainedDoFs::add_user_constraints(
  const unsigned int               level,
  const AffineConstraints<double> &constraints_on_level)
{
  AssertIndexRange(level, user_constraints.size());

  // Get the relevant DoFs from level_constraints if
  // the user constraint matrix has not been initialized
  if (user_constraints[level].get_local_lines().size() == 0)
    user_constraints[level].reinit(level_constraints[level].get_local_lines());

  user_constraints[level].merge(
    constraints_on_level,
    AffineConstraints<double>::MergeConflictBehavior::right_object_wins);
  user_constraints[level].close();
}


inline void
MGConstrainedDoFs::clear()
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
  user_constraints.clear();
}


inline bool
MGConstrainedDoFs::is_boundary_index(const unsigned int            level,
                                     const types::global_dof_index index) const
{
  if (boundary_indices.size() == 0)
    return false;

  AssertIndexRange(level, boundary_indices.size());
  return boundary_indices[level].is_element(index);
}

inline bool
MGConstrainedDoFs::at_refinement_edge(const unsigned int            level,
                                      const types::global_dof_index index) const
{
  AssertIndexRange(level, refinement_edge_indices.size());

  return refinement_edge_indices[level].is_element(index);
}

inline bool
MGConstrainedDoFs::is_interface_matrix_entry(
  const unsigned int            level,
  const types::global_dof_index i,
  const types::global_dof_index j) const
{
  const IndexSet &interface_dofs_on_level =
    this->get_refinement_edge_indices(level);

  return interface_dofs_on_level.is_element(i)     // at_refinement_edge(i)
         && !interface_dofs_on_level.is_element(j) // !at_refinement_edge(j)
         && !this->is_boundary_index(level, i)     // !on_boundary(i)
         && !this->is_boundary_index(level, j);    // !on_boundary(j)
}



inline const IndexSet &
MGConstrainedDoFs::get_boundary_indices(const unsigned int level) const
{
  AssertIndexRange(level, boundary_indices.size());
  return boundary_indices[level];
}



inline const IndexSet &
MGConstrainedDoFs::get_refinement_edge_indices(unsigned int level) const
{
  AssertIndexRange(level, refinement_edge_indices.size());
  return refinement_edge_indices[level];
}



inline bool
MGConstrainedDoFs::have_boundary_indices() const
{
  return boundary_indices.size() != 0;
}



inline const AffineConstraints<double> &
MGConstrainedDoFs::get_level_constraints(const unsigned int level) const
{
  AssertIndexRange(level, level_constraints.size());
  return level_constraints[level];
}



inline const AffineConstraints<double> &
MGConstrainedDoFs::get_user_constraint_matrix(const unsigned int level) const
{
  AssertIndexRange(level, user_constraints.size());
  return user_constraints[level];
}



DEAL_II_NAMESPACE_CLOSE

#endif


