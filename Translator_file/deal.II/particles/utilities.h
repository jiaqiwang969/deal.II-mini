//include/deal.II-translator/particles/utilities_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_particles_utilities
#define dealii_particles_utilities

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/particles/particle_handler.h>


DEAL_II_NAMESPACE_OPEN



namespace Particles
{
  /**
   * 一个命名空间，用于提供处理ParticleHandler对象及其与DoFHandler对象耦合的工具的函数。
   *
   */
  namespace Utilities
  {
    /**
     * 为粒子创建一个插值稀疏模式。        给出一个代表域
     * $\Omega$ 的三角形，一个 $\Omega$
     * 中的粒子处理程序，以及一个标量有限元空间 $V(\Omega)
     * = \text{span}\{v_j\}_{j=0}^n$ ，计算组装矩阵\f[ M_{i,j}
     * \dealcoloneq v_j(x_i) , \f]所需的稀疏模式，其中 $V(\Omega)$
     * 是与`空间_dh`相关的有限元空间，索引`i`是由位置为`x_i`的粒子id给出的。
     * 在矢量值有限元空间的情况下，必须进行插值的分量可以用分量掩码来选择。只支持原始的有限元空间。
     * 当选择一个以上的分量时，产生的稀疏度将等于`particle_handler.n_global_particles()
     * mask.n_selected_components()`乘以`space_dh.n_dofs()`，而相应的矩阵条目由\f[
     * M_{(i*n_comps+k),j} \dealcoloneq v_j(x_i) \cdot e_{comp_j}, \f] ]
     * 其中`comp_j`是矢量值基函数`v_j`的唯一非零分量（等于`fe.system_to_component_index(j).first`），`k`对应于其在掩码选定分量中的索引，而
     * $e_{comp_j}$ 是`comp_j`方向上的单位矢量。        稀疏度
     * "是通过定位粒子处理程序中索引为`i'的粒子相对于嵌入三角形的位置来填充的
     * $\Omega$
     * ，并按照在掩码中选择的顺序，将其与组件掩码 @p
     * space_comps, 中指定的所有局部自由度耦合 @p space_comps.
     * ，如果一个粒子不在 $\Omega$
     * 内，它将被忽略，稀疏度的相应行将为空。
     * 可以用 @p constraints
     * 参数提供AffineConstraints类所支持的形式的约束。方法
     * AffineConstraints::add_entries_local_to_global()
     * 被用来填充最终的稀疏度模式。
     *
     */
    template <int dim,
              int spacedim,
              typename SparsityType,
              typename number = double>
    void
    create_interpolation_sparsity_pattern(
      const DoFHandler<dim, spacedim> &                space_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      SparsityType &                                   sparsity,
      const AffineConstraints<number> &                constraints =
        AffineConstraints<number>(),
      const ComponentMask &space_comps = ComponentMask());

    /**
     * 为粒子创建一个插值矩阵。
     * 给出一个代表域的三角形 $\Omega$ ，一个 $\Omega$
     * 中的粒子处理程序，和一个标量有限元空间 $V(\Omega) =
     * \text{span}\{v_j\}_{j=0}^n$ ，计算矩阵\f[ M_{ij} \dealcoloneq
     * v_j(x_i) , \f]，其中 $V(\Omega)$
     * 是与`空间_dh`相关的有限元空间，索引`i`是由位置为`x_i`的粒子id给出。
     * 在矢量值有限元空间的情况下，必须进行插值的分量可以用分量掩码来选择。只支持原始的有限元空间。
     * 当选择一个以上的分量时，产生的稀疏度将等于`particle_handler.n_global_particles()
     * mask.n_selected_components()`乘以`space_dh.n_dofs()`，而相应的矩阵条目由\f[
     * M_{(i*n_comps+k),j} \dealcoloneq v_j(x_i) \cdot e_{comp_j}, \f] ]
     * 其中`comp_j`是矢量值基函数`v_j`的唯一非零分量（等于`fe.system_to_component_index(j).first`），`k`对应于它在掩码选定分量中的索引，
     * $e_{comp_j}$ 是`comp_j`方向上的单位矢量。
     * 矩阵通过定位粒子处理程序中索引为`i`的粒子相对于嵌入三角的位置来填充
     * $\Omega$ ，并按照在掩码中选择的顺序将其与组件掩码
     * @p space_comps, 中指定的所有局部自由度耦合 @p space_comps.
     * 。 如果一个粒子不在 $\Omega$
     * 内，它将被忽略，矩阵的相应行将为零。        可以用
     * @p constraints
     * 参数提供AffineConstraints类所支持的形式的约束。方法
     * AffineConstraints::distribute_local_to_global()
     * 用于分配矩阵的条目以尊重给定的约束。
     *
     */
    template <int dim, int spacedim, typename MatrixType>
    void
    create_interpolation_matrix(
      const DoFHandler<dim, spacedim> &                space_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      MatrixType &                                     matrix,
      const AffineConstraints<typename MatrixType::value_type> &constraints =
        AffineConstraints<typename MatrixType::value_type>(),
      const ComponentMask &space_comps = ComponentMask());

    /**
     * 给定一个DoFHandler和一个ParticleHandler，在粒子的位置插值一个矢量场。结果存储在一个输出向量中，其大小与本地拥有的粒子数和活动部件数相对应。
     * @param[in]  particle_handler
     * 粒子处理程序，其粒子作为插值点。          @param[in]
     * field_vector
     * 要插值的场的矢量。这个向量必须与提供的dof_handler相一致。
     * @param[in,out]  interpolated_field
     * 在粒子位置的场的内插值。矢量的大小必须是n_locally_owned_particles乘以n_components
     * @param[in]  field_comps
     * 一个可选的组件掩码，决定哪一个矢量场的子集被插值。
     *
     */
    template <int dim,
              int spacedim,
              typename InputVectorType,
              typename OutputVectorType>
    void
    interpolate_field_on_particles(
      const DoFHandler<dim, spacedim> &                field_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      const InputVectorType &                          field_vector,
      OutputVectorType &                               interpolated_field,
      const ComponentMask &field_comps = ComponentMask())
    {
      if (particle_handler.n_locally_owned_particles() == 0)
        {
          interpolated_field.compress(VectorOperation::add);
          return; // nothing else to do here
        }

      const auto &tria     = field_dh.get_triangulation();
      const auto &fe       = field_dh.get_fe();
      auto        particle = particle_handler.begin();

      // Take care of components
      const ComponentMask comps =
        (field_comps.size() == 0 ? ComponentMask(fe.n_components(), true) :
                                   field_comps);
      AssertDimension(comps.size(), fe.n_components());
      const auto n_comps = comps.n_selected_components();

      AssertDimension(field_vector.size(), field_dh.n_dofs());
      AssertDimension(interpolated_field.size(),
                      particle_handler.get_next_free_particle_index() *
                        n_comps);

      // Global to local indices
      std::vector<unsigned int> space_gtl(fe.n_components(),
                                          numbers::invalid_unsigned_int);
      for (unsigned int i = 0, j = 0; i < space_gtl.size(); ++i)
        if (comps[i])
          space_gtl[i] = j++;

      std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());

      while (particle != particle_handler.end())
        {
          const auto &cell = particle->get_surrounding_cell(tria);
          const auto &dh_cell =
            typename DoFHandler<dim, spacedim>::cell_iterator(*cell, &field_dh);
          dh_cell->get_dof_indices(dof_indices);
          const auto pic = particle_handler.particles_in_cell(cell);

          Assert(pic.begin() == particle, ExcInternalError());
          for (unsigned int i = 0; particle != pic.end(); ++particle, ++i)
            {
              const Point<dim> reference_location =
                particle->get_reference_location();

              const auto id = particle->get_id();

              for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
                {
                  const auto comp_j =
                    space_gtl[fe.system_to_component_index(j).first];
                  if (comp_j != numbers::invalid_unsigned_int)
                    interpolated_field[id * n_comps + comp_j] +=
                      fe.shape_value(j, reference_location) *
                      field_vector(dof_indices[j]);
                }
            }
        }
      interpolated_field.compress(VectorOperation::add);
    }

  } // namespace Utilities
} // namespace Particles
DEAL_II_NAMESPACE_CLOSE

#endif


