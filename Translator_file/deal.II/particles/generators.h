//include/deal.II-translator/particles/generators_0.txt
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

#ifndef dealii_particles_particle_generator_h
#define dealii_particles_particle_generator_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

#include <random>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * 一个命名空间，包含所有与粒子生成有关的类。
   *
   */
  namespace Generators
  {
    /**
     * 一个在指定的 @p particle_reference_locations.
     * 的每个单元中生成粒子的函数 被添加到 @p
     * particle_handler 对象中的粒子总数是 @p triangulation
     * 的本地拥有的单元数乘以 @p particle_reference_locations.
     * 中的位置数 一个可选的 @p mapping 参数可用于从 @p
     * particle_reference_locations 映射到真实的粒子位置。
     * @param  triangulation 与 @p particle_handler. 相关的三角形。
     * @param  particle_reference_locations
     * 单元格中的一个位置矢量。
     * 粒子将在每个单元格中的这些位置生成。          @param
     * particle_handler
     * 粒子处理程序，它将负责生成的粒子的所有权。
     * @param  mapping
     * 一个可选的映射对象，用于将单元格中的参考位置映射到三角形的真实单元。如果没有提供映射，则假定有MappingQ1。
     *
     */
    template <int dim, int spacedim = dim>
    void
    regular_reference_locations(
      const Triangulation<dim, spacedim> &triangulation,
      const std::vector<Point<dim>> &     particle_reference_locations,
      ParticleHandler<dim, spacedim> &    particle_handler,
      const Mapping<dim, spacedim> &      mapping =
        (ReferenceCells::get_hypercube<dim>()
           .template get_default_linear_mapping<dim, spacedim>()));

    /**
     * 一个在单元格 @p cell
     * 中的随机位置生成一个粒子的函数，其索引为 @p id.
     * 。该函数希望有一个随机数发生器，以避免为每个粒子生成和销毁一个发生器的昂贵费用，并且可以选择考虑到单元格的映射。在该函数中实现的算法在
     * @cite GLHPW2018  中有描述。简而言之，该算法在 @p cell.
     * 的边界框内生成随机位置，然后反转映射以检查生成的粒子是否在细胞本身内。这确保了该算法即使对于非线性映射和扭曲的单元格也能产生统计学上的随机位置。然而，如果界线盒和单元格体积之间的比率变得非常大的话
     *
     * - 即细胞变得强烈变形，例如一个铅笔形状的细胞位于域的对角线上
     *
     * - 那么该算法可能会变得非常低效。    因此，在抛出错误信息之前，它只尝试找到一个固定次数的单元格位置。          @param[in]  cell 生成粒子的单元。          @param[in]  id 将被分配给新粒子的粒子索引。          @param[in,out]  random_number_generator 将用于创建粒子的随机数发生器。          @param[in]  mapping 一个可选的映射对象，用于将单元格中的参考位置映射到真实单元格。如果没有提供映射，则假定有MappingQ1。
     *
     */
    template <int dim, int spacedim = dim>
    Particle<dim, spacedim>
    random_particle_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
      const types::particle_index                                        id,
      std::mt19937 &                random_number_generator,
      const Mapping<dim, spacedim> &mapping =
        (ReferenceCells::get_hypercube<dim>()
           .template get_default_linear_mapping<dim, spacedim>()));

    /**
     * 一个根据所提供的概率密度函数 @p
     * probability_density_function.
     * 在域中随机生成粒子密度的函数，添加到 @p
     * particle_handler 对象中的粒子总数为 @p n_particles_to_create.
     * 一个可选的 @p mapping 参数可用于从 @p
     * particle_reference_locations 映射到真实粒子位置。
     * 该函数可以通过计算每个单元的概率密度函数的积分并相应地创建粒子来确定每个单元的粒子数（如果选项
     * @p random_cell_selection
     * 设置为false），或者它可以根据概率密度函数和单元大小随机选择单元（如果选项
     * @p random_cell_selection
     * 设置为true）。在这两种情况下，细胞内单个粒子的位置都是随机计算的。
     * 该函数中实现的算法在  @cite GLHPW2018  中有描述。
     * @param[in]  triangulation 与  @p particle_handler.   @param[in]
     * probability_density_function
     * 与非负值有关的函数，决定在这个位置生成粒子的概率密度。该函数不需要被归一化。
     * @param[in]  random_cell_selection
     * 一个bool，决定如何计算每个单元的粒子数（见上面的描述）。
     * @param[in]  n_particles_to_create
     * 这个函数将创建的粒子的数量。          @param[in,out]
     * particle_handler（粒子处理程序）
     * 将对生成的粒子拥有所有权的粒子处理程序。
     * @param[in]  mapping
     * 一个可选的映射对象，用于将单元格中的参考位置映射到三角形的实际单元。如果没有提供映射，则假定是MappingQ1。
     * @param[in]  random_number_seed
     * 一个可选的种子，决定随机数发生器的初始状态。使用相同的数字来获得可重复的粒子分布，或者使用一个变化的数字（例如基于系统时间）来为每次对该函数的调用产生不同的粒子分布。
     *
     */
    template <int dim, int spacedim = dim>
    void
    probabilistic_locations(
      const Triangulation<dim, spacedim> &triangulation,
      const Function<spacedim> &          probability_density_function,
      const bool                          random_cell_selection,
      const types::particle_index         n_particles_to_create,
      ParticleHandler<dim, spacedim> &    particle_handler,
      const Mapping<dim, spacedim> &      mapping =
        (ReferenceCells::get_hypercube<dim>()
           .template get_default_linear_mapping<dim, spacedim>()),
      const unsigned int random_number_seed = 5432);


    /**
     * 一个在DoFHandler的支持点位置生成粒子的函数，可能是基于与用于构建ParticleHandler的三角图不同的三角图。
     * 被添加到 @p particle_handler
     * 对象中的粒子总数是被传递的DoFHandler中位于三角形内且其组件位于ComponentMask内的dofs数量。
     * 这个函数使用insert_global_particles，因此可能会引起相当大的mpi通信开销。
     * 这个函数在  step-70  中使用。          @param[in]  dof_handler
     * 一个DOF处理程序，可能存在于另一个三角形上，用于确定粒子的位置。
     * @param[in]  global_bounding_boxes
     * 一个包含所有处理器的所有边界盒的向量。这个向量可以通过首先使用
     * 'GridTools::compute_mesh_predicate_bounding_box()'
     * 来建立，然后使用 'Utilities::MPI::all_gather().
     * 来收集所有的边界框。  @param[in,out]  particle_handler
     * 粒子处理程序，将负责生成的粒子的所有权。生成的粒子将被附加到当前由粒子处理程序拥有的粒子上。
     * @param[in]  mapping
     * 一个可选的映射对象，用于映射DOF位置。如果没有提供映射，则假定是MappingQ1。
     * @param[in]  components
     * 决定使用哪个支持点的dof_handler的子集来生成粒子的组件掩码。
     * @param[in]  properties 每个要插入的粒子的可选属性向量。
     *
     */
    template <int dim, int spacedim = dim>
    void
    dof_support_points(
      const DoFHandler<dim, spacedim> &dof_handler,
      const std::vector<std::vector<BoundingBox<spacedim>>>
        &                             global_bounding_boxes,
      ParticleHandler<dim, spacedim> &particle_handler,
      const Mapping<dim, spacedim> &  mapping =
        (ReferenceCells::get_hypercube<dim>()
           .template get_default_linear_mapping<dim, spacedim>()),
      const ComponentMask &                   components = ComponentMask(),
      const std::vector<std::vector<double>> &properties = {});

    /**
     * 一个在Triangulation的正交点位置生成粒子的函数。这个三角图可以不同于用于构造ParticleHandler的三角图。
     * 被添加到 @p particle_handler
     * 对象中的粒子总数是单元格的数量乘以粒子参考位置的数量，这些参考位置通常是用正交点构建的。
     * 这个函数使用insert_global_particles，因此可能引起相当大的mpi通信开销。
     * @param[in]  triangulation
     * 可能是不匹配的三角形，用于将粒子插入域中。
     * @param[in]  正交法
     * 一个正交法，其参考位置用于在单元中插入粒子。
     * @param[in]  global_bounding_boxes
     * 一个包含所有处理器的所有边界盒的向量。这个向量可以通过首先使用
     * 'GridTools::compute_mesh_predicate_bounding_box()'
     * 来建立，并使用 'Utilities::MPI::all_gather.
     * 来收集所有的边界框。  @param[in,out]  particle_handler
     * 粒子处理程序，它将负责生成的粒子的所有权。
     * @param[in]  mapping
     * 一个可选的映射对象，用于映射正交位置。如果没有提供映射，则假定有MappingQ1。
     * @param[in]  properties
     * 一个可选的矢量，用于插入每个粒子的属性矢量。
     *
     */
    template <int dim, int spacedim = dim>
    void
    quadrature_points(
      const Triangulation<dim, spacedim> &triangulation,
      const Quadrature<dim> &             quadrature,
      const std::vector<std::vector<BoundingBox<spacedim>>>
        &                             global_bounding_boxes,
      ParticleHandler<dim, spacedim> &particle_handler,
      const Mapping<dim, spacedim> &  mapping =
        (ReferenceCells::get_hypercube<dim>()
           .template get_default_linear_mapping<dim, spacedim>()),
      const std::vector<std::vector<double>> &properties = {});
  } // namespace Generators
} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif


