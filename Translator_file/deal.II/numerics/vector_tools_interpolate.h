//include/deal.II-translator/numerics/vector_tools_interpolate_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_vector_tools_interpolate_h
#define dealii_vector_tools_interpolate_h

#include <deal.II/base/config.h>

#include <deal.II/fe/component_mask.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
class DoFHandler;
template <typename number>
class FullMatrix;
template <int dim, typename Number>
class Function;
template <class MeshType>
class InterGridMap;
template <int dim, int spacedim>
class Mapping;

namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
}

namespace VectorTools
{
  /**
   * @name  内插和投影
   *
   */
  //@{

  /**
   * 计算 @p function
   * 在支持点的内插到由给定DoFHandler参数初始化的Triangulation和FiniteElement对象描述的有限元空间。假设
   * @p function 的分量数与 @p dof.
   * 使用的有限元的分量数相匹配。注意，你可能要在事后用空间
   * @p dof
   * 的悬空节点调用<tt>hanging_nodes.distribution(vec)</tt>，以便使结果再次连续。
   * 更多信息请参见该命名空间的一般文档。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  interpolate(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType &                                               vec,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 与上述相同，但在hp-context中。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  interpolate(
    const hp::MappingCollection<dim, spacedim> &               mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType &                                               vec,
    const ComponentMask &component_mask = ComponentMask());


  /**
   * 用<tt>mapping=MappingQGeneric  @<dim,spacedim@>(1)</tt>. 调用上面的
   * @p interpolate() 函数。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  interpolate(
    const DoFHandler<dim, spacedim> &                          dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType &                                               vec,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 对不同的有限元空间进行插值。从 @p
   * dof_1 代表的FE空间到FE空间 @p dof_2. 上的向量 @p data_2
   * 执行向量 @p data_1 的插值（假设是重影，见 @ref
   * GlossGhostedVector ），每个单元上的插值由矩阵 @p transfer.
   * 表示，迄今为止，弯曲的边界被忽略了。
   * 请注意，你可能要在之后调用<tt>hanging_nodes.distribution(data_2)</tt>与空间
   * @p dof_2 中的悬挂节点，以使结果再次连续。
   * @note
   * 这个模板的实例化提供给一些矢量类型（见命名空间的一般文档），但只有InVector和OutVector的同一个矢量。其他组合必须通过手工实例化。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolate(const DoFHandler<dim, spacedim> &dof_1,
              const DoFHandler<dim, spacedim> &dof_2,
              const FullMatrix<double> &       transfer,
              const InVector &                 data_1,
              OutVector &                      data_2);

  /**
   * 这个函数是该系列中第一个interpolate()函数的一种概括或修改。它将一组函数插值到由DoFHandler参数定义的有限元空间上，根据每个单元的材料ID（见 @ref
   * GlossMaterialId ）来决定在每个单元上使用哪个函数。
   * @param[in]  mapping
   * 用来确定要评估函数的支持点位置的映射。    @param[in]
   * dof_handler
   * 用Triangulation和FiniteElement对象初始化的DoFHandler，它定义了有限元空间。
   * @param[in]  function_map A  std::map
   * 反映了那些应该被内插的单元上的材料ID与要内插到有限元空间上的函数之间的对应关系。
   * @param[out]  dst 全局有限元向量，保存内插值的输出。
   * @param[in]  component_mask 应被内插的组件的掩码。
   * @note  如果算法遇到一个单元，其材料ID没有列在给定的
   * @p function_map, 中，那么 @p dst
   * 将不会在输出向量的各自自由度中被更新。例如，如果
   * @p dst
   * 被初始化为零，那么在调用此函数后，那些对应于遗漏的材料ID的零将仍然保留在
   * @p dst 中。
   * @note
   * 位于不同材料id的单元之间的面的自由度将由在本函数中实现的各单元循环中最后调用的单元获得其值。由于单元格的顺序在某种程度上是任意的，你无法控制它。然而，如果你想控制单元格被访问的顺序，让我们看一下下面的例子。让
   * @p u 是一个感兴趣的变量，它被一些CG有限元所近似。让
   * @p 0,  @p 1 和 @p 2
   * 为三角形上的单元的材料id。让0：0.0，1：1.0，2：2.0是你想传递给这个函数的整个
   * @p function_map ，其中 @p key 是材料ID， @p value 是 @p u. 的值
   * 通过使用整个 @p function_map
   * ，你并不真正知道哪些值会被分配给面的DoF。另一方面，如果你把整个
   * @p
   * function_map分成三个较小的独立对象0：0.0和1：1.0以及2：2.0，并对这个函数进行三次不同的调用，分别传递这些对象（顺序取决于你想在单元间得到什么），那么后面的每次调用将重写前一次的单元间
   * @p  dofs。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  interpolate_based_on_material_id(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::map<types::material_id,
                   const Function<spacedim, typename VectorType::value_type> *>
      &                  function_map,
    VectorType &         dst,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 计算 @p dof1-function  @p u1 到 @p dof2-function  @p u2,
   * 的插值，其中 @p dof1 和 @p dof2
   * 代表具有共同粗网格的不同三角形。
   * dof1和dof2需要具有相同的有限元离散化。
   * 请注意，对于有悬挂节点的网格上的连续元素（即局部细化网格），这个函数并不能给出预期的输出。
   * 事实上，由于局部单元插值的原因，所产生的输出矢量不一定尊重悬空节点的连续性要求。
   * 对于这种情况（有悬挂节点的网格上的连续元素），请使用带有额外AffineConstraints参数的interpolate_to_different_mesh函数，见下文，或者通过调用悬挂节点约束对象的
   * @p AffineConstraints::distribute 函数使场符合要求。
   * @note 这个函数与 parallel::distributed::Triangulation,
   * 一起工作，但只有在两个网格的平行分区相同的情况下（见
   * parallel::distributed::Triangulation<dim>::no_automatic_repartitioning
   * 标志）。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  interpolate_to_different_mesh(const DoFHandler<dim, spacedim> &dof1,
                                const VectorType &               u1,
                                const DoFHandler<dim, spacedim> &dof2,
                                VectorType &                     u2);

  /**
   * 计算 @p dof1-function  @p u1 到 @p dof2-function  @p u2,
   * 的插值，其中 @p dof1 和 @p dof2
   * 代表具有共同粗网格的不同三角形。Dof1和Dof2需要具有相同的有限元离散化。
   * @p constraints  是对应于 @p
   * dof2的悬挂节点约束对象。当插值到具有悬挂节点的网格（局部细化网格）上的连续元素时，这个对象特别重要。
   * 没有它
   *
   * - 由于单元格内插的原因
   *
   * - 产生的输出矢量不一定尊重悬挂节点的连续性要求。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  interpolate_to_different_mesh(
    const DoFHandler<dim, spacedim> &                         dof1,
    const VectorType &                                        u1,
    const DoFHandler<dim, spacedim> &                         dof2,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType &                                              u2);

  /**
   * 与上述函数相同，但直接接受InterGridMap对象作为参数。对于同时插值几个向量很有用。
   * @p intergridmap 必须通过 InterGridMap::make_mapping
   * 从源DoFHandler指向目的DoFHandler进行初始化。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  interpolate_to_different_mesh(
    const InterGridMap<DoFHandler<dim, spacedim>> &           intergridmap,
    const VectorType &                                        u1,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType &                                              u2);

  //@}

  /**
   * 几何插值
   *
   */
  //@{
  /**
   * 给定一个至少包含一个间隔向量场的DoFHandler，这个函数在FE_Q()有限元的支持点上插值三角图，该有限元的度数与所需组件的度数相同。
   * 曲线流形得到尊重，产生的VectorType将是几何上一致的。产生的映射保证在FE_Q()有限元的支持点上是插值的，其程度与所需组件的程度相同。
   * 如果底层有限元是FE_Q(1)^spacedim，那么产生的 @p VectorType
   * 是三角结构顶点的有限元场表示。
   * 可选的ComponentMask参数可以用来指定FiniteElement的哪些组件来描述几何。如果在构造时没有指定掩码，那么将使用默认的构造掩码，这将被解释为假设FiniteElement的第一个`spacedim'分量来代表问题的几何形状。
   * 这个函数只对指定组件是原始的FiniteElements实现。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  get_position_vector(const DoFHandler<dim, spacedim> &dh,
                      VectorType &                     vector,
                      const ComponentMask &            mask = ComponentMask());

  /**
   * 与上述函数类似，但也以 @p mapping 为参数。
   * 如果映射的度数低于DoFHandler @p dh,
   * 中有限元素的度数，这将在流形指定的真实几何之间引入一个额外的近似值，但更重要的是它允许为不保留顶点位置的映射（如欧拉映射）填充位置向量。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  get_position_vector(const Mapping<dim, spacedim> &   mapping,
                      const DoFHandler<dim, spacedim> &dh,
                      VectorType &                     vector,
                      const ComponentMask &            mask = ComponentMask());

  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_interpolate_h


