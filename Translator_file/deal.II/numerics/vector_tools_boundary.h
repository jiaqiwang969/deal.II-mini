//include/deal.II-translator/numerics/vector_tools_boundary_0.txt
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


#ifndef dealii_vector_tools_boundary_h
#define dealii_vector_tools_boundary_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/hp/mapping_collection.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
namespace hp
{
  template <int dim>
  class QCollection;
} // namespace hp

namespace VectorTools
{
  /**
   * @name  内插和投影
   *
   */
  //@{

  /**
   * 计算对应于施加Dirichlet边界条件的解决方案的约束。 这个函数通过在边界周围插值，创建一个受迪里希特边界条件约束的自由度地图，以及要分配给它们的相应数值。对于边界上的每个自由度，如果它的索引已经存在于
   * @p boundary_values 中，那么它的边界值将被覆盖，否则将在
   * @p boundary_values.
   * 中插入一个具有适当索引和边界值的新条目。 参数 @p
   * function_map
   * 提供了一个由该函数处理的边界指标列表和相应的边界值函数。
   * 该图的键对应于面的编号 @p boundary_id 。
   * numbers::internal_face_boundary_id
   * 是这个键的非法值，因为它是为内部面保留的。关于如何用非空地图使用这个参数的例子，请看
   * step-16 的教程程序。    最后一个参数 @p component_mask
   * 中的标志表示有限元空间的哪些部分应被插值。如果保留默认值（即一个空数组），所有组件都被插值。如果它与默认值不同，则假定条目数等于边界函数和有限元中的分量数，给定边界函数中的那些分量将被用于在分量掩码中设置的相应标志。另见
   * @ref GlossComponentMask  。
   * 举个例子，假设你正在求解2d中的斯托克斯方程，变量为
   * $(u,v,p)$
   * ，你只想插值速度的边界值，那么分量掩码应该对应于
   * <code>(true,true,false)</code>  。
   * @note  无论是否指定了分量掩码， @p function_map
   * 中函数的分量数必须与 @p dof.
   * 使用的有限元的分量数一致。 ]
   * 换句话说，对于上面的例子，你需要提供一个有3个分量（两个速度和压力）的Function对象，尽管你只对其中的前两个感兴趣。
   * interpolate_boundary_values()然后会调用这个函数，在每个插值点获得3个值的向量，但只取前两个，舍弃第三个。换句话说，你可以自由地在
   * Function::vector_value,
   * 返回的向量的第三个分量中返回你喜欢的东西，但Function对象必须说明它有3个分量。
   * 如果使用的有限元的形状函数在一个以上的分量中是非零的（用deal.II的话说：它们是非原始的），那么这些分量目前不能用于内插边界值。
   * 因此，与这些非正则形状函数的分量相对应的分量掩码中的元素必须是
   * @p false.  更多信息请参见该命名空间的一般文档。
   * @note  当求解带有边界条件的偏微分方程
   * $u|_{\partial\Omega}=g$
   * （或边界的部分*）时，那么这个边界条件一般不能用
   * $u_h|_{\partial\Omega}=g$
   * 形式的有限元准确满足。这是因为函数 $g$
   * 一般不是多项式，而 $u_h|_{\partial\Omega}$
   * 在位于边界的网格的每个面上都是*项多项式。换句话说，一般来说，不可能施加*这样的边界条件；然而，可以*做的是施加@f[
   * u_h|_{\partial\Omega}=I_h^{\partial\Omega} g, @f]，其中
   * $I_h^{\partial\Omega} g$
   * 是一个函数，在位于边界的有限元空间的每个节点上等于
   * $g$ ，并且在两者之间是片断多项式。换句话说，
   * $I_h^{\partial\Omega}$ 是一个内插算子*， $I_h^{\partial\Omega} g$
   * 是内插的边界值
   *
   * --因此而得名。使用 $I_h^{\partial\Omega} g$ 而不是 $g$
   * 作为边界值会带来一个额外的误差（与使用正交引入一个额外的误差相比，能够准确计算弱形式的积分）。在大多数情况下，这个额外的误差与有限元方法中的其他误差项相同，尽管在测量
   * $L^2$
   * 准则中的误差时有一些微妙的差别。关于一些细节，请参见
   * @cite Bartels2004  。
   * @note  使用内插法的一个替代方法，@f[
   * u_h|_{\partial\Omega}=I_h^{\partial\Omega} g @f]是使用边界值 $g$
   * 对边界上的有限元空间的投影*：@f[
   * u_h|_{\partial\Omega}=\Pi_h^{\partial\Omega} g. @f]
   * 投影可以使用project_boundary_values()函数。使用投影可能有一些理论上的优势（参见
   * @cite Bartels2004
   * ），但有一个实际的缺点，即计算投影比计算插值要昂贵得多，因为后者可以一次完成一个面，而投影则需要解决整个边界上的问题。另一方面，插值只适用于
   * "节点式
   * "有限元空间（如FE_Q，但不包括FE_Q_Hierarchical），而投影则始终适用。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 像之前的函数一样，但采取一个映射集合，以配合具有hp能力的DoFHandler对象。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 和上面的函数一样，但只取一对边界指标和相应的边界函数。同样的评论适用于前面的函数，特别是关于组件掩码的使用和函数对象的要求大小。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &             mapping,
    const DoFHandler<dim, spacedim> &          dof,
    const types::boundary_id                   boundary_component,
    const Function<spacedim, number> &         boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 像之前的函数一样，但要取一个映射集合，以配合具有hp-capabilities的DoFHandler对象。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const types::boundary_id                    boundary_component,
    const Function<spacedim, number> &          boundary_function,
    std::map<types::global_dof_index, number> & boundary_values,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 调用另一个interpolate_boundary_values()函数，见上文，使用<tt>mapping=MappingQGeneric  @<dim,spacedim@>(1)</tt>.  与前一个函数的注释相同，特别是关于组件掩码的使用和函数对象的要求大小。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim> &          dof,
    const types::boundary_id                   boundary_component,
    const Function<spacedim, number> &         boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());


  /**
   * 调用另一个interpolate_boundary_values()函数，见上文，<tt>mapping=MappingQGeneric
   * @<dim,spacedim@>(1)</tt>.
   * 与前一个函数的评论相同，特别是关于使用分量掩码和函数对象的要求大小。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());


  /**
   * 将由于Dirichlet边界条件引起的（代数）约束插入AffineConstraints
   * @p constraints.
   * 这个函数识别受Dirichlet边界条件约束的自由度，将它们添加到
   * @p constraints
   * 中受约束的DoF列表中，并将各自的不均匀性设置为围绕边界插值的值。如果这个例程遇到一个已经被约束的DoF（例如被悬挂的节点约束，见下文，或任何其他类型的约束，例如来自周期性边界条件的约束），约束的旧设置（条目被约束的DoF，不均匀性）被保留，不会发生什么。
   * @note
   * 当在一个AffineConstraints对象中结合具有悬挂节点约束和边界条件的自适应细化网格时，悬挂节点约束应该总是首先被设置，然后是边界条件，因为边界条件不会在已经被约束的自由度的第二次操作中被设置。这可以确保离散化保持所需的一致性。参见
   * @ref constraints 模块中关于冲突约束的讨论。
   * 这个函数与上面的函数基本等同，只是它将其结果放入AffineConstraint对象中，而不是
   * `std::map`.  更多评论见上面的函数。
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask &      component_mask = ComponentMask());

  /**
   * 和前面的函数一样，但要取一个映射集合，以配合具有hp-capabilities的DoFHandler对象。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask &      component_mask = ComponentMask());

  /**
   * 与上述函数相同，但只取一对边界指标和相应的边界函数。同样的评论适用于前面的函数，特别是关于组件掩码的使用和函数对象的要求大小。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &    mapping,
    const DoFHandler<dim, spacedim> & dof,
    const types::boundary_id          boundary_component,
    const Function<spacedim, number> &boundary_function,
    AffineConstraints<number> &       constraints,
    const ComponentMask &             component_mask = ComponentMask());

  /**
   * 像之前的函数一样，但要取一个映射集合，以配合具有hp-capabilities的DoFHandler对象。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const types::boundary_id                    boundary_component,
    const Function<spacedim, number> &          boundary_function,
    AffineConstraints<number> &                 constraints,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * 调用另一个interpolate_boundary_values()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim,spacedim@>(1)</tt>.
   * 与前一个函数的注释相同，特别是关于组件掩码的使用和函数对象的要求大小。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim> & dof,
    const types::boundary_id          boundary_component,
    const Function<spacedim, number> &boundary_function,
    AffineConstraints<number> &       constraints,
    const ComponentMask &             component_mask = ComponentMask());


  /**
   * 调用另一个interpolate_boundary_values()函数，见上文，<tt>mapping=MappingQGeneric
   * @<dim,spacedim@>(1)</tt>.
   * 与前一个函数的注释相同，特别是关于组件掩码的使用和函数对象的要求大小。
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask &      component_mask = ComponentMask());


  /**
   * 投射一个函数或一组函数到域的边界。
   * 换句话说，计算以下问题的解决方案：找到 $u_h \in V_h$
   * （其中 $V_h$
   * 是由该函数的DoFHandler参数所代表的有限元空间），以便
   * @f{align*}{
   * \int_{\Gamma} \varphi_i u_h
   * = \sum_{k \in {\cal K}} \int_{\Gamma_k} \varphi_i f_k,
   * \qquad \forall \varphi_i \in V_h
   * @f}
   * 其中 $\Gamma = \bigcup_{k \in {\cal K}} \Gamma_k$ ,  $\Gamma_k \subset
   * \partial\Omega$ ,  $\cal K$ 是指数集， $f_k$
   * 是该函数的函数图参数 @p boundary_values
   * 所代表的相应边界函数，积分用正交法求值。这个问题在内部有一个非唯一的解决方案，但对于边界部分的自由度，
   * $\Gamma$ ，它是定义良好的，我们对其进行积分。
   * $u_h|_\Gamma$
   * 的值，即这个函数沿边界的自由度的节点值，就是这个函数所计算的。
   * 如果这个函数与 $H_{div}$
   * 符合的有限元空间一起使用，则会计算出一个不同的问题的解决方案，即。找到
   * $\vec{u}_h \in V_h \subset H(\text{div}; \Omega)$ ，以便
   * @f{align*}{
   * \int_{\Gamma} (\vec{\varphi}_i \cdot \vec{n}) (\vec{u}_h \cdot \vec{n})
   * = \sum_{k \in {\cal K}} \int_{\Gamma_k} (\vec{\varphi}_i \cdot \vec{n})
   * (\vec{f}_k \cdot \vec{n}),
   * \qquad \forall \vec{\varphi_i} \in V_h,
   * @f}
   * 其中 $\vec{n}$ 是一个外向法向量。    如果与 $H_\text{curl}$ 符合的元素一起使用，这个函数会抛出一个异常，所以应该使用project_boundary_values_curl_conforming_l2()来代替。      @param[in]  映射 将用于沿边界整合所需的转换的映射。    @param[in]  dof 描述有限元空间和自由度编号的DoFHandler。    @param[in]  boundary_functions 从边界指标到函数指针的映射，这些函数描述了边界上标有该边界指标的那些部分的期望值（见 @ref GlossBoundaryIndicator  "边界指标"
   * ）。
   * 投影只发生在边界的那些部分，这些部分的指标在这个地图中被代表。
   * @param[in]  q
   * 用于计算质量矩阵和投影的右手边所需的积分的面正交。
   * @param[out]  boundary_values
   * 这个函数的结果。它是一个包含边界上所有自由度指数的地图（如
   * @p boundary_functions)
   * 中的边界部分所涵盖的，以及该自由度的计算dof值。对于边界上的每个自由度，如果它的索引已经存在于
   * @p boundary_values 中，那么它的边界值将被覆盖，否则将在
   * @p boundary_values.
   * 中插入一个具有适当索引和边界值的新条目。 ]
   * component_mapping
   * 有时，将一个矢量值函数投射到有限元空间的一部分是很方便的（例如，将一个
   * <code>dim</code> 分量的函数投射到斯托克斯问题的
   * <code>dim+1</code>
   * 分量DoFHandler的速度分量）。为了允许这一点，这个参数允许分量被重新映射。如果矢量不是空的，它必须为
   * @p dof.
   * 中使用的有限元的每个矢量分量有一个条目，这个条目是
   * @p boundary_functions 中的分量编号，应该用于 @p dof.
   * 中的这个分量。 默认情况下，不应用重新映射。
   * @note
   * 使用投影*而不是边界值的插值*，在实践中没有什么区别。也就是说，计算投影的计算成本要高得多，因为它需要解决一个耦合边界上所有未知数的问题，而内插法则是一次只计算一个面。另一方面，插值只适用于
   * "节点型
   * "有限元空间（如FE_Q，但不包括FE_Q_Hierarchical），而投影则始终适用。(关于一些更多的理论考虑，见上面第一个interpolate_boundary_values()函数的文档)。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_functions,
    const Quadrature<dim - 1> &                q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * 调用project_boundary_values()函数，见上文，<tt>mapping=MappingQGeneric
   * @<dim,spacedim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_function,
    const Quadrature<dim - 1> &                q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * 和上面一样，但有hp-capabilities。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_functions,
    const hp::QCollection<dim - 1> &           q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * 调用project_boundary_values()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim,spacedim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_function,
    const hp::QCollection<dim - 1> &           q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * 将一个函数投影到域的边界，使用给定的面的正交公式。该函数识别受Dirichlet边界条件约束的自由度，将其添加到
   * @p constraints
   * 中的受约束自由度列表中，并将各自的不均匀性设置为投影操作产生的值。如果这个例程遇到一个已经被约束的DoF（例如被悬挂的节点约束，见下文，或任何其他类型的约束，例如来自周期性边界条件的约束），约束的旧设置（条目被约束的DoF，不均匀性）被保留，不会发生什么。
   * @note
   * 当在一个AffineConstraints对象中结合具有悬挂节点约束和边界条件的自适应细化网格时，悬挂节点约束应该总是首先被设置，然后是边界条件，因为边界条件不会在已经被约束的自由度的第二次操作中被设置。这可以确保离散化保持所需的一致性。参见
   * @ref constraints 模块中关于冲突约束的讨论。    如果 @p
   * component_mapping 为空，则假定 @p boundary_function 的分量数与
   * @p dof. 所使用的有限元的分量数一致。
   * 在1d中，投影等于插值。因此，interpolate_boundary_values被调用。
   * @arg   @p component_mapping:  如果 @p boundary_functions 和 @p dof
   * 中的组件不重合，这个向量允许它们被重新映射。如果这个向量不是空的，它必须有一个条目代表
   * @p 中的每个组件。这个条目是 @p boundary_functions
   * 中的分量编号，应该用于 @p dof. 中的这个分量。
   * 默认情况下，不应用重映射。
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        boundary_functions,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping = {});

  /**
   * 调用project_boundary_values()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim,spacedim@>(1)</tt>.  。
   * @ingroup constraints
   *
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        boundary_function,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping = {});

  /**
   * 这个函数是project_boundary_values_curl_conforming函数的一个更新版本。其目的是修复在使用前一个函数时与非矩形几何体（即具有非矩形面的元素）一起使用的一个问题。所用的L2投影方法取自PD
   * Ledger, K Morgan and O Hassan的论文 "Electromagnetic Scattering
   * simulation using an H (curl) conforming hp-finite element method in three
   * dimensions" ( Int. J. Num. Meth. Fluids, Volume 53, Issue 8, pages
   * 1267-1296) 。    该函数将计算对应于
   * $\vec{n}\times\vec{E}=\vec{n}\times\vec{F}$
   * 形式的Dirichlet边界条件的约束，即 $\vec{E}$ 和 $f$
   * 的切向分量应相吻合。    <h4>Computing constraints</h4>
   * 为了计算约束条件，我们使用基于上述论文的投影方法。在二维中，对于基于边缘的形状函数，无论有限元的顺序如何，这都是在单一阶段完成的。在三维中，这是在两个阶段完成的，首先是边，然后是面。
   * 对于每个单元，每个边 $e$ 通过解决线性系统 $Ax=b$
   * 进行投影，其中 $x$ 是边上自由度的约束向量， $A_{ij} =
   * \int_{e} (\vec{s}_{i}\cdot\vec{t})(\vec{s}_{j}\cdot\vec{t}) dS$  $b_{i} =
   * \int_{e} (\vec{s}_{i}\cdot\vec{t})(\vec{F}\cdot\vec{t}) dS$ 是
   * $\vec{s}_{i}$ 形状函数， $\vec{t}$ 是切线向量。
   * 一旦所有的边缘约束， $x$
   * 被计算出来，我们可以用类似的方式计算面的约束，同时考虑到边缘的残差。
   * 对于单元格上的每个面， $f$ ，我们解决线性系统 $By=c$
   * ，其中 $y$ 是面的自由度约束矢量， $B_{ij} = \int_{f}
   * (\vec{n} \times \vec{s}_{i}) \cdot (\vec{n} \times \vec{s}_{j}) dS$
   * $c_{i} = \int_{f} (\vec{n} \times \vec{r}) \cdot (\vec{n} \times
   * \vec{s}_i) dS$ 和 $\vec{r} = \vec{F}
   *
   * - \sum_{e
   * \in f} \sum{i \in e} x_{i}\vec{s}_i$ ，是边缘残差。
   * 然后在解决方案 $x$ 和 $y$ 中给出所产生的约束。
   * 如果AffineConstraints  @p constraints
   * 中包含数值或其他约束，那么如果要使用的边界部分的节点已经在约束列表中，那么新的约束将被添加或覆盖旧的约束。这可以通过使用不均匀约束来处理。请注意，当结合自适应网格和这种约束时，应该首先设置Dirichlet条件，然后通过悬挂节点约束完成，以确保离散化保持一致。参见
   * @ref constraints  模块中关于冲突约束的讨论。
   * <h4>Arguments to this function</h4>
   * 该函数明确用于FE_Nedelec元素，或包含FE_Nedelec元素的FESystem元素。如果在调用该函数时使用任何其他有限元，则会产生一个异常。用户在使用此函数时必须保证FESystem元素的正确设置，因为在这种情况下不可能进行检查。
   * 这个函数的第二个参数表示有限元的第一个矢量分量，对应于你希望约束的矢量函数。例如，如果我们正在求解三维麦克斯韦方程，并且有分量
   * $(E_x,E_y,E_z,B_x,B_y,B_z)$ ，我们想要边界条件
   * $\vec{n}\times\vec{B}=\vec{n}\times\vec{f}$ ，那么 @p
   * 第一_向量_分量将是3。在这个例子中， @p boundary_function
   * 必须返回6个分量，前3个对应于 $\vec{E}$ ，后3个对应于
   * $\vec{B}$  。隐含地假定矢量正好有 <code>dim</code>
   * 个分量，这些分量的排序方式与我们通常对坐标方向的排序方式相同，即
   * $x$  -， $y$  -，最后是 $z$  -分量。    参数 @p
   * boundary_component 对应于面的编号 @p 边界_id。
   * numbers::internal_face_boundary_id
   * 是一个非法值，因为它是为内部面保留的。
   * 最后一个参数被表示为计算边界点的法向量 $\vec{n}$ 。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim, typename number>
  void
  project_boundary_values_curl_conforming_l2(
    const DoFHandler<dim, dim> & dof_handler,
    const unsigned int           first_vector_component,
    const Function<dim, number> &boundary_function,
    const types::boundary_id     boundary_component,
    AffineConstraints<number> &  constraints,
    const Mapping<dim> &         mapping);


  /**
   * project_boundary_values_curl_conforming_l2（上文）的hp-namespace版本。
   * @ingroup constraints
   *
   */
  template <int dim, typename number>
  void
  project_boundary_values_curl_conforming_l2(
    const DoFHandler<dim, dim> &           dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, number> &          boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<number> &            constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection =
      hp::StaticMappingQ1<dim>::mapping_collection);


  /**
   * 计算对应于 $\vec{n}^T\vec{u}=\vec{n}^T\vec{f}$
   * 形式的边界条件的约束，即解 $u$ 和给定 $f$
   * 的法线分量应重合。函数 $f$ 由 @p boundary_function
   * 给出，所得到的约束被添加到 @p
   * 约束中，用于具有边界指标的面 @p boundary_component.
   * 这个函数是明确写给FE_RaviartThomas元素使用的。因此，如果它与其他有限元一起被调用，就会抛出一个异常。
   * 如果AffineConstraints对象 @p constraints
   * 之前包含数值或其他约束，那么如果要使用的边界部分的一个节点已经在约束列表中，那么就会添加新的约束或覆盖旧的约束。这可以通过使用不均匀约束来处理。请注意，当结合自适应网格和这种约束时，应该首先设置Dirichlet条件，然后通过悬挂节点约束完成，以确保离散化保持一致。参见
   * @ref constraints  模块中关于冲突约束的讨论。    参数 @p
   * first_vector_component
   * 表示有限元中的第一个矢量分量，对应于你要约束的矢量函数
   * $\vec{u}$ 。隐含地假定矢量正好有 <code>dim</code>
   * 个分量，这些分量的排序方式与我们通常对坐标方向的排序方式相同，即
   * $x$  -， $y$  -，最后是 $z$  -分量。    参数 @p
   * boundary_component 对应的是应用边界条件的面的 @p boundary_id
   * 。    numbers::internal_face_boundary_id
   * 是一个非法值，因为它被保留给内部面。 @p mapping
   * 用于计算边界点的法向量 $\vec{n}$ 。    <h4>Computing
   * constraints</h4>
   * 为了计算约束条件，我们在位于边界的每个面上使用Brezzi,
   * Fortin (Mixed and Hybrid Finite Element Methods, Springer,
   * 1991)中提出的插值运算。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  template <int dim>
  void
  project_boundary_values_div_conforming(
    const DoFHandler<dim, dim> & dof_handler,
    const unsigned int           first_vector_component,
    const Function<dim, double> &boundary_function,
    const types::boundary_id     boundary_component,
    AffineConstraints<double> &  constraints,
    const Mapping<dim> &         mapping);

  /**
   * 与上述hp-namespace的情况相同。
   * @ingroup constraints   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇表条目"
   *
   */
  template <int dim>
  void
  project_boundary_values_div_conforming(
    const DoFHandler<dim, dim> &           dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, double> &          boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<double> &            constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection =
      hp::StaticMappingQ1<dim>::mapping_collection);

  // @}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_boundary_h


