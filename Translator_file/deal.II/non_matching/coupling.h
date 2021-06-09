//include/deal.II-translator/non_matching/coupling_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_non_matching_coupling
#define dealii_non_matching_coupling

#include <deal.II/base/config.h>

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/affine_constraints.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个函数命名空间，提供处理两个没有对齐要求的网格的工具。
 * 通常这些函数允许对两个网格之间的实空间交点进行计算，例如表面积分和耦合矩阵的构造。
 *
 *
 */
namespace NonMatching
{
  /**
   * 为非匹配的、重叠的网格创建一个耦合稀疏模式。
   * 给出两个非匹配三角形，代表域 $\Omega$ 和 $B$  ，与 $B
   * \subseteq \Omega$ ，以及两个有限元空间 $V(\Omega) =
   * \text{span}\{v_i\}_{i=0}^n$ 和 $Q(B) = \text{span}\{w_j\}_{j=0}^m$
   * ，计算集合矩阵\f[ M_{ij} \dealcoloneq \int_{B} v_i(x) w_j(x) dx,
   * \quad i \in [0,n), j \in [0,m), \f]所需的稀疏模式，其中
   * $V(\Omega)$
   * ]是与传递给该函数的`space_dh'相关的有限元空间（如果在`space_comps'中指定，则为其一部分），而
   * $Q(B)$
   * 是与传递给该函数的`immersed_dh'相关的有限元空间（如果在`immersed_comps'中指定，则为其一部分）。
   * 稀疏度
   * "是通过定位正交点（由参考正交`四边形'获得）的位置来填补的，这些正交点定义在
   * $B$ 的元素上，相对于嵌入的三角形 $\Omega$
   * 。对于每个重叠的单元，对应于`space_dh`中的`space_comps`和`immersed_dh`中的`immersed_comps`的条目被添加到稀疏模式。
   * `space_comps'和`immersed_comps'掩码被假定为以同样的方式排序：`space_comps'的第一个分量将与`immersed_comps'的第一个分量耦合，第二个分量与第二个分量耦合，以此类推。如果两个掩码中的一个比另一个有更多的非零值，那么多余的成分将被忽略。
   * 如果域 $B$ 不在 $\Omega$
   * 之内，计算正交点位置的算法将抛出一个异常。特别要注意的是，这个函数只对`dim1`低于或等于`dim0`有感觉。一个静态断言可以保证实际情况是这样的。
   * 对于这两个空间，可以指定一个自定义的Mapping，默认为StaticMappingQ1。
   * 这个函数也可以并行工作，前提是浸入式三角形是
   * parallel::shared::Triangulation<dim1,spacedim>. 如果你使用浸入式
   * parallel::distributed::Triangulation<dim1,spacedim>.
   * ，会抛出一个异常，请看教程程序 step-60
   * 中关于如何使用这个函数的例子。
   *
   */
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number = double>
  void
  create_coupling_sparsity_pattern(
    const DoFHandler<dim0, spacedim> &space_dh,
    const DoFHandler<dim1, spacedim> &immersed_dh,
    const Quadrature<dim1> &          quad,
    Sparsity &                        sparsity,
    const AffineConstraints<number> & constraints = AffineConstraints<number>(),
    const ComponentMask &             space_comps = ComponentMask(),
    const ComponentMask &             immersed_comps = ComponentMask(),
    const Mapping<dim0, spacedim> &   space_mapping =
      StaticMappingQ1<dim0, spacedim>::mapping,
    const Mapping<dim1, spacedim> &immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * 与上述相同，但需要一个额外的 GridTools::Cache
   * 对象，而不是在内部创建一个。在这个版本的函数中，不能指定参数
   * @p space_mapping，因为它取自 @p cache 参数。
   *
   */
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename number = double>
  void
  create_coupling_sparsity_pattern(
    const GridTools::Cache<dim0, spacedim> &cache,
    const DoFHandler<dim0, spacedim> &      space_dh,
    const DoFHandler<dim1, spacedim> &      immersed_dh,
    const Quadrature<dim1> &                quad,
    Sparsity &                              sparsity,
    const AffineConstraints<number> &constraints = AffineConstraints<number>(),
    const ComponentMask &            space_comps = ComponentMask(),
    const ComponentMask &            immersed_comps = ComponentMask(),
    const Mapping<dim1, spacedim> &  immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping);


  /**
   * 为非匹配、重叠的网格创建耦合质量矩阵。
   * 给出两个非匹配三角形，代表域 $\Omega$ 和 $B$  ，与 $B
   * \subseteq \Omega$ ，以及两个有限元空间 $V(\Omega) =
   * \text{span}\{v_i\}_{i=0}^n$ 和 $Q(B) = \text{span}\{w_j\}_{j=0}^m$
   * ，计算耦合矩阵\f[ M_{ij} \dealcoloneq \int_{B} v_i(x) w_j(x) dx,
   * \quad i \in [0,n), j \in [0,m), \f]，其中 $V(\Omega)$
   * ]是与传递给该函数的`space_dh'相关的有限元空间（如果在`space_comps'中指定，则为其一部分），而
   * $Q(B)$
   * 是与传递给该函数的`immersed_dh'相关的有限元空间（如果在`immersed_comps`中指定，则为其一部分）。
   * 相应的稀疏性模式可以通过调用make_coupling_sparsity_pattern函数计算出来。矩阵的元素是通过定位定义在
   * $B$ 元素上的正交点相对于嵌入三角形 $\Omega$
   * 的位置来计算的。    假设 "space_comps "和 "immersed_comps
   * "掩码以同样的方式排序："space_comps "的第一个分量将与
   * "immersed_comps
   * "的第一个分量耦合，第二个分量与第二个分量耦合，以此类推。如果两个掩码中的一个比另一个有更多的非零项非零值，那么多余的分量将被忽略掉。
   * 如果域 $B$ 不在 $\Omega$
   * 之内，计算正交点位置的算法将抛出一个异常。特别要注意的是，这个函数只对`dim1`低于或等于`dim0`有意义。静态断言可以保证实际情况是这样的。
   * 对于这两个空间，可以指定一个自定义的Mapping，默认为StaticMappingQ1。
   * 这个函数也可以并行工作，前提是浸入式三角形是
   * parallel::shared::Triangulation<dim1,spacedim>. 如果你使用浸入式
   * parallel::distributed::Triangulation<dim1,spacedim>.
   * ，会抛出一个异常，请看教程程序 step-60
   * 中关于如何使用这个函数的例子。
   *
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const DoFHandler<dim0, spacedim> &                    space_dh,
    const DoFHandler<dim1, spacedim> &                    immersed_dh,
    const Quadrature<dim1> &                              quad,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints =
      AffineConstraints<typename Matrix::value_type>(),
    const ComponentMask &          space_comps    = ComponentMask(),
    const ComponentMask &          immersed_comps = ComponentMask(),
    const Mapping<dim0, spacedim> &space_mapping =
      StaticMappingQ1<dim0, spacedim>::mapping,
    const Mapping<dim1, spacedim> &immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * 与上述相同，但需要一个额外的 GridTools::Cache
   * 对象，而不是在内部创建一个。在这个版本的函数中，不能指定参数
   * @p space_mapping，因为它取自 @p cache 参数。
   *
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    const GridTools::Cache<dim0, spacedim> &              cache,
    const DoFHandler<dim0, spacedim> &                    space_dh,
    const DoFHandler<dim1, spacedim> &                    immersed_dh,
    const Quadrature<dim1> &                              quad,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints =
      AffineConstraints<typename Matrix::value_type>(),
    const ComponentMask &          space_comps    = ComponentMask(),
    const ComponentMask &          immersed_comps = ComponentMask(),
    const Mapping<dim1, spacedim> &immersed_mapping =
      StaticMappingQ1<dim1, spacedim>::mapping);

  /**
   * 为非匹配的独立网格创建一个耦合的稀疏模式，使用半径为ε的紧凑支持的卷积核。
   * 给出两个非匹配三角形，代表域 $\Omega^0$ 和 $\Omega^1$
   * ，都嵌入在 $\mathbb{R}^d$ 中，以及两个有限元空间
   * $V^0(\Omega^0) = \text{span}\{v_i\}_{i=0}^n$ 和 $V^1(\Omega^1) =
   * \text{span}\{w_\alpha\}_{\alpha=0}^m$ ，计算组装矩阵\f[
   * M_{i\alpha} \dealcoloneq \int_{\Omega^0} \int_{\Omega^1} v_i(x)
   * K^{\epsilon}(x-y) w_\alpha(y) dx \ dy, \quad i \in [0,n), \alpha \in
   * [0,m), \f]所需的稀疏模式，其中 $V^0(\Omega^0)$ 是与 @p dh0
   * 相关的有限元空间
   * ]传递给这个函数的有限元空间（如果在 @p comps0),
   * 中指定，则是其中一部分），而 $V^1(\Omega^1)$
   * 是传递给这个函数的 @p dh1 相关的有限元空间（如果在
   * @p comps1), 中指定，则是其中一部分）， $K^\epsilon$
   * 是一个从CutOffFunctionBase导出的函数，其紧凑支持包含在半径为
   * $\epsilon$ 的球中。    假设 @p comps0 和 @p comps1
   * 掩码以同样的方式排序： @p comps0 的第一个分量将与 @p
   * comps1,
   * 的第一个分量耦合，第二个分量与第二个分量耦合，以此类推。如果两个掩码中的一个比另一个有更多的有效成分，那么多余的成分将被忽略。
   * 对于这两个空间，可以指定一个自定义的Mapping，这两个空间的默认值是StaticMappingQ1。
   * 如果至少有一个三角形是
   * parallel::shared::Triangulation<dim1,spacedim>.
   * 类型的，这个函数也会并行工作，如果两个三角形都是
   * parallel::distributed::Triangulation<dim1,spacedim>.
   * 类型的，就会出现异常。 ]
   * 如果epsilon被设置为零，那么我们假定内核是狄拉克三角分布，并且调用被转发到这个名字空间中的同名方法，该方法不接受epsilon作为输入（但需要正交公式
   * @p quad
   * ）。在这种情况下，对两个空间需要有更多的限制条件。参见其他create_coupling_sparsity_pattern()函数的文档。
   *
   */
  template <int dim0,
            int dim1,
            int spacedim,
            typename Sparsity,
            typename Number = double>
  void
  create_coupling_sparsity_pattern(
    const double &                          epsilon,
    const GridTools::Cache<dim0, spacedim> &cache0,
    const GridTools::Cache<dim1, spacedim> &cache1,
    const DoFHandler<dim0, spacedim> &      dh0,
    const DoFHandler<dim1, spacedim> &      dh1,
    const Quadrature<dim1> &                quad,
    Sparsity &                              sparsity,
    const AffineConstraints<Number> &constraints0 = AffineConstraints<Number>(),
    const ComponentMask &            comps0       = ComponentMask(),
    const ComponentMask &            comps1       = ComponentMask());

  /**
   * 为非匹配的独立网格创建一个耦合质量矩阵，使用具有紧凑支持的卷积核。
   * 给出两个非匹配三角形，代表域 $\Omega^0$ 和 $\Omega^1$
   * ，都嵌入在 $\mathbb{R}^d$ 中，以及两个有限元空间
   * $V^0(\Omega^0) = \text{span}\{v_i\}_{i=0}^n$ 和 $V^1(\Omega^1) =
   * \text{span}\{w_\alpha\}_{\alpha=0}^m$ ，计算矩阵\f[ M_{i\alpha}
   * \dealcoloneq \int_{\Omega^0} \int_{\Omega^1} v_i(x) K^{\epsilon}(x-y)
   * w_\alpha(y) dx \ dy, \quad i \in [0,n), \alpha \in [0,m), \f] 其中
   * $V^0(\Omega^0)$ 是与 @p dh0  ]
   * 传递给这个函数的有限元空间（如果在  @p comps0),
   * 中指定，则为其中一部分），而  $V^1(\Omega^1)$
   * 是传递给这个函数的  @p dh1  相关的有限元空间（如果在
   * @p comps1),  中指定，则为其中一部分），而  $K^\epsilon$
   * 是一个从 CutOffFunctionBase
   * 派生的函数，紧凑支持包含在半径为  $\epsilon$
   * 的球中。
   * 相应的稀疏性模式可以通过调用make_coupling_sparsity_pattern()函数来计算。
   * 假设 @p comps0 和 @p comps1 掩码以同样的方式排序： @p
   * comps0 的第一个分量将与 @p comps1,
   * 的第一个分量耦合，第二个分量与第二个分量耦合，以此类推。如果两个掩码中的一个比另一个有更多的有效成分，那么多余的成分将被忽略。
   * 对于这两个空间，可以指定一个自定义的Mapping，这两个空间的默认值是StaticMappingQ1。
   * 这个函数也可以并行工作，前提是两个三角形中的一个是
   * parallel::shared::Triangulation<dim1,spacedim>.
   * 类型，如果两个三角形都是
   * parallel::distributed::Triangulation<dim1,spacedim>.
   * 类型，则会出现异常。 参数 @p epsilon
   * 用于设置用于计算卷积的截断函数的大小。如果epsilon被设置为0，那么我们假定内核是狄拉克德尔塔分布，并且调用被转发到这个名字空间中的同名方法，该方法不接受epsilon作为输入。
   *
   */
  template <int dim0, int dim1, int spacedim, typename Matrix>
  void
  create_coupling_mass_matrix(
    Functions::CutOffFunctionBase<spacedim> &             kernel,
    const double &                                        epsilon,
    const GridTools::Cache<dim0, spacedim> &              cache0,
    const GridTools::Cache<dim1, spacedim> &              cache1,
    const DoFHandler<dim0, spacedim> &                    dh0,
    const DoFHandler<dim1, spacedim> &                    dh1,
    const Quadrature<dim0> &                              quadrature0,
    const Quadrature<dim1> &                              quadrature1,
    Matrix &                                              matrix,
    const AffineConstraints<typename Matrix::value_type> &constraints0 =
      AffineConstraints<typename Matrix::value_type>(),
    const ComponentMask &comps0 = ComponentMask(),
    const ComponentMask &comps1 = ComponentMask());
} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif


