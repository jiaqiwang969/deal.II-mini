//include/deal.II-translator/numerics/vector_tools_project_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_vector_tools_project_h
#define dealii_vector_tools_project_h


#include <deal.II/base/config.h>

#include <functional>
#include <memory>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim, typename number, typename VectorizedArrayType>
class MatrixFree;
template <int dim>
class Quadrature;
template <int dim>
class QGauss;
template <typename Number, std::size_t width>
class VectorizedArray;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
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
   * 计算 @p function
   * 对有限元空间的投影。换句话说，给定一个函数
   * $f(\mathbf x)$ ，当前函数计算一个有限元函数 $f_h(\mathbf
   * x)=\sum_j F_j \varphi_j(\mathbf x)$ ，其特征是节点值 $F$
   * 的（输出）矢量满足方程
   * @f{align*}{
   * (\varphi_i, f_h)_\Omega = (\varphi_i,f)_\Omega
   * @f}
   * 对于所有测试函数  $\varphi_i$
   * 。这需要解决一个涉及质量矩阵的线性系统，因为上面的方程等同于线性系统
   * @f{align*}{
   * \sum_j (\varphi_i, \varphi_j)_\Omega F_j = (\varphi_i,f)_\Omega
   * @f}
   * 这也可以写成  $MF = \Phi$  与  $M_{ij} = (\varphi_i,
   * \varphi_j)_\Omega$  和  $\Phi_i = (\varphi_i,f)_\Omega$  。
   * 默认情况下， $f_h$
   * 的边界值不需要也没有强加，但是这个函数有一些可选的参数，允许强加零边界值，或者在第一步，以类似于上面的方式将
   * $f$
   * 的边界值投影到网格边界的有限元空间上，然后用这些值作为
   * $f_h$
   * 的强加边界值。这个函数的参数排序是这样的：如果你不想先投影到边界，你不需要给出第二个正交公式（类型为`正交<dim-1>`，用于计算矩阵和边界值投影的右手边），但如果你想这样做，你必须这样做。
   * 如果满足以下条件，则使用MatrixFree实现。
   *
   *
   *
   *
   *
   * -  @p enforce_zero_boundary  是假的。
   *
   *
   *
   *
   *
   * -  @p project_to_boundary_first  是假的。
   *
   *
   *
   *
   * - FiniteElement是由MatrixFree类支持的。
   *
   *
   * - FiniteElement的组件少于5个
   *
   * - FiniteElement的度数小于9。
   *
   *
   *
   * - dim==spacedim 在这种情况下，这个函数使用给定的正交公式执行数值正交，用于右手边的积分 $\Phi_i$ ，而QGauss(fe_degree+2)对象被用于质量算子。因此，你应该确保给定的正交公式对于创建右手边是足够精确的。    否则，只支持串行三角计算，质量矩阵使用 MatrixTools::create_mass_matrix. 组装，然后给定的正交规则用于矩阵和右手边。  因此，你应该确保给定的正交公式也足以用来创建质量矩阵。特别是，正交公式的度数必须足够高，以确保质量矩阵是可逆的。例如，如果你使用的是FE_Q(k)元素，那么矩阵项 $M_{ij}$ 的积分在每个变量中都是多项式度数 $2k$ ，你需要一个在每个坐标方向都有 $k+1$ 点的高斯正交公式来确保 $M$ 是可倒的。    更多信息请参见该命名空间的一般文档。    在1d中，边界正交公式的默认值是一个无效的对象，因为在1d中不会发生边界上的积分。      @param[in]  mapping 要使用的映射对象。    @param[in]  dof 描述投射到的有限元空间的DoFHandler，对应于  @p vec.   @param[in]  约束 在组装质量矩阵时使用的约束，通常在你有悬挂节点时需要。    @param[in]  正交 要用于组装质量矩阵的正交公式。    @param[in]  function 投射到有限元空间的函数。    @param[out]  vec 投射的函数将被存储在输出向量中。这个向量必须已经被初始化，并且不能有鬼魂元素。    @param[in]  enforce_zero_boundary 如果为真， @p vec  将有零边界条件。    @param[in]  q_boundary 如果 @p project_to_boundary_first 为真，将使用正交规则。    @param[in]  project_to_boundary_first 如果为真，在投射函数的内部之前执行对边界的投射。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  project(const Mapping<dim, spacedim> &                            mapping,
          const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const Quadrature<dim> &                                   quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec,
          const bool                 enforce_zero_boundary     = false,
          const Quadrature<dim - 1> &q_boundary                = (dim > 1 ?
                                                     QGauss<dim - 1>(2) :
                                                     Quadrature<dim - 1>(0)),
          const bool                 project_to_boundary_first = false);

  /**
   * 调用上面的project()函数，使用<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  project(const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const Quadrature<dim> &                                   quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec,
          const bool                 enforce_zero_boundary     = false,
          const Quadrature<dim - 1> &q_boundary                = (dim > 1 ?
                                                     QGauss<dim - 1>(2) :
                                                     Quadrature<dim - 1>(0)),
          const bool                 project_to_boundary_first = false);

  /**
   * 和上面一样，但有hp-capabilities。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  project(const hp::MappingCollection<dim, spacedim> &              mapping,
          const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const hp::QCollection<dim> &                              quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec,
          const bool                      enforce_zero_boundary = false,
          const hp::QCollection<dim - 1> &q_boundary = hp::QCollection<dim - 1>(
            dim > 1 ? QGauss<dim - 1>(2) : Quadrature<dim - 1>(0)),
          const bool project_to_boundary_first = false);

  /**
   * 调用上面的project()函数，使用一个 $Q_1$
   * 映射对象的集合，即使用
   * hp::StaticMappingQ1::mapping_collection. 。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  project(const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const hp::QCollection<dim> &                              quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec,
          const bool                      enforce_zero_boundary = false,
          const hp::QCollection<dim - 1> &q_boundary = hp::QCollection<dim - 1>(
            dim > 1 ? QGauss<dim - 1>(2) : Quadrature<dim - 1>(0)),
          const bool project_to_boundary_first = false);

  /**
   * 与上述标量值正交数据的投影相同。
   * 用户提供的函数应该根据单元格迭代器和正交数在正交点返回一个值，当然也应该与提供的
   * @p quadrature 对象一致，该对象将被用来组装右手边。
   * 这个函数可以和lambdas一起使用。
   * @code
   * VectorTools::project
   * (mapping,
   * dof_handler,
   * constraints,
   * quadrature_formula,
   * [&] (const typename DoFHandler<dim>::active_cell_iterator & cell,
   *     const unsigned int q)
   *
   * -> double
   * {
   *  return qp_data.get_data(cell)[q]->density;
   * },
   * field);
   * @endcode
   * 其中 <code>qp_data</code>
   * 是一个CellDataStorage对象，它存储正交点数据。
   *
   */
  template <int dim, typename VectorType, int spacedim>
  void
  project(const Mapping<dim, spacedim> &                            mapping,
          const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const Quadrature<dim> &                                   quadrature,
          const std::function<typename VectorType::value_type(
            const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
            const unsigned int)> &                                  func,
          VectorType &                                              vec_result);

  /**
   * 与上述标量值MatrixFree正交数据的投影相同。
   * 用户提供的函数 @p func
   * 应根据单元格号和正交点号在正交点返回一个VectorizedArray值，并应与
   * @p n_q_points_1d. 一致。这个函数可以与lambdas一起使用。
   * @code
   * VectorTools::project
   * (matrix_free_data,
   * constraints,
   * 3,
   * [&] (const unsigned int cell,
   *     const unsigned int q)
   *
   * -> VectorizedArray<double>
   * {
   *  return qp_data(cell,q);
   * },
   * field);
   * @endcode
   * 其中 <code>qp_data</code> 是一个类型为Table<2,
   * VectorizedArray<double>>的对象，它存储正交点数据。      @p
   * fe_component 允许额外指定使用 @p data
   * 的哪个组件，如果它是用 <code>std::vector<const
   * DoFHandler<dim>*></code>构造的。它将在内部用于FEEvaluation对象的构造器。
   *
   */
  template <int dim, typename VectorType>
  void
  project(
    std::shared_ptr<
      const MatrixFree<dim,
                       typename VectorType::value_type,
                       VectorizedArray<typename VectorType::value_type>>> data,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    const unsigned int                                        n_q_points_1d,
    const std::function<VectorizedArray<typename VectorType::value_type>(
      const unsigned int,
      const unsigned int)> &                                  func,
    VectorType &                                              vec_result,
    const unsigned int                                        fe_component = 0);

  /**
   * 与上述相同，但对于<code>n_q_points_1d =
   * matrix_free.get_dof_handler().get_fe().degree+1</code>。
   *
   */
  template <int dim, typename VectorType>
  void
  project(
    std::shared_ptr<
      const MatrixFree<dim,
                       typename VectorType::value_type,
                       VectorizedArray<typename VectorType::value_type>>> data,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    const std::function<VectorizedArray<typename VectorType::value_type>(
      const unsigned int,
      const unsigned int)> &                                  func,
    VectorType &                                              vec_result,
    const unsigned int                                        fe_component = 0);

  // @}

} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_project_h


