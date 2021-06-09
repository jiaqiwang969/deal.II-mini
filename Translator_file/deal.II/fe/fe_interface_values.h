//include/deal.II-translator/fe/fe_interface_values_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#ifndef dealii_fe_interface_values_h
#define dealii_fe_interface_values_h

#include <deal.II/base/config.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/q_collection.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
template <int dim, int spacedim>
class FEInterfaceValues;
#endif

/**
 * 使用提取器访问FEInterfaceValues得到的视图的命名空间。
 *
 *
 */
namespace FEInterfaceViews
{
  /**
   * 视图的基类。
   *
   */
  template <int dim, int spacedim = dim>
  class Base
  {
  public:
    /**
     * 构造函数。
     *
     */
    Base(const FEInterfaceValues<dim, spacedim> &fe_interface);

  protected:
    /**
     * 存储一个指向FEInterfaceValues实例的指针。
     *
     */
    const FEInterfaceValues<dim, spacedim> *fe_interface;
  };



  /**
   * FEInterfaceValues的标量变量的视图。
   *
   */
  template <int dim, int spacedim = dim>
  class Scalar : public Base<dim, spacedim>
  {
  public:
    /**
     * 这是返回数值的类型。
     *
     */
    using value_type = double;

    /**
     * 这是返回梯度的类型，例如从average_gradient()返回。
     *
     */
    using gradient_type =
      typename FEValuesViews::Scalar<dim, spacedim>::gradient_type;

    /**
     * 这是为hesians返回的类型，例如从jump_hessian()返回。
     *
     */
    using hessian_type =
      typename FEValuesViews::Scalar<dim, spacedim>::hessian_type;

    /**
     * 这是返回三阶导数的类型，例如从jump_hessian()返回。
     *
     */
    using third_derivative_type =
      typename FEValuesViews::Scalar<dim, spacedim>::third_derivative_type;

    /**
     * 表示一个单一标量分量的对象的构造函数
     *
     */
    Scalar(const FEInterfaceValues<dim, spacedim> &fe_interface,
           const unsigned int                      component);

    /**
     * 返回该视图所选分量的正交点 @p q_point
     * 中具有界面dof索引 @p interface_dof_index
     * 的形状函数的值。        参数 @p here_or_there
     * 在上游值和下游值之间进行选择，上游值由该正交点的法向量方向定义。如果
     * @p here_or_there
     * 为真，则使用界面的第一个单元的形状函数。
     * 换句话说，当从界面的两个单元之一接近给定的正交点时，该函数返回形状函数值的极限。
     * @note
     * 这个函数通常用于根据一个方向来挑选上游或下游的值。这可以通过使用
     * <code>(direction normal)>0</code>
     * 作为该函数的第一个参数来实现。
     *
     */
    value_type
    value(const bool         here_or_there,
          const unsigned int interface_dof_index,
          const unsigned int q_point) const;

    /**
     * 返回该视图所选部件的正交点 @p q_point 中形状函数 @p
     * interface_dof_index 在界面上的跳变 $\jump{u}=u_1
     *
     * - u_2$ 。
     *
     */
    value_type
    jump(const unsigned int interface_dof_index,
         const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point 的形状函数 @p
     * interface_dof_index 的界面上的平均值
     * $\average{u}=\frac{1}{2}(u_1 + u_2)$ 。
     *
     */
    value_type
    average(const unsigned int interface_dof_index,
            const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point 的形状函数 @p
     * interface_dof_index 在界面上的梯度 $\average{\nabla u}$
     * 的平均值。
     *
     */
    gradient_type
    average_gradient(const unsigned int interface_dof_index,
                     const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point 的形状函数 @p
     * interface_dof_index 在界面上的梯度 $\jump{nabla u}$ 的跳跃。
     *
     */
    gradient_type
    jump_gradient(const unsigned int interface_dof_index,
                  const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point处形状函数 @p
     * interface_dof_index 界面上的Hessian $\average{\nabla^2 u} =
     * \frac{1}{2}\nabla^2 u_{\text{cell0}} + \frac{1}{2} \nabla^2
     * u_{\text{cell1}}$ 的平均值。
     *
     */
    hessian_type
    average_hessian(const unsigned int interface_dof_index,
                    const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point 处的形状函数
     * @p interface_dof_index的界面上的梯度 $\jump{\nabla u}=\nabla
     * u_{\text{cell0}}
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * - \nabla u_{\text{cell1}}$ 的跳变。
     *
     */
    hessian_type
    jump_hessian(const unsigned int interface_dof_index,
                 const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point 处形状函数 @p
     * interface_dof_index 的界面上的第三导数 $\jump{\nabla^3 u} =
     * \nabla^3 u_{\text{cell0}}
     *
     * - \nabla^3 u_{\text{cell1}}$ 的跳变。
     *
     */
    third_derivative_type
    jump_3rd_derivative(const unsigned int interface_dof_index,
                        const unsigned int q_point) const;

  private:
    /**
     * 该视图的提取器。
     *
     */
    const FEValuesExtractors::Scalar extractor;
  };



  /**
   * FEInterfaceValues的矢量值变量的视图。
   *
   */
  template <int dim, int spacedim = dim>
  class Vector : public Base<dim, spacedim>
  {
  public:
    /**
     * 这是返回数值的类型。
     *
     */
    using value_type =
      typename FEValuesViews::Vector<dim, spacedim>::value_type;

    /**
     * 这是返回梯度的类型，例如从average_gradient()返回。
     *
     */
    using gradient_type =
      typename FEValuesViews::Vector<dim, spacedim>::gradient_type;

    /**
     * 这个类所代表的视图的二阶导数的类型的别名。这里，对于一组
     * <code>dim</code> 分量的有限元，Hessian是一个
     * <code>Tensor@<3,dim@></code>  。
     *
     */
    using hessian_type =
      typename FEValuesViews::Vector<dim, spacedim>::hessian_type;

    /**
     * 该类代表的视图的第三导数类型的别名。这里，对于有限元的一组
     * <code>dim</code> 分量，第三导数是一个
     * <code>Tensor@<4,dim@></code>  。
     *
     */
    using third_derivative_type =
      typename FEValuesViews::Vector<dim, spacedim>::third_derivative_type;

    /**
     * 代表一个矢量分量的对象的构造函数
     *
     */
    Vector(const FEInterfaceValues<dim, spacedim> &fe_interface,
           const unsigned int                      first_vector_component);

    /**
     * 返回该视图选择的具有界面dof索引 @p interface_dof_index
     * 的向量分量在正交点 @p q_point. 的值 参数 @p here_or_there
     * 在该正交点的法向量方向定义的上游值和下游值之间选择。如果
     * @p here_or_there
     * 为真，则使用界面的第一个单元的形状函数。
     * 换句话说，当从界面的两个单元之一接近给定正交点时，该函数返回形状函数的极限值。
     * @note
     * 这个函数通常用于根据一个方向来挑选上游或下游的值。这可以通过使用
     * <code>(direction normal)>0</code>
     * 作为该函数的第一个参数来实现。
     *
     */
    value_type
    value(const bool         here_or_there,
          const unsigned int interface_dof_index,
          const unsigned int q_point) const;

    /**
     * 返回正交点 @p q_point. 中形状函数 @p interface_dof_index
     * 的界面上的跳跃向量 $[\mathbf{u}]=\mathbf{u_1}
     *
     * - \mathbf{u_2}$ 。
     *
     */
    value_type
    jump(const unsigned int interface_dof_index,
         const unsigned int q_point) const;

    /**
     * 返回正交点 @p q_point. 的形状函数 @p interface_dof_index
     * 的界面上的平均矢量
     * $\average{\mathbf{u}}=\frac{1}{2}(\matbf{u_1} + \mathbf{u_2})$ 。
     *
     */
    value_type
    average(const unsigned int interface_dof_index,
            const unsigned int q_point) const;

    /**
     * 返回正交点 @p q_point. 的形状函数 @p interface_dof_index
     * 在界面上的梯度（等级2的张量） $\average{\nabla
     * \mathbf{u}}$ 的平均值。
     *
     */
    gradient_type
    average_gradient(const unsigned int interface_dof_index,
                     const unsigned int q_point) const;

    /**
     * 返回形状函数 @p interface_dof_index 在正交点 @p q_point.
     * 的界面上的梯度（等级为2的张量） $\jump{\nabla
     * \mathbf{u}}$ 的跳跃。
     *
     */
    gradient_type
    jump_gradient(const unsigned int interface_dof_index,
                  const unsigned int q_point) const;

    /**
     * 返回该视图所选组件的正交点 @p q_point处形状函数 @p
     * interface_dof_index 的界面上的Hessian $\average{\nabla^2 u} =
     * \frac{1}{2}\nabla^2 u_{\text{cell0}} + \frac{1}{2} \nabla^2
     * u_{\text{cell1}}$ 的平均值。
     *
     */
    hessian_type
    average_hessian(const unsigned int interface_dof_index,
                    const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point 处的形状函数
     * @p interface_dof_index的界面上的梯度 $\jump{\nabla u}=\nabla
     * u_{\text{cell0}}
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * - \nabla u_{\text{cell1}}$ 的跳变。
     *
     */
    hessian_type
    jump_hessian(const unsigned int interface_dof_index,
                 const unsigned int q_point) const;

    /**
     * 返回此视图所选组件的正交点 @p q_point 处形状函数 @p
     * interface_dof_index 的界面上的第三导数 $\jump{\nabla^3 u} =
     * \nabla^3 u_{\text{cell0}}
     *
     * - \nabla^3 u_{\text{cell1}}$ 的跳变。
     *
     */
    third_derivative_type
    jump_3rd_derivative(const unsigned int interface_dof_index,
                        const unsigned int q_point) const;

  private:
    /**
     * 此视图的提取器。
     *
     */
    const FEValuesExtractors::Vector extractor;
  };
} // namespace FEInterfaceViews



/**
 * FEInterfaceValues是一个数据结构，用于访问和组合网格中两个单元之间的界面上的有限元数据。
 * 它提供了一种方法来访问两个相邻单元之间的界面上的平均数、跳跃项和用于非连续加尔金方法的类似操作。这允许以类似于FEValues对单元、FEFaceValues对面的方式计算典型的网格相关的线性或双线性形式。在文献中，相邻单元之间的面被称为
 * "内部界面 "或 "面"。
 * 在内部，该类为两个FEFaceValues对象（或使用自适应细化时的FESubfaceValues）提供了一个抽象。该类引入了一个新的
 * "界面dof索引"，它在两个FEFaceValues对象的dof索引的联盟上行走。辅助函数允许在新的
 * "界面dof索引 "和相应的
 * "单元格索引"（0代表第一个单元格，1代表第二个单元格）以及该单元格中的
 * "dof索引 "之间进行转换。 该类是在 MeshWorker::mesh_loop().
 * 中使用的，其目的是作为MeshWorker和LocalIntegrators的低级替代品，与手动组装面状术语相比，它是一个更高级的抽象。
 *
 *
 */
template <int dim, int spacedim = dim>
class FEInterfaceValues
{
public:
  /**
   * 正交点的数量。
   *
   */
  const unsigned int n_quadrature_points;

  /**
   * 用一个有限元素（FiniteElement）构建FEInterfaceValues（面的两边都一样）。FEFaceValues对象将以给定的
   * @p mapping,   @p quadrature, 和 @p update_flags. 进行初始化。
   *
   */
  FEInterfaceValues(const Mapping<dim, spacedim> &      mapping,
                    const FiniteElement<dim, spacedim> &fe,
                    const Quadrature<dim - 1> &         quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * 与上述相同，但取一个正交规则的集合，以便不同的正交规则可以分配给不同的面。
   *
   */
  FEInterfaceValues(const Mapping<dim, spacedim> &      mapping,
                    const FiniteElement<dim, spacedim> &fe,
                    const hp::QCollection<dim - 1> &    quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * 用一个FiniteElement和一个Q1 Mapping构造FEInterfaceValues。
   * 见上面的构造函数。
   *
   */
  FEInterfaceValues(const FiniteElement<dim, spacedim> &fe,
                    const Quadrature<dim - 1> &         quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * 重新初始化这个对象，以用于由两个相邻单元的两个面给出的新接口。`cell`和`cell_neighbor`单元格将通过`cell_index`零和一在此调用后在所有需要识别界面相邻的两个单元格的地方被提及。
   * 使用  numbers::invalid_unsigned_int  表示  @p sub_face_no  或  @p
   * sub_face_no_neighbor 表示你想在整个面而不是子面工作。
   * 参数（包括它们的顺序）与 @p face_worker 中的
   * MeshWorker::mesh_loop(). 参数相同  @param[in]  cell
   * 与界面相邻的第一个单元的迭代器。    @param[in]  face_no
   * 一个整数，识别接口在第一个单元的哪个面上。
   * @param[in]  sub_face_no
   * 一个整数，标识界面所对应的面（由前面两个参数标识）的子面（子）。如果等于
   * numbers::invalid_unsigned_int, ，那么该界面被认为是整个面。
   * @param[in]  cell_neighbor
   * 一个与界面相邻的第二个单元格的迭代器。这个迭代器的类型不一定等于`cell`，但必须可以转换为`cell`。这允许对`cell`使用活动单元迭代器，对`cell_neighbor`使用`cell->neighbor(f)`，因为`cell->neighbor(f)`的返回类型只是一个单元迭代器（不一定是一个活动单元迭代器）。
   * @param[in]  face_no_neighbor
   * 和`face_no`一样，只是针对邻近的单元。    @param[in]
   * sub_face_no_neighbor
   * 和`sub_face_no`一样，只是针对相邻的单元。
   *
   */
  template <class CellIteratorType>
  void
  reinit(const CellIteratorType &                         cell,
         const unsigned int                               face_no,
         const unsigned int                               sub_face_no,
         const typename identity<CellIteratorType>::type &cell_neighbor,
         const unsigned int                               face_no_neighbor,
         const unsigned int                               sub_face_no_neighbor);

  /**
   * 重新初始化此对象，使其用于由单元格的单个面 @p
   * face_no 给出的界面 @p cell.
   * 这对于在域的边界上使用FEInterfaceValues很有用。
   * 因此，像jump()这样的成员将假定 "另一边
   * "的值为零。请注意，不需要sub_face_number，因为一个边界面不能与更细的单元相邻。
   * 调用此函数at_boundary()后将返回true。
   *
   */
  template <class CellIteratorType>
  void
  reinit(const CellIteratorType &cell, const unsigned int face_no);

  /**
   * 返回对界面中指定单元的FEFaceValues或FESubfaceValues对象的引用。
   * @p cell_index
   * 是0或1，对应于interface_dof_to_cell_and_dof_index()返回的单元格索引。
   *
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_face_values(const unsigned int cell_index) const;

  /**
   * 对所选映射对象的常量引用。
   *
   */
  const Mapping<dim, spacedim> &
  get_mapping() const;

  /**
   * 对所选有限元对象的常数引用。
   *
   */
  const FiniteElement<dim, spacedim> &
  get_fe() const;

  /**
   * 返回对使用中的正交对象的引用。
   *
   */
  const Quadrature<dim - 1> &
  get_quadrature() const;

  /**
   * 返回设置的更新标志。
   *
   */
  UpdateFlags
  get_update_flags() const;

  /**
   * @name  查询给定接口信息的函数  
     * @{ 
   *
   */

  /**
   * 返回当前界面是一个边界面还是一个有两个相邻单元的内部面。
   * 详见相应的reinit()函数。
   *
   */
  bool
  at_boundary() const;

  /**
   * 映射的正交权重。这个值等于映射的表面元素乘以正交点的权重。
   * 你可以把这个函数返回的数量看作是我们在这里通过四分法实现的积分中的表面元素
   * $ds$ 。      @dealiiRequiresUpdateFlags{update_JxW_values}
   *
   */
  double
  JxW(const unsigned int quadrature_point) const;

  /**
   * 返回每个正交点的JxW值的向量。
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   *
   */
  const std::vector<double> &
  get_JxW_values() const;

  /**
   * 返回每个正交点中界面的法向量。
   * 返回值与get_fe_face_values(0).get_normal_vectors()相同，因此，从这个界面的第一个单元的角度来看，是外部法向量。
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   *
   */
  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const;

  /**
   * 返回实空间中正交点的引用。
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   *
   */
  const std::vector<Point<spacedim>> &
  get_quadrature_points() const;

  /**
   * 返回当前界面上的DoFs（或形状函数）的数量。
   * @note
   * 这个数字只有在调用reinit()后才能得到，而且在调用reinit()后会发生变化。例如，在一个边界界面上，它等于单个FEFaceValues对象的道夫数，而对于一个DG元素的内部界面，它是这个数字的两倍。对于一个连续元素，它略小，因为界面上的两个单元共享一些道夫。
   *
   */
  unsigned
  n_current_interface_dofs() const;

  /**
   * 返回联合DoF指数的集合。这包括两个单元的指数。
   * 如果调用 reinit
   * 时有一个活动单元的迭代器，该指数基于活动指数（由
   * `DoFCellAccessor::get_dof_indices()`
   * 返回），如果是水平单元（即，如果 is_level_cell() 返回
   * true ），则返回 mg dof 指数。
   * @note
   * 这个函数只有在调用reinit()后才能使用，并且在调用reinit()后会发生变化。
   *
   */
  std::vector<types::global_dof_index>
  get_interface_dof_indices() const;

  /**
   * 将一个界面DoF指数转换成两个单元的相应本地DoF指数。如果一个界面DoF只在其中一个单元上起作用，另一个索引将是
   * numbers::invalid_unsigned_int.
   * 对于不连续的有限元，每个界面Dof将正好对应一个DoF索引。
   * @note
   * 该函数仅在调用reinit()后可用，并且可以从一次调用reinit()到下一次调用而改变。
   *
   */
  std::array<unsigned int, 2>
  interface_dof_to_dof_indices(const unsigned int interface_dof_index) const;

  /**
   * 返回给定正交点的法线。
   * 法线指向从这个界面的第一个单元看到的向外的方向。
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   *
   */
  Tensor<1, spacedim>
  normal(const unsigned int q_point_index) const;

  /**
   * @}
   *
   */

  /**
   * @name  用于评估形状函数数据的函数  @{
   *
   */

  /**
   * 返回正交点 @p q_point. 中具有接口道夫索引 @p
   * interface_dof_index 的形状函数值的分量 @p component  参数 @p
   * here_or_there 在0单元（这里为 @p true) ）和1单元（那里为
   * @p false). ）的值之间选择。你也可以把它解释为
   * "上游"（ @p true) 和 "下游"（ @p false)
   * 由这个正交点的法向量的方向定义。如果 @p here_or_there
   * 为真，则使用界面的第一个单元的形状函数。
   * 换句话说，当从界面的两个单元之一接近给定正交点时，该函数返回形状函数的极限值。
   * @note
   * 这个函数通常用于根据一个方向来挑选上游或下游的值。这可以通过使用
   * <code>(direction normal)>0</code>
   * 作为该函数的第一个参数来实现。
   *
   */
  double
  shape_value(const bool         here_or_there,
              const unsigned int interface_dof_index,
              const unsigned int q_point,
              const unsigned int component = 0) const;

  /**
   * 返回形状函数 @p interface_dof_index 在分量 @p component.
   * 的正交点 @p q_point 的界面上的跳变 $\jump{u}=u_{\text{cell0}}
   *
   * - u_{\text{cell1}}$
   * 注意，可以用不同的方式定义跳变（"那里 "的值减去
   * "这里
   * "的值，或者用其他方式；两者在有限元文献中都使用）。这里的定义使用
   * "这里的值减去那里的值"，从第一个单元格看。
   * 如果这是一个边界面（at_boundary()返回true），那么
   * $\jump{u}=u_{\text{cell0}}$  .
   *
   */
  double
  jump(const unsigned int interface_dof_index,
       const unsigned int q_point,
       const unsigned int component = 0) const;

  /**
   * 返回组件 @p component. 的正交点 @p q_point 处的形状函数
   * $\average{u}=\frac{1}{2}u_{\text{cell0}} + \frac{1}{2}u_{\text{cell1}}$
   * 在界面上的平均值
   * 如果这是一个边界面（at_boundary()返回true），那么
   * $\average{u}=u_{\text{cell0}}$  .
   *
   */
  double
  average(const unsigned int interface_dof_index,
          const unsigned int q_point,
          const unsigned int component = 0) const;

  /**
   * 返回组件 @p component. 的正交点 @p q_point处的形状函数
   * $\average{\nabla u} = \frac{1}{2}\nabla u_{\text{cell0}} + \frac{1}{2}
   * \nabla u_{\text{cell1}}$
   * 在界面上的梯度平均值，如果这是一个边界面（at_boundary()返回true），那么
   * $\average{\nabla u}=\nabla u_{\text{cell0}}$  .
   *
   */
  Tensor<1, spacedim>
  average_gradient(const unsigned int interface_dof_index,
                   const unsigned int q_point,
                   const unsigned int component = 0) const;

  /**
   * 返回形状函数 @p interface_dof_index 在组件 @p
   * q_point的正交点 @p component. 的界面上的Hessian
   * $\average{\nabla^2 u} = \frac{1}{2}\nabla^2 u_{\text{cell0}} +
   * \frac{1}{2} \nabla^2 u_{\text{cell1}}$ 的平均值
   * 如果这是一个边界面（at_boundary()返回true），那么
   * $\average{\nabla^2 u}=\nabla^2 u_{\text{cell0}}$  .
   *
   */
  Tensor<2, spacedim>
  average_hessian(const unsigned int interface_dof_index,
                  const unsigned int q_point,
                  const unsigned int component = 0) const;

  /**
   * 返回形状函数 @p interface_dof_index在组件 @p 的正交点 @p
   * q_point 的界面上的梯度跳变 $\jump{\nabla u}=\nabla
   * u_{\text{cell0}}
   *
   * - \nabla u_{\text{cell1}}$ 。
   * 如果这是一个边界面（at_boundary()返回true），那么
   * $\jump{\nabla u}=\nabla u_{\text{cell0}}$  .
   *
   */
  Tensor<1, spacedim>
  jump_gradient(const unsigned int interface_dof_index,
                const unsigned int q_point,
                const unsigned int component = 0) const;

  /**
   * 返回形状函数 @p interface_dof_index 在组件 @p component.
   * 的正交点 @p q_point 的界面上的Hessian $\jump{\nabla^2 u} =
   * \nabla^2 u_{\text{cell0}}
   *
   * - \nabla^2 u_{\text{cell1}}$ 的跳变
   * 如果这是一个边界面（at_boundary()返回true），那么
   * $\jump{\nabla^2 u} = \nabla^2 u_{\text{cell0}}$  .
   *
   */
  Tensor<2, spacedim>
  jump_hessian(const unsigned int interface_dof_index,
               const unsigned int q_point,
               const unsigned int component = 0) const;

  /**
   * 返回形状函数 @p interface_dof_index 在组件 @p component.
   * 的正交点 @p q_point 的界面上的第三导数 $\jump{\nabla^3 u} =
   * \nabla^3 u_{\text{cell0}}
   *
   * - \nabla^3 u_{\text{cell1}}$
   * 的跳变，如果这是一个边界面（at_boundary()返回真），那么
   * $\jump{\nabla^3 u} = \nabla^3 u_{\text{cell0}}$  .
   *
   */
  Tensor<3, spacedim>
  jump_3rd_derivative(const unsigned int interface_dof_index,
                      const unsigned int q_point,
                      const unsigned int component = 0) const;

  /**
   * 创建一个当前FEInterfaceValues对象的视图，该视图代表了可能是矢量值的有限元中的一个特定标量分量。
   * 视图的概念在命名空间FEValuesViews的文档中解释。
   *
   */
  const FEInterfaceViews::Scalar<dim, spacedim>
  operator[](const FEValuesExtractors::Scalar &scalar) const;

  /**
   * 为当前的FEInterfaceValues对象创建一个视图，该视图表示矢量值有限元的一组
   * <code>dim</code>
   * 标量分量（即一个矢量）。视图的概念在命名空间FEValuesViews的文档中有所解释。
   *
   */
  const FEInterfaceViews::Vector<dim, spacedim>
  operator[](const FEValuesExtractors::Vector &vector) const;

  /**
   * @}
   *
   */

private:
  /**
   * 当前界面的DoF指数列表，在reinit()中填写。
   *
   */
  std::vector<types::global_dof_index> interface_dof_indices;

  /**
   * 从界面DoF到FeFaceValues对象的两个局部DoF指数的映射。如果一个界面DoF只在其中一个单元上活动，另一个单元将有
   * numbers::invalid_unsigned_int. 。
   *
   */
  std::vector<std::array<unsigned int, 2>> dofmap;

  /**
   * 当前单元的FEFaceValues对象。
   *
   */
  FEFaceValues<dim, spacedim> internal_fe_face_values;

  /**
   * 当前单元格的FEFaceValues对象，如果该单元格被精炼。
   *
   */
  FESubfaceValues<dim, spacedim> internal_fe_subface_values;

  /**
   * 邻近单元格的FEFaceValues对象。
   *
   */
  FEFaceValues<dim, spacedim> internal_fe_face_values_neighbor;

  /**
   * 如果该单元格被细化，则为相邻单元格的FEFaceValues对象。
   *
   */
  FESubfaceValues<dim, spacedim> internal_fe_subface_values_neighbor;

  /**
   * 指向internal_fe_face_values或internal_fe_subface_values的指针，分别在reinit()中确定。
   *
   */
  FEFaceValuesBase<dim, spacedim> *fe_face_values;

  /**
   * 指向internal_fe_face_values_neighbor、internal_fe_subface_values_neighbor或nullptr的指针，分别在reinit()中确定。
   *
   */
  FEFaceValuesBase<dim, spacedim> *fe_face_values_neighbor;

  /* 让视图类成为该类的朋友，因为它们访问内部数据。 
*
*/
  template <int, int>
  friend class FEInterfaceViews::Scalar;
  template <int, int>
  friend class FEInterfaceViews::Vector;
};



#ifndef DOXYGEN

 /*---------------------- Inline functions ---------------------*/ 

template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const Mapping<dim, spacedim> &      mapping,
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1> &         quadrature,
  const UpdateFlags                   update_flags)
  : n_quadrature_points(quadrature.size())
  , internal_fe_face_values(mapping, fe, quadrature, update_flags)
  , internal_fe_subface_values(mapping, fe, quadrature, update_flags)
  , internal_fe_face_values_neighbor(mapping, fe, quadrature, update_flags)
  , internal_fe_subface_values_neighbor(mapping, fe, quadrature, update_flags)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const Mapping<dim, spacedim> &      mapping,
  const FiniteElement<dim, spacedim> &fe,
  const hp::QCollection<dim - 1> &    quadrature,
  const UpdateFlags                   update_flags)
  : n_quadrature_points(quadrature.max_n_quadrature_points())
  , internal_fe_face_values(mapping, fe, quadrature, update_flags)
  , internal_fe_subface_values(mapping, fe, quadrature, update_flags)
  , internal_fe_face_values_neighbor(mapping, fe, quadrature[0], update_flags)
  , internal_fe_subface_values_neighbor(mapping,
                                        fe,
                                        quadrature[0],
                                        update_flags)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1> &         quadrature,
  const UpdateFlags                   update_flags)
  : n_quadrature_points(quadrature.size())
  , internal_fe_face_values(
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe,
      quadrature,
      update_flags)
  , internal_fe_subface_values(
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe,
      quadrature,
      update_flags)
  , internal_fe_face_values_neighbor(
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe,
      quadrature,
      update_flags)
  , internal_fe_subface_values_neighbor(
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe,
      quadrature,
      update_flags)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{}



template <int dim, int spacedim>
template <class CellIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(
  const CellIteratorType &                         cell,
  const unsigned int                               face_no,
  const unsigned int                               sub_face_no,
  const typename identity<CellIteratorType>::type &cell_neighbor,
  const unsigned int                               face_no_neighbor,
  const unsigned int                               sub_face_no_neighbor)
{
  if (sub_face_no == numbers::invalid_unsigned_int)
    {
      internal_fe_face_values.reinit(cell, face_no);
      fe_face_values = &internal_fe_face_values;
    }
  else
    {
      internal_fe_subface_values.reinit(cell, face_no, sub_face_no);
      fe_face_values = &internal_fe_subface_values;
    }
  if (sub_face_no_neighbor == numbers::invalid_unsigned_int)
    {
      internal_fe_face_values_neighbor.reinit(cell_neighbor, face_no_neighbor);
      fe_face_values_neighbor = &internal_fe_face_values_neighbor;
    }
  else
    {
      internal_fe_subface_values_neighbor.reinit(cell_neighbor,
                                                 face_no_neighbor,
                                                 sub_face_no_neighbor);
      fe_face_values_neighbor = &internal_fe_subface_values_neighbor;
    }

  AssertDimension(fe_face_values->n_quadrature_points,
                  fe_face_values_neighbor->n_quadrature_points);

  const_cast<unsigned int &>(this->n_quadrature_points) =
    fe_face_values->n_quadrature_points;

  // Set up dof mapping and remove duplicates (for continuous elements).
  {
    // Get dof indices first:
    std::vector<types::global_dof_index> v(
      fe_face_values->get_fe().n_dofs_per_cell());
    cell->get_active_or_mg_dof_indices(v);
    std::vector<types::global_dof_index> v2(
      fe_face_values_neighbor->get_fe().n_dofs_per_cell());
    cell_neighbor->get_active_or_mg_dof_indices(v2);

    // Fill a map from the global dof index to the left and right
    // local index.
    std::map<types::global_dof_index, std::pair<unsigned int, unsigned int>>
                                          tempmap;
    std::pair<unsigned int, unsigned int> invalid_entry(
      numbers::invalid_unsigned_int, numbers::invalid_unsigned_int);

    for (unsigned int i = 0; i < v.size(); ++i)
      {
        // If not already existing, add an invalid entry:
        auto result = tempmap.insert(std::make_pair(v[i], invalid_entry));
        result.first->second.first = i;
      }

    for (unsigned int i = 0; i < v2.size(); ++i)
      {
        // If not already existing, add an invalid entry:
        auto result = tempmap.insert(std::make_pair(v2[i], invalid_entry));
        result.first->second.second = i;
      }

    // Transfer from the map to the sorted std::vectors.
    dofmap.resize(tempmap.size());
    interface_dof_indices.resize(tempmap.size());
    unsigned int idx = 0;
    for (auto &x : tempmap)
      {
        interface_dof_indices[idx] = x.first;
        dofmap[idx]                = {{x.second.first, x.second.second}};
        ++idx;
      }
  }
}



template <int dim, int spacedim>
template <class CellIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(const CellIteratorType &cell,
                                         const unsigned int      face_no)
{
  internal_fe_face_values.reinit(cell, face_no);
  fe_face_values          = &internal_fe_face_values;
  fe_face_values_neighbor = nullptr;

  interface_dof_indices.resize(fe_face_values->get_fe().n_dofs_per_cell());
  cell->get_active_or_mg_dof_indices(interface_dof_indices);

  dofmap.resize(interface_dof_indices.size());

  for (unsigned int i = 0; i < interface_dof_indices.size(); ++i)
    {
      dofmap[i] = {{i, numbers::invalid_unsigned_int}};
    }
}



template <int dim, int spacedim>
inline double
FEInterfaceValues<dim, spacedim>::JxW(const unsigned int q) const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->JxW(q);
}



template <int dim, int spacedim>
const std::vector<double> &
FEInterfaceValues<dim, spacedim>::get_JxW_values() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_JxW_values();
}



template <int dim, int spacedim>
const std::vector<Tensor<1, spacedim>> &
FEInterfaceValues<dim, spacedim>::get_normal_vectors() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_normal_vectors();
}



template <int dim, int spacedim>
const Mapping<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_mapping() const
{
  return internal_fe_face_values.get_mapping();
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe() const
{
  return internal_fe_face_values.get_fe();
}



template <int dim, int spacedim>
const Quadrature<dim - 1> &
FEInterfaceValues<dim, spacedim>::get_quadrature() const
{
  return internal_fe_face_values.get_quadrature();
}



template <int dim, int spacedim>
const std::vector<Point<spacedim>> &
FEInterfaceValues<dim, spacedim>::get_quadrature_points() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_quadrature_points();
}



template <int dim, int spacedim>
UpdateFlags
FEInterfaceValues<dim, spacedim>::get_update_flags() const
{
  return internal_fe_face_values.get_update_flags();
}



template <int dim, int spacedim>
unsigned
FEInterfaceValues<dim, spacedim>::n_current_interface_dofs() const
{
  Assert(
    interface_dof_indices.size() > 0,
    ExcMessage(
      "n_current_interface_dofs() is only available after a call to reinit()."));
  return interface_dof_indices.size();
}



template <int dim, int spacedim>
bool
FEInterfaceValues<dim, spacedim>::at_boundary() const
{
  return fe_face_values_neighbor == nullptr;
}



template <int dim, int spacedim>
std::vector<types::global_dof_index>
FEInterfaceValues<dim, spacedim>::get_interface_dof_indices() const
{
  return interface_dof_indices;
}



template <int dim, int spacedim>
std::array<unsigned int, 2>
FEInterfaceValues<dim, spacedim>::interface_dof_to_dof_indices(
  const unsigned int interface_dof_index) const
{
  AssertIndexRange(interface_dof_index, n_current_interface_dofs());
  return dofmap[interface_dof_index];
}



template <int dim, int spacedim>
const FEFaceValuesBase<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe_face_values(
  const unsigned int cell_index) const
{
  AssertIndexRange(cell_index, 2);
  Assert(
    cell_index == 0 || !at_boundary(),
    ExcMessage(
      "You are on a boundary, so you can only ask for the first FEFaceValues object."));

  return (cell_index == 0) ? *fe_face_values : *fe_face_values_neighbor;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::normal(const unsigned int q_point_index) const
{
  return fe_face_values->normal_vector(q_point_index);
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::shape_value(
  const bool         here_or_there,
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
    return get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                       q_point,
                                                       component);
  if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
    return get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                       q_point,
                                                       component);

  return 0.0;
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::jump(const unsigned int interface_dof_index,
                                       const unsigned int q_point,
                                       const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  double value = 0.0;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                         q_point,
                                                         component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                         q_point,
                                                         component);
  return value;
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::average(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                       q_point,
                                                       component);

  double value = 0.0;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                               q_point,
                                                               component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                               q_point,
                                                               component);

  return value;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::average_gradient(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                      q_point,
                                                      component);

  Tensor<1, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                              q_point,
                                                              component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_grad_component(dof_pair[1],
                                                              q_point,
                                                              component);

  return value;
}



template <int dim, int spacedim>
Tensor<2, spacedim>
FEInterfaceValues<dim, spacedim>::average_hessian(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                         q_point,
                                                         component);

  Tensor<2, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                                 q_point,
                                                                 component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_hessian_component(dof_pair[1],
                                                                 q_point,
                                                                 component);

  return value;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::jump_gradient(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                      q_point,
                                                      component);

  Tensor<1, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                        q_point,
                                                        component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_grad_component(dof_pair[1],
                                                        q_point,
                                                        component);

  return value;
}



template <int dim, int spacedim>
Tensor<2, spacedim>
FEInterfaceValues<dim, spacedim>::jump_hessian(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                         q_point,
                                                         component);

  Tensor<2, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                           q_point,
                                                           component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_hessian_component(dof_pair[1],
                                                           q_point,
                                                           component);

  return value;
}



template <int dim, int spacedim>
Tensor<3, spacedim>
FEInterfaceValues<dim, spacedim>::jump_3rd_derivative(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_3rd_derivative_component(dof_pair[0],
                                                                q_point,
                                                                component);

  Tensor<3, spacedim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_3rd_derivative_component(dof_pair[0],
                                                                  q_point,
                                                                  component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_3rd_derivative_component(dof_pair[1],
                                                                  q_point,
                                                                  component);

  return value;
}



 /*------------ Inline functions: FEInterfaceValues------------*/ 
template <int dim, int spacedim>
inline const FEInterfaceViews::Scalar<dim, spacedim>
  FEInterfaceValues<dim, spacedim>::
  operator[](const FEValuesExtractors::Scalar &scalar) const
{
  AssertIndexRange(scalar.component, this->get_fe().n_components());
  return FEInterfaceViews::Scalar<dim, spacedim>(*this, scalar.component);
}



template <int dim, int spacedim>
inline const FEInterfaceViews::Vector<dim, spacedim>
  FEInterfaceValues<dim, spacedim>::
  operator[](const FEValuesExtractors::Vector &vector) const
{
  const FiniteElement<dim, spacedim> &fe = this->get_fe();
  const unsigned int                  n_vectors =
    (fe.n_components() >= Tensor<1, spacedim>::n_independent_components ?
       fe.n_components() - Tensor<1, spacedim>::n_independent_components + 1 :
       0);
  (void)n_vectors;
  AssertIndexRange(vector.first_vector_component, n_vectors);
  return FEInterfaceViews::Vector<dim, spacedim>(*this,
                                                 vector.first_vector_component);
}



namespace FEInterfaceViews
{
  template <int dim, int spacedim>
  Base<dim, spacedim>::Base(
    const FEInterfaceValues<dim, spacedim> &fe_interface)
    : fe_interface(&fe_interface)
  {}



  template <int dim, int spacedim>
  Scalar<dim, spacedim>::Scalar(
    const FEInterfaceValues<dim, spacedim> &fe_interface,
    const unsigned int                      component)
    : Base<dim, spacedim>(fe_interface)
    , extractor(component)
  {}



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::value(const bool         here_or_there,
                               const unsigned int interface_dof_index,
                               const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
        dof_pair[1], q_point);

    return 0.0;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::jump(const unsigned int interface_dof_index,
                              const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    value_type value = 0.0;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::value_type
  Scalar<dim, spacedim>::average(const unsigned int interface_dof_index,
                                 const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    value_type value = 0.0;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value +=
        0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
                dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::gradient_type
  Scalar<dim, spacedim>::average_gradient(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .gradient(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::gradient_type
  Scalar<dim, spacedim>::jump_gradient(const unsigned int interface_dof_index,
                                       const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].gradient(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::hessian_type
  Scalar<dim, spacedim>::average_hessian(const unsigned int interface_dof_index,
                                         const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .hessian(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::third_derivative_type
  Scalar<dim, spacedim>::jump_3rd_derivative(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor]
        .third_derivative(dof_pair[0], q_point);

    third_derivative_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].third_derivative(
          dof_pair[0], q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -= (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                 .third_derivative(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Scalar<dim, spacedim>::hessian_type
  Scalar<dim, spacedim>::jump_hessian(const unsigned int interface_dof_index,
                                      const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].hessian(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  Vector<dim, spacedim>::Vector(
    const FEInterfaceValues<dim, spacedim> &fe_interface,
    const unsigned int                      first_vector_component)
    : Base<dim, spacedim>(fe_interface)
    , extractor(first_vector_component)
  {}



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::value(const bool         here_or_there,
                               const unsigned int interface_dof_index,
                               const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
      return (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
        dof_pair[1], q_point);

    return value_type();
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::jump(const unsigned int interface_dof_index,
                              const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    value_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::value_type
  Vector<dim, spacedim>::average(const unsigned int interface_dof_index,
                                 const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].value(
        dof_pair[0], q_point);

    value_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].value(dof_pair[0],
                                                                 q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value +=
        0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor].value(
                dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::gradient_type
  Vector<dim, spacedim>::average_gradient(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .gradient(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::gradient_type
  Vector<dim, spacedim>::jump_gradient(const unsigned int interface_dof_index,
                                       const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].gradient(
        dof_pair[0], q_point);

    gradient_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].gradient(dof_pair[0],
                                                                    q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].gradient(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::average_hessian(const unsigned int interface_dof_index,
                                         const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        0.5 *
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value += 0.5 * (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                       .hessian(dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::hessian_type
  Vector<dim, spacedim>::jump_hessian(const unsigned int interface_dof_index,
                                      const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor].hessian(
        dof_pair[0], q_point);

    hessian_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].hessian(dof_pair[0],
                                                                   q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -=
        (*(this->fe_interface->fe_face_values_neighbor))[extractor].hessian(
          dof_pair[1], q_point);

    return value;
  }



  template <int dim, int spacedim>
  typename Vector<dim, spacedim>::third_derivative_type
  Vector<dim, spacedim>::jump_3rd_derivative(
    const unsigned int interface_dof_index,
    const unsigned int q_point) const
  {
    const auto dof_pair = this->fe_interface->dofmap[interface_dof_index];

    if (this->fe_interface->at_boundary())
      return (*(this->fe_interface->fe_face_values))[extractor]
        .third_derivative(dof_pair[0], q_point);

    third_derivative_type value;

    if (dof_pair[0] != numbers::invalid_unsigned_int)
      value +=
        (*(this->fe_interface->fe_face_values))[extractor].third_derivative(
          dof_pair[0], q_point);

    if (dof_pair[1] != numbers::invalid_unsigned_int)
      value -= (*(this->fe_interface->fe_face_values_neighbor))[extractor]
                 .third_derivative(dof_pair[1], q_point);

    return value;
  }
} // namespace FEInterfaceViews

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


