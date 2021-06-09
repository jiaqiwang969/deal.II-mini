//include/deal.II-translator/fe/mapping_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_h
#define dealii_mapping_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>

#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/q_collection.h>

#include <array>
#include <cmath>
#include <memory>

DEAL_II_NAMESPACE_OPEN

template <typename ElementType, typename MemorySpaceType>
class ArrayView;
template <int dim>
class Quadrature;
template <int dim, int spacedim>
class FEValues;
template <int dim, int spacedim>
class FEValuesBase;
template <int dim, int spacedim>
class FEValues;
template <int dim, int spacedim>
class FEFaceValues;
template <int dim, int spacedim>
class FESubfaceValues;


/**
 * 用于 Mapping::transform() 函数的变换种类。
 * 特殊的有限元可能需要从参考单元到实际网格单元的特殊Mapping。为了做到最灵活，这个枚举为任意变换提供了一个可扩展的接口。尽管如此，这些必须在继承类的transform()函数中实现，才能发挥作用。
 *
 *
 * @ingroup mapping
 *
 *
 */
enum MappingKind
{
  /**
   * 无映射，即形状函数不从参考单元映射，而是直接在实空间单元上定义。
   *
   */
  mapping_none = 0x0000,

  /**
   * 协变映射（详见 Mapping::transform() ）。
   *
   */
  mapping_covariant = 0x0001,

  /**
   * 等变量映射（详见 Mapping::transform() ）。
   *
   */
  mapping_contravariant = 0x0002,

  /**
   * 协变向量场梯度的映射（详见 Mapping::transform() ）。
   *
   */
  mapping_covariant_gradient = 0x0003,

  /**
   * 不变向量场的梯度映射（详见 Mapping::transform() ）。
   *
   */
  mapping_contravariant_gradient = 0x0004,

  /**
   * 通常用于Hdiv元素的皮奥拉变换。Piola变换是H<sup>div</sup>中矢量值元素的标准变换。
   * 它相当于一个以体积元素的逆值为尺度的禁忌变换。
   *
   */
  mapping_piola = 0x0100,

  /**
   * 对应于mapping_piola变换的矢量场的梯度变换（详见
   * Mapping::transform() ）。
   *
   */
  mapping_piola_gradient = 0x0101,

  /**
   * 用于Nedelec元素的映射。
   * 卷曲元素被映射为协变向量。尽管如此，我们还是引入了一个单独的映射种类，这样我们就可以对向量和其梯度使用相同的标志（详见
   * Mapping::transform() ）。
   *
   */
  mapping_nedelec = 0x0200,

  /**
   * 用于Raviart-Thomas元素的映射。
   *
   */
  mapping_raviart_thomas = 0x0300,

  /**
   * 用于BDM元素的映射。
   *
   */
  mapping_bdm = mapping_raviart_thomas,

  /**
   * 2-forms和三阶张量的映射。
   * 这些是典型的应用于转换到参考单元的豫章的映射。
   * 协变向量场的赫斯的映射（详见 Mapping::transform() ）。
   *
   */
  mapping_covariant_hessian,

  /**
   * 忌变向量场的犹豫值的映射（详见 Mapping::transform() ）。
   *
   */
  mapping_contravariant_hessian,

  /**
   * 皮奥拉（piola）矢量场的赫斯的映射（详见
   * Mapping::transform() ）。
   *
   */
  mapping_piola_hessian
};


/**
 * @short  映射类的抽象基类。
 * 该类声明了用于描述从参考（单位）单元到实空间单元的映射功能的接口，以及用于填写使用FEValues、FEFaceValues和FESubfaceValues类所需的信息。这些接口的具体实现是在派生类中提供的。
 * <h3>Mathematics of the mapping</h3> 映射是一种转换 $\mathbf x =
 * \mathbf F_K(\hat{\mathbf  x})$ ，它将参考单元 $[0,1]^\text{dim}$
 * 中的点 $\hat{\mathbf x}$ 映射到实际网格单元 $K\subset{\mathbb
 * R}^\text{spacedim}$ 中的点 $\mathbf x$
 * 。这种映射的许多应用都需要这种映射的雅各布，
 * $J(\hat{\mathbf x}) =
 * \hat\nabla {\mathbf F}_K(\hat{\mathbf  x})$  。例如，如果dim=spacedim=2，我们有 @f[
 * J(\hat{\mathbf  x}) = \left(\begin{matrix}
 * \frac{\partial x}{\partial \hat x} & \frac{\partial x}{\partial \hat y}
 * \\
 * \frac{\partial y}{\partial \hat x} & \frac{\partial y}{\partial \hat y}
 * \end{matrix}\right)
 * @f] 。
 * <h4>%Mapping of scalar functions</h4>
 * *标量有限元的形状函数通常是在参考单元上定义的，然后根据规则@f[
 * \varphi(\mathbf x) = \varphi\bigl(\mathbf F_K(\hat{\mathbf  x})\bigr)
 * = \hat \varphi(\hat{\mathbf  x}).
 * @f]简单地进行映射。
 *
 *  <h4>%Mapping of integrals</h4>
 * 使用简单的变量变化，标量函数在一个单元上的积分  $K$  可以表示为在参考单元上的积分  $\hat K$  。具体来说，体积形式 $d\hat x$ 被转换为@f[
 * \int_K u(\mathbf x)\,dx = \int_{\hat K} \hat
 * u(\hat{\mathbf  x}) \left|\text{det}J(\hat{\mathbf  x})\right|
 * \,d\hat x.
 * @f]。
 * 在这种积分被正交近似的表达中，这就导致了形式为@f[
 * \int_K u(\mathbf x)\,dx
 * \approx
 * \sum_{q}
 * \hat u(\hat{\mathbf  x}_q)
 * \underbrace{\left|\text{det}J(\hat{\mathbf  x}_q)\right| w_q}_{=:
 * \text{JxW}_q}.
 * @f]的项。 这里，每个正交点的权重 $\text{JxW}_q$ （其中<i>JxW</i>象征着<i>Jacobian times Quadrature Weights</i>）在原始积分中扮演了 $dx$ 的角色。因此，它们出现在所有计算正交积分的代码中，并通过 FEValues::JxW(). 访问。
 * @todo  记录了在二维-1的情况下会发生什么。
 *
 *  <h4>%Mapping of vector fields, differential forms and gradients of vector
 * fields</h4> 矢量场或微分形式（标量函数的梯度） $\mathbf v$
 * ，以及矢量场的梯度 $\mathbf T$ 的变换遵循一般形式
 * @f[
 * \mathbf v(\mathbf x) = \mathbf A(\hat{\mathbf  x})
 * \hat{\mathbf  v}(\hat{\mathbf  x}),
 * \qquad
 * \mathbf T(\mathbf x) = \mathbf A(\hat{\mathbf  x})
 * \hat{\mathbf  T}(\hat{\mathbf  x}) \mathbf B(\hat{\mathbf  x}).
 * @f] 微分形式<b>A</b>和<b>B</b>是由被转换的对象的种类决定。这些转换是通过transform()函数进行的，被转换的对象的类型由其MappingKind参数指定。关于可能的选择，请看那里的文档。
 * <h4>Derivatives of the mapping</h4>
 * 一些应用需要映射的导数，其中一阶导数是映射的Jacobian，
 * $J_{iJ}(\hat{\mathbf x})=\frac{\partial x_i}{\partial \hat x_J}$
 * ，如上所述。映射的高阶导数也有类似的定义，例如，雅各布导数
 * $\hat H_{iJK}(\hat{\mathbf  x}) = \frac{\partial^2 x_i}{\partial \hat x_J
 * \partial \hat x_K}$  ，以及雅各布二阶导数  $\hat
 * K_{iJKL}(\hat{\mathbf  x}) = \frac{\partial^3 x_i}{\partial \hat x_J
 * \partial \hat x_K \partial \hat x_L}$  。定义高阶导数的 "前推
 * "版本也很有用：雅各布前推导数， $H_{ijk}(\hat{\mathbf x}) =
 * \frac{\partial^2 x_i}{\partial \hat x_J \partial \hat
 * x_K}(J_{jJ})^{-1}(J_{kK})^{-1}$  ，以及雅各布后推导数，
 * $K_{ijkl}(\hat{\mathbf  x}) = \frac{\partial^3 x_i}{\partial \hat x_J
 * \partial \hat x_K \partial \hat
 * x_L}(J_{jJ})^{-1}(J_{kK})^{-1}(J_{lL})^{-1}$
 * 。这些向前推的版本可以用来计算定义在参考单元上的函数相对于实际单元坐标的高阶导数。例如，相对于实细胞坐标的雅各布导数由以下公式给出。
 * @f[
 * \frac{\partial}{\partial x_j}\left[J_{iJ}(\hat{\mathbf  x})\right] =
 * H_{ikn}(\hat{\mathbf  x})J_{nJ}(\hat{\mathbf  x}),
 * @f]，而相对于实细胞坐标的雅各布反导也同样由以下公式给出。@f[
 * \frac{\partial}{\partial x_j}\left[\left(J_{iJ}(\hat{\mathbf
 * x})\right)^{-1}\right] =
 *
 * -H_{nik}(\hat{\mathbf  x})\left(J_{nJ}(\hat{\mathbf x})\right)^{-1}. @f]
 * 以类似的方式，在参考单元上定义的函数的高阶导数，相对于实数单元坐标，可以使用雅各布式推前高阶导数来定义。例如，Jacobian
 * pushed-forward导数相对于实际单元坐标的导数由以下公式给出。
 * @f[
 * \frac{\partial}{\partial x_l}\left[H_{ijk}(\hat{\mathbf  x})\right] =
 * K_{ijkl}(\hat{\mathbf  x})
 *
 * -H_{mjl}(\hat{\mathbf  x})H_{imk}(\hat{\mathbf
 * x})-H_{mkl}(\hat{\mathbf  x})H_{imj}(\hat{\mathbf  x}).
 * @f]
 * <h3>References</h3>
 * 关于微分几何和有限元的一般出版物是调查报告 <ul>   <li>  Douglas N. Arnold, Richard S. Falk, and Ragnar Winther. <i>Finite
 * element exterior calculus: from Hodge theory to numerical stability.</i>
 * Bull. Amer. Math. Soc. (N.S.), 47:281-354, 2010. <a
 * href="http://dx.doi.org/10.1090/S0273-0979-10-01278-4">DOI:
 * 10.1090/S0273-0979-10-01278-4</a>.   </ul>
 * 皮奥拉变换的描述来自休斯顿大学Ronald H. W. Hoppe的<a
 * href="http://www.math.uh.edu/~rohop/spring_05/downloads/">lecture
 * notes</a>，第七章。
 *
 *
 * @ingroup mapping
 *
 *
 */
template <int dim, int spacedim = dim>
class Mapping : public Subscriptor
{
public:
  /**
   * 虚拟解构器。
   *
   */
  virtual ~Mapping() override = default;

  /**
   * 返回一个指向当前对象副本的指针。然后，这个副本的调用者将拥有它的所有权。
   * 这个函数在这个基类中被声明为抽象的虚函数，派生类将不得不实现它。
   * 这个函数主要由 hp::MappingCollection 类使用。
   *
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const = 0;

  /**
   * 返回一个单元格的映射顶点。
   * 大多数时候，这些值将仅仅是由 <code>cell-@>vertex(v)</code>
   * 返回的顶点 <code>v</code>
   * 的坐标，即由三角法存储的信息。
   * 然而，也有增加位移或选择完全不同位置的映射，例如MappingQEulerian,
   * MappingQ1Eulerian, 或MappingFEField。
   * 这个函数的默认实现只是返回三角形所存储的信息，即：
   * <code>cell-@>vertex(v)</code>  .
   *
   */
  virtual boost::container::small_vector<Point<spacedim>,
                                         GeometryInfo<dim>::vertices_per_cell>
  get_vertices(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

  /**
   * 返回一个单元格的映射中心。
   * 如果你使用的是保留顶点位置的(bi-,tri-)线性映射，这个函数只是返回同样由`cell->center()`产生的值。然而，也有一些映射会增加位移或选择完全不同的位置，例如MappingQEulerian、MappingQ1Eulerian或MappingFEField，以及基于高阶多项式的映射，对于这些映射，中心可能不会与顶点位置的平均值重合。
   * 默认情况下，该函数返回参考单元中心的前推。如果参数
   * @p map_center_of_reference_cell
   * 被设置为false，那么返回值将是由get_vertices()方法返回的顶点位置的平均值。
   * @param[in]  cell 你想计算中心的单元格  @param[in]
   * map_center_of_reference_cell
   * 一个标志，用于将计算单元格中心的算法从应用于参考单元格中心的transform_unit_to_real_cell()转换为计算顶点平均数。
   *
   */
  virtual Point<spacedim>
  get_center(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
             const bool map_center_of_reference_cell = true) const;

  /**
   * 返回映射的单元格的边界盒。
   * 如果你使用的是保留顶点位置的(bi-,tri-)线性映射，这个函数简单地返回同样由`cell->bounding_box()`产生的值。然而，也有一些映射会增加位移或选择完全不同的位置，例如MappingQEulerian、MappingQ1Eulerian或MappingFEField。
   * 对于线性映射，该函数返回包含单元格所有顶点的边界框，如get_vertices()方法所返回的。对于通过支持点定义的高阶映射，边界盒只保证包含所有支持点，而且一般来说，它只是真正边界盒的近似值，可能更大。
   * @param[in] 单元格 你想计算边界框的单元格
   *
   */
  virtual BoundingBox<spacedim>
  get_bounding_box(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

  /**
   * 返回映射是否保留了顶点位置。换句话说，这个函数返回参考单元格顶点的映射位置（由
   * GeometryInfo::unit_cell_vertex()) 给出）是否等于
   * <code>cell-@>vertex()</code>
   * 的结果（即由三角法存储的信息）。
   * 例如，派生类中的实现对MappingQ、MappingQGeneric、MappingCartesian返回
   * @p true
   * ，但对MappingQEulerian、MappingQ1Eulerian、MappingFEField返回 @p
   * false 。
   *
   */
  virtual bool
  preserves_vertex_locations() const = 0;

  /**
   * 返回这个Mapping实例是否与 @p reference_cell.
   * 中的单元格类型兼容。
   *
   */
  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const = 0;

  /**
   * @name  引用单元格和实数单元格之间的映射点  
     * @{ 
   *
   */

  /**
   * 将单元格上的点 @p p 映射到实数单元格上的相应点 @p
   * cell.   @param  单元格 迭代器到将用于定义映射的单元。
   * @param  p 参考单元格上一个点的位置。    @return
   * 参考点的位置，使用由当前实现映射的派生类所定义的映射映射到实空间，以及第一个参数所确定的单元格的坐标。
   *
   */
  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &                                          p) const = 0;

  /**
   * 将实数 @p p 上的点 @p cell
   * 映射到单元格上的相应点，并返回其坐标。这个函数提供了transform_unit_to_real_cell()所提供的映射的逆映射。
   * 在一维的情况下，本函数返回实点 @p p 在 @p cell.
   * 所标识的曲线或曲面上的法线投影。
   * @note
   * 如果要计算反映射的点位于单元格边界之外，从参考（单位）单元格坐标到实数单元格坐标系的多项式映射并不总是可逆的。在这种情况下，当前函数可能无法计算参考单元上的一个点，该点在映射下的图像等于给定的点
   * @p p.  如果是这种情况，该函数会抛出一个
   * Mapping::ExcTransformationFailed  类型的异常。因此，给定的点
   * @p p
   * 是否位于单元格之外可以通过检查返回的参考坐标是否位于参考单元格之内或之外来确定（例如，使用
   * GeometryInfo::is_inside_unit_cell()) 或上述异常是否被抛出。
   * @param  cell 将用于定义映射的单元的迭代器。    @param  p
   * 给定单元格上的一个点的位置。    @return
   * 点的参考单元位置，当映射到实空间时等于第二个参数给出的坐标。这个映射使用由当前实现映射的派生类所定义的映射，以及第一个参数所确定的单元格的坐标。
   *
   */
  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &                                     p) const = 0;

  /**
   * 将多个点从真实点位置映射到参考位置的点。其功能基本上与在所有点上循环并为每个点单独调用
   * Mapping::transform_real_to_unit_cell()
   * 函数相同，但对于某些实现了更专业的版本的映射，如MappingQGeneric，其速度会更快。行为上的唯一区别是，这个函数永远不会抛出ExcTransformationFailed()异常。如果对`real_points[i]`转换失败，返回的`unit_points[i]`包含
   * std::numeric_limits<double>::infinity() 作为第一个条目。
   *
   */
  virtual void
  transform_points_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<spacedim>> &                    real_points,
    const ArrayView<Point<dim>> &unit_points) const;

  /**
   * 将实数 @p cell 上的点 @p p
   * 转换为参考单元上的对应点，然后将此点投射到给定面数
   * @p face_no.
   * 的面的坐标系中的一个(dim-1)维的点，理想情况下，点 @p
   * p 靠近面 @p face_no,
   * ，但技术上单元中的任何点都可以被投影。
   * 当dim=1时，这个函数没有物理意义，所以在这种情况下它会抛出一个异常。
   *
   */
  Point<dim - 1>
  project_real_point_to_unit_point_on_face(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Point<spacedim> &                                     p) const;

  /**
   * @}
   *
   */


  /**
   * @name  异常情况  @{
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException0(ExcInvalidData);


  /**
   * 计算实空间点和参考空间点之间的映射失败，通常是因为给定的点位于反向映射不唯一的单元之外。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(
    ExcTransformationFailed,
    "Computing the mapping between a real space point and a point in reference "
    "space failed, typically because the given point lies outside the cell "
    "where the inverse mapping is not unique.");

  /**
   * deal.II假设雅各布行列式为正。当单元格的几何形状在映射的图像下被扭曲时，映射变得无效，并抛出这个异常。
   * @ingroup Exceptions
   *
   */
  DeclException3(ExcDistortedMappedCell,
                 Point<spacedim>,
                 double,
                 int,
                 << "The image of the mapping applied to cell with center ["
                 << arg1 << "] is distorted. The cell geometry or the "
                 << "mapping are invalid, giving a non-positive volume "
                 << "fraction of " << arg2 << " in quadrature point " << arg3
                 << ".");

  /**
   * @}
   *
   */

  /**
   * @name  与FEValues的接口  
     * @{ 
   *
   */

public:
  /**
   * 用于映射对象内部数据的基类。内部机制是，在构建FEValues对象时，它要求将要使用的映射和有限元类为自己的目的分配内存，在其中可以存储只需要计算一次的数据。例如，大多数有限元将在这个对象中存储正交点的形状函数值，因为它们不会在单元之间变化，只需要计算一次。对于希望只在正交点评估一次用于映射的形状函数的映射类也是如此。
   * 由于使用不同正交规则的不同FEValues对象可能同时访问同一个映射对象，因此有必要为每个FEValues对象创建一个这样的对象。FEValues通过调用
   * Mapping::get_data(),
   * 或在现实中调用派生类中相应函数的实现来做到这一点。由
   * Mapping::get_data()
   * 创建的对象的所有权随后被转移到FEValues对象中，但每次要求它计算具体单元的信息时，这个对象的引用就会被传递给映射对象。当
   * FEValues::reinit()
   * （或FEFaceValues和FESubfaceValues中的相应类）调用
   * Mapping::fill_fe_values() （以及类似地通过
   * Mapping::fill_fe_face_values() 和 Mapping::fill_fe_subface_values()).
   * 该类的目的是让映射对象存储可以在开始时在参考单元上计算一次的信息，并在以后计算具体单元的信息时访问它。因此，交给
   * Mapping::fill_fe_values() 的对象被标记为 <code>const</code>
   * ，因为假设在使用这些信息的时候，不需要再次修改。然而，从Mapping派生出来的类也可以将这类对象用于其他两个目的。
   *
   *
   *
   *
   *
   * - 为在 Mapping::fill_fe_values() 和类似函数中进行的计算提供划痕空间。一些派生类希望使用从头开始的数组，如果每次调用这个函数时都要分配这些数组，只是在函数结束时再去分配，那就太浪费时间了。相反，可以把这块内存作为当前类的成员变量分配一次，然后在 Mapping::fill_fe_values(). 中简单地使用它。
   *
   *
   *
   *
   *
   * - 在调用 Mapping::fill_fe_values(), 后， FEValues::reinit() 调用 FiniteElement::fill_fe_values() ，其中有限元计算形状函数的值、梯度等，使用与这里描述的机制类似的开始时计算的信息（见 FiniteElement::InternalDataBase)  作为其工作的一部分， FiniteElement::fill_fe_values() 的一些实现需要转换形状函数数据，它们通过调用 Mapping::transform(). 来实现。对后者的调用也会收到对 Mapping::InternalDataBase 对象的引用。由于 Mapping::transform() 可能会在每个单元上被多次调用，有时值得派生类在 Mapping::fill_fe_values() 中只计算一次某些信息，并在 Mapping::transform(). 中重复使用，这些信息也可以存储在派生映射类从InternalDataBase派生的类中。    在这两种情况下，被传递的InternalDataBase对象都是 "道德上的约束"，也就是说，任何外部观察者都无法判断 Mapping::transform() 的抓取数组或一些中间数据是否被 Mapping::fill_fe_values() 所修改。因此，InternalDataBase对象总是作为 <code>const</code> 对象被传递。因此，想要利用上述两种额外用途的派生类需要将他们想要用于这些目的的成员变量标记为 <code>mutable</code> ，以允许它们被修改，尽管周围的对象被标记为 <code>const</code>  。
   *
   */
  class InternalDataBase
  {
  public:
    /**
     * 构造函数。设置update_flags为 @p update_default ， @p
     * first_cell 为 @p true.  。
     *
     */
    InternalDataBase();

    /**
     * 禁止复制构造。
     *
     */
    InternalDataBase(const InternalDataBase &) = delete;

    /**
     * 派生类的虚拟解构器
     *
     */
    virtual ~InternalDataBase() = default;

    /**
     * 一组更新标志，指定Mapping接口的实现需要在每个单元或面计算的信息种类，即在
     * Mapping::fill_fe_values() 和朋友圈。        这组标志被
     * Mapping::get_data(),  Mapping::get_face_data(), 或
     * Mapping::get_subface_data(),
     * 的实现存储在这里，是传递给那些需要对每个单元进行重新计算的函数的更新标志的子集。
     * (对应于在调用 Mapping::get_data()
     * 时已经可以一次性计算的信息的标志子集。
     *
     * - 或该接口的实现
     *
     * - 不需要存储在这里，因为它已经被处理过了)。
     *
     */
    UpdateFlags update_each;

    /**
     * 返回这个对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const;
  };


protected:
  /**
   * 给定一组更新标志，计算哪些其他的量<i>also</i>需要被计算，以满足给定标志的请求。
   * 然后返回原始标志集和刚刚计算的标志的组合。
   * 举个例子，如果 @p update_flags
   * 包含update_JxW_values（即雅各布式的行列式和正交公式提供的权重的乘积），一个映射可能需要计算完整的雅各布式矩阵，以便计算其行列式。然后他们将不仅返回update_JxW_values，而且还返回update_jacobians。在计算JxW值的派生类中，内部实际上不是这样做的
   *
   * - 他们设置了update_contravariant_transformation来代替，由此也可以计算出行列式。
   *
   * - 但这并不影响这个例子的启发性）。)     关于这个函数和FEValues之间的互动的广泛讨论可以在 @ref FE_vs_Mapping_vs_FEValues 文档模块中找到。      @see  UpdateFlags
   *
   */
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const = 0;

  /**
   * 创建并返回一个指向对象的指针，映射可以在该对象中存储数据，这些数据只需要计算一次，但在映射应用于具体单元时都可以使用（例如，在各种transform()函数中，以及构成映射与FEValues类接口的fill_fe_values()、fill_fe_face_values()和fill_fe_subface_values()中）。
   * 派生类将返回指向从 Mapping::InternalDataBase
   * 派生的类型的对象的指针（更多信息见那里），并且可能已经预先计算了一些信息（根据未来对映射的要求，由更新标志指定）和给定的正交对象。随后对transform()或fill_fe_values()和friends的调用将收到这里创建的对象（具有相同的更新标志集和相同的正交对象）。因此，派生类可以在其get_data()函数中预先计算一些信息，并将其存储在内部数据对象中。
   * 映射类不会跟踪由该函数创建的对象。因此，所有权将归调用者所有。
   * 关于这个函数和FEValues之间的互动的广泛讨论可以在 @ref
   * FE_vs_Mapping_vs_FEValues 文档模块中找到。      @param
   * update_flags
   * 一组标志，定义了在未来调用transform()或fill_fe_values()函数组时对映射类的期望。这组标志可能包含映射不知道如何处理的标志（例如，对于事实上由有限元类计算的信息，如
   * UpdateFlags::update_values).
   * 派生类将需要存储这些标志，或者至少是需要映射在fill_fe_values()中执行任何操作的标志子集，在
   * InternalDataBase::update_each.   @param  quadrature
   * 必须计算映射信息的正交对象。这包括正交点的位置和权重。
   * @return
   * 一个指向新创建的InternalDataBase类型（或派生类）对象的指针。该对象的所有权转移给调用函数。
   * @note
   * C++允许派生类中的虚拟函数可以返回不是InternalDataBase类型的对象的指针，但实际上是指向InternalDataBase的类<i>derived</i>的对象的指针。(这个特性被称为
   * "共变返回类型"。)这在某些情况下是很有用的，因为在派生类中的调用将立即使用返回的对象，知道它的真实（派生）类型。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_data(const UpdateFlags      update_flags,
           const Quadrature<dim> &quadrature) const = 0;

  /**
   * 像get_data()一样，但是为以后调用transform()或fill_fe_face_values()做准备，这些调用需要关于从参考面到具体单元面的映射信息。
   * @param  update_flags
   * 一组标志，定义在未来调用transform()或fill_fe_values()函数组时对映射类的期望。这组标志可能包含映射不知道如何处理的标志（例如，对于事实上由有限元类计算的信息，如
   * UpdateFlags::update_values).
   * 派生类将需要存储这些标志，或者至少是需要映射在fill_fe_values()中执行任何操作的标志子集，在
   * InternalDataBase::update_each.   @param  quadrature
   * 需要计算映射信息的正交对象。这包括正交点的位置和权重。
   * @return
   * 一个指向新创建的InternalDataBase类型（或派生类）对象的指针。该对象的所有权转移给调用函数。
   * @note
   * C++允许派生类中的虚拟函数可以返回不是InternalDataBase类型的对象的指针，但实际上是指向InternalDataBase的类<i>derived</i>的对象的指针。(这个特性被称为
   * "共变返回类型"。)这在某些情况下是很有用的，因为在派生类中的调用将立即使用返回的对象，知道它的真实（派生）类型。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_face_data(const UpdateFlags               update_flags,
                const hp::QCollection<dim - 1> &quadrature) const;

  /**
   * @deprecated  使用带有 hp::QCollection 参数的版本。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_face_data(const UpdateFlags          update_flags,
                const Quadrature<dim - 1> &quadrature) const;

  /**
   * 像get_data()和get_face_data()一样，但是为以后调用transform()或fill_fe_subface_values()做准备，这些调用将需要关于从参考面到具体单元的面的子（即子面）的映射信息。
   * @param  update_flags
   * 一组标志，定义在未来调用transform()或fill_fe_values()函数组时对映射类的期望。这组标志可能包含映射不知道如何处理的标志（例如，对于事实上由有限元类计算的信息，如
   * UpdateFlags::update_values).
   * 派生类将需要存储这些标志，或者至少是需要映射在fill_fe_values()中执行任何操作的标志子集，在
   * InternalDataBase::update_each.   @param  quadrature
   * 必须计算映射信息的正交对象。这包括正交点的位置和权重。
   * @return
   * 一个指向新创建的InternalDataBase类型（或派生类）对象的指针。该对象的所有权转移给调用函数。
   * @note
   * C++允许派生类中的虚拟函数可以返回不是InternalDataBase类型的对象的指针，但实际上是指向InternalDataBase的类<i>derived</i>的对象的指针。(这个特性被称为
   * "共变返回类型"。)这在某些情况下是很有用的，在这些情况下，调用是在派生类中，并且将立即使用返回的对象，知道它的真实（派生）类型。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_subface_data(const UpdateFlags          update_flags,
                   const Quadrature<dim - 1> &quadrature) const = 0;

  /**
   * 计算从参考单元到此函数的第一个参数所指示的实数单元的映射信息。派生类将不得不根据它们所代表的映射类型来实现这个函数。它被
   * FEValues::reinit().
   * 调用。从概念上讲，这个函数代表了从参考坐标
   * $\mathbf\in [0,1]^d$  到实空间坐标  $\mathbf x$  的映射  $K$
   * 的应用。它的目的是计算以下种类的数据。
   *
   *
   *
   *
   *
   *
   * - 从应用映射本身产生的数据，例如，计算实数单元上正交点的位置 $\mathbf x_q = \mathbf F_K(\hat{\mathbf x}_q)$ ，对FEValues的用户直接有用，例如在装配过程中。
   *
   *
   *
   *
   *
   * - 数据是有限元实现在真实单元上计算其形状函数所必需的。为此， FEValues::reinit() 函数在当前函数后调用 FiniteElement::fill_fe_values() ，该函数的输出作为 FiniteElement::fill_fe_values(). 的输入。这里需要计算的信息的例子是映射的Jacobian， $\hat\nabla \mathbf F_K(\hat{\mathbf x})$ 或其逆，例如，将参考单元上的形状函数的梯度转换为实单元上形状函数的梯度。    这个函数计算出来的信息被用来填充这个函数的输出参数的各个成员变量。该结构中的哪些成员变量应该被填充，由存储在传递给该函数的 Mapping::InternalDataBase 对象中的更新标志决定。    关于此函数和FEValues之间的互动的广泛讨论可以在 @ref FE_vs_Mapping_vs_FEValues 文档模块中找到。      @param[in]  cell 三角形中的单元格，本函数要计算从参考单元格到的映射。    @param[in]  cell_similarity 作为第一个参数的单元格是否是最近一次调用此函数的单元格的简单平移、旋转等。这个信息是通过匹配前一个单元和当前单元之间的顶点（由三角结构存储）简单计算出来的。这里传递的值可能被这个函数的实现所修改，然后应该被返回（见关于这个函数的返回值的讨论）。    @param[in]  quadrature 对当前评估中使用的正交公式的引用。这个正交对象与创建 @p internal_data 对象时使用的对象相同。该对象既用于映射正交点的位置，也用于计算每个正交点的JxW值（涉及正交权重）。    @param[in]  internal_data 一个对先前由get_data()创建的对象的引用，可用于存储映射在参考单元上可以计算一次的信息。参见 Mapping::InternalDataBase 类的文档，以了解这些对象的用途的广泛描述。    @param[out]  output_data 对成员变量应被计算的对象的引用。并非所有这个参数的成员都需要被填充；哪些成员需要被填充是由存储在 @p internal_data 对象内的更新标志决定的。    @return  这个函数的 @p cell_similarity 参数的一个更新值。当 FEValues::reinit() 调用 FiniteElement::fill_fe_values(). 时，返回的值将被用于相应的参数。在大多数情况下，派生类只想返回为 @p cell_similarity传递的值。然而，这个函数的实现可能会降低细胞相似度的级别。例如，对于那些不仅考虑到单元格顶点的位置（如Triangulation所报告的），而且还考虑到映射的其他特定信息的类，就是这种情况。目的是 FEValues::reinit() 可以只根据单元格的顶点来计算一个单元格是否与前一个单元格相似，而映射也可以考虑位移场（例如，在MappingQ1Eulerian和MappingFEField类中）。在这种情况下，映射可能会得出结论，先前计算的单元格相似度过于乐观，并通过返回一个不那么乐观的单元格相似度值，使其在随后的使用中无效 FiniteElement::fill_fe_values()  。
   * @note  FEValues确保这个函数总是用同一对 @p internal_data 和
   * @p output_data
   * 对象调用。换句话说，如果这个函数的实现知道它在之前的调用中已经把一个数据写入了输出参数，那么在以后的调用中，如果实现知道这是同一个值，就没有必要再把它复制到那里。
   *
   */
  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const = 0;

  /**
   * 这个函数等同于 Mapping::fill_fe_values(),
   * ，但用于单元格的面。有关其目的的广泛讨论，请参见那里。它被
   * FEFaceValues::reinit().  @param[in]
   * 单元格所调用，该函数要计算从参考单元格到的映射。
   * @param[in]  face_no 请求提供信息的给定单元的面的编号。
   * @param[in]  quadrature
   * 当前评估中使用的正交公式的引用。此正交对象与创建
   * @p internal_data
   * 对象时使用的对象相同。该对象既用于映射正交点的位置，也用于计算每个正交点的JxW值（涉及正交权重）。
   * @param[in]  internal_data
   * 一个对先前由get_data()创建的对象的引用，可用于存储映射在参考单元上可以计算一次的信息。参见
   * Mapping::InternalDataBase
   * 类的文档，以了解这些对象的用途的广泛描述。
   * @param[out]  output_data
   * 对成员变量应被计算的对象的引用。并非这个参数的所有成员都需要被填充；哪些成员需要被填充是由存储在
   * @p internal_data 对象内的更新标志决定的。
   *
   */
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * @deprecated  使用带有 hp::QCollection 参数的版本。
   *
   */
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * 这个函数等同于 Mapping::fill_fe_values(),
   * ，但适用于单元格的子面（即面的子女）。关于其目的的广泛讨论见那里。它被
   * FESubfaceValues::reinit().   @param[in]  cell
   * 三角形中的单元格，这个函数要为其计算从参考单元格到的映射。
   * @param[in]  face_no 请求提供信息的给定单元的面的编号。
   * @param[in]  subface_no
   * 请求提供信息的给定单元的面的子的编号。    @param[in]
   * quadrature
   * 对当前评估中使用的正交公式的引用。这个正交对象与创建
   * @p internal_data
   * 对象时使用的对象相同。该对象既用于映射正交点的位置，也用于计算每个正交点的JxW值（涉及正交权重）。
   * @param[in]  internal_data
   * 一个对先前由get_data()创建的对象的引用，可用于存储映射在参考单元上可以计算一次的信息。参见
   * Mapping::InternalDataBase
   * 类的文档，以了解这些对象的用途的广泛描述。
   * @param[out]  output_data
   * 对成员变量应被计算的对象的引用。并非所有这个参数的成员都需要被填充；哪些成员需要被填充是由存储在
   * @p internal_data 对象内的更新标志决定的。
   *
   */
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const = 0;

  /**
   * @}
   *
   */

public:
  /**
   * @name  将张量从参考坐标转换为实数坐标的函数  
     * @{ 
   *
   */

  /**
   * 根据所选的MappingKind对矢量或1-差分形式的场进行变换。
   * @note  通常情况下，这个函数被一个有限元调用，填充FEValues对象。对于这个有限元，应该有一个别名MappingKind，如 @p mapping_bdm,   @p mapping_nedelec,  等。这个别名应该优先于使用下面的种类。    目前由派生类实现的映射种类是。    <ul>   <li>   @p mapping_contravariant:  通过Jacobian将参考单元上的矢量场映射到物理单元。  @f[
   * \mathbf u(\mathbf x) = J(\hat{\mathbf  x})\hat{\mathbf  u}(\hat{\mathbf
   * x}).
   * @f]
   * 在物理学中，这通常被称为反变量变换。在数学上，它是一个矢量场的前推。
   * <li>   @p mapping_covariant:
   * 将参考单元上的一形场映射到物理单元上的一形场。理论上，这将指的是DerivativeForm<1,dim,1>，但我们将这种类型与Tensor<1,dim>进行规范性的识别）。在数学上，它是微分形式的回拉@f[
   * \mathbf u(\mathbf x) = J(\hat{\mathbf  x})(J(\hat{\mathbf  x})^{T}
   * J(\hat{\mathbf  x}))^{-1}\hat{\mathbf u}(\hat{\mathbf  x}).
   * @f]标量可微分函数的梯度是这样转化的。    在dim=spacedim的情况下，前面的公式简化为@f[
   * \mathbf u(\mathbf x) = J(\hat{\mathbf  x})^{-T}\hat{\mathbf
   * u}(\hat{\mathbf  x})
   * @f]，因为我们假设映射 $\mathbf F_K$
   * 总是可逆的，因此其雅各布 $J$ 是一个可逆矩阵。
   * <li>   @p mapping_piola:
   * 参考单元上的<i>dim-1</i>形式的场也由矢量场表示，但同样以不同的方式变换，即通过皮奥拉变换@f[
   * \mathbf u(\mathbf x) = \frac{1}{\text{det}\;J(\hat{\mathbf x})}
   * J(\hat{\mathbf x}) \hat{\mathbf  u}(\hat{\mathbf x}).
   * @f]  </ul>   @param[in]  输入 一个应该被映射的输入对象的数组（或数组的一部分）。    @param[in]  kind 要应用的映射的种类。    @param[in]  internal 一个指向 Mapping::InternalDataBase 类型的对象的指针，该对象包含先前由映射存储的信息。指向的对象是由get_data()、get_face_data()或get_subface_data()函数创建的，在调用当前函数之前，将作为对当前单元格的fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()调用的一部分而被更新。换句话说，这个对象也代表了与哪个单元格有关的变换应该被应用。    @param[out]  输出 一个数组（或数组的一部分），转换后的对象应该被放入其中。(注意，数组视图是 @p 常数，但它所指向的张量不是。)
   *
   */
  virtual void
  transform(const ArrayView<const Tensor<1, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<1, spacedim>> &output) const = 0;

  /**
   * 将一个微分形式的场从参考单元转换到物理单元。 认为 $\mathbf{T} = \nabla \mathbf u$ 和 $\hat{\mathbf  T} = \hat \nabla \hat{\mathbf  u}$ 是有用的， $\mathbf u$ 是一个矢量场。 目前由派生类实现的映射种类有。    <ul>   <li>   @p mapping_covariant:  将参考单元上的形式域映射到物理单元上的形式域。在数学上，它是微分形式的回拉@f[
   * \mathbf T(\mathbf x) = \hat{\mathbf  T}(\hat{\mathbf  x}) J(\hat{\mathbf
   * x})(J(\hat{\mathbf  x})^{T} J(\hat{\mathbf  x}))^{-1}.
   * @f]间隔向量值微分函数的雅各布斯是这样转换的。    在dim=spacedim的情况下，前面的公式简化为@f[
   * \mathbf T(\mathbf x) = \hat{\mathbf  u}(\hat{\mathbf  x}) J(\hat{\mathbf
   * x})^{-1}.
   * @f]  </ul>  。
   * @note
   * 如果把这个变换变成一个模板函数，其等级在<code>DerivativeForm
   * @<1,  dim, rank  @></code>.
   * 中会更合理，可惜C++不允许模板化虚拟函数。这就是为什么我们在这个函数transform()上面使用mapping_covariant()时，将<code>DerivativeForm
   * @<1,  dim, 1  @></code> 标识为 <code>Tensor@<1,dim@></code> 。
   * @param[in]  input
   * 应该被映射的输入对象的一个数组（或数组的一部分）。
   * @param[in]  kind 要应用的映射的种类。    @param[in]  internal
   * 一个指向 Mapping::InternalDataBase
   * 类型的对象的指针，该对象包含先前由映射存储的信息。指向的对象是由get_data()、get_face_data()或get_subface_data()函数创建的，在调用当前函数之前，将作为对当前单元格的fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()调用的一部分被更新。换句话说，这个对象也代表了与哪个单元格有关的变换应该被应用。
   * @param[out]  输出
   * 一个数组（或数组的一部分），转换后的对象应该被放入其中。(注意，数组视图是
   * @p 常数，但它所指向的张量不是。)
   *
   */
  virtual void
  transform(const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const = 0;

  /**
   * 将一个张量场从参考单元转换到物理单元。  这些张量通常是参考单元中已经从物理单元拉回来的矢量场的雅各布系数。 目前由派生类实现的映射种类有。    <ul>   <li>   @p mapping_contravariant_gradient: 它假设 $\mathbf u(\mathbf x)
   * = J \hat{\mathbf  u}$ 这样@f[
   * \mathbf T(\mathbf x) =
   * J(\hat{\mathbf  x}) \hat{\mathbf  T}(\hat{\mathbf  x})
   * J(\hat{\mathbf  x})^{-1}.
   * @f]  <li>   @p mapping_covariant_gradient: 它假设 $\mathbf u(\mathbf x) =
   * J^{-T} \hat{\mathbf  u}$ 这样@f[
   * \mathbf T(\mathbf x) =
   * J(\hat{\mathbf  x})^{-T} \hat{\mathbf  T}(\hat{\mathbf  x})
   * J(\hat{\mathbf  x})^{-1}.
   * @f]  <li>   @p mapping_piola_gradient: 它假设 $\mathbf u(\mathbf x) =
   * \frac{1}{\text{det}\;J(\hat{\mathbf x})} J(\hat{\mathbf x}) \hat{\mathbf
   * u}(\hat{\mathbf x})$ 这样 @f[
   * \mathbf T(\mathbf x) =
   * \frac{1}{\text{det}\;J(\hat{\mathbf x})}
   * J(\hat{\mathbf  x}) \hat{\mathbf  T}(\hat{\mathbf  x})
   * J(\hat{\mathbf  x})^{-1}.
   * @f] ]  </ul>   @todo  mapping_covariant_gradient、mapping_contravariant_gradient和mapping_piola_gradient的公式只对线性映射而言是真的。例如，如果映射是双线性的（或具有高阶多项式程度），那么就会有一个与  $J$  的导数相关的缺失项。      @param[in]  输入 一个应该被映射的输入对象的数组（或数组的一部分）。    @param[in]  kind 要应用的映射的种类。    @param[in]  internal 一个指向 Mapping::InternalDataBase 类型的对象的指针，该对象包含先前由映射存储的信息。指向的对象是由get_data()、get_face_data()或get_subface_data()函数创建的，在调用当前函数之前，将作为对当前单元格的fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()调用的一部分被更新。换句话说，这个对象也代表了与哪个单元格有关的变换应该被应用。    @param[out]  output 一个数组（或数组的一部分），转换后的对象应该被放入其中。(注意，数组视图是 @p 常数，但它所指向的张量不是。)
   *
   */
  virtual void
  transform(const ArrayView<const Tensor<2, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const = 0;

  /**
   * 将一个张量场从参考单元转换到物理单元。  这种张量在大多数情况下是参考单元中的向量场的 hessians，这些向量场已经从物理单元拉回来。    目前由派生类实现的映射种类有。    <ul>   <li>   @p mapping_covariant_gradient:  将参考单元上的形式场映射到物理单元上的形式场。在数学上，它是微分形式@f[
   * \mathbf T_{ijk}(\mathbf x) = \hat{\mathbf  T}_{iJK}(\hat{\mathbf  x})
   * J_{jJ}^{\dagger} J_{kK}^{\dagger}@f]的回拉，其中@f[ J^{\dagger} = J(\hat{\mathbf  x})(J(\hat{\mathbf  x})^{T}
   * J(\hat{\mathbf  x}))^{-1}.
   * @f]  </ul> 间隔向量值可微函数的Hessians是这样转化的（在减去导数与雅各布梯度的乘积后）。    在dim=spacedim的情况下，前面的公式简化为 @f[J^{\dagger} = J^{-1}@f]  @param[in]  input 应该被映射的输入对象的数组（或数组的一部分）。    @param[in]  kind 要应用的映射的种类。    @param[in]  internal 一个指向类型为 Mapping::InternalDataBase 的对象的指针，该对象包含先前由映射存储的信息。指向的对象是由get_data()、get_face_data()或get_subface_data()函数创建的，在调用当前函数之前，将作为对当前单元格的fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()调用的一部分被更新。换句话说，这个对象也代表了与哪个单元格有关的变换应该被应用。    @param[out]  output 一个数组（或数组的一部分），转换后的对象应该被放入其中。(注意，数组视图是 @p 常数，但它指向的张量不是。)
   *
   */
  virtual void
  transform(const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const = 0;

  /**
   * 将一个3差分形式的场从参考单元转换到物理单元。
   * 认为 $\mathbf{T}_{ijk} = D^2_{jk} \mathbf u_i$ 和 $\mathbf{\hat
   * T}_{IJK} = \hat D^2_{JK} \mathbf{\hat
   * u}_I$ 很有用， $\mathbf u_i$ 是一个矢量场。    目前由派生类实现的映射种类是。    <ul>   <li>   @p mapping_contravariant_hessian:  它假定 $\mathbf u_i(\mathbf x)
   * = J_{iI} \hat{\mathbf  u}_I$ 以便@f[
   * \mathbf T_{ijk}(\mathbf x) =
   * J_{iI}(\hat{\mathbf  x}) \hat{\mathbf  T}_{IJK}(\hat{\mathbf  x})
   * J_{jJ}(\hat{\mathbf  x})^{-1} J_{kK}(\hat{\mathbf  x})^{-1}.
   * @f]  <li>   @p mapping_covariant_hessian:  它假定 $\mathbf u_i(\mathbf x) =
   * J_{iI}^{-T} \hat{\mathbf  u}_I$  以便@f[
   * \mathbf T_{ijk}(\mathbf x) =
   * J_iI(\hat{\mathbf  x})^{-1} \hat{\mathbf  T}_{IJK}(\hat{\mathbf  x})
   * J_{jJ}(\hat{\mathbf  x})^{-1} J_{kK}(\hat{\mathbf  x})^{-1}.
   * @f]  <li>   @p mapping_piola_hessian:  ] 它假定  $\mathbf u_i(\mathbf x) =
   * \frac{1}{\text{det}\;J(\hat{\mathbf x})} J_{iI}(\hat{\mathbf x})
   * \hat{\mathbf u}(\hat{\mathbf x})$  这样 @f[
   * \mathbf T_{ijk}(\mathbf x) =
   * \frac{1}{\text{det}\;J(\hat{\mathbf x})}
   * J_{iI}(\hat{\mathbf  x}) \hat{\mathbf  T}_{IJK}(\hat{\mathbf  x})
   * J_{jJ}(\hat{\mathbf  x})^{-1} J_{kK}(\hat{\mathbf  x})^{-1}.
   * @f]  </ul>   @param[in]  输入 一个应该被映射的输入对象的数组（或数组的一部分）。    @param[in]  kind 要应用的映射的种类。    @param[in]  internal 一个指向 Mapping::InternalDataBase 类型的对象的指针，该对象包含先前由映射存储的信息。指向的对象是由get_data()、get_face_data()或get_subface_data()函数创建的，在调用当前函数之前，将作为对当前单元格的fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()调用的一部分被更新。换句话说，这个对象也代表了与哪个单元格有关的变换应该被应用。    @param[out]  output 一个数组（或数组的一部分），转换后的对象应该被放入其中。
   *
   */
  virtual void
  transform(const ArrayView<const Tensor<3, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const = 0;

  /**
   * @}
   *
   */


  // Give class @p FEValues access to the private <tt>get_...data</tt> and
  // <tt>fill_fe_...values</tt> functions.
  friend class FEValuesBase<dim, spacedim>;
  friend class FEValues<dim, spacedim>;
  friend class FEFaceValues<dim, spacedim>;
  friend class FESubfaceValues<dim, spacedim>;
};


/**
 * 返回一个适用于给定三角形的默认线性映射。在内部，这个函数为给定的三角结构所使用的参考单元调用上述函数，假设三角结构只使用单一的单元类型。如果三角剖分使用混合单元格类型，那么这个函数将触发一个异常。
 *
 *
 */
template <int dim, int spacedim>
const Mapping<dim, spacedim> &
get_default_linear_mapping(const Triangulation<dim, spacedim> &triangulation);


DEAL_II_NAMESPACE_CLOSE

#endif


