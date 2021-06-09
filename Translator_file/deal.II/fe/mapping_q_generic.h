//include/deal.II-translator/fe/mapping_q_generic_0.txt
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

#ifndef dealii_mapping_q_generic_h
#define dealii_mapping_q_generic_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <array>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int, int>
class MappingQ;

template <int, int>
class MappingQCache;


 /*!@addtogroup mapping */ 
 /*@{*/ 


/**
 * 这个类实现了多项式映射 $Q_p$ 的功能，其多项式程度 $p$
 * 将被用于网格的所有单元。MappingQ1和MappingQ类对这一行为进行了轻微的专业化。
 * 这个类的名字很糟糕。它确实应该被称为MappingQ，因为它在三角形的所有单元上持续使用
 * $Q_p$
 * 映射。然而，当我们为映射重写整个类的层次结构时，MappingQ这个名字已经被使用了。人们可能会争论说，应该总是使用MappingQGeneric而不是现有的MappingQ类（除非在构造对象时明确指定，否则它只使用程度为
 * $p$  <i>on cells at the boundary of the
 * domain</i>的映射）。另一方面，在很多情况下有很好的理由使用MappingQ：在很多情况下，曲线域只提供了关于边界处的边缘到底是如何形成的信息，但我们对内部的边缘一无所知。因此，在没有其他信息的情况下，我们只能假设内部边缘是直线，在这种情况下，内部单元也可以被视为双线性四边形或三线性六面体。(在
 * step-1 中已经展示了这样的网格的例子，但在 step-6 的
 * "结果 "部分也有讨论)
 * 。由于双线/三线映射的计算成本明显低于高阶映射，在这种情况下，只在域的边界单元上使用高阶映射是有利的。
 *
 * --也就是MappingQ的行为。当然，MappingQGeneric也对内部单元使用双线性映射，只要它不知道内部边缘的曲率，但它以昂贵的方式实现这一点：作为一般的
 * $Q_p$
 * 映射，其中映射支持点只是<i>happen</i>，沿着线性或双线性的边缘或面排列。
 * 有一些特殊情况值得考虑。
 *
 *
 *
 * - 如果你真的想对所有的单元格使用高阶映射，你可以使用当前的类来做，但这只有在你能真正提供关于网格内部边缘和面应该如何弯曲的信息时才有意义。这通常是通过将一个Manifold与内部单元和边缘关联来实现的。一个简单的例子在 step-6 的 "结果 "部分讨论；关于流形的完整讨论在 step-53 中提供。
 *
 *
 *
 * - 如果你正在处理描述嵌入更高空间维度的（弯曲的）流形的网格，即，如果dim!=spacedim，那么每个单元都位于域的边界，你很可能已经为所有单元附加了一个流形对象，然后也可以被映射类用于高阶映射。
 * <h4>Behavior along curved boundaries and with different manifolds</h4>
 * 如上所述，人们往往只知道一个表面的流形描述，而不知道计算域的内部。在这种情况下，一个FlatManifold对象将被分配给内部实体，它描述了一个通常的平面坐标系，其中高阶映射的附加点被准确地按照双/三线性映射放置。当与边界上的非平面流形结合时，例如一个圆凸入一个正方形单元的内部，这两个流形描述一般来说是不相容的。例如，仅通过单元格顶点定义的平坦流形会把内部点放在离边界沿直线的某个小距离epsilon处，因此一般是在圆的凹陷部分之外。如果MappingQ的多项式程度足够高，从参考单元到这样一个单元的转换一般会包含靠近边界的倒置区域。
 * 为了避免这种情况，该类应用了一种算法，利用所谓的转折插值使这种转换变得平滑，这种插值本质上是沿着周围实体的描述之间的线性混合。在计算附加点的算法中，即compute_mapping_support_points()方法，单元格的所有实体都是分层次通过的，从线开始到四边形，最后是六边形。层次结构中更高的对象上的点是从与该对象相关的流形中获得的，同时考虑到之前由与低维对象相关的流形计算的所有点，而不仅仅是顶点。如果只给一条线分配了一个弯曲的边界，但相邻的四边形却在一个平面流形上，那么四边形上的平面流形在插值四边形内的附加点的位置时将考虑到变形线上的点，从而总是导致一个定义明确的变换。
 * 本类中使用的插值方案确保曲线描述可以在单层元素内过度到平面描述，保持有限元插值的整体最佳收敛率。然而，如果随着网格的细化，弯曲边界和平坦内域之间的过渡被分散在更大的范围内，人们往往会得到更好的解质量。这是由特殊流形TransfiniteInterpolationManifold提供的。
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingQGeneric : public Mapping<dim, spacedim>
{
public:
  /**
   * 构造函数。   @p polynomial_degree
   * 表示用于从参考单元映射到实际单元的多项式程度。
   *
   */
  MappingQGeneric(const unsigned int polynomial_degree);

  /**
   * 复制构造函数。
   *
   */
  MappingQGeneric(const MappingQGeneric<dim, spacedim> &mapping);

  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * 返回映射的程度，即传递给构造函数的值。
   *
   */
  unsigned int
  get_degree() const;

  /**
   * 总是返回 @p true
   * ，因为这个类中函数的默认实现保留了顶点位置。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

  // for documentation, see the Mapping base class
  virtual BoundingBox<spacedim>
  get_bounding_box(const typename Triangulation<dim, spacedim>::cell_iterator
                     &cell) const override;

  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const override;

  /**
   * @name  参考单元和实数单元之间的映射点 
     * @{ 
   *
   */

  // for documentation, see the Mapping base class
  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &p) const override;

  // for documentation, see the Mapping base class
  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &p) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform_points_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<spacedim>> &                    real_points,
    const ArrayView<Point<dim>> &unit_points) const override;

  /**
   * @}
   *
   */

  /**
   * @name  将张量从参考坐标转换为实坐标的函数  @{
   *
   */

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<1, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<1, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<2, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<3, dim>> &                  input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  /**
   *
   */

  /**
   * @name  与FEValues和朋友的接口  
     * @{ 
   *
   */

  /**
   * 多项式映射的内部数据的存储。见 Mapping::InternalDataBase
   * 的广泛描述。
   * 对于当前的类，InternalData类存储了对象创建时（在get_data()中）计算一次的数据，以及类希望从调用fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()之间存储的数据，直到以后可能从有限元调用转化()等函数。后一类的成员变量被标记为
   * "可变"。
   *
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 构造函数。参数表示该对象所对应的映射的多项式程度。
     *
     */
    InternalData(const unsigned int polynomial_degree);

    /**
     * 根据给定的参数，初始化对象中与单元格数据相关的成员变量。
     * 该函数还调用compute_shape_function_values()来实际设置与映射形状函数的值和导数有关的成员变量。
     *
     */
    void
    initialize(const UpdateFlags      update_flags,
               const Quadrature<dim> &quadrature,
               const unsigned int     n_original_q_points);

    /**
     * 根据给定的参数，初始化对象中与单元格和面的数据有关的成员变量。为了初始化单元格数据，本函数调用initialize()。
     *
     */
    void
    initialize_face(const UpdateFlags      update_flags,
                    const Quadrature<dim> &quadrature,
                    const unsigned int     n_original_q_points);

    /**
     * 计算用于映射的形状函数的值和/或导数。
     * 哪些值、导数或高阶导数被计算是由哪些成员数组有非零大小决定的。它们通常被initialize()和initialize_face()函数设置为适当的大小，这些函数确实在内部调用这个函数。然而，用手来调整大小，然后直接调用这个函数是可能的（有时也很有用）。一个例子是在牛顿迭代中，我们更新了一个正交点的位置（例如，在
     * MappingQ::transform_real_to_uni_cell())
     * ，需要重新计算这个位置的映射和它的导数，但已经正确调整了所有内部数组的大小。
     *
     */
    void
    compute_shape_function_values(const std::vector<Point<dim>> &unit_points);

    /**
     * 正交点的形状函数。形状函数是按张量积顺序排列的，所以必须对顶点重新排序以获得变换。
     *
     */
    const double &
    shape(const unsigned int qpoint, const unsigned int shape_nr) const;

    /**
     * 正交点的形状函数。见上文。
     *
     */
    double &
    shape(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的梯度。见上文。
     *
     */
    const Tensor<1, dim> &
    derivative(const unsigned int qpoint, const unsigned int shape_nr) const;

    /**
     * 形状函数在正交点的梯度。见上文。
     *
     */
    Tensor<1, dim> &
    derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的二阶导数。见上文。
     *
     */
    const Tensor<2, dim> &
    second_derivative(const unsigned int qpoint,
                      const unsigned int shape_nr) const;

    /**
     * 形状函数在正交点的二阶导数。见上文。
     *
     */
    Tensor<2, dim> &
    second_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的三次导数。见上文。
     *
     */
    const Tensor<3, dim> &
    third_derivative(const unsigned int qpoint,
                     const unsigned int shape_nr) const;

    /**
     * 形状函数在正交点的三阶导数。见上文。
     *
     */
    Tensor<3, dim> &
    third_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 形状函数在正交点的第四次导数。见上文。
     *
     */
    const Tensor<4, dim> &
    fourth_derivative(const unsigned int qpoint,
                      const unsigned int shape_nr) const;

    /**
     * 形状函数在正交点的四次导数。见上文。
     *
     */
    Tensor<4, dim> &
    fourth_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * 返回这个对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 形状函数的值。通过函数访问  @p shape.  计算一次。
     *
     */
    std::vector<double> shape_values;

    /**
     * 形状函数导数的值。通过函数访问  @p derivative.
     * 计算一次。
     *
     */
    std::vector<Tensor<1, dim>> shape_derivatives;

    /**
     * 形状函数二次导数的值。通过函数 @p
     * second_derivative访问。        计算一次。
     *
     */
    std::vector<Tensor<2, dim>> shape_second_derivatives;

    /**
     * 形状函数第三导数的值。通过函数 @p
     * second_derivative访问。        计算一次。
     *
     */
    std::vector<Tensor<3, dim>> shape_third_derivatives;

    /**
     * 形状函数第四导数的值。通过函数 @p
     * second_derivative访问。        计算一次。
     *
     */
    std::vector<Tensor<4, dim>> shape_fourth_derivatives;

    /**
     * 单位切向量。用于计算边界形式和法向量。
     * 这个数组有`(dim-1)  GeometryInfo::faces_per_cell`
     * 条目。第一个 GeometryInfo::faces_per_cell
     * 包含每个面的第一个切向的向量；第二组
     * GeometryInfo::faces_per_cell
     * 条目包含第二个切向的向量（只有在3D中，因为每个面有两个切向），等等。
     * 填充一次。
     *
     */
    std::array<std::vector<Tensor<1, dim>>,
               GeometryInfo<dim>::faces_per_cell *(dim - 1)>
      unit_tangentials;

    /**
     * 映射的多项式程度。由于这里的对象也被MappingQ使用（稍作调整），我们需要存储这个。
     *
     */
    const unsigned int polynomial_degree;

    /**
     * 形状函数的数量。如果这是一个Q1映射，那么它就是简单的每个单元格的顶点数量。然而，由于派生类也使用这个类（例如Mapping_Q()类），形状函数的数量也可能不同。
     * 一般来说，它是 $(p+1)^\text{dim}$  ，其中 $p$
     * 是映射的多项式程度。
     *
     */
    const unsigned int n_shape_functions;

    /* 默认的线支持点。是在计算形状函数值时使用的。        正交点的数量取决于该类的度数，它与FE_Q<1>(this->degree)的自由度数相匹配。   
*
*/
    QGaussLobatto<1> line_support_points;

    /**
     * 如果给定的正交规则代表一个张量乘积，我们需要存储1d正交点上的1d多项式的值。这就是这个变量的作用。
     *
     */
    internal::MatrixFreeFunctions::ShapeInfo<VectorizedArray<double>>
      shape_info;

    /**
     * 如果给定的正交规则代表一个张量积，我们需要在这个对象中存储临时数据。
     *
     */
    mutable AlignedVector<VectorizedArray<double>> scratch;

    /**
     * 如果给定的正交规则代表一个张量积，那么在映射的支持点上的值将存储在这个对象中。
     *
     */
    mutable AlignedVector<VectorizedArray<double>> values_dofs;

    /**
     * 如果给定的正交规则代表一个张量乘积，那么正交点的值将存储在此对象中。
     *
     */
    mutable AlignedVector<VectorizedArray<double>> values_quad;

    /**
     * 如果给定的正交规则代表一个张量乘积，那么正交点的梯度将存储在这个对象中。
     *
     */
    mutable AlignedVector<VectorizedArray<double>> gradients_quad;

    /**
     * 如果给定的正交规则代表张量乘积，则正交点的斜率将存储在此对象中。
     *
     */
    mutable AlignedVector<VectorizedArray<double>> hessians_quad;

    /**
     * 表示给定的正交对象是否是张量积。
     *
     */
    bool tensor_product_quadrature;

    /**
     * 每个正交点的协变的张量。    存储的矩阵是Jacobian
     * G^{-1}，其中G = Jacobian^{t}的
     * Jacobian，是地图的第一基本形式；如果dim=spacedim，则还原为Jacobian矩阵的转置，其本身被存储在该结构的
     * @p contravariant 域中。        在每个单元格上计算。
     *
     */
    mutable std::vector<DerivativeForm<1, dim, spacedim>> covariant;

    /**
     * 每个正交点上的禁忌变换的张量。不变矩阵是变换的雅各布系数，即
     * $J_{ij}=dx_i/d\hat x_j$  。        在每个单元上计算。
     *
     */
    mutable std::vector<DerivativeForm<1, dim, spacedim>> contravariant;

    /**
     * 供内部使用的辅助向量。
     *
     */
    mutable std::vector<std::vector<Tensor<1, spacedim>>> aux;

    /**
     * 在 @p
     * cell_of_current_support_points上存储映射形状函数的支持点。
     *
     */
    mutable std::vector<Point<spacedim>> mapping_support_points;

    /**
     * 存储 @p mapping_support_points 的单元格。
     *
     */
    mutable typename Triangulation<dim, spacedim>::cell_iterator
      cell_of_current_support_points;

    /**
     * 每个正交点中的雅各布系数的行列式。如果#update_volume_elements就会被填满。
     *
     */
    mutable std::vector<double> volume_elements;
  };


  // documentation can be found in Mapping::requires_update_flags()
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  // documentation can be found in Mapping::get_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_data(const UpdateFlags, const Quadrature<dim> &quadrature) const override;

  using Mapping<dim, spacedim>::get_face_data;

  // documentation can be found in Mapping::get_face_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_face_data(const UpdateFlags               flags,
                const hp::QCollection<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::get_subface_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_subface_data(const UpdateFlags          flags,
                   const Quadrature<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::fill_fe_values()
  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  using Mapping<dim, spacedim>::fill_fe_face_values;

  // documentation can be found in Mapping::fill_fe_face_values()
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  // documentation can be found in Mapping::fill_fe_subface_values()
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;


  /**
   * 相对于其他的fill_fe_values()和fill_fe_face_values()函数依赖InternalDataBase的预计算信息，这个函数在传入当前函数的单元格和点上选择灵活的评估路径。
   * @param[in]  cell 要评估映射的单元  @param[in]  unit_points
   * 参考坐标中的点，应该在这里计算变换（Jacobians，位置）。
   * @param[in]  update_flags 应该被计算的信息种类。
   * @param[out]  output_data
   * 一个包含评估量的结构，例如在给定单元上应用映射及其底层流形后产生的雅各布系数。
   *
   */
  void
  fill_mapping_data_for_generic_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<dim>> &                         unit_points,
    const UpdateFlags                                           update_flags,
    dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * @}
   *
   */

protected:
  /**
   * 用作单元格映射的形状函数的多项式的程度。
   *
   */
  const unsigned int polynomial_degree;

  /* 默认的线支持点。在计算线和四边形上的支持点在实空间中的位置时，这些支持点是由Manifold<dim,spacedim>类需要的。    点的数量取决于这个类的程度，它与FE_Q<1>(this->degree)的自由度数量相匹配。 
*
*/
  const std::vector<Point<1>> line_support_points;

  /* 从线支持点定义为拉格朗日多项式的一维多项式。这些用于点评估，与FE_Q<1>(this->degree)的多项式空间相匹配。 
*
*/
  const std::vector<Polynomials::Polynomial<double>> polynomials_1d;

  /* 在扩展与映射支持点（以分层数字形式出现）的张量积时，使用的从词法到分层排序的编号。 
*
*/
  const std::vector<unsigned int> renumber_lexicographic_to_hierarchic;

  /* 参考坐标中的支持点。这些用于构建计算_mapping_support_points()的输出的近似值，而不是通过InternalData提供的FEValues接口来评估映射时。    点的数量取决于这个类的程度，它与FE_Q<dim>(this->degree)的自由度数量相匹配。 
*
*/
  const std::vector<Point<dim>> unit_cell_support_points;

  /**
   * 一个权重表的向量，我们将物体（直线、四边形、六边形）周边的支持点的位置与之相乘，得到内部支持点的位置。
   * 进入该表的方法是 @p [structdim-1],
   * ，即用0来访问直线上的支持点权重（即Gauss-Lobatto正交的内部点），用1来访问从周长到四边形内部的支持点权重，用2来访问从周长到六角形内部的支持点权重。
   * 该表本身包含有多少列，就有多少个特定对象的周边点（2代表直线，
   * <code>4 + 4*(degree-1)</code> 代表四边形， <code>8 + 12*(degree-1)
   * + 6*(degree-1)*(degree-1)</code>
   * 代表六边形）和多少行，就有多少个严格意义上的内部点。
   * 该表的定义见 "映射 "报告的公式（8）。
   *
   */
  const std::vector<Table<2, double>>
    support_point_weights_perimeter_to_interior;

  /**
   * 一个权重表，我们将单元格的顶点位置与之相乘，得到所有额外支持点的位置，包括线、四边形和六边形（根据情况）。这个数据结构是在我们一次性填充所有支持点时使用的，如果一个单元的所有子实体都连接着同一个流形，就会出现这种情况。这样一来，我们就可以避免为映射转换数据时的一些开销。
   * 该表的行数与单元格的顶点数相同（一维为2，二维为4，三维为8），行数与映射中的额外支持点数相同，即：<code>(degree+1)^dim
   *
   * - 2^dim</code>。
   *
   */
  const Table<2, double> support_point_weights_cell;

  /**
   * 返回该映射的支持点的位置。例如，对于 $Q_1$
   * 映射来说，这些是顶点，而对于高阶多项式映射来说，它们是顶点加上边、面和单元格内部的点，这些点是根据域和其边界的Manifold描述而放置。然而，其他类可以用不同的方式覆盖这个函数。特别是，MappingQ1Eulerian类正是这样做的，它不从当前单元的几何形状计算支持点，而是在单元的几何形状之外评估一个外部给定的位移场。
   * 这个函数的默认实现适用于大多数情况。它从底层流形中获取单元边界上的支持点的位置。然后使用低维实体（线、四边形）的内插法计算内部支持点（即二维的四边形支持点，三维的六边形支持点），以使转换尽可能平滑，而不会因为支持点的放置而在单元内引入额外的边界层。
   * 该函数从顶点（取自给定的单元）出发，通过线的支撑点（调用add_line_support_points()函数）和四边形面上的支撑点（在三维中，调用add_quad_support_points()函数）来工作。
   * 然后，它添加内部支持点，这些支持点是通过使用权重从周围的点插值计算出来的，如果dim<spacedim，它就向底层流形询问内部点的位置。
   *
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

  /**
   * 通过牛顿迭代将实单元上的点 @p p
   * 转换为单位单元上的相应点 @p cell 。
   *
   */
  Point<dim>
  transform_real_to_unit_cell_internal(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &                                     p,
    const Point<dim> &initial_p_unit) const;

  /**
   * 将位于给定单元格边界线上的所有形状函数的支持点追加到矢量
   * @p a. 中 位于线的顶点上的点不包括在内。
   * 该函数使用线的底层流形对象（如果没有设置，则使用单元格的底层流形对象）来确定请求的点的位置。这个函数通常由compute_mapping_support_points()函数调用。
   * 这个函数是虚拟的，以便让派生类选择形状函数支持点的方式与本类不同，本类选择的点是边界上的插值点。
   *
   */
  virtual void
  add_line_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim>> &                              a) const;

  /**
   * 将位于给定单元格的边界面（3D中的四边形）上的所有形状函数的支持点附加到矢量上
   * @p a.
   * 这个函数只定义于<tt>dim=3</tt>。位于四边形的顶点或线上的点不包括在内。
   * 该函数使用四边形的底层流形对象（如果没有设置，则使用单元格的底层流形对象）来确定所请求的点的位置。这个函数通常由compute_mapping_support_points()调用。
   * 这个函数是虚拟的，以便允许派生类以不同于本类的方式选择形状函数支持点，本类选择的点是边界上的插值点。
   *
   */
  virtual void
  add_quad_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim>> &                              a) const;

  // Make MappingQ a friend since it needs to call the fill_fe_values()
  // functions on its MappingQGeneric(1) sub-object.
  template <int, int>
  friend class MappingQ;

  // Make MappingQCache a friend since it needs to call the
  // compute_mapping_support_points() function.
  template <int, int>
  friend class MappingQCache;
};



 /*@}*/ 

 /*----------------------------------------------------------------------*/ 

#ifndef DOXYGEN

template <int dim, int spacedim>
inline const double &
MappingQGeneric<dim, spacedim>::InternalData::shape(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr, shape_values.size());
  return shape_values[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline double &
MappingQGeneric<dim, spacedim>::InternalData::shape(const unsigned int qpoint,
                                                    const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr, shape_values.size());
  return shape_values[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<1, dim> &
MappingQGeneric<dim, spacedim>::InternalData::derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_derivatives.size());
  return shape_derivatives[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline Tensor<1, dim> &
MappingQGeneric<dim, spacedim>::InternalData::derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_derivatives.size());
  return shape_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<2, dim> &
MappingQGeneric<dim, spacedim>::InternalData::second_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_second_derivatives.size());
  return shape_second_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<2, dim> &
MappingQGeneric<dim, spacedim>::InternalData::second_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_second_derivatives.size());
  return shape_second_derivatives[qpoint * n_shape_functions + shape_nr];
}

template <int dim, int spacedim>
inline const Tensor<3, dim> &
MappingQGeneric<dim, spacedim>::InternalData::third_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_third_derivatives.size());
  return shape_third_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<3, dim> &
MappingQGeneric<dim, spacedim>::InternalData::third_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_third_derivatives.size());
  return shape_third_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<4, dim> &
MappingQGeneric<dim, spacedim>::InternalData::fourth_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_fourth_derivatives.size());
  return shape_fourth_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<4, dim> &
MappingQGeneric<dim, spacedim>::InternalData::fourth_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_fourth_derivatives.size());
  return shape_fourth_derivatives[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline bool
MappingQGeneric<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}

#endif // DOXYGEN

 /* -------------- declaration of explicit specializations ------------- */ 


DEAL_II_NAMESPACE_CLOSE

#endif


