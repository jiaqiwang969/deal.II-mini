//include/deal.II-translator/fe/mapping_q_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_q_h
#define dealii_mapping_q_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup mapping */ 
 /*@{*/ 

/**
 * 一个在域的边界单元（或者，如果在构造函数中要求，对所有单元）上实现度数为
 * $p$ 的多项式映射的类，以及对内部单元的线性映射。
 * 该类实际上命名不当，因为（除非在构造对象时明确指定，见下文），它实际上不使用度数为
 * $p$
 * <i>everywhere</i>的映射，而只在边界的单元上使用。这与MappingQGeneric类形成对比，后者确实在任何地方都使用度数为
 * $Q_p$ 的多项式映射 $p$
 * 。当前类的重点是，在很多情况下，曲线域只提供了关于边界处的边缘到底是如何形成的信息，但我们对内部的边缘一无所知。因此，在没有其他信息的情况下，我们只能假设内部边缘是直线，在这种情况下，内部单元也可以被视为双线性四边形或三线性六面体。(在
 * step-1 中已经展示了这样的网格的例子，但在 step-6 的
 * "结果 "部分也有讨论)
 * 。因为双线/三线映射的计算成本明显低于高阶映射，在这种情况下，只在域的边界单元上使用高阶映射是有利的。这个类正好实现了这种行为。
 * 有一些特殊情况值得考虑。
 *
 *
 *
 * - 如果你想对所有单元使用高阶映射，你可以通过设置构造函数的第二个参数为true来实现。这只有在你能实际提供关于网格内部边缘和面应该如何弯曲的信息时才有意义。这通常是通过将Manifold与内部单元格和边缘关联来实现的。一个简单的例子在 step-6 的 "结果 "部分讨论；关于流形的完整讨论在 step-53 中提供。
 *
 *
 *
 * - 如果你将true作为该类的第二个参数，那么实际上完全等同于立即生成一个MappingQGeneric对象。
 *
 *
 *
 * - 如果提供的多项式程度是1，这个类也完全等同于MappingQGeneric。这是因为在这种情况下，在域的内部和边界的单元上使用的映射不能被区分。
 *
 *
 *
 * - 如果你正在处理嵌入更高空间维度的网格，也就是说，如果dim!=spacedim，那么每个单元都被认为是在域的边界，因此对所有单元都使用高阶映射；同样，这个类也相当于立即使用MappingQGeneric。
 * <h4>Behavior along curved boundaries and with different manifolds</h4>
 * 关于不同流形混合情况下的映射行为和收敛率，请参考MappingQGeneric的相关章节。
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingQ : public Mapping<dim, spacedim>
{
public:
  /**
   * 构造函数。   @p polynomial_degree
   * 表示用于映射单元格边界的多项式的程度。
   * 第二个参数决定了高阶映射是否也应该用于内部单元。如果它的值是
   * <code>false</code>
   * （默认值），那么在内部就使用低阶映射。这对大多数情况来说是足够的，因为高阶映射只是用来更好地接近边界。在这种情况下，内部以直线为界的单元是可以接受的。然而，在有些情况下，人们也希望在内部使用高阶映射。MappingQEulerian类就是这样一种情况。
   * 如果 @p dim 不等于 @p spacedim, ， @p use_mapping_q_on_all_cells
   * 的值就会被忽略，也就是说，如果我们考虑的是嵌入高维空间的表面的网格。
   *
   */
  MappingQ(const unsigned int polynomial_degree,
           const bool         use_mapping_q_on_all_cells = false);

  /**
   * 复制构造器。
   *
   */
  MappingQ(const MappingQ<dim, spacedim> &mapping);

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
   * 将单元格上的点 @p p 转换为实单元格 @p cell 上的点 @p
   * p_real ，并返回 @p p_real.  。
   *
   */
  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &p) const override;

  /**
   * 将实细胞上的点 @p p 转换到单元格 @p cell 上的点 @p p_unit
   * ，并返回 @p p_unit.  使用牛顿迭代和 @p
   * transform_unit_to_real_cell 函数。
   * 在一维情况下，该函数返回实点 @p p 在 @p cell.
   * 标识的曲线或曲面上的法线投影。
   * @note
   * 如果要计算反映射的点位于单元格边界之外，从参考（单位）单元格坐标到实数单元格坐标系的多项式映射并不总是可逆的。
   * 在这种情况下，当前函数可能无法计算出参考单元上的一个点，该点在映射下的图像等于给定的点
   * @p p.  如果是这种情况，那么这个函数会抛出一个
   * Mapping::ExcTransformationFailed  类型的异常。
   * 因此，给定的点 @p p
   * 是否位于单元格之外，可以通过检查返回的参考坐标是否位于参考单元格之内或之外来确定（例如，使用
   * GeometryInfo::is_inside_unit_cell) 或是否抛出了上述的异常。
   *
   */
  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &p) const override;

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
   * 返回一个指向当前对象副本的指针。然后，这个副本的调用者就拥有了它的所有权。
   *
   */

  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;


  /**
   * @name  与FEValues的接口 
     * @{ 
   *
   */

protected:
  /**
   * 该映射的内部数据的存储。见 Mapping::InternalDataBase
   * 的广泛描述。
   * 这包括在创建对象时（在get_data()中）计算一次的数据，以及该类希望从调用fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()之间存储的数据，直到以后可能从有限元调用转化()等函数。后一类的成员变量被标记为
   * "可变"。
   * 当前的类使用了与MappingQGeneric类基本相同的字段进行存储。因此，它继承自
   * MappingQGeneric::InternalData, 而不是 Mapping::InternalDataBase.
   * 。与 MappingQGeneric::InternalData
   * 的主要区别是，MappingQ会根据我们所在的单元格在 $Q_1$
   * 和 $Q_p$
   * 映射之间切换，所以内部数据对象也需要存储一个指向与
   * $Q_1$ 映射有关的InternalData对象的指针。
   *
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 构造函数。
     *
     */
    InternalData();


    /**
     * 返回这个对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 由<tt>fill_fe_[[sub]face]_values</tt>函数设置的标志。
     * 如果这个标志是 @p true
     * ，我们就在一个内部单元上，并且使用 @p
     * mapping_q1_data。
     *
     */
    mutable bool use_mapping_q1_on_current_cell;

    /**
     * 一个指向结构的指针，用于存储纯 $Q_1$
     * 映射的信息，默认情况下，在所有内部单元上使用。
     *
     */
    std::unique_ptr<typename MappingQGeneric<dim, spacedim>::InternalData>
      mapping_q1_data;

    /**
     * 一个指针，用于存储完整的 $Q_p$
     * 映射的信息，默认情况下，该映射用于所有边界单元。
     *
     */
    std::unique_ptr<typename MappingQGeneric<dim, spacedim>::InternalData>
      mapping_qp_data;
  };

protected:
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
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  using Mapping<dim, spacedim>::fill_fe_face_values;

  // documentation can be found in Mapping::fill_fe_face_values()
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  // documentation can be found in Mapping::fill_fe_subface_values()
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase &   internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  /**
   * @}
   *
   */

protected:
  /**
   * 将在域的边界的所有单元上使用的单元的多项式程度，如果指定的话，也可以在任何地方使用。
   *
   */
  const unsigned int polynomial_degree;

  /**
   * 如果这个标志设置为 @p true ，那么 @p MappingQ
   * 将用于所有单元，而不仅仅是用于边界单元。
   *
   */
  const bool use_mapping_q_on_all_cells;

  /**
   * 指向一个Q1映射的指针。除非在调用构造函数时设置了use_mapping_q_on_all_cells，否则该映射将用于内部单元。该映射也用于transform_real_to_unit_cell()中的任何单元，以便在我们使用完整的映射进行更昂贵的牛顿迭代之前，对点的位置进行一个廉价的初始猜测。
   * @note
   * MappingQEulerian将这个指针重设为MappingQ1Eulerian类型的对象，以确保Q1映射也知道欧拉位移的适当移位和变换。这也意味着我们确实需要在这里存储我们自己的Q1映射，而不是简单地求助于
   * StaticMappingQ1::mapping.  。
   * @note
   * 如果当前对象使用的多项式程度是1，那么qp_mapping和q1_mapping变量就会指向同一个底层对象。
   *
   */
  std::shared_ptr<const MappingQGeneric<dim, spacedim>> q1_mapping;

  /**
   * 指向一个Q_p映射的指针。除非在调用构造函数时设置了use_mapping_q_on_all_cells（在这种情况下，它用于所有单元），否则该映射将用于边界单元。
   * @note
   * MappingQEulerian和MappingC1将这个指针重置为他们自己实现的对象，以确保Q_p映射也知道欧拉位移的适当移动和转换（欧拉情况）和支持点的适当选择（C1情况）。
   * @note
   * 如果用于当前对象的多项式程度为1，那么qp_mapping和q1_mapping变量就会指向同一个底层对象。
   *
   */
  std::shared_ptr<const MappingQGeneric<dim, spacedim>> qp_mapping;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


