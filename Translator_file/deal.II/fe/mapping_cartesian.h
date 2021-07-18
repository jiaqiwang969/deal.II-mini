//include/deal.II-translator/fe/mapping_cartesian_0.txt
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

#ifndef dealii_mapping_cartesian_h
#define dealii_mapping_cartesian_h


#include <deal.II/base/config.h>

#include <deal.II/base/qprojector.h>

#include <deal.II/fe/mapping.h>

#include <cmath>


DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup mapping */ 
 /*@{*/ 

/**
 * 一个提供从参考单元到单元的映射的类，这些单元是轴平行的，即具有矩形（在2d中）或盒子（在3d中）的形状，边缘平行于坐标方向。因此，该类提供的功能等同于例如MappingQ为这类单元格提供的功能。然而，对单元格形状的了解使得这个类的效率大大提升。
 * 具体来说，该映射是针对那些从参考坐标到实际单元的映射是沿坐标方向的缩放的单元。每个单元上从参考坐标
 * $\hat {\mathbf x}$ 到实数坐标 $\mathbf x$ 的变换形式为
 *
 * @f{align*}{
 * {\mathbf x}(\hat {\mathbf x})
 * =
 * \begin{pmatrix}
 *   h_x & 0 \\
 *   0 & h_y
 * \end{pmatrix}
 * \hat{\mathbf x}
 * + {\mathbf v}_0
 * @f}
 * 在2d中，和
 * @f{align*}{
 * {\mathbf x}(\hat {\mathbf x})
 * =
 * \begin{pmatrix}
 *   h_x & 0 & 0 \\
 *   0 & h_y & 0 \\
 *   0 & 0 & h_z
 * \end{pmatrix}
 * \hat{\mathbf x}
 * + {\mathbf v}_0
 * @f}
 * 在3D中， ${\mathbf v}_0$ 是左下角的顶点， $h_x,h_y,h_z$
 * 是单元格沿轴线的延伸。
 * 这个类是为了提高效率，它没有做大量的错误检查。如果你将这种映射应用于不符合上述要求的单元格，你会得到奇怪的结果。
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingCartesian : public Mapping<dim, spacedim>
{
public:
  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * 返回 @p true 因为MappingCartesian保留了顶点位置。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

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
   * @}
   *
   */


private:
  /**
   * @name  与FEValues的接口  
     * @{ 
   *
   */

  /**
   * 存储映射的内部数据。见 Mapping::InternalDataBase
   * 的广泛描述。
   * 这包括创建对象时计算一次的数据（在get_data()中），以及该类希望从调用fill_fe_values()、fill_fe_face_values()或fill_fe_subface_values()之间存储的数据，直到以后可能从有限元调用转化()等函数。后一类的成员变量被标记为
   * "可变"。
   *
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 构造函数。
     *
     */
    InternalData(const Quadrature<dim> &quadrature);

    /**
     * 返回这个对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 我们在坐标方向上看到的最后一个单元的延伸，即<i>h<sub>x</sub></i>,
     * <i>h<sub>y</sub></i>, <i>h<sub>z</sub></i>。
     *
     */
    mutable Tensor<1, dim> cell_extents;

    /**
     * 体积元素
     *
     */
    mutable double volume_element;

    /**
     * 所有正交点的矢量。特别是，所有面的所有点。
     *
     */
    std::vector<Point<dim>> quadrature_points;
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

  /**
   * 用传入的单元格的尺寸更新传入的InternalData对象的cell_extents字段。
   *
   */
  void
  update_cell_extents(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const InternalData &                                        data) const;

  /**
   * 如果传入的InternalData对象的UpdateFlags显示应该更新正交点，则计算这些正交点。
   * 从fill_fe_values调用。
   *
   */
  void
  maybe_update_cell_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const InternalData &                                        data,
    std::vector<Point<dim>> &quadrature_points) const;

  /**
   * 如果传入的InternalData对象的UpdateFlags说它们应该被更新，则计算正交点。
   * 从fill_fe_face_values调用。
   *
   */
  void
  maybe_update_face_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const InternalData &                                        data,
    std::vector<Point<dim>> &quadrature_points) const;

  /**
   * 计算正交点，如果传入的InternalData对象的UpdateFlags说它们应该被更新。
   * 从fill_fe_subface_values调用。
   *
   */
  void
  maybe_update_subface_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const InternalData &                                        data,
    std::vector<Point<dim>> &quadrature_points) const;

  /**
   * 将InternalData中的正交点转换为实空间，在每个方向上用cell_extends缩放单位坐标。
   * 从各种 maybe_update_*_quadrature_points 函数中调用。
   *
   */
  void
  transform_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const InternalData &                                        data,
    const typename QProjector<dim>::DataSetDescriptor &         offset,
    std::vector<Point<dim>> &quadrature_points) const;

  /**
   * 如果传入的InternalData对象的UpdateFlags说它们应该被更新，则计算法向量。
   *
   */
  void
  maybe_update_normal_vectors(
    const unsigned int           face_no,
    const InternalData &         data,
    std::vector<Tensor<1, dim>> &normal_vectors) const;

  /**
   * 由于这个映射的Jacobian是常数，所以Jacobian的所有导数都是相同的零。如果相应的更新标志说它们应该被更新，那么就把这些量填成零。
   *
   */
  void
  maybe_update_jacobian_derivatives(
    const InternalData &             data,
    const CellSimilarity::Similarity cell_similarity,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


