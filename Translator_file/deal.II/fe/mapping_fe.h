//include/deal.II-translator/fe/mapping_fe_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_fe_h
#define dealii_mapping_fe_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_iterator.h>

#include <array>
#include <cmath>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup mapping */ 
 /*@{*/ 


/**
 * 这个类在三角形的所有单元上一致使用用户提供的有限元来实现多项式映射。
 * 如果人们用与离散化相同的FiniteElement初始化这个类，就会得到一个等参量映射。
 * 如果用FE_Q(degree)对象来初始化这个类，那么这个类就等同于MappingQGeneric(degree)。请注意，这里没有增加利用有限元的张量结构的优化。
 *
 *
 * @note  目前，只对张量程度==1和n_components==1的元素实现。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingFE : public Mapping<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  explicit MappingFE(const FiniteElement<dim, spacedim> &fe);

  /**
   * 复制构造函数。
   *
   */
  MappingFE(const MappingFE<dim, spacedim> &mapping);

  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * 返回映射的程度，即传递给构造函数的有限元的程度。
   *
   */
  unsigned int
  get_degree() const;

  // for documentation, see the Mapping base class
  virtual BoundingBox<spacedim>
  get_bounding_box(const typename Triangulation<dim, spacedim>::cell_iterator
                     &cell) const override;


  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const override;

  /**
   * 总是返回 @p true
   * ，因为这个类中的函数的默认实现保留了顶点位置。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

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
   *
   */

  /**
   * @name  与FEValues的接口  
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
     * 构造函数。
     *
     */
    InternalData(const FiniteElement<dim, spacedim> &fe);

    /**
     * 根据给定的参数，初始化对象中与单元格数据相关的成员变量。
     * 该函数还调用compute_shape_function_values()来实际设置与映射形状函数的值和导数相关的成员变量。
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
     *
     */
    void
    compute_shape_function_values(const std::vector<Point<dim>> &unit_points);


    /**
     * 正交点的形状函数。形状函数是按张量积顺序排列的，因此必须对顶点重新排序以获得变换。
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
     * 填充一次。
     *
     */
    std::array<std::vector<Tensor<1, dim>>,
               GeometryInfo<dim>::faces_per_cell *(dim - 1)>
      unit_tangentials;

    /**
     * 底层有限元。
     *
     */
    const FiniteElement<dim, spacedim> &fe;

    /**
     * 映射的多项式程度。
     *
     */
    const unsigned int polynomial_degree;

    /**
     * 形状函数的数量。
     *
     */
    const unsigned int n_shape_functions;

    /**
     * 每个正交点上的协方变换的张量。
     * 存储的矩阵是Jacobian G^{-1}，其中G = Jacobian^{t}
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
     * 在  @p
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

    /**
     * 投射的正交权重。
     *
     */
    mutable std::vector<double> quadrature_weights;
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
   * @}
   *
   */

protected:
  const std::unique_ptr<FiniteElement<dim, spacedim>> fe;

  /**
   * 用作单元格映射的形状函数的多项式的程度。
   *
   */
  const unsigned int polynomial_degree;

  /**
   * 返回映射的支持点的位置。
   *
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

private:
  Table<2, double> mapping_support_point_weights;
};



 /*@}*/ 

 /*----------------------------------------------------------------------*/ 

#ifndef DOXYGEN

template <int dim, int spacedim>
inline const double &
MappingFE<dim, spacedim>::InternalData::shape(const unsigned int qpoint,
                                              const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr, shape_values.size());
  return shape_values[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline double &
MappingFE<dim, spacedim>::InternalData::shape(const unsigned int qpoint,
                                              const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr, shape_values.size());
  return shape_values[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<1, dim> &
MappingFE<dim, spacedim>::InternalData::derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_derivatives.size());
  return shape_derivatives[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline Tensor<1, dim> &
MappingFE<dim, spacedim>::InternalData::derivative(const unsigned int qpoint,
                                                   const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_derivatives.size());
  return shape_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<2, dim> &
MappingFE<dim, spacedim>::InternalData::second_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_second_derivatives.size());
  return shape_second_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<2, dim> &
MappingFE<dim, spacedim>::InternalData::second_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_second_derivatives.size());
  return shape_second_derivatives[qpoint * n_shape_functions + shape_nr];
}

template <int dim, int spacedim>
inline const Tensor<3, dim> &
MappingFE<dim, spacedim>::InternalData::third_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_third_derivatives.size());
  return shape_third_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<3, dim> &
MappingFE<dim, spacedim>::InternalData::third_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_third_derivatives.size());
  return shape_third_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<4, dim> &
MappingFE<dim, spacedim>::InternalData::fourth_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_fourth_derivatives.size());
  return shape_fourth_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<4, dim> &
MappingFE<dim, spacedim>::InternalData::fourth_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_fourth_derivatives.size());
  return shape_fourth_derivatives[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline bool
MappingFE<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}



#endif // DOXYGEN

 /* -------------- declaration of explicit specializations ------------- */ 


DEAL_II_NAMESPACE_CLOSE

#endif


