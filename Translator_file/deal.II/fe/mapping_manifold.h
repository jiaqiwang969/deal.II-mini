//include/deal.II-translator/fe/mapping_manifold_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_manifold_h
#define dealii_mapping_manifold_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/mapping.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int, int>
class MappingQ;


 /*!@addtogroup mapping */ 
 /*@{*/ 


/**
 * 这个类实现了符合Manifold的映射的功能。这个映射通过利用来自底层Manifold对象的几何信息来计算参考单元和真实单元之间的转换。
 * 使用该映射计算的正交点位于准确的几何对象上，使用该类计算的正切和法向量是与底层几何体相切和相法。这与MappingQ类不同，后者使用某个阶次的多项式来逼近几何体，然后使用逼近的表面来计算法线和切线。
 * @warning
 * 由于数学上的原因，我们不可能在一个由SphericalManifold描述的几何体上使用这个类：更多信息见该类的注释。
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingManifold : public Mapping<dim, spacedim>
{
public:
  /**
   * 构造函数。
   *
   */
  MappingManifold() = default;

  /**
   * 复制构造函数。
   *
   */
  MappingManifold(const MappingManifold<dim, spacedim> &mapping);

  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * 总是返回 @p true
   * ，因为这个类假定顶点总是位于底层的漫游上。
   *
   */
  virtual bool
  preserves_vertex_locations() const override;

  virtual bool
  is_compatible_with(const ReferenceCell &cell_type) const override;

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
   * @name  将张量从参考坐标转换为实数坐标的函数  @{
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

  /**
   * @name  与FEValues的接口  
     * @{ 
   *
   */

public:
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
    InternalData() = default;

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
     * 计算与Manifold对象相关的权重，在计算正交点的位置时需要传递这些权重。
     *
     */
    void
    compute_manifold_quadrature_weights(const Quadrature<dim> &quadrature);

    /**
     * 内部存储顶点。
     *
     */
    void
    store_vertices(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

    /**
     * 返回这个对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 当前的单元格顶点。        计算每个。
     *
     */
    mutable std::vector<Point<spacedim>> vertices;

    /**
     * 当前的单元格。        计算每一个。
     *
     */
    mutable typename Triangulation<dim, spacedim>::cell_iterator cell;

    /**
     * 参考单元上的实际正交。        计算一次。
     *
     */
    Quadrature<dim> quad;


    /**
     * 用于流形正交公式的正交权重的值。
     * 流形类有一个函数 (Manifold::get_new_point())
     * ，根据流形上一些周围点的加权平均数返回新点。对于每个正交点，我们用一个正交公式来调用这个函数，这个公式是用当前单元格的顶点和在正交点本身评估的FE_Q(1)有限元的基函数值构建的。虽然每个单元的顶点都在变化，但每个正交点的权重可以计算一次。我们将这些信息存储在以下变量中，其中第一个索引贯穿正交点，第二个索引贯穿顶点索引。
     * 计算一次。
     *
     */
    std::vector<std::vector<double>> cell_manifold_quadrature_weights;


    /**
     * 用于 Manifold::get_new_point(). 的权重向量
     * 对于每个点（单元格的内部），我们计算每个顶点对这个点的权重。如果该点位于一个顶点，那么这个顶点的权重为1，其他所有顶点的权重为0。如果该点位于一个单元格的内部，那么每个顶点的权重只是在该点评估的与每个顶点相关的
     * $d$ -线性形状函数。        这个数组的大小为
     * GeometryInfo<dim>::vertices_per_cell,
     * ，但它不能被转换成一个固定大小的数组，因为它被用作
     * Manifold::get_new_point() 的输入，该数组希望看到一个
     * std::vector<double> 的权重。
     *
     */
    mutable std::vector<double> vertex_weights;

    /**
     * 单位切线向量。用于计算边界形式和法向量。
     * 这个数组有`(dim-1)  GeometryInfo::faces_per_cell` 项。第一组
     * GeometryInfo::faces_per_cell
     * 包含每个面的第一个切向的向量；第二组
     * GeometryInfo<dim>::faces_per_cell
     * 条目包含第二个切向的向量（仅在三维中，因为每个面有两个切向），等等。
     * 填充一次。
     *
     */
    std::array<std::vector<Tensor<1, dim>>,
               GeometryInfo<dim>::faces_per_cell *(dim - 1)>
      unit_tangentials;

    /**
     * 每个正交点的协方差变换的张量。
     * 存储的矩阵是Jacobian G^{-1}，其中G=Jacobian^{t}
     * Jacobian，是地图的第一基本形式；如果dim=spacedim，则还原为Jacobian矩阵的转置，其本身存储在该结构的
     * @p contravariant 域。        在每个单元格上计算。
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
     * 在每个正交点的雅各布系数的行列式。如果#update_volume_elements就会被填满。
     *
     */
    mutable std::vector<double> volume_elements;

    /**
     * 一个指向使用中的Manifold的指针。        每次更新。
     *
     */
    mutable SmartPointer<const Manifold<dim, spacedim>> manifold;
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
};



 /*@}*/ 

 /*----------------------------------------------------------------------*/ 

#ifndef DOXYGEN

template <int dim, int spacedim>
inline void
MappingManifold<dim, spacedim>::InternalData::store_vertices(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  vertices.resize(GeometryInfo<dim>::vertices_per_cell);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    vertices[i] = cell->vertex(i);
  this->cell = cell;
}


template <int dim, int spacedim>
inline void
MappingManifold<dim, spacedim>::InternalData::
  compute_manifold_quadrature_weights(const Quadrature<dim> &quad)
{
  cell_manifold_quadrature_weights.resize(
    quad.size(), std::vector<double>(GeometryInfo<dim>::vertices_per_cell));
  for (unsigned int q = 0; q < quad.size(); ++q)
    {
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        {
          cell_manifold_quadrature_weights[q][i] =
            GeometryInfo<dim>::d_linear_shape_function(quad.point(q), i);
        }
    }
}



template <int dim, int spacedim>
inline bool
MappingManifold<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}


template <int dim, int spacedim>
bool
MappingManifold<dim, spacedim>::is_compatible_with(
  const ReferenceCell &cell_type) const
{
  if (cell_type.get_dimension() != dim)
    return false; // TODO: or is this an error?

  if (cell_type.is_hyper_cube())
    return true;

  return false;
}



#endif // DOXYGEN

 /* -------------- declaration of explicit specializations ------------- */ 


DEAL_II_NAMESPACE_CLOSE

#endif


