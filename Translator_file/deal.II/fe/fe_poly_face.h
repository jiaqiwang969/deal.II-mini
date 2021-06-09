//include/deal.II-translator/fe/fe_poly_face_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
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

#ifndef dealii_fe_poly_face_h
#define dealii_fe_poly_face_h


#include <deal.II/base/config.h>

#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup febase */ 
 /*@{*/ 

/**
 * @warning  这个类还没有得到充分的测试!
 * 这个类给出了一个统一的框架来实现只位于网格面上的FiniteElement类。它们是基于多项式空间，如TensorProductPolynomials或PolynomialSpace类。
 * 每个实现以下功能的类都可以作为模板参数PolynomialType。
 *
 *
 * @code
 * double compute_value (const unsigned int i,
 *                     const Point<dim> &p) const;
 * @endcode
 * 示例类是TensorProductPolynomials、PolynomialSpace或PolynomialsP。
 * 这个类不是一个完全实现的FiniteElement类。相反，有几个在FiniteElement类中声明的纯虚拟函数不能由这个类来实现，而是留给派生类来实现。
 *
 *
 */
template <class PolynomialType,
          int dim      = PolynomialType::dimension + 1,
          int spacedim = dim>
class FE_PolyFace : public FiniteElement<dim, spacedim>
{
public:
  /**
   * 构造函数。
   *
   */
  FE_PolyFace(const PolynomialType &        poly_space,
              const FiniteElementData<dim> &fe_data,
              const std::vector<bool> &     restriction_is_additive_flags);

  /**
   * 返回该有限元的多项式程度，即传递给构造函数的值。
   *
   */
  unsigned int
  get_degree() const;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

protected:
  /*注意：以下函数的定义被内联到类的声明中，因为我们在MS Visual Studio中否则会遇到编译器错误。 
*
*/


  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    (void)update_flags;
    (void)mapping;
    (void)quadrature;
    (void)output_data;
    return std::make_unique<InternalData>();
  }

  using FiniteElement<dim, spacedim>::get_face_data;

  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags               update_flags,
    const Mapping<dim, spacedim> &  mapping,
    const hp::QCollection<dim - 1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    (void)mapping;
    (void)output_data;
    AssertDimension(quadrature.size(), 1);

    // generate a new data object and
    // initialize some fields
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
          data_ptr   = std::make_unique<InternalData>();
    auto &data       = dynamic_cast<InternalData &>(*data_ptr);
    data.update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature[0].size();

    // some scratch arrays
    std::vector<double>             values(0);
    std::vector<Tensor<1, dim - 1>> grads(0);
    std::vector<Tensor<2, dim - 1>> grad_grads(0);
    std::vector<Tensor<3, dim - 1>> empty_vector_of_3rd_order_tensors;
    std::vector<Tensor<4, dim - 1>> empty_vector_of_4th_order_tensors;

    // initialize fields only if really
    // necessary. otherwise, don't
    // allocate memory
    if (data.update_each & update_values)
      {
        values.resize(poly_space.n());
        data.shape_values.resize(poly_space.n(),
                                 std::vector<double>(n_q_points));
        for (unsigned int i = 0; i < n_q_points; ++i)
          {
            poly_space.evaluate(quadrature[0].point(i),
                                values,
                                grads,
                                grad_grads,
                                empty_vector_of_3rd_order_tensors,
                                empty_vector_of_4th_order_tensors);

            for (unsigned int k = 0; k < poly_space.n(); ++k)
              data.shape_values[k][i] = values[k];
          }
      }
    // No derivatives of this element
    // are implemented.
    if (data.update_each & update_gradients ||
        data.update_each & update_hessians)
      {
        Assert(false, ExcNotImplemented());
      }

    return data_ptr;
  }

  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    return get_face_data(
      update_flags,
      mapping,
      hp::QCollection<dim - 1>(QProjector<dim - 1>::project_to_all_children(
        ReferenceCells::get_hypercube<dim - 1>(), quadrature)),
      output_data);
  }

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * 独立于细胞的数据字段。
   * 关于这个类的一般用途的信息，请看基类的文档。
   *
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 包含一个面上的正交点的形状函数值的数组。
     * 每个形状函数都有一行，包含每个正交点的值。
     * 在这个数组中，我们将形状函数的值存储在单元格一个面上的正交点。由于这些值在转换到实际单元时不会改变，我们只需要在访问具体单元时将它们复制过来。
     * 特别是，我们可以简单地将同一组数值复制到每个面。
     *
     */
    std::vector<std::vector<double>> shape_values;
  };

  /**
   * 多项式空间。其类型由模板参数PolynomialType给出。
   *
   */
  PolynomialType poly_space;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


