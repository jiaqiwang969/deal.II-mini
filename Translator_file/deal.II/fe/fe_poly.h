//include/deal.II-translator/fe/fe_poly_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
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

#ifndef dealii_fe_poly_h
#define dealii_fe_poly_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/scalar_polynomials_base.h>

#include <deal.II/fe/fe.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup febase */ 
 /*@{*/ 

/**
 * 这个类给出了一个统一的框架，用于实现基于标量多项式空间的FiniteElement类，如TensorProductPolynomials或PolynomialSpace类。这个类在FE_PolyTensor类中有一个对应的张量值有限元的类。
 * 每一个拥有以下公共成员变量和函数的类都可以作为模板参数
 * @p PolynomialType.
 *
 *
 * @code
 * static const unsigned int dimension;
 *
 * void evaluate (const Point<dim>            &unit_point,
 *               std::vector<double>         &values,
 *               std::vector<Tensor<1,dim> > &grads,
 *               std::vector<Tensor<2,dim> > &grad_grads,
 *               std::vector<Tensor<3,dim> > &third_derivatives,
 *               std::vector<Tensor<4,dim> > &fourth_derivatives) const;
 *
 * double compute_value (const unsigned int i,
 *                      const Point<dim> &p) const;
 *
 * template <int order>
 * Tensor<order,dim> compute_derivative (const unsigned int i,
 *                                      const Point<dim> &p) const;
 * @endcode
 * 示例类是TensorProductPolynomials、PolynomialSpace或PolynomialsP。
 * 这个类不是一个完全实现的FiniteElement类。相反，有几个在FiniteElement和FiniteElement类中声明的纯虚拟函数不能被这个类实现，而是留待派生类实现。
 * @todo  由于spacedim !=
 * dim的几乎所有函数都是专用的，这个类需要清理。
 *
 *
 */

template <int dim, int spacedim = dim>
class FE_Poly : public FiniteElement<dim, spacedim>
{
public:
  /**
   * 构造函数。
   *
   */
  FE_Poly(const ScalarPolynomialsBase<dim> &poly_space,
          const FiniteElementData<dim> &    fe_data,
          const std::vector<bool> &         restriction_is_additive_flags,
          const std::vector<ComponentMask> &nonzero_components);

  /**
   * 复制构造函数。
   *
   */
  FE_Poly(const FE_Poly &fe);

  /**
   * 返回该有限元的多项式程度，即传递给构造函数的值。
   *
   */
  unsigned int
  get_degree() const;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * 返回底层多项式空间。
   *
   */
  const ScalarPolynomialsBase<dim> &
  get_poly_space() const;

  /**
   * 返回底层多项式空间的编号与基函数的lexicographic排序相比。返回
   * PolynomialType::get_numbering(). 。
   * @note
   * 这个类的一些实现不支持这个函数，因为对它们来说，基函数的lexicographic排序是不可能的。这方面的例子有。FE_SimplexP,
   * FE_WedgeP, 和 FE_PyramidP.
   *
   */
  std::vector<unsigned int>
  get_poly_space_numbering() const;

  /**
   * 返回底层多项式空间的反编号。返回
   * PolynomialType::get_numbering_inverse().
   * @note  参见get_poly_space_numbering()的说明。
   *
   */
  std::vector<unsigned int>
  get_poly_space_numbering_inverse() const;

  /**
   * 返回在点<tt>p</tt>处的<tt>i</tt>的形状函数的值。关于这个函数的语义，请看FiniteElement基类的更多信息。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回<tt>i</tt>第1个形状函数在点<tt>p</tt>处的<tt>分量</tt>的值。关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素是标量，返回值与调用不带<tt>_component</tt>后缀的函数相同，前提是指定的分量为零。
   *
   */
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const override;

  /**
   * 返回<tt>i</tt>第1个形状函数在点<tt>p</tt>的梯度。关于这个函数的语义，请看FiniteElement基类的更多信息。
   *
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回<tt>i</tt>第1个形状函数在点<tt>p</tt>处的<tt>分量</tt>向量分量的梯度。关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素是标量，返回值与调用不带<tt>_component</tt>后缀的函数相同，前提是指定的分量为零。
   *
   */
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const override;

  /**
   * 返回单元格上<tt>p</tt>点的<tt>i</tt>th形状函数的二阶导数张量。关于这个函数的语义，请看FiniteElement基类的更多信息。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回<tt>i</tt>第1个形状函数的<tt>分量</tt>在<tt>p</tt>点的二阶导数。关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素是标量，返回值与调用不带<tt>_component</tt>后缀的函数相同，前提是指定的分量为零。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const override;

  /**
   * 返回单元格上<tt>p</tt>点的<tt>i</tt>th形状函数的三阶导数的张量。关于这个函数的语义，请看FiniteElement基类的更多信息。
   *
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * 返回<tt>i</tt>第1个形状函数在<tt>p</tt>点的<tt>分量</tt>的3次导数。关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素是标量，返回值与调用不带<tt>_component</tt>后缀的函数相同，前提是指定的分量为零。
   *
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * 返回单元格上<tt>p</tt>点的<tt>i</tt>第4个形状函数的四阶导数的张量。关于这个函数的语义，请看FiniteElement基类的更多信息。
   *
   */
  virtual Tensor<4, dim>
  shape_4th_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * 返回<tt>i</tt>第1个形状函数在<tt>p</tt>点的<tt>分量</tt>向量分量的四阶导数。关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素是标量，返回值与调用不带<tt>_component</tt>后缀的函数相同，前提是指定的分量为零。
   *
   */
  virtual Tensor<4, dim>
  shape_4th_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * 返回这个对象的内存消耗估计值（以字节为单位）。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

protected:
  /*注意：以下函数的定义被内联到类声明中，因为我们在MS Visual Studio中会遇到编译器错误。 
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
    (void)mapping;

    // generate a new data object and
    // initialize some fields
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
          data_ptr   = std::make_unique<InternalData>();
    auto &data       = dynamic_cast<InternalData &>(*data_ptr);
    data.update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature.size();

    // initialize some scratch arrays. we need them for the underlying
    // polynomial to put the values and derivatives of shape functions
    // to put there, depending on what the user requested
    std::vector<double> values(
      update_flags & update_values ? this->n_dofs_per_cell() : 0);
    std::vector<Tensor<1, dim>> grads(
      update_flags & update_gradients ? this->n_dofs_per_cell() : 0);
    std::vector<Tensor<2, dim>> grad_grads(
      update_flags & update_hessians ? this->n_dofs_per_cell() : 0);
    std::vector<Tensor<3, dim>> third_derivatives(
      update_flags & update_3rd_derivatives ? this->n_dofs_per_cell() : 0);
    std::vector<Tensor<4, dim>>
      fourth_derivatives; // won't be needed, so leave empty

    // now also initialize fields the fields of this class's own
    // temporary storage, depending on what we need for the given
    // update flags.
    //
    // there is one exception from the rule: if we are dealing with
    // cells (i.e., if this function is not called via
    // get_(sub)face_data()), then we can already store things in the
    // final location where FEValues::reinit() later wants to see
    // things. we then don't need the intermediate space. we determine
    // whether we are on a cell by asking whether the number of
    // elements in the output array equals the number of quadrature
    // points (yes, it's a cell) or not (because in that case the
    // number of quadrature points we use here equals the number of
    // quadrature points summed over *all* faces or subfaces, whereas
    // the number of output slots equals the number of quadrature
    // points on only *one* face)
    if ((update_flags & update_values) &&
        !((output_data.shape_values.n_rows() > 0) &&
          (output_data.shape_values.n_cols() == n_q_points)))
      data.shape_values.reinit(this->n_dofs_per_cell(), n_q_points);

    if (update_flags & update_gradients)
      data.shape_gradients.reinit(this->n_dofs_per_cell(), n_q_points);

    if (update_flags & update_hessians)
      data.shape_hessians.reinit(this->n_dofs_per_cell(), n_q_points);

    if (update_flags & update_3rd_derivatives)
      data.shape_3rd_derivatives.reinit(this->n_dofs_per_cell(), n_q_points);

    // next already fill those fields of which we have information by
    // now. note that the shape gradients are only those on the unit
    // cell, and need to be transformed when visiting an actual cell
    if (update_flags & (update_values | update_gradients | update_hessians |
                        update_3rd_derivatives))
      for (unsigned int i = 0; i < n_q_points; ++i)
        {
          poly_space->evaluate(quadrature.point(i),
                               values,
                               grads,
                               grad_grads,
                               third_derivatives,
                               fourth_derivatives);

          // the values of shape functions at quadrature points don't change.
          // consequently, write these values right into the output array if
          // we can, i.e., if the output array has the correct size. this is
          // the case on cells. on faces, we already precompute data on *all*
          // faces and subfaces, but we later on copy only a portion of it
          // into the output object; in that case, copy the data from all
          // faces into the scratch object
          if (update_flags & update_values)
            if (output_data.shape_values.n_rows() > 0)
              {
                if (output_data.shape_values.n_cols() == n_q_points)
                  for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
                    output_data.shape_values[k][i] = values[k];
                else
                  for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
                    data.shape_values[k][i] = values[k];
              }

          // for everything else, derivatives need to be transformed,
          // so we write them into our scratch space and only later
          // copy stuff into where FEValues wants it
          if (update_flags & update_gradients)
            for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
              data.shape_gradients[k][i] = grads[k];

          if (update_flags & update_hessians)
            for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
              data.shape_hessians[k][i] = grad_grads[k];

          if (update_flags & update_3rd_derivatives)
            for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
              data.shape_3rd_derivatives[k][i] = third_derivatives[k];
        }
    return data_ptr;
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
     * 带有正交点的形状函数值的数组。每个形状函数都有一行，包含每个正交点的值。
     * 在这个数组中，我们将形状函数的值存储在单元格的正交点上。由于这些值在转换到实际单元时不会改变，我们只需要在访问具体单元时将它们复制过来。
     *
     */
    Table<2, double> shape_values;

    /**
     * 包含正交点的形状函数梯度的数组。每个形状函数都有一行，包含每个正交点的值。
     * 我们将梯度存储在单元格的正交点上。然后我们只需要在访问实际单元格时应用转换（这是一个矩阵-向量乘法）。
     *
     */
    Table<2, Tensor<1, dim>> shape_gradients;

    /**
     * 包含正交点的形状函数豫备数的数组。每个形状函数都有一行，包含每个正交点的值。
     * 我们在单元格的正交点上存储豫备值。然后，我们只需要在访问实际单元格时应用转换。
     *
     */
    Table<2, Tensor<2, dim>> shape_hessians;

    /**
     * 包含正交点的形状函数三阶导数的数组。每个形状函数都有一行，包含每个正交点的值。
     * 我们将三阶导数存储在单元格的正交点上。然后，我们只需要在访问实际单元格时应用转换。
     *
     */
    Table<2, Tensor<3, dim>> shape_3rd_derivatives;
  };

  /**
   * 通过减去对应于Jacobian推动的前向梯度的项来修正形状Hessians。    在修正之前，Hessians将由@f[
   * D_{ijk} = \frac{d^2\phi_i}{d \hat x_J d \hat x_K} (J_{jJ})^{-1}
   * (J_{kK})^{-1},
   * @f]给出，其中 $J_{iI}=\frac{d x_i}{d \hat x_I}$  。在校正之后，正确的黑森斯将由@f[
   * \frac{d^2 \phi_i}{d x_j d x_k} = D_{ijk}
   *
   * - H_{mjk} \frac{d \phi_i}{d x_m}, @f]给出，其中 $H_{ijk}$
   * 是雅各布式推前导数。
   *
   */
  void
  correct_hessians(
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &                mapping_data,
    const unsigned int n_q_points) const;

  /**
   * 通过减去对应于Jacobian推前梯度和第二导数的项来修正形状的第三导数。    在修正之前，第三导数将由@f[
   * D_{ijkl} = \frac{d^3\phi_i}{d \hat x_J d \hat x_K d \hat x_L} (J_{jJ})^{-1}
   * (J_{kK})^{-1} (J_{lL})^{-1},
   * @f]给出，其中 $J_{iI}=\frac{d x_i}{d \hat x_I}$  。修正后，正确的第三导数将由@f[
   * \frac{d^3\phi_i}{d x_j d x_k d x_l} = D_{ijkl}
   *
   * - H_{mjl} \frac{d^2 \phi_i}{d x_k d x_m}
   *
   *
   *
   *
   *
   *
   *
   * - H_{mkl} \frac{d^2 \phi_i}{d x_j d x_m}
   *
   * - H_{mjk} \frac{d^2 \phi_i}{d x_l d x_m}
   *
   *
   *
   *
   *
   *
   *
   * - K_{mjkl} \frac{d \phi_i}{d x_m}, @f]给出，其中 $H_{ijk}$
   * 是雅各布式推前导数， $K_{ijkl}$
   * 是雅各布式推前二导数。
   *
   */
  void
  correct_third_derivatives(
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &                mapping_data,
    const unsigned int n_q_points) const;


  /**
   * 多项式空间。
   *
   */
  const std::unique_ptr<ScalarPolynomialsBase<dim>> poly_space;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


