//include/deal.II-translator/fe/fe_poly_tensor_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_fe_poly_tensor_h
#define dealii_fe_poly_tensor_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/fe/fe.h>

#include <deal.II/lac/full_matrix.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类提供了一个统一的框架，用于实现基于张量值的多项式空间的FiniteElement类，如PolynomialsBDM和PolynomialsRaviartThomas。在这一点上，它是张量值的FE_Poly类的等同物。
 * 本质上，这个类所要求的是，派生类向它描述一个（矢量值）多项式空间，其中每个多项式都正好有
 * @p dim
 * 个矢量分量。提供这种实现的类都是从TensorPolynomialsBase类派生出来的，这些派生类型中的一个对象需要提供给这个类的构造函数。
 *
 *  <h3>Deriving classes</h3>
 * 这个类不是一个完全实现的FiniteElement类，但是实现了一些基于向量估值的多项式类的向量估值元素的共同特征。这里特别缺少的是关于节点值的拓扑位置的信息（即一个自由度在逻辑上是在顶点、边缘、面，还是在单元的内部
 *
 * --这些信息决定了跨单元界面的相关形状函数的连续性属性），而派生类需要提供这些信息。
 * 同样地，在许多情况下，节点函数取决于网格单元的形状，因为它们评估面的法向或切向分量。为了允许一系列的转换，引入了#mapping_kind这个变量。它需要在派生类的构造函数中设置。
 * 任何派生类都必须决定要使用的多项式空间。
 * 这个多项式空间应该简单地实现为一组矢量值的多项式，如PolynomialsBDM和PolynomialsRaviartThomas。
 * 为了便于实现，多项式空间选择哪种基础对当前的类并不重要。
 *
 * - 正如接下来所描述的，这个类处理了从多项式空间模板参数所选择的基础到我们内部要用于有限元计算的基础的转换。
 *
 *  <h4>Determining the correct basis</h4>
 * 在大多数情况下，描述多项式空间的类所使用的基，
 * $\{\tilde\varphi_j(\hat{\mathbf x})\}$
 * ，并不符合我们想用于有限元描述的基，
 * $\{\varphi_j(\hat{\mathbf x})\}$
 * 。相反，我们需要将有限元形状函数表达为多项式空间所提供的基础的线性组合。
 *
 * @f{align*}{
 * \varphi_j = \sum_k c_{jk} \tilde\varphi_j.
 * @f}
 * 这些展开系数 $c_{jk}$
 * 通常在派生类的构造函数中计算。为了方便，这个类一开始（除非另有告知，见下文），假定形状函数应该正是多项式空间所提供的。在派生类的构造函数中，我们通常会有如下形式的代码
 *
 * @code
 * // Now compute the inverse node matrix, generating the correct
 * // basis functions from the raw ones. For a discussion of what
 * // exactly happens here, see FETools::compute_node_matrix.
 * const FullMatrix<double> M = FETools::compute_node_matrix(*this);
 * this->inverse_node_matrix.reinit(n_dofs, n_dofs);
 * this->inverse_node_matrix.invert(M);
 * // From now on, the shape functions provided by FiniteElement::shape_value
 * // and similar functions will be the correct ones, not
 * // the raw shape functions from the polynomial space anymore.
 * @endcode
 * FETools::compute_node_matrix()
 * 函数更详细地解释了它到底计算什么，以及如何计算；无论如何，结果是
 * @p inverse_node_matrix 现在包含了扩展系数 $c_{jk}$
 * ，而且这块代码现在将矩阵设置为非零大小的事实向当前类的函数表明，当要求形状函数的值或导数时，它应该从那时起使用扩展的基础，
 * $\{\varphi_j(\hat{\mathbf x})\}$ ，而不再是原始，"原始 "基础
 * $\{\tilde\varphi_j(\hat{\mathbf x})\}$  。
 * 为了使这个方案奏效，必须确保在调用
 * FETools::compute_node_matrix() 时， @p inverse_node_matrix
 * 的大小为零；因此，对这个函数的调用不能被内联到最后一行
 *
 * - 调用的结果确实需要存储在临时对象 @p M. 中。
 *
 *  <h4>Setting the transformation</h4>
 * 在大多数情况下，矢量值的基函数从参考单元映射到实际网格单元时必须进行转换。这些转换可以从MappingKind集合中选择，并存储在#mapping_kind中。因此，每个构造函数都应该包含这样一行。
 *
 * @code
 * this->mapping_kind = {mapping_none};
 * @endcode
 * （在不需要映射的情况下）或使用MappingKind中定义的任何值，以适合你正在实现的元素。如果每个形状函数可以通过不同的映射来实现，那么
 * @p mapping_kind
 * 可以是一个元素数量与形状函数数量相同的向量。
 * @see  TensorPolynomialsBase
 * @ingroup febase
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_PolyTensor : public FiniteElement<dim, spacedim>
{
public:
  /**
   * 构造函数。这个构造函数通过 TensorPolynomialsBase::clone()
   * 函数对多项式对象进行深度拷贝，并存储一个指向该拷贝的指针。因此，调用网站可以简单地传递一个临时对象作为第一个参数。
   *
   */
  FE_PolyTensor(const TensorPolynomialsBase<dim> &polynomials,
                const FiniteElementData<dim> &    fe_data,
                const std::vector<bool> &         restriction_is_additive_flags,
                const std::vector<ComponentMask> &nonzero_components);


  /**
   * 拷贝构造函数。
   *
   */
  FE_PolyTensor(const FE_PolyTensor &fe);

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * 计算形状函数 @p i 在给定正交点 @p p. 的（标量）值。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  // documentation inherited from the base class
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const override;

  /**
   * 计算（标量）形状函数 @p i 在给定正交点的梯度  @p p.
   * 由于该类代表的元素是矢量值，没有这样的标量值，因此该函数抛出了一个异常。
   *
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  // documentation inherited from the base class
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const override;

  /**
   * 计算（标量）形状函数 @p i 在给定正交点的Hessian  @p p.
   * 由于该类所代表的元素是矢量值，没有这样的标量值，因此该函数抛出一个异常。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  // documentation inherited from the base class
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const override;

protected:
  /**
   * 用来将形状函数从参考单元映射到网格单元的映射类型。如果这个向量的长度为1，所有的形状函数都将采用相同的映射。如果向量的大小等于每个单元的有限元度数，那么每个形状函数将根据向量中的相应条目进行映射。
   *
   */
  std::vector<MappingKind> mapping_kind;

  /**
   * 返回一个布尔值，当有限元使用单一映射时为真，当有限元使用多个映射时为假。
   *
   */
  bool
  single_mapping_kind() const;

  /**
   * 对于三维中非标准面的方向，面（四边形）上的道夫必须被置换，以便与正确的形状函数相结合，另外还可以改变符号。给定一个四边形上的局部dof
   * @p index
   * ，如果该面有非标准的面朝向、面朝上或面朝下的旋转，则返回经过处理的形状函数的符号。在二维和一维中，没有必要进行包络，因此在这种情况下它没有任何作用。
   * permutation本身由界面类FiniteElement<dim>中实现的adjust_quad_dof_index_for_face_orientation返回。
   *
   */
  bool
  adjust_quad_dof_sign_for_face_orientation(const unsigned int index,
                                            const unsigned int face_no,
                                            const bool         face_orientation,
                                            const bool         face_flip,
                                            const bool face_rotation) const;

  /**
   *
   */
  std::vector<Table<2, bool>> adjust_quad_dof_sign_for_face_orientation_table;

  /**
   * 返回有限元的MappingKind  @p i  。
   *
   */
  MappingKind
  get_mapping_kind(const unsigned int i) const;

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
    (void)mapping;
    (void)output_data;
    // generate a new data object and
    // initialize some fields
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
          data_ptr   = std::make_unique<InternalData>();
    auto &data       = dynamic_cast<InternalData &>(*data_ptr);
    data.update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature.size();

    // some scratch arrays
    std::vector<Tensor<1, dim>> values(0);
    std::vector<Tensor<2, dim>> grads(0);
    std::vector<Tensor<3, dim>> grad_grads(0);
    std::vector<Tensor<4, dim>> third_derivatives(0);
    std::vector<Tensor<5, dim>> fourth_derivatives(0);

    if (update_flags & (update_values | update_gradients | update_hessians))
      data.dof_sign_change.resize(this->dofs_per_cell);

    // initialize fields only if really
    // necessary. otherwise, don't
    // allocate memory

    const bool update_transformed_shape_values =
      std::any_of(this->mapping_kind.begin(),
                  this->mapping_kind.end(),
                  [](const MappingKind t) { return t != mapping_none; });

    const bool update_transformed_shape_grads =
      std::any_of(this->mapping_kind.begin(),
                  this->mapping_kind.end(),
                  [](const MappingKind t) {
                    return (t == mapping_raviart_thomas || t == mapping_piola ||
                            t == mapping_nedelec || t == mapping_contravariant);
                  });

    const bool update_transformed_shape_hessian_tensors =
      update_transformed_shape_values;

    if (update_flags & update_values)
      {
        values.resize(this->n_dofs_per_cell());
        data.shape_values.reinit(this->n_dofs_per_cell(), n_q_points);
        if (update_transformed_shape_values)
          data.transformed_shape_values.resize(n_q_points);
      }

    if (update_flags & update_gradients)
      {
        grads.resize(this->n_dofs_per_cell());
        data.shape_grads.reinit(this->n_dofs_per_cell(), n_q_points);
        data.transformed_shape_grads.resize(n_q_points);

        if (update_transformed_shape_grads)
          data.untransformed_shape_grads.resize(n_q_points);
      }

    if (update_flags & update_hessians)
      {
        grad_grads.resize(this->n_dofs_per_cell());
        data.shape_grad_grads.reinit(this->n_dofs_per_cell(), n_q_points);
        data.transformed_shape_hessians.resize(n_q_points);
        if (update_transformed_shape_hessian_tensors)
          data.untransformed_shape_hessian_tensors.resize(n_q_points);
      }

    // Compute shape function values
    // and derivatives and hessians on
    // the reference cell.
    // Make sure, that for the
    // node values N_i holds
    // N_i(v_j)=\delta_ij for all basis
    // functions v_j
    if (update_flags & (update_values | update_gradients))
      for (unsigned int k = 0; k < n_q_points; ++k)
        {
          poly_space->evaluate(quadrature.point(k),
                               values,
                               grads,
                               grad_grads,
                               third_derivatives,
                               fourth_derivatives);

          if (update_flags & update_values)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  data.shape_values[i][k] = values[i];
              else
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  {
                    Tensor<1, dim> add_values;
                    for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
                      add_values += inverse_node_matrix(j, i) * values[j];
                    data.shape_values[i][k] = add_values;
                  }
            }

          if (update_flags & update_gradients)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  data.shape_grads[i][k] = grads[i];
              else
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  {
                    Tensor<2, dim> add_grads;
                    for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
                      add_grads += inverse_node_matrix(j, i) * grads[j];
                    data.shape_grads[i][k] = add_grads;
                  }
            }

          if (update_flags & update_hessians)
            {
              if (inverse_node_matrix.n_cols() == 0)
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  data.shape_grad_grads[i][k] = grad_grads[i];
              else
                for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                  {
                    Tensor<3, dim> add_grad_grads;
                    for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
                      add_grad_grads +=
                        inverse_node_matrix(j, i) * grad_grads[j];
                    data.shape_grad_grads[i][k] = add_grad_grads;
                  }
            }
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
   * FE_PolyTensor的独立于细胞的数据字段。在参考单元格上存储形状函数的值和它们的导数，以便以后使用。
   * 所有表格的组织方式是，正交点<i>k</i>的形状函数<i>i</i>的值可以通过索引<i>(i,k)</i>来访问。
   *
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 包含正交点的形状函数值的数组。每个形状函数都有一行，包含每个正交点的值。
     *
     */
    Table<2, Tensor<1, dim>> shape_values;

    /**
     * 以正交点为单位的形状函数梯度数组。每个形状函数有一行，包含每个正交点的值。
     *
     */
    Table<2, DerivativeForm<1, dim, spacedim>> shape_grads;

    /**
     * 以正交点为单位的形状函数豫备数组。每个形状函数有一行，包含每个正交点的值。
     *
     */
    Table<2, DerivativeForm<2, dim, spacedim>> shape_grad_grads;

    /**
     * 用于中间计算的抓取数组
     *
     */
    mutable std::vector<double>              dof_sign_change;
    mutable std::vector<Tensor<1, spacedim>> transformed_shape_values;
    // for shape_gradient computations
    mutable std::vector<Tensor<2, spacedim>> transformed_shape_grads;
    mutable std::vector<Tensor<2, dim>>      untransformed_shape_grads;
    // for shape_hessian computations
    mutable std::vector<Tensor<3, spacedim>> transformed_shape_hessians;
    mutable std::vector<Tensor<3, dim>> untransformed_shape_hessian_tensors;
  };



  /**
   * 传递给构造函数的对象的副本，描述多项式空间。
   *
   */
  const std::unique_ptr<const TensorPolynomialsBase<dim>> poly_space;

  /**
   * 应用于多项式<i>p<sub>j</sub></i>的节点值<i>N<sub>i</sub></i>的矩阵<i>a<sub>ij</sub></i>的倒数。这个矩阵用于将#poly_space中提供的
   * "原始
   * "基础中的多项式转换为参考单元上的节点函数的对偶基础。
   * 这个对象不是由FE_PolyTensor填充的，而是一个派生类允许重组基函数的机会。如果它留空，则使用#poly_space中的基础。
   *
   */
  FullMatrix<double> inverse_node_matrix;

  /**
   * 一个mutex，用来保护对下面的变量的访问。
   *
   */
  mutable std::mutex cache_mutex;

  /**
   * 如果一个形状函数是在一个点上计算的，我们必须计算所有的形状函数来应用#inverse_node_matrix。为了避免过多的开销，我们对点和函数值进行缓存，以便下次评估。
   *
   */
  mutable Point<dim> cached_point;

  /**
   * 调用shape_value_component()后缓存的形状函数值。
   *
   */
  mutable std::vector<Tensor<1, dim>> cached_values;

  /**
   * 调用shape_grad_component()后缓存的形状函数梯度。
   *
   */
  mutable std::vector<Tensor<2, dim>> cached_grads;

  /**
   * 在调用shape_grad_grad_component()后缓存形状函数的二阶导数。
   *
   */
  mutable std::vector<Tensor<3, dim>> cached_grad_grads;
};

DEAL_II_NAMESPACE_CLOSE

#endif


