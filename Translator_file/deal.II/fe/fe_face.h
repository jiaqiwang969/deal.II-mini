//include/deal.II-translator/fe/fe_face_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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

#ifndef dealii_fe_face_h
#define dealii_fe_face_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_poly_face.h>

DEAL_II_NAMESPACE_OPEN


/**
 * 一个有限元，在每个面上是张量积多项式，在单元内部未定义。面上的基函数是基于(dim-1)维高斯-洛巴托正交规则的支持点的拉格朗日多项式。对于1度和2度的元素，这些多项式因此对应于通常的等距点上的拉格朗日多项式。
 * 虽然名字没有透露，但该元素在单元面相交的位置是不连续的。特别是，这个有限元是FE_RaviartThomas在面的跟踪空间，在混合方法中发挥作用，例如与FE_DGQ元结合使用。它的使用在
 * step-51 教程程序中进行了演示。
 *
 *
 * @note
 * 由于这个元素只定义在面，所以只有FEFaceValues和FESubfaceValues会提供有用的信息。另一方面，如果你将这个元素与FEValues一起用于单元格积分，那么形状函数的值和导数将有无效的值，不可能产生任何有用的东西。为了使这个元素作为FES系统的一部分的使用更加简单，使用（单元）FEValues对象不会直接失败，但是组合元素中对应于FE_FaceQ的形状函数的那些分量将有上述的无效值。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_FaceQ
  : public FE_PolyFace<TensorProductPolynomials<dim - 1>, dim, spacedim>
{
public:
  /**
   * 度数为<tt>p</tt>的张量积多项式的构造器。使用这个构造函数创建的形状函数对应于每个坐标方向上的拉格朗日多项式。
   *
   */
  FE_FaceQ(const unsigned int p);

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。这个类返回<tt>FE_FaceQ<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>被适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * 在FiniteElement类中实现相应的函数。
   * 由于当前元素是插值的，所以节点值正好是支持点的值。此外，由于当前元素是标量的，支持点的值需要是长度为1的向量。
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  /**
   * 返回从一个元素的一个面插值到邻近元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * 如果形状函数 @p shape_index
   * 在面的某处有非零的函数值，这个函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * @name  支持hp的函数  
     * @{ 
   *
   */

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表， @p fe_other,
   * 是对代表该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是相等的，这两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   * 只有在dim==1的情况下，这样的约束集才是非空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线的自由度。
   * 只有在dim==2的情况下，这种约束的集合才是非空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理四边形上的自由度。
   * 只有在dim==3的情况下，这种约束的集合才是非空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * 返回该元素是否以新的方式实现其悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * @}
   *
   */

  /**
   * 返回一个元素的常量模式列表。对于这个元素，它只是简单地返回一行，所有条目都设置为真。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

private:
  /**
   * 返回每个顶点、线、四边形、六边形的道夫向量。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg);
};



/**
 * FE_FaceQ对一维的特殊化。在这种情况下，有限元只包括一个单元的两个面（=顶点）中的一个自由度，而不考虑程度。然而，这个元素仍然在其构造函数中接受一个度，同时也返回这个度。这样一来，在一维中也可以用痕量元素进行不依赖维度的编程（尽管在一维中根本没有计算上的好处）。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int spacedim>
class FE_FaceQ<1, spacedim> : public FiniteElement<1, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  FE_FaceQ(const unsigned int p);

  virtual std::unique_ptr<FiniteElement<1, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_FaceQ<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>由适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * 返回从一个元素的面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<1, spacedim> &source,
                                FullMatrix<double> &              matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<1, spacedim> &source,
    const unsigned int                subface,
    FullMatrix<double> &              matrix,
    const unsigned int                face_no = 0) const override;

  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 返回该元素是否以新的方式实现其悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs中的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是相等的，这两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   * 只有在dim==1的情况下，这样的约束集才是非空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<1, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线的自由度。
   * 只有在dim==2的情况下，这种约束的集合才是非空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<1, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理四边形上的自由度。
   * 只有在dim==3的情况下，这种约束的集合才是非空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<1, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * 返回一个元素的常数模式列表。对于这个元素，它只是返回一行，所有条目都设置为真。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

protected:
  /*注意：以下函数的定义被内联到类声明中，因为我们在MS Visual Studio中会遇到编译器错误。 
*
*/


  virtual std::unique_ptr<typename FiniteElement<1, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags  /*update_flags*/ ,
    const Mapping<1, spacedim> &  /*mapping*/ ,
    const Quadrature<1> &  /*quadrature*/ ,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &  /*output_data*/ ) const override
  {
    return std::make_unique<
      typename FiniteElement<1, spacedim>::InternalDataBase>();
  }

  using FiniteElement<1, spacedim>::get_face_data;

  std::unique_ptr<typename FiniteElement<1, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags update_flags,
    const Mapping<1, spacedim> &  /*mapping*/ ,
    const hp::QCollection<0> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &  /*output_data*/ ) const override
  {
    AssertDimension(quadrature.size(), 1);

    // generate a new data object and initialize some fields
    auto data_ptr =
      std::make_unique<typename FiniteElement<1, spacedim>::InternalDataBase>();
    data_ptr->update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature[0].size();
    AssertDimension(n_q_points, 1);
    (void)n_q_points;

    // No derivatives of this element are implemented.
    if (data_ptr->update_each & update_gradients ||
        data_ptr->update_each & update_hessians)
      {
        Assert(false, ExcNotImplemented());
      }

    return data_ptr;
  }

  std::unique_ptr<typename FiniteElement<1, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags           update_flags,
    const Mapping<1, spacedim> &mapping,
    const Quadrature<0> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override
  {
    return get_face_data(update_flags,
                         mapping,
                         hp::QCollection<0>(quadrature),
                         output_data);
  }

  virtual void
  fill_fe_values(
    const typename Triangulation<1, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                          cell_similarity,
    const Quadrature<1> &                                     quadrature,
    const Mapping<1, spacedim> &                              mapping,
    const typename Mapping<1, spacedim>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<1,
                                                                       spacedim>
      &                                                          mapping_data,
    const typename FiniteElement<1, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<1, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<1, spacedim>::cell_iterator &cell,
    const unsigned int                                        face_no,
    const hp::QCollection<0> &                                quadrature,
    const Mapping<1, spacedim> &                              mapping,
    const typename Mapping<1, spacedim>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<1,
                                                                       spacedim>
      &                                                          mapping_data,
    const typename FiniteElement<1, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<1, spacedim>::cell_iterator &cell,
    const unsigned int                                        face_no,
    const unsigned int                                        sub_no,
    const Quadrature<0> &                                     quadrature,
    const Mapping<1, spacedim> &                              mapping,
    const typename Mapping<1, spacedim>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<1,
                                                                       spacedim>
      &                                                          mapping_data,
    const typename FiniteElement<1, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                       spacedim>
      &output_data) const override;

private:
  /**
   * 返回每个顶点、线条、四边形、六边形的道夫向量。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg);
};



/**
 * 一个有限元，它是每个面的完全多项式的Legendre元（即相当于FE_DGP在单元上的面），在单元的内部未定义。面上的基函数来自
 * Polynomials::Legendre. 。
 * 虽然名字没有说明，但该元素在单元格面相交的位置是不连续的。该元素适用于混合方法，例如与FE_DGP元素结合使用。在
 * step-51 教程程序中可以找到一个混合方法的例子。
 *
 *
 * @note
 * 由于这个元素只定义在面，所以只有FEFaceValues和FESubfaceValues会提供有用的信息。另一方面，如果你将这个元素与FEValues一起用于单元格积分，那么形状函数的值和导数将有无效的值，不可能产生任何有用的东西。为了使这个元素作为FES系统的一部分的使用更加简单，使用（单元）FEValues对象不会直接失败，但是组合元素中对应FE_FaceP的形状函数的那些分量将有上述的无效值。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_FaceP : public FE_PolyFace<PolynomialSpace<dim - 1>, dim, spacedim>
{
public:
  /**
   * 度数为<tt>p</tt>的多项式的完整基础的构造函数。使用这个构造函数创建的形状函数对应于每个坐标方向的Legendre多项式。
   *
   */
  FE_FaceP(unsigned int p);

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_FaceP<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>由适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 返回这个元素是否以新的方式实现其悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * 返回一个元素的常量模式的列表。对于这个元素，每个面的第一个条目是真的，其他的都是假的（因为常数函数是由Legendre多项式的第一基函数表示的）。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

private:
  /**
   * 返回每个顶点、直线、四边形、六边形的道夫向量。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg);
};



/**
 * 一维的FE_FaceP，即在元素顶点上有自由度。更多信息请参见通用模板的文档。
 *
 *
 */
template <int spacedim>
class FE_FaceP<1, spacedim> : public FE_FaceQ<1, spacedim>
{
public:
  /**
   * 构造函数。
   *
   */
  FE_FaceP(const unsigned int p);

  /**
   * 返回元素的名称
   *
   */
  std::string
  get_name() const override;
};


DEAL_II_NAMESPACE_CLOSE

#endif


