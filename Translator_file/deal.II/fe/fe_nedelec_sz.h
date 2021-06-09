//include/deal.II-translator/fe/fe_nedelec_sz_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

#ifndef dealii_fe_nedelec_sz_h
#define dealii_fe_nedelec_sz_h

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/polynomials_integrated_legendre_sz.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 这个类代表了S.
 * Zaglmayr博士论文中描述的H<sup>curl</sup>-conforming
 * N&eacute;d&eacute;lec元素的实现，<b>High Order Finite Element Methods
 * for Electromagnetic Field Computation</b>，Johannes Kepler
 * Universit&auml;t Linz,
 * 2006。它的使用范围与FE_Nedelec类描述的顶部所述相同。
 * 这个元素克服了传统N&eacute;d&eacute;lec元素中存在的符号冲突问题，这些问题来自于基函数中使用的边缘和面的参数化。因此，这个元素应该为一般的四边形和六面体元素提供一致的结果，对于这些元素，从所有相邻单元看到的边和面的相对方向往往难以确定。
 * 该元素解决符号冲突问题的方法是为局部边和面分配一个全局定义的方向。本地边的方向总是被选择为：定义边的第一个顶点是具有最高全局顶点编号的顶点，而第二个边的顶点是具有最低全局顶点编号的顶点。
 * 同样地，面的方向总是被选择为：第一个顶点被选择为构成面的四个顶点中全局顶点编号最高的那个。然后，第三个顶点被选择为与第一个顶点相对的几何图形，第二和第四个顶点被决定为第二个顶点的全局顶点编号高于第四个顶点。
 * 请注意，这个元素目前不支持不符合要求的网格。
 * 关于这个元素的进一步细节，包括一些基准测试，可以在R.
 * Kynch, P. Ledger: <b>Resolving the sign conflict problem for hp-hexahedral
 * N&eacute;d&eacute;lec elements with application to eddy current
 * problems</b>, Computers & Structures 181, 41-54,
 * 2017（见https://doi.org/10.1016/j.compstruc.2016.05.021）。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_NedelecSZ : public FiniteElement<dim, dim>
{
public:
  static_assert(dim == spacedim,
                "FE_NedelecSZ is only implemented for dim==spacedim!");

  /**
   *
   */
  FE_NedelecSZ(const unsigned int order);

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  /**
   * 这个元素是矢量值，所以这个函数会抛出一个异常。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 没有实现。
   *
   */
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const override;

  /**
   * 这个元素是向量值的，所以这个函数会抛出一个异常。
   *
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 没有实现。
   *
   */
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const override;

  /**
   * 这个元素是向量值的，所以这个函数会抛出一个异常。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 没有实现。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const override;

protected:
  /**
   * 用来将形状函数从参考单元映射到网格单元的映射种类。
   *
   */
  MappingKind mapping_kind;

  virtual std::unique_ptr<
    typename dealii::FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * 计算由第一个参数表示的单元上的形状函数信息。请注意，这个函数必须重新计算与单元格相关的自由度，因此目前还不是线程安全的。
   *
   */
  virtual void
  fill_fe_values(
    const typename Triangulation<dim, dim>::cell_iterator &cell,
    const CellSimilarity::Similarity                       cell_similarity,
    const Quadrature<dim> &                                quadrature,
    const Mapping<dim, dim> &                              mapping,
    const typename Mapping<dim, dim>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, dim>
      &                                                       mapping_data,
    const typename FiniteElement<dim, dim>::InternalDataBase &fedata,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
      &data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  /**
   * 计算由前两个参数表示的单元格和面的形状函数的信息。请注意，这个函数必须重新计算与单元格相关的自由度，因此目前不是线程安全的。
   *
   */
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, dim>::cell_iterator &cell,
    const unsigned int                                     face_no,
    const hp::QCollection<dim - 1> &                       quadrature,
    const Mapping<dim, dim> &                              mapping,
    const typename Mapping<dim, dim>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, dim>
      &                                                       mapping_data,
    const typename FiniteElement<dim, dim>::InternalDataBase &fedata,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
      &data) const override;

  /**
   * 没有实现。
   *
   */
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, dim>::cell_iterator &cell,
    const unsigned int                                     face_no,
    const unsigned int                                     sub_no,
    const Quadrature<dim - 1> &                            quadrature,
    const Mapping<dim, dim> &                              mapping,
    const typename Mapping<dim, dim>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, dim>
      &                                                       mapping_data,
    const typename FiniteElement<dim, dim>::InternalDataBase &fedata,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
      &data) const override;

  /**
   * 衍生的内部数据，用于存储与单元无关的数据。  请注意，由于这个元素的性质，一些有用的预计算量被存储起来，用于计算与细胞相关的形状函数。    被存储的主要数量与边缘和面的参数化相关。这些是。    <ul>   <li>   $\lambda_{i}$ 。
   *
   * - 三线函数，在 $i$ 第1个顶点等于1，在所有其他顶点等于0。 </li>   <li>   $\sigma_{i}$ 。
   *
   * - 与第 $i$ 个顶点相关的线性函数。 </li>   </ul>  这些函数的定义，以及边和面的参数化和边和面的扩展参数，可以在Zaglmayr论文的第82页找到。全局定义的边缘和面的方向的详细定义可以在第67页找到。
   *
   */
  class InternalData : public FiniteElement<dim, dim>::InternalDataBase
  {
  public:
    /**
     * 参考元素上的形状函数的存储。我们只预先计算基于单元的DoFs，因为基于边缘和面的DoFs取决于单元。
     * 由于基于单元的DoFs，这个变量被声明为可变的。
     *
     */
    mutable std::vector<std::vector<Tensor<1, dim>>> shape_values;

    /**
     * 存储参考元素上的形状函数梯度。我们只预先计算基于单元的DoFs，因为基于边缘和面的DoFs取决于单元。
     * 由于基于单元的DoFs，这个变量被声明为可改变的。
     *
     */
    mutable std::vector<std::vector<DerivativeForm<1, dim, dim>>> shape_grads;

    /**
     * 存储顶点之间所有可能的边缘参数化。在计算基于边和面的DoF时需要这些参数，而这些参数是依赖于单元格的。
     * 从顶点i开始到顶点 $j$ 结束的边E的参数化由 $\sigma_{E}
     * = \sigma_{i}
     *
     * - \sigma{j}$ 给出。 sigma_imj_values[q][i][j]存储由顶点 $i$ 和
     * $j$ 在第q个正交点连接的边参数化值。
     * 请注意，并非所有的 $i$ 和 $j$
     * 组合都会在六面体单元上产生有效的边，但它们是以这种方式计算的，以便用于非标准的边和面的方向。
     *
     */
    std::vector<std::vector<std::vector<double>>> sigma_imj_values;

    /**
     * 存储顶点之间所有可能的边缘参数化的梯度。在计算基于边缘和面的DoF时需要这些梯度，这些梯度与细胞有关。请注意，梯度的组成部分是恒定的。
     * 从顶点 $i$ 开始到顶点 $j$ 结束的边的参数化由
     * $\sigma_{E} = \sigma_{i}
     *
     * - \sigma{j}$ 给出。 sigma_imj_grads[i][j][d]存储由顶点 $i$ 和
     * $j$ 连接的边参数化的梯度在组件 $d$ 中。
     * 请注意，边缘参数化的梯度在边缘上是恒定的，所以我们不需要在每个正交点上存储它。
     *
     */
    std::vector<std::vector<std::vector<double>>> sigma_imj_grads;

    /**
     * 储存在正交点的边缘参数化的值。edge_sigma_values[m][q]
     * 存储边m上第q个正交点的边参数化值。这些值随着物理单元的边的方向而变化，因此在用于计算时必须考虑到
     * "符号"。
     *
     */
    std::vector<std::vector<double>> edge_sigma_values;

    /**
     * 储存正交点的边缘参数化梯度。
     * edge_sigma_grads[m][d]存储边m上组件d的边参数化梯度。这些值随物理单元的边的方向变化，因此在用于计算时必须考虑
     * "符号"。
     *
     */
    std::vector<std::vector<double>> edge_sigma_grads;

    /**
     * 储存正交点的边缘扩展参数。这些参数是为12条边存储的，这样全局顶点编号将遵循
     * "标准 "deal.II单元定义的顺序。        从顶点 $i$
     * 开始到顶点 $j$ 结束的一条边的扩展参数由 $\lambda_{E} =
     * \lambda_{i} + \lambda_{j}$ 给出。
     * 请注意，在这个定义下， $\lambda_{E}$
     * 的值不随边的方向变化。 edge_lambda_values[m][q]存储边 $q$
     * 上的第1个正交点的边扩展参数值。
     *
     */
    std::vector<std::vector<double>> edge_lambda_values;

    /**
     * 存储2D中边缘扩展参数的梯度。在这种情况下，它们是常数。edge_lambda_grads_2d[m][d]存储边
     * $d$ 上组件 $m$
     * 的边延伸参数的梯度，以便全局顶点编号*遵循 "标准
     * "deal.II单元定义的顺序。
     *
     */
    std::vector<std::vector<double>> edge_lambda_grads_2d;

    /**
     * 存储3D中边缘扩展参数的梯度。在这种情况下，它们是非恒定的。edge_lambda_grads_3d[m][q][d]
     * 存储了边m上第 $q$ 个正交点的组件 $d$
     * 的边延伸参数的梯度。
     *
     */
    std::vector<std::vector<std::vector<double>>> edge_lambda_grads_3d;

    /**
     * 存储三维边缘扩展参数的二阶导数，这些参数在整个单元中是恒定的。这些存储在12条边上，以便全局顶点编号*遵循
     * "标准
     * "deal.II单元定义的顺序。edge_lambda_gradgrads_3d[m][d1][d2]
     * 存储边上的边延伸参数相对于分量d1和d2的二阶导数
     * $m$  。
     *
     */
    std::vector<std::vector<std::vector<double>>> edge_lambda_gradgrads_3d;

    /**
     * 存储面的扩展参数。这些参数是为6个面而存储的，这样全局顶点编号将遵循
     * "标准 "deal.II单元所定义的顺序。        由顶点v1, v2, v3,
     * v4定义的面的扩展参数F由 $\lambda_{F} = \lambda_{v1} +
     * \lambda_{v2} + \lambda_{v3} + \lambda_{v4}$ 给出。
     * 请注意，在这个定义下， $\lambda_{F}$
     * 的值不会随着面的方向而改变。 face_lambda_values[m][q]
     * 存储面 $q$ 的第1个正交点的面扩展参数值。
     *
     */
    std::vector<std::vector<double>> face_lambda_values;

    /**
     * 存储面扩展参数的梯度。face_lambda_grads[m][d] 存储面 $m$
     * 上 $d$ 组件的面扩展参数的梯度。
     *
     */
    std::vector<std::vector<double>> face_lambda_grads;
  };

private:
  /**
   * 内部函数返回一个 "每个对象的道夫
   * "的向量，返回的向量的分量指的是。  0 = 顶点 1 = 边缘
   * 2 = 面（在2D中是一个单元） 3 = 单元
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * 内部存储所有需要的集成Legendre多项式。
   *
   */
  std::vector<Polynomials::Polynomial<double>> IntegratedLegendrePolynomials;

  /**
   * 内部函数用于填充内部的集成Legendre多项式数组。
   *
   */
  void
  create_polynomials(const unsigned int degree);

  /**
   * 返回基础集中的DoF数量。
   *
   */
  unsigned int
  compute_num_dofs(const unsigned int degree) const;

  /**
   * 在给定的InternalData对象上填充与细胞相关的基于边缘的形状函数。
   *
   */
  void
  fill_edge_values(const typename Triangulation<dim, dim>::cell_iterator &cell,
                   const Quadrature<dim> &quadrature,
                   const InternalData &   fedata) const;

  /**
   * 在给定的InternalData对象上填充依赖细胞的基于面的形状函数。
   *
   */
  void
  fill_face_values(const typename Triangulation<dim, dim>::cell_iterator &cell,
                   const Quadrature<dim> &quadrature,
                   const InternalData &   fedata) const;
};



 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


