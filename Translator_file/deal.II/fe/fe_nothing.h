//include/deal.II-translator/fe/fe_nothing_0.txt
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

#ifndef dealii_fe_nothing_h
#define dealii_fe_nothing_h

#include <deal.II/base/config.h>

#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 定义一个自由度为零的有限元空间，因此，只能表示一个单一的函数：零函数。
 * 这个类很有用（在hp方法的背景下），可以代表三角结构中不应该分配自由度的空单元，或者描述一个由零扩展到域中不需要的部分的场。因此，一个三角剖分可以分为两个区域：一个是使用正常元素的活动区域，另一个是使用FE_Nothing元素的非活动区域。因此，DoFHandler将不给FE_Nothing单元分配自由度，这个子区域因此被隐式地从计算中删除。
 * step-10 和 step-46 展示了这个元素的使用案例。在论文  @cite
 * Cangiani2012  中也介绍了这个元素的一个有趣的应用。
 *
 *  <h3>FE_Nothing as seen as a function space</h3>
 * 有限元通常最好被解释为形成一个[函数空间](https://en.wikipedia.org/wiki/Function_space)，即形成一个[矢量空间]的函数集合(https://en.wikipedia.org/wiki/Vector_space)。我们确实可以从这个角度来解释FE_Nothing：它对应于函数空间
 * $V_h=\{0\}$
 * ，即到处为零的函数集合。(构造函数可以接受一个参数，如果该参数大于1，则将空间扩展为一个有多个分量的向量值函数的空间，所有的分量都等于零)。事实上，这是一个矢量空间，因为矢量空间中元素的每一个线性组合也是矢量空间中的一个元素，正如单元素零的每一个倍数一样。很明显，函数空间没有自由度，因此该类的名称是。
 *
 *  <h3>FE_Nothing in combination with other elements</h3> 在诸如 step-46
 * 的情况下，人们在对解变量不感兴趣的单元上使用FE_Nothing。例如，在流体结构相互作用问题中，流体速度只定义在域的流体部分的单元上。然后在域的固体部分的单元上使用FE_Nothing来描述速度的有限元空间。换句话说，从概念上讲，速度无处不在，但在域的那些不感兴趣的部分，它是完全为零的，不会占用那里的任何自由度。
 * 问题是，在对解感兴趣的区域（使用 "正常
 * "有限元）和不感兴趣的区域（使用FE_Nothing）之间的界面上会发生什么。该界面上的解应该是零
 *
 * - 即我们考虑一个 "连续的 "有限元场，在使用FE_Nothing的那个区域刚好为零
 *
 * - 或者说，对界面上的连续性没有要求。在deal.II语言中，这是由函数 FiniteElement::compare_for_domination() 返回的内容编码的。如果FE_Nothing "占优势"，那么解在界面上必须为零；如果没有，那么就没有要求，可以认为FE_Nothing是一个函数空间，一般来说是不连续的（即在单元界面上没有任何形式的连续性要求），但在每个单元上都等于零。
 * 一个构造参数表示该元素是否应该被认为是支配性的。默认情况下，它不被视为主导，即FE_Nothing被视为一个不连续的元素。
 *
 *  <h3>FE_Nothing in the context of hanging nodes</h3>
 * 请注意，当引入FE_Nothing元素时，必须注意所得到的网格拓扑结构是否继续有意义。在处理悬挂节点约束时尤其如此，因为该库对这些约束的性质做了一些基本假设。以下的几何形状是可以接受的。
 *
 * @code
 * +---------+----+----+
 * |         | 0  |    |
 * |    1    +----+----+
 * |         | 0  |    |
 * +---------+----+----+
 * @endcode
 *
 *
 * @code
 * +---------+----+----+
 * |         | 1  |    |
 * |    0    +----+----+
 * |         | 1  |    |
 * +---------+----+----+
 * @endcode
 * 这里，0表示一个FE_Nothing单元，1表示一些其他的元素类型。在这些情况下，该库不难计算出必要的悬挂节点约束（即没有约束）。然而，以下的几何形状是不能接受的（至少在目前的实现中）。
 *
 * @code
 * +---------+----+----+
 * |         | 0  |    |
 * |    1    +----+----+
 * |         | 1  |    |
 * +---------+----+----+
 * @endcode
 * 区别在于子面的混合性质，这种情况我们目前还没有实现。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_Nothing : public FiniteElement<dim, spacedim>
{
public:
  /**
   * 构造器。      @param[in]  type 指定参考单元格的类型。
   * @param[in]  n_components
   * 表示赋予这个有限元的向量分量的数量。默认为1。
   * @param[in]  dominate
   * 决定FE_Nothing在compare_for_domination()中是否会支配任何其他FE(默认为`false')。
   * 因此在例如 $Q_1$
   * 与FE_Nothing相遇的界面上，我们将强制这两个函数的轨迹是相同的。因为FE_Nothing编码的空间在任何地方都是零，这意味着
   * $Q_1$
   * 字段在这个界面将被强制变成零。也请看这个类的一般文档中的讨论。
   *
   */
  FE_Nothing(const ReferenceCell &type,
             const unsigned int   n_components = 1,
             const bool           dominate     = false);

  /**
   * 与上述相同，但为超立方体参考单元类型。
   *
   */
  FE_Nothing(const unsigned int n_components = 1, const bool dominate = false);

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。名称为<tt>FE_Nothing  @<dim,spacedim@>(type,  n_components, dominating)</tt>其中<tt>dim</tt>, <tt>spacedim</tt>, <tt>type</tt>, 和<tt>n_components</tt>都是由构造函数或类型签名指定的，有以下例外。    <ol>   <li>  如果<tt>spacedim == dim</tt>，那么该字段不被打印。 </li>   <li>  如果<tt>type</tt>是一个超立方体，那么该字段不被打印。 </li>   <li>  如果<tt>n_components == 1</tt>，则该字段不被打印。 </li>   <li>  如果<tt>dominate == false</tt>那么该字段不被打印。 </li>   </ol>   </ol> .
   *
   */
  virtual std::string
  get_name() const override;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * 返回 @p ith 形状函数在 @p p.  @p p
   * 点的值，该点是参考元素上的一个点。因为当前元素没有自由度，这个函数在实践中显然不应该被调用。
   * 因此，这个函数真正做的是触发一个异常。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

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
   * 准备内部数据结构并填入独立于单元格的值。返回一个指向对象的指针，然后该函数的调用者必须承担该对象的所有权（包括在不再需要时进行销毁）。
   * 在当前情况下，这个函数只是返回一个默认的指针，因为这个元素不存在有意义的数据。
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
      &output_data) const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination()
   * 在当前情况下，如果构造函数 @p dominate
   * 中的第二个参数为真，则该元素被认为是主导的。当这个参数是假的，并且
   * @p fe_other
   * 也是FE_Nothing()的类型，任何一个元素都可以占主导地位。否则就没有_要求了。
   * 也请看这个类的一般文档中的讨论。
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;



  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 返回从给定的有限元到当前有限元的插值矩阵。由于当前的有限元没有自由度，插值矩阵必然是空的。
   *
   */
  virtual void
  get_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source_fe,
    FullMatrix<double> &                interpolation_matrix) const override;

  /**
   * 返回从一个元素的一个面插值到邻近元素的面的矩阵。矩阵的大小是<tt>source.#dofs_per_face</tt>乘以<tt>this->#dofs_per_face</tt>。
   * 由于当前的有限元没有自由度，插值矩阵必然是空的。
   *
   */

  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source_fe,
                                FullMatrix<double> &interpolation_matrix,
                                const unsigned int  face_no = 0) const override;


  /**
   * 返回从一个元素的面插值到邻近元素的子面的矩阵。矩阵的大小是<tt>source.#dofs_per_face</tt>乘以<tt>this->#dofs_per_face</tt>。
   * 由于当前的有限元没有自由度，插值矩阵必然是空的。
   *
   */

  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source_fe,
    const unsigned int                  index,
    FullMatrix<double> &                interpolation_matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * @return  如果FE支配任何其他的，则为真。
   *
   */
  bool
  is_dominating() const;

private:
  /**
   * 如果为真，在compare_for_domination()中，这个元素将支配除它自己之外的任何其他元素。这是因为与其他任何有限元空间相比，只包含零函数的空间肯定更小(因而也更占优势)。
   *
   */
  const bool dominate;
};


 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


