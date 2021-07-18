//include/deal.II-translator/fe/fe_simplex_p_0.txt
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

#ifndef dealii_fe_fe_p_h
#define dealii_fe_fe_p_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/**
 * FE_SimplexP和FE_SimplexDGP的基类。
 *
 *
 * @note  只在2D和3D中实现。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_SimplexPoly : public dealii::FE_Poly<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  FE_SimplexPoly(const unsigned int                                degree,
                 const std::vector<unsigned int> &                 dpo_vector,
                 const typename FiniteElementData<dim>::Conformity conformity);

  /**
   * 返回一个元素的恒定模式的列表。对于这个元素，该列表由所有组件的真实参数组成。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * @copydoc   dealii::FiniteElement::get_prolongation_matrix() 。
   * @note  只对 RefinementCase::isotropic_refinement. 实施。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * @copydoc   dealii::FiniteElement::get_restriction_matrix()  *
   * dealii::FiniteElement::get_restriction_matrix() 。
   * @note  只对 RefinementCase::isotropic_refinement. 实施
   *
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * @copydoc   dealii::FiniteElement::get_face_interpolation_matrix()
   *
   */
  void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source_fe,
                                FullMatrix<double> &interpolation_matrix,
                                const unsigned int  face_no) const override;



  /**
   * @copydoc   dealii::FiniteElement::hp_constraints_are_implemented() *
   */
  void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &x_source_fe,
    const unsigned int                  subface,
    FullMatrix<double> &                interpolation_matrix,
    const unsigned int                  face_no) const override;

  /**
   *
   */
  bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc
   * dealii::FiniteElement::convert_generalized_support_point_values_to_dof_values()
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  mutable Threads::Mutex mutex;
};



/**
 * 实现标量拉格朗日有限元 $P_k$ ，得到连续的、程度为 $k$
 * 的分片多项式的有限元空间。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_SimplexP : public FE_SimplexPoly<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  FE_SimplexP(const unsigned int degree);

  /**
   * @copydoc   dealii::FiniteElement::clone()
   *
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_SimplexP<dim>(degree)</tt>，
   * @p dim 和 @p degree 用适当的值替换。
   *
   */
  std::string
  get_name() const override;

  /**
   * @copydoc   dealii::FiniteElement::compare_for_domination() 。
   *
   */
  FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim) const override;

  /**
   * @copydoc   dealii::FiniteElement::hp_vertex_dof_identities()
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;


  /**
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;
};



/**
 * 实现标量不连续拉格朗日有限元 $P_k$ ，有时表示为
 * $P_{-k}$ ，得到不连续的、程度为 $k$
 * 的分片多项式的有限元空间。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_SimplexDGP : public FE_SimplexPoly<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  FE_SimplexDGP(const unsigned int degree);

  /**
   * @copydoc   dealii::FiniteElement::clone()  构建函数。
   *
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_SimplexDGP<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 被适当的值替换。
   *
   */
  std::string
  get_name() const override;

  /**
   * @copydoc   dealii::FiniteElement::compare_for_domination() 。
   *
   */
  FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim) const override;

  /**
   * @copydoc   dealii::FiniteElement::hp_vertex_dof_identities()
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc   dealii::FiniteElement::hp_line_dof_identities()
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;
};

DEAL_II_NAMESPACE_CLOSE

#endif


