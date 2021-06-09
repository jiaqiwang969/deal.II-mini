//include/deal.II-translator/fe/fe_wedge_p_0.txt
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

#ifndef dealii_fe_fe_p_wedge_h
#define dealii_fe_fe_p_wedge_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_wedge.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/**
 * FE_WedgeP 和 FE_WedgeDGP 的基类。
 *
 *
 * @note 只为3D实现。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_WedgePoly : public dealii::FE_Poly<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  FE_WedgePoly(const unsigned int                                degree,
               const internal::GenericDoFsPerObject &            dpos,
               const typename FiniteElementData<dim>::Conformity conformity);
};

/**
 * 实现楔形上的标量拉格朗日有限元，得到连续的、程度为
 * $k$ 的分片多项式的有限元空间。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_WedgeP : public FE_WedgePoly<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  FE_WedgeP(const unsigned int degree);

  /**
   * @copydoc   dealii::FiniteElement::clone()
   * dealii::FiniteElement::clone() .
   *
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_WedgeP<dim>(degree)</tt>，
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
   * @copydoc   dealii::FiniteElement::hp_line_dof_identities()
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * @copydoc   dealii::FiniteElement::hp_quad_dof_identities()
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;
};

/**
 * 在楔形上实现标量拉格朗日有限元，得到不连续的、程度为
 * $k$ 的分片多项式的有限元空间。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_WedgeDGP : public FE_WedgePoly<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  FE_WedgeDGP(const unsigned int degree);

  /**
   * @copydoc   dealii::FiniteElement::clone()  构建函数。
   *
   */
  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_WedgeDGP<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 由适当的值代替。
   *
   */
  std::string
  get_name() const override;
};

DEAL_II_NAMESPACE_CLOSE

#endif


