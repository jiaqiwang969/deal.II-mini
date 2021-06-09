//include/deal.II-translator/fe/mapping_c1_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#ifndef dealii_mapping_c1_h
#define dealii_mapping_c1_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup mapping */ 
 /*@{*/ 

/**
 * 使用边界的C1（连续可微）立方映射的映射类。这个类是建立在MappingQ之上的，它简单地确定了边界的立方体映射的不同插值点。MappingQ选择它们是为了对边界进行插值，而这个类选择它们是为了使离散的边界是全局连续可微的。
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingC1 : public MappingQ<dim, spacedim>
{
public:
  /**
   * 构造函数。将固定度 @p 3
   * 传给基类，因为一个立方体映射足以生成边界的连续映射。
   *
   */
  MappingC1();

  /**
   * 返回一个指向当前对象副本的指针。然后，这个副本的调用者就拥有了它的所有权。
   *
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

protected:
  /**
   * 一个派生自MappingQGeneric的类，它为通用映射提供了边界对象的支持点，从而使相应的Q3映射最终成为C1。
   *
   */
  class MappingC1Generic : public MappingQGeneric<dim, spacedim>
  {
  public:
    /**
     * 构造函数。
     *
     */
    MappingC1Generic();

    /**
     * 对于<tt>dim=2,3</tt>。将所有位于边界线上的形状函数的支持点追加到向量中
     * @p a.  位于线上但在顶点的点不包括在内。
     * 这个函数选择各自的点不是为了插值边界（就像基类那样），而是为了使产生的立方体映射是一个连续的映射。
     *
     */
    virtual void
    add_line_support_points(
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim>> &                         a) const override;

    /**
     * 对于<tt>dim=3</tt>。将所有位于边界面（3D中的四边形）上的形状函数的支持点追加到向量中
     * @p a.  位于线上但在顶点上的点不包括在内。
     * 这个函数选择各自的点不是为了插值边界（就像基类那样），而是为了使产生的立方体映射是一个连续的映射。
     *
     */
    virtual void
    add_quad_support_points(
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim>> &                         a) const override;
  };
};

 /*@}*/ 

 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN

template <>
void
MappingC1<1>::MappingC1Generic::add_line_support_points(
  const Triangulation<1>::cell_iterator &,
  std::vector<Point<1>> &) const;
template <>
void
MappingC1<2>::MappingC1Generic::add_line_support_points(
  const Triangulation<2>::cell_iterator &cell,
  std::vector<Point<2>> &                a) const;

template <>
void
MappingC1<1>::MappingC1Generic::add_quad_support_points(
  const Triangulation<1>::cell_iterator &,
  std::vector<Point<1>> &) const;
template <>
void
MappingC1<2>::MappingC1Generic::add_quad_support_points(
  const Triangulation<2>::cell_iterator &,
  std::vector<Point<2>> &) const;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


