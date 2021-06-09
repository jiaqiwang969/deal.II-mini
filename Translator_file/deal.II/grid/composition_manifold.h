//include/deal.II-translator/grid/composition_manifold_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_composition_manifold_h
#define dealii_composition_manifold_h


 /*----------------------------   composition_manifold.h     ------------*/ 

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * CompositionManifold。
 * 取两个ChartManifold对象，并使其组成。CompositionManifold对象是一个ChartManifold，从第一个ChartManifold的图到第二个ChartManifold的嵌入空间。如果第一个ChartManifold是周期性的，那么产生的ChartManifold也是周期性的，具有相同的周期性。第二个ChartManifold的周期性是不允许的，如果第二个Manifold是周期性的，构造函数将抛出一个异常。
 * 这个类只适用于dim <= chartdim <= intermediate_spacedim <=
 * spacedim。如果你试图实例化任何不同的东西，在ChartManifold类中的一个违反这个条件的类将抛出一个异常。
 * 给定ChartManifold F和ChartManifold
 * G，这个类代表F之后的G的组成。 模板参数有以下含义。
 * @tparam  dim 产生的ChartManifold的尺寸  @tparam  spacedim
 * 产生的ChartManifold的空间尺寸  @tparam  chartdim
 * 产生的ChartManifold的图表尺寸  @tparam  ] intermediate_dim
 * 第一个ChartManifold的空间尺寸  @tparam  dim1
 * 第一个ChartManifold的尺寸，它也与第二个ChartManifold的图表尺寸重合
 * @tparam  dim2 第二个ChartManifold的尺寸
 *
 *
 * @ingroup manifold
 *
 */
template <int dim,
          int spacedim         = dim,
          int chartdim         = dim,
          int intermediate_dim = dim,
          int dim1             = dim,
          int dim2             = dim>
class CompositionManifold : public ChartManifold<dim, spacedim, chartdim>
{
public:
  /**
   * 构建两个给定流形的组合。
   *
   */
  CompositionManifold(const ChartManifold<dim1, intermediate_dim, chartdim> &F,
                      const ChartManifold<dim2, spacedim, intermediate_dim> &G);

  /**
   * 制作该流形的克隆。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * 将spacedim中的给定点拉回欧几里得图表维空间。这个函数调用G的pull_back()函数，然后再调用F的pull_back()函数。
   *
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const;


  /**
   * 将chartdim维度的点向前推到一个欧几里得的间隔点。该函数首先调用F的push_forward()，然后调用G的push_forward()。
   *
   */
  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const;

  /**
   * 返回F之后的G的组成的导数。
   *
   */
  virtual DerivativeForm<1, chartdim, spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const;

private:
  /**
   * 第一个ChartManifold。
   *
   */
  SmartPointer<
    const ChartManifold<dim1, intermediate_dim, chartdim>,
    CompositionManifold<dim, spacedim, chartdim, dim1, dim2, intermediate_dim>>
    F;


  /**
   * 第二个ChartManifold。
   *
   */
  SmartPointer<
    const ChartManifold<dim2, spacedim, intermediate_dim>,
    CompositionManifold<dim, spacedim, chartdim, dim1, dim2, intermediate_dim>>
    G;
};


 /*------------------Template Implementations------------------------*/ 

template <int dim,
          int spacedim,
          int chartdim,
          int intermediate_dim,
          int dim1,
          int dim2>
CompositionManifold<dim, spacedim, chartdim, intermediate_dim, dim1, dim2>::
  CompositionManifold(const ChartManifold<dim1, intermediate_dim, chartdim> &F,
                      const ChartManifold<dim2, spacedim, intermediate_dim> &G)
  : ChartManifold<dim, spacedim, chartdim>(F.get_periodicity())
  , F(&F)
  , G(&G)
{
  // We don't know what to do with a periodicity in the second manifold, so
  // throw an assertion if the second manifold is periodic
  Assert(G.get_periodicity().norm() == 0.0,
         ExcMessage("The second manifold cannot be periodic."));
}


template <int dim,
          int spacedim,
          int chartdim,
          int intermediate_dim,
          int dim1,
          int dim2>
std::unique_ptr<Manifold<dim, spacedim>>
CompositionManifold<dim, spacedim, chartdim, intermediate_dim, dim1, dim2>::
  clone() const
{
  return std::make_unique<
    CompositionManifold<dim, spacedim, chartdim, intermediate_dim, dim1, dim2>>(
    *F, *G);
}


template <int dim,
          int spacedim,
          int chartdim,
          int intermediate_dim,
          int dim1,
          int dim2>
Point<chartdim>
CompositionManifold<dim, spacedim, chartdim, intermediate_dim, dim1, dim2>::
  pull_back(const Point<spacedim> &space_point) const
{
  return F->pull_back(G->pull_back(space_point));
}



template <int dim,
          int spacedim,
          int chartdim,
          int intermediate_dim,
          int dim1,
          int dim2>
Point<spacedim>
CompositionManifold<dim, spacedim, chartdim, intermediate_dim, dim1, dim2>::
  push_forward(const Point<chartdim> &chart_point) const
{
  return G->push_forward(F->push_forward(chart_point));
}



template <int dim,
          int spacedim,
          int chartdim,
          int intermediate_dim,
          int dim1,
          int dim2>
DerivativeForm<1, chartdim, spacedim>
CompositionManifold<dim, spacedim, chartdim, intermediate_dim, dim1, dim2>::
  push_forward_gradient(const Point<chartdim> &chart_point) const
{
  DerivativeForm<1, chartdim, intermediate_dim> DF =
    F->push_forward_gradient(chart_point);

  DerivativeForm<1, intermediate_dim, spacedim> DG =
    G->push_forward_gradient(F->push_forward(chart_point));

  DerivativeForm<1, chartdim, spacedim> DF_DG;

  for (unsigned int d = 0; d < spacedim; ++d)
    for (unsigned int c = 0; c < chartdim; ++c)
      for (unsigned int s = 0; s < intermediate_dim; ++s)
        DF_DG[d][c] += DG[d][s] * DF[s][c];

  return DF_DG;
}


DEAL_II_NAMESPACE_CLOSE

#endif


