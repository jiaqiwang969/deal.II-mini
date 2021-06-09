//include/deal.II-translator/boost_adaptors/bounding_box_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_boost_adaptor_bounding_box_h
#define dealii_boost_adaptor_bounding_box_h

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <deal.II/boost_adaptors/point.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


namespace boost
{
  namespace geometry
  {
    namespace traits
    {
      /**
       * 用于 dealii::BoundingBox. 的标签适配器。
       *
       */
      template <int dim, class Number>
      struct tag<dealii::BoundingBox<dim, Number>>
      {
        using type = box_tag;
      };

      /**
       * 用于 dealii::BoundingBox. 的点类型适配器
       *
       */
      template <int dim, class Number>
      struct point_type<dealii::BoundingBox<dim, Number>>
      {
        using type = dealii::Point<dim, Number>;
      };

      /**
       * 访问一个 dealii::BoundingBox. 的左下角的D-th坐标
       *
       */
      template <int dim, class Number, std::size_t D>
      struct indexed_access<dealii::BoundingBox<dim, Number>, min_corner, D>
      {
        /**
         * 获取 dealii::BoundingBox. 左下角的D-th坐标的函数。
         *
         */
        static inline double
        get(dealii::BoundingBox<dim, Number> const &box)
        {
          return box.get_boundary_points().first[D];
        }

        /**
         * 为 dealii::BoundingBox.
         * 的左下角的D-th坐标设置的函数。
         *
         */
        static inline void
        set(dealii::BoundingBox<dim, Number> &box, Number value)
        {
          box.get_boundary_points().first[D] = value;
        }
      };

      /**
       * 访问一个 dealii::BoundingBox. 的右上角的D-th坐标。
       *
       */
      template <int dim, class Number, std::size_t D>
      struct indexed_access<dealii::BoundingBox<dim, Number>, max_corner, D>
      {
        /**
         * 获取 dealii::BoundingBox. 右上角D-th坐标的函数。
         *
         */
        static inline double
        get(dealii::BoundingBox<dim, Number> const &box)
        {
          return box.get_boundary_points().second[D];
        }

        /**
         * 为 dealii::BoundingBox.
         * 的右上角的D-th坐标设置的函数。
         *
         */
        static inline void
        set(dealii::BoundingBox<dim, Number> &box, Number value)
        {
          box.get_boundary_points().second[D] = value;
        }
      };
    } // namespace traits
  }   // namespace geometry
} // namespace boost

#endif


