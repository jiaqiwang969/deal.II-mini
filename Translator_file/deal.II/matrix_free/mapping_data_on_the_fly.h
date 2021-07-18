//include/deal.II-translator/matrix_free/mapping_data_on_the_fly_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


#ifndef dealii_matrix_free_mapping_data_on_the_fly_h
#define dealii_matrix_free_mapping_data_on_the_fly_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/matrix_free/mapping_info.h>
#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * 这个类使用标准的deal.II信息，以FEEvaluation和朋友们可以使用的形式，提供评估的映射信息，以便进行矢量化访问。由于DoFHandler/Triangulation单元迭代器没有对单元进行矢量化，所以FEEvaluation的矢量化模型的接口是使用同一元素的
     * @p  VectorizedArray::n_array_element
     * 副本。因此，这个接口主要适用于对同一个单元的几个运算符进行评估，例如，在组装单元矩阵时。
     * 与deal.II中的Mapping类不同，这个类实际上并不提供可用于评估几何体的边界描述，而是以FEEvaluation可访问的形式提供来自给定deal.II映射的评估的几何体（如传递给这个类的构造函数）。
     *
     */
    template <int dim, typename Number, typename VectorizedArrayType>
    class MappingDataOnTheFly
    {
      static_assert(
        std::is_same<Number, typename VectorizedArrayType::value_type>::value,
        "Type of Number and of VectorizedArrayType do not match.");

    public:
      /**
       * 构造函数，类似于FEValues。因为这个类只评估几何，所以不需要指定有限元素，最简单的元素FE_Nothing在内部被用于底层的FEValues对象。
       *
       */
      MappingDataOnTheFly(const Mapping<dim> & mapping,
                          const Quadrature<1> &quadrature,
                          const UpdateFlags    update_flags);

      /**
       * 构造函数。这个构造函数与另一个构造函数等价，只是它使对象隐含地使用
       * $Q_1$ 映射（即MappingQGeneric(1)类型的对象）。
       *
       */
      MappingDataOnTheFly(const Quadrature<1> &quadrature,
                          const UpdateFlags    update_flags);

      /**
       * 用给定的单元格迭代器进行初始化。
       *
       */
      void
      reinit(typename dealii::Triangulation<dim>::cell_iterator cell);

      /**
       * 返回reinit()是否被调用过至少一次，也就是说，一个单元格被设置过。
       *
       */
      bool
      is_initialized() const;

      /**
       * 返回当前单元格的三角形迭代器。
       *
       */
      typename dealii::Triangulation<dim>::cell_iterator
      get_cell() const;

      /**
       * 返回一个对底层FEValues对象的引用，该对象评估了某些数量（只有与映射有关的，如Jacobian或映射的正交点可以访问，因为实际上没有使用有限元数据）。
       *
       */
      const dealii::FEValues<dim> &
      get_fe_values() const;

      /**
       * 返回对MappingInfoStorage类型的底层存储域的引用，其格式与MappingInfo中的数据域相同。这确保了与MappingInfo类中预计算的数据字段的兼容性。
       *
       */
      const MappingInfoStorage<dim, dim, Number, VectorizedArrayType> &
      get_data_storage() const;

      /**
       * 返回这个对象底层的一维正交的引用。
       *
       */
      const Quadrature<1> &
      get_quadrature() const;

    private:
      /**
       * 一个单元格迭代器，以备我们在飞行中生成数据，能够检查我们是否需要重新生成存储在这个类中的信息。
       *
       */
      typename dealii::Triangulation<dim>::cell_iterator present_cell;

      /**
       * 为初始化FEValues对象所需的假有限元对象。
       *
       */
      FE_Nothing<dim> fe_dummy;

      /**
       * 一个执行（标量）评估的底层FEValues对象。
       *
       */
      dealii::FEValues<dim> fe_values;

      /**
       * 获取用于重新初始化形状信息的一维正交公式。
       *
       */
      const Quadrature<1> quadrature_1d;

      /**
       * 为单个单元创建的存储部分，并与MappingInfo类比持有。
       *
       */
      MappingInfoStorage<dim, dim, Number, VectorizedArrayType>
        mapping_info_storage;
    };


     /*-------------------------- Inline functions ---------------------------*/ 

    template <int dim, typename Number, typename VectorizedArrayType>
    inline MappingDataOnTheFly<dim, Number, VectorizedArrayType>::
      MappingDataOnTheFly(const Mapping<dim> & mapping,
                          const Quadrature<1> &quadrature,
                          const UpdateFlags    update_flags)
      : fe_values(
          mapping,
          fe_dummy,
          Quadrature<dim>(quadrature),
          MappingInfo<dim, Number, VectorizedArrayType>::compute_update_flags(
            update_flags))
      , quadrature_1d(quadrature)
    {
      mapping_info_storage.descriptor.resize(1);
      mapping_info_storage.descriptor[0].initialize(quadrature);
      mapping_info_storage.data_index_offsets.resize(1);
      mapping_info_storage.JxW_values.resize(fe_values.n_quadrature_points);
      mapping_info_storage.jacobians[0].resize(fe_values.n_quadrature_points);
      if (update_flags & update_quadrature_points)
        {
          mapping_info_storage.quadrature_point_offsets.resize(1, 0);
          mapping_info_storage.quadrature_points.resize(
            fe_values.n_quadrature_points);
        }
      if (fe_values.get_update_flags() & update_normal_vectors)
        {
          mapping_info_storage.normal_vectors.resize(
            fe_values.n_quadrature_points);
          mapping_info_storage.normals_times_jacobians[0].resize(
            fe_values.n_quadrature_points);
        }
      Assert(!(fe_values.get_update_flags() & update_jacobian_grads),
             ExcNotImplemented());
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline MappingDataOnTheFly<dim, Number, VectorizedArrayType>::
      MappingDataOnTheFly(const Quadrature<1> &quadrature,
                          const UpdateFlags    update_flags)
      : MappingDataOnTheFly(::dealii::StaticMappingQ1<dim, dim>::mapping,
                            quadrature,
                            update_flags)
    {}



    template <int dim, typename Number, typename VectorizedArrayType>
    inline void
    MappingDataOnTheFly<dim, Number, VectorizedArrayType>::reinit(
      typename dealii::Triangulation<dim>::cell_iterator cell)
    {
      if (present_cell == cell)
        return;
      present_cell = cell;
      fe_values.reinit(present_cell);
      for (unsigned int q = 0; q < fe_values.get_quadrature().size(); ++q)
        {
          if (fe_values.get_update_flags() & update_JxW_values)
            mapping_info_storage.JxW_values[q] = fe_values.JxW(q);
          if (fe_values.get_update_flags() & update_jacobians)
            {
              Tensor<2, dim> jac = fe_values.jacobian(q);
              jac                = invert(transpose(jac));
              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  mapping_info_storage.jacobians[0][q][d][e] = jac[d][e];
            }
          if (fe_values.get_update_flags() & update_quadrature_points)
            for (unsigned int d = 0; d < dim; ++d)
              mapping_info_storage.quadrature_points[q][d] =
                fe_values.quadrature_point(q)[d];
          if (fe_values.get_update_flags() & update_normal_vectors)
            {
              for (unsigned int d = 0; d < dim; ++d)
                mapping_info_storage.normal_vectors[q][d] =
                  fe_values.normal_vector(q)[d];
              mapping_info_storage.normals_times_jacobians[0][q] =
                mapping_info_storage.normal_vectors[q] *
                mapping_info_storage.jacobians[0][q];
            }
        }
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline bool
    MappingDataOnTheFly<dim, Number, VectorizedArrayType>::is_initialized()
      const
    {
      return present_cell !=
             typename dealii::Triangulation<dim>::cell_iterator();
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline typename dealii::Triangulation<dim>::cell_iterator
    MappingDataOnTheFly<dim, Number, VectorizedArrayType>::get_cell() const
    {
      return fe_values.get_cell();
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline const dealii::FEValues<dim> &
    MappingDataOnTheFly<dim, Number, VectorizedArrayType>::get_fe_values() const
    {
      return fe_values;
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline const MappingInfoStorage<dim, dim, Number, VectorizedArrayType> &
    MappingDataOnTheFly<dim, Number, VectorizedArrayType>::get_data_storage()
      const
    {
      return mapping_info_storage;
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline const Quadrature<1> &
    MappingDataOnTheFly<dim, Number, VectorizedArrayType>::get_quadrature()
      const
    {
      return quadrature_1d;
    }


  } // end of namespace MatrixFreeFunctions
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif


