//include/deal.II-translator/integrators/divergence_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_integrators_divergence_h
#define dealii_integrators_divergence_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/integrators/grad_div.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/meshworker/dof_info.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
  /**
   * @brief Local integrators related to the divergence operator and its
   * 追踪。
   * @ingroup Integrators
   *
   */
  namespace Divergence
  {
    /**
     * 用于发散的单元矩阵。导数是在试验函数上。    \f[
     * \int_Z v\nabla \cdot \mathbf u \,dx \f]
     * 这是强发散算子，试验空间至少应该是<b>H</b><sup>div</sup>。试验函数可以是不连续的。
     *
     */
    template <int dim>
    void
    cell_matrix(FullMatrix<double> &     M,
                const FEValuesBase<dim> &fe,
                const FEValuesBase<dim> &fetest,
                double                   factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int i = 0; i < t_dofs; ++i)
            {
              const double vv = fetest.shape_value(i, k);
              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int j = 0; j < n_dofs; ++j)
                  {
                    const double du = fe.shape_grad_component(j, k, d)[d];
                    M(i, j) += dx * du * vv;
                  }
            }
        }
    }

    /**
     * 强形式的发散算子的残差。\f[ \int_Z v\nabla \cdot \mathbf u
     * \,dx
     * \f]这是强发散算子，试验空间至少应该是<b>H</b><sup>div</sup>。试验函数可能是不连续的。
     * 函数cell_matrix()是该函数相对于测试函数的Frechet导数。
     *
     */
    template <int dim, typename number>
    void
    cell_residual(Vector<number> &                                    result,
                  const FEValuesBase<dim> &                           fetest,
                  const ArrayView<const std::vector<Tensor<1, dim>>> &input,
                  const double factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertVectorVectorDimension(input, dim, fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert(result.size() == t_dofs,
             ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k = 0; k < fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) += dx * input[d][k][d] * fetest.shape_value(i, k);
        }
    }


    /**
     * 弱形式的发散算子的残差。\f[
     *
     * - \int_Z
     * \nabla v \cdot \mathbf u \,dx
     * \f]这是弱发散算子，试验空间至少应是<b>H</b><sup>1</sup>。试验函数可能是不连续的。
     * @todo
     * 验证：函数cell_matrix()是该函数相对于试验函数的Frechet导数。
     *
     */
    template <int dim, typename number>
    void
    cell_residual(Vector<number> &                            result,
                  const FEValuesBase<dim> &                   fetest,
                  const ArrayView<const std::vector<double>> &input,
                  const double                                factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertVectorVectorDimension(input, dim, fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert(result.size() == t_dofs,
             ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k = 0; k < fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) -= dx * input[d][k] * fetest.shape_grad(i, k)[d];
        }
    }


    /**
     * 梯度的单元格矩阵。该导数是关于试验函数的。\f[
     * \int_Z \nabla u \cdot \mathbf v\,dx
     * \f]这是强梯度，试验空间至少要在<i>H</i><sup>1</sup>。试验函数可以是不连续的。
     *
     */
    template <int dim>
    void
    gradient_matrix(FullMatrix<double> &     M,
                    const FEValuesBase<dim> &fe,
                    const FEValuesBase<dim> &fetest,
                    double                   factor = 1.)
    {
      const unsigned int t_dofs = fetest.dofs_per_cell;
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertDimension(fe.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int i = 0; i < t_dofs; ++i)
              {
                const double vv = fetest.shape_value_component(i, k, d);
                for (unsigned int j = 0; j < n_dofs; ++j)
                  {
                    const Tensor<1, dim> &Du = fe.shape_grad(j, k);
                    M(i, j) += dx * vv * Du[d];
                  }
              }
        }
    }

    /**
     * 强形式的梯度算子的残差。\f[ \int_Z \mathbf v\cdot\nabla u
     * \,dx
     * \f]这是强梯度算子，试验空间至少应该是<b>H</b><sup>1</sup>。试验函数可能是不连续的。
     * 函数gradient_matrix()是该函数相对于试验函数的Frechet导数。
     *
     */
    template <int dim, typename number>
    void
    gradient_residual(Vector<number> &                   result,
                      const FEValuesBase<dim> &          fetest,
                      const std::vector<Tensor<1, dim>> &input,
                      const double                       factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertDimension(input.size(), fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert(result.size() == t_dofs,
             ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k = 0; k < fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) +=
                dx * input[k][d] * fetest.shape_value_component(i, k, d);
        }
    }

    /**
     * 梯度算子的残差为弱形式。\f[
     *
     * -\int_Z
     * \nabla\cdot \mathbf v u \,dx
     * \f]这是弱梯度算子，试验空间至少应该是<b>H</b><sup>div</sup>。试验函数可能是不连续的。
     * @todo
     * 验证：函数gradient_matrix()是该函数相对于试验函数的Frechet导数。
     *
     */
    template <int dim, typename number>
    void
    gradient_residual(Vector<number> &           result,
                      const FEValuesBase<dim> &  fetest,
                      const std::vector<double> &input,
                      const double               factor = 1.)
    {
      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertDimension(input.size(), fetest.n_quadrature_points);
      const unsigned int t_dofs = fetest.dofs_per_cell;
      Assert(result.size() == t_dofs,
             ExcDimensionMismatch(result.size(), t_dofs));

      for (unsigned int k = 0; k < fetest.n_quadrature_points; ++k)
        {
          const double dx = factor * fetest.JxW(k);

          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) -=
                dx * input[k] * fetest.shape_grad_component(i, k, d)[d];
        }
    }

    /**
     * 发散算子的踪迹，即矢量值试验空间和测试空间的法向分量的乘积。    @f[ \int_F (\mathbf u\cdot \mathbf n) v \,ds @f]
     *
     */
    template <int dim>
    void
    u_dot_n_matrix(FullMatrix<double> &     M,
                   const FEValuesBase<dim> &fe,
                   const FEValuesBase<dim> &fetest,
                   double                   factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const Tensor<1, dim> ndx = factor * fe.JxW(k) * fe.normal_vector(k);
          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              for (unsigned int d = 0; d < dim; ++d)
                M(i, j) += ndx[d] * fe.shape_value_component(j, k, d) *
                           fetest.shape_value(i, k);
        }
    }

    /**
     * 发散算子的轨迹，即向量值试验空间的法线分量与试验空间的乘积。    @f[
     * \int_F (\mathbf u\cdot \mathbf n) v \,ds
     * @f]
     *
     */
    template <int dim, typename number>
    void
    u_dot_n_residual(Vector<number> &                            result,
                     const FEValuesBase<dim> &                   fe,
                     const FEValuesBase<dim> &                   fetest,
                     const ArrayView<const std::vector<double>> &data,
                     double                                      factor = 1.)
    {
      const unsigned int t_dofs = fetest.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(fetest.get_fe().n_components(), 1);
      AssertDimension(result.size(), t_dofs);
      AssertVectorVectorDimension(data, dim, fe.n_quadrature_points);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const Tensor<1, dim> ndx = factor * fe.normal_vector(k) * fe.JxW(k);

          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) += ndx[d] * fetest.shape_value(i, k) * data[d][k];
        }
    }

    /**
     * 梯度算子的轨迹，即矢量值试验空间的法线分量与试验空间的乘积。    @f[
     * \int_F u (\mathbf v\cdot \mathbf n) \,ds
     * @f]
     *
     */
    template <int dim, typename number>
    void
    u_times_n_residual(Vector<number> &           result,
                       const FEValuesBase<dim> &  fetest,
                       const std::vector<double> &data,
                       double                     factor = 1.)
    {
      const unsigned int t_dofs = fetest.dofs_per_cell;

      AssertDimension(fetest.get_fe().n_components(), dim);
      AssertDimension(result.size(), t_dofs);
      AssertDimension(data.size(), fetest.n_quadrature_points);

      for (unsigned int k = 0; k < fetest.n_quadrature_points; ++k)
        {
          const Tensor<1, dim> ndx =
            factor * fetest.normal_vector(k) * fetest.JxW(k);

          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) +=
                ndx[d] * fetest.shape_value_component(i, k, d) * data[k];
        }
    }

    /**
     * 发散算子的踪迹，即矢量值试验函数的法线分量的跳跃与试验函数的平均值的乘积。    @f[
     * \int_F (\mathbf u_1\cdot \mathbf n_1 + \mathbf u_2 \cdot \mathbf n_2)
     * \frac{v_1+v_2}{2} \,ds
     * @f]
     *
     */
    template <int dim>
    void
    u_dot_n_matrix(FullMatrix<double> &     M11,
                   FullMatrix<double> &     M12,
                   FullMatrix<double> &     M21,
                   FullMatrix<double> &     M22,
                   const FEValuesBase<dim> &fe1,
                   const FEValuesBase<dim> &fe2,
                   const FEValuesBase<dim> &fetest1,
                   const FEValuesBase<dim> &fetest2,
                   double                   factor = 1.)
    {
      const unsigned int n_dofs = fe1.dofs_per_cell;
      const unsigned int t_dofs = fetest1.dofs_per_cell;

      AssertDimension(fe1.get_fe().n_components(), dim);
      AssertDimension(fe2.get_fe().n_components(), dim);
      AssertDimension(fetest1.get_fe().n_components(), 1);
      AssertDimension(fetest2.get_fe().n_components(), 1);
      AssertDimension(M11.m(), t_dofs);
      AssertDimension(M11.n(), n_dofs);
      AssertDimension(M12.m(), t_dofs);
      AssertDimension(M12.n(), n_dofs);
      AssertDimension(M21.m(), t_dofs);
      AssertDimension(M21.n(), n_dofs);
      AssertDimension(M22.m(), t_dofs);
      AssertDimension(M22.n(), n_dofs);

      for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
        {
          const double dx = factor * fe1.JxW(k);
          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              for (unsigned int d = 0; d < dim; ++d)
                {
                  const double un1 = fe1.shape_value_component(j, k, d) *
                                     fe1.normal_vector(k)[d];
                  const double un2 = -fe2.shape_value_component(j, k, d) *
                                     fe1.normal_vector(k)[d];
                  const double v1 = fetest1.shape_value(i, k);
                  const double v2 = fetest2.shape_value(i, k);

                  M11(i, j) += .5 * dx * un1 * v1;
                  M12(i, j) += .5 * dx * un2 * v1;
                  M21(i, j) += .5 * dx * un1 * v2;
                  M22(i, j) += .5 * dx * un2 * v2;
                }
        }
    }

    /**
     * 正态分量的跳变 @f[
     * \int_F
     * (\mathbf u_1\cdot \mathbf n_1 + \mathbf u_2 \cdot \mathbf n_2)
     * (\mathbf v_1\cdot \mathbf n_1 + \mathbf v_2 \cdot \mathbf n_2)
     * \,ds
     * @f] 。
     *
     */
    template <int dim>
    void
    u_dot_n_jump_matrix(FullMatrix<double> &     M11,
                        FullMatrix<double> &     M12,
                        FullMatrix<double> &     M21,
                        FullMatrix<double> &     M22,
                        const FEValuesBase<dim> &fe1,
                        const FEValuesBase<dim> &fe2,
                        double                   factor = 1.)
    {
      const unsigned int n_dofs = fe1.dofs_per_cell;

      AssertDimension(fe1.get_fe().n_components(), dim);
      AssertDimension(fe2.get_fe().n_components(), dim);
      AssertDimension(M11.m(), n_dofs);
      AssertDimension(M11.n(), n_dofs);
      AssertDimension(M12.m(), n_dofs);
      AssertDimension(M12.n(), n_dofs);
      AssertDimension(M21.m(), n_dofs);
      AssertDimension(M21.n(), n_dofs);
      AssertDimension(M22.m(), n_dofs);
      AssertDimension(M22.n(), n_dofs);

      for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
        {
          const double dx = factor * fe1.JxW(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              for (unsigned int d = 0; d < dim; ++d)
                {
                  const double un1 = fe1.shape_value_component(j, k, d) *
                                     fe1.normal_vector(k)[d];
                  const double un2 = -fe2.shape_value_component(j, k, d) *
                                     fe1.normal_vector(k)[d];
                  const double vn1 = fe1.shape_value_component(i, k, d) *
                                     fe1.normal_vector(k)[d];
                  const double vn2 = -fe2.shape_value_component(i, k, d) *
                                     fe1.normal_vector(k)[d];

                  M11(i, j) += dx * un1 * vn1;
                  M12(i, j) += dx * un2 * vn1;
                  M21(i, j) += dx * un1 * vn2;
                  M22(i, j) += dx * un2 * vn2;
                }
        }
    }

    /**
     * 由FEValuesBase对象确定的正交集上的发散的<i>L</i><sup>2</sup>-norm。
     * 该向量预计由长度等于正交点数量的dim向量组成。有限元的分量数量必须等于空间维度。
     *
     */
    template <int dim>
    double
    norm(const FEValuesBase<dim> &                           fe,
         const ArrayView<const std::vector<Tensor<1, dim>>> &Du)
    {
      AssertDimension(fe.get_fe().n_components(), dim);
      AssertVectorVectorDimension(Du, dim, fe.n_quadrature_points);

      double result = 0;
      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          double div = Du[0][k][0];
          for (unsigned int d = 1; d < dim; ++d)
            div += Du[d][k][d];
          result += div * div * fe.JxW(k);
        }
      return result;
    }

  } // namespace Divergence
} // namespace LocalIntegrators


DEAL_II_NAMESPACE_CLOSE

#endif


