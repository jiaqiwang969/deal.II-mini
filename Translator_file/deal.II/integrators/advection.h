//include/deal.II-translator/integrators/advection_0.txt
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

#ifndef dealii_integrators_advection_h
#define dealii_integrators_advection_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/meshworker/dof_info.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
  /**
   * @brief Local integrators related to advection along a vector field and
   * 其DG公式
   * 所有的平流算子都依赖于一个平流速度，在下面的公式中用<b>w</b>表示。它在参数列表中被表示为<tt>velocity</tt>。
   * 函数cell_matrix()和upwind_value_matrix()都是以弱形式取方程，也就是说，方向导数是在测试函数上。
   * @ingroup Integrators
   *
   */
  namespace Advection
  {
    /**
     * 沿着<b>w</b>方向的平流以弱形式进行，导数在测试函数\f[
     * m_{ij} = \int_Z u_j\,(\mathbf w \cdot \nabla) v_i
     * \, dx. \f]上。
     * <tt>fe</tt>中的FiniteElement可以是标量或矢量值的。在后一种情况下，平移算子将分别应用于每个分量。
     * @param  M：作为结果得到的平移矩阵  @param  fe:
     * 描述本地试验函数空间的FEValues对象。必须设置#update_values和#update_gradients，以及#update_JxW_values。
     * @param
     * fetest。描述本地测试函数空间的FEValues对象。必须设置#update_values和#update_gradients。
     * @param
     * velocity（速度）。平流速度，一个维度为<tt>dim</tt>的矢量。每个分量可以包含一个长度为1的矢量，在这种情况下，假定速度是恒定的；如果速度不是恒定的，则包含一个与正交点一样多的条目的矢量。
     * @param 因子是一个可选的结果的乘法因子。
     *
     */
    template <int dim>
    void
    cell_matrix(FullMatrix<double> &                        M,
                const FEValuesBase<dim> &                   fe,
                const FEValuesBase<dim> &                   fetest,
                const ArrayView<const std::vector<double>> &velocity,
                const double                                factor = 1.)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int t_dofs       = fetest.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();

      AssertDimension(velocity.size(), dim);
      // If the size of the
      // velocity vectors is one,
      // then do not increment
      // between quadrature points.
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;

      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      AssertDimension(M.n(), n_dofs);
      AssertDimension(M.m(), t_dofs);

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double       dx     = factor * fe.JxW(k);
          const unsigned int vindex = k * v_increment;

          for (unsigned j = 0; j < n_dofs; ++j)
            for (unsigned i = 0; i < t_dofs; ++i)
              for (unsigned int c = 0; c < n_components; ++c)
                {
                  double wgradv =
                    velocity[0][vindex] * fe.shape_grad_component(i, k, c)[0];
                  for (unsigned int d = 1; d < dim; ++d)
                    wgradv +=
                      velocity[d][vindex] * fe.shape_grad_component(i, k, c)[d];
                  M(i, j) -= dx * wgradv * fe.shape_value_component(j, k, c);
                }
        }
    }



    /**
     * Scalar advection residual operator in strong form \f[ r_i = \int_Z
     * (\mathbf w \cdot \nabla)u\, v_i \, dx. \f]
     * （警告）这不是与cell_matrix()一致的残差，而是与它的转置一致。
     *
     */
    template <int dim>
    inline void
    cell_residual(Vector<double> &                            result,
                  const FEValuesBase<dim> &                   fe,
                  const std::vector<Tensor<1, dim>> &         input,
                  const ArrayView<const std::vector<double>> &velocity,
                  double                                      factor = 1.)
    {
      const unsigned int nq     = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      Assert(input.size() == nq, ExcDimensionMismatch(input.size(), nq));
      Assert(result.size() == n_dofs,
             ExcDimensionMismatch(result.size(), n_dofs));

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k = 0; k < nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned i = 0; i < n_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) += dx * input[k][d] * fe.shape_value(i, k) *
                           velocity[d][k * v_increment];
        }
    }



    /**
     * Vector-valued advection residual operator in strong form \f[ r_i =
     * \int_Z \bigl((\mathbf w \cdot \nabla) \mathbf u\bigr) \cdot\mathbf v_i
     * \, dx. \f] （警告
     * 这不是与cell_matrix()一致的残差，而是与它的转置一致。
     *
     */
    template <int dim>
    inline void
    cell_residual(Vector<double> &                                    result,
                  const FEValuesBase<dim> &                           fe,
                  const ArrayView<const std::vector<Tensor<1, dim>>> &input,
                  const ArrayView<const std::vector<double>> &        velocity,
                  double factor = 1.)
    {
      const unsigned int nq     = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_comp = fe.get_fe().n_components();

      AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
      Assert(result.size() == n_dofs,
             ExcDimensionMismatch(result.size(), n_dofs));

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k = 0; k < nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned i = 0; i < n_dofs; ++i)
            for (unsigned int c = 0; c < n_comp; ++c)
              for (unsigned int d = 0; d < dim; ++d)
                result(i) += dx * input[c][k][d] *
                             fe.shape_value_component(i, k, c) *
                             velocity[d][k * v_increment];
        }
    }



    /**
     * 弱形式的标量平流残差算子 \f[ r_i = \int_Z  (\mathbf w \cdot
     * \nabla)v\, u_i \, dx. \f] 。
     *
     */
    template <int dim>
    inline void
    cell_residual(Vector<double> &                            result,
                  const FEValuesBase<dim> &                   fe,
                  const std::vector<double> &                 input,
                  const ArrayView<const std::vector<double>> &velocity,
                  double                                      factor = 1.)
    {
      const unsigned int nq     = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      Assert(input.size() == nq, ExcDimensionMismatch(input.size(), nq));
      Assert(result.size() == n_dofs,
             ExcDimensionMismatch(result.size(), n_dofs));

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k = 0; k < nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned i = 0; i < n_dofs; ++i)
            for (unsigned int d = 0; d < dim; ++d)
              result(i) -= dx * input[k] * fe.shape_grad(i, k)[d] *
                           velocity[d][k * v_increment];
        }
    }



    /**
     * 弱形式的矢量值平流残差算子 \f[ r_i = \int_Z \bigl((\mathbf
     * w \cdot \nabla) \mathbf v\bigr) \cdot\mathbf u_i \, dx. \f] 。
     *
     */
    template <int dim>
    inline void
    cell_residual(Vector<double> &                            result,
                  const FEValuesBase<dim> &                   fe,
                  const ArrayView<const std::vector<double>> &input,
                  const ArrayView<const std::vector<double>> &velocity,
                  double                                      factor = 1.)
    {
      const unsigned int nq     = fe.n_quadrature_points;
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_comp = fe.get_fe().n_components();

      AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
      Assert(result.size() == n_dofs,
             ExcDimensionMismatch(result.size(), n_dofs));

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k = 0; k < nq; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned i = 0; i < n_dofs; ++i)
            for (unsigned int c = 0; c < n_comp; ++c)
              for (unsigned int d = 0; d < dim; ++d)
                result(i) -= dx * input[c][k] *
                             fe.shape_grad_component(i, k, c)[d] *
                             velocity[d][k * v_increment];
        }
    }



    /**
     * 弱平流算子的边界上风通量。这是在外流边界的试验函数值，其他为零： @f[
     * a_{ij} = \int_{\partial\Omega}
     * [\mathbf w\cdot\mathbf n]_+
     * u_i v_j \, ds
     * @f] <tt>速度</tt>以ArrayView形式提供，具有<tt>dim</tt>向量，每个速度分量一个。如果平流速度是常数，每个向量必须只有一个条目，或者每个正交点都有一个条目。        有限元可以有几个分量，在这种情况下，每个分量的平流速度是相同的。
     *
     */
    template <int dim>
    void
    upwind_value_matrix(FullMatrix<double> &                        M,
                        const FEValuesBase<dim> &                   fe,
                        const FEValuesBase<dim> &                   fetest,
                        const ArrayView<const std::vector<double>> &velocity,
                        double                                      factor = 1.)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int t_dofs       = fetest.dofs_per_cell;
      unsigned int       n_components = fe.get_fe().n_components();
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);

          double nv = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            nv += fe.normal_vector(k)[d] * velocity[d][k * v_increment];

          if (nv > 0)
            {
              for (unsigned i = 0; i < t_dofs; ++i)
                for (unsigned j = 0; j < n_dofs; ++j)
                  {
                    if (fe.get_fe().is_primitive())
                      M(i, j) +=
                        dx * nv * fe.shape_value(i, k) * fe.shape_value(j, k);
                    else
                      for (unsigned int c = 0; c < n_components; ++c)
                        M(i, j) += dx * nv *
                                   fetest.shape_value_component(i, k, c) *
                                   fe.shape_value_component(j, k, c);
                  }
            }
        }
    }



    /**
     * 标量情况。弱平流算子的边界处上风通量的残差。这是在流出边界的试验函数值和流入边界的入流边界条件值：@f[
     * a_{ij} = \int_{\partial\Omega}
     * (\mathbf w\cdot\mathbf n)
     * \widehat u v_j \, ds
     * @f]这里，数值通量 $\widehat u$ 是面的上风值，即有限元函数，其值在流出边界的参数`input`中给出。在流入边界，它是参数`data'中的非均质边界值。        <tt>速度</tt>以ArrayView形式提供，有<tt>dim</tt>向量，每个速度分量一个。如果平流速度是常数，每个向量必须只有一个条目，或者每个正交点都有一个条目。        有限元可以有多个分量，在这种情况下，每个分量的平流速度是相同的。
     *
     */
    template <int dim>
    inline void
    upwind_value_residual(Vector<double> &                            result,
                          const FEValuesBase<dim> &                   fe,
                          const std::vector<double> &                 input,
                          const std::vector<double> &                 data,
                          const ArrayView<const std::vector<double>> &velocity,
                          double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(input.size(), fe.n_quadrature_points);
      AssertDimension(data.size(), fe.n_quadrature_points);

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }


      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);

          double nv = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            nv += fe.normal_vector(k)[d] * velocity[d][k * v_increment];

          // Always use the upwind value
          const double val = (nv > 0.) ? input[k] : -data[k];

          for (unsigned i = 0; i < n_dofs; ++i)
            {
              const double v = fe.shape_value(i, k);
              result(i) += dx * nv * val * v;
            }
        }
    }



    /**
     * 矢量值的情况。边界上的上风通量的残差，用于弱平流算子。这是在流出边界的试验函数的值和流入边界的入流边界条件的值：@f[
     * a_{ij} = \int_{\partial\Omega}
     * (\mathbf w\cdot\mathbf n)
     * \widehat u v_j \, ds
     * @f]这里，数值通量 $\widehat u$ 是面的上风值，即有限元函数，其值在流出边界的参数`input`中给出。在流入边界，它是参数`data'中的非均质边界值。        <tt>速度</tt>以ArrayView形式提供，有<tt>dim</tt>向量，每个速度分量一个。如果平流速度是常数，每个向量必须只有一个条目，或者每个正交点都有一个条目。        有限元可以有几个分量，在这种情况下，每个分量的平流速度是相同的。
     *
     */
    template <int dim>
    inline void
    upwind_value_residual(Vector<double> &                            result,
                          const FEValuesBase<dim> &                   fe,
                          const ArrayView<const std::vector<double>> &input,
                          const ArrayView<const std::vector<double>> &data,
                          const ArrayView<const std::vector<double>> &velocity,
                          double factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int n_comp = fe.get_fe().n_components();

      AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
      AssertVectorVectorDimension(data, n_comp, fe.n_quadrature_points);

      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe.n_quadrature_points);
        }


      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);

          double nv = 0.;
          for (unsigned int d = 0; d < dim; ++d)
            nv += fe.normal_vector(k)[d] * velocity[d][k * v_increment];

          std::vector<double> val(n_comp);

          for (unsigned int d = 0; d < n_comp; ++d)
            {
              val[d] = (nv > 0.) ? input[d][k] : -data[d][k];
              for (unsigned i = 0; i < n_dofs; ++i)
                {
                  const double v = fe.shape_value_component(i, k, d);
                  result(i) += dx * nv * val[d] * v;
                }
            }
        }
    }



    /**
     * 弱平流算子的内部上风通量。矩阵条目对应于试验函数的上风值，乘以试验函数的跳变 @f[
     * a_{ij} = \int_F \left|\mathbf w
     * \cdot \mathbf n\right|
     * u^\uparrow
     * (v^\uparrow-v^\downarrow)
     * \,ds
     * @f] <tt>速度</tt>作为ArrayView提供，具有<tt>dim</tt>向量，每个速度分量一个。如果平流速度是常数，每个向量必须只有一个条目，或者每个正交点都有一个条目。        有限元可以有几个分量，在这种情况下，每个分量的平流方式相同。
     *
     */
    template <int dim>
    void
    upwind_value_matrix(FullMatrix<double> &                        M11,
                        FullMatrix<double> &                        M12,
                        FullMatrix<double> &                        M21,
                        FullMatrix<double> &                        M22,
                        const FEValuesBase<dim> &                   fe1,
                        const FEValuesBase<dim> &                   fe2,
                        const FEValuesBase<dim> &                   fetest1,
                        const FEValuesBase<dim> &                   fetest2,
                        const ArrayView<const std::vector<double>> &velocity,
                        const double                                factor = 1.)
    {
      const unsigned int n1 = fe1.dofs_per_cell;
      // Multiply the quadrature point
      // index below with this factor to
      // have simpler data for constant
      // velocities.
      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe1.n_quadrature_points);
        }

      for (unsigned k = 0; k < fe1.n_quadrature_points; ++k)
        {
          double nbeta = fe1.normal_vector(k)[0] * velocity[0][k * v_increment];
          for (unsigned int d = 1; d < dim; ++d)
            nbeta += fe1.normal_vector(k)[d] * velocity[d][k * v_increment];
          const double        dx_nbeta = factor * std::abs(nbeta) * fe1.JxW(k);
          FullMatrix<double> &M1       = nbeta > 0. ? M11 : M22;
          FullMatrix<double> &M2       = nbeta > 0. ? M21 : M12;
          const FEValuesBase<dim> &fe  = nbeta > 0. ? fe1 : fe2;
          const FEValuesBase<dim> &fetest  = nbeta > 0. ? fetest1 : fetest2;
          const FEValuesBase<dim> &fetestn = nbeta > 0. ? fetest2 : fetest1;
          for (unsigned i = 0; i < n1; ++i)
            for (unsigned j = 0; j < n1; ++j)
              {
                if (fe1.get_fe().is_primitive())
                  {
                    M1(i, j) += dx_nbeta * fe.shape_value(j, k) *
                                fetest.shape_value(i, k);
                    M2(i, j) -= dx_nbeta * fe.shape_value(j, k) *
                                fetestn.shape_value(i, k);
                  }
                else
                  {
                    for (unsigned int d = 0; d < fe1.get_fe().n_components();
                         ++d)
                      {
                        M1(i, j) += dx_nbeta *
                                    fe.shape_value_component(j, k, d) *
                                    fetest.shape_value_component(i, k, d);
                        M2(i, j) -= dx_nbeta *
                                    fe.shape_value_component(j, k, d) *
                                    fetestn.shape_value_component(i, k, d);
                      }
                  }
              }
        }
    }



    /**
     * 标量情况。在内部的上风通量为弱平流算子。    矩阵条目对应于试验函数的上风值，乘以试验函数的跳变 @f[
     * a_{ij} = \int_F \left|\mathbf w
     * \cdot \mathbf n\right|
     * u^\uparrow
     * (v^\uparrow-v^\downarrow)
     * \,ds
     * @f] <tt>速度</tt>作为ArrayView提供，具有<tt>dim</tt>向量，每个速度分量一个。如果平流速度是常数，每个向量必须只有一个条目，或者每个正交点都有一个条目。        有限元可以有几个分量，在这种情况下，每个分量的平流方式相同。
     *
     */
    template <int dim>
    void
    upwind_face_residual(Vector<double> &                            result1,
                         Vector<double> &                            result2,
                         const FEValuesBase<dim> &                   fe1,
                         const FEValuesBase<dim> &                   fe2,
                         const std::vector<double> &                 input1,
                         const std::vector<double> &                 input2,
                         const ArrayView<const std::vector<double>> &velocity,
                         const double factor = 1.)
    {
      Assert(fe1.get_fe().n_components() == 1,
             ExcDimensionMismatch(fe1.get_fe().n_components(), 1));
      Assert(fe2.get_fe().n_components() == 1,
             ExcDimensionMismatch(fe2.get_fe().n_components(), 1));

      const unsigned int n1 = fe1.dofs_per_cell;
      // Multiply the quadrature point
      // index below with this factor to
      // have simpler data for constant
      // velocities.
      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe1.n_quadrature_points);
        }

      for (unsigned k = 0; k < fe1.n_quadrature_points; ++k)
        {
          double nbeta = fe1.normal_vector(k)[0] * velocity[0][k * v_increment];
          for (unsigned int d = 1; d < dim; ++d)
            nbeta += fe1.normal_vector(k)[d] * velocity[d][k * v_increment];
          const double dx_nbeta = factor * nbeta * fe1.JxW(k);

          for (unsigned i = 0; i < n1; ++i)
            {
              const double v1 = fe1.shape_value(i, k);
              const double v2 = fe2.shape_value(i, k);
              const double u1 = input1[k];
              const double u2 = input2[k];
              if (nbeta > 0)
                {
                  result1(i) += dx_nbeta * u1 * v1;
                  result2(i) -= dx_nbeta * u1 * v2;
                }
              else
                {
                  result1(i) += dx_nbeta * u2 * v1;
                  result2(i) -= dx_nbeta * u2 * v2;
                }
            }
        }
    }



    /**
     * 矢量值的情况。在内部的上风通量为弱平流算子。矩阵条目对应于试验函数的上风值，乘以试验函数的跳变 @f[
     * a_{ij} = \int_F \left|\mathbf w
     * \cdot \mathbf n\right|
     * u^\uparrow
     * (v^\uparrow-v^\downarrow)
     * \,ds
     * @f] <tt>速度</tt>作为ArrayView提供，具有<tt>dim</tt>向量，每个速度分量一个。如果平流速度是常数，每个向量必须只有一个条目，或者每个正交点都有一个条目。        有限元可以有几个分量，在这种情况下，每个分量的平流方式相同。
     *
     */
    template <int dim>
    void
    upwind_face_residual(Vector<double> &                            result1,
                         Vector<double> &                            result2,
                         const FEValuesBase<dim> &                   fe1,
                         const FEValuesBase<dim> &                   fe2,
                         const ArrayView<const std::vector<double>> &input1,
                         const ArrayView<const std::vector<double>> &input2,
                         const ArrayView<const std::vector<double>> &velocity,
                         const double factor = 1.)
    {
      const unsigned int n_comp = fe1.get_fe().n_components();
      const unsigned int n1     = fe1.dofs_per_cell;
      AssertVectorVectorDimension(input1, n_comp, fe1.n_quadrature_points);
      AssertVectorVectorDimension(input2, n_comp, fe2.n_quadrature_points);

      // Multiply the quadrature point
      // index below with this factor to
      // have simpler data for constant
      // velocities.
      AssertDimension(velocity.size(), dim);
      const unsigned int v_increment = (velocity[0].size() == 1) ? 0 : 1;
      if (v_increment == 1)
        {
          AssertVectorVectorDimension(velocity, dim, fe1.n_quadrature_points);
        }

      for (unsigned k = 0; k < fe1.n_quadrature_points; ++k)
        {
          double nbeta = fe1.normal_vector(k)[0] * velocity[0][k * v_increment];
          for (unsigned int d = 1; d < dim; ++d)
            nbeta += fe1.normal_vector(k)[d] * velocity[d][k * v_increment];
          const double dx_nbeta = factor * nbeta * fe1.JxW(k);

          for (unsigned i = 0; i < n1; ++i)
            for (unsigned int d = 0; d < n_comp; ++d)
              {
                const double v1 = fe1.shape_value_component(i, k, d);
                const double v2 = fe2.shape_value_component(i, k, d);
                const double u1 = input1[d][k];
                const double u2 = input2[d][k];
                if (nbeta > 0)
                  {
                    result1(i) += dx_nbeta * u1 * v1;
                    result2(i) -= dx_nbeta * u1 * v2;
                  }
                else
                  {
                    result1(i) += dx_nbeta * u2 * v1;
                    result2(i) -= dx_nbeta * u2 * v2;
                  }
              }
        }
    }

  } // namespace Advection
} // namespace LocalIntegrators


DEAL_II_NAMESPACE_CLOSE

#endif


