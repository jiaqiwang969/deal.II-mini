//include/deal.II-translator/integrators/maxwell_0.txt
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

#ifndef dealii_integrators_maxwell_h
#define dealii_integrators_maxwell_h


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
   * @brief Local integrators related to curl operators and their traces.
   * 我们对卷积算子采用以下约定。首先，在三个空间维度上@f[
   * \nabla\times \mathbf u = \begin{pmatrix} \partial_2 u_3
   *
   * - \partial_3 u_2 \\ \partial_3 u_1
   *
   * - \partial_1 u_3 \\ \partial_1 u_2
   *
   * - \partial_2 u_1 \end{pmatrix}.
   * @f] 在两个空间维度上，通过将一个向量<b>u</b>扩展到 $(u_1, u_2, 0)^T$ 和一个标量<i>p</i>扩展到 $(0,0,p)^T$ 得到curl。  计算非零分量，我们得到一个向量函数的标量卷度和一个标量函数的向量卷度。目前的实现交换了符号，我们有。  @f[
   * \nabla \times \mathbf u = \partial_1 u_2
   *
   * - \partial_2 u_1, \qquad \nabla \times p = \begin{pmatrix} \partial_2 p
   * \\
   *
   * -\partial_1 p \end{pmatrix} @f]
   * @ingroup Integrators
   *
   */
  namespace Maxwell
  {
    /**
     * 辅助函数。给出<tt>dim</tt>二阶导数的张量，计算向量函数的curl的curl。在二维和三维中的结果是。    @f[
     * \nabla\times\nabla\times \mathbf u = \begin{pmatrix}
     * \partial_1\partial_2 u_2
     *
     * - \partial_2^2 u_1 \\
     * \partial_1\partial_2 u_1
     *
     * - \partial_1^2 u_2
     * \end{pmatrix}
     *
     * \nabla\times\nabla\times \mathbf u = \begin{pmatrix}
     * \partial_1\partial_2 u_2 + \partial_1\partial_3 u_3
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * - (\partial_2^2+\partial_3^2) u_1 \\
     * \partial_2\partial_3 u_3 + \partial_2\partial_1 u_1
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * - (\partial_3^2+\partial_1^2) u_2 \\
     * \partial_3\partial_1 u_1 + \partial_3\partial_2 u_2
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * - (\partial_1^2+\partial_2^2) u_3
     * \end{pmatrix}
     * @f]
     * @note
     * 第三个张量参数在二维中不使用，例如，可以重复前面的一个。
     *
     */
    template <int dim>
    Tensor<1, dim>
    curl_curl(const Tensor<2, dim> &h0,
              const Tensor<2, dim> &h1,
              const Tensor<2, dim> &h2)
    {
      Tensor<1, dim> result;
      switch (dim)
        {
          case 2:
            result[0] = h1[0][1] - h0[1][1];
            result[1] = h0[0][1] - h1[0][0];
            break;
          case 3:
            result[0] = h1[0][1] + h2[0][2] - h0[1][1] - h0[2][2];
            result[1] = h2[1][2] + h0[1][0] - h1[2][2] - h1[0][0];
            result[2] = h0[2][0] + h1[2][1] - h2[0][0] - h2[1][1];
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
      return result;
    }

    /**
     * 辅助函数。给出第一导数的<tt>dim</tt>张量和法向量，计算切向卷曲 @f[
     * \mathbf n \times \nabla \times u.
     * @f] 。
     * @note
     * 第三个张量参数在二维中不使用，例如，可以重复前面的一个。
     *
     */
    template <int dim>
    Tensor<1, dim>
    tangential_curl(const Tensor<1, dim> &g0,
                    const Tensor<1, dim> &g1,
                    const Tensor<1, dim> &g2,
                    const Tensor<1, dim> &normal)
    {
      Tensor<1, dim> result;

      switch (dim)
        {
          case 2:
            result[0] = normal[1] * (g1[0] - g0[1]);
            result[1] = -normal[0] * (g1[0] - g0[1]);
            break;
          case 3:
            result[0] =
              normal[2] * (g2[1] - g0[2]) + normal[1] * (g1[0] - g0[1]);
            result[1] =
              normal[0] * (g0[2] - g1[0]) + normal[2] * (g2[1] - g1[2]);
            result[2] =
              normal[1] * (g1[0] - g2[1]) + normal[0] * (g0[2] - g2[0]);
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
      return result;
    }

    /**
     * 弱形式的卷曲算子@f[
     * \int_Z \nabla\times u \cdot
     * \nabla \times v \,dx
     * @f]。
     *
     */
    template <int dim>
    void
    curl_curl_matrix(FullMatrix<double> &     M,
                     const FEValuesBase<dim> &fe,
                     const double             factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      // Depending on the dimension,
      // the cross product is either
      // a scalar (2d) or a vector
      // (3d). Accordingly, in the
      // latter case we have to sum
      // up three bilinear forms, but
      // in 2d, we don't. Thus, we
      // need to adapt the loop over
      // all dimensions
      const unsigned int d_max = (dim == 2) ? 1 : dim;

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = factor * fe.JxW(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              for (unsigned int d = 0; d < d_max; ++d)
                {
                  const unsigned int d1 = (d + 1) % dim;
                  const unsigned int d2 = (d + 2) % dim;

                  const double cv = fe.shape_grad_component(i, k, d2)[d1] -
                                    fe.shape_grad_component(i, k, d1)[d2];
                  const double cu = fe.shape_grad_component(j, k, d2)[d1] -
                                    fe.shape_grad_component(j, k, d1)[d2];

                  M(i, j) += dx * cu * cv;
                }
        }
    }

    /**
     * 卷曲算子的矩阵 @f[
     * \int_Z \nabla \times u \cdot v \,dx.
     * @f] 这是3D中的标准卷曲算子和2D中的标量卷曲。矢量curl算子可以通过交换测试和试用函数得到。
     *
     */
    template <int dim>
    void
    curl_matrix(FullMatrix<double> &     M,
                const FEValuesBase<dim> &fe,
                const FEValuesBase<dim> &fetest,
                double                   factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      const unsigned int t_dofs = fetest.dofs_per_cell;
      AssertDimension(fe.get_fe().n_components(), dim);
      // There should be the right number of components (3 in 3D, otherwise 1)
      // for the curl.
      AssertDimension(fetest.get_fe().n_components(), (dim == 3) ? dim : 1);
      AssertDimension(M.m(), t_dofs);
      AssertDimension(M.n(), n_dofs);

      const unsigned int d_max = (dim == 2) ? 1 : dim;

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int i = 0; i < t_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              for (unsigned int d = 0; d < d_max; ++d)
                {
                  const unsigned int d1 = (d + 1) % dim;
                  const unsigned int d2 = (d + 2) % dim;

                  const double vv = fetest.shape_value_component(i, k, d);
                  const double cu = fe.shape_grad_component(j, k, d2)[d1] -
                                    fe.shape_grad_component(j, k, d1)[d2];
                  M(i, j) += dx * cu * vv;
                }
        }
    }

    /**
     * 麦克斯韦系统中切向分量的尼采型弱边界条件的矩阵。        @f[
     * \int_F \biggl( 2\gamma
     * (u\times n) (v\times n)
     *
     * -
     * (u\times n)(\nu \nabla\times
     * v)
     *
     * - (v\times
     * n)(\nu \nabla\times u)
     * \biggr)
     * @f]
     *
     */
    template <int dim>
    void
    nitsche_curl_matrix(FullMatrix<double> &     M,
                        const FEValuesBase<dim> &fe,
                        const unsigned int       face_no,
                        double                   penalty,
                        double                   factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      // Depending on the
      // dimension, the cross
      // product is either a scalar
      // (2d) or a vector
      // (3d). Accordingly, in the
      // latter case we have to sum
      // up three bilinear forms,
      // but in 2d, we don't. Thus,
      // we need to adapt the loop
      // over all dimensions
      const unsigned int d_max = (dim == 2) ? 1 : dim;

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double         dx = factor * fe.JxW(k);
          const Tensor<1, dim> n  = fe.normal_vector(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              if (fe.get_fe().has_support_on_face(i, face_no) &&
                  fe.get_fe().has_support_on_face(j, face_no))
                {
                  for (unsigned int d = 0; d < d_max; ++d)
                    {
                      const unsigned int d1 = (d + 1) % dim;
                      const unsigned int d2 = (d + 2) % dim;

                      const double cv = fe.shape_grad_component(i, k, d2)[d1] -
                                        fe.shape_grad_component(i, k, d1)[d2];
                      const double cu = fe.shape_grad_component(j, k, d2)[d1] -
                                        fe.shape_grad_component(j, k, d1)[d2];
                      const double v =
                        fe.shape_value_component(i, k, d1) * n[d2] -
                        fe.shape_value_component(i, k, d2) * n[d1];
                      const double u =
                        fe.shape_value_component(j, k, d1) * n[d2] -
                        fe.shape_value_component(j, k, d2) * n[d1];

                      M(i, j) += dx * (2. * penalty * u * v - cv * u - cu * v);
                    }
                }
        }
    }
    /**
     * 两个切向迹线的乘积，@f[
     * \int_F (u\times n)(v\times n)
     * \, ds.
     * @f]。
     *
     */
    template <int dim>
    void
    tangential_trace_matrix(FullMatrix<double> &     M,
                            const FEValuesBase<dim> &fe,
                            double                   factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;

      AssertDimension(fe.get_fe().n_components(), dim);
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);

      // Depending on the
      // dimension, the cross
      // product is either a scalar
      // (2d) or a vector
      // (3d). Accordingly, in the
      // latter case we have to sum
      // up three bilinear forms,
      // but in 2d, we don't. Thus,
      // we need to adapt the loop
      // over all dimensions
      const unsigned int d_max = (dim == 2) ? 1 : dim;

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double         dx = factor * fe.JxW(k);
          const Tensor<1, dim> n  = fe.normal_vector(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              for (unsigned int d = 0; d < d_max; ++d)
                {
                  const unsigned int d1 = (d + 1) % dim;
                  const unsigned int d2 = (d + 2) % dim;

                  const double v = fe.shape_value_component(i, k, d1) * n(d2) -
                                   fe.shape_value_component(i, k, d2) * n(d1);
                  const double u = fe.shape_value_component(j, k, d1) * n(d2) -
                                   fe.shape_value_component(j, k, d2) * n(d1);

                  M(i, j) += dx * u * v;
                }
        }
    }

    /**
     * 麦克斯韦系统的内部惩罚通量。        @f[
     * \int_F \biggl( \gamma
     * \{u\times n\}\{v\times n\}
     *
     * -
     * \{u\times n\}\{\nu \nabla\times
     * v\}- \{v\times
     * n\}\{\nu \nabla\times u\}
     * \biggr)\;dx
     * @f]
     *
     */
    template <int dim>
    inline void
    ip_curl_matrix(FullMatrix<double> &     M11,
                   FullMatrix<double> &     M12,
                   FullMatrix<double> &     M21,
                   FullMatrix<double> &     M22,
                   const FEValuesBase<dim> &fe1,
                   const FEValuesBase<dim> &fe2,
                   const double             pen,
                   const double             factor1 = 1.,
                   const double             factor2 = -1.)
    {
      const unsigned int n_dofs = fe1.n_dofs_per_cell();

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

      const double nu1     = factor1;
      const double nu2     = (factor2 < 0) ? factor1 : factor2;
      const double penalty = .5 * pen * (nu1 + nu2);

      // Depending on the
      // dimension, the cross
      // product is either a scalar
      // (2d) or a vector
      // (3d). Accordingly, in the
      // latter case we have to sum
      // up three bilinear forms,
      // but in 2d, we don't. Thus,
      // we need to adapt the loop
      // over all dimensions
      const unsigned int d_max = (dim == 2) ? 1 : dim;

      for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
        {
          const double         dx = fe1.JxW(k);
          const Tensor<1, dim> n  = fe1.normal_vector(k);
          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int j = 0; j < n_dofs; ++j)
              for (unsigned int d = 0; d < d_max; ++d)
                {
                  const unsigned int d1 = (d + 1) % dim;
                  const unsigned int d2 = (d + 2) % dim;
                  // curl u, curl v
                  const double cv1 =
                    nu1 * fe1.shape_grad_component(i, k, d2)[d1] -
                    fe1.shape_grad_component(i, k, d1)[d2];
                  const double cv2 =
                    nu2 * fe2.shape_grad_component(i, k, d2)[d1] -
                    fe2.shape_grad_component(i, k, d1)[d2];
                  const double cu1 =
                    nu1 * fe1.shape_grad_component(j, k, d2)[d1] -
                    fe1.shape_grad_component(j, k, d1)[d2];
                  const double cu2 =
                    nu2 * fe2.shape_grad_component(j, k, d2)[d1] -
                    fe2.shape_grad_component(j, k, d1)[d2];

                  // u x n, v x n
                  const double u1 =
                    fe1.shape_value_component(j, k, d1) * n(d2) -
                    fe1.shape_value_component(j, k, d2) * n(d1);
                  const double u2 =
                    -fe2.shape_value_component(j, k, d1) * n(d2) +
                    fe2.shape_value_component(j, k, d2) * n(d1);
                  const double v1 =
                    fe1.shape_value_component(i, k, d1) * n(d2) -
                    fe1.shape_value_component(i, k, d2) * n(d1);
                  const double v2 =
                    -fe2.shape_value_component(i, k, d1) * n(d2) +
                    fe2.shape_value_component(i, k, d2) * n(d1);

                  M11(i, j) +=
                    .5 * dx * (2. * penalty * u1 * v1 - cv1 * u1 - cu1 * v1);
                  M12(i, j) +=
                    .5 * dx * (2. * penalty * v1 * u2 - cv1 * u2 - cu2 * v1);
                  M21(i, j) +=
                    .5 * dx * (2. * penalty * u1 * v2 - cv2 * u1 - cu1 * v2);
                  M22(i, j) +=
                    .5 * dx * (2. * penalty * u2 * v2 - cv2 * u2 - cu2 * v2);
                }
        }
    }


  } // namespace Maxwell
} // namespace LocalIntegrators


DEAL_II_NAMESPACE_CLOSE

#endif


