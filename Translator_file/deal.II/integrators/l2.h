//include/deal.II-translator/integrators/l2_0.txt
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

#ifndef dealii_integrators_l2_h
#define dealii_integrators_l2_h


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
   * @brief Local integrators related to <i>L<sup>2</sup></i>-inner products.
   * @ingroup Integrators
   *
   */
  namespace L2
  {
    /**
     * 标量或矢量值有限元的质量矩阵。\f[ \int_Z uv\,dx \quad
     * \text{or} \quad \int_Z \mathbf u\cdot \mathbf v\,dx \f]
     * 同样地，这个术语可以用在面，它计算积分 \f[ \int_F
     * uv\,ds \quad \text{or} \quad \int_F \mathbf u\cdot
     * \mathbf v\,ds \f]  @param  M 作为结果得到的质量矩阵。
     * @param  fe
     * 描述局部试验函数空间的FEValues对象。必须设置#update_values和#update_JxW_values。
     * @param  factor 一个常数，与质量矩阵相乘。
     *
     */
    template <int dim>
    void
    mass_matrix(FullMatrix<double> &     M,
                const FEValuesBase<dim> &fe,
                const double             factor = 1.)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * factor;
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              double Mii = 0.0;
              for (unsigned int d = 0; d < n_components; ++d)
                Mii += dx * fe.shape_value_component(i, k, d) *
                       fe.shape_value_component(i, k, d);

              M(i, i) += Mii;

              for (unsigned int j = i + 1; j < n_dofs; ++j)
                {
                  double Mij = 0.0;
                  for (unsigned int d = 0; d < n_components; ++d)
                    Mij += dx * fe.shape_value_component(j, k, d) *
                           fe.shape_value_component(i, k, d);

                  M(i, j) += Mij;
                  M(j, i) += Mij;
                }
            }
        }
    }

    /**
     * 标量或矢量值有限元的加权质量矩阵。    \f[ \int_Z
     * \omega(x) uv\,dx \quad \text{or} \quad \int_Z \omega(x) \mathbf u\cdot
     * \mathbf v\,dx \f]
     * 同样地，这个术语可以用在面，它计算积分 \f[ \int_F
     * \omega(x) uv\,ds \quad \text{or} \quad \int_F
     * \omega(x) \mathbf u\cdot \mathbf v\,ds \f]  @param  M
     * 作为结果得到的加权质量矩阵。      @param  fe
     * 描述局部试验函数空间的FEValues对象。必须设置#update_values和#update_JxW_values。
     * @param  weights 在有限元中正交点评估的权重， $\omega(x)$
     * （大小必须等于元素中正交点的数量）。
     *
     */
    template <int dim>
    void
    weighted_mass_matrix(FullMatrix<double> &       M,
                         const FEValuesBase<dim> &  fe,
                         const std::vector<double> &weights)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int n_components = fe.get_fe().n_components();
      AssertDimension(M.m(), n_dofs);
      AssertDimension(M.n(), n_dofs);
      AssertDimension(weights.size(), fe.n_quadrature_points);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double dx = fe.JxW(k) * weights[k];
          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              double Mii = 0.0;
              for (unsigned int d = 0; d < n_components; ++d)
                Mii += dx * fe.shape_value_component(i, k, d) *
                       fe.shape_value_component(i, k, d);

              M(i, i) += Mii;

              for (unsigned int j = i + 1; j < n_dofs; ++j)
                {
                  double Mij = 0.0;
                  for (unsigned int d = 0; d < n_components; ++d)
                    Mij += dx * fe.shape_value_component(j, k, d) *
                           fe.shape_value_component(i, k, d);

                  M(i, j) += Mij;
                  M(j, i) += Mij;
                }
            }
        }
    }

    /**
     * <i>L<sup>2</sup></i>-标量函数的内积。        \f[ \int_Z fv\,dx
     * \quad \text{or} \quad \int_F fv\,ds \f]  @param  result
     * 作为结果得到的向量。      @param  fe
     * 描述本地试用函数空间的FEValues对象。必须设置#update_values和#update_JxW_values。
     * @param  input 在有限元中正交点评估的 $f$
     * 的表示（大小必须等于元素中正交点的数量）。
     * @param 因子 一个常数，与结果相乘。
     *
     */
    template <int dim, typename number>
    void
    L2(Vector<number> &           result,
       const FEValuesBase<dim> &  fe,
       const std::vector<double> &input,
       const double               factor = 1.)
    {
      const unsigned int n_dofs = fe.dofs_per_cell;
      AssertDimension(result.size(), n_dofs);
      AssertDimension(fe.get_fe().n_components(), 1);
      AssertDimension(input.size(), fe.n_quadrature_points);

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        for (unsigned int i = 0; i < n_dofs; ++i)
          result(i) += fe.JxW(k) * factor * input[k] * fe.shape_value(i, k);
    }

    /**
     * <i>L<sup>2</sup></i>-内积
     * 为矢量值的右手边的一个片断。\f[ \int_Z \mathbf f\cdot
     * \mathbf v\,dx \quad \text{or}
     * \quad \int_F \mathbf f\cdot \mathbf v\,ds \f]  @param  result
     * 作为结果得到的向量。      @param  fe
     * 描述本地试验函数空间的FEValues对象。必须设置#update_values和#update_JxW_values。
     * @param  input 在有限元中的正交点上评估的 $\mathbf f$
     * 的矢量值表示（每个分量的大小必须等于元素中正交点的数量）。
     * @param  因子 一个常数，与结果相乘。
     *
     */
    template <int dim, typename number>
    void
    L2(Vector<number> &                            result,
       const FEValuesBase<dim> &                   fe,
       const ArrayView<const std::vector<double>> &input,
       const double                                factor = 1.)
    {
      const unsigned int n_dofs       = fe.dofs_per_cell;
      const unsigned int n_components = input.size();

      AssertDimension(result.size(), n_dofs);
      AssertDimension(input.size(), fe.get_fe().n_components());

      for (unsigned int k = 0; k < fe.n_quadrature_points; ++k)
        for (unsigned int i = 0; i < n_dofs; ++i)
          for (unsigned int d = 0; d < n_components; ++d)
            result(i) += fe.JxW(k) * factor *
                         fe.shape_value_component(i, k, d) * input[d][k];
    }

    /**
     * 标量或矢量值有限元的两个单元之间的跳跃矩阵。注意，因子
     * $\gamma$ 可以用来实现加权跳变。\f[ \int_F [\gamma u][\gamma
     * v]\,ds \quad \text{or}
     * \int_F [\gamma \mathbf u]\cdot [\gamma \mathbf v]\,ds
     * \f]使用适当的权重，这个词可以用来惩罚<i>H<sup>1</sup></i>中的违反一致性。
     * 请注意，对于后面的参数，外部矩阵指的是单元间的通量，而内部矩阵指的是单元内的条目耦合。
     * @param  M11 作为结果得到的第一个单元的内部矩阵。
     * @param  M12 作为结果得到的第一个单元的外部矩阵。
     * @param  M21 作为结果获得的第二个单元的外部矩阵。
     * @param  M22 作为结果获得的第二个单元的内部矩阵。
     * @param  fe1
     * 描述第一个单元的局部试验函数空间的FEValues对象。必须设置#update_values和#update_JxW_values。
     * @param  fe2
     * 描述第二个单元的本地试用函数空间的FEValues对象。必须设置#update_values和#update_JxW_values。
     * @param  factor1
     * 一个常数，与第一个单元的形状函数相乘。      @param
     * factor2 一个常数，用于乘以第二个单元的形状函数。
     *
     */
    template <int dim>
    void
    jump_matrix(FullMatrix<double> &     M11,
                FullMatrix<double> &     M12,
                FullMatrix<double> &     M21,
                FullMatrix<double> &     M22,
                const FEValuesBase<dim> &fe1,
                const FEValuesBase<dim> &fe2,
                const double             factor1 = 1.,
                const double             factor2 = 1.)
    {
      const unsigned int n1_dofs      = fe1.n_dofs_per_cell();
      const unsigned int n2_dofs      = fe2.n_dofs_per_cell();
      const unsigned int n_components = fe1.get_fe().n_components();

      Assert(n1_dofs == n2_dofs, ExcNotImplemented());
      (void)n2_dofs;
      AssertDimension(n_components, fe2.get_fe().n_components());
      AssertDimension(M11.m(), n1_dofs);
      AssertDimension(M12.m(), n1_dofs);
      AssertDimension(M21.m(), n2_dofs);
      AssertDimension(M22.m(), n2_dofs);
      AssertDimension(M11.n(), n1_dofs);
      AssertDimension(M12.n(), n2_dofs);
      AssertDimension(M21.n(), n1_dofs);
      AssertDimension(M22.n(), n2_dofs);

      for (unsigned int k = 0; k < fe1.n_quadrature_points; ++k)
        {
          const double dx = fe1.JxW(k);

          for (unsigned int i = 0; i < n1_dofs; ++i)
            for (unsigned int j = 0; j < n1_dofs; ++j)
              for (unsigned int d = 0; d < n_components; ++d)
                {
                  const double u1 =
                    factor1 * fe1.shape_value_component(j, k, d);
                  const double u2 =
                    -factor2 * fe2.shape_value_component(j, k, d);
                  const double v1 =
                    factor1 * fe1.shape_value_component(i, k, d);
                  const double v2 =
                    -factor2 * fe2.shape_value_component(i, k, d);

                  M11(i, j) += dx * u1 * v1;
                  M12(i, j) += dx * u2 * v1;
                  M21(i, j) += dx * u1 * v2;
                  M22(i, j) += dx * u2 * v2;
                }
        }
    }
  } // namespace L2
} // namespace LocalIntegrators

DEAL_II_NAMESPACE_CLOSE

#endif


