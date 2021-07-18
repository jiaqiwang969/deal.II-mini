//include/deal.II-translator/matrix_free/evaluation_selector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_selector_h
#define dealii_matrix_free_evaluation_selector_h

#include <deal.II/base/config.h>

#include <deal.II/matrix_free/evaluation_kernels.h>

DEAL_II_NAMESPACE_OPEN



/**
 * 这个类根据模板参数和shape_info变量选择合适的评估策略，为一些内部函数提供捷径。
 *
 *
 */
template <int dim, int fe_degree, int n_q_points_1d, typename Number>
struct SelectEvaluator
{
  /**
   * 为evaluate函数选择合适的评估策略，即这调用
   * internal::FEEvaluationImpl::evaluate(),
   * internal::FEEvaluationImplCollocation::evaluate() 或
   * internal::FEEvaluationImplTransformToCollocation::evaluate()
   * 的适当模板参数。
   *
   */
  static void
  evaluate(const unsigned int                     n_components,
           const EvaluationFlags::EvaluationFlags evaluation_flag,
           const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
           Number *values_dofs_actual,
           Number *values_quad,
           Number *gradients_quad,
           Number *hessians_quad,
           Number *scratch_data);

  /**
   * 为积分函数选择合适的评估策略，即调用
   * internal::FEEvaluationImpl::integrate(),
   * internal::FEEvaluationImplCollocation::integrate() 或
   * internal::FEEvaluationImplTransformToCollocation::integrate()
   * ，并带有合适的模板参数。
   *
   */
  static void
  integrate(const unsigned int                     n_components,
            const EvaluationFlags::EvaluationFlags integration_flag,
            const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
            Number *   values_dofs_actual,
            Number *   values_quad,
            Number *   gradients_quad,
            Number *   scratch_data,
            const bool sum_into_values_array = false);
};

//----------------------Implementation for SelectEvaluator---------------------
#ifndef DOXYGEN

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
inline void
SelectEvaluator<dim, fe_degree, n_q_points_1d, Number>::evaluate(
  const unsigned int                                      n_components,
  const EvaluationFlags::EvaluationFlags                  evaluation_flag,
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
  Number *                                                values_dofs_actual,
  Number *                                                values_quad,
  Number *                                                gradients_quad,
  Number *                                                hessians_quad,
  Number *                                                scratch_data)
{
  Assert(fe_degree >= 0 && n_q_points_1d > 0, ExcInternalError());

  internal::FEEvaluationImplEvaluateSelector<dim, Number>::
    template run<fe_degree, n_q_points_1d>(n_components,
                                           evaluation_flag,
                                           shape_info,
                                           values_dofs_actual,
                                           values_quad,
                                           gradients_quad,
                                           hessians_quad,
                                           scratch_data);
}



template <int dim, int fe_degree, int n_q_points_1d, typename Number>
inline void
SelectEvaluator<dim, fe_degree, n_q_points_1d, Number>::integrate(
  const unsigned int                                      n_components,
  const EvaluationFlags::EvaluationFlags                  integration_flag,
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
  Number *                                                values_dofs_actual,
  Number *                                                values_quad,
  Number *                                                gradients_quad,
  Number *                                                scratch_data,
  const bool                                              sum_into_values_array)
{
  Assert(fe_degree >= 0 && n_q_points_1d > 0, ExcInternalError());

  internal::FEEvaluationImplIntegrateSelector<dim, Number>::
    template run<fe_degree, n_q_points_1d>(n_components,
                                           integration_flag,
                                           shape_info,
                                           values_dofs_actual,
                                           values_quad,
                                           gradients_quad,
                                           scratch_data,
                                           sum_into_values_array);
}
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


