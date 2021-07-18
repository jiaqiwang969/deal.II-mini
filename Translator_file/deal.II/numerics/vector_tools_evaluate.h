//include/deal.II-translator/numerics/vector_tools_evaluate_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


#ifndef dealii_vector_tools_evaluation_h
#define dealii_vector_tools_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  /**
   * 用于point_values()和point_gradients()的标志的命名空间。
   *
   */
  namespace EvaluationFlags
  {
    /**
     * 用于point_values()和point_gradients()的标志。
     *
     */
    enum EvaluationFlags
    {
      /**
       * 计算平均数。
       *
       */
      avg = 0,
      /**
       * 计算最大值。
       * @note 只对标量值有效。
       *
       */
      max = 1,
      /**
       * 计算最小值。
       * @note  仅适用于标量值。
       *
       */
      min = 2,
      /**
       * 取任何值。
       *
       */
      insert = 3
    };
  } // namespace EvaluationFlags

  /**
   * 给定一个（分布式）解决方案向量  @p vector,  评估由  @p
   * evaluation_points.  指定的（任意甚至是远程）点的值
   * @warning
   * 这是一个集体调用，需要由通信器的所有处理器执行。
   *
   */
  template <int n_components, int dim, int spacedim, typename VectorType>
  std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  point_values(
    const Mapping<dim> &                                  mapping,
    const DoFHandler<dim, spacedim> &                     dof_handler,
    const VectorType &                                    vector,
    const std::vector<Point<spacedim>> &                  evaluation_points,
    Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const EvaluationFlags::EvaluationFlags flags = EvaluationFlags::avg);

  /**
   * 给定一个（分布式）解决方案向量 @p vector, ，评估由 @p
   * cache 指定的点的值，这些点可能是由上述函数设置的。
   * @note
   * 细化/粗化/重新分区导致缓存无效，因此必须再次调用上述函数。
   * @warning
   * 这是一个集体调用，需要由通信器中的所有处理器来执行。
   *
   */
  template <int n_components, int dim, int spacedim, typename VectorType>
  std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  point_values(
    const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const DoFHandler<dim, spacedim> &                           dof_handler,
    const VectorType &                                          vector,
    const EvaluationFlags::EvaluationFlags flags = EvaluationFlags::avg);

  /**
   * 给定一个（分布式）解决方案向量  @p vector,  评估由  @p
   * evaluation_points.
   * 指定的（任意的甚至是遥远的）点的梯度。
   *
   */
  template <int n_components, int dim, int spacedim, typename VectorType>
  std::vector<typename FEPointEvaluation<n_components, dim>::gradient_type>
  point_gradients(
    const Mapping<dim> &                                  mapping,
    const DoFHandler<dim, spacedim> &                     dof_handler,
    const VectorType &                                    vector,
    const std::vector<Point<spacedim>> &                  evaluation_points,
    Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const EvaluationFlags::EvaluationFlags flags = EvaluationFlags::avg);

  /**
   * 给定一个（分布式）解向量 @p vector, ，评估由 @p cache
   * 指定的点的梯度，这可能是由上述函数设置的。
   * @note
   * 细化/粗化/重新划分导致缓存无效，因此必须再次调用上述函数。
   *
   */
  template <int n_components, int dim, int spacedim, typename VectorType>
  std::vector<typename FEPointEvaluation<n_components, dim>::gradient_type>
  point_gradients(
    const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const DoFHandler<dim, spacedim> &                           dof_handler,
    const VectorType &                                          vector,
    const EvaluationFlags::EvaluationFlags flags = EvaluationFlags::avg);



  // inlined functions


#ifndef DOXYGEN
  template <int n_components, int dim, int spacedim, typename VectorType>
  inline std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  point_values(const Mapping<dim> &                mapping,
               const DoFHandler<dim, spacedim> &   dof_handler,
               const VectorType &                  vector,
               const std::vector<Point<spacedim>> &evaluation_points,
               Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
               const EvaluationFlags::EvaluationFlags                flags)
  {
    cache.reinit(evaluation_points, dof_handler.get_triangulation(), mapping);

    return point_values<n_components>(cache, dof_handler, vector, flags);
  }



  template <int n_components, int dim, int spacedim, typename VectorType>
  inline std::vector<
    typename FEPointEvaluation<n_components, dim>::gradient_type>
  point_gradients(const Mapping<dim> &                mapping,
                  const DoFHandler<dim, spacedim> &   dof_handler,
                  const VectorType &                  vector,
                  const std::vector<Point<spacedim>> &evaluation_points,
                  Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
                  const EvaluationFlags::EvaluationFlags                flags)
  {
    cache.reinit(evaluation_points, dof_handler.get_triangulation(), mapping);

    return point_gradients<n_components>(cache, dof_handler, vector, flags);
  }



  namespace internal
  {
    /**
     * 对标量进行还原。
     *
     */
    template <typename T>
    T
    reduce(const EvaluationFlags::EvaluationFlags &flags,
           const ArrayView<const T> &              values)
    {
      switch (flags)
        {
          case EvaluationFlags::avg:
            {
              return std::accumulate(values.begin(), values.end(), T{}) /
                     (T(1.0) * values.size());
            }
          case EvaluationFlags::max:
            return *std::max_element(values.begin(), values.end());
          case EvaluationFlags::min:
            return *std::min_element(values.begin(), values.end());
          case EvaluationFlags::insert:
            return values[0];
          default:
            Assert(false, ExcNotImplemented());
            return values[0];
        }
    }



    /**
     * 对张量进行还原。
     *
     */
    template <int rank, int dim, typename Number>
    Tensor<rank, dim, Number>
    reduce(const EvaluationFlags::EvaluationFlags &          flags,
           const ArrayView<const Tensor<rank, dim, Number>> &values)
    {
      switch (flags)
        {
          case EvaluationFlags::avg:
            {
              return std::accumulate(values.begin(),
                                     values.end(),
                                     Tensor<rank, dim, Number>{}) /
                     (Number(1.0) * values.size());
            }
          case EvaluationFlags::insert:
            return values[0];
          default:
            Assert(false, ExcNotImplemented());
            return values[0];
        }
    }



    template <int n_components,
              int dim,
              int spacedim,
              typename VectorType,
              typename value_type>
    inline std::vector<value_type>
    evaluate_at_points(
      const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
      const DoFHandler<dim, spacedim> &                           dof_handler,
      const VectorType &                                          vector,
      const EvaluationFlags::EvaluationFlags                      flags,
      const UpdateFlags                                           update_flags,
      const dealii::EvaluationFlags::EvaluationFlags evaluation_flags,
      const std::function<
        value_type(const FEPointEvaluation<n_components, dim> &,
                   const unsigned int &)> process_quadrature_point)
    {
      Assert(cache.is_ready(),
             ExcMessage(
               "Utilities::MPI::RemotePointEvaluation is not ready yet! "
               "Please call Utilities::MPI::RemotePointEvaluation::reinit() "
               "yourself or another function that does this for you."));

      Assert(
        &dof_handler.get_triangulation() == &cache.get_triangulation(),
        ExcMessage(
          "The provided Utilities::MPI::RemotePointEvaluation and DoFHandler "
          "object have been set up with different Triangulation objects, "
          "a scenario not supported!"));

      // evaluate values at points if possible
      const auto evaluation_point_results = [&]() {
        // helper function for accessing the global vector and interpolating
        // the results onto the points
        const auto evaluation_function = [&](auto &      values,
                                             const auto &cell_data) {
          std::vector<typename VectorType::value_type> solution_values;

          std::vector<std::unique_ptr<FEPointEvaluation<n_components, dim>>>
            evaluators(dof_handler.get_fe_collection().size());

          const auto get_evaluator = [&](const unsigned int active_fe_index)
            -> FEPointEvaluation<n_components, dim> & {
            if (evaluators[active_fe_index] == nullptr)
              evaluators[active_fe_index] =
                std::make_unique<FEPointEvaluation<n_components, dim>>(
                  cache.get_mapping(),
                  dof_handler.get_fe(active_fe_index),
                  update_flags);

            return *evaluators[active_fe_index];
          };

          for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
            {
              typename DoFHandler<dim>::active_cell_iterator cell = {
                &cache.get_triangulation(),
                cell_data.cells[i].first,
                cell_data.cells[i].second,
                &dof_handler};

              const ArrayView<const Point<dim>> unit_points(
                cell_data.reference_point_values.data() +
                  cell_data.reference_point_ptrs[i],
                cell_data.reference_point_ptrs[i + 1] -
                  cell_data.reference_point_ptrs[i]);

              solution_values.resize(
                dof_handler.get_fe(cell->active_fe_index()).n_dofs_per_cell());
              cell->get_dof_values(vector,
                                   solution_values.begin(),
                                   solution_values.end());

              auto &evaluator = get_evaluator(cell->active_fe_index());

              evaluator.reinit(cell, unit_points);
              evaluator.evaluate(solution_values, evaluation_flags);

              for (unsigned int q = 0; q < unit_points.size(); ++q)
                values[q + cell_data.reference_point_ptrs[i]] =
                  process_quadrature_point(evaluator, q);
            }
        };

        std::vector<value_type> evaluation_point_results;
        std::vector<value_type> buffer;

        cache.template evaluate_and_process<value_type>(
          evaluation_point_results, buffer, evaluation_function);

        return evaluation_point_results;
      }();

      if (cache.is_map_unique())
        {
          // each point has exactly one result (unique map)
          return evaluation_point_results;
        }
      else
        {
          // map is not unique (multiple or no results): postprocessing is
          // needed
          std::vector<value_type> unique_evaluation_point_results(
            cache.get_point_ptrs().size() - 1);

          const auto &ptr = cache.get_point_ptrs();

          for (unsigned int i = 0; i < ptr.size() - 1; ++i)
            {
              const auto n_entries = ptr[i + 1] - ptr[i];
              if (n_entries == 0)
                continue;

              unique_evaluation_point_results[i] =
                reduce(flags,
                       ArrayView<const value_type>(
                         evaluation_point_results.data() + ptr[i], n_entries));
            }

          return unique_evaluation_point_results;
        }
    }
  } // namespace internal

  template <int n_components, int dim, int spacedim, typename VectorType>
  inline std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  point_values(
    const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const DoFHandler<dim, spacedim> &                           dof_handler,
    const VectorType &                                          vector,
    const EvaluationFlags::EvaluationFlags                      flags)
  {
    return internal::evaluate_at_points<
      n_components,
      dim,
      spacedim,
      VectorType,
      typename FEPointEvaluation<n_components, dim>::value_type>(
      cache,
      dof_handler,
      vector,
      flags,
      update_values,
      dealii::EvaluationFlags::values,
      [](const auto &evaluator, const auto &q) {
        return evaluator.get_value(q);
      });
  }

  template <int n_components, int dim, int spacedim, typename VectorType>
  inline std::vector<
    typename FEPointEvaluation<n_components, dim>::gradient_type>
  point_gradients(
    const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const DoFHandler<dim, spacedim> &                           dof_handler,
    const VectorType &                                          vector,
    const EvaluationFlags::EvaluationFlags                      flags)
  {
    return internal::evaluate_at_points<
      n_components,
      dim,
      spacedim,
      VectorType,
      typename FEPointEvaluation<n_components, dim>::gradient_type>(
      cache,
      dof_handler,
      vector,
      flags,
      update_gradients,
      dealii::EvaluationFlags::gradients,
      [](const auto &evaluator, const auto &q) {
        return evaluator.get_gradient(q);
      });
  }

#endif
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_boundary_h


