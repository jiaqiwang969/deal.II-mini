//include/deal.II-translator/matrix_free/fe_point_evaluation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_fe_point_evaluation_h
#define dealii_fe_point_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FEPointEvaluation
  {
    /**
     * 用于区分FlexibleEvaluator类所使用的不同数量组件的值和梯度类型的结构。
     *
     */
    template <int dim, int n_components, typename Number>
    struct EvaluatorTypeTraits
    {
      using value_type    = Tensor<1, n_components, Number>;
      using gradient_type = Tensor<1, n_components, Tensor<1, dim, Number>>;

      static void
      read_value(const Number       vector_entry,
                 const unsigned int component,
                 value_type &       result)
      {
        AssertIndexRange(component, n_components);
        result[component] = vector_entry;
      }

      static void
      write_value(Number &           vector_entry,
                  const unsigned int component,
                  const value_type & result)
      {
        AssertIndexRange(component, n_components);
        vector_entry = result[component];
      }

      static void
      set_gradient(
        const Tensor<1, dim, Tensor<1, n_components, VectorizedArray<Number>>>
          &                value,
        const unsigned int vector_lane,
        gradient_type &    result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i][vector_lane];
      }

      static void get_gradient(
        Tensor<1, dim, Tensor<1, n_components, VectorizedArray<Number>>> &value,
        const unsigned int   vector_lane,
        const gradient_type &result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i][vector_lane] = result[i][d];
      }

      static void
      set_value(const Tensor<1, n_components, VectorizedArray<Number>> &value,
                const unsigned int vector_lane,
                value_type &       result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          result[i] = value[i][vector_lane];
      }

      static void
        get_value(Tensor<1, n_components, VectorizedArray<Number>> &value,
                  const unsigned int                                vector_lane,
                  const value_type &                                result)
      {
        for (unsigned int i = 0; i < n_components; ++i)
          value[i][vector_lane] = result[i];
      }

      template <typename Number2>
      static Number2 &access(Tensor<1, n_components, Number2> &value,
                             const unsigned int                component)
      {
        return value[component];
      }

      template <typename Number2>
      static const Number2 &
      access(const Tensor<1, n_components, Number2> &value,
             const unsigned int                      component)
      {
        return value[component];
      }
    };

    template <int dim, typename Number>
    struct EvaluatorTypeTraits<dim, 1, Number>
    {
      using value_type    = Number;
      using gradient_type = Tensor<1, dim, Number>;

      static void
      read_value(const Number vector_entry,
                 const unsigned int,
                 value_type &result)
      {
        result = vector_entry;
      }

      static void
      write_value(Number &vector_entry,
                  const unsigned int,
                  const value_type &result)
      {
        vector_entry = result;
      }

      static void
      set_gradient(const Tensor<1, dim, VectorizedArray<Number>> &value,
                   const unsigned int                             vector_lane,
                   gradient_type &                                result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          result[d] = value[d][vector_lane];
      }

      static void get_gradient(Tensor<1, dim, VectorizedArray<Number>> &value,
                               const unsigned int   vector_lane,
                               const gradient_type &result)
      {
        for (unsigned int d = 0; d < dim; ++d)
          value[d][vector_lane] = result[d];
      }

      static void
      set_value(const VectorizedArray<Number> &value,
                const unsigned int             vector_lane,
                value_type &                   result)
      {
        result = value[vector_lane];
      }

      static void
      get_value(VectorizedArray<Number> &value,
                const unsigned int       vector_lane,
                const value_type &       result)
      {
        value[vector_lane] = result;
      }

      template <typename Number2>
      static Number2 &
      access(Number2 &value, const unsigned int)
      {
        return value;
      }

      template <typename Number2>
      static const Number2 &
      access(const Number2 &value, const unsigned int)
      {
        return value;
      }
    };

    template <int dim, typename Number>
    struct EvaluatorTypeTraits<dim, dim, Number>
    {
      using value_type    = Tensor<1, dim, Number>;
      using gradient_type = Tensor<2, dim, Number>;

      static void
      read_value(const Number       vector_entry,
                 const unsigned int component,
                 value_type &       result)
      {
        result[component] = vector_entry;
      }

      static void
      write_value(Number &           vector_entry,
                  const unsigned int component,
                  const value_type & result)
      {
        vector_entry = result[component];
      }

      static void
      set_gradient(
        const Tensor<1, dim, Tensor<1, dim, VectorizedArray<Number>>> &value,
        const unsigned int vector_lane,
        gradient_type &    result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            result[i][d] = value[d][i][vector_lane];
      }

      static void get_gradient(
        Tensor<1, dim, Tensor<1, dim, VectorizedArray<Number>>> &value,
        const unsigned int                                       vector_lane,
        const gradient_type &                                    result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int d = 0; d < dim; ++d)
            value[d][i][vector_lane] = result[i][d];
      }

      static void
      set_value(const Tensor<1, dim, VectorizedArray<Number>> &value,
                const unsigned int                             vector_lane,
                value_type &                                   result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          result[i] = value[i][vector_lane];
      }

      static void get_value(Tensor<1, dim, VectorizedArray<Number>> &value,
                            const unsigned int vector_lane,
                            const value_type & result)
      {
        for (unsigned int i = 0; i < dim; ++i)
          value[i][vector_lane] = result[i];
      }

      static Number &
      access(value_type &value, const unsigned int component)
      {
        return value[component];
      }

      static const Number &
      access(const value_type &value, const unsigned int component)
      {
        return value[component];
      }

      static Tensor<1, dim> &
      access(gradient_type &value, const unsigned int component)
      {
        return value[component];
      }

      static const Tensor<1, dim> &
      access(const gradient_type &value, const unsigned int component)
      {
        return value[component];
      }
    };

    template <typename Number>
    struct EvaluatorTypeTraits<1, 1, Number>
    {
      using value_type    = Number;
      using gradient_type = Tensor<1, 1, Number>;

      static void
      read_value(const Number vector_entry,
                 const unsigned int,
                 value_type &result)
      {
        result = vector_entry;
      }

      static void
      write_value(Number &vector_entry,
                  const unsigned int,
                  const value_type &result)
      {
        vector_entry = result;
      }

      static void
      set_gradient(const Tensor<1, 1, VectorizedArray<Number>> &value,
                   const unsigned int                           vector_lane,
                   gradient_type &                              result)
      {
        result[0] = value[0][vector_lane];
      }

      static void get_gradient(Tensor<1, 1, VectorizedArray<Number>> &value,
                               const unsigned int   vector_lane,
                               const gradient_type &result)
      {
        value[0][vector_lane] = result[0];
      }

      static void
      set_value(const VectorizedArray<Number> &value,
                const unsigned int             vector_lane,
                value_type &                   result)
      {
        result = value[vector_lane];
      }

      static void
      get_value(VectorizedArray<Number> &value,
                const unsigned int       vector_lane,
                const value_type &       result)
      {
        value[vector_lane] = result;
      }

      template <typename Number2>
      static Number2 &
      access(Number2 &value, const unsigned int)
      {
        return value;
      }

      template <typename Number2>
      static const Number2 &
      access(const Number2 &value, const unsigned int)
      {
        return value;
      }
    };

    template <int dim, int spacedim>
    bool
    is_fast_path_supported(const FiniteElement<dim, spacedim> &fe,
                           const unsigned int base_element_number);

    template <int dim, int spacedim>
    std::vector<Polynomials::Polynomial<double>>
    get_polynomial_space(const FiniteElement<dim, spacedim> &fe);
  } // namespace FEPointEvaluation
} // namespace internal



/**
 * 该类提供了一个接口，用于评估任意参考点位置上的单元格的内插解数值和梯度。这些点可以从一个单元格到另一个单元格，在数量和位置方面都可以改变。两个典型的用例是在非匹配网格上的评估和粒子模拟。
 * 该类的使用与FEValues或FEEvaluation类似。该类首先通过调用
 * `FEPointEvaluation::reinit(cell,
 * unit_points)`初始化为一个单元，与其他概念的主要区别是需要传递参考坐标的底层点。然后，在调用evaluation()或integration()时，用户可以计算给定点的信息。最终，访问函数get_value()或get_gradient()允许在一个特定的点索引上查询这些信息。
 * 该功能类似于在每个单元的`unit_points`上分别创建一个带有正交对象的FEValues对象，然后调用 FEValues::get_function_values
 * 或 FEValues::get_function_gradients,
 * ，对于一些元素和映射，这就是内部实际发生的情况。然而，对于映射和有限元素的特定组合，有一种更有效的实现方式，可以避免FEValues的内存分配和其他昂贵的启动成本。目前，该功能专门用于从MappingQGeneric派生的映射，以及与
 * @ref matrixfree
 * 模块一起工作的具有张量积结构的有限元。在这些情况下，该类所隐含的成本与使用
 * `FEValues::reinit(cell)` 和 `FEValues::get_function_gradients`.
 * 相似（有时甚至更低）。
 *
 *
 */
template <int n_components,
          int dim,
          int spacedim    = dim,
          typename Number = double>
class FEPointEvaluation
{
public:
  using value_type = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, n_components, Number>::value_type;
  using gradient_type = typename internal::FEPointEvaluation::
    EvaluatorTypeTraits<dim, n_components, Number>::gradient_type;

  /**
   * 构建器。      @param  mapping
   * 描述传递给evaluation()函数的单元格的实际几何形状的Mapping类。
   * @param  fe
   * 用于评估的FiniteElement对象，通常在所有要评估的单元上都是一样的。
   * @param  update_flags
   * 指定在调用reinit()时由映射计算的数量。在evaluation()或integration()期间，这些数据被查询以产生所需的结果（例如，有限元解的梯度）。
   * @param  first_selected_component
   * 对于多分量的FiniteElement对象，该参数允许从该参数开始选择`n_components`分量范围。
   *
   */
  FEPointEvaluation(const Mapping<dim> &      mapping,
                    const FiniteElement<dim> &fe,
                    const UpdateFlags         update_flags,
                    const unsigned int        first_selected_component = 0);

  /**
   * 设置给定单元的映射信息，例如，如果要求函数的梯度，通过计算给定点的映射的雅各布系数。
   * @param[in]  cell 当前单元的迭代器  @param[in]  unit_points
   * 当前单元参考位置的点的列表，在evaluate()和integration()函数中，FiniteElement对象应被评估/整合。
   *
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const ArrayView<const Point<dim>> &unit_points);

  /**
   * 这个函数在传递给reinit()的单元和`unit_points`上插值有限元解，由`solution_values`表示。
   * @param[in]  solution_values
   * 这个数组应该包含由`cell->get_dof_values(global_vector,
   * solution_values)'返回的元素上的未知值。      @param[in]
   * evaluation_flags 指明哪些量应该在点上被评估的标志。
   *
   */
  void
  evaluate(const ArrayView<const Number> &         solution_values,
           const EvaluationFlags::EvaluationFlags &evaluation_flags);

  /**
   * 这个函数将之前submit_value()或submit_gradient()调用所传递的量乘以测试函数的值或梯度，并对所有给定的点进行求和。这类似于测试函数的双线性形式的积分，不同的是这个公式不包括`JxW`因子。这就使得该类方法可以自然地将点信息（如粒子）嵌入到有限元公式中。当然，通过给submit_value()的数据乘以一个`JxW`信息，积分也可以用这个类来表示。
   * @param[out]  solution_values
   * 这个数组将包含积分的结果，可以用来在`cell->set_dof_values(solution_values,
   * global_vector)`或`cell->distribute_local_to_global(solution_values,
   * global_vector)`时使用。注意，对于多分量系统，只有部分分量被本类选择，本类未触及的条目将被清零。
   * @param[in]  integration_flags
   * 指明哪些量应该在点上进行积分的标志。
   *
   */
  void
  integrate(const ArrayView<Number> &               solution_values,
            const EvaluationFlags::EvaluationFlags &integration_flags);

  /**
   * 在调用 FEPointEvaluation::evaluate() 并设置了
   * EvaluationFlags::value 后，返回正交点编号 @p point_index
   * 的值，或者调用 FEPointEvaluation::submit_value().
   * 后，返回已经存储在那里的值，如果该对象是矢量值的，将给出一个矢量值的返回参数。
   *
   */
  const value_type &
  get_value(const unsigned int point_index) const;

  /**
   * 写一个值到包含有point_index组件的点上的值的字段。与通过get_value()访问同一字段。如果在调用设置了
   * EvaluationFlags::values 的函数 FEPointEvaluation::integrate()
   * 之前应用，这指定了由当前单元格上的所有基函数测试并整合的值。
   *
   */
  void
  submit_value(const value_type &value, const unsigned int point_index);

  /**
   * 在调用 FEPointEvaluation::evaluate() 并设置
   * EvaluationFlags::gradient
   * 后，返回索引为`point_index`的点的实坐标梯度，或者调用
   * FEPointEvaluation::submit_gradient(). 后存储在那里的梯度。
   * 实坐标梯度是通过获取单位梯度（也可以通过get_unit_gradient()获取）并应用映射的逆雅各布系数获得。如果对象是矢量值的，则会给出一个矢量值的返回参数。
   *
   */
  const gradient_type &
  get_gradient(const unsigned int point_index) const;

  /**
   * 在调用 FEPointEvaluation::evaluate() 并设置
   * EvaluationFlags::gradient
   * 后，返回索引为`point_index`的点的单位坐标梯度，或者在调用
   * FEPointEvaluation::submit_gradient().
   * 后，返回已经存储在那里的梯度
   * 如果对象是矢量值的，将给出一个矢量值的返回参数。请注意，当矢量化被启用时，来自几个点的值会被分组在一起。
   *
   */
  const gradient_type &
  get_unit_gradient(const unsigned int point_index) const;

  /**
   * 写一个贡献值，这个贡献值被梯度测试到包含给定`point_index`的点上的值的区域。与通过get_gradient()访问的字段相同。如果在调用函数
   * FEPointEvaluation::integrate(EvaluationFlags::gradients)
   * 之前应用，这将指定由当前单元上的所有基函数梯度测试的内容，并对其进行整合。
   *
   */
  void
  submit_gradient(const gradient_type &, const unsigned int point_index);

  /**
   * 返回当前单元格上给定点指数的变换的雅各布系数。前提是。这个类需要用包含`update_jacobian`的UpdateFlags构建。
   *
   */
  DerivativeForm<1, dim, spacedim>
  jacobian(const unsigned int point_index) const;

  /**
   * 返回当前单元格上给定点指数的变形的雅各布系数的逆值。前提是。这个类需要用包含`update_inverse_jacobian`或`update_gradients`的UpdateFlags构建。
   *
   */
  DerivativeForm<1, spacedim, dim>
  inverse_jacobian(const unsigned int point_index) const;

  /**
   * 返回传给reinit()的点中，给定的点的实坐标位置。
   *
   */
  Point<spacedim>
  real_point(const unsigned int point_index) const;

  /**
   * 返回给定的点索引的单位/参考坐标位置，即传递给reinit()函数的各个点。
   *
   */
  Point<dim>
  unit_point(const unsigned int point_index) const;

private:
  /**
   * 指向传递给构造函数的Mapping对象的指针。
   *
   */
  SmartPointer<const Mapping<dim, spacedim>> mapping;

  /**
   * 指向MappingQGeneric类的指针，使该类的快速路径。
   *
   */
  const MappingQGeneric<dim, spacedim> *mapping_q_generic;

  /**
   * 指向传递给构造函数的FiniteElement对象的指针。
   *
   */
  SmartPointer<const FiniteElement<dim>> fe;

  /**
   * 描述用于该类的快速路径的张量积元素的一维多项式基础，使用张量积评估器。
   *
   */
  std::vector<Polynomials::Polynomial<double>> poly;

  /**
   * 存储多项式是否是线性的，节点在0和1。
   *
   */
  bool polynomials_are_hat_functions;

  /**
   * 在FiniteElement类所隐含的未知数的未知数和用于张量代码路径的词典式编号之间重新编号。
   *
   */
  std::vector<unsigned int> renumber;

  /**
   * 临时数组，用于存储传递给evaluate()函数的`solution_values'，其格式与张量积评估器兼容。对于矢量值的设置，这个数组使用`张量<1,
   * n_components>`类型来收集特定基函数的未知数。
   *
   */
  std::vector<value_type> solution_renumbered;

  /**
   * 临时数组用于存储在`integrate()`过程中计算的`solution_values`的矢量版本，其格式与张量积评估器兼容。对于矢量值的设置，这个数组使用`Tensor<1,
   * n_components, VectorizedArray<Number>>格式。
   *
   */
  AlignedVector<typename internal::FEPointEvaluation::EvaluatorTypeTraits<
    dim,
    n_components,
    VectorizedArray<Number>>::value_type>
    solution_renumbered_vectorized;

  /**
   * 临时数组来存储各点的数值。
   *
   */
  std::vector<value_type> values;

  /**
   * 临时数组用于存储各点的单位坐标的梯度。
   *
   */
  std::vector<gradient_type> unit_gradients;

  /**
   * 临时数组，用于存储各点的实坐标梯度。
   *
   */
  std::vector<gradient_type> gradients;

  /**
   * 每个组件的未知数，即唯一的基函数的数量，对于所选择的FiniteElement（或基础元素）。
   *
   */
  unsigned int dofs_per_component;

  /**
   * 对于复杂的FiniteElement对象，这个变量告诉我们哪些未知数在所选组件中实际带有自由度。
   *
   */
  std::vector<std::array<bool, n_components>> nonzero_shape_function_component;

  /**
   * 评估所需的更新标志。
   *
   */
  UpdateFlags update_flags;

  /**
   * 快速评估路径中特定于映射的更新标志。
   *
   */
  UpdateFlags update_flags_mapping;

  /**
   * 慢速评估路径下的FEValues对象。
   *
   */
  std::shared_ptr<FEValues<dim, spacedim>> fe_values;

  /**
   * 用于存储快速评估路径的映射所计算的临时数据的数组。
   *
   */
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    mapping_data;

  /**
   * 在reinit()中指定的参考点。
   *
   */
  std::vector<Point<dim>> unit_points;
};

// ----------------------- template and inline function ----------------------


template <int n_components, int dim, int spacedim, typename Number>
FEPointEvaluation<n_components, dim, spacedim, Number>::FEPointEvaluation(
  const Mapping<dim> &      mapping,
  const FiniteElement<dim> &fe,
  const UpdateFlags         update_flags,
  const unsigned int        first_selected_component)
  : mapping(&mapping)
  , mapping_q_generic(
      dynamic_cast<const MappingQGeneric<dim, spacedim> *>(&mapping))
  , fe(&fe)
  , update_flags(update_flags)
  , update_flags_mapping(update_default)
{
  bool         same_base_element   = true;
  unsigned int base_element_number = 0;
  unsigned int component           = 0;
  for (; base_element_number < fe.n_base_elements(); ++base_element_number)
    if (component + fe.element_multiplicity(base_element_number) >
        first_selected_component)
      {
        if (first_selected_component + n_components >
            component + fe.element_multiplicity(base_element_number))
          same_base_element = false;
        break;
      }
    else
      component += fe.element_multiplicity(base_element_number);
  if (mapping_q_generic != nullptr &&
      internal::FEPointEvaluation::is_fast_path_supported(
        fe, base_element_number) &&
      same_base_element)
    {
      internal::MatrixFreeFunctions::ShapeInfo<double> shape_info;

      shape_info.reinit(QMidpoint<1>(), fe, base_element_number);
      renumber           = shape_info.lexicographic_numbering;
      dofs_per_component = shape_info.dofs_per_component_on_cell;
      poly               = internal::FEPointEvaluation::get_polynomial_space(
        fe.base_element(base_element_number));

      polynomials_are_hat_functions =
        (poly.size() == 2 && poly[0].value(0.) == 1. &&
         poly[0].value(1.) == 0. && poly[1].value(0.) == 0. &&
         poly[1].value(1.) == 1.);
    }
  else
    {
      nonzero_shape_function_component.resize(fe.n_dofs_per_cell());
      for (unsigned int d = 0; d < n_components; ++d)
        {
          const unsigned int component = first_selected_component + d;
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
            {
              const bool is_primitive = fe.is_primitive() || fe.is_primitive(i);
              if (is_primitive)
                nonzero_shape_function_component[i][d] =
                  (component == fe.system_to_component_index(i).first);
              else
                nonzero_shape_function_component[i][d] =
                  (fe.get_nonzero_components(i)[component] == true);
            }
        }
    }

  // translate update flags
  if (update_flags & update_jacobians)
    update_flags_mapping |= update_jacobians;
  if (update_flags & update_gradients ||
      update_flags & update_inverse_jacobians)
    update_flags_mapping |= update_inverse_jacobians;
  if (update_flags & update_quadrature_points)
    update_flags_mapping |= update_quadrature_points;
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const ArrayView<const Point<dim>> &                         unit_points)
{
  this->unit_points.resize(unit_points.size());
  std::copy(unit_points.begin(), unit_points.end(), this->unit_points.begin());

  if (!poly.empty())
    mapping_q_generic->fill_mapping_data_for_generic_points(
      cell, unit_points, update_flags_mapping, mapping_data);
  else
    {
      fe_values = std::make_shared<FEValues<dim, spacedim>>(
        *mapping,
        *fe,
        Quadrature<dim>(
          std::vector<Point<dim>>(unit_points.begin(), unit_points.end())),
        update_flags | update_flags_mapping);
      fe_values->reinit(cell);
      mapping_data.initialize(unit_points.size(), update_flags_mapping);
      if (update_flags_mapping & update_jacobians)
        for (unsigned int q = 0; q < unit_points.size(); ++q)
          mapping_data.jacobians[q] = fe_values->jacobian(q);
      if (update_flags_mapping & update_inverse_jacobians)
        for (unsigned int q = 0; q < unit_points.size(); ++q)
          mapping_data.inverse_jacobians[q] = fe_values->inverse_jacobian(q);
      if (update_flags_mapping & update_quadrature_points)
        for (unsigned int q = 0; q < unit_points.size(); ++q)
          mapping_data.quadrature_points[q] = fe_values->quadrature_point(q);
    }

  if (update_flags & update_values)
    values.resize(unit_points.size(), numbers::signaling_nan<value_type>());
  if (update_flags & update_gradients)
    gradients.resize(unit_points.size(),
                     numbers::signaling_nan<gradient_type>());
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::evaluate(
  const ArrayView<const Number> &         solution_values,
  const EvaluationFlags::EvaluationFlags &evaluation_flag)
{
  if (unit_points.empty())
    return;

  AssertDimension(solution_values.size(), fe->dofs_per_cell);
  if (((evaluation_flag & EvaluationFlags::values) ||
       (evaluation_flag & EvaluationFlags::gradients)) &&
      !poly.empty())
    {
      // fast path with tensor product evaluation
      if (solution_renumbered.size() != dofs_per_component)
        solution_renumbered.resize(dofs_per_component);
      for (unsigned int comp = 0; comp < n_components; ++comp)
        for (unsigned int i = 0; i < dofs_per_component; ++i)
          internal::FEPointEvaluation::
            EvaluatorTypeTraits<dim, n_components, Number>::read_value(
              solution_values[renumber[comp * dofs_per_component + i]],
              comp,
              solution_renumbered[i]);

      // unit gradients are currently only implemented with the fast tensor
      // path
      unit_gradients.resize(unit_points.size(),
                            numbers::signaling_nan<gradient_type>());

      const std::size_t n_points = unit_points.size();
      const std::size_t n_lanes  = VectorizedArray<Number>::size();
      for (unsigned int i = 0; i < n_points; i += n_lanes)
        {
          // convert to vectorized format
          Point<dim, VectorizedArray<Number>> vectorized_points;
          for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              vectorized_points[d][j] = unit_points[i + j][d];

          // compute
          const auto val_and_grad =
            internal::evaluate_tensor_product_value_and_gradient(
              poly,
              solution_renumbered,
              vectorized_points,
              polynomials_are_hat_functions);

          // convert back to standard format
          if (evaluation_flag & EvaluationFlags::values)
            for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
              internal::FEPointEvaluation::
                EvaluatorTypeTraits<dim, n_components, Number>::set_value(
                  val_and_grad.first, j, values[i + j]);
          if (evaluation_flag & EvaluationFlags::gradients)
            {
              Assert(update_flags & update_gradients ||
                       update_flags & update_inverse_jacobians,
                     ExcNotInitialized());
              for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
                {
                  Assert(update_flags_mapping & update_inverse_jacobians,
                         ExcNotInitialized());
                  internal::FEPointEvaluation::EvaluatorTypeTraits<
                    dim,
                    n_components,
                    Number>::set_gradient(val_and_grad.second,
                                          j,
                                          unit_gradients[i + j]);
                  gradients[i + j] = apply_transformation(
                    mapping_data.inverse_jacobians[i + j].transpose(),
                    unit_gradients[i + j]);
                }
            }
        }
    }
  else if ((evaluation_flag & EvaluationFlags::values) ||
           (evaluation_flag & EvaluationFlags::gradients))
    {
      // slow path with FEValues
      Assert(fe_values.get() != nullptr,
             ExcMessage(
               "Not initialized. Please call FEPointEvaluation::reinit()!"));

      if (evaluation_flag & EvaluationFlags::values)
        {
          values.resize(unit_points.size());
          std::fill(values.begin(), values.end(), value_type());
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              const Number value = solution_values[i];
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        values[q], d) += fe_values->shape_value(i, q) * value;
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        values[q], d) +=
                      fe_values->shape_value_component(i, q, d) * value;
            }
        }

      if (evaluation_flag & EvaluationFlags::gradients)
        {
          gradients.resize(unit_points.size());
          std::fill(gradients.begin(), gradients.end(), gradient_type());
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              const Number value = solution_values[i];
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        gradients[q], d) += fe_values->shape_grad(i, q) * value;
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    internal::FEPointEvaluation::
                      EvaluatorTypeTraits<dim, n_components, Number>::access(
                        gradients[q], d) +=
                      fe_values->shape_grad_component(i, q, d) * value;
            }
        }
    }
}



template <int n_components, int dim, int spacedim, typename Number>
void
FEPointEvaluation<n_components, dim, spacedim, Number>::integrate(
  const ArrayView<Number> &               solution_values,
  const EvaluationFlags::EvaluationFlags &integration_flags)
{
  if (unit_points.size() == 0) // no evaluation points provided
    {
      std::fill(solution_values.begin(), solution_values.end(), 0.0);
      return;
    }

  AssertDimension(solution_values.size(), fe->dofs_per_cell);
  if (((integration_flags & EvaluationFlags::values) ||
       (integration_flags & EvaluationFlags::gradients)) &&
      !poly.empty())
    {
      // fast path with tensor product integration

      if (integration_flags & EvaluationFlags::values)
        AssertIndexRange(unit_points.size(), values.size() + 1);
      if (integration_flags & EvaluationFlags::gradients)
        AssertIndexRange(unit_points.size(), gradients.size() + 1);

      if (solution_renumbered_vectorized.size() != dofs_per_component)
        solution_renumbered_vectorized.resize(dofs_per_component);
      // zero content
      solution_renumbered_vectorized.fill(
        typename internal::FEPointEvaluation::EvaluatorTypeTraits<
          dim,
          n_components,
          VectorizedArray<Number>>::value_type());

      const std::size_t n_points = unit_points.size();
      const std::size_t n_lanes  = VectorizedArray<Number>::size();
      for (unsigned int i = 0; i < n_points; i += n_lanes)
        {
          // convert to vectorized format
          Point<dim, VectorizedArray<Number>> vectorized_points;
          for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
            for (unsigned int d = 0; d < dim; ++d)
              vectorized_points[d][j] = unit_points[i + j][d];

          typename internal::ProductTypeNoPoint<value_type,
                                                VectorizedArray<Number>>::type
            value = {};
          Tensor<1,
                 dim,
                 typename internal::ProductTypeNoPoint<
                   value_type,
                   VectorizedArray<Number>>::type>
            gradient;

          if (integration_flags & EvaluationFlags::values)
            for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
              internal::FEPointEvaluation::
                EvaluatorTypeTraits<dim, n_components, Number>::get_value(
                  value, j, values[i + j]);
          if (integration_flags & EvaluationFlags::gradients)
            for (unsigned int j = 0; j < n_lanes && i + j < n_points; ++j)
              {
                Assert(update_flags_mapping & update_inverse_jacobians,
                       ExcNotInitialized());
                gradients[i + j] =
                  apply_transformation(mapping_data.inverse_jacobians[i + j],
                                       gradients[i + j]);
                internal::FEPointEvaluation::
                  EvaluatorTypeTraits<dim, n_components, Number>::get_gradient(
                    gradient, j, gradients[i + j]);
              }

          // compute
          internal::integrate_add_tensor_product_value_and_gradient(
            poly,
            value,
            gradient,
            vectorized_points,
            solution_renumbered_vectorized);
        }

      // add between the lanes and write into the result
      std::fill(solution_values.begin(), solution_values.end(), Number());
      for (unsigned int comp = 0; comp < n_components; ++comp)
        for (unsigned int i = 0; i < dofs_per_component; ++i)
          {
            VectorizedArray<Number> result;
            internal::FEPointEvaluation::
              EvaluatorTypeTraits<dim, n_components, VectorizedArray<Number>>::
                write_value(result, comp, solution_renumbered_vectorized[i]);
            for (unsigned int lane = n_lanes / 2; lane > 0; lane /= 2)
              for (unsigned int j = 0; j < lane; ++j)
                result[j] += result[lane + j];
            solution_values[renumber[comp * dofs_per_component + i]] =
              result[0];
          }
    }
  else if ((integration_flags & EvaluationFlags::values) ||
           (integration_flags & EvaluationFlags::gradients))
    {
      // slow path with FEValues

      Assert(fe_values.get() != nullptr,
             ExcMessage(
               "Not initialized. Please call FEPointEvaluation::reinit()!"));
      std::fill(solution_values.begin(), solution_values.end(), 0.0);

      if (integration_flags & EvaluationFlags::values)
        {
          AssertIndexRange(unit_points.size(), values.size() + 1);
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values->shape_value(i, q) *
                      internal::FEPointEvaluation::
                        EvaluatorTypeTraits<dim, n_components, Number>::access(
                          values[q], d);
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values->shape_value_component(i, q, d) *
                      internal::FEPointEvaluation::
                        EvaluatorTypeTraits<dim, n_components, Number>::access(
                          values[q], d);
            }
        }

      if (integration_flags & EvaluationFlags::gradients)
        {
          AssertIndexRange(unit_points.size(), gradients.size() + 1);
          for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
            {
              for (unsigned int d = 0; d < n_components; ++d)
                if (nonzero_shape_function_component[i][d] &&
                    (fe->is_primitive(i) || fe->is_primitive()))
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values->shape_grad(i, q) *
                      internal::FEPointEvaluation::
                        EvaluatorTypeTraits<dim, n_components, Number>::access(
                          gradients[q], d);
                else if (nonzero_shape_function_component[i][d])
                  for (unsigned int q = 0; q < unit_points.size(); ++q)
                    solution_values[i] +=
                      fe_values->shape_grad_component(i, q, d) *
                      internal::FEPointEvaluation::
                        EvaluatorTypeTraits<dim, n_components, Number>::access(
                          gradients[q], d);
            }
        }
    }
}



template <int n_components, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluation<n_components, dim, spacedim, Number>::
  value_type &
  FEPointEvaluation<n_components, dim, spacedim, Number>::get_value(
    const unsigned int point_index) const
{
  AssertIndexRange(point_index, values.size());
  return values[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluation<n_components, dim, spacedim, Number>::
  gradient_type &
  FEPointEvaluation<n_components, dim, spacedim, Number>::get_gradient(
    const unsigned int point_index) const
{
  AssertIndexRange(point_index, gradients.size());
  return gradients[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline const typename FEPointEvaluation<n_components, dim, spacedim, Number>::
  gradient_type &
  FEPointEvaluation<n_components, dim, spacedim, Number>::get_unit_gradient(
    const unsigned int point_index) const
{
  Assert(!poly.empty(),
         ExcMessage("Unit gradients are currently only implemented for tensor "
                    "product finite elements combined with MappingQGeneric "
                    "mappings"));
  AssertIndexRange(point_index, unit_gradients.size());
  return unit_gradients[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::submit_value(
  const value_type & value,
  const unsigned int point_index)
{
  AssertIndexRange(point_index, unit_points.size());
  values[point_index] = value;
}



template <int n_components, int dim, int spacedim, typename Number>
inline void
FEPointEvaluation<n_components, dim, spacedim, Number>::submit_gradient(
  const gradient_type &gradient,
  const unsigned int   point_index)
{
  AssertIndexRange(point_index, unit_points.size());
  gradients[point_index] = gradient;
}



template <int n_components, int dim, int spacedim, typename Number>
inline DerivativeForm<1, dim, spacedim>
FEPointEvaluation<n_components, dim, spacedim, Number>::jacobian(
  const unsigned int point_index) const
{
  Assert(update_flags_mapping & update_jacobians, ExcNotInitialized());
  AssertIndexRange(point_index, mapping_data.jacobians.size());
  return mapping_data.jacobians[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline DerivativeForm<1, spacedim, dim>
FEPointEvaluation<n_components, dim, spacedim, Number>::inverse_jacobian(
  const unsigned int point_index) const
{
  Assert(update_flags_mapping & update_inverse_jacobians ||
           update_flags_mapping & update_gradients,
         ExcNotInitialized());
  AssertIndexRange(point_index, mapping_data.inverse_jacobians.size());
  return mapping_data.inverse_jacobians[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline Point<spacedim>
FEPointEvaluation<n_components, dim, spacedim, Number>::real_point(
  const unsigned int point_index) const
{
  Assert(update_flags_mapping & update_quadrature_points, ExcNotInitialized());
  AssertIndexRange(point_index, mapping_data.quadrature_points.size());
  return mapping_data.quadrature_points[point_index];
}



template <int n_components, int dim, int spacedim, typename Number>
inline Point<dim>
FEPointEvaluation<n_components, dim, spacedim, Number>::unit_point(
  const unsigned int point_index) const
{
  AssertIndexRange(point_index, unit_points.size());
  return unit_points[point_index];
}

DEAL_II_NAMESPACE_CLOSE

#endif


