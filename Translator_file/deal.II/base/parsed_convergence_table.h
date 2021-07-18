//include/deal.II-translator/base/parsed_convergence_table_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_base_parsed_convergence_table_h
#define dealii_base_parsed_convergence_table_h

#include <deal.II/base/config.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_integrate_difference.h>

DEAL_II_NAMESPACE_OPEN

/**
 *
 * @brief The ParsedConvergenceTable class
 * 这个类简化了收敛表的构造，从参数文件中读取生成表格的选项。它提供了一系列的方法，可以用来计算给定参考精确解的误差，或者两个数值解之间的差异，或者任何其他自定义的误差计算，通过
 * std::function 对象给出。
 * 这个类的一个使用例子是这样给出的
 *
 * @code
 * ParsedConvergenceTable table;
 *
 * ParameterHandler prm;
 * table.add_parameters(prm);
 *
 * for (unsigned int i = 0; i < n_cycles; ++i)
 * {
 *   ... // do some computations
 *   table.error_from_exact(dof_handler, solution, exact_solution);
 * }
 * table.output_table(std::cout);
 * @endcode
 *
 * 上面的代码构造了一个ParsedConvergenceTable，它适用于标量问题，并将产生一个包含误差的`H1_norm`、`L2_norm`和`Linfty_norm`规范的误差表。
 * 每当调用error_from_exact()或difference()方法时，这个类的实例会检查其参数，计算构造时给出的参数所指定的所有规范，可能通过参数文件修改，计算所有使用add_extra_column()方法指定的额外列条目，并写入收敛表的一行。
 * 一旦你完成了计算，对output_table()的调用将在提供的流上生成一个格式化的收敛表，并将其发送到参数文件中指定的文件（如果有的话）。
 * 只要稍加修改，同样的代码就可以用来估计混合或多物理场问题的误差，例如。
 *
 * @code
 * using namespace VectorTools;
 * ParsedConvergenceTable table({"u,u,p"},{{H1_norm, L2_norm}, {L2_norm}});
 *
 * ParameterHandler prm;
 * table.add_parameters(prm);
 *
 * for (unsigned int i = 0; i < n_cycles; ++i)
 * {
 *   ... // do some computations
 *   table.error_from_exact(dof_handler, solution, exact_solution);
 * }
 * table.output_table(std::cout);
 * @endcode
 *
 * 上面的代码假设你正在解决一个有三个分量的斯托克斯问题。两个分量为矢量速度场`u`，一个分量为压力场`p`，并将产生一个误差表，其中包括速度场（前两个分量）的误差`H1`和`L2`规范以及压力场的`L2`误差。
 * 你也可以调用`table.output_table()`，不加参数，只把表写到参数文件中指定的文件。
 * 通过调用方法add_parameters()传递一个ParameterHandler对象，以下选项将被定义在给定的ParameterHandler对象中（在ParameterHandler对象的当前级别，即你用
 * ParamterHandler::enter_subsection()
 * 方法输入的任何级别），并可以在运行时通过参数文件进行修改。
 *
 * @code
 * set Enable computation of the errors = true
 * set Error file name                  =
 * set Error precision                  = 3
 * set Exponent for p-norms             = 2
 * set Extra columns                    = dofs, cells
 * set List of error norms to compute   = Linfty_norm, L2_norm, H1_norm
 * set Rate key                         = dofs
 * set Rate mode                        = reduction_rate_log2
 * @endcode
 *  在使用这个类时，请引证
 *
 * @code{.bib}
 * @article{SartoriGiulianiBardelloni-2018-a,
 * Author = {Sartori, Alberto and Giuliani, Nicola and
 *          Bardelloni, Mauro and Heltai, Luca},
 * Journal = {SoftwareX},
 * Pages = {318--327},
 * Title = {{deal2lkit: A toolkit library for high performance
 *          programming in deal.II}},
 * Doi = {10.1016/j.softx.2018.09.004},
 * Volume = {7},
 * Year = {2018}}
 * @endcode
 *
 *
 *
 */
class ParsedConvergenceTable
{
public:
  /**
   * ParsedConvergenceTable对象的最小构造函数。
   * 组件的数量必须与用于计算误差的有限元空间的组件数量一致。如果一个分量的名称重复，那么它将被解释为一个矢量场，重复的分量的误差将被归为一组。
   * 矢量 @p list_of_error_norms
   * 的大小必须与唯一分量名称的数量相匹配，并且可以包含零个或多个逗号分隔的标识符，用于计算每个分量的规范（关于可用的选项，请参见
   * VectorTools::NormType 的文档）。    例如，下面的构造函数
   * @code
   * using namespace VectorTools;
   * ParsedConvergenceTable table({"u", "v", "v"},
   *                            {{Linfty_norm}, {L2_norm, H1_norm}});
   * @endcode
   * 将产生（如果参数文件未被触动）一个类似于以下的表格
   * @code
   * cells dofs u_Linfty_norm    v_L2_norm      v_H1_norm
   * 4     9    1.183e-01
   *
   * -    5.156e-02
   *
   * -    2.615e-01
   *
   * -
   * 16    25   3.291e-02 2.50 1.333e-02 2.65 1.272e-01 1.41
   * 64    81   8.449e-03 2.31 3.360e-03 2.34 6.313e-02 1.19
   * 256   289  2.126e-03 2.17 8.418e-04 2.18 3.150e-02 1.09
   * 1024  1089 5.325e-04 2.09 2.106e-04 2.09 1.574e-02 1.05
   * @endcode
   * 请参阅其他构造函数，了解你可以改变的所有参数的文件。
   * @param  component_names 指定组件的名称；  @param
   * list_of_error_norms
   * 指定为每个独特的组件名称计算什么错误规范。
   *
   */
  ParsedConvergenceTable(
    const std::vector<std::string> &                    component_names = {"u"},
    const std::vector<std::set<VectorTools::NormType>> &list_of_error_norms = {
      {VectorTools::H1_norm, VectorTools::L2_norm, VectorTools::Linfty_norm}});

  /**
   * ParsedConvergenceTable的完整构造函数。      @param
   * component_names
   * 组件的名称。重复的连续名称被解释为矢量值域的组件；
   * @param  list_of_error_norms
   * 为每个独特的组件名称指定要计算的误差准则；  @param
   * exponent 在p-norms中使用的指数；  @param  extra_columns
   * 要添加的额外列。这些可以是 "单元格 "或 "道夫"；
   * @param  rate_key
   * 指定额外的列，我们将通过它来计算错误率。这个键可以是
   * "cell "或 "dofs
   * "中的一个，或者，如果你通过add_extra_column()方法向表中添加额外的列，它可以是你添加的额外列之一；
   * @param  rate_mode 指定计算错误率时要使用的比率模式。
   * 这可能是 "reduction_rate"，"reduction_rate_log2"，或者
   * "none"。参见 ConvergenceTable::RateMode
   * 的文档，了解每种模式的行为方式；  @param  error_file_name
   * 错误输出文件的名称（扩展名为txt, gpl, tex, 或
   * org）。如果与空字符串不同，比output_table()也会以从其扩展名推断出的格式写入该文件；
   * @param  precision 写入错误时要使用多少位数；  @param
   * compute_error
   * 控制是否启用填充表。这个标志可以用来在运行时禁用任何错误计算；你用这个构造函数指定的参数可以通过调用add_parameters()方法写入ParameterHandler对象中。一旦你调用add_parameters()方法，以下选项将被定义在给定的ParameterHandler对象中，并且这个类的实例的参数将跟随你在运行时对ParameterHandler对象的修改。
   * @code
   * # Listing of Parameters
   * #
   *
   * ---------------------
   * # When set to false, no computations are performed.
   * set Enable computation of the errors = true
   *
   * # Set this to a filename with extension .txt, .gpl, .org, or .tex to enable
   * # writing the convergence table to a file.
   * set Error file name                  =
   *
   * # Number of digits to use when printing the error.
   * set Error precision                  = 3
   *
   * # Extra columns to add to the table. Available options are dofs and cells.
   * set Extra columns                    = dofs, cells
   *
   * # The exponent to use when computing p-norms.
   * set Exponent for p-norms             = 2
   *
   * # Each component is separated by a semicolon and each norm by a comma. See
   * # the documentation of VectorTools::NormType for a list of implemented
   * # norms. If you want to skip a component, leave its entry empty.
   * set List of error norms to compute   = Linfty_norm, L2_norm, H1_norm
   *
   * # Key to use when computing convergence rates. If this is set to a
   * # column that is not present, or to none, then no error rates are computed.
   * set Rate key                         = dofs
   *
   * # What type of error rate to compute. Available options are
   * # reduction_rate_log2, reduction_rate, and none.
   * set Rate mode                        = reduction_rate_log2
   * @endcode
   *
   */
  ParsedConvergenceTable(
    const std::vector<std::string> &                    component_names,
    const std::vector<std::set<VectorTools::NormType>> &list_of_error_norms,
    const double                                        exponent,
    const std::set<std::string> &                       extra_columns,
    const std::string &                                 rate_key,
    const std::string &                                 rate_mode,
    const std::string &                                 error_file_name,
    const unsigned int                                  precision,
    const bool                                          compute_error);

  /**
   * 将这个类中的所有参数附加到参数处理程序的条目上  @p
   * prm.  每当 @p prm
   * 的内容发生变化，这个类的参数就会被更新。
   *
   */
  void
  add_parameters(ParameterHandler &prm);

  /**
   * 在错误表中添加一行，包含 @p solution 和 @p exact
   * 函数之间的错误，在参数文件中指定的规范（s）。
   * 如果你在这个调用中指定了一个 @p weight
   * 函数，那么这将用于计算加权误差。权重函数可以是一个标量函数（将用于所有组件），也可以是一个向量函数。
   * 当它是一个矢量函数时，如果分量的数量与底层有限元空间的分量数量不一致，将触发一个断言。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  error_from_exact(const DoFHandler<dim, spacedim> &vspace,
                   const VectorType &               solution,
                   const Function<spacedim> &       exact,
                   const Function<spacedim> *       weight = nullptr);

  /**
   * 和上面一样，有不同的映射。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  error_from_exact(const Mapping<dim, spacedim> &   mapping,
                   const DoFHandler<dim, spacedim> &vspace,
                   const VectorType &               solution,
                   const Function<spacedim> &       exact,
                   const Function<spacedim> *       weight = nullptr);

  /**
   * 通过在调用error_from_exact()或difference()时调用函数 @p
   * custom_function, ，在表中添加一个额外的列（名称为 @p
   * column_name) ）。
   * 你可以根据你的需要多次调用这个方法。如果 @p
   * column_name
   * 已经在之前的调用中使用过，那么以相同的名字调用这个方法将覆盖你之前指定的任何函数。如果你在这个调用中使用了一个lambda函数，请确保在lambda函数内部使用的变量在调用error_from_exact()或difference()之前一直有效。
   * 确保在第一次调用error_from_exact()或difference()之前添加所有额外列。在你已经开始填充收敛表之后再向收敛表添加额外的列将会触发一个异常。
   * 例如，这个方法可以用来计算时间步长的误差，例如。
   * @code
   * using namespace VectorTools;
   * ParsedConvergenceTable table({"u"}, {{L2_norm}});
   *
   * double dt = .5;
   * auto dt_function = [&]() {
   *      return dt;
   * };
   *
   * table.add_extra_column("dt", dt_function, false);
   *
   * for (unsigned int i = 0; i < n_cycles; ++i)
   * {
   *   // ... compute solution at the current dt
   *
   *   table.error_from_exact(dof_handler, solution, exact_solution);
   *   dt /= 2.0;
   * }
   * table.output_table(std::cout);
   * @endcode
   * 将产生一个类似于以下的表格
   * @code
   *  dt        u_L2_norm
   * 5.000e-1    5.156e-02
   *
   * -
   * 2.500e-2    1.333e-02 2.65
   * 1.250e-2    3.360e-03 2.34
   * 6.250e-3    8.418e-04 2.18
   * @endcode
   * 只要你使用以下参数文件（这里只显示非默认条目）。
   * @code
   * set Extra columns                  =
   * set List of error norms to compute = L2_norm
   * set Rate key                       = dt
   * @endcode
   * @param  column_name 要添加的列的名称；  @param  custom_function
   * %将被调用以填充给定条目的函数。你需要确保这个函数的范围在调用error_from_exact()或difference()之前是有效的；
   * @param  compute_rate
   * 如果设置为true，那么这个列将被包含在计算错误率的列的列表中。如果你想计算与这一列有关的错误率，你可能想把它设置为false。在这种情况下，你还应该在参数文件中指定
   * @p column_name 作为比率键。
   *
   */
  void
  add_extra_column(const std::string &            column_name,
                   const std::function<double()> &custom_function,
                   const bool                     compute_rate = true);

  /**
   * 同一矢量空间中两个解决方案之间的差异。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  difference(const DoFHandler<dim, spacedim> &,
             const VectorType &,
             const VectorType &,
             const Function<spacedim> *weight = nullptr);

  /**
   * 与上述相同，有一个非默认的映射。
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  difference(const Mapping<dim, spacedim> &mapping,
             const DoFHandler<dim, spacedim> &,
             const VectorType &,
             const VectorType &,
             const Function<spacedim> *weight = nullptr);

  /**
   * 将错误表写入 @p out
   * 流（文本格式），以及（可能）写入参数中指定的文件流（格式由文件名扩展名推断）。
   *
   */
  void
  output_table(std::ostream &out);

  /**
   * 将错误表写到参数中指定的文件流中。
   * 如果参数文件中的 "错误文件名
   * "选项被设置为空字符串，则不写输出。
   *
   */
  void
  output_table();

private:
  /**
   * 向输出表添加速率。
   *
   */
  void
  prepare_table_for_output();

  /**
   * 解决方案组件的名称。
   *
   */
  const std::vector<std::string> component_names;

  /**
   * 与上述相同，但只包含一次重复的成分名称。
   *
   */
  const std::vector<std::string> unique_component_names;

  /**
   * 每个独特组件名称的掩码。
   *
   */
  const std::vector<ComponentMask> unique_component_masks;

  /**
   * 向表中添加行时要调用的额外方法。
   *
   */
  std::map<std::string, std::pair<std::function<double()>, bool>>
    extra_column_functions;

  /**
   * 计算每个组件的错误类型。
   *
   */
  std::vector<std::set<VectorTools::NormType>> norms_per_unique_component;

  /**
   * 在p-norm类型中使用的指数。
   *
   */
  double exponent;

  /**
   * 实际的表格
   *
   */
  ConvergenceTable table;

  /**
   * 要添加到表中的额外列。
   *
   */
  std::set<std::string> extra_columns;

  /**
   * 我们计算收敛率所涉及的列的名称。
   *
   */
  std::string rate_key;

  /**
   * 收敛率模式。参见 ConvergenceTable::RateMode 的文件。
   *
   */
  std::string rate_mode;

  /**
   * 用于输出表格的精度。
   *
   */
  unsigned int precision;

  /**
   * 写入文件时要使用的文件名。
   *
   */
  std::string error_file_name;

  /**
   * 计算错误。如果这是错误的，所有进行错误计算的方法都被禁用，不做任何事情。
   *
   */
  bool compute_error;
};



#ifndef DOXYGEN
// ============================================================
// Template functions
// ============================================================
template <int dim, int spacedim, typename VectorType>
void
ParsedConvergenceTable::difference(const DoFHandler<dim, spacedim> &dh,
                                   const VectorType &               solution1,
                                   const VectorType &               solution2,
                                   const Function<spacedim> *       weight)
{
  AssertThrow(solution1.size() == solution2.size(),
              ExcDimensionMismatch(solution1.size(), solution2.size()));
  VectorType solution(solution1);
  solution -= solution2;
  error_from_exact(
    dh,
    solution,
    Functions::ConstantFunction<spacedim>(0, component_names.size()),
    weight);
}



template <int dim, int spacedim, typename VectorType>
void
ParsedConvergenceTable::difference(const Mapping<dim, spacedim> &   mapping,
                                   const DoFHandler<dim, spacedim> &dh,
                                   const VectorType &               solution1,
                                   const VectorType &               solution2,
                                   const Function<spacedim> *       weight)
{
  AssertThrow(solution1.size() == solution2.size(),
              ExcDimensionMismatch(solution1.size(), solution2.size()));
  VectorType solution(solution1);
  solution -= solution2;
  error_from_exact(
    mapping,
    dh,
    solution,
    Functions::ConstantFunction<spacedim>(0, component_names.size()),
    weight);
}



template <int dim, int spacedim, typename VectorType>
void
ParsedConvergenceTable::error_from_exact(const DoFHandler<dim, spacedim> &dh,
                                         const VectorType &        solution,
                                         const Function<spacedim> &exact,
                                         const Function<spacedim> *weight)
{
  error_from_exact(get_default_linear_mapping(dh.get_triangulation()),
                   dh,
                   solution,
                   exact,
                   weight);
}



template <int dim, int spacedim, typename VectorType>
void
ParsedConvergenceTable::error_from_exact(const Mapping<dim, spacedim> &mapping,
                                         const DoFHandler<dim, spacedim> &dh,
                                         const VectorType &        solution,
                                         const Function<spacedim> &exact,
                                         const Function<spacedim> *weight)
{
  const auto n_components = component_names.size();

  if (compute_error)
    {
      AssertDimension(exact.n_components, n_components);
      AssertDimension(dh.get_fe().n_components(), n_components);

      const types::global_cell_index n_active_cells =
        dh.get_triangulation().n_global_active_cells();
      const unsigned int n_dofs = dh.n_dofs();

      for (const auto &col : extra_columns)
        if (col == "cells")
          {
            table.add_value("cells", n_active_cells);
            table.set_tex_caption("cells", "\\# cells");
            table.set_tex_format("cells", "r");
          }
        else if (col == "dofs")
          {
            table.add_value("dofs", n_dofs);
            table.set_tex_caption("dofs", "\\# dofs");
            table.set_tex_format("dofs", "r");
          }

      // A vector of zero std::functions with n_components components
      const std::vector<std::function<double(const Point<spacedim> &)>>
        zero_components(n_components,
                        [](const Point<spacedim> &) { return 0.0; });

      // The default weight function, with n_components components
      std::vector<std::function<double(const Point<spacedim> &)>>
        weight_components(n_components,
                          [](const Point<spacedim> &) { return 1.0; });

      if (weight != nullptr)
        {
          if (weight->n_components == 1)
            {
              for (auto &f : weight_components)
                f = [&](const Point<spacedim> &p) { return weight->value(p); };
            }
          else
            {
              AssertDimension(weight->n_components, n_components);
              for (unsigned int i = 0; i < n_components; ++i)
                weight_components[i] = [&](const Point<spacedim> &p) {
                  return weight->value(p, i);
                };
            }
        }

      for (unsigned int i = 0; i < norms_per_unique_component.size(); ++i)
        {
          std::map<VectorTools::NormType, double> errors;

          const auto &norms = norms_per_unique_component[i];
          const auto &mask  = unique_component_masks[i];

          // Simple case first
          if (norms.empty())
            continue;

          auto components_expr = zero_components;
          for (unsigned int i = 0; i < n_components; ++i)
            if (mask[i] == true)
              components_expr[i] = weight_components[i];

          FunctionFromFunctionObjects<spacedim> select_component(
            components_expr);

          Vector<float> difference_per_cell(
            dh.get_triangulation().n_global_active_cells());

          QGauss<dim> q_gauss((dh.get_fe().degree + 1) * 2);

          for (const auto &norm : norms)
            {
              difference_per_cell = 0;
              VectorTools::integrate_difference(mapping,
                                                dh,
                                                solution,
                                                exact,
                                                difference_per_cell,
                                                q_gauss,
                                                norm,
                                                &select_component,
                                                exponent);

              errors[norm] = VectorTools::compute_global_error(
                dh.get_triangulation(), difference_per_cell, norm, exponent);

              std::string name = unique_component_names[i] + "_" +
                                 Patterns::Tools::to_string(norm);
              std::string latex_name = "$\\| " + unique_component_names[i] +
                                       " - " + unique_component_names[i] +
                                       "_h \\|_{" +
                                       Patterns::Tools::to_string(norm) + "}$";

              table.add_value(name, errors[norm]);
              table.set_precision(name, precision);
              table.set_scientific(name, true);
              table.set_tex_caption(name, latex_name);
            }
        }

      for (const auto &extra_col : extra_column_functions)
        {
          const double custom_error = extra_col.second.first();

          std::string name = extra_col.first;
          table.add_value(name, custom_error);
          table.set_precision(name, precision);
          table.set_scientific(name, true);
        }
    }
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif


