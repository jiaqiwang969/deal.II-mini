//include/deal.II-translator/base/convergence_table_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_convergence_table_h
#define dealii_convergence_table_h


#include <deal.II/base/config.h>

#include <deal.II/base/table_handler.h>

DEAL_II_NAMESPACE_OPEN


/**
 * ConvergenceTable类是TableHandler类的一个应用，它存储了一些收敛数据，比如cg方法的残差，或者一些离散解的评估<i>L<sup>2</sup></i>-错误等，并评估收敛率或顺序。
 * 已经实现的#RateMode是#reduction_rate，其中收敛率是下面两行的商，以及#reduction_rate_log2，它评估收敛的顺序。这些标准评估对全局细化很有用，对于局部细化来说，这可能不是一个合适的方法，因为收敛率的设置应该与单元数或DoF数有关。这些非标准方法的实现是留给用户的。
 * 例如，可以通过调用`add_value("n cells",
 * n_cells)`将单元格的数量添加到表中。DoF的数量可以通过调用`add_value("n
 * dofs",
 * n_dofs)`添加到表中。当然，我们也可以通过调用add_value()和其他参数来添加更多种类的信息。在任何情况下，在输出表格之前，可以调用函数evaluate_convergence_rates()和evaluate_all_convergence_rates()。
 * 对于如何评估同一RateMode中多列的收敛率，有两种可能性。  <ol>   <li>  对所有需要的列调用evaluate_convergence_rates()  <li>  对所有不需要这种评估的列调用omit_column_from_convergence_rate_evaluation() ，然后evaluate_all_convergence_rates() 来评估所有没有被标记为忽略的列的收敛率。  </ol>
 * 在 step-7 和 step-13
 * 的示例程序中也可以找到关于这个类的详细讨论。它也被用在
 * step-74  中。
 *
 *
 * @ingroup textoutput
 *
 *
 */
class ConvergenceTable : public TableHandler
{
public:
  /**
   * 构建器。
   *
   */
  ConvergenceTable() = default;

  /**
   * 相对于行的比率。
   *
   */
  enum RateMode
  {
    /**
     * 不做任何事情。
     *
     */
    none,
    /**
     * 前一行和这一行的值的商数。
     *
     */
    reduction_rate,
    /**
     * 以2为底的#reduction_rate的对数，代表网格大小减半时的收敛顺序，例如，从h到h/2。
     *
     */
    reduction_rate_log2
  };

  /**
   * 评估数据列<tt>data_column_key</tt>的收敛率，由于#RateMode与参考列<tt>reference_column_key</tt>的关系。要确保数据列和参考数据列的表项的值类型是一个数字，即double,
   * float, (unsigned) int，等等。
   * 由于这个类没有关于参考列与值列所依据的空间维度的信息，所以需要作为最后一个参数传递给这个方法。<i>default
   * dimension for the reference
   * column</i>是2，这对于2D中的单元格数量是合适的。如果参考列是
   * $1/h$
   * ，记得在三维中工作时也要将维数设置为1，以获得正确的速率。
   * 新的速率列和数据列将被合并为一个超级列。
   * 超级列的文本标题将（默认）与数据列的标题相同。这可以通过使用基类TableHandler的<tt>set_tex_supercaption
   * (..)</tt>函数来改变。    这个方法的行为方式是这样的。
   * 如果RateMode是reduction_rate，那么计算出来的输出是 $
   * \frac{e_{n-1}/k_{n-1}}{e_n/k_n}, $ ，其中 $k$
   * 是参考列（没有维度依赖！）。
   * 如果RateMode是reduction_rate_log2，那么计算的输出结果是 $
   * dim \frac{\log |e_{n-1}/e_{n}|}{\log |k_n/k_{n-1}|} $  。
   * 这很有用，例如，如果我们使用自由度的数量作为参考键，或者更好的是使用单元格的数量。
   * 假设二维的误差与 $ C (1/\sqrt{k})^r $
   * 成正比，那么这个方法将产生 $r$
   * 的结果。对于一般的维度，如这个函数的最后一个参数所描述的，公式需要是
   * $ C (1/\sqrt[dim]{k})^r $  。
   * @note
   * 因为这个函数是在几行已经被填满之后才向表格添加列的，所以它关闭了TableHandler基类的自动填充模式。如果你打算用自动填充的方式进一步添加数据，你必须在调用这个函数后重新启用它。
   *
   */
  void
  evaluate_convergence_rates(const std::string &data_column_key,
                             const std::string &reference_column_key,
                             const RateMode     rate_mode,
                             const unsigned int dim = 2);


  /**
   * 评估数据列<tt>data_column_key</tt>由于#RateMode的收敛率。
   * 要确保数据列的表项的值类型是一个数字，即双数、浮点数、（无符号）int，等等。
   * 新的速率列和数据列将被合并为一个超级列。
   * 超级列的文本标题将（默认）与数据列的标题相同。这可以通过使用基类TableHandler的set_tex_supercaption()函数来改变。
   * @note
   * 由于这个函数是在几行已经被填满之后才向表添加列的，所以它关闭了TableHandler基类的自动填充模式。如果你打算用自动填充的方式进一步添加数据，你必须在调用这个函数后重新启用它。
   *
   */
  void
  evaluate_convergence_rates(const std::string &data_column_key,
                             const RateMode     rate_mode);

  /**
   * 在评估 "所有
   * "列的收敛率时省略此列<tt>key</tt>（不是超级列！）（见下面两个函数）。
   * Column::flag==1 是为从收敛率评估中省略该列而保留的。
   *
   */
  void
  omit_column_from_convergence_rate_evaluation(const std::string &key);

  /**
   * 评估由于<tt>rate_mode</tt>与参考列<tt>reference_column_key</tt>有关的收敛率。这个函数评估所有列的速率，除了要省略的列（见前面的函数）和以前评估过的速率列之外。
   * 这个函数允许评估一个表中几乎所有列的收敛率，而不需要为每一列单独调用evaluate_convergence_rates()。
   * 例子。像<tt>n cells</tt>或<tt>n
   * dofs</tt>列可能希望在评估收敛率时被省略。因此他们应该通过调用
   * omit_column_from_convergence_rate_evaluation()来省略。
   *
   */
  void
  evaluate_all_convergence_rates(const std::string &reference_column_key,
                                 const RateMode     rate_mode);

  /**
   * 评估由于<tt>rate_mode</tt>引起的收敛率。这个函数评估所有列的速率，除了将被省略的列（见前一个函数）和先前评估的速率列之外。
   * 这个函数允许评估一个表中几乎所有列的收敛率，而不需要为每一列单独调用evaluate_convergence_rates()。
   * 例子。像<tt>n cells</tt>或<tt>n
   * dofs</tt>列可能希望在评估收敛率时被省略。因此他们应该通过调用
   * omit_column_from_convergence_rate_evaluation()来省略。
   *
   */
  void
  evaluate_all_convergence_rates(const RateMode rate_mode);

  /**
   * @addtogroup  Exceptions
   *  
     * @{ 
   *
   */

  /**
   * Exceptions
   *
   */
  DeclException1(ExcRateColumnAlreadyExists,
                 std::string,
                 << "Rate column <" << arg1 << "> does already exist.");
  //@}
};


DEAL_II_NAMESPACE_CLOSE

#endif


