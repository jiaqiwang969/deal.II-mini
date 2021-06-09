//include/deal.II-translator/base/auto_derivative_function_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#ifndef dealii_auto_derivative_function_h
#define dealii_auto_derivative_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类通过采用数字差商自动计算一个函数的梯度。这只是在用户函数没有自己提供梯度函数的情况下。
 * 下面是一个用户定义函数的例子，它只重载并实现了value()函数，但没有实现gradient()函数。如果梯度()函数被调用，那么由AutoDerivativeFunction实现的梯度函数将被调用，其中后者的函数采用了数字差商。
 *
 *
 * @code
 * class UserFunction: public AutoDerivativeFunction
 * {
 * // access to one component at one point
 * double value (const Point<dim>   &p,
 *               const unsigned int component = 0) const override
 * {
 *   // Implementation ....
 * };
 * };
 *
 * UserFunction user_function;
 *
 * // gradient by employing difference quotients.
 * Tensor<1,dim> grad=user_function.gradient(some_point);
 * @endcode
 *
 * 如果用户也重载并实现了梯度函数，那么，当然会调用用户的梯度函数。
 * 注意，上面解释的value()和gradient()函数的用法也适用于value_list()和gradient_list()函数，以及这些函数的向量值版本，例如vector_value(),
 * vector_gradient(), vector_value_list() 和 vector_gradient_list()。
 * gradient()和gradient_list()函数使用了 Function::value()
 * 函数。vector_gradient()和vector_gradient_list()使用的是
 * Function::vector_value()
 * 函数。请确保用户定义的函数分别实现value()函数和vector_value()函数。
 * 此外要注意，这个类的对象确实<b>not</b>代表了一个函数的导数，像FunctionDerivative，通过调用value()函数给出一个方向性的导数。事实上，这个类（AutoDerivativeFunction类）可以替代Function类作为用户定义类的基类。这个类实现了自动计算数字差商的梯度()函数，并作为基函数类和用户定义函数类之间的中间类。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim>
class AutoDerivativeFunction : public Function<dim>
{
public:
  /**
   * 差额公式的名称。
   *
   */
  enum DifferenceFormula
  {
    /**
     * 二阶的对称欧拉公式。    @f[
     * u'(t) \approx
     * \frac{u(t+h)
     *
     * -
     * u(t-h)}{2h}.
     * @f]
     *
     */
    Euler,
    /**
     * 一阶的上风欧拉公式。    @f[
     * u'(t) \approx
     * \frac{u(t)
     *
     * -
     * u(t-h)}{h}.
     * @f]
     *
     */
    UpwindEuler,
    /**
     * 四阶方案@f[
     * u'(t) \approx
     * \frac{u(t-2h)
     *
     * - 8u(t-h)
     * +  8u(t+h)
     *
     * - u(t+2h)}{12h}.
     * @f]。
     *
     */
    FourthOrder
  };

  /**
   * 构造函数。取差分步长<tt>h</tt>。在这里选择一个合适的值是用户的责任。<tt>h</tt>的选择应该考虑到绝对值以及函数的局部变化量。设置<tt>h=1e-6</tt>可能是绝对值为1左右的函数的一个好选择，而且变化不大。
   * <tt>h</tt>可以在以后使用set_h()函数来改变。
   * 设置DifferenceFormula
   * <tt>formula</tt>为set_formula()函数的默认<tt>Euler</tt>公式。通过调用set_formula()函数改变这个预设公式。
   *
   */
  AutoDerivativeFunction(const double       h,
                         const unsigned int n_components = 1,
                         const double       initial_time = 0.0);

  /**
   * 虚拟析构器；在这种情况下绝对必要。
   *
   */
  virtual ~AutoDerivativeFunction() override = default;

  /**
   * 选择差异公式。参见枚举#DifferenceFormula以了解可用的选择。
   *
   */
  void
  set_formula(const DifferenceFormula formula = Euler);

  /**
   * 取差值的步长<tt>h</tt>。在这里选择一个合适的值是用户的责任。<tt>h</tt>的选择应该考虑到函数的绝对值以及局部变化的量。设置<tt>h=1e-6</tt>可能是绝对值为1左右的函数的一个好选择，而且变化不大。
   *
   */
  void
  set_h(const double h);

  /**
   * 返回函数的指定分量在给定点的梯度。
   * 使用预设的#DifferenceFormula来计算数字差商。
   *
   */
  virtual Tensor<1, dim>
  gradient(const Point<dim> & p,
           const unsigned int component = 0) const override;

  /**
   * 返回函数的所有分量在给定点的梯度。
   * 使用预设的#DifferenceFormula来计算数字差商。
   *
   */
  virtual void
  vector_gradient(const Point<dim> &           p,
                  std::vector<Tensor<1, dim>> &gradients) const override;

  /**
   * 将<tt>gradients</tt>设为函数的指定分量在<tt>点</tt>的梯度。
   * 假设<tt>gradients</tt>已经有合适的大小，即与<tt>points</tt>数组的大小相同。
   * 使用预设的#DifferenceFormula来计算数字差商。
   *
   */
  virtual void
  gradient_list(const std::vector<Point<dim>> &points,
                std::vector<Tensor<1, dim>> &  gradients,
                const unsigned int             component = 0) const override;

  /**
   * 将<tt>gradients</tt>设置为函数在<tt>points</tt>处的梯度，适用于所有组件。假设<tt>gradients</tt>已经有正确的大小，即与<tt>points</tt>数组的大小相同。
   * <tt>gradients</tt>的外循环是在列表中的点上，内循环是在函数的不同成分上。
   * 使用预设的#DifferenceFormula来计算数字差商。
   *
   */
  virtual void
  vector_gradient_list(
    const std::vector<Point<dim>> &           points,
    std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;

  /**
   * 返回一个最小顺序为<tt>ord</tt>的#DifferenceFormula。
   *
   */
  static DifferenceFormula
  get_formula_of_order(const unsigned int ord);


private:
  /**
   * 差分公式的步骤大小。由set_h()函数设定。
   *
   */
  double h;

  /**
   * 包括按<tt>h</tt>缩放的单位向量。
   *
   */
  std::vector<Tensor<1, dim>> ht;

  /**
   * 差分公式。由set_formula()函数设置。
   *
   */
  DifferenceFormula formula;
};


DEAL_II_NAMESPACE_CLOSE

#endif


