//include/deal.II-translator/base/tensor_function_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_tensor_function_h
#define dealii_tensor_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_time.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类是一个张量值函数的模型。该类的接口与函数类的接口基本相同，只是不支持有多个分量的向量值函数，但返回类型总是张量值的。对这种类型的对象进行评估的返回值总是整个张量，而对于<tt>Function</tt>类，可以只要求一个特定的分量，或者使用<tt>vector_value</tt>函数，但它不返回值，而是将其写入其第二个参数提供的地址中。类的不同行为的原因是，在张量值函数的情况下，编译器预先知道参数的大小，这样就可以在堆栈中为返回值分配正确的内存量；另一方面，对于矢量值函数，编译器不知道其大小，所以内存必须在堆中分配，导致相对昂贵的复制操作。因此，我们可以认为这个类是<tt>Function</tt>类的特化，对于这个类，大小是已知的。另一个好处是可以返回任意等级的张量，而不仅仅是向量，因为它们的大小也可以简单地确定。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int rank, int dim, typename Number = double>
class TensorFunction
  : public FunctionTime<typename numbers::NumberTraits<Number>::real_type>,
    public Subscriptor
{
public:
  /**
   * <tt>value</tt>函数的返回类型的别名。
   *
   */
  using value_type = Tensor<rank, dim, Number>;

  /**
   * <tt>gradient</tt>函数的返回类型的别名。
   *
   */
  using gradient_type = Tensor<rank + 1, dim, Number>;

  /**
   * 用于表示时间的标量值实数类型。
   *
   */
  using time_type = typename FunctionTime<
    typename numbers::NumberTraits<Number>::real_type>::time_type;

  /**
   * 构造函数。可以为时间变量取一个初始值，默认为零。
   *
   */
  TensorFunction(const time_type initial_time = time_type(0.0));

  /**
   * 虚拟析构器；在这种情况下是绝对必要的，因为类通常不是通过它们的真实类型来使用，而是通过指向这个基类的指针。
   *
   */
  virtual ~TensorFunction() override = default;

  /**
   * 返回函数在给定点的值。
   *
   */
  virtual value_type
  value(const Point<dim> &p) const;

  /**
   * 将<tt>values</tt>设为函数在<tt>点</tt>的点值。
   * 假设<tt>values</tt>已经有合适的大小，即与<tt>points</tt>数组的大小相同。
   *
   */
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<value_type> &      values) const;

  /**
   * 返回函数在给定点的梯度。
   *
   */
  virtual gradient_type
  gradient(const Point<dim> &p) const;

  /**
   * 设置<tt>gradients</tt>为函数在<tt>points</tt>的梯度。
   * 假设<tt>values</tt>已经有合适的大小，即与<tt>points</tt>数组的大小相同。
   *
   */
  virtual void
  gradient_list(const std::vector<Point<dim>> &points,
                std::vector<gradient_type> &   gradients) const;
};



/**
 * 提供一个张量值的函数，它总是返回一个恒定的张量值。很明显，这个函数的所有导数都是零。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int rank, int dim, typename Number = double>
class ConstantTensorFunction : public TensorFunction<rank, dim, Number>
{
public:
  /**
   * 用于表示时间的标量值实数类型。
   *
   */
  using time_type = typename TensorFunction<rank, dim, Number>::time_type;

  /**
   * 构造函数；将常数张量值作为参数。参考值在内部被复制。
   * 可以指定时间变量的初始值，否则默认为零。
   *
   */
  ConstantTensorFunction(const dealii::Tensor<rank, dim, Number> &value,
                         const time_type initial_time = 0.0);

  virtual ~ConstantTensorFunction() override = default;

  virtual typename dealii::TensorFunction<rank, dim, Number>::value_type
  value(const Point<dim> &p) const override;

  virtual void
  value_list(
    const std::vector<Point<dim>> &points,
    std::vector<typename dealii::TensorFunction<rank, dim, Number>::value_type>
      &values) const override;

  virtual typename dealii::TensorFunction<rank, dim, Number>::gradient_type
  gradient(const Point<dim> &p) const override;

  virtual void
  gradient_list(
    const std::vector<Point<dim>> &points,
    std::vector<
      typename dealii::TensorFunction<rank, dim, Number>::gradient_type>
      &gradients) const override;

private:
  const dealii::Tensor<rank, dim, Number> _value;
};



/**
 * 提供一个张量值的函数，它总是返回0。很明显，这个函数的所有派生都是零。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int rank, int dim, typename Number = double>
class ZeroTensorFunction : public ConstantTensorFunction<rank, dim, Number>
{
public:
  /**
   * 用于表示时间的标量值实数类型。
   *
   */
  using time_type =
    typename ConstantTensorFunction<rank, dim, Number>::time_type;

  /**
   * 构造函数。
   * 可以指定时间变量的初始值，否则默认为零。
   *
   */
  ZeroTensorFunction(const time_type initial_time = 0.0);
};


DEAL_II_NAMESPACE_CLOSE

#endif


