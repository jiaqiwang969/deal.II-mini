//include/deal.II-translator/base/signaling_nan_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_signaling_nan_h
#define dealii_signaling_nan_h

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <limits>


DEAL_II_NAMESPACE_OPEN

namespace numbers
{
  namespace internal
  {
    /**
     * 一个命名空间，用于实现创建信令NaN对象的函数。这就是
     * numbers::signaling_nan() 函数调用的地方。
     *
     */
    namespace SignalingNaN
    {
      /**
       * 一个通用模板，用于知道如何用信号NaN来初始化 @p T
       * 类型的对象，以表示无效的值。
       * 这个类的真正实现发生在对模板参数 @p T.
       * 的特定值的（部分）专门化中。
       *
       */
      template <typename T>
      struct NaNInitializer;


      /**
       * 一般NaNInitializer类的特殊化，提供了一个函数，返回一个等于无效信号NaN的
       * @p float 值。
       *
       */
      template <>
      struct NaNInitializer<float>
      {
        static float
        invalid_element()
        {
          return std::numeric_limits<float>::signaling_NaN();
        }
      };


      /**
       * 一般NaNInitializer类的特殊化，它提供了一个函数，返回一个等于无效信号NaN的
       * @p double 值。
       *
       */
      template <>
      struct NaNInitializer<double>
      {
        static double
        invalid_element()
        {
          return std::numeric_limits<double>::signaling_NaN();
        }
      };


      /**
       * 一般NaNInitializer类的特殊化，它提供了一个函数，返回一个Tensor<1,dim>值，其组成部分是无效的信号NaN值。
       *
       */
      template <int dim, typename T>
      struct NaNInitializer<Tensor<1, dim, T>>
      {
        static Tensor<1, dim, T>
        invalid_element()
        {
          Tensor<1, dim, T> nan_tensor;

          for (unsigned int i = 0; i < dim; ++i)
            nan_tensor[i] = NaNInitializer<T>::invalid_element();

          return nan_tensor;
        }
      };



      /**
       * 一般NaNInitializer类的特殊化，它提供了一个函数来返回一个Tensor<rank,dim>值，其组件是无效的信号NaN值。
       *
       */
      template <int rank, int dim, typename T>
      struct NaNInitializer<Tensor<rank, dim, T>>
      {
        static Tensor<rank, dim, T>
        invalid_element()
        {
          Tensor<rank, dim, T> nan_tensor;

          // recursively initialize sub-tensors with invalid elements
          for (unsigned int i = 0; i < dim; ++i)
            nan_tensor[i] =
              NaNInitializer<Tensor<rank - 1, dim, T>>::invalid_element();

          return nan_tensor;
        }
      };



      /**
       * 一般NaNInitializer类的特殊化，它提供了一个函数来返回一个Tensor<rank,dim>值，其组件是无效的信号NaN值。
       *
       */
      template <int dim, typename T>
      struct NaNInitializer<Point<dim, T>>
      {
        static Point<dim, T>
        invalid_element()
        {
          Point<dim, T> nan_point;

          for (unsigned int i = 0; i < dim; ++i)
            nan_point[i] = NaNInitializer<T>::invalid_element();

          return nan_point;
        }
      };



      /**
       * 一般NaNInitializer类的特殊化，它提供了一个函数，返回一个SymmetricTensor<rank,dim>值，其组件是无效的信号NaN值。
       *
       */
      template <int rank, int dim, typename T>
      struct NaNInitializer<SymmetricTensor<rank, dim, T>>
      {
        static SymmetricTensor<rank, dim, T>
        invalid_element()
        {
          // initialize symmetric tensors via the unrolled list of elements
          T initializers
            [SymmetricTensor<rank, dim, T>::n_independent_components];
          for (unsigned int i = 0;
               i < SymmetricTensor<rank, dim, T>::n_independent_components;
               ++i)
            initializers[i] = NaNInitializer<T>::invalid_element();

          return SymmetricTensor<rank, dim, T>(initializers);
        }
      };



      /**
       * 一般NaNInitializer类的特殊化，它提供了一个函数，返回一个DerivativeForm<order,dim,spacedim>值，其组件是无效的信号NaN值。
       *
       */
      template <int order, int dim, int spacedim, typename T>
      struct NaNInitializer<DerivativeForm<order, dim, spacedim, T>>
      {
        static DerivativeForm<order, dim, spacedim, T>
        invalid_element()
        {
          DerivativeForm<order, dim, spacedim, T> form;

          // recursively initialize sub-tensors with invalid elements
          for (unsigned int i = 0; i < spacedim; ++i)
            form[i] = NaNInitializer<Tensor<order, dim, T>>::invalid_element();

          return form;
        }
      };
    } // namespace SignalingNaN
  }   // namespace internal



  /**
   * 提供一个类型为 @p T
   * 的对象，里面充满了信号NaN，在计算中使用时会引起异常。这些对象的内容是一个
   * "信号NaN"（"NaN "代表 "不是一个数字"，而 "信号
   * "意味着至少在支持这个的平台上，任何使用它们的算术运算都会终止程序）。这类对象的目的是把它们作为无效对象和数组的标记，这些对象和数组需要在某个时候被初始化为有效值，当这种后来的初始化没有在第一次使用之前发生时，就会触发错误。一个例子是像这样的代码。
   * @code
   * double x = numbers::signaling_nan<double>();
   * if (some condition)
   * {
   *   ...much code computing a,b,c...
   *   x = f(a,b,c);
   * }
   * else
   * {
   *   ...more code...
   *   // bug: we forgot to assign a value to 'x'.
   * }
   *
   * return std::sin(x);
   * @endcode
   * 错误在于`else`分支忘记向`x`变量写入一个值。如果你的平台支持信号NaN，那么这个错误在上面的最后一行会变得很明显，因为程序会被一个浮点异常终止。处理器意识到代码正试图对仍然存储在`x`中的信令NaN进行操作，并终止了程序，从而便于找到问题所在。如果在第一行将`x`初始化为零（或者完全不初始化），这就不是一个容易发现的错误。在这种情况下，最后一行对
   * `std::sin` 的调用将简单地计算 "某某
   * "的正弦，如果`某某条件==假`，但这种无效的结果对调用者来说可能并不明显，需要大量的调试才能发现，因为下游的计算将是错误的，没有任何指示*它们为什么是错误的。
   * @tparam  T
   * 返回无效对象的类型。这个类型可以是标量，也可以是Tensor、SymmetricTensor或DerivativeForm类型。如果
   * internal::SignalingNaN::NaNInitializer
   * 类对该类型有相应的专门化，可以支持其他类型。
   * @note  因为 @p T
   * 类型不作为函数参数使用，编译器不能从参数的类型中推断出它。因此，你必须明确地提供它。例如，这一行
   * @code
   *   Tensor<1,dim> tensor = Utilities::signaling_nan<Tensor<1,dim> >();
   * @endcode
   * 用无效的值初始化一个张量。
   *
   */
  template <class T>
  T
  signaling_nan()
  {
    // dispatch to the classes in the internal namespace because there
    // we can do partial specializations, which is not possible for
    // template functions such as the current one
    return internal::SignalingNaN::NaNInitializer<T>::invalid_element();
  }
} // namespace numbers


DEAL_II_NAMESPACE_CLOSE

#endif


