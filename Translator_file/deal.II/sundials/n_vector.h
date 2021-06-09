//include/deal.II-translator/sundials/n_vector_0.txt
//-----------------------------------------------------------
//
//    Copyright (C) 2020 - 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#ifndef dealii_sundials_n_vector_h
#define dealii_sundials_n_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  include <sundials/sundials_nvector.h>

#  include <functional>
#  include <memory>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
namespace SUNDIALS
{
  namespace internal
  {
    template <typename VectorType>
    class NVectorView;
  }
} // namespace SUNDIALS
#  endif

namespace SUNDIALS
{
  namespace internal
  {
    /**
     * 创建一个给定的NVectorView  @p vector.
     * 这个调用的目的是用来作为
     * @code
     * auto view = make_nvector_view(vector);
     * @endcode
     * 产生的对象`view`必须被保留在周围，只要任何其他对象将使用内部查看的N_Vector。
     * @tparam  VectorType
     * 被观察的向量的类型。这个参数可以自动推导，并将尊重一个潜在的const-qualifier。
     * @param  矢量 要查看的矢量。      @return   @p vector.
     * 的NVectorView  @related  NVectorView
     *
     */
    template <typename VectorType>
    NVectorView<VectorType>
    make_nvector_view(VectorType &vector);

    /**
     * 检索连接到N_Vector的底层向量  @p v.
     * 只有当底层向量不是常量时，这个调用才会成功。在这种情况下，请使用unwrap_nvector_const()。
     * @note
     * 用户必须确保在调用此函数时询问正确的VectorType，并且没有类型安全检查。
     * @tparam  VectorType 存储在 @p v 中的向量的类型  @param  v
     * 要解包的向量  @return  存储在 @p v 中的向量
     *
     */
    template <typename VectorType>
    VectorType *
    unwrap_nvector(N_Vector v);

    /**
     * 检索连接到N_Vector  @p v
     * 的底层向量，作为一个常数指针。
     * @note
     * 用户必须确保在调用此函数时询问正确的VectorType，并且没有类型安全检查。
     * @tparam  VectorType 存储在 @p v 中的向量的类型  @param  v
     * 要解包的向量  @return  存储在 @p v 中的向量
     *
     */
    template <typename VectorType>
    const VectorType *
    unwrap_nvector_const(N_Vector v);

    /**
     * 一个指向向量的视图，只要需要N_Vector就可以使用。
     * 这个类的对象最好是通过make_nvector_view()创建，因为
     * @code
     * auto view = make_nvector_view(vector);
     * @endcode
     * 产生的N_Vector是实际矢量的视图，而不是拥有内存。另外，不能对生成的N_Vector调用N_VDestroy()，因为这将导致在析构器中出现双倍删除。
     * @note
     * SUNDIALS永远不会在它自己没有创建的向量上调用N_VDestroy()，因此上述约束没有限制用户。
     * @tparam  VectorType 储存在 @p v 的向量的类型。
     *
     */
    template <typename VectorType>
    class NVectorView
    {
    public:
      /**
       * 默认构造函数。
       * 该对象实际上没有查看任何东西，需要用operator=(NVectorView
       * &&)来分配。
       *
       */
      NVectorView() = default;

      /**
       * 构造函数。创建 @p vector. 的视图。
       *
       */
      NVectorView(VectorType &vector);

      /**
       * 移动赋值。
       *
       */
      NVectorView(NVectorView &&) noexcept = default;

      /**
       * 移动构造器。
       *
       */
      NVectorView &
      operator=(NVectorView &&) noexcept = default;

      /**
       * 明确删除复制的ctor。这个类只有移动。
       *
       */
      NVectorView(const NVectorView &) = delete;

      /**
       * 明确删除复制赋值。这个类是只允许移动的。
       *
       */
      NVectorView &
      operator=(const NVectorView &) = delete;

      /**
       * 解构器。
       * @note  这将不会破坏被查看的向量。
       *
       */
      ~NVectorView() = default;

      /**
       * 隐式转换为N_Vector。这个操作符使NVectorView看起来像一个实际的N_Vector，它可以直接作为许多SUNDIALS函数的参数使用。
       *
       */
      operator N_Vector() const;

      /**
       * 访问这个对象所查看的N_Vector。
       *
       */
      N_Vector operator->() const;

    private:
      /**
       * 实际指向该类所查看的向量的指针。
       *
       */
      std::unique_ptr<_generic_N_Vector, std::function<void(N_Vector)>>
        vector_ptr;
    };
  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE

#endif
#endif


