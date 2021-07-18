//include/deal.II-translator/base/array_view_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
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

#ifndef dealii_array_view_h
#define dealii_array_view_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
template <int N, typename T>
class Table;

template <typename number>
class LAPACKFullMatrix;


/**
 * 一个表示 @p ElementType
 * 类型的内存位置窗口的类，并将其作为一个数组来呈现。从本质上讲，这个类只不过是一个指向第一个位置的指针和一个代表数组长度的整数元素而已。内存仍然属于分配它的人，因为这个类并没有接管所有权。
 * 使用这个类的好处是，你不需要传递成对的指针，而且
 * <code>operator[]</code>
 * 会检查你下标这个数组视图的索引是否有效。注意，只有当底层数据存储在CPU内存中时，才允许访问元素。
 * 这个类可以处理对非常量和常量内存位置的视图。如果你想表示一个常数的视图，那么这个类的模板参数类型也需要是
 * @p const 。下面的代码片断给出了一个例子。
 *
 * @code
 * std::vector<int> array = get_data(); // a writable array
 * ArrayView<int> view (&array[5], 5); // a view of elements 5..9 (inclusive)
 * view[2] = 42; // array[7] is set to 42
 * ArrayView<const int> const_view (&array[5], 5); // same view, but read-only
 * int element_7 = const_view[2]; // set element_7 to 42
 * const_view[2] = 42; // this line won't compile; can't write into this view
 * @endcode
 * 在任何一种情况下，访问一个视图的元素都不会改变ArrayView对象本身，因此
 * ArrayView::operator[] 是一个 @p
 * 的const函数。这对应于这样一个概念：视图只是代表了一个，嗯，由别人拥有的内存的
 * "视图"。因此，访问视图中的元素会改变其他对象所管理的内存，但不会改变视图本身，这使得我们可以将
 * ArrayView::operator[] 变成 @p const 成员函数。这与 std::vector,
 * 相反， std::vector, 管理着它所指向的内存，因此改变
 * std::vector 的一个元素会改变 std::vector 对象本身
 *
 * - 因此， std::vector::operator[] 就不是 @p const. 。
 *
 *
 * @note  这个类与
 * [`std::span`](https://en.cppreference.com/w/cpp/container/span),
 * 类似，但后者只在C++20中开始使用。
 *
 *
 * @ingroup data
 *
 */
template <typename ElementType, typename MemorySpaceType = MemorySpace::Host>
class ArrayView
{
public:
  /**
   * 一个别名，表示这个容器类的 "value_type"，即它 "存储
   * "或指向的元素的类型。
   *
   */
  using value_type = ElementType;

  /**
   * 指向数组的迭代器的一个别名。
   *
   */
  using iterator = value_type *;

  /**
   * 指向数组的常数迭代器的别名。
   *
   */
  using const_iterator = const ElementType *;

  /**
   * 默认的构造函数。
   *
   */
  ArrayView();

  /**
   * 构造函数。      @param[in]  starting_element
   * 指向该对象所代表的数组的第一个元素的指针。
   * @param[in]  n_elements
   * 这个对象应该代表的内存块的长度（元素）。
   * @note
   * 由这些参数构建的对象不知道它所指向的对象到底有多大。因此，每当你调用
   * ArrayView::operator[],
   * 时，数组视图可以检查给定的索引是否在视图的范围内，但它不能检查视图是否确实是分配该范围的底层对象的有效元素范围的子集。换句话说，你需要确保这个构造函数的两个参数所指定的视图范围实际上是它所指向的数组元素的一个子集。做到这一点的适当方法是使用make_array_view()函数。
   *
   */
  ArrayView(value_type *starting_element, const std::size_t n_elements);

  /**
   * 指向非 @p const
   * 元素的数组视图的复制构造函数。如果当前对象将指向非
   * @p const
   * 元素，那么这就是一个直接的复制构造函数。另一方面，如果当前类型的
   * @p ElementType 模板参数是一个 @p const
   * 限定类型，那么当前构造函数是一个转换构造函数，将非
   * @p const 视图转换为 @p const 视图，类似于将非 @p const
   * 指针转换为 @p const 指针。
   *
   */
  ArrayView(const ArrayView<typename std::remove_cv<value_type>::type,
                            MemorySpaceType> &view);

  /**
   * 一个构造函数，从一个value_type对象自动创建一个视图。这样创建的视图的长度为1。
   *
   */
  explicit ArrayView(value_type &element);

  /**
   * 一个构造函数，从一个 std::vector
   * 对象中自动创建一个视图。
   * 该视图包含了给定向量的所有元素。
   * 当调用一个以ArrayView对象为参数的函数，并传入一个
   * std::vector. 时，这个隐式转换构造函数特别有用。
   * @note  这个构造函数需要一个 @p const
   * 向量的引用作为参数。    它只能用于初始化指向 @p const
   * 内存位置的ArrayView对象，例如 <code>ArrayView@<const
   * double@></code>  。    不能用这样的参数初始化指向非 @p
   * const 内存的ArrayView对象，如 <code>ArrayView@<double@></code>
   * 。
   *
   */
  ArrayView(
    const std::vector<typename std::remove_cv<value_type>::type> &vector);

  /**
   * 一个构造函数，从一个 std::vector
   * 对象中自动创建一个视图。
   * 该视图包含了给定向量的所有元素。
   * 当调用一个以ArrayView对象为参数的函数，并传入一个
   * std::vector. 时，这个隐式转换构造函数就特别有用。
   * @note  这个构造函数需要一个非 @p const
   * 向量的引用作为参数。它可以用来初始化ArrayView对象，该对象指向
   * @p const 内存位置，例如 <code>ArrayView@<const double@></code>
   * ，或者指向非 @p const 内存，例如
   * <code>ArrayView@<double@></code>  。
   *
   */
  ArrayView(std::vector<typename std::remove_cv<value_type>::type> &vector);

  /**
   * 一个构造函数，为一个给定的C-style数组自动创建一个视图。
   * 这个构造函数可以如下使用。
   * @code
   * ArrayView<const int>
   * get_data_table ()
   * {
   *   static const int my_data[7] = { 1, 1, 2, 3, 5, 8, 13 };
   *   return {my_data};
   * }
   * @endcode
   * 这样返回的对象是一个数组的视图，其大小被正确推导出来。
   *
   */
  template <std::size_t N>
  ArrayView(value_type (&array)[N]);

  /**
   * 一个构造函数可以从一个 std::array
   * 对象中自动创建一个视图。
   * 该视图包含了给定向量的所有元素。
   * 当调用一个以ArrayView对象为参数的函数，并传入一个
   * std::array. 时，这个隐式转换构造函数特别有用。
   *
   */
  template <std::size_t N>
  ArrayView(
    const std::array<typename std::remove_cv<value_type>::type, N> &vector);

  /**
   * 一个构造函数，可以从一个 std::array
   * 对象自动创建一个视图。
   * 该视图包含了给定向量的所有元素。
   * 当调用一个以ArrayView对象为参数的函数，并传入一个
   * std::array. 时，这个隐式转换构造函数就特别有用。
   *
   */
  template <std::size_t N>
  ArrayView(std::array<typename std::remove_cv<value_type>::type, N> &vector);

  /**
   * 重新初始化一个视图。      @param[in]  starting_element
   * 指向该对象所代表的数组的第一个元素的指针。
   * @param[in]  n_elements
   * 这个对象应该代表的内存块的长度（以元素计）。
   * @note
   * 由这些参数构建的对象不知道它所指向的对象到底有多大。因此，每当你调用
   * ArrayView::operator[],
   * 时，数组视图可以检查给定的索引是否在视图的范围内，但它不能检查视图是否确实是分配该范围的底层对象的有效元素范围的子集。换句话说，你需要确保这个构造函数的两个参数所指定的视图范围实际上是它所指向的数组元素的一个子集。做到这一点的适当方法是使用make_array_view()函数。
   *
   */
  void
  reinit(value_type *starting_element, const std::size_t n_elements);

  /**
   * 比较两个相同类型的ArrayView对象。如果两个对象具有相同的大小和相同的起始指针，则认为它们是相等的。
   * 这个版本总是与const value_type进行比较。
   *
   */
  bool
  operator==(
    const ArrayView<const value_type, MemorySpaceType> &other_view) const;

  /**
   * 比较两个相同类型的ArrayView对象。如果两个对象具有相同的大小和相同的起始指针，则认为它们是相等的。
   * 这个版本总是与非const value_type进行比较。
   *
   */
  bool
  operator==(const ArrayView<typename std::remove_cv<value_type>::type,
                             MemorySpaceType> &other_view) const;

  /**
   * 比较两个相同类型的ArrayView对象。如果两个对象具有相同的大小和相同的起始指针，则认为它们是相等的。
   * 这个版本总是与const value_type进行比较。
   *
   */
  bool
  operator!=(
    const ArrayView<const value_type, MemorySpaceType> &other_view) const;

  /**
   * 比较两个相同类型的ArrayView对象。如果两个对象具有相同的大小和相同的起始指针，则认为它们是相等的。
   * 这个版本总是与非const value_type进行比较。
   *
   */
  bool
  operator!=(const ArrayView<typename std::remove_cv<value_type>::type,
                             MemorySpaceType> &other_view) const;

  /**
   * 返回这个对象所代表的内存视图的大小（以元素为单位）。
   *
   */
  std::size_t
  size() const;

  /**
   * 返回一个指向作为元素存储的底层数组的指针。
   * 如果该容器是空的，将返回一个nullptr。
   *
   */
  value_type *
  data() const noexcept;

  /**
   * 返回一个指向数组视图开头的迭代器。
   *
   */
  iterator
  begin() const;

  /**
   * 返回一个指向数组视图结束后的迭代器。
   *
   */
  iterator
  end() const;

  /**
   * 返回一个指向数组视图开始的常数迭代器。
   *
   */
  const_iterator
  cbegin() const;

  /**
   * 返回一个指向数组视图结束后的常数迭代器。
   *
   */
  const_iterator
  cend() const;

  /**
   * 返回一个指向当前对象所代表的范围内第 $i$
   * 个元素的引用。    这个函数被标记为 @p const
   * ，因为它不改变 <em> 视图对象  </em>
   * 。但是它可能会返回一个非  @p const
   * 内存位置的引用，这取决于类的模板类型是否是  @p
   * const。
   * 只有当底层数据确实存储在CPU内存中时，才允许调用这个函数。
   *
   */
  value_type &operator[](const std::size_t i) const;

private:
  /**
   * 指向该对象所代表的内存中的位置范围的第一个元素的指针。
   *
   */
  value_type *starting_element;

  /**
   * 这个对象所代表的数组的长度。
   *
   */
  std::size_t n_elements;

  friend class ArrayView<const ElementType, MemorySpaceType>;
};



//---------------------------------------------------------------------------


namespace internal
{
  namespace ArrayViewHelper
  {
    template <typename MemorySpaceType>
    inline bool
    is_in_correct_memory_space(const void *const ptr)
    {
#ifndef DEAL_II_COMPILER_CUDA_AWARE
      (void)ptr;
      static_assert(std::is_same<MemorySpaceType, MemorySpace::Host>::value,
                    "If the compiler doesn't understand CUDA code, "
                    "the only possible memory space is 'MemorySpace::Host'!");
      return true;
#else
      cudaPointerAttributes attributes;
      const cudaError_t cuda_error = cudaPointerGetAttributes(&attributes, ptr);
      if (cuda_error != cudaErrorInvalidValue)
        {
          AssertCuda(cuda_error);
          if (std::is_same<MemorySpaceType, MemorySpace::Host>::value)
            return (attributes.type == cudaMemoryTypeHost) ||
                   (attributes.type == cudaMemoryTypeUnregistered);
          else
            return attributes.type == cudaMemoryTypeDevice;
        }
      else
        {
          // ignore and reset the error since host pointers produce an error
          cudaGetLastError();
          return std::is_same<MemorySpaceType, MemorySpace::Host>::value;
        }
#endif
    }
  } // namespace ArrayViewHelper
} // namespace internal



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView()
  : starting_element(nullptr)
  , n_elements(0)
{}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  value_type *      starting_element,
  const std::size_t n_elements)
  : starting_element(starting_element)
  , n_elements(n_elements)
{
  Assert(
    n_elements == 0 ||
      internal::ArrayViewHelper::is_in_correct_memory_space<MemorySpaceType>(
        starting_element),
    ExcMessage("The memory space indicated by the template parameter "
               "and the one derived from the pointer value do not match!"));
}



template <typename ElementType, typename MemorySpaceType>
inline void
ArrayView<ElementType, MemorySpaceType>::reinit(value_type *starting_element,
                                                const std::size_t n_elements)
{
  *this = ArrayView(starting_element, n_elements);
}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(ElementType &element)
  : starting_element(&element)
  , n_elements(1)
{}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  const ArrayView<typename std::remove_cv<value_type>::type, MemorySpaceType>
    &view)
  : starting_element(view.starting_element)
  , n_elements(view.n_elements)
{}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  const std::vector<typename std::remove_cv<value_type>::type> &vector)
  : // use delegating constructor
  ArrayView(vector.data(), vector.size())
{
  // the following static_assert is not strictly necessary because,
  // if we got a const std::vector reference argument but ElementType
  // is not itself const, then the call to the forwarding constructor
  // above will already have failed: vector.data() will have returned
  // a const pointer, but we need a non-const pointer.
  //
  // nevertheless, leave the static_assert in since it provides a
  // more descriptive error message that will simply come after the first
  // error produced above
  static_assert(std::is_const<value_type>::value == true,
                "This constructor may only be called if the ArrayView "
                "object has a const value_type. In other words, you can "
                "only create an ArrayView to const values from a const "
                "std::vector.");
}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  std::vector<typename std::remove_cv<value_type>::type> &vector)
  : // use delegating constructor
  ArrayView(vector.data(), vector.size())
{}



template <typename ElementType, typename MemorySpaceType>
template <std::size_t N>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  const std::array<typename std::remove_cv<value_type>::type, N> &vector)
  : // use delegating constructor
  ArrayView(vector.data(), vector.size())
{
  // the following static_assert is not strictly necessary because,
  // if we got a const std::array reference argument but ElementType
  // is not itself const, then the call to the forwarding constructor
  // above will already have failed: vector.data() will have returned
  // a const pointer, but we need a non-const pointer.
  //
  // nevertheless, leave the static_assert in since it provides a
  // more descriptive error message that will simply come after the first
  // error produced above
  static_assert(std::is_const<value_type>::value == true,
                "This constructor may only be called if the ArrayView "
                "object has a const value_type. In other words, you can "
                "only create an ArrayView to const values from a const "
                "std::array.");
}



template <typename ElementType, typename MemorySpaceType>
template <std::size_t N>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  ElementType (&array)[N])
  : ArrayView(&array[0], N)
{}



template <typename ElementType, typename MemorySpaceType>
template <std::size_t N>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  std::array<typename std::remove_cv<value_type>::type, N> &vector)
  : // use delegating constructor
  ArrayView(vector.data(), vector.size())
{}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::
operator==(const ArrayView<const value_type, MemorySpaceType> &other_view) const
{
  return (other_view.data() == starting_element) &&
         (other_view.size() == n_elements);
}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::
operator==(const ArrayView<typename std::remove_cv<value_type>::type,
                           MemorySpaceType> &other_view) const
{
  return (other_view.data() == starting_element) &&
         (other_view.size() == n_elements);
}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::
operator!=(const ArrayView<const value_type, MemorySpaceType> &other_view) const
{
  return !(*this == other_view);
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::value_type *
ArrayView<ElementType, MemorySpaceType>::data() const noexcept
{
  if (n_elements == 0)
    return nullptr;
  else
    return starting_element;
}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::
operator!=(const ArrayView<typename std::remove_cv<value_type>::type,
                           MemorySpaceType> &other_view) const
{
  return !(*this == other_view);
}



template <typename ElementType, typename MemorySpaceType>
inline std::size_t
ArrayView<ElementType, MemorySpaceType>::size() const
{
  return n_elements;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::iterator
ArrayView<ElementType, MemorySpaceType>::begin() const
{
  return starting_element;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::iterator
ArrayView<ElementType, MemorySpaceType>::end() const
{
  return starting_element + n_elements;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::const_iterator
ArrayView<ElementType, MemorySpaceType>::cbegin() const
{
  return starting_element;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::const_iterator
ArrayView<ElementType, MemorySpaceType>::cend() const
{
  return starting_element + n_elements;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::value_type &
  ArrayView<ElementType, MemorySpaceType>::operator[](const std::size_t i) const
{
  AssertIndexRange(i, n_elements);
  Assert(
    (std::is_same<MemorySpaceType, MemorySpace::Host>::value),
    ExcMessage(
      "Accessing elements is only allowed if the data is stored in CPU memory!"));

  return *(starting_element + i);
}



#ifndef DOXYGEN
namespace internal
{
  namespace ArrayViewHelper
  {
    /**
     * 返回在给定的迭代器范围内，通过取消迭代器得到的对象是否在内存中形成一个连续的范围。
     *
     */
    template <class Iterator>
    bool
    is_contiguous(const Iterator &first, const Iterator &last)
    {
      const auto n = std::distance(first, last);
      for (typename std::decay<decltype(n)>::type i = 0; i < n; ++i)
        if (std::addressof(*(std::next(first, i))) !=
            std::next(std::addressof(*first), i))
          return false;
      return true;
    }


    /**
     * 返回在给定的迭代器范围内通过解引用迭代器得到的对象是否在内存中形成一个连续的范围。
     * 这种对（  @p const  或非  @p const)
     * 指针的特殊化无条件地返回  @p true
     * ，因为指针所指向的对象是连续的这一事实已经嵌入到C++的内存模型中。
     *
     */
    template <class T>
    constexpr bool
    is_contiguous(T *, T *)
    {
      return true;
    }
  } // namespace ArrayViewHelper
} // namespace internal
#endif



/**
 * 创建一个ArrayView，它需要一对迭代器作为参数。ArrayView的类型是从迭代器的值类型推断出来的（例如，从两个常数迭代器创建的视图将有一个常数类型）。
 * @warning  迭代器 @p begin 和 @p end
 * 必须绑定（以通常的半开放方式）一个连续的内存范围内的值。这个函数的目的是用于像
 * <code>boost::container::small_vector</code> or <code>std::vector</code>
 * 这样的容器中的迭代器，而不会在例如
 * <code>boost::container::stable_vector</code> or <code>std::deque</code>
 * 中正确工作。在调试模式下，我们检查所提供的迭代器是否确实代表连续的内存。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename Iterator, typename MemorySpaceType = MemorySpace::Host>
ArrayView<typename std::remove_reference<
            typename std::iterator_traits<Iterator>::reference>::type,
          MemorySpaceType>
make_array_view(const Iterator begin, const Iterator end)
{
  static_assert(
    std::is_same<typename std::iterator_traits<Iterator>::iterator_category,
                 typename std::random_access_iterator_tag>::value,
    "The provided iterator should be a random access iterator.");
  Assert(begin <= end,
         ExcMessage(
           "The beginning of the array view should be before the end."));
  Assert(internal::ArrayViewHelper::is_contiguous(begin, end),
         ExcMessage("The provided range isn't contiguous in memory!"));
  // the reference type, not the value type, knows the constness of the iterator
  return ArrayView<typename std::remove_reference<
                     typename std::iterator_traits<Iterator>::reference>::type,
                   MemorySpaceType>(std::addressof(*begin), end - begin);
}



/**
 * 从一对指针创建一个视图。  <code>ElementType</code>
 * 可以是const-qualified的。
 * @warning  指针 @p begin 和 @p end
 * 必须绑定（以通常的半开放方式）一个连续的内存中的值范围。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType, typename MemorySpaceType = MemorySpace::Host>
ArrayView<ElementType, MemorySpaceType>
make_array_view(ElementType *const begin, ElementType *const end)
{
  Assert(begin <= end,
         ExcMessage(
           "The beginning of the array view should be before the end."));
  return ArrayView<ElementType, MemorySpaceType>(begin, end - begin);
}



/**
 * 从一个ArrayView本身创建一个视图。 这个函数用于 @p const
 * 对ArrayView类型对象的引用。它的存在只是为了兼容的目的。
 * @param[in]  array_view 我们希望复制的ArrayView。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename Number, typename MemorySpaceType>
inline ArrayView<const Number, MemorySpaceType>
make_array_view(const ArrayView<Number, MemorySpaceType> &array_view)
{
  return make_array_view(array_view.cbegin(), array_view.cend());
}



/**
 * 从ArrayView本身创建一个视图。 这个函数用于非  @p const
 * ArrayView类型的对象的引用。它的存在只是为了兼容的目的。
 * @param[in]  array_view 我们希望复制的ArrayView。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename Number, typename MemorySpaceType>
inline ArrayView<Number, MemorySpaceType>
make_array_view(ArrayView<Number, MemorySpaceType> &array_view)
{
  return make_array_view(array_view.begin(), array_view.end());
}



/**
 * 为整个Tensor对象创建一个视图。这相当于用第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于 @p const
 * 对Tensor类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对
 * @p const 对象集合的一个视图。
 * @param[in]  张量
 * 我们希望有一个数组视图对象的张量。数组视图对应于
 * <em> 整个 </em>
 * 对象，但条目在数组中的呈现顺序是一个实现细节，不应该被依赖。
 * @relatesalso ArrayView
 *
 *
 */
template <int rank, int dim, typename Number>
inline ArrayView<const Number>
make_array_view(const Tensor<rank, dim, Number> &tensor)
{
  return make_array_view(tensor.begin_raw(), tensor.end_raw());
}



/**
 * 为整个Tensor对象创建一个视图。这相当于用一个指向第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于非  @p const
 * Tensor类型对象的引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的视图。
 * @param[in]  Tensor
 * 我们希望有一个数组视图对象的Tensor。数组视图对应于
 * <em> 整个 </em>
 * 对象，但条目在数组中的呈现顺序是一个实现细节，不应该被依赖。
 * @relatesalso  ArrayView
 *
 *
 */
template <int rank, int dim, typename Number>
inline ArrayView<Number>
make_array_view(Tensor<rank, dim, Number> &tensor)
{
  return make_array_view(tensor.begin_raw(), tensor.end_raw());
}



/**
 * 为整个SymmetricTensor对象创建一个视图。这相当于用第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于 @p const
 * 对SymmetricTensor类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对
 * @p const 对象集合的一个视图。
 * @param[in]  tensor
 * 我们希望有一个数组视图对象的SymmetricTensor。数组视图对应于
 * <em> 整个 </em>
 * 对象，但条目在数组中的呈现顺序是一个实现细节，不应该被依赖。
 * @relatesalso ArrayView
 *
 *
 */
template <int rank, int dim, typename Number>
inline ArrayView<const Number>
make_array_view(const SymmetricTensor<rank, dim, Number> &tensor)
{
  return make_array_view(tensor.begin_raw(), tensor.end_raw());
}



/**
 * 为整个SymmetricTensor对象创建一个视图。这相当于用第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于对SymmetricTensor类型的对象的非 @p const
 * 引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的视图。
 * @param[in]  tensor
 * 我们希望有一个数组视图对象的SymmetricTensor。数组视图对应于
 * <em> 整个 </em>
 * 对象，但条目在数组中的呈现顺序是一个实现细节，不应该被依赖。
 * @relatesalso ArrayView
 *
 *
 */
template <int rank, int dim, typename Number>
inline ArrayView<Number>
make_array_view(SymmetricTensor<rank, dim, Number> &tensor)
{
  return make_array_view(tensor.begin_raw(), tensor.end_raw());
}



/**
 * 为整个C风格的数组创建一个视图。这相当于用一个指向第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 产生的ArrayView是否可写，取决于ElementType是否为常量类型。
 * @param[in]  array
 * 我们希望有一个ArrayView对象的C型数组。ArrayView对应的是
 * <em> 整个 </em> 向量。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType, int N>
inline ArrayView<ElementType> make_array_view(ElementType (&array)[N])
{
  return ArrayView<ElementType>(array, N);
}



/**
 * 为整个Vector对象创建一个视图。这相当于用一个指向第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于非 @p const
 * 对Vector类型对象的引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的视图。
 * @param[in]  vector
 * 我们希望有一个数组视图对象的Vector。数组视图对应于
 * <em>  整个 </em>  Vector。
 * @relatesalso ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(Vector<ElementType> &vector)
{
  return ArrayView<ElementType>(vector.begin(), vector.size());
}



/**
 * 创建一个视图到整个Vector对象。这相当于用一个指向第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于 @p const
 * 对Vector类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对
 * @p const 对象集合的一个视图。
 * @param[in]  vector
 * 我们希望有一个数组视图对象的Vector。该数组视图对应于
 * <em> 整个 </em> 向量。
 * @relatesalso ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Vector<ElementType> &vector)
{
  return ArrayView<const ElementType>(vector.begin(), vector.size());
}



/**
 * 创建一个视图到整个 std::vector
 * 对象。这相当于用一个指向第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于非  @p const
 * 矢量类型的对象的引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的视图。
 * @param[in]  矢量
 * 我们希望有一个阵列视图对象的矢量。该数组视图对应于
 * <em> 整个 </em> 向量。
 * @relatesalso ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(std::vector<ElementType> &vector)
{
  return ArrayView<ElementType>(vector.data(), vector.size());
}



/**
 * 创建一个视图到整个 std::vector
 * 对象。这相当于用第一个元素的指针和给定参数的大小初始化一个ArrayView对象。
 * 这个函数用于 @p const
 * 对矢量类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对一组
 * @p const 对象的视图。
 * @param[in]  矢量
 * 我们希望有一个数组视图对象的矢量。该数组视图对应于
 * <em> 整个 </em> 向量。
 * @relatesalso ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const std::vector<ElementType> &vector)
{
  return ArrayView<const ElementType>(vector.data(), vector.size());
}



/**
 * 创建一个视图到 std::vector
 * 对象的一部分。这相当于用一个指向 @p starting_index-
 * 第三元素的指针和 @p size_of_view
 * 作为视图的长度初始化ArrayView对象。
 * 这个函数用于对矢量类型对象的非 @p const
 * 引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的视图。
 * @param[in]  矢量 我们希望有一个数组视图对象的矢量。
 * @param[in]  starting_index
 * 将成为该视图一部分的向量第一个元素的索引。  @param[in]
 * size_of_view 新ArrayView中的元素数量。 @pre   <code>starting_index
 * + size_of_view <= vector.size()</code>
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(std::vector<ElementType> &vector,
                const std::size_t         starting_index,
                const std::size_t         size_of_view)
{
  Assert(starting_index + size_of_view <= vector.size(),
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of the given vector."));
  return ArrayView<ElementType>(&vector[starting_index], size_of_view);
}



/**
 * 创建一个视图到 std::vector
 * 对象的一部分。这相当于用一个指向 @p starting_index-
 * 第1个元素的指针和 @p size_of_view
 * 作为视图的长度初始化ArrayView对象。 这个函数用于 @p
 * const
 * 对矢量类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对一组
 * @p const 对象的视图。
 * @param[in]  矢量 我们希望有一个数组视图对象的矢量。
 * @param[in]  starting_index
 * 将成为该视图一部分的向量的第一个元素的索引。
 * @param[in]  size_of_view 新ArrayView中的元素数量。 @pre
 * <code>starting_index + size_of_view <= vector.size()</code>
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const std::vector<ElementType> &vector,
                const std::size_t               starting_index,
                const std::size_t               size_of_view)
{
  Assert(starting_index + size_of_view <= vector.size(),
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of the given vector."));
  return ArrayView<const ElementType>(&vector[starting_index], size_of_view);
}



/**
 * 为Table<2>对象的整个行创建一个视图。这相当于用一个指向给定行的第一个元素的指针初始化一个ArrayView对象，并将该行的长度作为视图的长度。
 * 这个函数用于非 @p const
 * 对Table类型对象的引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的视图。
 * @param[in]  table
 * 我们希望有一个数组视图对象的表。数组视图对应于  <em>
 * 整个  </em>  行。  @param[in]  row
 * 该视图所对应的表的行的索引。
 * @relatesalso ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<ElementType>
  make_array_view(Table<2, ElementType> &                         table,
                  const typename Table<2, ElementType>::size_type row)
{
  AssertIndexRange(row, table.size()[0]);
  return ArrayView<ElementType>(&table[row][0], table.size()[1]);
}



/**
 * 为整个Table<2>对象创建一个视图。这相当于用一个指向给定表的第一个元素的指针初始化一个ArrayView对象，并将表的条目数作为视图的长度。
 * 这个函数用于非 @p const
 * 的Table类型对象的引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的一个视图。
 * @param[in]  table
 * 我们希望有一个数组视图对象的表。数组视图对应于 <em>
 * 整个 </em>
 * 表，但条目在数组中的呈现顺序是一个实现细节，不应该被依赖。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<ElementType> make_array_view(Table<2, ElementType> &table)
{
  return ArrayView<ElementType>(&table[0][0], table.n_elements());
}



/**
 * 为整个Table<2>对象创建一个视图。这相当于用一个指向给定表的第一个元素的指针初始化一个ArrayView对象，并将表项的数量作为视图的长度。
 * 这个函数用于 @p const
 * 对Table类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对一组
 * @p const 对象的视图。
 * @param[in]  table
 * 我们希望有一个数组视图对象的表。数组视图对应于 <em>
 * 整个 </em>
 * 表，但条目在数组中的呈现顺序是一个实现细节，不应该被依赖。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Table<2, ElementType> &table)
{
  return ArrayView<const ElementType>(&table[0][0], table.n_elements());
}


/**
 * 为整个LAPACKFullMatrix对象创建一个视图。这相当于用一个指向给定对象的第一个元素的指针初始化一个ArrayView对象，并将条目数作为视图的长度。
 * 这个函数用于 @p non-const
 * 对LAPACKFullMatrix类型对象的引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组
 * @p non-const 对象的视图。
 * @param[in]  矩阵
 * 我们希望有一个阵列视图对象的LAPACKFullMatrix。数组视图对应于
 * <em> 整个 </em>
 * 对象，但条目在数组中的呈现顺序是一个实现细节，不应该被依赖。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(LAPACKFullMatrix<ElementType> &matrix)
{
  return ArrayView<ElementType>(&matrix(0, 0), matrix.n_elements());
}



/**
 * 为整个LAPACKFullMatrix对象创建一个视图。这相当于用一个指向给定对象的第一个元素的指针初始化一个ArrayView对象，并将条目数作为视图的长度。
 * 这个函数用于 @p const
 * 对LAPACKFullMatrix类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对
 * @p const 对象集合的一个视图。
 * @param[in]  矩阵
 * 我们希望有一个阵列视图对象的LAPACKFullMatrix。阵列视图对应于
 * <em> 整个 </em>
 * 对象，但条目在阵列中的呈现顺序是一个实现细节，不应依赖。
 * @relatesalso ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const LAPACKFullMatrix<ElementType> &matrix)
{
  return ArrayView<const ElementType>(&matrix(0, 0), matrix.n_elements());
}



/**
 * 为Table<2>对象的整个行创建一个视图。这相当于用一个指向给定行的第一个元素的指针初始化一个ArrayView对象，并将该行的长度作为视图的长度。
 * 这个函数用于 @p const
 * 对Table类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对一组
 * @p const 对象的视图。
 * @param[in]  table
 * 我们希望有一个数组视图对象的表。数组视图对应于一个
 * <em> 整个 </em> 行。  @param[in]  row
 * 该视图所对应的表的行的索引。
 * @relatesalso ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Table<2, ElementType> &                   table,
                const typename Table<2, ElementType>::size_type row)
{
  AssertIndexRange(row, table.size()[0]);
  return ArrayView<const ElementType>(&table[row][0], table.size()[1]);
}



/**
 * 为Table<2>对象的一行（部分）创建一个视图。
 * 这个函数用于对Table类型对象的非  @p const
 * 引用。这种对象包含可以被写入的元素。因此，这个函数的返回类型是对一组可写对象的视图。
 * @param[in]  table
 * 我们希望有一个数组视图对象的表。数组视图对应于  <em>
 * 整个  </em>  行。  @param[in]  row
 * 该视图所对应的表的行的索引。  @param[in]  starting_column
 * 与该视图的第一个元素相对应的表内的列的索引。
 * @param[in]  size_of_view
 * 这个视图应该有多少个元素。这对应于该视图应该对应的当前行中的列数。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<ElementType> make_array_view(
  Table<2, ElementType> &                         table,
  const typename Table<2, ElementType>::size_type row,
  const typename Table<2, ElementType>::size_type starting_column,
  const std::size_t                               size_of_view)
{
  AssertIndexRange(row, table.size()[0]);
  AssertIndexRange(starting_column, table.size()[1]);
  Assert(starting_column + size_of_view <= table.size()[1],
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of a column of the given table."));
  return ArrayView<ElementType>(&table[row][starting_column], size_of_view);
}



/**
 * 为Table<2>对象的某一行创建一个视图（一部分）。
 * 这个函数用于 @p const
 * 对Table类型对象的引用，因为它们包含不可变的元素。因此，这个函数的返回类型是对一组
 * @p const 对象的视图。
 * @param[in]  table
 * 我们希望有一个数组视图对象的表。数组视图对应于一个
 * <em> 整个 </em> 行。  @param[in]  row
 * 该视图所对应的表的行的索引。  @param[in]  starting_column
 * 与该视图的第一个元素相对应的表内的列的索引。
 * @param[in]  size_of_view
 * 这个视图应该有多少个元素。这对应于该视图应该对应的当前行中的列数。
 * @relatesalso  ArrayView
 *
 *
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Table<2, ElementType> &                   table,
                const typename Table<2, ElementType>::size_type row,
                const typename Table<2, ElementType>::size_type starting_column,
                const std::size_t                               size_of_view)
{
  AssertIndexRange(row, table.size()[0]);
  AssertIndexRange(starting_column, table.size()[1]);
  Assert(starting_column + size_of_view <= table.size()[1],
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of a column of the given table."));
  return ArrayView<const ElementType>(&table[row][starting_column],
                                      size_of_view);
}



/* 创建一个不允许修改其指向的容器的视图。如果传入的对象已经不是`const'，并且一个函数在签名中要求视图为常量内存，这就很有用。
* 这个函数返回一个`ArrayView<const T>`类型的对象，其中`T`是容器的元素类型。
*  @relatesalso  ArrayView

* 
*
*/
template <typename Container>
inline auto
make_const_array_view(const Container &container)
  -> decltype(make_array_view(container))
{
  return make_array_view(container);
}


DEAL_II_NAMESPACE_CLOSE

#endif


