//include/deal.II-translator/A-headers/c++_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


/**
 *    @defgroup CPP11 deal.II and the C++11 standard
 * 从9.0版本开始，deal.II需要一个至少支持<a
 * href="http://en.wikipedia.org/wiki/C%2B%2B11">C++11</a>的编译器。作为其中的一部分，deal.II内部实现中的许多地方现在都在使用C++11中才引入的功能。
 * 也就是说，deal.II也有一些函数和类，使其使用C++11的功能更容易。
 * 一个例子是对C++11<a
 * href="http://en.wikipedia.org/wiki/C++11#Range-based_for_loop">range-based
 * for
 * loops</a>的支持。基于deal.II的代码通常有许多类似的循环。
 * @code
 * Triangulation<dim> triangulation;
 * ...
 * typename Triangulation<dim>::active_cell_iterator
 *   cell = triangulation.begin_active(),
 *   endc = triangulation.end();
 * for (; cell!=endc; ++cell)
 *   cell->set_refine_flag();
 * @endcode
 * 使用C++11的基于范围的for循环，你现在可以这样写。
 * @code
 * Triangulation<dim> triangulation;
 * ...
 * for (auto &cell : triangulation.active_cell_iterators())
 *   cell->set_refine_flag();
 * @endcode
 * 这依赖于诸如 Triangulation::active_cell_iterators(),
 * 和DoF处理类中的等价函数， DoFHandler::active_cell_iterators(),
 * hp::DoFHandler::active_cell_iterators().
 * 这些函数的变体为所有单元格（不仅仅是活动单元格）和单个层次上的单元格提供迭代器范围。
 * 库中还有许多其他的函数，允许习惯性的使用基于范围的for循环。例子有
 * GeometryInfo::face_indices(),   GeometryInfo::vertex_indices(),
 * FEValuesBase::quadrature_point_indices(), 等许多。
 * C++11还引入了[constexpr](https://en.cppreference.com/w/cpp/language/constexpr)变量和函数的概念。定义为
 * "constexpr
 * "的变量是在程序编译过程中计算出来的常量值，因此其初始化的运行时间成本为零。此外，`constexpr`常量有正确定义的生命期，完全避免了所谓的
 * "静态初始化顺序惨败"。%函数可以被标记为
 * "constexpr"，表明如果它们的输入参数是常量表达式，它们可以产生编译时常量返回值。此外，至少有一个`constexpr`构造函数的类可以被初始化为`constexpr`。
 * 作为一个例子，由于构造函数 Tensor::Tensor(const  array_type
 * &)是`constexpr`，我们可以在编译时用一个数组初始化一个张量，如：。
 * @code
 * constexpr double[2][2] entries = {{1., 0.}, {0., 1.}};
 * constexpr Tensor<2, 2> A(entries);
 * @endcode
 * 在这里，A的内容没有被存储在堆栈中。相反，它们在编译时被初始化并插入到可执行程序的`.data`部分。程序可以在运行时使用这些值，而无需花费时间进行初始化。初始化张量可以用一句话来简化。
 * @code
 * constexpr Tensor<2, 2> A({{1., 0.}, {0., 1.}});
 * @endcode
 * 一些函数如determinant()被指定为`constexpr`，但它们需要具有C++14能力的编译器。因此，这个函数在内部被声明为。
 * @code
 * template <int dim, typename Number>
 * DEAL_II_CONSTEXPR Number determinant(const Tensor<2, dim, Number> &t);
 * @endcode
 * 如果有一个支持C++14的编译器，宏 @ref
 * DEAL_II_CONSTEXPR
 * 简化为`constexpr`。否则，对于旧的编译器，它完全忽略了DEAL_II_CONSTEXPR。因此，在较新的编译器中，用户可以编写
 * @code
 * constexpr double det_A = determinant(A);
 * @endcode
 * 假设`A`是用`constexpr`指定器声明的。这个例子显示了使用`constexpr`的性能收益，因为这里我们在编译时进行了一个具有
 * $O(\text{dim}^3)$ 复杂性的操作，避免了任何运行时的成本。
 *
 */



/**
 * deal.II目前只需要一个符合C++11的编译器，但有一些来自C++14标准的函数和类，在编译器只支持C++11的情况下也很容易提供。
 * 这些都收集在当前的命名空间中。 最明显的例子是<a
 * href="https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique">`std::make_unique`</a>函数，它可以说是一个疏忽，因为它没有被包含在C++11中（鉴于C++11中有<a
 * href="https://en.cppreference.com/w/cpp/memory/shared_ptr/make_shared">`std::make_shared`</a>）。
 * 在这个命名空间中还有其他一些小的补充，使我们在这一点上已经可以使用C++14的特性，尽管我们不需要一个符合C++14的编译器。
 *
 *
 * @note
 * 如果使用的编译器确实支持C++14，那么这个命名空间的内容只是从命名空间`std`导入的类和函数。也就是说，我们回到了编译器提供的内容，而不是我们自己的实现。
 *
 */
namespace std_cxx14
{}



/**
 * deal.II目前只需要一个符合C++11的编译器，但有一些来自C++17标准的函数和类，在编译器只支持C++11的情况下也很容易提供。
 * 这些都收集在当前的命名空间。 最显著的例子是<a
 * href="https://en.cppreference.com/w/cpp/utility/optional">`std::optional`</a>类，它从C++17标准开始被引入到C++中。
 * 在这个命名空间中还有其他一些小的补充，使我们在这一点上已经可以使用C++17的特性，尽管我们不需要一个符合C++17的编译器。
 *
 *
 * @note
 * 如果使用的编译器确实支持C++17，那么这个命名空间的内容只是从命名空间`std`导入的类和函数。也就是说，我们回到了编译器提供的东西，而不是我们自己的实现。
 *
 */
namespace std_cxx17
{}



/**
 * deal.II目前只需要一个符合C++11的编译器，但有一些来自C++20标准的函数和类，在编译器只支持C++11的情况下也很容易提供。
 * 这些都收集在当前的命名空间中。 一个例子是<a
 * href="https://en.cppreference.com/w/cpp/ranges/iota_view"
 * >`std::ranges::iota_view`</a>类，它从C++20标准开始被引入到C++中。它被用作
 * GeometryInfo::face_indices(),  GeometryInfo::vertex_indices(), 和
 * FEValuesBase::quadrature_point_indices()
 * 等函数的返回类型，以支持基于范围的for循环（参见 @ref
 * CPP11
 * 中基于范围的for循环的例子，以及上述函数的文档）。
 * 在这个命名空间中还有其他一些小的补充，使我们在这一点上已经可以使用C++20的特性，尽管我们不需要一个兼容C++20的编译器。
 *
 *
 * @note
 * 如果使用的编译器确实支持C++20，那么这个命名空间的内容只是从命名空间`std`导入的类和函数。也就是说，我们回到了编译器提供的内容，而不是我们自己的实现。
 *
 */
namespace std_cxx20
{}


