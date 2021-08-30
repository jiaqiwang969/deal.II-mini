//include/deal.II-translator/A-headers/instantiations_0.txt
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


/**
 * @page Instantiations Template instantiations
 * 复杂的类和函数模板的实例化在编译时间和磁盘空间方面都很昂贵。因此，我们尽可能地将模板的声明和实现分开，并确保只有在必要时才由编译器读取其实现。
 * <tt>deal.II</tt>中的模板类可以分为三类，取决于可能的不同实例的数量。这三类将在下文中讨论。
 *
 *
 * @section  Inst1 已知和固定数量的实例化
 * 这些是具有模板参数的类，其值非常容易预测。典型的原型是
 * @code
 * template <int dim> class Function;
 * @endcode
 *
 * 在这里，我们在设计库的时候有少量的实例 (  <code>dim =
 * 1,2,3</code>  )
 * 已知。因此，这个类的成员函数被定义在源文件目录下的<tt>.cc</tt>文件中，我们在源文件中明确地将这些已知值实例化为模板。
 * 从应用程序的角度来看，你实际看到的只是模板的声明。成员函数的实际实例化发生在库中，并且是在你编译库的时候完成的，而不是在你编译应用程序代码的时候。
 * 对于这些类，为新的参数增加实例化需要改变库。然而，这很少需要，当然，除非你不满足于只在1d、2d或3d中计算。
 *
 *
 * @subsection  Inst1a 可用实例
 * 如果模板参数是<tt>dim</tt>，如果没有其他信息，可用的实例是为<tt>dim=1,2,3</tt>。
 * 还有一些类的情况（不取决于空间维度），只支持一定的、少量的模板参数，并在库中提供显式实例。特别是，这包括所有在标量底层存储值类型上进行模板化的线性代数类：我们只支持
 * <code>double</code>, <code>float</code>  ，以及在某些情况下支持
 * <code>std::complex@<double@></code>  和
 * <code>std::complex@<float@></code>  。
 *
 *
 * @section  Inst2 一些实例化，其中大部分是已知的
 * 这些是通常有少量实例的类模板，但可能需要额外的实例。因此，在库中提供了一组最可能的参数的实例化，但在一个特殊的头文件中提供了模板的实现，以便在有人想对一个未预见的参数进行实例化时，可以访问它。
 * 这方面的典型例子是一些线性代数类，它们将向量类型作为模板参数。例如，它们将在库中被实例化为
 * <code>Vector&lt;double&gt;</code> ,  <code>Vector&lt;float&gt;</code>,
 * <code>BlockVector&lt;double&gt;</code> , 和
 * <code>BlockVector&lt;float&gt;</code>
 * 。然而，它们也可以与其他矢量类型一起使用，如
 * <code>long double</code> 和 <code>std::complex@<long double@></code>
 * ，只要它们满足某些接口，包括不属于库的但可能在应用程序中定义的矢量类型。在这种情况下，应用程序可以通过手工实例化这些模板，如下一节所述。
 *
 *
 * @subsection  Inst2c 创建新实例
 * 从你的源文件中选择一个来提供所需的实例。假设你想让头文件<tt>XXXX</tt>中定义的类模板<tt>XX.h</tt>，用模板参数<tt>long
 * double</tt>来实例化。那么，你的文件应该包含以下几行
 * @code
 *                 // Include class template declaration
 * #include <xxxx.h>
 *                 // Include member definitions
 * #include <xxxx.templates.h>
 *
 * ...
 *
 * template class XXXX<long double>;
 * @endcode
 *
 *
 * @subsection  Inst2p提供的实例
 * 与 @ref
 * Inst1
 * 节中的类一样，库中提供的实例通常以类似于这样的形式列在该类的文档中。
 * @verbatim
 * Template Instantiations: some  (<p1>a,b,c<p2>)
 * @endverbatim
 *
 *
 * @section  Inst3 许多未知的实例
 * 这些类，不存在合理的预先确定的实例集。因此，所有的成员定义都包含在头文件中，并在需要的地方被实例化。
 * 一个例子是SmartPointer类模板，它几乎可以与任何模板参数一起使用。
 *
 */


