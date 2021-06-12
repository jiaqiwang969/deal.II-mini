//include/deal.II-translator/A-headers/laoperators_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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
 *    @defgroup LAOperators Linear Operators <h3>Linear Operator</h3>
 * deal.II包括支持以一种非常普遍的方式来描述线性变换。这是用LinearOperator类来完成的，就像 @ref ConceptMatrixType "MatrixType概念 "
 * 一样，它为<i>applying</i>向量上的线性操作定义了一个最小的接口。
 *
 * @code
 * std::function<void(Range &, const Domain &)> vmult;
 * std::function<void(Range &, const Domain &)> vmult_add;
 * std::function<void(Domain &, const Range &)> Tvmult;
 * std::function<void(Domain &, const Range &)> Tvmult_add;
 * @endcode
 *
 * LinearOperator类的最大优势在于它为复杂的矩阵-向量操作提供了语法糖。作为一个例子，考虑操作
 * $(A+k\,B)\,C$  ，其中  $A$  、  $B$  和  $C$
 * 表示（可能不同的）稀疏矩阵对象。为了构造一个LinearOperator
 * <code>op</code>
 * ，当应用于一个向量时执行上述计算，我们可以写道。
 * @code
 * #include <deal.II/lac/linear_operator_tools.h>
 *
 * dealii::SparseMatrix<double> A, B, C;
 * double k;
 * // Setup and assembly...
 *
 * const auto op_a = linear_operator(A);
 * const auto op_b = linear_operator(B);
 * const auto op_c = linear_operator(C);
 *
 * const auto op = (op_a + k op_b) op_c;
 * @endcode
 * 现在， <code>op</code>
 * 可以作为一个矩阵对象用于进一步的计算。
 * linear_operator()函数可以用来将一个普通的矩阵或预处理对象包装成一个LinearOperator。线性算子可以用transpose_operator()进行转置，或者通过使用inverse_operator()与迭代求解器一起进行反转。
 * 对于LinearOperator类型的对象，所有的向量空间操作，即加减法、标量乘法和组合（兼容的线性算子）都被实现。
 * @code
 * dealii::LinearOperator<> op_a, op_b;
 * double k;
 *
 * // vector space addition, subtraction and scalar multiplication
 * op_a + op_b;
 * op_a
 *
 * - op_b;
 * k op_a;
 * op_a k;
 *
 * // in-place variants
 * op_a += op_b;
 * op_a
 *
 * -= op_b;
 * op_a= k;
 *
 * // operator composition
 * op_a op_b;
 * op_a= op_b; // If op_b is an endomorphism of the domain space of op_a
 * @endcode
 *
 * block_operator()和block_diagonal_operator()提供了对单个线性算子的进一步封装，使其成为封锁的线性算子变体。
 * step-20  教程中有一个关于LinearOperator类的详细使用例子。
 *
 *
 * @note  如下所述，当使用LinearOperator作为 <code>res =
 * op_a*x</code>
 * 时，会在幕后生成PackagedOperation类实例。因此，用户程序必须包括这两个类的头文件才能编译成功。为了更容易决定在什么情况下包含哪些头文件，并防止隐藏的与模板相关的编译器错误，所有与LinearOperator相关的头文件都被归入了`<deal.II/lac/linear_operator_tools.h>`头文件。
 * <h3>Packaged Operation</h3> 通过 <code>operator*</code>
 * 将LinearOperator对象应用于一个向量，会产生一个PackagedOperation对象来存储这个计算。
 * PackagedOperation类允许对涉及向量和线性运算符的表达式进行懒惰的评估。这是通过存储计算表达式来实现的，只有当对象被隐含地转换为向量对象，或者
 * PackagedOperation::apply() （或 PackagedOperation::apply_add())
 * 被手动调用时才执行计算。这就避免了不必要的中间结果的临时存储。
 * 作为一个例子，考虑多个向量的相加。
 * @code
 * #include <deal.II/lac/linear_operator_tools.h>
 *
 * dealii::Vector<double> a, b, c, d;
 * // ..
 * dealii::Vector<double> result = a + b
 *
 * - c + d;
 * @endcode
 * 转换PackagedOperation <code>a + b
 *
 * - c + d</code>到一个向量的结果是，代码相当于以下代码
 * @code
 * dealii::Vector<double> a, b, c, d;
 * // ..
 * dealii::Vector<double> result = a;
 * result += b;
 * result
 *
 * -= c;
 * result += d;
 * @endcode
 * 这避免了任何中间存储。作为第二个例子（涉及一个LinearOperator对象），考虑计算一个残差
 * $b-Ax$  。
 *
 * @code
 * dealii::SparseMatrix<double> A;
 * dealii::Vector<double> b, x;
 * // ..
 * const auto op_a = linear_operator(A);
 *
 * dealii::Vector<double> residual =  b
 *
 * - op_a x;
 * @endcode
 * 这里，表达式<code>b
 *
 * - op_a x</code>的结果又是一个PackagedOperation类型的对象，它存储了应该使用两个向量和线性运算符执行的<i>sequence of
 * operations</i>。将表达式转换为矢量（就像这里发生的对矢量
 * <code>residual</code>
 * 的赋值一样），执行计算（见下面的注释）。
 *
 *
 * @note
 * 计算表达式的懒惰评估必然涉及对底层向量和矩阵对象的引用。例如，创建一个
 * <code>residual_expr</code> 对象
 * @code
 * auto residual_expr =  b
 *
 * - op_a x;
 * @endcode
 * 存储残差的计算表达式，引用向量  <code>b</code> and matrix
 * <code>A</code>
 * 。在这一点上，它不进行任何计算。特别是，如果
 * <code>b</code> 或 <code>A</code> 被改变<b>after</b>，
 * <code>residual_expr</code>
 * 的创建，每一个后续的表达式的评估都用新的值进行。
 * @code
 * auto residual_expr =  b
 *
 * - op_a x;
 * residual_expr.apply(tmp);  // tmp is a Vector<double>
 *
 * // modify b, or A
 *
 * residual_expr.apply(tmp2); // tmp2 is a Vector<double>
 *
 * // tmp and tmp2 are different
 * @endcode
 * 因此，作为一种保障，如果你想马上计算一个表达式的结果，总是明确地在左边使用一个向量类型（而不是
 * <code>auto</code>  ）。
 * @code
 * Vector<double> residual =  b
 *
 * - op_a x; // computes the residual at this point
 * @endcode
 *
 * @note   step-20
 * 教程中有一个PackagedOperation类的详细使用例子。
 *
 *
 * @note
 * LinearOperator的许多用例导致中间表达式需要一个PackagedOperation。为了一次性包含所有必要的头文件，可以考虑使用
 * @code
 * #include <deal.II/lac/linear_operator_tools.h>
 * @endcode
 *
 * @ingroup LAC
 * @ingroup MATRICES
 *
 */


