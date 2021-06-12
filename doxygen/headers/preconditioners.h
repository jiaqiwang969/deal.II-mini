//include/deal.II-translator/A-headers/preconditioners_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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
 *    @defgroup Preconditioners Preconditioners and Relaxation Operators
 * <h3>Preconditioners</h3>
 * 先决条件是用来加速线性系统的迭代求解的。典型的前置条件是Jacobi、Gauss-Seidel或SSOR，但该库也支持更复杂的前置条件，如Vanka或不完全LU分解（ILU）。此外，稀疏直接求解器在可用的情况下也可以作为预处理器使用。
 * 广义上讲，预处理器是一种运算器，它与矩阵相乘以改善条件。其想法是，经过预处理的系统<i>P<sup>-1</sup>Ax
 * = P<sup>-1</sup>b</i>比原始系统<i>Ax =
 * b</i>更容易解决。这到底意味着什么，取决于矩阵的结构，在此不能作一般性讨论。对于对称的正定矩阵<i>A</i>和<i>P</i>，这意味着<i>P<sup>-1</sup>A</i>的光谱条件数（最大和最小特征值的商）比<i>A</i>的要小得多。
 * 在最简单的例子中，Richardson迭代，在SolverRichardson中实现，预处理迭代看起来像@f[
 * x^{k+1} = x^k
 *
 * - P^{-1} \bigl(A x^k
 *
 * - b\bigr).
 * @f]。 因此，预处理相当于对残差应用一个线性算子，因此，预处理器<i>P<sup>-1</sup></i>的动作被实现为<tt>vmult()</tt>。deal.II中需要预处理器的模板用 @ref ConceptPreconditionerType "预处理器类型概念 "来表示。在实践中，我们通常可以将任何定义了 <code>vmult()</code> and <code>Tvmult()</code> 的矩阵类对象作为一个预处理程序。本模块中的所有预处理类都实现了这个接口。
 * 当用于Krylov空间方法时，由该方法决定是用<i>A</i>替换<i>P<sup>-1</sup>A</i>的乘法（例如SolverBicgstab），还是做更复杂的事情。例如，SolverCG使用<i>P<sup>-1</sup></i>来定义内积，这就是为什么它需要一个对称的正定算子<i>P</i>的原因。
 * <h3>Relaxation methods</h3> 许多预处理程序依赖于加法拆分<i>A
 * = P
 *
 * - N</i>为两个矩阵。在这种情况下，上述Richardson方法的迭代步骤可以简化为@f[
 * x^{k+1} = P^{-1} \bigl(N x^k + b\bigr),
 * @f]，从而完全避免了与<i>A</i>的乘法。我们把以这种方式将前一个迭代<i>x<sup>k</sup></i>映射到下一个迭代的算子称为放松算子。它们的通用接口由 @ref ConceptRelaxationType "松弛类型概念 "给出。本模块中名字以<tt>Relaxation</tt>开头的类实现了这个接口，还有预处理器PreconditionJacobi、PreconditionSOR、PreconditionBlockJacobi、PreconditionBlockSOR和PreconditionBlockSSOR。
 * <h3>The interface</h3>
 * 在这一节中，我们讨论了预设条件器通常必须提供的接口，以便在deal.II库内工作。
 * <h4>Initialization</h4>
 * 为了能够存储在容器中，所有的预处理程序都有一个没有参数的构造函数。由于这通常会产生一个无用的对象，所有的预处理程序都有一个函数
 * @code
 * void initialize (...)
 * @endcode
 *
 * 这个函数接收要预处理的矩阵以及其他所需的参数，并设置预处理程序的内部结构。
 * <h4>Relaxation methods</h4>
 *
 * @ingroup LAC
 * @ingroup Matrices
 *
 */


