//include/deal.II-translator/lac/schur_complement_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

#ifndef dealii_schur_complement_h
#define dealii_schur_complement_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/vector_memory.h>


DEAL_II_NAMESPACE_OPEN

/**
 * @name  创建一个与舒尔补码有关的LinearOperator
 *
 *
 */
//@{

/**
 * @relatesalso  LinearOperator
 * 返回一个LinearOperator，执行与Schur补数相关的操作。有两个额外的辅助函数，condense_schur_rhs()和postprocess_schur_solution()，为了用这个操作符在线性代数中执行任何有用的任务，很可能需要使用这些函数。
 * 我们以如下方式构建Schur补码的定义。
 * 考虑一个一般的线性方程组，它可以被分解成两个主要的方程组。@f{eqnarray*}{
 * \mathbf{K}\mathbf{d} = \mathbf{f}
 * \quad \Rightarrow\quad
 * \left(\begin{array}{cc}
 *  A & B \\ C & D
 * \end{array}\right)
 * \left(\begin{array}{cc}
 *  x \\ y
 * \end{array}\right)
 * =
 * \left(\begin{array}{cc}
 *  f \\ g
 * \end{array}\right),
 * @f}
 * 其中 $ A,B,C,D $ 代表矩阵 $ \mathbf{K} $
 * 的一般子块，同样地， $ \mathbf{d},\mathbf{f} $
 * 的一般子向量由 $ x,y,f,g $ 给出。
 * 这等同于以下两个声明。@f{eqnarray*}{
 * (1) \quad Ax + By &=& f \\
 * (2) \quad Cx + Dy &=& g \quad .
 * @f}
 *
 * 假设 $ A,D $ 既是正方形又是可倒置的，那么我们可以进行两种可能的替换，@f{eqnarray*}{
 * (3) \quad x &=& A^{-1}(f
 *
 * - By) \quad \text{from} \quad (1) \\
 * (4) \quad y &=& D^{-1}(g
 *
 * - Cx) \quad \text{from} \quad (2) ,
 * @f} 。
 * 这相当于对这个方程组进行块状高斯消除。
 * 就目前的实现而言，我们选择将（3）代入（2）@f{eqnarray*}{
 * C \: A^{-1}(f
 *
 * - By) + Dy &=& g \\
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -C \: A^{-1} \: By + Dy &=& g
 *
 * - C \: A^{-1} \: f \quad .
 * @f}。
 * 这导致了结果@f[
 * (5) \quad (D
 *
 * - C\: A^{-1} \:B)y  = g
 *
 * - C \: A^{-1} f
 *     \quad \Rightarrow \quad Sy = g'
 * @f]，其中 $ S = (D
 *
 * - C\: A^{-1} \:B) $ 是舒尔补数，而修改后的右手边矢量 $ g' =
 * g
 *
 * - C \: A^{-1} f $ 则来自凝结步骤。请注意，对于 $ S $
 * 的这种选择，子矩阵 $ D $
 * 不需要是可逆的，因此可能是空矩阵。理想情况下 $ A $
 * 应该是条件良好的。
 * 所以对于任何任意矢量 $ a $ ，舒尔补码进行以下操作。@f[
 * (6) \quad Sa = (D
 *
 * - C \: A^{-1} \: B)a
 * @f]
 * 一个典型的解决线性系统(1),(2)所需的步骤集是。1. 定义逆矩阵  @p A_inv  （使用inverse_operator()）。2. 2. 定义舒尔补码 $ S $ （使用schur_complement()）。3. 定义迭代逆矩阵 $ S^{-1} $ ，使(6)成立。有必要使用带有预处理程序的求解器来计算 $ S $ 的近似逆运算，因为我们从未直接计算 $ S $ ，而是计算其运算的结果。为了实现这一点，可以再次使用inverse_operator()与我们刚刚构建的Schur补码结合起来。请注意， $ S $ 和它的预处理程序都是在与 $ D $ 相同的空间内运行的。4. 4. 使用 condense_schur_rhs() 对(5)的RHS进行预处理步骤。   @f[
 *    g' = g
 *
 * - C \: A^{-1} \: f
 *  @f] 5. 求解（5）中的 $ y $ 。   @f[
 * y =  S^{-1} g'
 *  @f] 6. 使用postprocess_schur_solution()执行（3）中的后处理步骤。   @f[
 * x =  A^{-1} (f
 *
 * - By) @f]
 * 下面给出了该算子在全耦合系统中的典型用法说明。
 *
 * @code
 * #include<deal.II/lac/schur_complement.h>
 *
 * // Given BlockMatrix K and BlockVectors d,F
 *
 * // Decomposition of tangent matrix
 * const auto A = linear_operator(K.block(0,0));
 * const auto B = linear_operator(K.block(0,1));
 * const auto C = linear_operator(K.block(1,0));
 * const auto D = linear_operator(K.block(1,1));
 *
 * // Decomposition of solution vector
 * auto x = d.block(0);
 * auto y = d.block(1);
 *
 * // Decomposition of RHS vector
 * auto f = F.block(0);
 * auto g = F.block(1);
 *
 * // Construction of inverse of Schur complement
 * const auto prec_A = PreconditionSelector<...>(A);
 * const auto A_inv = inverse_operator<...>(A,prec_A);
 * const auto S = schur_complement(A_inv,B,C,D);
 *
 * // D and S operate on same space
 * const auto S_prec = PreconditionSelector<...>(D);
 * const auto S_inv = inverse_operator<...>(S,...,prec_S);
 *
 * // Solve reduced block system
 * // PackagedOperation that represents the condensed form of g
 * auto rhs = condense_schur_rhs (A_inv,C,f,g);
 *
 * // Solve for y
 * y = S_inv rhs;
 *
 * // Compute x using resolved solution y
 * x = postprocess_schur_solution (A_inv,B,y,f);
 * @endcode
 *
 * 在上面的例子中， $ S $ 的预处理器被定义为 $ D $
 * 的预处理器，这是有效的，因为它们在同一空间内运行。然而，如果
 * $ D $ 和 $ S $
 * 太不相似，那么这可能导致大量的求解器迭代，因为 $
 * \text{prec}(D) $ 不是 $ S^{-1} $ 的良好近似。
 * 在这种情况下，更好的预处理程序是为 $ S^{-1} $
 * 提供更具代表性的近似。一种方法显示在 step-22 中，其中
 * $ D $ 是空矩阵， $ S^{-1} $
 * 的预处理器是由这个空间的质量矩阵导出的。
 * 从另一个角度来看，类似的结果可以通过首先构造一个代表
 * $ S $ 的近似对象来实现，其中昂贵的操作，即 $ A^{-1} $
 * ，被近似。此后，我们构建近似的反算子 $ \tilde{S}^{-1} $
 * ，然后将其作为计算 $ S^{-1} $ 的前提条件。
 *
 * @code
 * // Construction of approximate inverse of Schur complement
 * const auto A_inv_approx = linear_operator(preconditioner_A);
 * const auto S_approx = schur_complement(A_inv_approx,B,C,D);
 *
 * // D and S_approx operate on same space
 * const auto S_approx_prec = PreconditionSelector<...>(D);
 *
 * // Inner solver: Typically limited to few iterations
 * //               using IterationNumberControl
 * auto S_inv_approx = inverse_operator(S_approx,...,S_approx_prec);
 *
 * // Construction of exact inverse of Schur complement
 * const auto S = schur_complement(A_inv,B,C,D);
 *
 * // Outer solver
 * const auto S_inv = inverse_operator(S,...,S_inv_approx);
 *
 * // Solve reduced block system
 * auto rhs = condense_schur_rhs (A_inv,C,f,g);
 *
 * // Solve for y
 * y = S_inv rhs;
 * x = postprocess_schur_solution (A_inv,B,y,f);
 * @endcode
 * 请注意，由于构建 @c  S_inv_approx和随后的 @c
 * S_inv，有一对嵌套的迭代求解器，可能会共同消耗大量的资源。因此，在构建迭代逆运算器的过程中，应该注意选择。我们可以考虑使用IterationNumberControl（或类似的机制）来限制内部求解器的迭代次数。这可以控制近似逆运算
 * $ \tilde{S}^{-1} $ 的精度，该运算仅作为 $ S^{-1} $
 * 的预处理程序发挥作用。此外， $ \tilde{S}^{-1} $
 * 的先决条件，在本例中是 $ \text{prec}(D) $
 * ，理想情况下应该是计算成本低的。
 * 然而，如果使用基于IterationNumberControl的迭代求解器作为预处理器，那么预处理操作就不是线性操作了。这里最好采用SolverFGMRES（灵活的GMRES）这样的灵活求解器作为外部求解器，以处理预处理器的可变行为。否则，迭代求解器可能会在预处理器的容忍度附近停滞不前，或者普遍表现得不稳定。另外，使用ReductionControl可以确保前置条件器总是以相同的容差进行求解，从而使其行为恒定。
 * 这种功能的更多例子可以在测试套件中找到，如
 * <code>tests/lac/schur_complement_01.cc</code>
 * 。使用schur_complement解决一个多组件问题（即  step-22
 * ），可以在  <code>tests/lac/schur_complement_03.cc</code>
 * 中找到。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range_1,
          typename Domain_1,
          typename Range_2,
          typename Domain_2,
          typename Payload>
LinearOperator<Range_2, Domain_2, Payload>
schur_complement(const LinearOperator<Domain_1, Range_1, Payload> &A_inv,
                 const LinearOperator<Range_1, Domain_2, Payload> &B,
                 const LinearOperator<Range_2, Domain_1, Payload> &C,
                 const LinearOperator<Range_2, Domain_2, Payload> &D)
{
  // We return the result of the compound LinearOperator
  // directly, so as to ensure that the underlying Payload
  // definition aligns with the operations expressed here.
  // All of the memory allocations etc. are taken care of
  // internally.
  if (D.is_null_operator == false)
    return D - C * A_inv * B;
  else
    return -1.0 * C * A_inv * B;
}

//@}


/**
 * @name  创建与舒尔补数有关的打包操作对象
 *
 */
//@{

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 对于方程组@f{eqnarray*}{
 * Ax + By &=& f \\
 * Cx + Dy &=& g \quad ,
 * @f}来说
 * 该操作对RHS子向量 @p g 进行预处理（凝结）步骤，以便舒尔补码可以用来解决这个方程组。更具体地说，它产生一个对象，代表子向量 @p g, 的浓缩形式，即@f[
 * g' = g
 *
 * - C \: A^{-1} \: f
 * @f] 。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range_1,
          typename Domain_1,
          typename Range_2,
          typename Payload>
PackagedOperation<Range_2>
condense_schur_rhs(const LinearOperator<Range_1, Domain_1, Payload> &A_inv,
                   const LinearOperator<Range_2, Domain_1, Payload> &C,
                   const Range_1 &                                   f,
                   const Range_2 &                                   g)
{
  // We return the result of the compound PackagedOperation
  // directly, so as to ensure that the underlying Payload
  // definition aligns with the operations expressed here.
  // All of the memory allocations etc. are taken care of
  // internally.
  return g - C * A_inv * f;
}

/**
 * @relatesalso PackagedOperation
 * 对于方程组@f{eqnarray*}{
 * Ax + By &=& f \\
 * Cx + Dy &=& g \quad ,
 * @f}来说
 * *该操作执行舒尔补码的后处理步骤，在已知子向量 @p y 的情况下求解第二个子向量 @p x ，其结果是@f[
 * x =  A^{-1}(f
 *
 * - By)
 * @f]
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range_1,
          typename Domain_1,
          typename Domain_2,
          typename Payload>
PackagedOperation<Domain_1>
postprocess_schur_solution(
  const LinearOperator<Range_1, Domain_1, Payload> &A_inv,
  const LinearOperator<Range_1, Domain_2, Payload> &B,
  const Domain_2 &                                  y,
  const Range_1 &                                   f)
{
  // We return the result of the compound PackagedOperation
  // directly, so as to ensure that the underlying Payload
  // definition aligns with the operations expressed here.
  // All of the memory allocations etc. are taken care of
  // internally.
  return A_inv * (f - B * y);
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif


