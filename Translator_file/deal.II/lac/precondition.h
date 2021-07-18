//include/deal.II-translator/lac/precondition_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_precondition_h
#define dealii_precondition_h

// This file contains simple preconditioners.

#include <deal.II/base/config.h>

#include <deal.II/base/cuda_size.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class SparseMatrix;
namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename, typename>
    class Vector;
    template <typename>
    class BlockVector;
  } // namespace distributed
} // namespace LinearAlgebra
#endif


/*!   @addtogroup Preconditioners  
     * @{ 

 
*
*/


/**
 * 没有预设条件。
 * 如果你想使用一个没有预处理的线性求解器，这个类可以帮助你。LAC中的所有求解器都需要一个预处理。因此，你必须使用这里提供的身份来避免预处理。它可以按以下方式使用。
 *
 *
 * @code
 * SolverControl           solver_control (1000, 1e-12);
 * SolverCG<>              cg (solver_control);
 * cg.solve (system_matrix, solution, system_rhs, PreconditionIdentity());
 * @endcode
 *
 * 参见 step-3 教程程序中的一个例子和其他解释。
 * 另外，也可以用IdentityMatrix类来预设条件，方法是这样的。
 *
 *
 */
class PreconditionIdentity : public Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 这个函数的出现只是为了提供一个预处理程序的接口，以交给平滑器。
   * 这没有什么作用。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。
     *
     */
    AdditionalData() = default;
  };

  /**
   * 构造函数，将域和范围的大小设置为它们的默认值。
   *
   */
  PreconditionIdentity();

  /**
   * 矩阵参数被忽略，这里只是为了与更复杂的预处理程序兼容。
   *
   */
  template <typename MatrixType>
  void
  initialize(const MatrixType &    matrix,
             const AdditionalData &additional_data = AdditionalData());

  /**
   * 应用预处理程序。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * 应用转置预处理程序。由于这是一个身份，这个函数与vmult()相同。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;

  /**
   * 应用预处理程序，在前一个值上添加。
   *
   */
  template <class VectorType>
  void
  vmult_add(VectorType &, const VectorType &) const;

  /**
   * 应用转置的预处理程序，添加。由于这是一个身份，这个函数与vmult_add()相同。
   *
   */
  template <class VectorType>
  void
  Tvmult_add(VectorType &, const VectorType &) const;

  /**
   * 这个函数的出现只是为了提供一个预处理程序的接口，以交给平滑器。
   * 这个函数什么都不做。
   *
   */
  void
  clear()
  {}

  /**
   * 返回共域（或范围）空间的维数。注意，矩阵的维度是
   * $m \times n$  。
   * @note
   * 只有在预处理程序被初始化的情况下才能调用这个函数。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维数。请注意，矩阵的维度是 $m \times n$
   * 。
   * @note
   * 只有在预处理程序被初始化的情况下才可以调用这个函数。
   *
   */
  size_type
  n() const;

private:
  /**
   * 范围空间的维度。
   *
   */
  size_type n_rows;

  /**
   * 域空间的维度。
   *
   */
  size_type n_columns;
};



/**
 * 用Richardson方法进行预处理。这个预处理只是用AdditionalData对象提供的常数松弛因子对向量进行缩放。
 * 在Krylov空间方法中，这个预处理器不应该有任何影响。使用SolverRichardson，两个松弛参数将只是相乘。不过，这个类在多网格平滑器对象（MGSmootherRelaxation）中还是很有用的。
 *
 *
 */
class PreconditionRichardson : public Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 理查德森预处理程序的参数。
   *
   */
  class AdditionalData
  {
  public:
    /**
     * 构造器。由于没有合理的默认参数，必须给出块的大小。
     *
     */
    AdditionalData(const double relaxation = 1.);

    /**
     * 松弛参数。
     *
     */
    double relaxation;
  };

  /**
   * 构造函数，将松弛参数、域和范围大小设置为它们的默认值。
   *
   */
  PreconditionRichardson();

  /**
   * 改变松弛参数。
   *
   */
  void
  initialize(const AdditionalData &parameters);

  /**
   * 以与其他预处理程序一致的方式改变松弛参数。矩阵参数被忽略，这里只是为了与更复杂的预处理程序兼容。
   *
   */
  template <typename MatrixType>
  void
  initialize(const MatrixType &matrix, const AdditionalData &parameters);

  /**
   * 应用预处理程序。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * 应用转置预处理程序。由于这是一个身份，这个函数与vmult()相同。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;
  /**
   * 应用预处理程序，在前一个值的基础上增加。
   *
   */
  template <class VectorType>
  void
  vmult_add(VectorType &, const VectorType &) const;

  /**
   * 应用转置的预处理程序，添加。由于这是一个身份，这个函数与vmult_add()相同。
   *
   */
  template <class VectorType>
  void
  Tvmult_add(VectorType &, const VectorType &) const;

  /**
   * 这个函数的出现只是为了提供一个预处理程序的接口，以交给平滑器。
   * 这个函数什么都不做。
   *
   */
  void
  clear()
  {}

  /**
   * 返回共域（或范围）空间的维数。注意，矩阵的维度是
   * $m \times n$  。
   * @note
   * 只有在预处理程序被初始化的情况下才能调用这个函数。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维数。请注意，矩阵的维度是 $m \times n$
   * 。
   * @note
   * 只有在预处理程序被初始化的情况下才可以调用这个函数。
   *
   */
  size_type
  n() const;

private:
  /**
   * 松弛参数与向量相乘。
   *
   */
  double relaxation;

  /**
   * 范围空间的维度。
   *
   */
  size_type n_rows;

  /**
   * 域空间的维度。
   *
   */
  size_type n_columns;
};



/**
 * 使用矩阵构建函数的预处理程序。
 * 该类形成了适合于LAC求解器类的预处理方法。由于许多预处理方法是基于矩阵条目的，因此这些方法必须作为底层矩阵实现的成员函数来实现。现在这个类的目的是允许从LAC求解器类中轻松访问这些成员函数。
 * 似乎所有内置的预处理程序都有一个松弛参数，所以请为这些程序使用PreconditionRelaxation。
 * 你通常不会想要创建这种类型的命名对象，尽管有可能。最常见的用法是这样的。
 *
 * @code
 * SolverGMRES<SparseMatrix<double> Vector<double>> gmres(control,memory,500);
 *
 * gmres.solve(
 * matrix, solution, right_hand_side,
 * PreconditionUseMatrix<SparseMatrix<double>,Vector<double> >(
 *   matrix, &SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
 * 这将创建一个未命名的对象，作为第四个参数传递给SolverGMRES类的求解器函数。它假设SparseMatrix类有一个函数<tt>precondition_Jacobi</tt>，以两个向量（源和目的）作为参数（实际上，没有这样的函数，现有的函数需要第三个参数，表示松弛参数；因此这个例子只是为了说明一般的想法）。
 * 请注意，由于默认的模板参数，上面的例子可以写得更短，如下。
 *
 * @code
 * ...
 * gmres.solve(
 * matrix, solution, right_hand_side,
 * PreconditionUseMatrix<>(
 *   matrix,&SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
 *
 *
 *
 */
template <typename MatrixType = SparseMatrix<double>,
          class VectorType    = Vector<double>>
class PreconditionUseMatrix : public Subscriptor
{
public:
  /**
   * 矩阵的预处理函数的类型。
   *
   */
  using function_ptr = void (MatrixType::*)(VectorType &,
                                            const VectorType &) const;

  /**
   * 构造函数。
   * 这个构造函数存储了一个对矩阵对象的引用，供以后使用，并选择一个预处理方法，这个方法必须是该矩阵的成员函数。
   *
   */
  PreconditionUseMatrix(const MatrixType &M, const function_ptr method);

  /**
   * 执行预处理。调用传递给该对象构造函数的函数，并在此给出两个参数。
   *
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

private:
  /**
   * 指向使用中的矩阵的指针。
   *
   */
  const MatrixType &matrix;

  /**
   * 指向预处理函数的指针。
   *
   */
  const function_ptr precondition;
};



/**
 * 其他预处理程序的基类。这里，只实现了一些常见的Jacobi、SOR和SSOR预处理函数。对于预处理，请参考派生类。
 *
 *
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionRelaxation : public Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = typename MatrixType::size_type;

  /**
   * 用于参数的类。
   *
   */
  class AdditionalData
  {
  public:
    /**
     * 构造函数。
     *
     */
    AdditionalData(const double relaxation = 1.);

    /**
     * 放松参数。
     *
     */
    double relaxation;
  };

  /**
   * 初始化矩阵和松弛参数。矩阵只是存储在预处理程序对象中。松弛参数应大于零，并小于2，因为数值上的原因。它的默认值为1。
   *
   */
  void
  initialize(const MatrixType &    A,
             const AdditionalData &parameters = AdditionalData());

  /**
   * 释放矩阵并重置其指针。
   *
   */
  void
  clear();

  /**
   * 返回共域（或范围）空间的维数。注意，矩阵的维度是
   * $m \times n$  。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维度。注意，矩阵的维度是 $m \times n$  .
   *
   */
  size_type
  n() const;

protected:
  /**
   * 指向矩阵对象的指针。
   *
   */
  SmartPointer<const MatrixType, PreconditionRelaxation<MatrixType>> A;

  /**
   * 松弛参数。
   *
   */
  double relaxation;
};



/**
 * 使用矩阵内置函数的雅可比预处理。 使用的<tt>MatrixType</tt>类需要有一个<tt>precondition_Jacobi(VectorType&, const VectorType&, double</tt>)函数。这个类满足了 @ref ConceptRelaxationType 的 "松弛概念"
 * 。
 *
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner:
 *
 * PreconditionJacobi<SparseMatrix<double> > precondition;
 * precondition.initialize(
 * A, PreconditionJacobi<SparseMatrix<double>>::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 *
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionJacobi : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * 基类AdditionalData的一个别名。
   *
   */
  using AdditionalData =
    typename PreconditionRelaxation<MatrixType>::AdditionalData;

  /**
   * 应用预处理程序。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * 应用转置的预处理程序。由于这是一个对称的预处理程序，这个函数与vmult()相同。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;

  /**
   * 执行预处理的Richardson迭代的一个步骤。
   *
   */
  template <class VectorType>
  void
  step(VectorType &x, const VectorType &rhs) const;

  /**
   * 执行预设条件的Richardson迭代的一个步骤。
   *
   */
  template <class VectorType>
  void
  Tstep(VectorType &x, const VectorType &rhs) const;
};


/**
 * 使用矩阵内置函数进行SOR预处理。
 * 假设矩阵<i>A = D + L + U</i>被分割成对角线<i>D</i>以及严格的下三角和上三角<i>L</i>和<i>U</i>，那么具有松弛参数<i>r</i>的SOR预处理器是@f[
 * P^{-1} = r (D+rL)^{-1}.
 * @f] 正是这个算子<i>P<sup>-1</sup></i>，由vmult（）通过正置换实现。类似地，Tvmult()实现了<i>r(D+rU)<sup>-1</sup></i>的操作。
 * SOR迭代本身可以直接写成@f[
 * x^{k+1} = x^k
 *
 * - r D^{-1} \bigl(L x^{k+1} + U x^k
 *
 * - b\bigr).
 * @f] 使用右手边的<i>b</i>和前一个迭代<i>x</i>，这是由step()实现的操作。
 * 使用的MatrixType类需要有<tt>precondition_SOR(VectorType&, const VectorType&, double)</tt>和<tt>precondition_TSOR(VectorType&, const VectorType&, double)</tt>函数。这个类满足了 @ref ConceptRelaxationType 的 "放松概念"
 * 。
 *
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner
 *
 * PreconditionSOR<SparseMatrix<double> > precondition;
 * precondition.initialize(
 * A, PreconditionSOR<SparseMatrix<double>>::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 *
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionSOR : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * 基类AdditionalData的一个别名。
   *
   */
  using AdditionalData =
    typename PreconditionRelaxation<MatrixType>::AdditionalData;

  /**
   * 应用预处理程序。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * 应用转置的预处理程序。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;

  /**
   * 执行预处理的Richardson迭代的一个步骤。
   *
   */
  template <class VectorType>
  void
  step(VectorType &x, const VectorType &rhs) const;

  /**
   * 执行一个转置的预处理Richardson迭代步骤。
   *
   */
  template <class VectorType>
  void
  Tstep(VectorType &x, const VectorType &rhs) const;
};



/**
 * 使用矩阵内置函数的SSOR预处理。 使用的<tt>MatrixType</tt>类需要有一个函数<tt>precondition_SSOR(VectorType&, const VectorType&, double)</tt>。这个类满足了 @ref ConceptRelaxationType 的 "放松概念"
 * 。
 *
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner
 *
 * PreconditionSSOR<SparseMatrix<double> > precondition;
 * precondition.initialize(
 * A, PreconditionSSOR<SparseMatrix<double>>::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 *
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionSSOR : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * 基类AdditionalData的一个别名。
   *
   */
  using AdditionalData =
    typename PreconditionRelaxation<MatrixType>::AdditionalData;

  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = typename MatrixType::size_type;

  /**
   * 对基类的一个别名。
   *
   */
  using BaseClass = PreconditionRelaxation<MatrixType>;


  /**
   * 初始化矩阵和放松参数。矩阵只是存储在预处理器对象中。松弛参数应该大于零，并小于2，因为数值上的原因。它的默认值为1。
   *
   */
  void
  initialize(const MatrixType &                        A,
             const typename BaseClass::AdditionalData &parameters =
               typename BaseClass::AdditionalData());

  /**
   * 应用预处理程序。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * 应用转置预处理程序。由于这是一个对称的预处理程序，这个函数与vmult()相同。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;


  /**
   * 执行预处理的Richardson迭代的一个步骤
   *
   */
  template <class VectorType>
  void
  step(VectorType &x, const VectorType &rhs) const;

  /**
   * 执行预设条件的Richardson迭代的一个步骤。
   *
   */
  template <class VectorType>
  void
  Tstep(VectorType &x, const VectorType &rhs) const;

private:
  /**
   * 一个数组，为每个矩阵行存储对角线后第一个位置的位置。
   *
   */
  std::vector<std::size_t> pos_right_of_diagonal;
};


/**
 * 使用矩阵内置函数的修正SOR预处理程序。
 * 使用的<tt>MatrixType</tt>类需要有<tt>PSOR(VectorType&, const
 * VectorType&, double)</tt>和<tt>TPSOR(VectorType&, const VectorType&,
 * double)</tt>函数。
 *
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * std::vector<unsigned int> permutation(x.size());
 * std::vector<unsigned int> inverse_permutation(x.size());
 *
 * //...fill permutation and its inverse with reasonable values
 *
 * // Define and initialize preconditioner
 *
 * PreconditionPSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, permutation, inverse_permutation, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionPSOR : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = typename MatrixType::size_type;

  /**
   * PreconditionPSOR的参数。
   *
   */
  class AdditionalData
  {
  public:
    /**
     * 构造函数。关于参数的描述，见下文。
     * 包容向量是作为一个参考来存储的。因此，必须保证向量的寿命超过预处理程序的寿命。
     * 松弛参数应该大于零，并且由于数值原因小于2。它的默认值是1。
     *
     */
    AdditionalData(
      const std::vector<size_type> &permutation,
      const std::vector<size_type> &inverse_permutation,
      const typename PreconditionRelaxation<MatrixType>::AdditionalData
        &parameters =
          typename PreconditionRelaxation<MatrixType>::AdditionalData());

    /**
     * 包容向量的存储。
     *
     */
    const std::vector<size_type> &permutation;
    /**
     * 储存反互换向量。
     *
     */
    const std::vector<size_type> &inverse_permutation;
    /**
     * 放松参数
     *
     */
    typename PreconditionRelaxation<MatrixType>::AdditionalData parameters;
  };

  /**
   * 初始化矩阵和松弛参数。矩阵只是存储在前置条件器对象中。
   * 扰动向量以指针形式存储。因此，必须保证该向量的寿命超过预调节器的寿命。
   * 松弛参数应该大于零，并且由于数值原因小于2。它的默认值是1。
   *
   */
  void
  initialize(const MatrixType &            A,
             const std::vector<size_type> &permutation,
             const std::vector<size_type> &inverse_permutation,
             const typename PreconditionRelaxation<MatrixType>::AdditionalData
               &parameters =
                 typename PreconditionRelaxation<MatrixType>::AdditionalData());

  /**
   * 初始化矩阵和松弛参数。矩阵只是存储在preconditioner对象中。
   * 关于可能的参数的更多细节，请参阅类的文档和
   * PreconditionPSOR::AdditionalData 类的文档。
   * 这个函数被调用后，预处理程序就可以使用了（使用派生类的
   * <code>vmult</code> 函数）。
   *
   */
  void
  initialize(const MatrixType &A, const AdditionalData &additional_data);

  /**
   * 应用预处理程序。
   *
   */
  template <class VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * 应用转置预处理程序。
   *
   */
  template <class VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;

private:
  /**
   * 储存互换向量。
   *
   */
  const std::vector<size_type> *permutation;
  /**
   * 储存反置换向量。
   *
   */
  const std::vector<size_type> *inverse_permutation;
};



/**
 * 用Chebyshev多项式对对称正定矩阵进行预处理。这个预处理是基于 @p PreconditionerType 类型的内部预处理的迭代，其系数被调整为最佳覆盖最大特征值 $\lambda_{\max{}}$ 到由可选参数 @p smoothing_range. 指定的特定较低特征值 $\lambda_{\min{}}$ 之间的特征值范围。该算法是基于以下三期递归。@f[
 * x^{n+1} = x^{n} + \rho_n \rho_{n-1} (x^{n}
 *
 * - x^{n-1}) +
 *   \frac{\rho_n}{\lambda_{\max{}}-\lambda_{\min{}}} P^{-1} (b-Ax^n).
 * @f] 其中参数 $\rho_0$ 被设置为 $\rho_0 = 2
 * \frac{\lambda_{\max{}}-\lambda_{\min{}}}{\lambda_{\max{}}+\lambda_{\min{}}}$
 * 的最大特征值 $\lambda_{\max{}}$ 并通过 $\rho_n =
 * \left(2\frac{\lambda_{\max{}}+\lambda_{\min{}}}
 * {\lambda_{\max{}}-\lambda_{\min{}}}
 *
 * - \rho_{n-1}\right)^{-1}$
 * 更新。切比雪夫多项式的构造是为了强烈阻尼
 * $\lambda_{\min{}}$ 和 $\lambda_{\max{}}$
 * 之间的特征值范围，例如，在
 * Utilities::LinearAlgebra::chebyshev_filter(). 中可以看到。
 * 预处理程序的典型用例是通过DiagonalMatrix指定的雅可比预处理程序，这也是预处理程序的默认值。请注意，如果程度变量被设置为1，切比雪夫迭代就对应于一个雅可比预处理器（或基础预处理器类型），其松弛参数根据指定的平滑范围。
 * 除了默认选择的点式雅可比预处理外，该类还允许更高级的预处理类型，例如DG方法中的迭代块-雅可比预处理。
 * 除了内部预处理对象，这种迭代不需要访问矩阵条目，这使得它成为无矩阵计算的理想成分。在这种情况下，这个类可以被用作多网格平滑器，它实际上是%并行的（假设矩阵-向量乘积是%并行的，内部预处理器是%并行的）。在
 * step-37 和 step-59 的教程程序中演示了它的使用。
 * <h4>Estimation of the eigenvalues</h4>
 * 切比雪夫方法依赖于对矩阵特征值的估计，该估计在第一次调用vmult()时被计算出来。该算法调用共轭梯度求解器（即Lanczos迭代），因此要求（预处理）矩阵系统的对称性和正确定性。特征值算法可以通过
 * PreconditionChebyshev::AdditionalData::eig_cg_n_iterations
 * 控制，指定应该进行多少次迭代。迭代从一个初始向量开始，这个初始向量取决于向量类型。对于具有快速元素访问的
 * dealii::Vector 或 dealii::LinearAlgebra::distributed::Vector,
 * 类，它是一个条目为`(-5.5,
 *
 * - .5,
 *
 * - .5,
 *
 *
 *
 * - .5, ..., 3.5, 4.5, 5.5)`，有适当的尾声，并进行调整，使其平均值始终为零，这对拉普拉斯的效果很好。这种设置在并行中是稳定的，因为对于不同数量的处理器，但未知数的排序相同，除了舍入误差外，将计算相同的初始向量，从而计算特征值分布。对于其他的矢量类型，初始矢量包含所有的1，按矢量的长度缩放，除了最开始的条目是0，再次触发高频内容。
 * 特征值的计算发生在第一次调用vmult()、Tvmult()、step()或Tstep()函数之一时，或者直接调用improve_eigenvalues()时。在后一种情况下，有必要提供一个临时向量，其布局与应用预处理程序时使用的源向量和目的向量相同。
 * 最小和最大特征值的估计值取自SolverCG（即使求解器没有在要求的迭代次数中收敛）。最后，最大特征值被乘以1.2的安全系数。
 * 由于特征值估计的成本，如果反复应用，例如在几何多网格求解器的平滑器中，该类最适合，而该平滑器又可用于解决几个线性系统。
 * <h4>Bypassing the eigenvalue computation</h4>
 * 在某些情况下，这一类的自动特征值计算可能会导致质量不好，或者在与不同的自由度枚举并行使用时可能不稳定，使计算强烈依赖于并行配置。可以通过设置
 * AdditionalData::eig_cg_n_iterations
 * 为零来绕过自动特征值计算，而提供变量
 * AdditionalData::max_eigenvalue
 * 来代替。最小特征值是通过`max_eigenvalue/smoothing_range`隐式指定的。
 * <h4>Using the PreconditionChebyshev as a solver</h4>
 * 如果范围<tt>[max_eigenvalue/smoothing_range,
 * max_eigenvalue]</tt>包含预处理矩阵系统的所有特征值，并且程度（即迭代次数）足够高，这个类也可以作为一个直接求解器使用。关于Chebyshev迭代的误差估计，可用于确定迭代次数，见Varga（2009）。
 * 为了使用切比雪夫作为求解器，将程度设置为
 * numbers::invalid_unsigned_int
 * 以强制自动计算达到给定目标公差所需的迭代次数。在这种情况下，目标公差从变量
 * PreconditionChebyshev::AdditionalData::smoothing_range
 * 中读取（它需要是一个小于1的数字，以强制进行任何迭代，显然）。
 * 关于该算法的细节，请参见第5.1节。
 *
 * @code{.bib}
 * @book{Varga2009,
 * Title     = {Matrix iterative analysis},
 * Author    = {Varga, R. S.},
 * Publisher = {Springer},
 * Address   = {Berlin},
 * Edition   = {2nd},
 * Year      = {2009},
 * }
 * @endcode
 *  <h4>Requirements on the templated classes</h4>
 * MatrixType类必须从Subscriptor派生出来，因为MatrixType的SmartPointer被保存在该类中。特别是，这意味着矩阵对象需要在PreconditionChebyshev的生命周期内持续存在。先决条件被保存在一个shared_ptr中，被复制到类的AdditionalData成员变量中，所以用于初始化的变量在调用initialize()后可以安全地被丢弃。矩阵和预处理程序都需要为矩阵-向量乘积提供
 * @p vmult() 函数，为访问（方形）矩阵的行数提供 @p m()
 * 函数。此外，矩阵必须提供<tt>el(i,i)</tt>方法，以便在预处理程序类型为DiagonalMatrix时访问矩阵对角线。尽管强烈建议在一个单独的预处理对象中传递反向对角线条目，以实现Jacobi方法（这是用MPI进行%并行计算时操作该类的唯一可能方式，因为没有关于本地存储的条目范围的知识，只需要从矩阵中获取），但有一个向后兼容的函数，可以在串行计算的情况下提取对角线。
 *
 *
 */
template <typename MatrixType         = SparseMatrix<double>,
          typename VectorType         = Vector<double>,
          typename PreconditionerType = DiagonalMatrix<VectorType>>
class PreconditionChebyshev : public Subscriptor
{
public:
  /**
   * 申报容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 标准化的数据结构，用于向预处理程序输送附加参数。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造器。
     *
     */
    AdditionalData(const unsigned int degree              = 1,
                   const double       smoothing_range     = 0.,
                   const unsigned int eig_cg_n_iterations = 8,
                   const double       eig_cg_residual     = 1e-2,
                   const double       max_eigenvalue      = 1);

    /**
     * 复制赋值运算符。
     *
     */
    AdditionalData &
    operator=(const AdditionalData &other_data);

    /**
     * 这决定了切比雪夫多项式的程度。多项式的度数给出了应用vmult()操作所要进行的矩阵-向量乘积的数量。度数1对应的是阻尼雅可比方法。
     * 如果度数被设置为 numbers::invalid_unsigned_int,
     * ，算法将根据主类讨论中提到的通常的切比雪夫误差公式自动确定必要的迭代次数。
     *
     */
    unsigned int degree;

    /**
     * 这设置了矩阵中最大的特征值和最小的特征值之间的处理范围。如果该参数被设置为一个小于1的数字，最大和最小的特征值的估计值将在内部计算。对于一个大于1的平滑范围，切比雪夫多项式将作用于区间
     * $[\lambda_\mathrm{max}/ \tt{smoothing\_range}, \lambda_\mathrm{max}]$
     * ，其中 $\lambda_\mathrm{max}$
     * 是对矩阵最大特征值的估计。选择<tt>smoothing_range</tt>在5和20之间是有用的，以防预处理程序被用作多网格的平滑器。
     *
     */
    double smoothing_range;

    /**
     * 为寻找最大特征值而进行的最大CG迭代次数。如果设置为零，则不进行计算。相反，用户必须通过变量
     * PreconditionChebyshev::AdditionalData::max_eigenvalue.
     * 提供一个最大特征值。
     *
     */
    unsigned int eig_cg_n_iterations;

    /**
     * 为寻找最大特征值而进行的CG迭代的容忍度。
     *
     */
    double eig_cg_residual;

    /**
     * 用于工作的最大特征值。只有在 @p
     * eig_cg_n_iterations被设置为0时才有效，否则该参数被忽略。
     *
     */
    double max_eigenvalue;

    /**
     * 将用于给定的运算器的约束条件。这个变量用于在创建初始猜测时将正确的条目清零。
     *
     */
    AffineConstraints<double> constraints;

    /**
     * 存储切比雪夫所包涵的预处理对象。
     *
     */
    std::shared_ptr<PreconditionerType> preconditioner;
  };


  /**
   * 构造函数。
   *
   */
  PreconditionChebyshev();

  /**
   * 初始化函数。接受用于形成预处理程序的矩阵，如果有额外的标志，则接受额外的标志。这个函数只有在输入矩阵有一个操作符<tt>el(i,i)</tt>用于访问对角线上的所有元素时才能工作。或者，可以在AdditionalData字段的帮助下提供对角线。
   * 这个函数在给定的迭代次数为正数的情况下，使用修改过的CG迭代法计算矩阵的特征值范围，并以其对角线为权重进行估算。
   *
   */
  void
  initialize(const MatrixType &    matrix,
             const AdditionalData &additional_data = AdditionalData());

  /**
   * 计算预处理程序对<tt>src</tt>的作用，将结果存储在<tt>dst</tt>。
   *
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * 计算转置的预处理程序对<tt>src</tt>的作用，将结果存入<tt>dst</tt>。
   *
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * 执行预处理Richardson迭代的一个步骤。
   *
   */
  void
  step(VectorType &dst, const VectorType &src) const;

  /**
   * 执行预设条件的Richardson迭代的一个步骤。
   *
   */
  void
  Tstep(VectorType &dst, const VectorType &src) const;

  /**
   * 重置预处理程序。
   *
   */
  void
  clear();

  /**
   * 返回码域（或范围）空间的维数。注意，矩阵的维度为
   * $m \times n$  。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维度。注意该矩阵的维度为 $m \times n$  .
   *
   */
  size_type
  n() const;

  /**
   * 一个包含由PreconditionChebychev类执行的特征值估计信息的结构。
   *
   */
  struct EigenvalueInformation
  {
    /**
     * 最小特征值的估计。
     *
     */
    double min_eigenvalue_estimate;
    /**
     * 对最大特征值的估计。
     *
     */
    double max_eigenvalue_estimate;
    /**
     * 执行的CG迭代次数或0。
     *
     */
    unsigned int cg_iterations;
    /**
     * 切比雪夫多项式的度数（使用 AdditionalData::degree
     * 设置或按该处所述估计）。
     *
     */
    unsigned int degree;
    /**
     * 用无效值初始化的构造函数。
     *
     */
    EigenvalueInformation()
      : min_eigenvalue_estimate{std::numeric_limits<double>::max()}
      , max_eigenvalue_estimate{std::numeric_limits<double>::lowest()}
      , cg_iterations{0}
      , degree{0}
    {}
  };

  /**
   * 计算预处理程序所需的特征值估计。
   * 如果用户没有调用该函数，则在第一次使用预处理程序时自动调用。向量
   * @p src 的布局用于创建内部临时向量，其内容并不重要。
   * 基于特征值的计算，初始化因子theta和delta。如果用户在AdditionalData中为最大的特征值设置了值，则不进行计算，而使用用户给出的信息。
   *
   */
  EigenvalueInformation
  estimate_eigenvalues(const VectorType &src) const;

private:
  /**
   * 一个指向底层矩阵的指针。
   *
   */
  SmartPointer<
    const MatrixType,
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>>
    matrix_ptr;

  /**
   * 用于<tt>vmult</tt>操作的内部向量。
   *
   */
  mutable VectorType solution_old;

  /**
   * 用于<tt>vmult</tt>操作的内部向量。
   *
   */
  mutable VectorType temp_vector1;

  /**
   * 用于<tt>vmult</tt>操作的内部向量。
   *
   */
  mutable VectorType temp_vector2;

  /**
   * 存储传递给初始化函数的额外数据，通过复制操作获得。
   *
   */
  AdditionalData data;

  /**
   * 所考虑的最大和最小的特征值的平均值。
   *
   */
  double theta;

  /**
   * 所考虑的最大和最小特征值之间的一半区间长度。
   *
   */
  double delta;

  /**
   * 存储预处理程序是否已被设置，特征值是否已被计算。
   *
   */
  bool eigenvalues_are_initialized;

  /**
   * 一个mutex，以避免不同线程的多个vmult()调用覆盖了临时向量。
   *
   */
  mutable Threads::Mutex mutex;
};



 /*@}*/ 
 /* ---------------------------------- Inline functions ------------------- */ 

#ifndef DOXYGEN

inline PreconditionIdentity::PreconditionIdentity()
  : n_rows(0)
  , n_columns(0)
{}

template <typename MatrixType>
inline void
PreconditionIdentity::initialize(const MatrixType &matrix,
                                 const PreconditionIdentity::AdditionalData &)
{
  n_rows    = matrix.m();
  n_columns = matrix.n();
}


template <class VectorType>
inline void
PreconditionIdentity::vmult(VectorType &dst, const VectorType &src) const
{
  dst = src;
}



template <class VectorType>
inline void
PreconditionIdentity::Tvmult(VectorType &dst, const VectorType &src) const
{
  dst = src;
}

template <class VectorType>
inline void
PreconditionIdentity::vmult_add(VectorType &dst, const VectorType &src) const
{
  dst += src;
}



template <class VectorType>
inline void
PreconditionIdentity::Tvmult_add(VectorType &dst, const VectorType &src) const
{
  dst += src;
}

inline PreconditionIdentity::size_type
PreconditionIdentity::m() const
{
  Assert(n_rows != 0, ExcNotInitialized());
  return n_rows;
}

inline PreconditionIdentity::size_type
PreconditionIdentity::n() const
{
  Assert(n_columns != 0, ExcNotInitialized());
  return n_columns;
}

//---------------------------------------------------------------------------

inline PreconditionRichardson::AdditionalData::AdditionalData(
  const double relaxation)
  : relaxation(relaxation)
{}


inline PreconditionRichardson::PreconditionRichardson()
  : relaxation(0)
  , n_rows(0)
  , n_columns(0)
{
  AdditionalData add_data;
  relaxation = add_data.relaxation;
}



inline void
PreconditionRichardson::initialize(
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
}



template <typename MatrixType>
inline void
PreconditionRichardson::initialize(
  const MatrixType &                            matrix,
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
  n_rows     = matrix.m();
  n_columns  = matrix.n();
}



template <class VectorType>
inline void
PreconditionRichardson::vmult(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.equ(relaxation, src);
}



template <class VectorType>
inline void
PreconditionRichardson::Tvmult(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.equ(relaxation, src);
}

template <class VectorType>
inline void
PreconditionRichardson::vmult_add(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.add(relaxation, src);
}



template <class VectorType>
inline void
PreconditionRichardson::Tvmult_add(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.add(relaxation, src);
}

inline PreconditionRichardson::size_type
PreconditionRichardson::m() const
{
  Assert(n_rows != 0, ExcNotInitialized());
  return n_rows;
}

inline PreconditionRichardson::size_type
PreconditionRichardson::n() const
{
  Assert(n_columns != 0, ExcNotInitialized());
  return n_columns;
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionRelaxation<MatrixType>::initialize(const MatrixType &    rA,
                                               const AdditionalData &parameters)
{
  A          = &rA;
  relaxation = parameters.relaxation;
}


template <typename MatrixType>
inline void
PreconditionRelaxation<MatrixType>::clear()
{
  A = nullptr;
}

template <typename MatrixType>
inline typename PreconditionRelaxation<MatrixType>::size_type
PreconditionRelaxation<MatrixType>::m() const
{
  Assert(A != nullptr, ExcNotInitialized());
  return A->m();
}

template <typename MatrixType>
inline typename PreconditionRelaxation<MatrixType>::size_type
PreconditionRelaxation<MatrixType>::n() const
{
  Assert(A != nullptr, ExcNotInitialized());
  return A->n();
}

//---------------------------------------------------------------------------

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::vmult(VectorType &      dst,
                                      const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_Jacobi(dst, src, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::Tvmult(VectorType &      dst,
                                       const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_Jacobi(dst, src, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::step(VectorType &      dst,
                                     const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->Jacobi_step(dst, src, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::Tstep(VectorType &      dst,
                                      const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  step(dst, src);
}



//---------------------------------------------------------------------------

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::vmult(VectorType &dst, const VectorType &src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_SOR(dst, src, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::Tvmult(VectorType &      dst,
                                    const VectorType &src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_TSOR(dst, src, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::step(VectorType &dst, const VectorType &src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->SOR_step(dst, src, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::Tstep(VectorType &dst, const VectorType &src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->TSOR_step(dst, src, this->relaxation);
}



//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionSSOR<MatrixType>::initialize(
  const MatrixType &                        rA,
  const typename BaseClass::AdditionalData &parameters)
{
  this->PreconditionRelaxation<MatrixType>::initialize(rA, parameters);

  // in case we have a SparseMatrix class, we can extract information about
  // the diagonal.
  const SparseMatrix<typename MatrixType::value_type> *mat =
    dynamic_cast<const SparseMatrix<typename MatrixType::value_type> *>(
      &*this->A);

  // calculate the positions first after the diagonal.
  if (mat != nullptr)
    {
      const size_type n = this->A->n();
      pos_right_of_diagonal.resize(n, static_cast<std::size_t>(-1));
      for (size_type row = 0; row < n; ++row)
        {
          // find the first element in this line which is on the right of the
          // diagonal.  we need to precondition with the elements on the left
          // only. note: the first entry in each line denotes the diagonal
          // element, which we need not check.
          typename SparseMatrix<typename MatrixType::value_type>::const_iterator
            it = mat->begin(row) + 1;
          for (; it < mat->end(row); ++it)
            if (it->column() > row)
              break;
          pos_right_of_diagonal[row] = it - mat->begin();
        }
    }
}


template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::vmult(VectorType &      dst,
                                    const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_SSOR(dst, src, this->relaxation, pos_right_of_diagonal);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::Tvmult(VectorType &      dst,
                                     const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_SSOR(dst, src, this->relaxation, pos_right_of_diagonal);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::step(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->SSOR_step(dst, src, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::Tstep(VectorType &      dst,
                                    const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  step(dst, src);
}



//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionPSOR<MatrixType>::initialize(
  const MatrixType &                                                 rA,
  const std::vector<size_type> &                                     p,
  const std::vector<size_type> &                                     ip,
  const typename PreconditionRelaxation<MatrixType>::AdditionalData &parameters)
{
  permutation         = &p;
  inverse_permutation = &ip;
  PreconditionRelaxation<MatrixType>::initialize(rA, parameters);
}


template <typename MatrixType>
inline void
PreconditionPSOR<MatrixType>::initialize(const MatrixType &    A,
                                         const AdditionalData &additional_data)
{
  initialize(A,
             additional_data.permutation,
             additional_data.inverse_permutation,
             additional_data.parameters);
}


template <typename MatrixType>
template <typename VectorType>
inline void
PreconditionPSOR<MatrixType>::vmult(VectorType &      dst,
                                    const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionPSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionPSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  dst = src;
  this->A->PSOR(dst, *permutation, *inverse_permutation, this->relaxation);
}



template <typename MatrixType>
template <class VectorType>
inline void
PreconditionPSOR<MatrixType>::Tvmult(VectorType &      dst,
                                     const VectorType &src) const
{
  static_assert(
    std::is_same<typename PreconditionPSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionPSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  dst = src;
  this->A->TPSOR(dst, *permutation, *inverse_permutation, this->relaxation);
}

template <typename MatrixType>
PreconditionPSOR<MatrixType>::AdditionalData::AdditionalData(
  const std::vector<size_type> &permutation,
  const std::vector<size_type> &inverse_permutation,
  const typename PreconditionRelaxation<MatrixType>::AdditionalData &parameters)
  : permutation(permutation)
  , inverse_permutation(inverse_permutation)
  , parameters(parameters)
{}


//---------------------------------------------------------------------------


template <typename MatrixType, class VectorType>
PreconditionUseMatrix<MatrixType, VectorType>::PreconditionUseMatrix(
  const MatrixType & M,
  const function_ptr method)
  : matrix(M)
  , precondition(method)
{}



template <typename MatrixType, class VectorType>
void
PreconditionUseMatrix<MatrixType, VectorType>::vmult(
  VectorType &      dst,
  const VectorType &src) const
{
  (matrix.*precondition)(dst, src);
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline PreconditionRelaxation<MatrixType>::AdditionalData::AdditionalData(
  const double relaxation)
  : relaxation(relaxation)
{}



//---------------------------------------------------------------------------

namespace internal
{
  namespace PreconditionChebyshevImplementation
  {
    // for deal.II vectors, perform updates for Chebyshev preconditioner all
    // at once to reduce memory transfer. Here, we select between general
    // vectors and deal.II vectors where we expand the loop over the (local)
    // size of the vector

    // generic part for non-deal.II vectors
    template <typename VectorType, typename PreconditionerType>
    inline void
    vector_updates(const VectorType &        rhs,
                   const PreconditionerType &preconditioner,
                   const unsigned int        iteration_index,
                   const double              factor1,
                   const double              factor2,
                   VectorType &              solution_old,
                   VectorType &              temp_vector1,
                   VectorType &              temp_vector2,
                   VectorType &              solution)
    {
      if (iteration_index == 0)
        {
          solution.equ(factor2, rhs);
          preconditioner.vmult(solution_old, solution);
        }
      else if (iteration_index == 1)
        {
          // compute t = P^{-1} * (b-A*x^{n})
          temp_vector1.sadd(-1.0, 1.0, rhs);
          preconditioner.vmult(solution_old, temp_vector1);

          // compute x^{n+1} = x^{n} + f_1 * x^{n} + f_2 * t
          solution_old.sadd(factor2, 1 + factor1, solution);
        }
      else
        {
          // compute t = P^{-1} * (b-A*x^{n})
          temp_vector1.sadd(-1.0, 1.0, rhs);
          preconditioner.vmult(temp_vector2, temp_vector1);

          // compute x^{n+1} = x^{n} + f_1 * (x^{n}-x^{n-1}) + f_2 * t
          solution_old.sadd(-factor1, factor2, temp_vector2);
          solution_old.add(1 + factor1, solution);
        }

      solution.swap(solution_old);
    }

    // worker routine for deal.II vectors. Because of vectorization, we need
    // to put the loop into an extra structure because the virtual function of
    // VectorUpdatesRange prevents the compiler from applying vectorization.
    template <typename Number>
    struct VectorUpdater
    {
      VectorUpdater(const Number *     rhs,
                    const Number *     matrix_diagonal_inverse,
                    const unsigned int iteration_index,
                    const Number       factor1,
                    const Number       factor2,
                    Number *           solution_old,
                    Number *           tmp_vector,
                    Number *           solution)
        : rhs(rhs)
        , matrix_diagonal_inverse(matrix_diagonal_inverse)
        , iteration_index(iteration_index)
        , factor1(factor1)
        , factor2(factor2)
        , solution_old(solution_old)
        , tmp_vector(tmp_vector)
        , solution(solution)
      {}

      void
      apply_to_subrange(const std::size_t begin, const std::size_t end) const
      {
        // To circumvent a bug in gcc
        // (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63945), we create
        // copies of the variables factor1 and factor2 and do not check based on
        // factor1.
        const Number factor1        = this->factor1;
        const Number factor1_plus_1 = 1. + this->factor1;
        const Number factor2        = this->factor2;
        if (iteration_index == 0)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (std::size_t i = begin; i < end; ++i)
              solution[i] = factor2 * matrix_diagonal_inverse[i] * rhs[i];
          }
        else if (iteration_index == 1)
          {
            // x^{n+1} = x^{n} + f_1 * x^{n} + f_2 * P^{-1} * (b-A*x^{n})
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (std::size_t i = begin; i < end; ++i)
              // for efficiency reason, write back to temp_vector that is
              // already read (avoid read-for-ownership)
              tmp_vector[i] =
                factor1_plus_1 * solution[i] +
                factor2 * matrix_diagonal_inverse[i] * (rhs[i] - tmp_vector[i]);
          }
        else
          {
            // x^{n+1} = x^{n} + f_1 * (x^{n}-x^{n-1})
            //           + f_2 * P^{-1} * (b-A*x^{n})
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (std::size_t i = begin; i < end; ++i)
              solution_old[i] =
                factor1_plus_1 * solution[i] - factor1 * solution_old[i] +
                factor2 * matrix_diagonal_inverse[i] * (rhs[i] - tmp_vector[i]);
          }
      }

      const Number *     rhs;
      const Number *     matrix_diagonal_inverse;
      const unsigned int iteration_index;
      const Number       factor1;
      const Number       factor2;
      mutable Number *   solution_old;
      mutable Number *   tmp_vector;
      mutable Number *   solution;
    };

    template <typename Number>
    struct VectorUpdatesRange : public ::dealii::parallel::ParallelForInteger
    {
      VectorUpdatesRange(const VectorUpdater<Number> &updater,
                         const std::size_t            size)
        : updater(updater)
      {
        if (size < internal::VectorImplementation::minimum_parallel_grain_size)
          VectorUpdatesRange::apply_to_subrange(0, size);
        else
          apply_parallel(
            0,
            size,
            internal::VectorImplementation::minimum_parallel_grain_size);
      }

      ~VectorUpdatesRange() override = default;

      virtual void
      apply_to_subrange(const std::size_t begin,
                        const std::size_t end) const override
      {
        updater.apply_to_subrange(begin, end);
      }

      const VectorUpdater<Number> &updater;
    };

    // selection for diagonal matrix around deal.II vector
    template <typename Number>
    inline void
    vector_updates(const ::dealii::Vector<Number> &                rhs,
                   const DiagonalMatrix<::dealii::Vector<Number>> &jacobi,
                   const unsigned int        iteration_index,
                   const double              factor1,
                   const double              factor2,
                   ::dealii::Vector<Number> &solution_old,
                   ::dealii::Vector<Number> &temp_vector1,
                   ::dealii::Vector<Number> &,
                   ::dealii::Vector<Number> &solution)
    {
      VectorUpdater<Number> upd(rhs.begin(),
                                jacobi.get_vector().begin(),
                                iteration_index,
                                factor1,
                                factor2,
                                solution_old.begin(),
                                temp_vector1.begin(),
                                solution.begin());
      VectorUpdatesRange<Number>(upd, rhs.size());

      // swap vectors x^{n+1}->x^{n}, given the updates in the function above
      if (iteration_index == 0)
        {
          // nothing to do here because we can immediately write into the
          // solution vector without remembering any of the other vectors
        }
      else if (iteration_index == 1)
        {
          solution.swap(temp_vector1);
          solution_old.swap(temp_vector1);
        }
      else
        solution.swap(solution_old);
    }

    // selection for diagonal matrix around parallel deal.II vector
    template <typename Number>
    inline void
    vector_updates(
      const LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &rhs,
      const DiagonalMatrix<
        LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>> &jacobi,
      const unsigned int iteration_index,
      const double       factor1,
      const double       factor2,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &solution_old,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &temp_vector1,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &solution)
    {
      VectorUpdater<Number> upd(rhs.begin(),
                                jacobi.get_vector().begin(),
                                iteration_index,
                                factor1,
                                factor2,
                                solution_old.begin(),
                                temp_vector1.begin(),
                                solution.begin());
      VectorUpdatesRange<Number>(upd, rhs.locally_owned_size());

      // swap vectors x^{n+1}->x^{n}, given the updates in the function above
      if (iteration_index == 0)
        {
          // nothing to do here because we can immediately write into the
          // solution vector without remembering any of the other vectors
        }
      else if (iteration_index == 1)
        {
          solution.swap(temp_vector1);
          solution_old.swap(temp_vector1);
        }
      else
        solution.swap(solution_old);
    }

    template <typename MatrixType, typename PreconditionerType>
    inline void
    initialize_preconditioner(
      const MatrixType &                   matrix,
      std::shared_ptr<PreconditionerType> &preconditioner)
    {
      (void)matrix;
      (void)preconditioner;
      AssertThrow(preconditioner.get() != nullptr, ExcNotInitialized());
    }

    template <typename MatrixType, typename VectorType>
    inline void
    initialize_preconditioner(
      const MatrixType &                           matrix,
      std::shared_ptr<DiagonalMatrix<VectorType>> &preconditioner)
    {
      if (preconditioner.get() == nullptr || preconditioner->m() != matrix.m())
        {
          if (preconditioner.get() == nullptr)
            preconditioner = std::make_shared<DiagonalMatrix<VectorType>>();

          Assert(
            preconditioner->m() == 0,
            ExcMessage(
              "Preconditioner appears to be initialized but not sized correctly"));

          // This part only works in serial
          if (preconditioner->m() != matrix.m())
            {
              preconditioner->get_vector().reinit(matrix.m());
              for (typename VectorType::size_type i = 0; i < matrix.m(); ++i)
                preconditioner->get_vector()(i) = 1. / matrix.el(i, i);
            }
        }
    }

    template <typename VectorType>
    void
    set_initial_guess(VectorType &vector)
    {
      vector = 1. / std::sqrt(static_cast<double>(vector.size()));
      if (vector.locally_owned_elements().is_element(0))
        vector(0) = 0.;
    }

    template <typename Number>
    void
    set_initial_guess(::dealii::Vector<Number> &vector)
    {
      // Choose a high-frequency mode consisting of numbers between 0 and 1
      // that is cheap to compute (cheaper than random numbers) but avoids
      // obviously re-occurring numbers in multi-component systems by choosing
      // a period of 11
      for (unsigned int i = 0; i < vector.size(); ++i)
        vector(i) = i % 11;

      const Number mean_value = vector.mean_value();
      vector.add(-mean_value);
    }

    template <typename Number>
    void
    set_initial_guess(
      ::dealii::LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &vector)
    {
      // Choose a high-frequency mode consisting of numbers between 0 and 1
      // that is cheap to compute (cheaper than random numbers) but avoids
      // obviously re-occurring numbers in multi-component systems by choosing
      // a period of 11.
      // Make initial guess robust with respect to number of processors
      // by operating on the global index.
      types::global_dof_index first_local_range = 0;
      if (!vector.locally_owned_elements().is_empty())
        first_local_range = vector.locally_owned_elements().nth_index_in_set(0);
      for (unsigned int i = 0; i < vector.locally_owned_size(); ++i)
        vector.local_element(i) = (i + first_local_range) % 11;

      const Number mean_value = vector.mean_value();
      vector.add(-mean_value);
    }

    template <typename Number>
    void
    set_initial_guess(
      ::dealii::LinearAlgebra::distributed::BlockVector<Number> &vector)
    {
      for (unsigned int block = 0; block < vector.n_blocks(); ++block)
        set_initial_guess(vector.block(block));
    }


#  ifdef DEAL_II_COMPILER_CUDA_AWARE
    template <typename Number>
    __global__ void
    set_initial_guess_kernel(const types::global_dof_index offset,
                             const unsigned int            locally_owned_size,
                             Number *                      values)

    {
      const unsigned int index = threadIdx.x + blockDim.x * blockIdx.x;
      if (index < locally_owned_size)
        values[index] = (index + offset) % 11;
    }

    template <typename Number>
    void
    set_initial_guess(
      ::dealii::LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA>
        &vector)
    {
      // Choose a high-frequency mode consisting of numbers between 0 and 1
      // that is cheap to compute (cheaper than random numbers) but avoids
      // obviously re-occurring numbers in multi-component systems by choosing
      // a period of 11.
      // Make initial guess robust with respect to number of processors
      // by operating on the global index.
      types::global_dof_index first_local_range = 0;
      if (!vector.locally_owned_elements().is_empty())
        first_local_range = vector.locally_owned_elements().nth_index_in_set(0);

      const auto n_local_elements = vector.locally_owned_size();
      const int  n_blocks =
        1 + (n_local_elements - 1) / CUDAWrappers::block_size;
      set_initial_guess_kernel<<<n_blocks, CUDAWrappers::block_size>>>(
        first_local_range, n_local_elements, vector.get_values());
      AssertCudaKernel();

      const Number mean_value = vector.mean_value();
      vector.add(-mean_value);
    }
#  endif // DEAL_II_COMPILER_CUDA_AWARE

    struct EigenvalueTracker
    {
    public:
      void
      slot(const std::vector<double> &eigenvalues)
      {
        values = eigenvalues;
      }

      std::vector<double> values;
    };
  } // namespace PreconditionChebyshevImplementation
} // namespace internal



template <typename MatrixType, class VectorType, typename PreconditionerType>
inline PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  AdditionalData::AdditionalData(const unsigned int degree,
                                 const double       smoothing_range,
                                 const unsigned int eig_cg_n_iterations,
                                 const double       eig_cg_residual,
                                 const double       max_eigenvalue)
  : degree(degree)
  , smoothing_range(smoothing_range)
  , eig_cg_n_iterations(eig_cg_n_iterations)
  , eig_cg_residual(eig_cg_residual)
  , max_eigenvalue(max_eigenvalue)
{}



template <typename MatrixType, class VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::AdditionalData &
                  PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  AdditionalData::operator=(const AdditionalData &other_data)
{
  degree              = other_data.degree;
  smoothing_range     = other_data.smoothing_range;
  eig_cg_n_iterations = other_data.eig_cg_n_iterations;
  eig_cg_residual     = other_data.eig_cg_residual;
  max_eigenvalue      = other_data.max_eigenvalue;
  preconditioner      = other_data.preconditioner;
  constraints.copy_from(other_data.constraints);

  return *this;
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  PreconditionChebyshev()
  : theta(1.)
  , delta(1.)
  , eigenvalues_are_initialized(false)
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionChebyshev and VectorType must have the same size_type.");
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::initialize(
  const MatrixType &    matrix,
  const AdditionalData &additional_data)
{
  matrix_ptr = &matrix;
  data       = additional_data;
  Assert(data.degree > 0,
         ExcMessage("The degree of the Chebyshev method must be positive."));
  internal::PreconditionChebyshevImplementation::initialize_preconditioner(
    matrix, data.preconditioner);
  eigenvalues_are_initialized = false;
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::clear()
{
  eigenvalues_are_initialized = false;
  theta = delta = 1.0;
  matrix_ptr    = nullptr;
  {
    VectorType empty_vector;
    solution_old.reinit(empty_vector);
    temp_vector1.reinit(empty_vector);
    temp_vector2.reinit(empty_vector);
  }
  data.preconditioner.reset();
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::EigenvalueInformation
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  estimate_eigenvalues(const VectorType &src) const
{
  Assert(eigenvalues_are_initialized == false, ExcInternalError());
  Assert(data.preconditioner.get() != nullptr, ExcNotInitialized());

  PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
    EigenvalueInformation info{};

  solution_old.reinit(src);
  temp_vector1.reinit(src, true);

  if (data.eig_cg_n_iterations > 0)
    {
      Assert(data.eig_cg_n_iterations > 2,
             ExcMessage(
               "Need to set at least two iterations to find eigenvalues."));

      // set a very strict tolerance to force at least two iterations
      ReductionControl control(
        data.eig_cg_n_iterations,
        std::sqrt(
          std::numeric_limits<typename VectorType::value_type>::epsilon()),
        1e-10,
        false,
        false);

      internal::PreconditionChebyshevImplementation::EigenvalueTracker
                           eigenvalue_tracker;
      SolverCG<VectorType> solver(control);
      solver.connect_eigenvalues_slot(
        [&eigenvalue_tracker](const std::vector<double> &eigenvalues) {
          eigenvalue_tracker.slot(eigenvalues);
        });

      // set an initial guess that contains some high-frequency parts (to the
      // extent possible without knowing the discretization and the numbering)
      // to trigger high eigenvalues according to the external function
      internal::PreconditionChebyshevImplementation::set_initial_guess(
        temp_vector1);
      data.constraints.set_zero(temp_vector1);

      try
        {
          solver.solve(*matrix_ptr,
                       solution_old,
                       temp_vector1,
                       *data.preconditioner);
        }
      catch (SolverControl::NoConvergence &)
        {}

      // read the eigenvalues from the attached eigenvalue tracker
      if (eigenvalue_tracker.values.empty())
        info.min_eigenvalue_estimate = info.max_eigenvalue_estimate = 1.;
      else
        {
          info.min_eigenvalue_estimate = eigenvalue_tracker.values.front();

          // include a safety factor since the CG method will in general not
          // be converged
          info.max_eigenvalue_estimate = 1.2 * eigenvalue_tracker.values.back();
        }

      info.cg_iterations = control.last_step();
    }
  else
    {
      info.max_eigenvalue_estimate = data.max_eigenvalue;
      info.min_eigenvalue_estimate = data.max_eigenvalue / data.smoothing_range;
    }

  const double alpha = (data.smoothing_range > 1. ?
                          info.max_eigenvalue_estimate / data.smoothing_range :
                          std::min(0.9 * info.max_eigenvalue_estimate,
                                   info.min_eigenvalue_estimate));

  // in case the user set the degree to invalid unsigned int, we have to
  // determine the number of necessary iterations from the Chebyshev error
  // estimate, given the target tolerance specified by smoothing_range. This
  // estimate is based on the error formula given in section 5.1 of
  // R. S. Varga, Matrix iterative analysis, 2nd ed., Springer, 2009
  if (data.degree == numbers::invalid_unsigned_int)
    {
      const double actual_range = info.max_eigenvalue_estimate / alpha;
      const double sigma        = (1. - std::sqrt(1. / actual_range)) /
                           (1. + std::sqrt(1. / actual_range));
      const double eps = data.smoothing_range;
      const_cast<
        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(
        this)
        ->data.degree =
        1 + static_cast<unsigned int>(
              std::log(1. / eps + std::sqrt(1. / eps / eps - 1.)) /
              std::log(1. / sigma));
    }

  info.degree = data.degree;

  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(this)
    ->delta = (info.max_eigenvalue_estimate - alpha) * 0.5;
  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(this)
    ->theta = (info.max_eigenvalue_estimate + alpha) * 0.5;

  // We do not need the second temporary vector in case we have a
  // DiagonalMatrix as preconditioner and use deal.II's own vectors
  using NumberType = typename VectorType::value_type;
  if (std::is_same<PreconditionerType, DiagonalMatrix<VectorType>>::value ==
        false ||
      (std::is_same<VectorType, dealii::Vector<NumberType>>::value == false &&
       ((std::is_same<VectorType,
                      LinearAlgebra::distributed::
                        Vector<NumberType, MemorySpace::Host>>::value ==
         false) ||
        (std::is_same<VectorType,
                      LinearAlgebra::distributed::
                        Vector<NumberType, MemorySpace::CUDA>>::value ==
         false))))
    temp_vector2.reinit(src, true);
  else
    {
      VectorType empty_vector;
      temp_vector2.reinit(empty_vector);
    }

  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(this)
    ->eigenvalues_are_initialized = true;

  return info;
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::vmult(
  VectorType &      solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  internal::PreconditionChebyshevImplementation::vector_updates(
    rhs,
    *data.preconditioner,
    0,
    0.,
    1. / theta,
    solution_old,
    temp_vector1,
    temp_vector2,
    solution);

  // if delta is zero, we do not need to iterate because the updates will be
  // zero
  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      matrix_ptr->vmult(temp_vector1, solution);
      const double rhokp   = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      internal::PreconditionChebyshevImplementation::vector_updates(
        rhs,
        *data.preconditioner,
        k + 1,
        factor1,
        factor2,
        solution_old,
        temp_vector1,
        temp_vector2,
        solution);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::Tvmult(
  VectorType &      solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  internal::PreconditionChebyshevImplementation::vector_updates(
    rhs,
    *data.preconditioner,
    0,
    0.,
    1. / theta,
    solution_old,
    temp_vector1,
    temp_vector2,
    solution);

  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      matrix_ptr->Tvmult(temp_vector1, solution);
      const double rhokp   = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      internal::PreconditionChebyshevImplementation::vector_updates(
        rhs,
        *data.preconditioner,
        k + 1,
        factor1,
        factor2,
        solution_old,
        temp_vector1,
        temp_vector2,
        solution);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::step(
  VectorType &      solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  matrix_ptr->vmult(temp_vector1, solution);
  internal::PreconditionChebyshevImplementation::vector_updates(
    rhs,
    *data.preconditioner,
    1,
    0.,
    1. / theta,
    solution_old,
    temp_vector1,
    temp_vector2,
    solution);

  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      matrix_ptr->vmult(temp_vector1, solution);
      const double rhokp   = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      internal::PreconditionChebyshevImplementation::vector_updates(
        rhs,
        *data.preconditioner,
        k + 2,
        factor1,
        factor2,
        solution_old,
        temp_vector1,
        temp_vector2,
        solution);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::Tstep(
  VectorType &      solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  matrix_ptr->Tvmult(temp_vector1, solution);
  internal::PreconditionChebyshevImplementation::vector_updates(
    rhs,
    *data.preconditioner,
    1,
    0.,
    1. / theta,
    solution_old,
    temp_vector1,
    temp_vector2,
    solution);

  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      matrix_ptr->Tvmult(temp_vector1, solution);
      const double rhokp   = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      internal::PreconditionChebyshevImplementation::vector_updates(
        rhs,
        *data.preconditioner,
        k + 2,
        factor1,
        factor2,
        solution_old,
        temp_vector1,
        temp_vector2,
        solution);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::size_type
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::m() const
{
  Assert(matrix_ptr != nullptr, ExcNotInitialized());
  return matrix_ptr->m();
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::size_type
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::n() const
{
  Assert(matrix_ptr != nullptr, ExcNotInitialized());
  return matrix_ptr->n();
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


