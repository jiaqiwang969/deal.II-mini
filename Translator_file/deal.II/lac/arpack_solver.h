//include/deal.II-translator/lac/arpack_solver_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_arpack_solver_h
#define dealii_arpack_solver_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/solver_control.h>

#include <cstring>


#ifdef DEAL_II_WITH_ARPACK

DEAL_II_NAMESPACE_OPEN


extern "C" void
dnaupd_(int *         ido,
        char *        bmat,
        unsigned int *n,
        char *        which,
        unsigned int *nev,
        const double *tol,
        double *      resid,
        int *         ncv,
        double *      v,
        int *         ldv,
        int *         iparam,
        int *         ipntr,
        double *      workd,
        double *      workl,
        int *         lworkl,
        int *         info);

extern "C" void
dsaupd_(int *         ido,
        char *        bmat,
        unsigned int *n,
        char *        which,
        unsigned int *nev,
        double *      tol,
        double *      resid,
        int *         ncv,
        double *      v,
        int *         ldv,
        int *         iparam,
        int *         ipntr,
        double *      workd,
        double *      workl,
        int *         lworkl,
        int *         info);

extern "C" void
dneupd_(int *         rvec,
        char *        howmany,
        int *         select,
        double *      d,
        double *      di,
        double *      z,
        int *         ldz,
        double *      sigmar,
        double *      sigmai,
        double *      workev,
        char *        bmat,
        unsigned int *n,
        char *        which,
        unsigned int *nev,
        double *      tol,
        double *      resid,
        int *         ncv,
        double *      v,
        int *         ldv,
        int *         iparam,
        int *         ipntr,
        double *      workd,
        double *      workl,
        int *         lworkl,
        int *         info);

extern "C" void
dseupd_(int *         rvec,
        char *        howmany,
        int *         select,
        double *      d,
        double *      z,
        int *         ldz,
        double *      sigmar,
        char *        bmat,
        unsigned int *n,
        char *        which,
        unsigned int *nev,
        double *      tol,
        double *      resid,
        int *         ncv,
        double *      v,
        int *         ldv,
        int *         iparam,
        int *         ipntr,
        double *      workd,
        double *      workl,
        int *         lworkl,
        int *         info);

/**
 * 使用ARPACK的接口。ARPACK是一个Fortran77的子程序集，旨在解决大规模的特征值问题。
 * 这里我们为ARPACK的例程 <code>dnaupd</code> and <code>dneupd</code>
 * 提供接口。如果运算符被指定为对称的，我们就使用ARPACK的对称接口
 * <code>dsaupd</code> and <code>dseupd</code> 代替。
 * 该软件包被设计用来计算一般n乘n矩阵A的几个特征值和相应的特征向量。它最适用于大的稀疏矩阵A。
 * 在这个类中，我们利用了应用于广义特征谱问题 $(A-\lambda
 * B)x=0$ 的方法，用于 $x\neq0$ ；其中 $A$ 是一个系统矩阵，
 * $B$ 是一个质量矩阵，而 $\lambda, x$
 * 分别是一组特征值和特征向量。
 * ArpackSolver可以通过以下方式用于带有串行对象的应用代码中。
 *
 * @code
 * SolverControl solver_control(1000, 1e-9);
 * ArpackSolver solver(solver_control);
 * solver.solve(A, B, OP, lambda, x, size_of_spectrum);
 * @endcode
 * 对于广义的特征值问题 $Ax=B\lambda x$  ，其中变量
 * <code>size_of_spectrum</code>
 * 告诉ARPACK要解决的特征向量/特征值对的数量。这里，
 * <code>lambda</code> 是一个包含计算的特征值的向量，
 * <code>x</code> 是一个包含计算的特征向量的向量，
 * <code>OP</code> 是对矩阵 <code>A</code>
 * 的逆运算。围绕零点的移位和反转变换被应用。
 * 通过AdditionalData，用户可以指定一些要设置的参数。
 * 关于ARPACK例程 <code>dsaupd</code> ,  <code>dseupd</code>,
 * <code>dnaupd</code> and <code>dneupd</code>
 * 如何工作以及如何适当设置参数的进一步信息，请看ARPACK手册。
 *
 *
 * @note  每当你使用AffineConstraints消除自由度时，你会产生虚假的特征值和特征向量。如果你确保消除的矩阵行的对角线都等于1，你会得到一个额外的特征值。但要注意deal.II中的一些函数将这些对角线设置为相当任意的（从特征值问题的角度来看）值。参见 @ref step_36  "  step-36  "
 * 的例子。
 *
 *
 */
class ArpackSolver : public Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;


  /**
   * 一个枚举，列出在solve()函数中计算哪些特征值的可能选择。
   *
   */
  enum WhichEigenvalues
  {
    /**
     * 代数上最大的特征值。
     *
     */
    algebraically_largest,
    /**
     * 代数上最小的特征值。
     *
     */
    algebraically_smallest,
    /**
     * 具有最大量级的特征值。
     *
     */
    largest_magnitude,
    /**
     * 具有最小量级的特征值。
     *
     */
    smallest_magnitude,
    /**
     * 具有最大实部的特征值。
     *
     */
    largest_real_part,
    /**
     * 具有最小实部的特征值。
     *
     */
    smallest_real_part,
    /**
     * 具有最大虚部的特征值。
     *
     */
    largest_imaginary_part,
    /**
     * 具有最小虚部的特征值。
     *
     */
    smallest_imaginary_part,
    /**
     * 从频谱的高端计算一半的特征值，另一半从低端计算。如果要求的特征向量的数量是奇数，那么额外的特征向量来自谱的高端。
     *
     */
    both_ends
  };

  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。默认情况下，将Arnoldi向量的数量（如果问题是对称的，则为Lanczos向量）设置为15。对于非对称问题，设置求解器寻找最大的特征值）。)
     *
     */
    explicit AdditionalData(
      const unsigned int     number_of_arnoldi_vectors = 15,
      const WhichEigenvalues eigenvalue_of_interest    = largest_magnitude,
      const bool             symmetric                 = false);

    /**
     * Arnoldi/Lanczos向量的数量。这个数字应该小于问题的大小，但要大于2倍的特征值数量（如果设置了n_eigenvalues）加1。
     *
     */
    const unsigned int number_of_arnoldi_vectors;

    /**
     * 指定感兴趣的特征值。
     *
     */
    const WhichEigenvalues eigenvalue_of_interest;

    /**
     * 指定问题是否是对称的。
     *
     */
    const bool symmetric;
  };

  /**
   * 访问控制收敛的对象。
   *
   */
  SolverControl &
  control() const;

  /**
   * 构造函数。
   *
   */
  ArpackSolver(SolverControl &       control,
               const AdditionalData &data = AdditionalData());

  /**
   * 设置初始矢量，用于构建Krylov空间。
   *
   */
  template <typename VectorType>
  void
  set_initial_vector(const VectorType &vec);

  /**
   * 设置移位 @p sigma ，用于移位和反转的光谱变换。
   * 如果这个函数没有被调用，则假定移位为零。
   *
   */
  void
  set_shift(const std::complex<double> sigma);

  /**
   * 通过调用ARPACK的 <code>dsaupd</code> and <code>dseupd</code> 或
   * <code>dnaupd</code> and <code>dneupd</code>
   * 函数来解决广义的特征直角问题 $A x=\lambda B x$ 。
   * 该函数返回一个长度为<i>n</i>的特征值向量和一个长度为<i>n</i>的特征向量（在对称情况下），以及长度为<i>n+1</i>的非对称情况下的向量。在对称情况下，所有的特征向量都是实数。在非对称情况下，复杂的特征值总是以复数共轭对的形式出现。因此，具有非零复数部分的特征值的特征向量是通过将实数和虚数部分放在连续的实值向量中来存储的。复数共轭特征值的特征向量不需要被存储，因为它只是存储的特征向量的复数共轭。因此，如果最后n个特征值有一个非零的虚部，Arpack总共需要n+1个实值向量来存储特征向量的实部和虚部。
   * @param  A
   * 我们要计算特征值的算子。实际上，这个参数完全没有使用。
   * @param  B
   * 基础空间的内积，通常是质量矩阵。对于受限问题，它可以是一个部分质量矩阵，例如，像斯托克斯问题的速度质量矩阵。只有它的函数
   * <code>vmult()</code> 被使用。      @param  inverse
   * 这是实际使用的可能移位的逆矩阵，而不是 <code>A</code>.
   * Only its function <code>vmult()</code> 被使用。      @param
   * eigenvalues 是一个复数的向量，其中的特征值被返回。
   * @param
   * 特征向量是一个<b>real</b>特征向量的向量，包含所有特征向量的实部和与复数共轭特征值对对应的特征向量的虚部。
   * 因此，其长度在对称情况下应为<i>n</i>，在非对称情况下应为<i>n+1</i>。在非对称情况下，存储方案会导致例如以下模式。假设前两个特征值是实数，第三和第四是一对复数共轭。询问三个特征对的结果是<i>[real(v1),real(v2),
   * real(v3),imag(v3)]</i>。请注意，如果我们在这个例子中要求四个特征对，我们会得到同样的模式，因为第四个特征向量只是第三个的复共轭。
   * @param  n_eigenvalues
   * 这个参数的目的并不清楚，但将其设置为
   * <code>eigenvalues</code> 或更大的规模是安全的。
   * 让它保持默认的0，它将在内部被重置为
   * <code>eigenvalues</code> 的大小。
   *
   */
  template <typename VectorType,
            typename MatrixType1,
            typename MatrixType2,
            typename INVERSE>
  void
  solve(const MatrixType1 &                A,
        const MatrixType2 &                B,
        const INVERSE &                    inverse,
        std::vector<std::complex<double>> &eigenvalues,
        std::vector<VectorType> &          eigenvectors,
        const unsigned int                 n_eigenvalues = 0);

protected:
  /**
   * 对控制迭代求解器收敛的对象的引用。
   *
   */
  SolverControl &solver_control;

  /**
   * 存储这个特定求解器的标志的副本。
   *
   */
  const AdditionalData additional_data;

  /**
   * 存储一个初始向量
   *
   */
  bool                initial_vector_provided;
  std::vector<double> resid;

  /**
   * 移位的实数部分
   *
   */
  double sigmar;

  /**
   * 移位的虚数部分
   *
   */
  double sigmai;


private:
  /**
   * 例外的情况。
   *
   */
  DeclException2(ArpackExcInvalidNumberofEigenvalues,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(ArpackExcInvalidEigenvectorSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvectors " << arg2);

  DeclException2(
    ArpackExcInvalidEigenvectorSizeNonsymmetric,
    int,
    int,
    << "To store the real and complex parts of " << arg1
    << " eigenvectors in real-valued vectors, their size (currently set to "
    << arg2 << ") should be greater than or equal to " << arg1 + 1);

  DeclException2(ArpackExcInvalidEigenvalueSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvalues " << arg2);

  DeclException2(ArpackExcInvalidNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(ArpackExcSmallNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is too small to obtain " << arg2 << " eigenvalues");

  DeclException1(ArpackExcArpackIdo,
                 int,
                 << "This ido " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(ArpackExcArpackMode,
                 int,
                 << "This mode " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodsaupd,
                 int,
                 << "Error with dsaupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodnaupd,
                 int,
                 << "Error with dnaupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodseupd,
                 int,
                 << "Error with dseupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodneupd,
                 int,
                 << "Error with dneupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfoMaxIt,
                 int,
                 << "Maximum number " << arg1 << " of iterations reached.");

  DeclExceptionMsg(ArpackExcArpackNoShifts,
                   "No shifts could be applied during implicit"
                   " Arnoldi update, try increasing the number of"
                   " Arnoldi vectors.");
};


inline ArpackSolver::AdditionalData::AdditionalData(
  const unsigned int     number_of_arnoldi_vectors,
  const WhichEigenvalues eigenvalue_of_interest,
  const bool             symmetric)
  : number_of_arnoldi_vectors(number_of_arnoldi_vectors)
  , eigenvalue_of_interest(eigenvalue_of_interest)
  , symmetric(symmetric)
{
  // Check for possible options for symmetric problems
  if (symmetric)
    {
      Assert(
        eigenvalue_of_interest != largest_real_part,
        ExcMessage(
          "'largest real part' can only be used for non-symmetric problems!"));
      Assert(
        eigenvalue_of_interest != smallest_real_part,
        ExcMessage(
          "'smallest real part' can only be used for non-symmetric problems!"));
      Assert(
        eigenvalue_of_interest != largest_imaginary_part,
        ExcMessage(
          "'largest imaginary part' can only be used for non-symmetric problems!"));
      Assert(
        eigenvalue_of_interest != smallest_imaginary_part,
        ExcMessage(
          "'smallest imaginary part' can only be used for non-symmetric problems!"));
    }
  // Check for possible options for asymmetric problems
  else
    {
      Assert(
        eigenvalue_of_interest != algebraically_largest,
        ExcMessage(
          "'largest algebraic part' can only be used for symmetric problems!"));
      Assert(
        eigenvalue_of_interest != algebraically_smallest,
        ExcMessage(
          "'smallest algebraic part' can only be used for symmetric problems!"));
      Assert(eigenvalue_of_interest != both_ends,
             ExcMessage(
               "'both ends' can only be used for symmetric problems!"));
    }
}


inline ArpackSolver::ArpackSolver(SolverControl &       control,
                                  const AdditionalData &data)
  : solver_control(control)
  , additional_data(data)
  , initial_vector_provided(false)
  , sigmar(0.0)
  , sigmai(0.0)
{}



inline void
ArpackSolver::set_shift(const std::complex<double> sigma)
{
  sigmar = sigma.real();
  sigmai = sigma.imag();
}



template <typename VectorType>
inline void
ArpackSolver::set_initial_vector(const VectorType &vec)
{
  initial_vector_provided = true;
  resid.resize(vec.size());
  for (size_type i = 0; i < vec.size(); ++i)
    resid[i] = vec[i];
}


template <typename VectorType,
          typename MatrixType1,
          typename MatrixType2,
          typename INVERSE>
inline void
ArpackSolver::solve(const MatrixType1 &  /*system_matrix*/ ,
                    const MatrixType2 &                mass_matrix,
                    const INVERSE &                    inverse,
                    std::vector<std::complex<double>> &eigenvalues,
                    std::vector<VectorType> &          eigenvectors,
                    const unsigned int                 n_eigenvalues)
{
  // Problem size
  unsigned int n = eigenvectors[0].size();

  // Number of eigenvalues
  const unsigned int nev_const =
    (n_eigenvalues == 0) ? eigenvalues.size() : n_eigenvalues;
  // nev for arpack, which might change by plus one during dneupd
  unsigned int nev = nev_const;

  // check input sizes
  if (additional_data.symmetric)
    {
      Assert(nev <= eigenvectors.size(),
             ArpackExcInvalidEigenvectorSize(nev, eigenvectors.size()));
    }
  else
    Assert(nev + 1 <= eigenvectors.size(),
           ArpackExcInvalidEigenvectorSizeNonsymmetric(nev,
                                                       eigenvectors.size()));

  Assert(nev <= eigenvalues.size(),
         ArpackExcInvalidEigenvalueSize(nev, eigenvalues.size()));

  // check large enough problem size
  Assert(nev < n, ArpackExcInvalidNumberofEigenvalues(nev, n));

  Assert(additional_data.number_of_arnoldi_vectors < n,
         ArpackExcInvalidNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, n));

  // check whether we have enough Arnoldi vectors
  Assert(additional_data.number_of_arnoldi_vectors > 2 * nev + 1,
         ArpackExcSmallNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, nev));

  // ARPACK mode for dsaupd/dnaupd, here only mode 3, i.e. shift-invert mode
  int mode = 3;

  // reverse communication parameter
  int ido = 0;

  // 'G' generalized eigenvalue problem 'I' standard eigenvalue problem
  char bmat[2] = "G";

  // Specify the eigenvalues of interest, possible parameters "LA" algebraically
  // largest "SA" algebraically smallest "LM" largest magnitude "SM" smallest
  // magnitude "LR" largest real part "SR" smallest real part "LI" largest
  // imaginary part "SI" smallest imaginary part "BE" both ends of spectrum
  // simultaneous.
  char which[3];
  switch (additional_data.eigenvalue_of_interest)
    {
      case algebraically_largest:
        std::strcpy(which, "LA");
        break;
      case algebraically_smallest:
        std::strcpy(which, "SA");
        break;
      case largest_magnitude:
        std::strcpy(which, "LM");
        break;
      case smallest_magnitude:
        std::strcpy(which, "SM");
        break;
      case largest_real_part:
        std::strcpy(which, "LR");
        break;
      case smallest_real_part:
        std::strcpy(which, "SR");
        break;
      case largest_imaginary_part:
        std::strcpy(which, "LI");
        break;
      case smallest_imaginary_part:
        std::strcpy(which, "SI");
        break;
      case both_ends:
        std::strcpy(which, "BE");
        break;
    }

  // tolerance for ARPACK
  double tol = control().tolerance();

  // if the starting vector is used it has to be in resid
  if (!initial_vector_provided || resid.size() != n)
    resid.resize(n, 1.);

  // number of Arnoldi basis vectors specified
  // in additional_data
  int ncv = additional_data.number_of_arnoldi_vectors;

  int                 ldv = n;
  std::vector<double> v(ldv * ncv, 0.0);

  // information to the routines
  std::vector<int> iparam(11, 0);

  iparam[0] = 1; // shift strategy

  // maximum number of iterations
  iparam[2] = control().max_steps();

  // Set the mode of dsaupd. 1 is exact shifting, 2 is user-supplied shifts,
  // 3 is shift-invert mode, 4 is buckling mode, 5 is Cayley mode.

  iparam[6] = mode;
  std::vector<int> ipntr(14, 0);

  // work arrays for ARPACK
  std::vector<double> workd(3 * n, 0.);
  int                 lworkl =
    additional_data.symmetric ? ncv * ncv + 8 * ncv : 3 * ncv * ncv + 6 * ncv;
  std::vector<double> workl(lworkl, 0.);

  // information out of the iteration
  int info = 1;

  while (ido != 99)
    {
      // call of ARPACK dsaupd/dnaupd routine
      if (additional_data.symmetric)
        dsaupd_(&ido,
                bmat,
                &n,
                which,
                &nev,
                &tol,
                resid.data(),
                &ncv,
                v.data(),
                &ldv,
                iparam.data(),
                ipntr.data(),
                workd.data(),
                workl.data(),
                &lworkl,
                &info);
      else
        dnaupd_(&ido,
                bmat,
                &n,
                which,
                &nev,
                &tol,
                resid.data(),
                &ncv,
                v.data(),
                &ldv,
                iparam.data(),
                ipntr.data(),
                workd.data(),
                workl.data(),
                &lworkl,
                &info);

      if (ido == 99)
        break;

      switch (mode)
        {
          case 3:
            {
              switch (ido)
                {
                  case -1:
                    {
                      VectorType src, dst, tmp;
                      src.reinit(eigenvectors[0]);
                      dst.reinit(src);
                      tmp.reinit(src);


                      for (size_type i = 0; i < src.size(); ++i)
                        src(i) = workd[ipntr[0] - 1 + i];

                      // multiplication with mass matrix M
                      mass_matrix.vmult(tmp, src);
                      // solving linear system
                      inverse.vmult(dst, tmp);

                      for (size_type i = 0; i < dst.size(); ++i)
                        workd[ipntr[1] - 1 + i] = dst(i);
                    }
                    break;

                  case 1:
                    {
                      VectorType src, dst, tmp, tmp2;
                      src.reinit(eigenvectors[0]);
                      dst.reinit(src);
                      tmp.reinit(src);
                      tmp2.reinit(src);

                      for (size_type i = 0; i < src.size(); ++i)
                        {
                          src(i) = workd[ipntr[2] - 1 + i];
                          tmp(i) = workd[ipntr[0] - 1 + i];
                        }
                      // solving linear system
                      inverse.vmult(dst, src);

                      for (size_type i = 0; i < dst.size(); ++i)
                        workd[ipntr[1] - 1 + i] = dst(i);
                    }
                    break;

                  case 2:
                    {
                      VectorType src, dst;
                      src.reinit(eigenvectors[0]);
                      dst.reinit(src);

                      for (size_type i = 0; i < src.size(); ++i)
                        src(i) = workd[ipntr[0] - 1 + i];

                      // Multiplication with mass matrix M
                      mass_matrix.vmult(dst, src);

                      for (size_type i = 0; i < dst.size(); ++i)
                        workd[ipntr[1] - 1 + i] = dst(i);
                    }
                    break;

                  default:
                    Assert(false, ArpackExcArpackIdo(ido));
                    break;
                }
            }
            break;
          default:
            Assert(false, ArpackExcArpackMode(mode));
            break;
        }
    }

  // Set number of used iterations in SolverControl
  control().check(iparam[2], 0.);

  if (info < 0)
    {
      if (additional_data.symmetric)
        {
          Assert(false, ArpackExcArpackInfodsaupd(info));
        }
      else
        Assert(false, ArpackExcArpackInfodnaupd(info));
    }
  else
    {
      // 1 - compute eigenvectors, 0 - only eigenvalues
      int rvec = 1;

      // which eigenvectors
      char howmany = 'A';

      std::vector<int> select(ncv, 1);

      int ldz = n;

      std::vector<double> eigenvalues_real(nev + 1, 0.);
      std::vector<double> eigenvalues_im(nev + 1, 0.);

      // call of ARPACK dseupd/dneupd routine
      if (additional_data.symmetric)
        {
          std::vector<double> z(ldz * nev, 0.);
          dseupd_(&rvec,
                  &howmany,
                  select.data(),
                  eigenvalues_real.data(),
                  z.data(),
                  &ldz,
                  &sigmar,
                  bmat,
                  &n,
                  which,
                  &nev,
                  &tol,
                  resid.data(),
                  &ncv,
                  v.data(),
                  &ldv,
                  iparam.data(),
                  ipntr.data(),
                  workd.data(),
                  workl.data(),
                  &lworkl,
                  &info);
        }
      else
        {
          std::vector<double> workev(3 * ncv, 0.);
          dneupd_(&rvec,
                  &howmany,
                  select.data(),
                  eigenvalues_real.data(),
                  eigenvalues_im.data(),
                  v.data(),
                  &ldz,
                  &sigmar,
                  &sigmai,
                  workev.data(),
                  bmat,
                  &n,
                  which,
                  &nev,
                  &tol,
                  resid.data(),
                  &ncv,
                  v.data(),
                  &ldv,
                  iparam.data(),
                  ipntr.data(),
                  workd.data(),
                  workl.data(),
                  &lworkl,
                  &info);
        }

      if (info == 1)
        {
          Assert(false, ArpackExcArpackInfoMaxIt(control().max_steps()));
        }
      else if (info == 3)
        {
          Assert(false, ArpackExcArpackNoShifts());
        }
      else if (info != 0)
        {
          if (additional_data.symmetric)
            {
              Assert(false, ArpackExcArpackInfodseupd(info));
            }
          else
            Assert(false, ArpackExcArpackInfodneupd(info));
        }

      for (unsigned int i = 0; i < nev; ++i)
        for (unsigned int j = 0; j < n; ++j)
          eigenvectors[i](j) = v[i * n + j];

      for (unsigned int i = 0; i < nev_const; ++i)
        eigenvalues[i] =
          std::complex<double>(eigenvalues_real[i], eigenvalues_im[i]);
    }
}


inline SolverControl &
ArpackSolver::control() const
{
  return solver_control;
}

DEAL_II_NAMESPACE_CLOSE


#endif
#endif


