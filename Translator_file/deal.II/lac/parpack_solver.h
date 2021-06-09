//include/deal.II-translator/lac/parpack_solver_0.txt
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

#ifndef dealii_parpack_solver_h
#define dealii_parpack_solver_h

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_operation.h>

#include <cstring>


#ifdef DEAL_II_ARPACK_WITH_PARPACK

DEAL_II_NAMESPACE_OPEN

extern "C"
{
  // http://www.mathkeisan.com/usersguide/man/pdnaupd.html
  void
  pdnaupd_(MPI_Fint *comm,
           int *     ido,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // http://www.mathkeisan.com/usersguide/man/pdsaupd.html
  void
  pdsaupd_(MPI_Fint *comm,
           int *     ido,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // http://www.mathkeisan.com/usersguide/man/pdneupd.html
  void
  pdneupd_(MPI_Fint *comm,
           int *     rvec,
           char *    howmany,
           int *     select,
           double *  d,
           double *  di,
           double *  z,
           int *     ldz,
           double *  sigmar,
           double *  sigmai,
           double *  workev,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // http://www.mathkeisan.com/usersguide/man/pdseupd.html
  void
  pdseupd_(MPI_Fint *comm,
           int *     rvec,
           char *    howmany,
           int *     select,
           double *  d,
           double *  z,
           int *     ldz,
           double *  sigmar,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // other resources:
  //    http://acts.nersc.gov/superlu/example5/pnslac.c.html
  //    https://github.com/phpisciuneri/tijo/blob/master/dvr_parpack.cpp
}

/**
 * 使用PARPACK的接口。PARPACK是一个Fortran77子程序的集合，旨在解决大规模的特征值问题。这里我们为PARPACK的例程
 * <code>pdneupd</code>, <code>pdseupd</code> 、 <code>pdnaupd</code>,
 * <code>pdsaupd</code> 提供接口。
 * 该软件包被设计用来计算一般n乘n矩阵A的几个特征值和相应的特征向量。它最适用于大的稀疏矩阵A。
 * 在这个类中，我们利用了应用于广义特征谱问题 $(A-\lambda
 * B)x=0$ 的方法，用于 $x\neq0$ ；其中 $A$ 是一个系统矩阵，
 * $B$ 是一个质量矩阵，而 $\lambda, x$
 * 分别是一组特征值和特征向量。
 * ArpackSolver可以通过以下方式在应用代码中使用。
 *
 * @code
 * SolverControl solver_control (1000, 1e-9);
 * const unsigned int num_arnoldi_vectors = 2*size_of_spectrum + 2;
 * PArpackSolver<V>::AdditionalData
 *   additional_data(num_arnoldi_vectors,
 *                   dealii::PArpackSolver<V>::largest_magnitude,
 *                   true);
 *
 *  PArpackSolver<V> eigensolver (solver_control,
 *                                mpi_communicator,
 *                                additional_data);
 *  eigensolver.set_shift(sigma);
 *  eigensolver.reinit(locally_owned_dofs);
 *  eigensolver.solve (A,
 *                     B,
 *                     OP,
 *                     lambda,
 *                     x,
 *                     size_of_spectrum);
 * @endcode
 * 对于广义的特征值问题 $Ax=B\lambda x$  ，其中变量
 * <code>size_of_spectrum</code>
 * 告诉PARPACK要解决的特征向量/特征值对的数量。这里，
 * <code>lambda</code> 是一个包含计算的特征值的向量，
 * <code>x</code> 是一个包含计算的特征向量的 <code>V</code>
 * 类型的对象的向量。
 * 目前，只有(P)Arpack的三种模式被实现。在模式3（默认）中，
 * <code>OP</code> 是对矩阵<code>A的逆运算。
 *
 * - sigma B</code>，其中 <code> sigma </code> 是一个移位值，默认设置为零。而在模式2中， <code>OP</code> 是 <code>M</code> 的逆运算。最后，模式1对应于没有谱系转换的标准特征值问题  $Ax=\lambda x$  。模式可以通过AdditionalData对象指定。请注意，对于移位和反转（模式=3），所寻求的特征对是应用谱系转换后的特征对。
 * <code>OP</code>  可以通过使用LinearOperator来指定。
 *
 * @code
 * const double shift = 5.0;
 * const auto op_A = linear_operator<vector_t>(A);
 * const auto op_B = linear_operator<vector_t>(B);
 * const auto op_shift = op_A
 *
 * - shift op_B;
 * SolverControl solver_control_lin (1000, 1e-10,false,false);
 *
 * SolverCG<vector_t> cg(solver_control_lin);
 * const auto op_shift_invert =
 *   inverse_operator(op_shift, cg, PreconditionIdentity ());
 * @endcode
 *
 * 该类旨在与MPI一起使用，可以在任意的向量和矩阵分布式类上工作。
 * 对称和非对称的 <code>A</code> 都支持。 关于PARPACK例程
 * <code>pdneupd</code> 、 <code>pdseupd</code>, <code>pdnaupd</code>,
 * <code>pdsaupd</code>
 * 如何工作以及如何适当设置参数的进一步信息，请查看PARPACK手册。
 *
 *
 */
template <typename VectorType>
class PArpackSolver : public Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 一个枚举，列出在solve()函数中计算哪些特征值的可能选择。注意，这与应用移位和反转（目前唯一支持的谱系变换）后的问题相对应。
   * 一个特定的选择是基于对称或非对称矩阵 <code>A</code>
   * 考虑的限制。
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
   * 标准化的数据结构，如果需要的话，可以将额外的数据输送到求解器中。
   *
   */
  struct AdditionalData
  {
    const unsigned int     number_of_arnoldi_vectors;
    const WhichEigenvalues eigenvalue_of_interest;
    const bool             symmetric;
    const int              mode;
    AdditionalData(
      const unsigned int     number_of_arnoldi_vectors = 15,
      const WhichEigenvalues eigenvalue_of_interest    = largest_magnitude,
      const bool             symmetric                 = false,
      const int              mode                      = 3);
  };

  /**
   * 对控制收敛的对象的访问。
   *
   */
  SolverControl &
  control() const;

  /**
   * 构造函数。
   *
   */
  PArpackSolver(SolverControl &       control,
                const MPI_Comm &      mpi_communicator,
                const AdditionalData &data = AdditionalData());

  /**
   * 初始化内部变量。
   *
   */
  void
  reinit(const IndexSet &locally_owned_dofs);

  /**
   * 在使用BlockVectors时初始化内部变量。    @p locally_owned_dofs
   * 用于设置问题的维度，而 @p partitioning
   * 用于调用所使用的deal.II blockvector的reinit。
   *
   */
  void
  reinit(const IndexSet &             locally_owned_dofs,
         const std::vector<IndexSet> &partitioning);

  /**
   * 从输入中初始化内部变量  @p distributed_vector.  。
   *
   */
  void
  reinit(const VectorType &distributed_vector);

  /**
   * 设置初始矢量，用于建立Krylov空间。
   *
   */
  void
  set_initial_vector(const VectorType &vec);

  /**
   * 设置移位 @p sigma ，用于移位和反转的光谱变换。
   * 如果这个函数没有被调用，则假定移位为零。
   * @note 只与 <code>mode=3</code>
   * 有关（关于什么是不同模式的定义，请看这个类的一般文档）。
   *
   */
  void
  set_shift(const std::complex<double> sigma);

  /**
   * 通过调用PARPACK的 <code>pd(n/s)eupd</code> and
   * <code>pd(n/s)aupd</code> 函数来解决广义的特征直角问题 $A
   * x=\lambda B x$ 。    在 <code>mode=3</code> 中， @p inverse
   * 应对应于 $[A-\sigma B]^{-1}$ ，而在 <code>mode=2</code>
   * 中应代表 $B^{-1}$ 。对于 <code>mode=1</code> ， @p B 和 @p
   * inverse 都被忽略了。
   *
   */
  template <typename MatrixType1, typename MatrixType2, typename INVERSE>
  void
  solve(const MatrixType1 &                A,
        const MatrixType2 &                B,
        const INVERSE &                    inverse,
        std::vector<std::complex<double>> &eigenvalues,
        std::vector<VectorType> &          eigenvectors,
        const unsigned int                 n_eigenvalues);

  /**
   * 与上述相同，但将特征向量作为指针。
   *
   */
  template <typename MatrixType1, typename MatrixType2, typename INVERSE>
  void
  solve(const MatrixType1 &                A,
        const MatrixType2 &                B,
        const INVERSE &                    inverse,
        std::vector<std::complex<double>> &eigenvalues,
        std::vector<VectorType *> &        eigenvectors,
        const unsigned int                 n_eigenvalues);

  /**
   * 以字节为单位返回该类的内存消耗。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 对控制迭代求解器收敛性的对象的引用。
   *
   */
  SolverControl &solver_control;

  /**
   * 存储这个特定求解器的标志的副本。
   *
   */
  const AdditionalData additional_data;

  // keep MPI communicator non-const as Arpack functions are not const either:

  /**
   * C++ MPI通信器。
   *
   */
  MPI_Comm mpi_communicator;

  /**
   * Fortran MPI通信器。
   *
   */
  MPI_Fint mpi_communicator_fortran;

  // C++98 guarantees that the elements of a vector are stored contiguously

  /**
   * 工作数组的长度workl。
   *
   */
  int lworkl;

  /**
   * 长度为lworkl的双精度工作数组
   *
   */
  std::vector<double> workl;

  /**
   * 长度为3*N的双精度工作数组
   *
   */
  std::vector<double> workd;

  /**
   * 本地自由度的数量。
   *
   */
  int nloc;

  /**
   * Additional_data中指定的Arnoldi基向量的数量
   *
   */
  int ncv;


  /**
   * 数组v的前导维度
   *
   */
  int ldv;

  /**
   * 双精度向量，大小为ldv的NCV。
   * 将包含最终的Arnoldi基向量集。
   *
   */
  std::vector<double> v;

  /**
   * 一个辅助标志，当提供初始向量时被设置为真。
   *
   */
  bool initial_vector_provided;

  /**
   * 初始残差向量，可能来自先前的运行。
   * 在输出时，它包含最终的残差向量。
   *
   */
  std::vector<double> resid;

  /**
   * 数组Z的前导尺寸等于nloc。
   *
   */
  int ldz;

  /**
   * 一个最小尺寸为nloc的NEV+1的向量。 Z包含了特征系统A*z =
   * lambda*B*z的B-正态里兹向量，对应于里兹值的近似值。
   *
   */
  std::vector<double> z;

  /**
   * Workev数组的大小。
   *
   */
  int lworkev;

  /**
   * 尺寸为3*NCV的双精度工作阵列。
   *
   */
  std::vector<double> workev;

  /**
   * 维度为NCV的向量。
   *
   */
  std::vector<int> select;

  /**
   * 在Arpack和deal.II之间使用的临时向量
   *
   */
  VectorType src, dst, tmp;

  /**
   * 局部自由度的索引。
   *
   */
  std::vector<types::global_dof_index> local_indices;

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
   * 初始化依赖于 @p locally_owned_dofs. 的内部变量
   * 该函数在reinit()函数中被调用。
   *
   */
  void
  internal_reinit(const IndexSet &locally_owned_dofs);

  /**
   * PArpackExcInfoPdnaupds.
   *
   */
  DeclException2(PArpackExcConvergedEigenvectors,
                 int,
                 int,
                 << arg1 << " eigenpairs were requested, but only " << arg2
                 << " converged");

  DeclException2(PArpackExcInvalidNumberofEigenvalues,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(PArpackExcInvalidEigenvectorSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvectors " << arg2);

  DeclException2(
    PArpackExcInvalidEigenvectorSizeNonsymmetric,
    int,
    int,
    << "To store the real and complex parts of " << arg1
    << " eigenvectors in real-valued vectors, their size (currently set to "
    << arg2 << ") should be greater than or equal to " << arg1 + 1);

  DeclException2(PArpackExcInvalidEigenvalueSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvalues " << arg2);

  DeclException2(PArpackExcInvalidNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(PArpackExcSmallNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is too small to obtain " << arg2 << " eigenvalues");

  DeclException1(PArpackExcIdo,
                 int,
                 << "This ido " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(PArpackExcMode,
                 int,
                 << "This mode " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(PArpackExcInfoPdnaupd,
                 int,
                 << "Error with Pdnaupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(PArpackExcInfoPdneupd,
                 int,
                 << "Error with Pdneupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(PArpackExcInfoMaxIt,
                 int,
                 << "Maximum number " << arg1 << " of iterations reached.");

  DeclException1(PArpackExcNoShifts,
                 int,
                 << "No shifts could be applied during implicit"
                 << " Arnoldi update, try increasing the number of"
                 << " Arnoldi vectors.");
};



template <typename VectorType>
std::size_t
PArpackSolver<VectorType>::memory_consumption() const
{
  return MemoryConsumption::memory_consumption(double()) *
           (workl.size() + workd.size() + v.size() + resid.size() + z.size() +
            workev.size()) +
         src.memory_consumption() + dst.memory_consumption() +
         tmp.memory_consumption() +
         MemoryConsumption::memory_consumption(types::global_dof_index()) *
           local_indices.size();
}



template <typename VectorType>
PArpackSolver<VectorType>::AdditionalData::AdditionalData(
  const unsigned int     number_of_arnoldi_vectors,
  const WhichEigenvalues eigenvalue_of_interest,
  const bool             symmetric,
  const int              mode)
  : number_of_arnoldi_vectors(number_of_arnoldi_vectors)
  , eigenvalue_of_interest(eigenvalue_of_interest)
  , symmetric(symmetric)
  , mode(mode)
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
  Assert(mode >= 1 && mode <= 3,
         ExcMessage("Currently, only modes 1, 2 and 3 are supported."));
}



template <typename VectorType>
PArpackSolver<VectorType>::PArpackSolver(SolverControl &       control,
                                         const MPI_Comm &      mpi_communicator,
                                         const AdditionalData &data)
  : solver_control(control)
  , additional_data(data)
  , mpi_communicator(mpi_communicator)
  , mpi_communicator_fortran(MPI_Comm_c2f(mpi_communicator))
  , lworkl(0)
  , nloc(0)
  , ncv(0)
  , ldv(0)
  , initial_vector_provided(false)
  , ldz(0)
  , lworkev(0)
  , sigmar(0.0)
  , sigmai(0.0)
{}



template <typename VectorType>
void
PArpackSolver<VectorType>::set_shift(const std::complex<double> sigma)
{
  sigmar = sigma.real();
  sigmai = sigma.imag();
}



template <typename VectorType>
void
PArpackSolver<VectorType>::set_initial_vector(const VectorType &vec)
{
  initial_vector_provided = true;
  Assert(resid.size() == local_indices.size(),
         ExcDimensionMismatch(resid.size(), local_indices.size()));
  vec.extract_subvector_to(local_indices.begin(),
                           local_indices.end(),
                           resid.data());
}



template <typename VectorType>
void
PArpackSolver<VectorType>::internal_reinit(const IndexSet &locally_owned_dofs)
{
  // store local indices to write to vectors
  locally_owned_dofs.fill_index_vector(local_indices);

  // scalars
  nloc = locally_owned_dofs.n_elements();
  ncv  = additional_data.number_of_arnoldi_vectors;

  AssertDimension(local_indices.size(), nloc);

  // vectors
  ldv = nloc;
  v.resize(ldv * ncv, 0.0);

  resid.resize(nloc, 1.0);

  // work arrays for ARPACK
  workd.resize(3 * nloc, 0.0);

  lworkl =
    additional_data.symmetric ? ncv * ncv + 8 * ncv : 3 * ncv * ncv + 6 * ncv;
  workl.resize(lworkl, 0.);

  ldz = nloc;
  z.resize(ldz * ncv, 0.); // TODO we actually need only ldz*nev

  // WORKEV  Double precision  work array of dimension 3*NCV.
  lworkev = additional_data.symmetric ? 0  /*not used in symmetric case*/ 
                                        :
                                        3 * ncv;
  workev.resize(lworkev, 0.);

  select.resize(ncv, 0);
}



template <typename VectorType>
void
PArpackSolver<VectorType>::reinit(const IndexSet &locally_owned_dofs)
{
  internal_reinit(locally_owned_dofs);

  // deal.II vectors:
  src.reinit(locally_owned_dofs, mpi_communicator);
  dst.reinit(locally_owned_dofs, mpi_communicator);
  tmp.reinit(locally_owned_dofs, mpi_communicator);
}



template <typename VectorType>
void
PArpackSolver<VectorType>::reinit(const VectorType &distributed_vector)
{
  internal_reinit(distributed_vector.locally_owned_elements());

  // deal.II vectors:
  src.reinit(distributed_vector);
  dst.reinit(distributed_vector);
  tmp.reinit(distributed_vector);
}



template <typename VectorType>
void
PArpackSolver<VectorType>::reinit(const IndexSet &locally_owned_dofs,
                                  const std::vector<IndexSet> &partitioning)
{
  internal_reinit(locally_owned_dofs);

  // deal.II vectors:
  src.reinit(partitioning, mpi_communicator);
  dst.reinit(partitioning, mpi_communicator);
  tmp.reinit(partitioning, mpi_communicator);
}



template <typename VectorType>
template <typename MatrixType1, typename MatrixType2, typename INVERSE>
void
PArpackSolver<VectorType>::solve(const MatrixType1 &                A,
                                 const MatrixType2 &                B,
                                 const INVERSE &                    inverse,
                                 std::vector<std::complex<double>> &eigenvalues,
                                 std::vector<VectorType> &eigenvectors,
                                 const unsigned int       n_eigenvalues)
{
  std::vector<VectorType *> eigenvectors_ptr(eigenvectors.size());
  for (unsigned int i = 0; i < eigenvectors.size(); ++i)
    eigenvectors_ptr[i] = &eigenvectors[i];
  solve(A, B, inverse, eigenvalues, eigenvectors_ptr, n_eigenvalues);
}



template <typename VectorType>
template <typename MatrixType1, typename MatrixType2, typename INVERSE>
void
PArpackSolver<VectorType>::solve(const MatrixType1 &system_matrix,
                                 const MatrixType2 &mass_matrix,
                                 const INVERSE &    inverse,
                                 std::vector<std::complex<double>> &eigenvalues,
                                 std::vector<VectorType *> &eigenvectors,
                                 const unsigned int         n_eigenvalues)
{
  if (additional_data.symmetric)
    {
      Assert(n_eigenvalues <= eigenvectors.size(),
             PArpackExcInvalidEigenvectorSize(n_eigenvalues,
                                              eigenvectors.size()));
    }
  else
    Assert(n_eigenvalues + 1 <= eigenvectors.size(),
           PArpackExcInvalidEigenvectorSizeNonsymmetric(n_eigenvalues,
                                                        eigenvectors.size()));

  Assert(n_eigenvalues <= eigenvalues.size(),
         PArpackExcInvalidEigenvalueSize(n_eigenvalues, eigenvalues.size()));


  // use eigenvectors to get the problem size so that it is possible to
  // employ LinearOperator for mass_matrix.
  Assert(n_eigenvalues < eigenvectors[0]->size(),
         PArpackExcInvalidNumberofEigenvalues(n_eigenvalues,
                                              eigenvectors[0]->size()));

  Assert(additional_data.number_of_arnoldi_vectors < eigenvectors[0]->size(),
         PArpackExcInvalidNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, eigenvectors[0]->size()));

  Assert(additional_data.number_of_arnoldi_vectors > 2 * n_eigenvalues + 1,
         PArpackExcSmallNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, n_eigenvalues));

  int mode = additional_data.mode;

  // reverse communication parameter
  // must be zero on the first call to pdnaupd
  int ido = 0;

  // 'G' generalized eigenvalue problem
  // 'I' standard eigenvalue problem
  char bmat[2];
  bmat[0] = (mode == 1) ? 'I' : 'G';
  bmat[1] = '\0';

  // Specify the eigenvalues of interest, possible parameters:
  // "LA" algebraically largest
  // "SA" algebraically smallest
  // "LM" largest magnitude
  // "SM" smallest magnitude
  // "LR" largest real part
  // "SR" smallest real part
  // "LI" largest imaginary part
  // "SI" smallest imaginary part
  // "BE" both ends of spectrum simultaneous
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

  // information to the routines
  std::vector<int> iparam(11, 0);

  iparam[0] = 1;
  // shift strategy: exact shifts with respect to the current Hessenberg matrix
  // H.

  // maximum number of iterations
  iparam[2] = control().max_steps();

  // Parpack currently works only for NB = 1
  iparam[3] = 1;

  // Sets the mode of dsaupd:
  // 1 is A*x=lambda*x, OP = A, B = I
  // 2 is A*x = lambda*M*x, OP = inv[M]*A, B = M
  // 3 is shift-invert mode, OP = inv[A-sigma*M]*M, B = M
  // 4 is buckling mode,
  // 5 is Cayley mode.

  iparam[6] = mode;
  std::vector<int> ipntr(14, 0);

  // information out of the iteration
  //  If INFO .EQ. 0, a random initial residual vector is used.
  //  If INFO .NE. 0, RESID contains the initial residual vector,
  //  possibly from a previous run.
  // Typical choices in this situation might be to use the final value
  // of the starting vector from the previous eigenvalue calculation
  int info = initial_vector_provided ? 1 : 0;

  // Number of eigenvalues of OP to be computed. 0 < NEV < N.
  int nev             = n_eigenvalues;
  int n_inside_arpack = nloc;

  // IDO = 99: done
  while (ido != 99)
    {
      // call of ARPACK pdnaupd routine
      if (additional_data.symmetric)
        pdsaupd_(&mpi_communicator_fortran,
                 &ido,
                 bmat,
                 &n_inside_arpack,
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
        pdnaupd_(&mpi_communicator_fortran,
                 &ido,
                 bmat,
                 &n_inside_arpack,
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

      AssertThrow(info == 0, PArpackExcInfoPdnaupd(info));

      // if we converge, we shall not modify anything in work arrays!
      if (ido == 99)
        break;

      // IPNTR(1) is the pointer into WORKD for X,
      // IPNTR(2) is the pointer into WORKD for Y.
      const int shift_x = ipntr[0] - 1;
      const int shift_y = ipntr[1] - 1;
      Assert(shift_x >= 0, dealii::ExcInternalError());
      Assert(shift_x + nloc <= static_cast<int>(workd.size()),
             dealii::ExcInternalError());
      Assert(shift_y >= 0, dealii::ExcInternalError());
      Assert(shift_y + nloc <= static_cast<int>(workd.size()),
             dealii::ExcInternalError());

      src = 0.;

      // switch based on both ido and mode
      if ((ido == -1) || (ido == 1 && mode < 3))
        // compute  Y = OP * X
        {
          src.add(nloc, local_indices.data(), workd.data() + shift_x);
          src.compress(VectorOperation::add);

          if (mode == 3)
            // OP = inv[K - sigma*M]*M
            {
              mass_matrix.vmult(tmp, src);
              inverse.vmult(dst, tmp);
            }
          else if (mode == 2)
            // OP = inv[M]*K
            {
              system_matrix.vmult(tmp, src);
              // store M*X in X
              tmp.extract_subvector_to(local_indices.begin(),
                                       local_indices.end(),
                                       workd.data() + shift_x);
              inverse.vmult(dst, tmp);
            }
          else if (mode == 1)
            {
              system_matrix.vmult(dst, src);
            }
          else
            AssertThrow(false, PArpackExcMode(mode));
        }
      else if (ido == 1 && mode >= 3)
        // compute  Y = OP * X for mode 3, 4 and 5, where
        // the vector B * X is already available in WORKD(ipntr(3)).
        {
          const int shift_b_x = ipntr[2] - 1;
          Assert(shift_b_x >= 0, dealii::ExcInternalError());
          Assert(shift_b_x + nloc <= static_cast<int>(workd.size()),
                 dealii::ExcInternalError());

          // B*X
          src.add(nloc, local_indices.data(), workd.data() + shift_b_x);
          src.compress(VectorOperation::add);

          // solving linear system
          Assert(mode == 3, ExcNotImplemented());
          inverse.vmult(dst, src);
        }
      else if (ido == 2)
        // compute  Y = B * X
        {
          src.add(nloc, local_indices.data(), workd.data() + shift_x);
          src.compress(VectorOperation::add);

          // Multiplication with mass matrix M
          if (mode == 1)
            {
              dst = src;
            }
          else
            // mode 2,3 and 5 have B=M
            {
              mass_matrix.vmult(dst, src);
            }
        }
      else
        AssertThrow(false, PArpackExcIdo(ido));
      // Note: IDO = 3 does not appear to be required for currently
      // implemented modes

      // store the result
      dst.extract_subvector_to(local_indices.begin(),
                               local_indices.end(),
                               workd.data() + shift_y);
    } // end of pd*aupd_ loop

  // 1 - compute eigenvectors,
  // 0 - only eigenvalues
  int rvec = 1;

  // which eigenvectors
  char howmany[4] = "All";

  std::vector<double> eigenvalues_real(n_eigenvalues + 1, 0.);
  std::vector<double> eigenvalues_im(n_eigenvalues + 1, 0.);

  // call of ARPACK pdneupd routine
  if (additional_data.symmetric)
    pdseupd_(&mpi_communicator_fortran,
             &rvec,
             howmany,
             select.data(),
             eigenvalues_real.data(),
             z.data(),
             &ldz,
             &sigmar,
             bmat,
             &n_inside_arpack,
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
    pdneupd_(&mpi_communicator_fortran,
             &rvec,
             howmany,
             select.data(),
             eigenvalues_real.data(),
             eigenvalues_im.data(),
             v.data(),
             &ldz,
             &sigmar,
             &sigmai,
             workev.data(),
             bmat,
             &n_inside_arpack,
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

  if (info == 1)
    {
      AssertThrow(false, PArpackExcInfoMaxIt(control().max_steps()));
    }
  else if (info == 3)
    {
      AssertThrow(false, PArpackExcNoShifts(1));
    }
  else if (info != 0)
    {
      AssertThrow(false, PArpackExcInfoPdneupd(info));
    }

  for (int i = 0; i < nev; ++i)
    {
      (*eigenvectors[i]) = 0.0;
      AssertIndexRange(i * nloc + nloc, v.size() + 1);

      eigenvectors[i]->add(nloc, local_indices.data(), &v[i * nloc]);
      eigenvectors[i]->compress(VectorOperation::add);
    }

  for (size_type i = 0; i < n_eigenvalues; ++i)
    eigenvalues[i] =
      std::complex<double>(eigenvalues_real[i], eigenvalues_im[i]);

  // Throw an error if the solver did not converge.
  AssertThrow(iparam[4] >= static_cast<int>(n_eigenvalues),
              PArpackExcConvergedEigenvectors(n_eigenvalues, iparam[4]));

  // both PDNAUPD and PDSAUPD compute eigenpairs of inv[A - sigma*M]*M
  // with respect to a semi-inner product defined by M.

  // resid likely contains residual with respect to M-norm.
  {
    tmp = 0.0;
    tmp.add(nloc, local_indices.data(), resid.data());
    solver_control.check(iparam[2], tmp.l2_norm());
  }
}



template <typename VectorType>
SolverControl &
PArpackSolver<VectorType>::control() const
{
  return solver_control;
}

DEAL_II_NAMESPACE_CLOSE


#endif
#endif


