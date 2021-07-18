//include/deal.II-translator/lac/solver_gmres_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_solver_gmres_h
#define dealii_solver_gmres_h



#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Solvers */ 
 /*@{*/ 

namespace internal
{
  /**
   * 一个用于GMRES求解器的辅助类的命名空间。
   *
   */
  namespace SolverGMRESImplementation
  {
    /**
     * 用于保存临时向量的类。
     * 一旦需要，这个类就会自动分配一个新的向量。
     * 未来的版本还应该能够自动转移向量，避免重新启动。
     *
     */

    template <typename VectorType>
    class TmpVectors
    {
    public:
      /**
       * 构造函数。准备一个长度为 @p  max_size的 @p VectorType
       * 数组。
       *
       */
      TmpVectors(const unsigned int max_size, VectorMemory<VectorType> &vmem);

      /**
       * 销毁器。删除所有分配的向量。
       *
       */
      ~TmpVectors() = default;

      /**
       * 获取向量编号  @p i.
       * 如果这个向量之前未被使用，会发生错误。
       *
       */
      VectorType &operator[](const unsigned int i) const;

      /**
       * 获取向量编号  @p i.  必要时分配它。
       * 如果必须分配一个向量，则使用 @p temp
       * 将其重新设定为适当的尺寸。
       *
       */
      VectorType &
      operator()(const unsigned int i, const VectorType &temp);

      /**
       * 返回数据向量的大小。它在求解器中用于存储Arnoldi向量。
       *
       */
      unsigned int
      size() const;


    private:
      /**
       * 获得向量的池子。
       *
       */
      VectorMemory<VectorType> &mem;

      /**
       * 用于存储向量的字段。
       *
       */
      std::vector<typename VectorMemory<VectorType>::Pointer> data;
    };
  } // namespace SolverGMRESImplementation
} // namespace internal

/**
 * 重启预设条件的直接广义最小残差法的实现。停止的标准是残差的规范。
 * AdditionalData结构包含使用的临时向量的数量。阿诺尔迪基的大小是这个数字减去3。此外，它允许你选择右预处理或左预处理。默认是左预处理。最后，它包括一个标志，表明是否使用默认残差作为停止标准。
 *
 *  <h3>Left versus right preconditioning</h3> @p AdditionalData
 * 允许你在左和右预处理之间进行选择。正如预期的那样，这将分别在对系统<i>P<sup>-1</sup>A</i>和<i>AP<sup>-1</sup></i>的求解之间切换。
 * 第二个结果是用来衡量收敛性的残差类型。对于左预处理，这是<b>preconditioned</b>的残差，而对于右预处理，它是未预处理系统的残差。
 * 可以选择使用标志 AdditionalData::use_default_residual.
 * 来覆盖这一行为。<tt>true</tt>值是指上一段描述的行为，而<tt>false</tt>则是恢复它。但是要注意，在这种情况下必须计算额外的残差，阻碍了求解器的整体性能。
 *
 *  <h3>The size of the Arnoldi basis</h3> 最大的基数大小是由
 * AdditionalData::max_n_tmp_vectors,
 * 控制的，它是这个数字减去2。如果迭代步数超过这个数字，所有的基向量都会被丢弃，并从目前得到的近似值重新开始迭代。
 * 请注意，GMRes的最小化特性只涉及到Arnoldi基所跨越的Krylov空间。因此，重新启动的GMRes不再是<b>not</b>的最小化。基长的选择是内存消耗和收敛速度之间的权衡，因为较长的基意味着在更大的空间内实现最小化。
 * 关于使用该类对矩阵和向量的要求，请参见Solver基类的文档。
 *
 *  <h3>Observing the progress of linear solver iterations</h3>
 * 这个类的solve()函数使用Solver基类中描述的机制来确定收敛性。这个机制也可以用来观察迭代的进度。
 *
 *  <h3>Eigenvalue and condition number estimates</h3>
 * 该类可以在求解过程中估计特征值和条件数。这是在内部迭代过程中通过创建海森堡矩阵来实现的。特征值被估计为海森堡矩阵的特征值，条件数被估计为海森堡矩阵的最大和最小奇异值的比率。估算值可以通过使用
 * @p  connect_condition_number_slot和 @p connect_eigenvalues_slot.
 * 连接一个函数作为槽来获得，然后这些槽将被解算器以估算值为参数调用。
 *
 *
 */
template <class VectorType = Vector<double>>
class SolverGMRES : public SolverBase<VectorType>
{
public:
  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。默认情况下，设置临时向量的数量为30，即每28次迭代做一次重启。同时从左边开始设置预设条件，停止准则的残差为默认残差，仅在必要时进行重新正交。
     *
     */
    explicit AdditionalData(const unsigned int max_n_tmp_vectors     = 30,
                            const bool         right_preconditioning = false,
                            const bool         use_default_residual  = true,
                            const bool force_re_orthogonalization    = false);

    /**
     * 临时向量的最大数量。这个参数控制Arnoldi基的大小，由于历史原因，这个参数是#max_n_tmp_vectors-2。SolverGMRES假设至少有三个临时向量，所以这个值必须大于或等于3。
     *
     */
    unsigned int max_n_tmp_vectors;

    /**
     * 右侧预处理的标志。
     * @note
     * 在左和右预处理之间的改变也会改变残差的评估方式。参见SolverGMRES中的相应章节。
     *
     */
    bool right_preconditioning;

    /**
     * 用于衡量收敛性的默认残差的标志。
     *
     */
    bool use_default_residual;

    /**
     * 在每一步中强制对正交基进行重新正交的标志。
     * 如果设置为false，求解器每5次迭代都会自动检查正交性的损失，只有在必要时才启用重新正交化。
     *
     */
    bool force_re_orthogonalization;
  };

  /**
   * 构造函数。
   *
   */
  SolverGMRES(SolverControl &           cn,
              VectorMemory<VectorType> &mem,
              const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverGMRES(SolverControl &cn, const AdditionalData &data = AdditionalData());

  /**
   * 复制构造函数被删除。
   *
   */
  SolverGMRES(const SolverGMRES<VectorType> &) = delete;

  /**
   * 求解x的线性系统 $Ax=b$ 。
   *
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

  /**
   * 连接一个槽来检索估计的条件数。如果every_iteration=true，则在每次外部迭代时调用，否则在迭代结束时调用一次（即因为已经达到收敛，或者因为检测到分歧）。
   *
   */
  boost::signals2::connection
  connect_condition_number_slot(const std::function<void(double)> &slot,
                                const bool every_iteration = false);

  /**
   * 连接一个槽来检索估计的特征值。如果every_iteration=true，则在每次外部迭代时调用，否则在迭代结束时调用一次（即，要么因为已经实现收敛，要么因为检测到分歧）。
   *
   */
  boost::signals2::connection
  connect_eigenvalues_slot(
    const std::function<void(const std::vector<std::complex<double>> &)> &slot,
    const bool every_iteration = false);

  /**
   * 连接一个槽，检索由初始矩阵在Krylov基上的投影得到的海森堡矩阵。如果every_iteration=true，则在每次外部迭代时调用，否则在迭代结束时调用一次（即因为已经达到收敛，或者因为检测到分歧）。
   *
   */
  boost::signals2::connection
  connect_hessenberg_slot(
    const std::function<void(const FullMatrix<double> &)> &slot,
    const bool every_iteration = true);

  /**
   * 连接一个槽，检索由Arnoldi算法生成的Krylov空间的基向量。当迭代完成时（即，因为已经达到收敛，或者因为检测到发散），一次性调用。
   *
   */
  boost::signals2::connection
  connect_krylov_space_slot(
    const std::function<
      void(const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)>
      &slot);


  /**
   * 连接一个槽，当向量被重新正交时，检索一个通知。
   *
   */
  boost::signals2::connection
  connect_re_orthogonalization_slot(const std::function<void(int)> &slot);


  DeclException1(ExcTooFewTmpVectors,
                 int,
                 << "The number of temporary vectors you gave (" << arg1
                 << ") is too small. It should be at least 10 for "
                 << "any results, and much more for reasonable ones.");

protected:
  /**
   * 包括tmp向量的最大数量。
   *
   */
  AdditionalData additional_data;

  /**
   * 用来检索估计条件数的信号。在所有迭代结束时调用一次。
   *
   */
  boost::signals2::signal<void(double)> condition_number_signal;

  /**
   * 用来检索估计条件数的信号。在每次外部迭代时被调用。
   *
   */
  boost::signals2::signal<void(double)> all_condition_numbers_signal;

  /**
   * 用来检索估计的特征值的信号。在所有迭代结束时调用一次。
   *
   */
  boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    eigenvalues_signal;

  /**
   * 用来检索估计的特征值的信号。在每次外部迭代时被调用。
   *
   */
  boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    all_eigenvalues_signal;

  /**
   * 用于检索海森堡矩阵的信号。在所有迭代结束时调用一次。
   *
   */
  boost::signals2::signal<void(const FullMatrix<double> &)> hessenberg_signal;

  /**
   * 用于检索海森堡矩阵的信号。在每次外部迭代时被调用。
   *
   */
  boost::signals2::signal<void(const FullMatrix<double> &)>
    all_hessenberg_signal;

  /**
   * 用来检索Krylov空间基向量的信号。当所有迭代结束时调用一次。
   *
   */
  boost::signals2::signal<void(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)>
    krylov_space_signal;

  /**
   * 当向量被重新正交化时，用于检索通知的信号。
   *
   */
  boost::signals2::signal<void(int)> re_orthogonalize_signal;

  /**
   * 实现对残差准则的计算。
   *
   */
  virtual double
  criterion();

  /**
   * 通过最后一列的基辅旋转，将上黑森伯格矩阵转化为三对角结构。
   *
   */
  void
  givens_rotation(Vector<double> &h,
                  Vector<double> &b,
                  Vector<double> &ci,
                  Vector<double> &si,
                  int             col) const;

  /**
   * 使用修改后的Gram-Schmidt算法将向量 @p vv
   * 与第一个参数给出的 @p dim （正交）向量正交。
   * 用于正交的因子存储在 @p h. 中。 布尔值 @p
   * re_orthogonalize指定是否应将修正的Gram-Schmidt算法应用两次。该算法每隔五步就会检查程序中的正交性损失，并在这种情况下将该标志设置为真。
   * 所有随后的迭代都使用重新正交。
   * 如果信号re_orthogonalize_signal被连接，则调用它。
   *
   */
  static double
  modified_gram_schmidt(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType>
      &                                       orthogonal_vectors,
    const unsigned int                        dim,
    const unsigned int                        accumulated_iterations,
    VectorType &                              vv,
    Vector<double> &                          h,
    bool &                                    re_orthogonalize,
    const boost::signals2::signal<void(int)> &re_orthogonalize_signal =
      boost::signals2::signal<void(int)>());

  /**
   * 从内部迭代过程中产生的海森堡矩阵H_orig中估计出特征值。使用这些估计值来计算条件数。以这些估计值为参数调用信号eigenvalues_signal和
   * cond_signal。
   *
   */
  static void
  compute_eigs_and_cond(
    const FullMatrix<double> &H_orig,
    const unsigned int        dim,
    const boost::signals2::signal<
      void(const std::vector<std::complex<double>> &)> &eigenvalues_signal,
    const boost::signals2::signal<void(const FullMatrix<double> &)>
      &                                          hessenberg_signal,
    const boost::signals2::signal<void(double)> &cond_signal);

  /**
   * 投射的系统矩阵
   *
   */
  FullMatrix<double> H;

  /**
   * 用于反转的辅助矩阵  @p H  。
   *
   */
  FullMatrix<double> H1;
};

/**
 * 实现具有灵活预处理的广义最小残差法（灵活的GMRES或FGMRES）。
 * 这个灵活版本的GMRES方法允许在每个迭代步骤中使用不同的预处理。因此，对于预处理器的不准确评估，它也更加稳健。一个重要的应用是在预处理程序中使用Krylov空间方法。与允许在左右预处理之间进行选择的SolverGMRES相反，这个求解器总是从右边应用预处理。
 * FGMRES在每个迭代步骤中需要两个向量，总共产生
 * <tt>2*SolverFGMRES::%AdditionalData::%max_basis_size+1</tt>
 * 个辅助向量。除此之外，与GMRES相比，FGMRES每次迭代所需的操作数量大致相同，只是在每次重启和solve()结束时少用一次预处理程序。
 * 更多细节见  @cite Saad1991  。
 *
 *
 */
template <class VectorType = Vector<double>>
class SolverFGMRES : public SolverBase<VectorType>
{
public:
  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   *
   */
  struct AdditionalData
  {
    /**
     * 构造函数。默认情况下，将最大基础尺寸设置为30。
     *
     */
    explicit AdditionalData(const unsigned int max_basis_size = 30)
      : max_basis_size(max_basis_size)
    {}

    /**
     * 最大基数大小。
     *
     */
    unsigned int max_basis_size;
  };

  /**
   * 构造函数。
   *
   */
  SolverFGMRES(SolverControl &           cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverFGMRES(SolverControl &       cn,
               const AdditionalData &data = AdditionalData());

  /**
   * 求解x的线性系统 $Ax=b$ 。
   *
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

private:
  /**
   * 额外的标志。
   *
   */
  AdditionalData additional_data;

  /**
   * 投射的系统矩阵
   *
   */
  FullMatrix<double> H;

  /**
   * 用于反转的辅助矩阵  @p H  。
   *
   */
  FullMatrix<double> H1;
};

 /*@}*/ 
 /* --------------------- Inline and template functions ------------------- */ 


#ifndef DOXYGEN
namespace internal
{
  namespace SolverGMRESImplementation
  {
    template <class VectorType>
    inline TmpVectors<VectorType>::TmpVectors(const unsigned int max_size,
                                              VectorMemory<VectorType> &vmem)
      : mem(vmem)
      , data(max_size)
    {}



    template <class VectorType>
    inline VectorType &TmpVectors<VectorType>::
                       operator[](const unsigned int i) const
    {
      AssertIndexRange(i, data.size());

      Assert(data[i] != nullptr, ExcNotInitialized());
      return *data[i];
    }



    template <class VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator()(const unsigned int i,
                                       const VectorType & temp)
    {
      AssertIndexRange(i, data.size());
      if (data[i] == nullptr)
        {
          data[i] = std::move(typename VectorMemory<VectorType>::Pointer(mem));
          data[i]->reinit(temp);
        }
      return *data[i];
    }



    template <class VectorType>
    unsigned int
    TmpVectors<VectorType>::size() const
    {
      return (data.size() > 0 ? data.size() - 1 : 0);
    }



    // A comparator for better printing eigenvalues
    inline bool
    complex_less_pred(const std::complex<double> &x,
                      const std::complex<double> &y)
    {
      return x.real() < y.real() ||
             (x.real() == y.real() && x.imag() < y.imag());
    }
  } // namespace SolverGMRESImplementation
} // namespace internal



template <class VectorType>
inline SolverGMRES<VectorType>::AdditionalData::AdditionalData(
  const unsigned int max_n_tmp_vectors,
  const bool         right_preconditioning,
  const bool         use_default_residual,
  const bool         force_re_orthogonalization)
  : max_n_tmp_vectors(max_n_tmp_vectors)
  , right_preconditioning(right_preconditioning)
  , use_default_residual(use_default_residual)
  , force_re_orthogonalization(force_re_orthogonalization)
{
  Assert(3 <= max_n_tmp_vectors,
         ExcMessage("SolverGMRES needs at least three "
                    "temporary vectors."));
}



template <class VectorType>
SolverGMRES<VectorType>::SolverGMRES(SolverControl &           cn,
                                     VectorMemory<VectorType> &mem,
                                     const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
SolverGMRES<VectorType>::SolverGMRES(SolverControl &       cn,
                                     const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <class VectorType>
inline void
SolverGMRES<VectorType>::givens_rotation(Vector<double> &h,
                                         Vector<double> &b,
                                         Vector<double> &ci,
                                         Vector<double> &si,
                                         int             col) const
{
  for (int i = 0; i < col; i++)
    {
      const double s     = si(i);
      const double c     = ci(i);
      const double dummy = h(i);
      h(i)               = c * dummy + s * h(i + 1);
      h(i + 1)           = -s * dummy + c * h(i + 1);
    };

  const double r = 1. / std::sqrt(h(col) * h(col) + h(col + 1) * h(col + 1));
  si(col)        = h(col + 1) * r;
  ci(col)        = h(col) * r;
  h(col)         = ci(col) * h(col) + si(col) * h(col + 1);
  b(col + 1)     = -si(col) * b(col);
  b(col) *= ci(col);
}



template <class VectorType>
inline double
SolverGMRES<VectorType>::modified_gram_schmidt(
  const internal::SolverGMRESImplementation::TmpVectors<VectorType>
    &                                       orthogonal_vectors,
  const unsigned int                        dim,
  const unsigned int                        accumulated_iterations,
  VectorType &                              vv,
  Vector<double> &                          h,
  bool &                                    reorthogonalize,
  const boost::signals2::signal<void(int)> &reorthogonalize_signal)
{
  Assert(dim > 0, ExcInternalError());
  const unsigned int inner_iteration = dim - 1;

  // need initial norm for detection of re-orthogonalization, see below
  double     norm_vv_start = 0;
  const bool consider_reorthogonalize =
    (reorthogonalize == false) && (inner_iteration % 5 == 4);
  if (consider_reorthogonalize)
    norm_vv_start = vv.l2_norm();

  // Orthogonalization
  h(0) = vv * orthogonal_vectors[0];
  for (unsigned int i = 1; i < dim; ++i)
    h(i) = vv.add_and_dot(-h(i - 1),
                          orthogonal_vectors[i - 1],
                          orthogonal_vectors[i]);
  double norm_vv =
    std::sqrt(vv.add_and_dot(-h(dim - 1), orthogonal_vectors[dim - 1], vv));

  // Re-orthogonalization if loss of orthogonality detected. For the test, use
  // a strategy discussed in C. T. Kelley, Iterative Methods for Linear and
  // Nonlinear Equations, SIAM, Philadelphia, 1995: Compare the norm of vv
  // after orthogonalization with its norm when starting the
  // orthogonalization. If vv became very small (here: less than the square
  // root of the machine precision times 10), it is almost in the span of the
  // previous vectors, which indicates loss of precision.
  if (consider_reorthogonalize)
    {
      if (norm_vv >
          10. * norm_vv_start *
            std::sqrt(
              std::numeric_limits<typename VectorType::value_type>::epsilon()))
        return norm_vv;

      else
        {
          reorthogonalize = true;
          if (!reorthogonalize_signal.empty())
            reorthogonalize_signal(accumulated_iterations);
        }
    }

  if (reorthogonalize == true)
    {
      double htmp = vv * orthogonal_vectors[0];
      h(0) += htmp;
      for (unsigned int i = 1; i < dim; ++i)
        {
          htmp = vv.add_and_dot(-htmp,
                                orthogonal_vectors[i - 1],
                                orthogonal_vectors[i]);
          h(i) += htmp;
        }
      norm_vv =
        std::sqrt(vv.add_and_dot(-htmp, orthogonal_vectors[dim - 1], vv));
    }

  return norm_vv;
}



template <class VectorType>
inline void
SolverGMRES<VectorType>::compute_eigs_and_cond(
  const FullMatrix<double> &H_orig,
  const unsigned int        dim,
  const boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    &eigenvalues_signal,
  const boost::signals2::signal<void(const FullMatrix<double> &)>
    &                                          hessenberg_signal,
  const boost::signals2::signal<void(double)> &cond_signal)
{
  // Avoid copying the Hessenberg matrix if it isn't needed.
  if ((!eigenvalues_signal.empty() || !hessenberg_signal.empty() ||
       !cond_signal.empty()) &&
      dim > 0)
    {
      LAPACKFullMatrix<double> mat(dim, dim);
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          mat(i, j) = H_orig(i, j);
      hessenberg_signal(H_orig);
      // Avoid computing eigenvalues if they are not needed.
      if (!eigenvalues_signal.empty())
        {
          // Copy mat so that we can compute svd below. Necessary since
          // compute_eigenvalues will leave mat in state
          // LAPACKSupport::unusable.
          LAPACKFullMatrix<double> mat_eig(mat);
          mat_eig.compute_eigenvalues();
          std::vector<std::complex<double>> eigenvalues(dim);
          for (unsigned int i = 0; i < mat_eig.n(); ++i)
            eigenvalues[i] = mat_eig.eigenvalue(i);
          // Sort eigenvalues for nicer output.
          std::sort(eigenvalues.begin(),
                    eigenvalues.end(),
                    internal::SolverGMRESImplementation::complex_less_pred);
          eigenvalues_signal(eigenvalues);
        }
      // Calculate condition number, avoid calculating the svd if a slot
      // isn't connected. Need at least a 2-by-2 matrix to do the estimate.
      if (!cond_signal.empty() && (mat.n() > 1))
        {
          mat.compute_svd();
          double condition_number =
            mat.singular_value(0) / mat.singular_value(mat.n() - 1);
          cond_signal(condition_number);
        }
    }
}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverGMRES<VectorType>::solve(const MatrixType &        A,
                               VectorType &              x,
                               const VectorType &        b,
                               const PreconditionerType &preconditioner)
{
  // TODO:[?] Check, why there are two different start residuals.
  // TODO:[GK] Make sure the parameter in the constructor means maximum basis
  // size

  LogStream::Prefix prefix("GMRES");

  // extra call to std::max to placate static analyzers: coverity rightfully
  // complains that data.max_n_tmp_vectors - 2 may overflow
  const unsigned int n_tmp_vectors =
    std::max(additional_data.max_n_tmp_vectors, 3u);

  // Generate an object where basis vectors are stored.
  internal::SolverGMRESImplementation::TmpVectors<VectorType> tmp_vectors(
    n_tmp_vectors, this->memory);

  // number of the present iteration; this
  // number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  const bool do_eigenvalues =
    !condition_number_signal.empty() || !all_condition_numbers_signal.empty() ||
    !eigenvalues_signal.empty() || !all_eigenvalues_signal.empty() ||
    !hessenberg_signal.empty() || !all_hessenberg_signal.empty();
  // for eigenvalue computation, need to collect the Hessenberg matrix (before
  // applying Givens rotations)
  FullMatrix<double> H_orig;
  if (do_eigenvalues)
    H_orig.reinit(n_tmp_vectors, n_tmp_vectors - 1);

  // matrix used for the orthogonalization process later
  H.reinit(n_tmp_vectors, n_tmp_vectors - 1);

  // some additional vectors, also used in the orthogonalization
  dealii::Vector<double> gamma(n_tmp_vectors), ci(n_tmp_vectors - 1),
    si(n_tmp_vectors - 1), h(n_tmp_vectors - 1);


  unsigned int dim = 0;

  SolverControl::State iteration_state = SolverControl::iterate;
  double               last_res        = -std::numeric_limits<double>::max();

  // switch to determine whether we want a left or a right preconditioner. at
  // present, left is default, but both ways are implemented
  const bool left_precondition = !additional_data.right_preconditioning;

  // Per default the left preconditioned GMRes uses the preconditioned
  // residual and the right preconditioned GMRes uses the unpreconditioned
  // residual as stopping criterion.
  const bool use_default_residual = additional_data.use_default_residual;

  // define two aliases
  VectorType &v = tmp_vectors(0, x);
  VectorType &p = tmp_vectors(n_tmp_vectors - 1, x);

  // Following vectors are needed when we are not using the default residuals
  // as stopping criterion
  typename VectorMemory<VectorType>::Pointer r;
  typename VectorMemory<VectorType>::Pointer x_;
  std::unique_ptr<dealii::Vector<double>>    gamma_;
  if (!use_default_residual)
    {
      r  = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      x_ = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      r->reinit(x);
      x_->reinit(x);

      gamma_ = std::make_unique<dealii::Vector<double>>(gamma.size());
    }

  bool re_orthogonalize = additional_data.force_re_orthogonalization;

  ///////////////////////////////////////////////////////////////////////////
  // outer iteration: loop until we either reach convergence or the maximum
  // number of iterations is exceeded. each cycle of this loop amounts to one
  // restart
  do
    {
      // reset this vector to the right size
      h.reinit(n_tmp_vectors - 1);

      if (left_precondition)
        {
          A.vmult(p, x);
          p.sadd(-1., 1., b);
          preconditioner.vmult(v, p);
        }
      else
        {
          A.vmult(v, x);
          v.sadd(-1., 1., b);
        };

      double rho = v.l2_norm();

      // check the residual here as well since it may be that we got the exact
      // (or an almost exact) solution vector at the outset. if we wouldn't
      // check here, the next scaling operation would produce garbage
      if (use_default_residual)
        {
          last_res = rho;
          iteration_state =
            this->iteration_status(accumulated_iterations, rho, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }
      else
        {
          deallog << "default_res=" << rho << std::endl;

          if (left_precondition)
            {
              A.vmult(*r, x);
              r->sadd(-1., 1., b);
            }
          else
            preconditioner.vmult(*r, v);

          double res = r->l2_norm();
          last_res   = res;
          iteration_state =
            this->iteration_status(accumulated_iterations, res, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }

      gamma(0) = rho;

      v *= 1. / rho;

      // inner iteration doing at most as many steps as there are temporary
      // vectors. the number of steps actually been done is propagated outside
      // through the @p dim variable
      for (unsigned int inner_iteration = 0;
           ((inner_iteration < n_tmp_vectors - 2) &&
            (iteration_state == SolverControl::iterate));
           ++inner_iteration)
        {
          ++accumulated_iterations;
          // yet another alias
          VectorType &vv = tmp_vectors(inner_iteration + 1, x);

          if (left_precondition)
            {
              A.vmult(p, tmp_vectors[inner_iteration]);
              preconditioner.vmult(vv, p);
            }
          else
            {
              preconditioner.vmult(p, tmp_vectors[inner_iteration]);
              A.vmult(vv, p);
            }

          dim = inner_iteration + 1;

          const double s         = modified_gram_schmidt(tmp_vectors,
                                                 dim,
                                                 accumulated_iterations,
                                                 vv,
                                                 h,
                                                 re_orthogonalize,
                                                 re_orthogonalize_signal);
          h(inner_iteration + 1) = s;

          // s=0 is a lucky breakdown, the solver will reach convergence,
          // but we must not divide by zero here.
          if (s != 0)
            vv *= 1. / s;

          // for eigenvalues, get the resulting coefficients from the
          // orthogonalization process
          if (do_eigenvalues)
            for (unsigned int i = 0; i < dim + 1; ++i)
              H_orig(i, inner_iteration) = h(i);

          //  Transformation into tridiagonal structure
          givens_rotation(h, gamma, ci, si, inner_iteration);

          //  append vector on matrix
          for (unsigned int i = 0; i < dim; ++i)
            H(i, inner_iteration) = h(i);

          //  default residual
          rho = std::fabs(gamma(dim));

          if (use_default_residual)
            {
              last_res = rho;
              iteration_state =
                this->iteration_status(accumulated_iterations, rho, x);
            }
          else
            {
              deallog << "default_res=" << rho << std::endl;

              dealii::Vector<double> h_(dim);
              *x_     = x;
              *gamma_ = gamma;
              H1.reinit(dim + 1, dim);

              for (unsigned int i = 0; i < dim + 1; ++i)
                for (unsigned int j = 0; j < dim; ++j)
                  H1(i, j) = H(i, j);

              H1.backward(h_, *gamma_);

              if (left_precondition)
                for (unsigned int i = 0; i < dim; ++i)
                  x_->add(h_(i), tmp_vectors[i]);
              else
                {
                  p = 0.;
                  for (unsigned int i = 0; i < dim; ++i)
                    p.add(h_(i), tmp_vectors[i]);
                  preconditioner.vmult(*r, p);
                  x_->add(1., *r);
                };
              A.vmult(*r, *x_);
              r->sadd(-1., 1., b);
              // Now *r contains the unpreconditioned residual!!
              if (left_precondition)
                {
                  const double res = r->l2_norm();
                  last_res         = res;

                  iteration_state =
                    this->iteration_status(accumulated_iterations, res, x);
                }
              else
                {
                  preconditioner.vmult(*x_, *r);
                  const double preconditioned_res = x_->l2_norm();
                  last_res                        = preconditioned_res;

                  iteration_state =
                    this->iteration_status(accumulated_iterations,
                                           preconditioned_res,
                                           x);
                }
            }
        };
      // end of inner iteration. now calculate the solution from the temporary
      // vectors
      h.reinit(dim);
      H1.reinit(dim + 1, dim);

      for (unsigned int i = 0; i < dim + 1; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          H1(i, j) = H(i, j);

      compute_eigs_and_cond(H_orig,
                            dim,
                            all_eigenvalues_signal,
                            all_hessenberg_signal,
                            condition_number_signal);

      H1.backward(h, gamma);

      if (left_precondition)
        for (unsigned int i = 0; i < dim; ++i)
          x.add(h(i), tmp_vectors[i]);
      else
        {
          p = 0.;
          for (unsigned int i = 0; i < dim; ++i)
            p.add(h(i), tmp_vectors[i]);
          preconditioner.vmult(v, p);
          x.add(1., v);
        };
      // end of outer iteration. restart if no convergence and the number of
      // iterations is not exceeded
    }
  while (iteration_state == SolverControl::iterate);

  compute_eigs_and_cond(H_orig,
                        dim,
                        eigenvalues_signal,
                        hessenberg_signal,
                        condition_number_signal);

  if (!krylov_space_signal.empty())
    krylov_space_signal(tmp_vectors);

  // in case of failure: throw exception
  AssertThrow(iteration_state == SolverControl::success,
              SolverControl::NoConvergence(accumulated_iterations, last_res));
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_condition_number_slot(
  const std::function<void(double)> &slot,
  const bool                         every_iteration)
{
  if (every_iteration)
    {
      return all_condition_numbers_signal.connect(slot);
    }
  else
    {
      return condition_number_signal.connect(slot);
    }
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_eigenvalues_slot(
  const std::function<void(const std::vector<std::complex<double>> &)> &slot,
  const bool every_iteration)
{
  if (every_iteration)
    {
      return all_eigenvalues_signal.connect(slot);
    }
  else
    {
      return eigenvalues_signal.connect(slot);
    }
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_hessenberg_slot(
  const std::function<void(const FullMatrix<double> &)> &slot,
  const bool                                             every_iteration)
{
  if (every_iteration)
    {
      return all_hessenberg_signal.connect(slot);
    }
  else
    {
      return hessenberg_signal.connect(slot);
    }
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_krylov_space_slot(
  const std::function<void(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)> &slot)
{
  return krylov_space_signal.connect(slot);
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_re_orthogonalization_slot(
  const std::function<void(int)> &slot)
{
  return re_orthogonalize_signal.connect(slot);
}



template <class VectorType>
double
SolverGMRES<VectorType>::criterion()
{
  // dummy implementation. this function is not needed for the present
  // implementation of gmres
  Assert(false, ExcInternalError());
  return 0;
}


//----------------------------------------------------------------------//

template <class VectorType>
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl &           cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl &       cn,
                                       const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverFGMRES<VectorType>::solve(const MatrixType &        A,
                                VectorType &              x,
                                const VectorType &        b,
                                const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("FGMRES");

  SolverControl::State iteration_state = SolverControl::iterate;

  const unsigned int basis_size = additional_data.max_basis_size;

  // Generate an object where basis vectors are stored.
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> v(
    basis_size, this->memory);
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> z(
    basis_size, this->memory);

  // number of the present iteration; this number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  // matrix used for the orthogonalization process later
  H.reinit(basis_size + 1, basis_size);

  // Vectors for projected system
  Vector<double> projected_rhs;
  Vector<double> y;

  // Iteration starts here
  double res = -std::numeric_limits<double>::max();

  typename VectorMemory<VectorType>::Pointer aux(this->memory);
  aux->reinit(x);
  do
    {
      A.vmult(*aux, x);
      aux->sadd(-1., 1., b);

      double beta     = aux->l2_norm();
      res             = beta;
      iteration_state = this->iteration_status(accumulated_iterations, res, x);
      if (iteration_state == SolverControl::success)
        break;

      H.reinit(basis_size + 1, basis_size);
      double a = beta;

      for (unsigned int j = 0; j < basis_size; ++j)
        {
          if (a != 0) // treat lucky breakdown
            v(j, x).equ(1. / a, *aux);
          else
            v(j, x) = 0.;


          preconditioner.vmult(z(j, x), v[j]);
          A.vmult(*aux, z[j]);

          // Gram-Schmidt
          H(0, j) = *aux * v[0];
          for (unsigned int i = 1; i <= j; ++i)
            H(i, j) = aux->add_and_dot(-H(i - 1, j), v[i - 1], v[i]);
          H(j + 1, j) = a = std::sqrt(aux->add_and_dot(-H(j, j), v[j], *aux));

          // Compute projected solution

          if (j > 0)
            {
              H1.reinit(j + 1, j);
              projected_rhs.reinit(j + 1);
              y.reinit(j);
              projected_rhs(0) = beta;
              H1.fill(H);

              // check convergence. note that the vector 'x' we pass to the
              // criterion is not the final solution we compute if we
              // decide to jump out of the iteration (we update 'x' again
              // right after the current loop)
              Householder<double> house(H1);
              res = house.least_squares(y, projected_rhs);
              iteration_state =
                this->iteration_status(++accumulated_iterations, res, x);
              if (iteration_state != SolverControl::iterate)
                break;
            }
        }

      // Update solution vector
      for (unsigned int j = 0; j < y.size(); ++j)
        x.add(y(j), z[j]);
    }
  while (iteration_state == SolverControl::iterate);

  // in case of failure: throw exception
  if (iteration_state != SolverControl::success)
    AssertThrow(false,
                SolverControl::NoConvergence(accumulated_iterations, res));
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


