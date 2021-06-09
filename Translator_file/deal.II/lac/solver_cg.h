//include/deal.II-translator/lac/solver_cg_0.txt
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

#ifndef dealii_solver_cg_h
#define dealii_solver_cg_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/tridiagonal_matrix.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#ifndef DOXYGEN
class PreconditionIdentity;
#endif


 /*!@addtogroup Solvers */ 
 /*@{*/ 

/**
 * 该类实现了预设条件的共轭梯度（CG）方法，可用于解决具有对称正定矩阵的线性系统。该类首先用于
 * step-3 和 step-4
 * ，但也用于其他许多教程程序。像所有其他求解器类一样，只要满足一定的要求，它可以在任何类型的向量和矩阵上工作（关于使用该类时对矩阵和向量的要求，见求解器基类的文档）。求解向量的类型必须作为模板参数传递，默认为
 * dealii::Vector<double>.  。
 *
 *
 * @note  这个版本的CG取自D.
 * Braess的《有限元》一书。它需要一个对称的预处理程序（也就是说，例如，SOR不是一个可能的选择）。
 *
 *  <h3>Eigenvalue computation</h3>
 * cg方法对原始预处理的线性系统进行正交投影，使其成为另一个维度较小的系统。此外，投影矩阵
 * @p T 是三对角的。由于投影是正交的， @p T
 * 的特征值与原始预处理矩阵 @p PA.
 * 的特征值相近。事实上，在 @p n 步之后，其中 @p n
 * 是原始系统的维度，两个矩阵的特征值相等。但是，即使是小数目的迭代步骤，
 * @p T 的条件数也是对 @p PA. 的条件数的良好估计。 在 @p m
 * 步之后，矩阵T_m可以写成系数 @p alpha 和 @p beta  ]
 * 为三对角矩阵，其对角线元素为<tt>1/alpha_0</tt>, <tt>1/alpha_1
 * + beta_0/alpha_0</tt>, ...
 * ，<tt>1/alpha_{m-1</tt>+beta_{m-2}/alpha_{m-2}}}和非对角线元素<tt>sqrt(beta_0)/alpha_0</tt>,
 * ...,
 * <tt>sqrt(beta_{m-2</tt>)/alpha_{m-2}。这个矩阵的特征值可以通过后处理计算出来。
 * @see  Y. Saad: "Iterative methods for Sparse Linear Systems", section
 * 6.7.3了解详情。
 * 系数、特征值和条件数（计算为最大特征值与最小特征值之比）可以通过使用函数
 * @p  connect_coefficients_slot,  @p connect_eigenvalues_slot  和  @p
 * connect_condition_number_slot之一将一个函数作为槽连接到求解器中来获得。然后，这些槽将从求解器中被调用，并以估计值为参数。
 * <h3>Observing the progress of linear solver iterations</h3>
 * 这个类的solve()函数使用Solver基类中描述的机制来确定收敛性。这个机制也可以用来观察迭代的进度。
 *
 *
 */
template <typename VectorType = Vector<double>>
class SolverCG : public SolverBase<VectorType>
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 标准化的数据结构，用于向求解器输送额外的数据。
   * 这里，它不存储任何东西，只是为了与其他求解器类保持一致而存在。
   *
   */
  struct AdditionalData
  {};

  /**
   * 构造函数。
   *
   */
  SolverCG(SolverControl &           cn,
           VectorMemory<VectorType> &mem,
           const AdditionalData &    data = AdditionalData());

  /**
   * 构造函数。使用一个GrowingVectorMemory类型的对象作为默认分配内存。
   *
   */
  SolverCG(SolverControl &cn, const AdditionalData &data = AdditionalData());

  /**
   * 虚拟解构器。
   *
   */
  virtual ~SolverCG() override = default;

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
   * 连接一个槽来检索CG系数。该槽将以α作为第一个参数，以β作为第二个参数被调用，其中α和β遵循Y.
   * Saad: "Iterative methods for Sparse Linear Systems", section 6.7.
   * 每次迭代调用一次
   *
   */
  boost::signals2::connection
  connect_coefficients_slot(
    const std::function<void(typename VectorType::value_type,
                             typename VectorType::value_type)> &slot);

  /**
   * 连接一个槽来检索估计的条件数。如果every_iteration=true，则在每次迭代时调用，否则在迭代结束时调用一次（即，因为已经达到收敛，或者因为检测到发散）。
   *
   */
  boost::signals2::connection
  connect_condition_number_slot(const std::function<void(double)> &slot,
                                const bool every_iteration = false);

  /**
   * 连接一个槽来检索估计的特征值。如果every_iteration=true，则在每次迭代时调用，否则在迭代结束时调用一次（即，要么因为已经达到收敛，要么因为检测到分歧）。
   *
   */
  boost::signals2::connection
  connect_eigenvalues_slot(
    const std::function<void(const std::vector<double> &)> &slot,
    const bool every_iteration = false);

protected:
  /**
   * 派生类的接口。这个函数在每个步骤中获得当前迭代向量、残差和更新向量。它可以用于收敛历史的图形输出。
   *
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

  /**
   * 估计对角线和非对角线的特征值。使用这些估计值来计算条件数。以这些估计值为参数调用信号eigenvalues_signal和
   * cond_signal。
   *
   */
  static void
  compute_eigs_and_cond(
    const std::vector<typename VectorType::value_type> &diagonal,
    const std::vector<typename VectorType::value_type> &offdiagonal,
    const boost::signals2::signal<void(const std::vector<double> &)>
      &                                          eigenvalues_signal,
    const boost::signals2::signal<void(double)> &cond_signal);

  /**
   * 附加参数。
   *
   */
  AdditionalData additional_data;

  /**
   * 用于检索CG系数的信号。在每次迭代时调用。
   *
   */
  boost::signals2::signal<void(typename VectorType::value_type,
                               typename VectorType::value_type)>
    coefficients_signal;

  /**
   * 用来检索估计条件数的信号。当所有迭代结束时调用一次。
   *
   */
  boost::signals2::signal<void(double)> condition_number_signal;

  /**
   * 用于检索估计条件数的信号。在每次迭代时被调用。
   *
   */
  boost::signals2::signal<void(double)> all_condition_numbers_signal;

  /**
   * 用于检索估计的特征值的信号。在所有迭代结束时调用一次。
   *
   */
  boost::signals2::signal<void(const std::vector<double> &)> eigenvalues_signal;

  /**
   * 用于检索估计的特征值的信号。在每次迭代时被调用。
   *
   */
  boost::signals2::signal<void(const std::vector<double> &)>
    all_eigenvalues_signal;
};

 /*@}*/ 

 /*------------------------- Implementation ----------------------------*/ 

#ifndef DOXYGEN

template <typename VectorType>
SolverCG<VectorType>::SolverCG(SolverControl &           cn,
                               VectorMemory<VectorType> &mem,
                               const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <typename VectorType>
SolverCG<VectorType>::SolverCG(SolverControl &cn, const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <typename VectorType>
void
SolverCG<VectorType>::print_vectors(const unsigned int,
                                    const VectorType &,
                                    const VectorType &,
                                    const VectorType &) const
{}



template <typename VectorType>
inline void
SolverCG<VectorType>::compute_eigs_and_cond(
  const std::vector<typename VectorType::value_type> &diagonal,
  const std::vector<typename VectorType::value_type> &offdiagonal,
  const boost::signals2::signal<void(const std::vector<double> &)>
    &                                          eigenvalues_signal,
  const boost::signals2::signal<void(double)> &cond_signal)
{
  // Avoid computing eigenvalues unless they are needed.
  if (!cond_signal.empty() || !eigenvalues_signal.empty())
    {
      TridiagonalMatrix<typename VectorType::value_type> T(diagonal.size(),
                                                           true);
      for (size_type i = 0; i < diagonal.size(); ++i)
        {
          T(i, i) = diagonal[i];
          if (i < diagonal.size() - 1)
            T(i, i + 1) = offdiagonal[i];
        }
      T.compute_eigenvalues();
      // Need two eigenvalues to estimate the condition number.
      if (diagonal.size() > 1)
        {
          auto condition_number = T.eigenvalue(T.n() - 1) / T.eigenvalue(0);
          // Condition number is real valued and nonnegative; simply take
          // the absolute value:
          cond_signal(std::abs(condition_number));
        }
      // Avoid copying the eigenvalues of T to a vector unless a signal is
      // connected.
      if (!eigenvalues_signal.empty())
        {
          std::vector<double> eigenvalues(T.n());
          for (unsigned int j = 0; j < T.n(); ++j)
            {
              // for a hermitian matrix, all eigenvalues are real-valued
              // and non-negative, simply return the absolute value:
              eigenvalues[j] = std::abs(T.eigenvalue(j));
            }
          eigenvalues_signal(eigenvalues);
        }
    }
}



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverCG<VectorType>::solve(const MatrixType &        A,
                            VectorType &              x,
                            const VectorType &        b,
                            const PreconditionerType &preconditioner)
{
  using number = typename VectorType::value_type;

  SolverControl::State conv = SolverControl::iterate;

  LogStream::Prefix prefix("cg");

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer g_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer d_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer h_pointer(this->memory);

  // define some aliases for simpler access
  VectorType &g = *g_pointer;
  VectorType &d = *d_pointer;
  VectorType &h = *h_pointer;

  // Should we build the matrix for eigenvalue computations?
  const bool do_eigenvalues =
    !condition_number_signal.empty() || !all_condition_numbers_signal.empty() ||
    !eigenvalues_signal.empty() || !all_eigenvalues_signal.empty();

  // vectors used for eigenvalue computations
  std::vector<typename VectorType::value_type> diagonal;
  std::vector<typename VectorType::value_type> offdiagonal;

  typename VectorType::value_type eigen_beta_alpha = 0;

  // resize the vectors, but do not set the values since they'd be overwritten
  // soon anyway.
  g.reinit(x, true);
  d.reinit(x, true);
  h.reinit(x, true);

  int    it        = 0;
  number gh        = number();
  number beta      = number();
  number alpha     = number();
  number old_alpha = number();

  // compute residual. if vector is zero, then short-circuit the full
  // computation
  if (!x.all_zero())
    {
      A.vmult(g, x);
      g.add(-1., b);
    }
  else
    g.equ(-1., b);

  double res = g.l2_norm();
  conv       = this->iteration_status(0, res, x);
  if (conv != SolverControl::iterate)
    return;

  while (conv == SolverControl::iterate)
    {
      it++;
      old_alpha = alpha;

      if (it > 1)
        {
          if (std::is_same<PreconditionerType, PreconditionIdentity>::value ==
              false)
            {
              preconditioner.vmult(h, g);
              beta = gh;
              Assert(std::abs(beta) != 0., ExcDivideByZero());
              gh   = g * h;
              beta = gh / beta;
              d.sadd(beta, -1., h);
            }
          else
            {
              beta = gh;
              gh   = res * res;
              beta = gh / beta;
              d.sadd(beta, -1., g);
            }
        }
      else
        {
          if (std::is_same<PreconditionerType, PreconditionIdentity>::value ==
              false)
            {
              preconditioner.vmult(h, g);
              d.equ(-1., h);
              gh = g * h;
            }
          else
            {
              d.equ(-1., g);
              gh = res * res;
            }
        }

      A.vmult(h, d);

      alpha = d * h;
      Assert(std::abs(alpha) != 0., ExcDivideByZero());
      alpha = gh / alpha;

      x.add(alpha, d);
      res = std::sqrt(std::abs(g.add_and_dot(alpha, h, g)));

      print_vectors(it, x, g, d);

      if (it > 1)
        {
          this->coefficients_signal(old_alpha, beta);
          // set up the vectors containing the diagonal and the off diagonal of
          // the projected matrix.
          if (do_eigenvalues)
            {
              diagonal.push_back(number(1.) / old_alpha + eigen_beta_alpha);
              eigen_beta_alpha = beta / old_alpha;
              offdiagonal.push_back(std::sqrt(beta) / old_alpha);
            }
          compute_eigs_and_cond(diagonal,
                                offdiagonal,
                                all_eigenvalues_signal,
                                all_condition_numbers_signal);
        }

      conv = this->iteration_status(it, res, x);
    }

  compute_eigs_and_cond(diagonal,
                        offdiagonal,
                        eigenvalues_signal,
                        condition_number_signal);

  // in case of failure: throw exception
  if (conv != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(it, res));
  // otherwise exit as normal
}



template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_coefficients_slot(
  const std::function<void(typename VectorType::value_type,
                           typename VectorType::value_type)> &slot)
{
  return coefficients_signal.connect(slot);
}



template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_condition_number_slot(
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



template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_eigenvalues_slot(
  const std::function<void(const std::vector<double> &)> &slot,
  const bool                                              every_iteration)
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



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


