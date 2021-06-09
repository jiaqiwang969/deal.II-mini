//include/deal.II-translator/lac/ginkgo_solver_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#ifndef dealii_ginkgo_solver_h
#  define dealii_ginkgo_solver_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_GINKGO

#    include <deal.II/lac/block_sparse_matrix.h>
#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/solver_control.h>
#    include <deal.II/lac/sparse_matrix.h>
#    include <deal.II/lac/vector.h>

#    include <ginkgo/ginkgo.hpp>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

namespace GinkgoWrappers
{
  /**
   * 该类构成了Ginkgo所有迭代求解器的基类。
   * 各个派生类只接受特定于它们的额外数据，并解决给定的线性系统。Ginkgo实现的全部求解器集合可在<a
   * Ginkgo href="https://ginkgo-project.github.io/ginkgo/doc/develop/">
   * documentation and manual pages</a>中找到。
   * @ingroup GinkgoWrappers
   *
   */
  template <typename ValueType, typename IndexType>
  class SolverBase
  {
  public:
    /**
     * 构造器。         @p exec_type
     * 定义了计算解决方案的范式。
     * 它是一个字符串，可以选择 "OMP"、"参考 "或 "Cuda"。
     * 各自的字符串创建各自的执行器，如下所述。
     * Ginkgo目前支持三种不同的执行器类型。        +
     * OmpExecutor指定数据应被存储并在支持OpenMP的设备（如主机CPU）上执行相关操作；
     * ``` auto omp =  gko::create<gko::OmpExecutor>();  ``` +
     * CudaExecutor指定数据应被存储并在NVIDIA
     * GPU加速器上执行操作； ```
     * if(gko::CudaExecutor::get_num_devices()  > 0 ){ auto cuda =
     * gko::create<gko::CudaExecutor>();  }     ``` +
     * ReferenceExecutor执行一个非优化的引用实现，可用于调试库。
     * ``` auto ref =  gko::create<gko::ReferenceExecutor>();  ```
     * 下面的代码片段演示了使用OpenMP执行器来创建一个求解器，它将使用OpenMP范式在CPU上求解系统。
     * ``` auto omp =  gko::create<gko::OmpExecutor>();  using cg =
     * gko::solver::Cg<>;  auto solver_gen =  cg::build()  .with_criteria(
     * gko::stop::Iteration::build().with_max_iters(20u).on(omp),
     * gko::stop::ResidualNormReduction<>::build()
     * ].with_reduction_factor(1e-6) .on(mp)) .on(mp); auto solver =
     * solver_gen->generate(system_matrix); solver->apply(lend(rhs),
     * lend(solution)); ``  @p solver_control  对象与其他
     * deal.II迭代解算器相同。
     *
     */
    SolverBase(SolverControl &solver_control, const std::string &exec_type);

    /**
     * 解构器。
     *
     */
    virtual ~SolverBase() = default;

    /**
     * 初始化矩阵并将其数据复制到Ginkgo的数据结构中。
     *
     */
    void
    initialize(const SparseMatrix<ValueType> &matrix);

    /**
     * 解决线性系统<tt>Ax=b</tt>。根据派生类提供的信息，选择Ginkgo的一个线性求解器。
     *
     */
    void
    apply(Vector<ValueType> &solution, const Vector<ValueType> &rhs);

    /**
     * 解决线性系统<tt>Ax=b</tt>。根据派生类提供的信息，选择Ginkgo的一个线性求解器。
     *
     */
    void
    solve(const SparseMatrix<ValueType> &matrix,
          Vector<ValueType> &            solution,
          const Vector<ValueType> &      rhs);

    /**
     * 访问控制收敛的对象。
     *
     */
    SolverControl &
    control() const;


  protected:
    /**
     * 对控制迭代求解器收敛性的对象的引用。
     *
     */
    SolverControl &solver_control;

    /**
     * Ginkgo生成的求解器工厂对象。
     *
     */
    std::shared_ptr<gko::LinOpFactory> solver_gen;

    /**
     * 残差准则对象，该对象根据solver_control成员中设置的容忍度来控制残差的减少。
     *
     */
    std::shared_ptr<gko::stop::ResidualNormReduction<>::Factory>
      residual_criterion;

    /**
     * Ginkgo收敛记录器，用于检查收敛情况和其他需要的求解器数据。
     *
     */
    std::shared_ptr<gko::log::Convergence<>> convergence_logger;

    /**
     * Ginkgo组合工厂对象用于创建一个组合停顿准则，以传递给求解器。
     *
     */
    std::shared_ptr<gko::stop::Combined::Factory> combined_factory;

    /**
     * Ginkgo中的执行范式。可以在 `gko::OmpExecutor`,
     * `gko::CudaExecutor` 和 `gko::ReferenceExecutor`
     * 之间选择，更多细节可以在Ginkgo的文档中找到。
     *
     */
    std::shared_ptr<gko::Executor> executor;

  private:
    /**
     * 用事件掩码初始化Ginkgo记录器对象。参照<a
     * href="https://github.com/ginkgo-project/ginkgo/blob/develop/include/ginkgo/core/log/logger.hpp">Ginkgo's
     * logging event masks.</a>。
     *
     */
    void
    initialize_ginkgo_log();

    /**
     * Ginkgo矩阵数据结构。第一个模板参数是用于存储矩阵的非零点阵列。第二个是用于行指针和列索引。
     * @todo  基于矩阵类型的模板化。
     *
     */
    std::shared_ptr<gko::matrix::Csr<ValueType, IndexType>> system_matrix;

    /**
     * 执行范式为一个字符串，由用户设置。可以在`omp'、`cuda'和`reference'之间选择，更多细节可以在Ginkgo的文档中找到。
     *
     */
    const std::string exec_type;
  };


  /**
   * 使用Ginkgo CG解算器的解算器接口的实现。
   * @ingroup GinkgoWrappers
   *
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverCG : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * 一个标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造函数。          @param[in,out]  solver_control
     * 然后，求解器控制对象被用来设置参数，并从求解线性系统的CG工厂设置CG求解器。
     * @param[in]  exec_type CG求解器的执行范式。          @param[in]
     * data 解算器所需的额外数据。
     *
     */
    SolverCG(SolverControl &       solver_control,
             const std::string &   exec_type,
             const AdditionalData &data = AdditionalData());

    /**
     * 构造器。          @param[in,out]  solver_control
     * 解算器控制对象然后用于设置参数，并从解算线性系统的CG工厂设置CG解算器。
     * @param[in]  exec_type CG求解器的执行范式。          @param[in]
     * preconditioner 解算器的预处理程序。          @param[in]  data
     * 解算器所需的额外数据。
     *
     */
    SolverCG(SolverControl &                           solver_control,
             const std::string &                       exec_type,
             const std::shared_ptr<gko::LinOpFactory> &preconditioner,
             const AdditionalData &                    data = AdditionalData());

  protected:
    /**
     * 存储该特定求解器的设置副本。
     *
     */
    const AdditionalData additional_data;
  };


  /**
   * 使用Ginkgo Bicgstab解算器的解算器接口的实现。
   * @ingroup GinkgoWrappers
   *
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverBicgstab : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * 一个标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造函数。          @param[in,out]  solver_control
     * 然后，求解器控制对象被用来设置参数，并从Bicgstab工厂设置Bicgstab求解器，从而求解线性系统。
     * @param[in]  exec_type Bicgstab求解器的执行范式。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverBicgstab(SolverControl &       solver_control,
                   const std::string &   exec_type,
                   const AdditionalData &data = AdditionalData());

    /**
     * 构造器。          @param[in,out]  solver_control
     * 解算器控制对象然后用于设置参数，并从Bicgstab工厂设置Bicgstab解算器，该解算器解决线性系统。
     * @param[in]  exec_type Bicgstab求解器的执行范式。
     * @param[in]  preconditioner 解算器的预处理程序。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverBicgstab(SolverControl &                           solver_control,
                   const std::string &                       exec_type,
                   const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                   const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储该特定求解器的设置副本。
     *
     */
    const AdditionalData additional_data;
  };

  /**
   * 使用Ginkgo CGS求解器的求解器接口的实现。
   * CGS或共轭梯度平方法是一种迭代型Krylov子空间方法，适用于一般系统。
   * @ingroup GinkgoWrappers
   *
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverCGS : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * 一个标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。          @param[in,out]  solver_control
     * 然后，求解器控制对象被用来设置参数，并从CGS工厂设置CGS求解器，从而求解线性系统。
     * @param[in]  exec_type CGS求解器的执行范式。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverCGS(SolverControl &       solver_control,
              const std::string &   exec_type,
              const AdditionalData &data = AdditionalData());

    /**
     * 构造器。          @param[in,out]  solver_control
     * 解算器控制对象，然后用于设置参数，并从CGS工厂设置CGS解算器，解算线性系统。
     * @param[in]  exec_type CGS求解器的执行范式。
     * @param[in]  preconditioner 解算器的预处理程序。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverCGS(SolverControl &                           solver_control,
              const std::string &                       exec_type,
              const std::shared_ptr<gko::LinOpFactory> &preconditioner,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储该特定求解器的设置副本。
     *
     */
    const AdditionalData additional_data;
  };

  /**
   * 使用Ginkgo FCG求解器的求解器接口的实现。
   * FCG或灵活共轭梯度法是一种迭代型Krylov子空间方法，适用于对称正定法。
   * 虽然这种方法对对称正定矩阵表现非常好，但一般来说，它不适合一般的矩阵。
   * 与基于Polack-Ribiere公式的标准CG相比，灵活的CG使用Fletcher-Reeves公式来创建跨越Krylov子空间的正态向量。这增加了每个Krylov求解器迭代的计算成本，但允许使用非恒定的预处理程序。
   * @ingroup GinkgoWrappers
   *
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverFCG : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * 一个标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。          @param[in,out]  solver_control
     * 然后，求解器控制对象被用来设置参数，并从FCG工厂设置FCG求解器，该求解器解决线性系统。
     * @param[in]  exec_type FCG求解器的执行范式。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverFCG(SolverControl &       solver_control,
              const std::string &   exec_type,
              const AdditionalData &data = AdditionalData());

    /**
     * 构造器。          @param[in,out]  solver_control
     * 解算器控制对象，然后用于设置参数，并从FCG工厂设置FCG解算器，解算线性系统。
     * @param[in]  exec_type FCG求解器的执行范式。
     * @param[in]  preconditioner 解算器的预处理程序。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverFCG(SolverControl &                           solver_control,
              const std::string &                       exec_type,
              const std::shared_ptr<gko::LinOpFactory> &preconditioner,
              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储该特定求解器的设置副本。
     *
     */
    const AdditionalData additional_data;
  };

  /**
   * 使用Ginkgo GMRES求解器的求解器接口的实现。
   * @ingroup GinkgoWrappers
   *
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverGMRES : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * 一个标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，设置临时向量的数量为30，即每30次迭代做一次重启。
       *
       */
      AdditionalData(const unsigned int restart_parameter = 30);

      /**
       * 临时向量的最大数量。
       *
       */
      unsigned int restart_parameter;
    };

    /**
     * 构造函数。          @param[in,out]  solver_control
     * 然后，求解器控制对象被用来设置参数，并从GMRES工厂设置GMRES求解器，求解线性系统。
     * @param[in]  exec_type GMRES求解器的执行范式。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverGMRES(SolverControl &       solver_control,
                const std::string &   exec_type,
                const AdditionalData &data = AdditionalData());

    /**
     * 构造器。          @param[in,out]  solver_control
     * 解算器控制对象然后用于设置参数，并从GMRES工厂设置GMRES解算器，该解算器解决线性系统。
     * @param[in]  exec_type GMRES求解器的执行范式。
     * @param[in]  preconditioner 解算器的预处理程序。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverGMRES(SolverControl &                           solver_control,
                const std::string &                       exec_type,
                const std::shared_ptr<gko::LinOpFactory> &preconditioner,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储该特定求解器的设置副本。
     *
     */
    const AdditionalData additional_data;
  };

  /**
   * 使用Ginkgo IR求解器的求解器接口的实现。
   * 迭代细化（IR）是一种迭代方法，它使用另一种粗略的方法，通过当前的残差对当前解的误差进行近似。
   * @ingroup GinkgoWrappers
   *
   */
  template <typename ValueType = double, typename IndexType = int32_t>
  class SolverIR : public SolverBase<ValueType, IndexType>
  {
  public:
    /**
     * 一个标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造函数。          @param[in,out]  solver_control
     * 然后，求解器控制对象被用来设置参数，并从解决线性系统的IR工厂设置IR求解器。
     * @param[in]  exec_type 红外求解器的执行范式。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverIR(SolverControl &       solver_control,
             const std::string &   exec_type,
             const AdditionalData &data = AdditionalData());

    /**
     * 构造器。          @param[in,out]  solver_control
     * 然后，求解器控制对象被用来设置参数，并从求解线性系统的IR工厂设置IR求解器。
     * @param[in]  exec_type 红外求解器的执行范式。
     * @param[in]  inner_solver 用于IR求解器的内部求解器。
     * @param[in]  data 解算器所需的额外数据。
     *
     */
    SolverIR(SolverControl &                           solver_control,
             const std::string &                       exec_type,
             const std::shared_ptr<gko::LinOpFactory> &inner_solver,
             const AdditionalData &                    data = AdditionalData());

  protected:
    /**
     * 存储该特定求解器的设置副本。
     *
     */
    const AdditionalData additional_data;
  };


} // namespace GinkgoWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_GINKGO

#endif
 /*----------------------------   ginkgo_solver.h ---------------------------*/ 


