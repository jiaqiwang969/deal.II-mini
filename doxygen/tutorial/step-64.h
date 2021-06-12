/**
@page step_64 The step-64 tutorial program
This tutorial depends on step-7, step-37.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The test case</a>
        <li><a href="#Movingdatatoandfromthedevice">Moving data to and from the device</a>
        <li><a href="#Matrixvectorproductimplementation">Matrix-vector product implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ClasscodeVaryingCoefficientFunctorcode">Class <code>VaryingCoefficientFunctor</code></a>
        <li><a href="#ClasscodeHelmholtzOperatorQuadcode">Class <code>HelmholtzOperatorQuad</code></a>
        <li><a href="#ClasscodeLocalHelmholtzOperatorcode">Class <code>LocalHelmholtzOperator</code></a>
        <li><a href="#ClasscodeHelmholtzOperatorcode">Class <code>HelmholtzOperator</code></a>
        <li><a href="#ClasscodeHelmholtzProblemcode">Class <code>HelmholtzProblem</code></a>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-64/doc/intro.dox

 <br> 

<i>
This program was contributed by Bruno Turcksin and Daniel Arndt, Oak Ridge National Laboratory.
</i>




<a name="Introduction"></a><h1>Introduction</h1>


这个例子展示了如何使用CUDA在GPU上实现超立方体上系数可变的亥姆霍兹方程的无矩阵方法。该线性系统将使用共轭梯度法进行求解，并通过MPI进行并行化。

在过去的几年里，一般的异构计算，特别是GPU，已经获得了很多的青睐。这是因为在给定的功率预算下，GPU比CPU提供更好的计算能力和内存带宽。在2019年初的架构中，对于PDE相关的任务，GPU的功率效率约为服务器CPU的2-3倍，宽<a
href="https://en.wikipedia.org/wiki/SIMD">SIMD</a>。GPU也是机器学习中最受欢迎的架构。另一方面，GPU并不容易编程。这个程序探索了deal.II的能力，看看这样的程序可以如何有效地实现。

虽然我们试图让CPU和GPU的无矩阵类的接口尽可能接近，但还是有一些区别。当在GPU上使用无矩阵框架时，人们必须编写一些CUDA代码。然而，其数量相当少，而且CUDA的使用仅限于几个关键词。




<a name="Thetestcase"></a><h3>The test case</h3>


在这个例子中，我们考虑亥姆霍兹问题@f{eqnarray*} - \nabla \cdot
\nabla u + a(\mathbf x) u &=&1,\\ u &=& 0 \quad \text{on } \partial \Omega @f} 。

其中 $a(\mathbf x)$ 是一个可变系数。

我们选择 $\Omega=[0,1]^3$ 和 $a(\mathbf x)=\frac{10}{0.05 +
2\|\mathbf x\|^2}$ 作为域。由于系数是围绕原点对称的，但域却不是，我们最终会得到一个非对称的解决方案。

如果你在本教程中读到这里，你就会知道这个问题的弱式表述是怎样的，以及原则上是怎样为它组建线性系统的。当然，在这个程序中，我们实际上不会形成矩阵，而只是表示它与之相乘时的作用。




<a name="Movingdatatoandfromthedevice"></a><h3>Moving data to and from the device</h3>


GPU（我们从现在开始用 "设备 "一词来指代GPU）有自己的内存，与CPU（我们从现在开始用 "主机 "一词）可访问的内存分开。设备上的正常计算可以分为三个独立的步骤。

-# 数据从主机移到设备上。

-#计算是在设备上完成的。

-# 结果从设备移回主机。

数据移动可以由用户代码显式完成，也可以使用UVM（统一虚拟内存）自动完成。在deal.II中，只支持第一种方法。虽然这意味着用户有额外的负担，但这可以更好地控制数据移动，更重要的是可以避免在主机而不是设备上错误地运行重要的内核。

deal.II中的数据移动是使用 LinearAlgebra::ReadWriteVector. 完成的，这些向量可以被看作是主机上的缓冲区，用于存储从设备接收的数据或向设备发送数据。有两种类型的向量可以在设备上使用。

-  LinearAlgebra::CUDAWrappers::Vector, ，它类似于更常见的Vector<Number>，和

-  LinearAlgebra::distributed::Vector<Number,   MemorySpace::CUDA>, 这是一个普通的 LinearAlgebra::distributed::Vector ，我们已经指定了要使用哪个内存空间。

如果没有指定内存空间，默认为 MemorySpace::Host. 。

接下来，我们展示如何使用 LinearAlgebra::CUDAWrappers::Vector: 将数据移入/移出设备。

@code
  unsigned int size = 10;
  LinearAlgebra::ReadWriteVector<double> rw_vector(size);


  ...do something with the rw_vector...


  // Move the data to the device:
  LinearAlgebra::CUDAWrappers::Vector<double> vector_dev(size);
  vector_dev.import(rw_vector, VectorOperations::insert);


  ...do some computations on the device...


  // Move the data back to the host:
  rw_vector.import(vector_dev, VectorOperations::insert);
@endcode

这里使用的两个向量类都只在一台机器上工作，也就是说，一个内存空间在主机上，一个在设备上。

但在有些情况下，人们希望在一些机器上的多个MPI进程之间运行并行计算，而每个机器都配备了GPU。在这种情况下，人们希望使用 `LinearAlgebra::distributed::Vector<Number,MemorySpace::CUDA>`, ，它是类似的，但`import()`阶段可能涉及MPI通信。

@code
  IndexSet locally_owned_dofs, locally_relevant_dofs;
  ...fill the two IndexSet objects...


  // Create the ReadWriteVector using an IndexSet instead of the size
  LinearAlgebra::ReadWriteVector<double> owned_rw_vector(locally_owned_dofs);


  ...do something with the rw_vector...


  // Move the data to the device:
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
    distributed_vector_dev(locally_owned_dofs, MPI_COMM_WORLD);
  distributed_vector_dev.import(owned_rw_vector, VectorOperations::insert);


  ...do something with the dev_vector...


  // Create a ReadWriteVector with a different IndexSet:
  LinearAlgebra::ReadWriteVector<double>
    relevant_rw_vector(locally_relevant_dofs);


  // Move the data to the host, possibly using MPI communication:
  relevant_rw_vector.import(distributed_vector_dev, VectorOperations::insert);
@endcode

`relevant_rw_vector`是一个存储向量所有元素的子集的对象。通常情况下，这些是 @ref GlossLocallyRelevantDof "本地相关的DoF"，这意味着它们在不同的MPI进程之间是重叠的。因此，一台机器上存储在该向量中的元素可能与该机器上的GPU存储的元素不一致，需要MPI通信来导入它们。

在所有这些情况下，在导入矢量时，可以插入数值（使用 VectorOperation::insert) 或添加到矢量的先前内容中（使用 VectorOperation::add).  ）。




<a name="Matrixvectorproductimplementation"></a><h3>Matrix-vector product implementation</h3>


在设备上评估无矩阵算子所需的代码与主机上的代码非常相似。然而，也有一些区别，主要是Step-37中的`local_apply()`函数和正交点的循环都需要封装在自己的函数中。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * First include the necessary files from the deal.II library known from the
 * previous tutorials.
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/tria.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_cg.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * The following ones include the data structures for the
 * implementation of matrix-free methods on GPU:
 * 
 * @code
 * #include <deal.II/base/cuda.h>
 * 
 * #include <deal.II/matrix_free/cuda_fe_evaluation.h>
 * #include <deal.II/matrix_free/cuda_matrix_free.h>
 * #include <deal.II/matrix_free/operators.h>
 * 
 * #include <fstream>
 * 
 * 
 * @endcode
 * 
 * As usual, we enclose everything into a namespace of its own:
 * 
 * @code
 * namespace Step64
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClasscodeVaryingCoefficientFunctorcode"></a> 
 * <h3>Class <code>VaryingCoefficientFunctor</code></h3>
 * 

 * 
 * Next, we define a class that implements the varying coefficients
 * we want to use in the Helmholtz operator. Later, we want to pass
 * an object of this type to a CUDAWrappers::MatrixFree
 * object that expects the class to have an `operator()` that fills the
 * values provided in the constructor for a given cell. This operator
 * needs to run on the device, so it needs to be marked as `__device__`
 * for the compiler.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class VaryingCoefficientFunctor
 *   {
 *   public:
 *     VaryingCoefficientFunctor(double *coefficient)
 *       : coef(coefficient)
 *     {}
 * 
 *     __device__ void operator()(
 *       const unsigned int                                          cell,
 *       const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data);
 * 
 * @endcode
 * 
 * Since CUDAWrappers::MatrixFree::Data doesn't know about the size of its
 * arrays, we need to store the number of quadrature points and the numbers
 * of degrees of freedom in this class to do necessary index conversions.
 * 
 * @code
 *     static const unsigned int n_dofs_1d = fe_degree + 1;
 *     static const unsigned int n_local_dofs =
 *       dealii::Utilities::pow(n_dofs_1d, dim);
 *     static const unsigned int n_q_points =
 *       dealii::Utilities::pow(n_dofs_1d, dim);
 * 
 *   private:
 *     double *coef;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * The following function implements this coefficient. Recall from
 * the introduction that we have defined it as $a(\mathbf
 * x)=\frac{10}{0.05 + 2\|\mathbf x\|^2}$
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   __device__ void VaryingCoefficientFunctor<dim, fe_degree>::operator()(
 *     const unsigned int                                          cell,
 *     const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data)
 *   {
 *     const unsigned int pos = CUDAWrappers::local_q_point_id<dim, double>(
 *       cell, gpu_data, n_dofs_1d, n_q_points);
 *     const Point<dim> q_point =
 *       CUDAWrappers::get_quadrature_point<dim, double>(cell,
 *                                                       gpu_data,
 *                                                       n_dofs_1d);
 * 
 *     double p_square = 0.;
 *     for (unsigned int i = 0; i < dim; ++i)
 *       {
 *         const double coord = q_point[i];
 *         p_square += coord * coord;
 *       }
 *     coef[pos] = 10. / (0.05 + 2. * p_square);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClasscodeHelmholtzOperatorQuadcode"></a> 
 * <h3>Class <code>HelmholtzOperatorQuad</code></h3>
 * 

 * 
 * The class `HelmholtzOperatorQuad` implements the evaluation of
 * the Helmholtz operator at each quadrature point. It uses a
 * similar mechanism as the MatrixFree framework introduced in
 * step-37. In contrast to there, the actual quadrature point
 * index is treated implicitly by converting the current thread
 * index. As before, the functions of this class need to run on
 * the device, so need to be marked as `__device__` for the
 * compiler.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class HelmholtzOperatorQuad
 *   {
 *   public:
 *     __device__ HelmholtzOperatorQuad(double coef)
 *       : coef(coef)
 *     {}
 * 
 *     __device__ void
 *     operator()(CUDAWrappers::FEEvaluation<dim, fe_degree> *fe_eval) const;
 * 
 *   private:
 *     double coef;
 *   };
 * 
 * 
 * @endcode
 * 
 * The Helmholtz problem we want to solve here reads in weak form as follows:
 * @f{eqnarray*}
 * (\nabla v, \nabla u)+ (v, a(\mathbf x) u) &=&(v,1) \quad \forall v.
 * @f}
 * If you have seen step-37, then it will be obvious that
 * the two terms on the left-hand side correspond to the two function calls
 * here:
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   __device__ void HelmholtzOperatorQuad<dim, fe_degree>::
 *                   operator()(CUDAWrappers::FEEvaluation<dim, fe_degree> *fe_eval) const
 *   {
 *     fe_eval->submit_value(coef * fe_eval->get_value());
 *     fe_eval->submit_gradient(fe_eval->get_gradient());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClasscodeLocalHelmholtzOperatorcode"></a> 
 * <h3>Class <code>LocalHelmholtzOperator</code></h3>
 * 

 * 
 * Finally, we need to define a class that implements the whole operator
 * evaluation that corresponds to a matrix-vector product in matrix-based
 * approaches.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class LocalHelmholtzOperator
 *   {
 *   public:
 *     LocalHelmholtzOperator(double *coefficient)
 *       : coef(coefficient)
 *     {}
 * 
 *     __device__ void operator()(
 *       const unsigned int                                          cell,
 *       const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
 *       CUDAWrappers::SharedData<dim, double> *                     shared_data,
 *       const double *                                              src,
 *       double *                                                    dst) const;
 * 
 * @endcode
 * 
 * Again, the CUDAWrappers::MatrixFree object doesn't know about the number
 * of degrees of freedom and the number of quadrature points so we need
 * to store these for index calculations in the call operator.
 * 
 * @code
 *     static const unsigned int n_dofs_1d    = fe_degree + 1;
 *     static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
 *     static const unsigned int n_q_points   = Utilities::pow(fe_degree + 1, dim);
 * 
 *   private:
 *     double *coef;
 *   };
 * 
 * 
 * @endcode
 * 
 * This is the call operator that performs the Helmholtz operator evaluation
 * on a given cell similar to the MatrixFree framework on the CPU.
 * In particular, we need access to both values and gradients of the source
 * vector and we write value and gradient information to the destination
 * vector.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   __device__ void LocalHelmholtzOperator<dim, fe_degree>::operator()(
 *     const unsigned int                                          cell,
 *     const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
 *     CUDAWrappers::SharedData<dim, double> *                     shared_data,
 *     const double *                                              src,
 *     double *                                                    dst) const
 *   {
 *     const unsigned int pos = CUDAWrappers::local_q_point_id<dim, double>(
 *       cell, gpu_data, n_dofs_1d, n_q_points);
 * 
 *     CUDAWrappers::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double>
 *       fe_eval(cell, gpu_data, shared_data);
 *     fe_eval.read_dof_values(src);
 *     fe_eval.evaluate(true, true);
 *     fe_eval.apply_for_each_quad_point(
 *       HelmholtzOperatorQuad<dim, fe_degree>(coef[pos]));
 *     fe_eval.integrate(true, true);
 *     fe_eval.distribute_local_to_global(dst);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClasscodeHelmholtzOperatorcode"></a> 
 * <h3>Class <code>HelmholtzOperator</code></h3>
 * 

 * 
 * The `HelmholtzOperator` class acts as wrapper for
 * `LocalHelmholtzOperator` defining an interface that can be used
 * with linear solvers like SolverCG. In particular, like every
 * class that implements the interface of a linear operator, it
 * needs to have a `vmult()` function that performs the action of
 * the linear operator on a source vector.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class HelmholtzOperator
 *   {
 *   public:
 *     HelmholtzOperator(const DoFHandler<dim> &          dof_handler,
 *                       const AffineConstraints<double> &constraints);
 * 
 *     void
 *     vmult(LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &dst,
 *           const LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
 *             &src) const;
 * 
 *     void initialize_dof_vector(
 *       LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &vec) const;
 * 
 *   private:
 *     CUDAWrappers::MatrixFree<dim, double>       mf_data;
 *     LinearAlgebra::CUDAWrappers::Vector<double> coef;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * The following is the implementation of the constructor of this
 * class. In the first part, we initialize the `mf_data` member
 * variable that is going to provide us with the necessary
 * information when evaluating the operator.
 *   

 * 
 * In the second half, we need to store the value of the coefficient
 * for each quadrature point in every active, locally owned cell.
 * We can ask the parallel triangulation for the number of active, locally
 * owned cells but only have a DoFHandler object at hand. Since
 * DoFHandler::get_triangulation() returns a Triangulation object, not a
 * parallel::TriangulationBase object, we have to downcast the return value.
 * This is safe to do here because we know that the triangulation is a
 * parallel:distributed::Triangulation object in fact.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   HelmholtzOperator<dim, fe_degree>::HelmholtzOperator(
 *     const DoFHandler<dim> &          dof_handler,
 *     const AffineConstraints<double> &constraints)
 *   {
 *     MappingQGeneric<dim> mapping(fe_degree);
 *     typename CUDAWrappers::MatrixFree<dim, double>::AdditionalData
 *       additional_data;
 *     additional_data.mapping_update_flags = update_values | update_gradients |
 *                                            update_JxW_values |
 *                                            update_quadrature_points;
 *     const QGauss<1> quad(fe_degree + 1);
 *     mf_data.reinit(mapping, dof_handler, constraints, quad, additional_data);
 * 
 * 
 *     const unsigned int n_owned_cells =
 *       dynamic_cast<const parallel::TriangulationBase<dim> *>(
 *         &dof_handler.get_triangulation())
 *         ->n_locally_owned_active_cells();
 *     coef.reinit(Utilities::pow(fe_degree + 1, dim) * n_owned_cells);
 * 
 *     const VaryingCoefficientFunctor<dim, fe_degree> functor(coef.get_values());
 *     mf_data.evaluate_coefficients(functor);
 *   }
 * 
 * 
 * @endcode
 * 
 * The key step then is to use all of the previous classes to loop over
 * all cells to perform the matrix-vector product. We implement this
 * in the next function.
 *   

 * 
 * When applying the Helmholtz operator, we have to be careful to handle
 * boundary conditions correctly. Since the local operator doesn't know about
 * constraints, we have to copy the correct values from the source to the
 * destination vector afterwards.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void HelmholtzOperator<dim, fe_degree>::vmult(
 *     LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &      dst,
 *     const LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &src)
 *     const
 *   {
 *     dst = 0.;
 *     LocalHelmholtzOperator<dim, fe_degree> helmholtz_operator(
 *       coef.get_values());
 *     mf_data.cell_loop(helmholtz_operator, src, dst);
 *     mf_data.copy_constrained_values(src, dst);
 *   }
 * 
 * 
 * 
 *   template <int dim, int fe_degree>
 *   void HelmholtzOperator<dim, fe_degree>::initialize_dof_vector(
 *     LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &vec) const
 *   {
 *     mf_data.initialize_dof_vector(vec);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClasscodeHelmholtzProblemcode"></a> 
 * <h3>Class <code>HelmholtzProblem</code></h3>
 * 

 * 
 * This is the main class of this program. It defines the usual
 * framework we use for tutorial programs. The only point worth
 * commenting on is the `solve()` function and the choice of vector
 * types.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class HelmholtzProblem
 *   {
 *   public:
 *     HelmholtzProblem();
 * 
 *     void run();
 * 
 *   private:
 *     void setup_system();
 * 
 *     void assemble_rhs();
 * 
 *     void solve();
 * 
 *     void output_results(const unsigned int cycle) const;
 * 
 *     MPI_Comm mpi_communicator;
 * 
 *     parallel::distributed::Triangulation<dim> triangulation;
 * 
 *     FE_Q<dim>       fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     IndexSet locally_owned_dofs;
 *     IndexSet locally_relevant_dofs;
 * 
 *     AffineConstraints<double>                          constraints;
 *     std::unique_ptr<HelmholtzOperator<dim, fe_degree>> system_matrix_dev;
 * 
 * @endcode
 * 
 * Since all the operations in the `solve()` function are executed on the
 * graphics card, it is necessary for the vectors used to store their values
 * on the GPU as well. LinearAlgebra::distributed::Vector can be told which
 * memory space to use. There is also LinearAlgebra::CUDAWrappers::Vector
 * that always uses GPU memory storage but doesn't work with MPI. It might
 * be worth noticing that the communication between different MPI processes
 * can be improved if the MPI implementation is CUDA-aware and the configure
 * flag `DEAL_II_MPI_WITH_CUDA_SUPPORT` is enabled. (The value of this
 * flag needs to be set at the time you call `cmake` when installing
 * deal.II.)
 *     

 * 
 * In addition, we also keep a solution vector with CPU storage such that we
 * can view and display the solution as usual.
 * 
 * @code
 *     LinearAlgebra::distributed::Vector<double, MemorySpace::Host>
 *                                                                   ghost_solution_host;
 *     LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> solution_dev;
 *     LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
 *       system_rhs_dev;
 * 
 *     ConditionalOStream pcout;
 *   };
 * 
 * 
 * @endcode
 * 
 * The implementation of all the remaining functions of this class apart from
 * `Helmholtzproblem::solve()` doesn't contain anything new and we won't
 * further comment much on the overall approach.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   HelmholtzProblem<dim, fe_degree>::HelmholtzProblem()
 *     : mpi_communicator(MPI_COMM_WORLD)
 *     , triangulation(mpi_communicator)
 *     , fe(fe_degree)
 *     , dof_handler(triangulation)
 *     , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {}
 * 
 * 
 * 
 *   template <int dim, int fe_degree>
 *   void HelmholtzProblem<dim, fe_degree>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     locally_owned_dofs = dof_handler.locally_owned_dofs();
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 *     system_rhs_dev.reinit(locally_owned_dofs, mpi_communicator);
 * 
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(),
 *                                              constraints);
 *     constraints.close();
 * 
 *     system_matrix_dev.reset(
 *       new HelmholtzOperator<dim, fe_degree>(dof_handler, constraints));
 * 
 *     ghost_solution_host.reinit(locally_owned_dofs,
 *                                locally_relevant_dofs,
 *                                mpi_communicator);
 *     system_matrix_dev->initialize_dof_vector(solution_dev);
 *     system_rhs_dev.reinit(solution_dev);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Unlike programs such as step-4 or step-6, we will not have to
 * assemble the whole linear system but only the right hand side
 * vector. This looks in essence like we did in step-4, for example,
 * but we have to pay attention to using the right constraints
 * object when copying local contributions into the global
 * vector. In particular, we need to make sure the entries that
 * correspond to boundary nodes are properly zeroed out. This is
 * necessary for CG to converge.  (Another solution would be to
 * modify the `vmult()` function above in such a way that we pretend
 * the source vector has zero entries by just not taking them into
 * account in matrix-vector products. But the approach used here is
 * simpler.)
 *   

 * 
 * At the end of the function, we can't directly copy the values
 * from the host to the device but need to use an intermediate
 * object of type LinearAlgebra::ReadWriteVector to construct the
 * correct communication pattern necessary.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void HelmholtzProblem<dim, fe_degree>::assemble_rhs()
 *   {
 *     LinearAlgebra::distributed::Vector<double, MemorySpace::Host>
 *                       system_rhs_host(locally_owned_dofs,
 *                       locally_relevant_dofs,
 *                       mpi_communicator);
 *     const QGauss<dim> quadrature_formula(fe_degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points |
 *                               update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     Vector<double> cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell_rhs = 0;
 * 
 *           fe_values.reinit(cell);
 * 
 *           for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
 *             {
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 cell_rhs(i) += (fe_values.shape_value(i, q_index) * 1.0 *
 *                                 fe_values.JxW(q_index));
 *             }
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           constraints.distribute_local_to_global(cell_rhs,
 *                                                  local_dof_indices,
 *                                                  system_rhs_host);
 *         }
 *     system_rhs_host.compress(VectorOperation::add);
 * 
 *     LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
 *     rw_vector.import(system_rhs_host, VectorOperation::insert);
 *     system_rhs_dev.import(rw_vector, VectorOperation::insert);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This solve() function finally contains the calls to the new classes
 * previously discussed. Here we don't use any preconditioner, i.e.,
 * precondition by the identity matrix, to focus just on the peculiarities of
 * the CUDAWrappers::MatrixFree framework. Of course, in a real application
 * the choice of a suitable preconditioner is crucial but we have at least the
 * same restrictions as in step-37 since matrix entries are computed on the
 * fly and not stored.
 *   

 * 
 * After solving the linear system in the first part of the function, we
 * copy the solution from the device to the host to be able to view its
 * values and display it in `output_results()`. This transfer works the same
 * as at the end of the previous function.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void HelmholtzProblem<dim, fe_degree>::solve()
 *   {
 *     PreconditionIdentity preconditioner;
 * 
 *     SolverControl solver_control(system_rhs_dev.size(),
 *                                  1e-12 * system_rhs_dev.l2_norm());
 *     SolverCG<LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>> cg(
 *       solver_control);
 *     cg.solve(*system_matrix_dev, solution_dev, system_rhs_dev, preconditioner);
 * 
 *     pcout << "  Solved in " << solver_control.last_step() << " iterations."
 *           << std::endl;
 * 
 *     LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
 *     rw_vector.import(solution_dev, VectorOperation::insert);
 *     ghost_solution_host.import(rw_vector, VectorOperation::insert);
 * 
 *     constraints.distribute(ghost_solution_host);
 * 
 *     ghost_solution_host.update_ghost_values();
 *   }
 * 
 * @endcode
 * 
 * The output results function is as usual since we have already copied the
 * values back from the GPU to the CPU.
 *   

 * 
 * While we're already doing something with the function, we might
 * as well compute the $L_2$ norm of the solution. We do this by
 * calling VectorTools::integrate_difference(). That function is
 * meant to compute the error by evaluating the difference between
 * the numerical solution (given by a vector of values for the
 * degrees of freedom) and an object representing the exact
 * solution. But we can easily compute the $L_2$ norm of the
 * solution by passing in a zero function instead. That is, instead
 * of evaluating the error $\|u_h-u\|_{L_2(\Omega)}$, we are just
 * evaluating $\|u_h-0\|_{L_2(\Omega)}=\|u_h\|_{L_2(\Omega)}$
 * instead.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void HelmholtzProblem<dim, fe_degree>::output_results(
 *     const unsigned int cycle) const
 *   {
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(ghost_solution_host, "solution");
 *     data_out.build_patches();
 * 
 *     DataOutBase::VtkFlags flags;
 *     flags.compression_level = DataOutBase::VtkFlags::best_speed;
 *     data_out.set_flags(flags);
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", cycle, mpi_communicator, 2);
 * 
 *     Vector<float> cellwise_norm(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(dof_handler,
 *                                       ghost_solution_host,
 *                                       Functions::ZeroFunction<dim>(),
 *                                       cellwise_norm,
 *                                       QGauss<dim>(fe.degree + 2),
 *                                       VectorTools::L2_norm);
 *     const double global_norm =
 *       VectorTools::compute_global_error(triangulation,
 *                                         cellwise_norm,
 *                                         VectorTools::L2_norm);
 *     pcout << "  solution norm: " << global_norm << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * There is nothing surprising in the `run()` function either. We simply
 * compute the solution on a series of (globally) refined meshes.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void HelmholtzProblem<dim, fe_degree>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 7 - dim; ++cycle)
 *       {
 *         pcout << "Cycle " << cycle << std::endl;
 * 
 *         if (cycle == 0)
 *           GridGenerator::hyper_cube(triangulation, 0., 1.);
 *         triangulation.refine_global(1);
 * 
 *         setup_system();
 * 
 *         pcout << "   Number of active cells:       "
 *               << triangulation.n_global_active_cells() << std::endl
 *               << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl;
 * 
 *         assemble_rhs();
 *         solve();
 *         output_results(cycle);
 *         pcout << std::endl;
 *       }
 *   }
 * } // namespace Step64
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main()</code> function</h3>
 * 

 * 
 * Finally for the `main()` function.  By default, all the MPI ranks
 * will try to access the device with number 0, which we assume to be
 * the GPU device associated with the CPU on which a particular MPI
 * rank runs. This works, but if we are running with MPI support it
 * may be that multiple MPI processes are running on the same machine
 * (for example, one per CPU core) and then they would all want to
 * access the same GPU on that machine. If there is only one GPU in
 * the machine, there is nothing we can do about it: All MPI ranks on
 * that machine need to share it. But if there are more than one GPU,
 * then it is better to address different graphic cards for different
 * processes. The choice below is based on the MPI process id by
 * assigning GPUs round robin to GPU ranks. (To work correctly, this
 * scheme assumes that the MPI ranks on one machine are
 * consecutive. If that were not the case, then the rank-GPU
 * association may just not be optimal.) To make this work, MPI needs
 * to be initialized before using this function.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace Step64;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
 * 
 *       int         n_devices       = 0;
 *       cudaError_t cuda_error_code = cudaGetDeviceCount(&n_devices);
 *       AssertCuda(cuda_error_code);
 *       const unsigned int my_mpi_id =
 *         Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
 *       const int device_id = my_mpi_id % n_devices;
 *       cuda_error_code     = cudaSetDevice(device_id);
 *       AssertCuda(cuda_error_code);
 * 
 *       HelmholtzProblem<3, 3> helmholtz_problem;
 *       helmholtz_problem.run();
 *     }
 *   catch (std::exception &exc)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Exception on processing: " << std::endl
 *                 << exc.what() << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       return 1;
 *     }
 *   catch (...)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Unknown exception!" << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       return 1;
 *     }
 * 
 *   return 0;
 * }
 * @endcode
examples/step-64/doc/results.dox



<a name="Results"></a><h1>Results</h1>


由于本教程的主要目的是演示如何使用 CUDAWrappers::MatrixFree 接口，而不是计算任何有用的东西本身，所以我们在这里只是显示预期的输出。

@code
Cycle 0
   Number of active cells:       8
   Number of degrees of freedom: 343
  Solved in 27 iterations.
  solution norm: 0.0205439


Cycle 1
   Number of active cells:       64
   Number of degrees of freedom: 2197
  Solved in 60 iterations.
  solution norm: 0.0205269


Cycle 2
   Number of active cells:       512
   Number of degrees of freedom: 15625
  Solved in 114 iterations.
  solution norm: 0.0205261


Cycle 3
   Number of active cells:       4096
   Number of degrees of freedom: 117649
  Solved in 227 iterations.
  solution norm: 0.0205261
@endcode



在这里，我们可以做两个观察。首先，数值解的准则在收敛，大概是收敛到精确（但未知）解的准则。其次，每次细化网格时，迭代次数大约增加一倍。这与CG迭代次数随矩阵条件数的平方根增长的预期一致；而且我们知道二阶微分运算的矩阵条件数的增长方式是 ${\cal O}(h^{-2})$  。这当然是相当低效的，因为一个最佳解算器的迭代次数与问题的大小无关。但是要有这样一个求解器，就需要使用比我们在这里使用的身份矩阵更好的预处理。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


目前，这个程序完全没有使用预处理程序。这主要是因为构建一个高效的无矩阵预处理程序是不容易的。  然而，只需要相应矩阵的对角线的简单选择是很好的选择，这些也可以用无矩阵的方式计算。另外，也许更好的是，我们可以扩展教程，使用类似步骤37的切比雪夫平滑器的多重网格。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-64.cu"
*/
