/**
@page step_76 The step-76 tutorial program
This tutorial depends on step-67.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#SharedmemoryandhybridparallelizationwithMPI30">Shared-memory and hybrid parallelization with MPI-3.0</a>
      <ul>
        <li><a href="#Motivation">Motivation</a>
        <li><a href="#BasicMPI30commands">Basic MPI-3.0 commands</a>
        <li><a href="#MPI30andLinearAlgebradistributedVector">MPI-3.0 and LinearAlgebra::distributed::Vector</a>
        <li><a href="#MPI30andMatrixFree">MPI-3.0 and MatrixFree</a>
      </ul>
        <li><a href="#Cellcentricloops">Cell-centric loops</a>
      <ul>
        <li><a href="#MotivationFCLvsCCL">Motivation: FCL vs. CCL</a>
        <li><a href="#CellcentricloopsandMatrixFree">Cell-centric loops and MatrixFree</a>
      </ul>
        <li><a href="#ProvidinglambdastoMatrixFreeloops">Providing lambdas to MatrixFree loops</a>
        <li><a href="#VectorizedArrayType">VectorizedArrayType</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Parametersandutilityfunctions">Parameters and utility functions</a>
        <li><a href="#EuleroperatorusingacellcentricloopandMPI30sharedmemory">Euler operator using a cell-centric loop and MPI-3.0 shared memory</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#ExtensiontothecompressibleNavierStokesequations">Extension to the compressible Navier-Stokes equations</a>
        <li><a href="#BlockGaussSeidellikepreconditioners">Block Gauss-Seidel-like preconditioners</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-76/doc/intro.dox



 <br> 

<i>
This program was contributed by Martin Kronbichler, Peter Munch, and David
Schneider. Many of the features shown here have been added to deal.II during
and for the development of the deal.II-based, efficient, matrix-free
finite-element library for high-dimensional partial differential equations
hyper.deal (see https://github.com/hyperdeal/hyperdeal). For more details and
for applications of the presented features in slightly different contexts
(high-dimensional advection equation and Vlasov-Poisson equations) see the release
paper @cite munch2020hyperdeal.


This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA) and
by the Bavarian government through the project "High-order matrix-free finite
element implementations with hybrid parallelization and improved data locality"
within the KONWIHR program.
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序求解流体力学的欧拉方程，使用显式时间积分器和无矩阵框架应用于空间的高阶非连续Galerkin离散化。这里使用的数值方法与step-67中使用的相同，但是，我们利用不同的高级MatrixFree技术来达到更高的吞吐量。

本教程的两个主要特点是。

- 使用MPI-3.0中的共享内存特性和

- 使用以单元为中心的循环，它只允许向全局向量写入一次，因此，是使用共享内存的理想选择。

我们在本教程中讨论的其他主题是模板参数VectorizedArrayType的用法和好处（而不是简单地使用VectorizedArray<Number>），以及向MatrixFree循环传递lambdas的可能性。

关于数字的细节，我们参考步骤-67的文件。我们在这里只集中讨论关键的差异。

<a name="SharedmemoryandhybridparallelizationwithMPI30"></a><h3>Shared-memory and hybrid parallelization with MPI-3.0</h3>


<a name="Motivation"></a><h4>Motivation</h4>


存在许多基于线程的共享内存库，如TBB、OpenMP或TaskFlow。将这些库集成到现有的MPI程序中，就可以使用共享内存。然而，这些库对程序员来说有一定的开销，因为所有可并行的代码部分都必须根据所使用的库进行查找和转换，包括当一些第三方数值库，如迭代求解器包，只依赖MPI时的困难。

考虑到一个纯粹的MPI并行化的有限元应用，我们可以发现，使用共享内存的主要时间和内存优势来自于访问同一计算节点上的进程所拥有的解决方案矢量的部分，而不需要进行明确的复制和缓冲。因此，MPI-3.0提供了基于所谓窗口的共享内存功能，进程可以直接访问同一共享内存域中的邻居的数据。

<a name="BasicMPI30commands"></a><h4>Basic MPI-3.0 commands</h4>


有几个相关的MPI-3.0命令值得详细讨论。一个新的MPI通信器 <code>comm_sm</code>  ，由通信器 <code>comm</code> 的进程组成，这些进程可以访问相同的共享内存，可以通过以下方式创建。

@code
MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &comm_sm);
@endcode



下面的代码片断显示了共享内存的简化分配例程，包括值类型  <code>T</code>  和大小  <code>local_size</code>  ，以及如何查询属于同一共享内存域的进程的数据指针。

@code
MPI_Win          win;         // window
T *              data_this;   // pointer to locally-owned data
std::vector<T *> data_others; // pointers to shared data


// configure shared memory
MPI_Info info;
MPI_Info_create(&info);
MPI_Info_set(info, "alloc_shared_noncontig", "true");


// allocate shared memory
MPI_Win_allocate_shared(local_size * sizeof(T), sizeof(T), info, comm_sm, &data_this, &win);


// get pointers to the shared data owned by the processes in the same sm domain
data_others.resize(size_sm);
int disp_unit = 0; // displacement size - an output parameter we don't need right now
MPI_Aint ssize = 0; // window size - an output parameter we don't need right now
for (int i = 0; i < size_sm; ++i)
  MPI_Win_shared_query(win, i, &ssize, &disp_unit, &data_others[i]);


Assert(data_this == data_others[rank_sm], ExcMessage("Something went wrong!"));
@endcode



一旦不再需要这些数据，窗口就必须被释放，这也释放了本地拥有的数据。

@code
MPI_Win_free(&win);
@endcode



<a name="MPI30andLinearAlgebradistributedVector"></a><h4>MPI-3.0 and LinearAlgebra::distributed::Vector</h4>


上一节提到的命令被整合到了 LinearAlgebra::distributed::Vector 中，如果为reinit()函数提供了可选的（第二个）通信器，那么这些命令就被用来分配共享内存。

例如，可以用一个分区器（包含全局通信器）和一个子通信器（包含同一计算节点上的进程）来设置一个向量。

@code
vec.reinit(partitioner, comm_sm);
@endcode



本地拥有的值和幽灵值可以像往常一样被处理。然而，现在用户也可以通过该函数读取共享内存邻居的值。

@code
const std::vector<ArrayView<const Number>> &
LinearAlgebra::distributed::Vector::shared_vector_data() const;
@endcode



<a name="MPI30andMatrixFree"></a><h4>MPI-3.0 and MatrixFree</h4>


虽然 LinearAlgebra::distributed::Vector 提供了分配共享内存和以协调的方式访问相邻进程的共享内存值的选项，但它实际上并没有利用共享内存的使用本身的好处。

然而，MatrixFree的基础设施确实如此。

- 一方面，在无矩阵循环 MatrixFree::loop(),  MatrixFree::cell_loop(), 和 MatrixFree::loop_cell_centric(), 中，只有需要更新的幽灵值 <em> 被 </em> 更新。来自共享内存邻居的幽灵值可以被直接访问，这使得缓冲，即把值复制到矢量的幽灵区域可能是多余的。   为了处理可能的竞赛条件，在MatrixFree中进行了必要的同步。在数值必须被缓冲的情况下，数值被直接从邻近的共享内存进程中复制，绕过了基于  <code>MPI_ISend</code>  和  <code>MPI_IRecv</code>  的更昂贵的MPI操作。

- 另一方面，像FEEvaluation和FEFaceEvaluation这样的类可以直接从共享内存中读取，所以在某些情况下确实不需要缓冲值。

为了能够使用MatrixFree的共享内存功能，MatrixFree必须通过提供用户创建的子通信器进行适当的配置。

@code
typename MatrixFree<dim, Number>::AdditionalData additional_data;


// set flags as usual (not shown)


additional_data.communicator_sm = comm_sm;


data.reinit(mapping, dof_handler, constraint, quadrature, additional_data);
@endcode






<a name="Cellcentricloops"></a><h3>Cell-centric loops</h3>


<a name="MotivationFCLvsCCL"></a><h4>Motivation: FCL vs. CCL</h4>


"以面为中心的循环"（简称FCL）在单独的循环中访问单元和面（内部和边界的）。因此，每个实体只被访问一次，单元之间的通量只被评估一次。如何在 MatrixFree::loop() 的帮助下，通过提供三个函数（一个用于细胞积分，一个用于内部，一个用于边界面）来执行以面为中心的循环，已经在步骤59和步骤67中提出。

与此相反，"以单元为中心的循环"（在hyper.deal发布的论文中简称CCL或ECL（代表以元素为中心的循环），处理一个单元并直接连续处理其所有面（即访问所有面两次）。在文献 @cite KronbichlerKormann2019 中，它们的好处对于现代CPU处理器架构来说已经很清楚了，尽管这种循环意味着通量必须被计算两次（对于内部面的每一面）。CCL有两个主要优点。

- 一方面，在CCL的情况下，解向量中的条目正好被写回主内存一次，而在FCL的情况下，尽管高速缓存有效地调度了单元和面环，但由于高速缓存容量的错过，至少有一次。

- 另一方面，由于解向量的每个条目只被访问一次，在CCL的情况下，访问解向量时不需要线程间的同步。在写入目标向量的过程中不存在竞赛条件，这使得CCL特别适用于共享内存并行化。

我们还应该注意到，尽管在CCL的情况下通量被计算了两次，但这并不自动转化为计算量的翻倍，因为已经内插到单元正交点的值可以用简单的一维内插法内插到一个面上。

<a name="CellcentricloopsandMatrixFree"></a><h4>Cell-centric loops and MatrixFree</h4>


对于以单元为中心的循环实现，可以使用函数 MatrixFree::loop_cell_centric() ，用户可以向其传递一个应该在每个单元上执行的函数。

为了得出一个适当的函数，可以在 MatrixFree::loop_cell_centric(), 中传递，原则上可以转换/合并以下三个函数，它们可以传递给 MatrixFree::loop(): 。

@code
matrix_free.template loop<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    // operation performed on cells


    FEEvaluation<dim, degree, degree + 1, 1, Number> phi(data);
    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, cell_evaluation_flags);


        // some operations on the cell quadrature points


        phi.integrate_scatter(cell_evaluation_flags, dst);
      }
  },
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    // operation performed inner faces


    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_m(data, /*is_interior_face=*/true);
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_p(data, /*is_interior_face=*/false);


    for (unsigned int face = range.first; face < range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, face_evaluation_flags);
        phi_p.reinit(face);
        phi_p.gather_evaluate(src, face_evaluation_flags);


        // some operations on the face quadrature points


        phi_m.integrate_scatter(face_evaluation_flags, dst);
        phi_p.integrate_scatter(face_evaluation_flags, dst);
      }
  },
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    // operation performed boundary faces


    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_m(data, /*is_interior_face=*/true);


    for (unsigned int face = range.first; face < range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, face_evaluation_flags);


        // some operations on the face quadrature points


        phi_m.integrate_scatter(face_evaluation_flags, dst);
      }
  },
  dst,
  src);
@endcode



以下列方式进行。

@code
matrix_free.template loop_cell_centric<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    FEEvaluation<dim, degree, degree + 1, 1, Number>     phi(data);
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_m(data, /*is_interior_face=*/true);
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_p(data, /*is_interior_face=*/false);


    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, cell_evaluation_flags);


        // some operations on the cell quadrature points


        phi.integrate_scatter(cell_evaluation_flags, dst);


        // loop over all faces of cell
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          {
            if (data.get_faces_by_cells_boundary_id(cell, face)[0] ==
                numbers::internal_face_boundary_id)
              {
                // internal face
                phi_m.reinit(cell, face);
                phi_m.gather_evaluate(src, face_evaluation_flags);
                phi_p.reinit(cell, face);
                phi_p.gather_evaluate(src, face_evaluation_flags);


                // some operations on the face quadrature points


                phi_m.integrate_scatter(face_evaluation_flags, dst);
              }
            else
              {
                // boundary face
                phi_m.reinit(cell, face);
                phi_m.gather_evaluate(src, face_evaluation_flags);


                // some operations on the face quadrature points


                phi_m.integrate_scatter(face_evaluation_flags, dst);
              }
          }
      }
  },
  dst,
  src);
@endcode



应该注意的是，FEFaceEvaluation现在是用两个数字初始化的，即单元号和本地面孔号。给出的例子只是强调了如何将以面为中心的循环转化为以单元为中心的循环，而且绝非高效，因为数据要从全局向量中多次读写，而且计算也是重复进行的。下面，我们将讨论针对这些问题的高级技术。

为了能够使用 MatrixFree::loop_cell_centric(), ，必须启用 MatrixFree::AdditionalData 的下列标志。

@code
typename MatrixFree<dim, Number>::AdditionalData additional_data;


// set flags as usual (not shown)


additional_data.hold_all_faces_to_owned_cells       = true;
additional_data.mapping_update_flags_faces_by_cells =
  additional_data.mapping_update_flags_inner_faces |
  additional_data.mapping_update_flags_boundary_faces;


data.reinit(mapping, dof_handler, constraint, quadrature, additional_data);
@endcode



特别是，这些标志使内部数据结构能够为所有单元格的面设置。

目前，deal.II中以单元为中心的循环只适用于均匀细化的网格，并且不应用任何约束条件（这是通常使用的DG的标准情况）。




<a name="ProvidinglambdastoMatrixFreeloops"></a><h3>Providing lambdas to MatrixFree loops</h3>


上面给出的例子已经使用了lambdas，它已经被提供给无矩阵循环。下面的简短例子介绍了如何在使用类和指向其方法之一的指针的版本和利用lambdas的变体之间转换函数。

在下面的代码中，一个类和它的一个方法的指针被传递给了 MatrixFree::loop(): ，这个方法应该被解释为单元格积分。

@code
void
local_apply_cell(const MatrixFree<dim, Number> &              data,
                 VectorType &                                 dst,
                 const VectorType &                           src,
                 const std::pair<unsigned int, unsigned int> &range) const
{
  FEEvaluation<dim, degree, degree + 1, 1, Number> phi(data);
  for (unsigned int cell = range.first; cell < range.second; ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(src, cell_evaluation_flags);


      // some operations on the quadrature points


      phi.integrate_scatter(cell_evaluation_flags, dst);
    }
}
@endcode



@code
matrix_free.cell_loop(&Operator::local_apply_cell, this, dst, src);
@endcode



然而，也可以通过lambda函数传递匿名函数，结果是一样的。

@code
matrix_free.template cell_loop<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    FEEvaluation<dim, degree, degree + 1, 1, Number> phi(data);
    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, cell_evaluation_flags);


        // some operations on the quadrature points


        phi.integrate_scatter(cell_evaluation_flags, dst);
      }
  },
  dst,
  src);
@endcode



<a name="VectorizedArrayType"></a><h3>VectorizedArrayType</h3>


VectorizedArray<Number>类是实现deal.II中无矩阵算法的高节点级性能的一个关键组件。它是一个围绕Number类型的 $n$ 条目的短向量的包装类，并通过内在函数将算术操作映射到适当的单指令/多数据（SIMD）概念。向量的长度可以通过 VectorizedArray::size() 查询，其底层数字类型可以通过 VectorizedArray::value_type. 查询。

在默认情况下（ <code>VectorizedArray<Number></code> ），向量长度在库的编译时被设置为与给定的处理器架构所支持的最高值相匹配。然而，也可以指定第二个可选的模板参数，如 <code>VectorizedArray<Number, size></code>, where <code>size</code> 明确控制特定指令集能力范围内的向量长度。下表列出了支持的向量长度的完整列表。

 <table align="center" class="doxtable">
  <tr>
   <th>double</th>
   <th>float</th>
   <th>ISA</th>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 1></code></td>
   <td><code>VectorizedArray<float, 1></code></td>
   <td>(auto-vectorization)</td>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 2></code></td>
   <td><code>VectorizedArray<float, 4></code></td>
   <td>SSE2/AltiVec</td>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 4></code></td>
   <td><code>VectorizedArray<float, 8></code></td>
   <td>AVX/AVX2</td>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 8></code></td>
   <td><code>VectorizedArray<float, 16></code></td>
   <td>AVX-512</td>
  </tr>
</table> 

这允许用户选择矢量长度/ISA，因此，在无矩阵算子评估中一次处理的单元数，可能会减少对缓存的压力，这对非常高的度数（和尺寸）来说是一个严重的问题。

减少填充通道数量的另一个可能的原因是为了简化调试：不用看例如8个单元，而是可以集中在一个单元上。

VectorizedArray的接口也能够被任何具有匹配接口的类型所替代。具体来说，这为deal.II准备了 <code>std::simd</code> 类，它计划成为C++23标准的一部分。下表比较了deal.II特定的SIMD类和相应的C++23类。


 <table align="center" class="doxtable">
  <tr>
   <th>VectorizedArray (deal.II)</th>
   <th>std::simd (C++23)</th>
  </tr>
  <tr>
   <td><code>VectorizedArray<Number></code></td>
   <td><code>std::experimental::native_simd<Number></code></td>
  </tr>
  <tr>
   <td><code>VectorizedArray<Number, size></code></td>
   <td><code>std::experimental::fixed_size_simd<Number, size></code></td>
  </tr>
</table> 


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Parametersandutilityfunctions"></a> 
 * <h3>Parameters and utility functions</h3>
 * 

 * 
 * 包括与 step-67 中相同的内容。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/time_stepping.h> 
 * #include <deal.II/base/timer.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/vectorization.h> 
 * 
 * #include <deal.II/distributed/tria.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * 
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_system.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/tria_accessor.h> 
 * #include <deal.II/grid/tria_iterator.h> 
 * 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/la_parallel_vector.h> 
 * 
 * #include <deal.II/matrix_free/fe_evaluation.h> 
 * #include <deal.II/matrix_free/matrix_free.h> 
 * #include <deal.II/matrix_free/operators.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <fstream> 
 * #include <iomanip> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 一个新的包含，用于根据其边界ID对单元格进行分类。
 * 

 * 
 * 
 * @code
 * #include <deal.II/matrix_free/tools.h> 
 * 
 * namespace Euler_DG 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 与  step-67  中的输入参数相同。
 * 

 * 
 * 
 * @code
 *   constexpr unsigned int testcase             = 1; 
 *   constexpr unsigned int dimension            = 2; 
 *   constexpr unsigned int n_global_refinements = 2; 
 *   constexpr unsigned int fe_degree            = 5; 
 *   constexpr unsigned int n_q_points_1d        = fe_degree + 2; 
 * 
 * @endcode
 * 
 * 这个参数指定共享内存组的大小。目前，只有1和 numbers::invalid_unsigned_int 的值是可能的，这导致了内存功能可以被关闭，或者所有访问同一共享内存域的进程被分组。
 * 

 * 
 * 
 * @code
 *   constexpr unsigned int group_size = numbers::invalid_unsigned_int; 
 * 
 *   using Number = double; 
 * 
 * @endcode
 * 
 * 这里，数据结构的类型被选择为矢量化。在默认情况下，使用VectorizedArray<Number>，也就是说，在给定的硬件上使用最高的指令集架构扩展，具有最大数量的向量通道。然而，人们可能会减少填充通道的数量，例如，通过编写 <code>using VectorizedArrayType = VectorizedArray<Number, 4></code> ，只处理4个单元。
 * 

 * 
 * 
 * @code
 *   using VectorizedArrayType = VectorizedArray<Number>; 
 * 
 * @endcode
 * 
 * 以下参数没有改变。
 * 

 * 
 * 
 * @code
 *   constexpr double gamma       = 1.4; 
 *   constexpr double final_time  = testcase == 0 ? 10 : 2.0; 
 *   constexpr double output_tick = testcase == 0 ? 1 : 0.05; 
 * 
 *   const double courant_number = 0.15 / std::pow(fe_degree, 1.5); 
 * 
 * @endcode
 * 
 * 指定对性能研究有用的最大时间步骤数。
 * 

 * 
 * 
 * @code
 *   constexpr unsigned int max_time_steps = numbers::invalid_unsigned_int; 
 * 
 * @endcode
 * 
 * 与Runge-Kutta有关的函数从 step-67 复制过来，并稍作修改，以尽量减少全局矢量访问。
 * 

 * 
 * 
 * @code
 *   enum LowStorageRungeKuttaScheme 
 *   { 
 *     stage_3_order_3, 
 *     stage_5_order_4, 
 *     stage_7_order_4, 
 *     stage_9_order_5, 
 *   }; 
 *   constexpr LowStorageRungeKuttaScheme lsrk_scheme = stage_5_order_4; 
 * 
 *   class LowStorageRungeKuttaIntegrator 
 *   { 
 *   public: 
 *     LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme) 
 *     { 
 *       TimeStepping::runge_kutta_method lsrk; 
 *       switch (scheme) 
 *         { 
 *           case stage_3_order_3: 
 *             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3; 
 *             break; 
 *           case stage_5_order_4: 
 *             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4; 
 *             break; 
 *           case stage_7_order_4: 
 *             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4; 
 *             break; 
 *           case stage_9_order_5: 
 *             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5; 
 *             break; 
 * 
 *           default: 
 *             AssertThrow(false, ExcNotImplemented()); 
 *         } 
 *       TimeStepping::LowStorageRungeKutta< 
 *         LinearAlgebra::distributed::Vector<Number>> 
 *                           rk_integrator(lsrk); 
 *       std::vector<double> ci; // not used 
 *       rk_integrator.get_coefficients(ai, bi, ci); 
 *     } 
 * 
 *     unsigned int n_stages() const 
 *     { 
 *       return bi.size(); 
 *     } 
 * 
 *     template <typename VectorType, typename Operator> 
 *     void perform_time_step(const Operator &pde_operator, 
 *                            const double    current_time, 
 *                            const double    time_step, 
 *                            VectorType &    solution, 
 *                            VectorType &    vec_ri, 
 *                            VectorType &    vec_ki) const 
 *     { 
 *       AssertDimension(ai.size() + 1, bi.size()); 
 * 
 *       vec_ki.swap(solution); 
 * 
 *       double sum_previous_bi = 0; 
 *       for (unsigned int stage = 0; stage < bi.size(); ++stage) 
 *         { 
 *           const double c_i = stage == 0 ? 0 : sum_previous_bi + ai[stage - 1]; 
 * 
 *           pde_operator.perform_stage(stage, 
 *                                      current_time + c_i * time_step, 
 *                                      bi[stage] * time_step, 
 *                                      (stage == bi.size() - 1 ? 
 *                                         0 : 
 *                                         ai[stage] * time_step), 
 *                                      (stage % 2 == 0 ? vec_ki : vec_ri), 
 *                                      (stage % 2 == 0 ? vec_ri : vec_ki), 
 *                                      solution); 
 * 
 *           if (stage > 0) 
 *             sum_previous_bi += bi[stage - 1]; 
 *         } 
 *     } 
 * 
 *   private: 
 *     std::vector<double> bi; 
 *     std::vector<double> ai; 
 *   }; 
 * 
 * @endcode
 * 
 * 来自  step-67  的欧拉特定实用函数。
 * 

 * 
 * 
 * @code
 *   enum EulerNumericalFlux 
 *   { 
 *     lax_friedrichs_modified, 
 *     harten_lax_vanleer, 
 *   }; 
 *   constexpr EulerNumericalFlux numerical_flux_type = lax_friedrichs_modified; 
 * 
 *   template <int dim> 
 *   class ExactSolution : public Function<dim> 
 *   { 
 *   public: 
 *     ExactSolution(const double time) 
 *       : Function<dim>(dim + 2, time) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double ExactSolution<dim>::value(const Point<dim> & x, 
 *                                    const unsigned int component) const 
 *   { 
 *     const double t = this->get_time(); 
 * 
 *     switch (testcase) 
 *       { 
 *         case 0: 
 *           { 
 *             Assert(dim == 2, ExcNotImplemented()); 
 *             const double beta = 5; 
 * 
 *             Point<dim> x0; 
 *             x0[0] = 5.; 
 *             const double radius_sqr = 
 *               (x - x0).norm_square() - 2. * (x[0] - x0[0]) * t + t * t; 
 *             const double factor = 
 *               beta / (numbers::PI * 2) * std::exp(1. - radius_sqr); 
 *             const double density_log = std::log2( 
 *               std::abs(1. - (gamma - 1.) / gamma * 0.25 * factor * factor)); 
 *             const double density = std::exp2(density_log * (1. / (gamma - 1.))); 
 *             const double u       = 1. - factor * (x[1] - x0[1]); 
 *             const double v       = factor * (x[0] - t - x0[0]); 
 * 
 *             if (component == 0) 
 *               return density; 
 *             else if (component == 1) 
 *               return density * u; 
 *             else if (component == 2) 
 *               return density * v; 
 *             else 
 *               { 
 *                 const double pressure = 
 *                   std::exp2(density_log * (gamma / (gamma - 1.))); 
 *                 return pressure / (gamma - 1.) + 
 *                        0.5 * (density * u * u + density * v * v); 
 *               } 
 *           } 
 * 
 *         case 1: 
 *           { 
 *             if (component == 0) 
 *               return 1.; 
 *             else if (component == 1) 
 *               return 0.4; 
 *             else if (component == dim + 1) 
 *               return 3.097857142857143; 
 *             else 
 *               return 0.; 
 *           } 
 * 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *           return 0.; 
 *       } 
 *   } 
 * 
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, dim, Number> 
 *     euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables) 
 *   { 
 *     const Number inverse_density = Number(1.) / conserved_variables[0]; 
 * 
 *     Tensor<1, dim, Number> velocity; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       velocity[d] = conserved_variables[1 + d] * inverse_density; 
 * 
 *     return velocity; 
 *   } 
 * 
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Number 
 *     euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables) 
 *   { 
 *     const Tensor<1, dim, Number> velocity = 
 *       euler_velocity<dim>(conserved_variables); 
 * 
 *     Number rho_u_dot_u = conserved_variables[1] * velocity[0]; 
 *     for (unsigned int d = 1; d < dim; ++d) 
 *       rho_u_dot_u += conserved_variables[1 + d] * velocity[d]; 
 * 
 *     return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u); 
 *   } 
 * 
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>> 
 *     euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables) 
 *   { 
 *     const Tensor<1, dim, Number> velocity = 
 *       euler_velocity<dim>(conserved_variables); 
 *     const Number pressure = euler_pressure<dim>(conserved_variables); 
 * 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>> flux; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       { 
 *         flux[0][d] = conserved_variables[1 + d]; 
 *         for (unsigned int e = 0; e < dim; ++e) 
 *           flux[e + 1][d] = conserved_variables[e + 1] * velocity[d]; 
 *         flux[d + 1][d] += pressure; 
 *         flux[dim + 1][d] = 
 *           velocity[d] * (conserved_variables[dim + 1] + pressure); 
 *       } 
 * 
 *     return flux; 
 *   } 
 * 
 *   template <int n_components, int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, n_components, Number> 
 *     operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix, 
 *               const Tensor<1, dim, Number> &                         vector) 
 *   { 
 *     Tensor<1, n_components, Number> result; 
 *     for (unsigned int d = 0; d < n_components; ++d) 
 *       result[d] = matrix[d] * vector; 
 *     return result; 
 *   } 
 * 
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, dim + 2, Number> 
 *     euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m, 
 *                          const Tensor<1, dim + 2, Number> &u_p, 
 *                          const Tensor<1, dim, Number> &    normal) 
 *   { 
 *     const auto velocity_m = euler_velocity<dim>(u_m); 
 *     const auto velocity_p = euler_velocity<dim>(u_p); 
 * 
 *     const auto pressure_m = euler_pressure<dim>(u_m); 
 *     const auto pressure_p = euler_pressure<dim>(u_p); 
 * 
 *     const auto flux_m = euler_flux<dim>(u_m); 
 *     const auto flux_p = euler_flux<dim>(u_p); 
 * 
 *     switch (numerical_flux_type) 
 *       { 
 *         case lax_friedrichs_modified: 
 *           { 
 *             const auto lambda = 
 *               0.5 * std::sqrt(std::max(velocity_p.norm_square() + 
 *                                          gamma * pressure_p * (1. / u_p[0]), 
 *                                        velocity_m.norm_square() + 
 *                                          gamma * pressure_m * (1. / u_m[0]))); 
 * 
 *             return 0.5 * (flux_m * normal + flux_p * normal) + 
 *                    0.5 * lambda * (u_m - u_p); 
 *           } 
 * 
 *         case harten_lax_vanleer: 
 *           { 
 *             const auto avg_velocity_normal = 
 *               0.5 * ((velocity_m + velocity_p) * normal); 
 *             const auto   avg_c = std::sqrt(std::abs( 
 *               0.5 * gamma * 
 *               (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0])))); 
 *             const Number s_pos = 
 *               std::max(Number(), avg_velocity_normal + avg_c); 
 *             const Number s_neg = 
 *               std::min(Number(), avg_velocity_normal - avg_c); 
 *             const Number inverse_s = Number(1.) / (s_pos - s_neg); 
 * 
 *             return inverse_s * 
 *                    ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) - 
 *                     s_pos * s_neg * (u_m - u_p)); 
 *           } 
 * 
 *         default: 
 *           { 
 *             Assert(false, ExcNotImplemented()); 
 *             return {}; 
 *           } 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 来自  step-67  的通用实用函数。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename VectorizedArrayType> 
 *   VectorizedArrayType 
 *   evaluate_function(const Function<dim> &                  function, 
 *                     const Point<dim, VectorizedArrayType> &p_vectorized, 
 *                     const unsigned int                     component) 
 *   { 
 *     VectorizedArrayType result; 
 *     for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v) 
 *       { 
 *         Point<dim> p; 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           p[d] = p_vectorized[d][v]; 
 *         result[v] = function.value(p, component); 
 *       } 
 *     return result; 
 *   } 
 * 
 *   template <int dim, typename VectorizedArrayType, int n_components = dim + 2> 
 *   Tensor<1, n_components, VectorizedArrayType> 
 *   evaluate_function(const Function<dim> &                  function, 
 *                     const Point<dim, VectorizedArrayType> &p_vectorized) 
 *   { 
 *     AssertDimension(function.n_components, n_components); 
 *     Tensor<1, n_components, VectorizedArrayType> result; 
 *     for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v) 
 *       { 
 *         Point<dim> p; 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           p[d] = p_vectorized[d][v]; 
 *         for (unsigned int d = 0; d < n_components; ++d) 
 *           result[d][v] = function.value(p, d); 
 *       } 
 *     return result; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="EuleroperatorusingacellcentricloopandMPI30sharedmemory"></a> 
 * <h3>Euler operator using a cell-centric loop and MPI-3.0 shared memory</h3>
 * 

 * 
 * 来自 step-67 的欧拉算子，有一些变化，详见下文。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   class EulerOperator 
 *   { 
 *   public: 
 *     static constexpr unsigned int n_quadrature_points_1d = n_points_1d; 
 * 
 *     EulerOperator(TimerOutput &timer_output); 
 * 
 *     ~EulerOperator(); 
 * 
 *     void reinit(const Mapping<dim> &   mapping, 
 *                 const DoFHandler<dim> &dof_handler); 
 * 
 *     void set_inflow_boundary(const types::boundary_id       boundary_id, 
 *                              std::unique_ptr<Function<dim>> inflow_function); 
 * 
 *     void set_subsonic_outflow_boundary( 
 *       const types::boundary_id       boundary_id, 
 *       std::unique_ptr<Function<dim>> outflow_energy); 
 * 
 *     void set_wall_boundary(const types::boundary_id boundary_id); 
 * 
 *     void set_body_force(std::unique_ptr<Function<dim>> body_force); 
 * 
 *     void 
 *     perform_stage(const unsigned int                                stage, 
 *                   const Number                                      cur_time, 
 *                   const Number                                      bi, 
 *                   const Number                                      ai, 
 *                   const LinearAlgebra::distributed::Vector<Number> &current_ri, 
 *                   LinearAlgebra::distributed::Vector<Number> &      vec_ki, 
 *                   LinearAlgebra::distributed::Vector<Number> &solution) const; 
 * 
 *     void project(const Function<dim> &                       function, 
 *                  LinearAlgebra::distributed::Vector<Number> &solution) const; 
 * 
 *     std::array<double, 3> compute_errors( 
 *       const Function<dim> &                             function, 
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const; 
 * 
 *     double compute_cell_transport_speed( 
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const; 
 * 
 *     void 
 *     initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const; 
 * 
 *   private: 
 * 
 * @endcode
 * 
 * 包含子通信器的SubCommunicatorWrapper实例，我们需要将其传递给 MatrixFree::reinit() ，以便能够利用MPI-3.0的共享内存功能。
 * 

 * 
 * 
 * @code
 *     MPI_Comm subcommunicator; 
 * 
 *     MatrixFree<dim, Number, VectorizedArrayType> data; 
 * 
 *     TimerOutput &timer; 
 * 
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>> 
 *       inflow_boundaries; 
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>> 
 *                                    subsonic_outflow_boundaries; 
 *     std::set<types::boundary_id>   wall_boundaries; 
 *     std::unique_ptr<Function<dim>> body_force; 
 *   }; 
 * 
 * @endcode
 * 
 * 新的构造函数，可以创建一个子通信器。用户可以通过全局参数group_size指定子通信器的大小。如果大小被设置为-1，一个共享内存域的所有MPI进程将被合并为一个组。指定的大小对于MatrixFree的共享内存能力的好处是决定性的，因此，设置为 <code>size</code> to <code>-1</code> 是一个合理的选择。通过设置 <code>1</code> ，用户明确地禁用了MatrixFree的MPI-3.0共享内存功能，而完全依赖MPI-2.0功能，如 <code>MPI_Isend</code> 和 <code>MPI_Irecv</code>  。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   EulerOperator<dim, degree, n_points_1d>::EulerOperator(TimerOutput &timer) 
 *     : timer(timer) 
 *   { 
 * #if DEAL_II_MPI_VERSION_GTE(3, 0) 
 *     if (group_size == 1) 
 *       { 
 *         this->subcommunicator = MPI_COMM_SELF; 
 *       } 
 *     else if (group_size == numbers::invalid_unsigned_int) 
 *       { 
 *         const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD); 
 * 
 *         MPI_Comm_split_type(MPI_COMM_WORLD, 
 *                             MPI_COMM_TYPE_SHARED, 
 *                             rank, 
 *                             MPI_INFO_NULL, 
 *                             &subcommunicator); 
 *       } 
 *     else 
 *       { 
 *         Assert(false, ExcNotImplemented()); 
 *       } 
 * #else 
 *     (void)subcommunicator; 
 *     (void)group_size; 
 *     this->subcommunicator = MPI_COMM_SELF; 
 * #endif 
 *   } 
 * 
 * @endcode
 * 
 * 新增负责释放子通信器的析构器。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   EulerOperator<dim, degree, n_points_1d>::~EulerOperator() 
 *   { 
 * #ifdef DEAL_II_WITH_MPI 
 *     if (this->subcommunicator != MPI_COMM_SELF) 
 *       MPI_Comm_free(&subcommunicator); 
 * #endif 
 *   } 
 * 
 * @endcode
 * 
 * 修改了 reinit() 函数，以设置 MatrixFree 中的内部数据结构，使其可以被以单元为中心的循环使用，并使用 MPI-3.0 的共享内存功能。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::reinit( 
 *     const Mapping<dim> &   mapping, 
 *     const DoFHandler<dim> &dof_handler) 
 *   { 
 *     const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler}; 
 *     const AffineConstraints<double>            dummy; 
 *     const std::vector<const AffineConstraints<double> *> constraints = {&dummy}; 
 *     const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d), 
 *                                                     QGauss<1>(fe_degree + 1)}; 
 * 
 *     typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData 
 *       additional_data; 
 *     additional_data.mapping_update_flags = 
 *       (update_gradients | update_JxW_values | update_quadrature_points | 
 *        update_values); 
 *     additional_data.mapping_update_flags_inner_faces = 
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors | 
 *        update_values); 
 *     additional_data.mapping_update_flags_boundary_faces = 
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors | 
 *        update_values); 
 *     additional_data.tasks_parallel_scheme = 
 *       MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData::none; 
 * 
 * @endcode
 * 
 * 对单元格进行分类，使所有车道的每个面都有相同的边界ID。这在严格意义上是没有必要的，但是，可以在 EulerOperator::perform_stage() 中写出更简单的代码，不需要屏蔽，因为可以保证所有分组的单元格（在一个VectorizedArray中）也要对面进行完全相同的操作。
 * 

 * 
 * 
 * @code
 *     MatrixFreeTools::categorize_by_boundary_ids(dof_handler.get_triangulation(), 
 *                                                 additional_data); 
 * 
 * @endcode
 * 
 * 通过提供子通信器在MatrixFree中启用MPI-3.0共享内存功能。
 * 

 * 
 * 
 * @code
 *     additional_data.communicator_sm = subcommunicator; 
 * 
 *     data.reinit( 
 *       mapping, dof_handlers, constraints, quadratures, additional_data); 
 *   } 
 * 
 * @endcode
 * 
 * 下面这个函数是做Runge--Kutta更新的整个阶段，并且是
 * 

 * 
 * - 旁边的设置稍作修改
 * 

 * 
 * 与 step-67 相比，我们没有依次执行平流步骤（使用 MatrixFree::loop()) 和反质量矩阵步骤（使用 MatrixFree::cell_loop()) ），而是在 MatrixFree::loop_cell_centric(). 中一次性评估所有内容。 这个函数期望在每个本地拥有的（宏）单元上执行一个单独的函数作为参数，这样我们就需要在该单元的所有面上循环，自行执行需要的积分步骤。
 * 

 * 
 * 以下函数在很大程度上包含了 step-67 中的以下函数的副本，所以这里跳过了与弱形式的评估有关的评论。
 * 

 * 
 * -  <code>EulerDG::EulerOperator::local_apply_cell</code>  
 * 

 * 
 * -  <code>EulerDG::EulerOperator::local_apply_face</code>  
 * 

 * 
 * -  <code>EulerDG::EulerOperator::local_apply_boundary_face</code>  
 * 

 * 
 * -  <code>EulerDG::EulerOperator::local_apply_inverse_mass_matrix</code>  
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::perform_stage( 
 *     const unsigned int                                stage, 
 *     const Number                                      current_time, 
 *     const Number                                      bi, 
 *     const Number                                      ai, 
 *     const LinearAlgebra::distributed::Vector<Number> &current_ri, 
 *     LinearAlgebra::distributed::Vector<Number> &      vec_ki, 
 *     LinearAlgebra::distributed::Vector<Number> &      solution) const 
 *   { 
 *     for (auto &i : inflow_boundaries) 
 *       i.second->set_time(current_time); 
 *     for (auto &i : subsonic_outflow_boundaries) 
 *       i.second->set_time(current_time); 
 * 
 * @endcode
 * 
 * 通过调用 MatrixFree::loop_cell_centric() 并提供一个包含单元、面和边界-面积分效果的lambda来运行以单元为中心的循环。
 * 

 * 
 * 
 * @code
 *     data.template loop_cell_centric<LinearAlgebra::distributed::Vector<Number>, 
 *                                     LinearAlgebra::distributed::Vector<Number>>( 
 *       [&](const auto &data, auto &dst, const auto &src, const auto cell_range) { 
 *         using FECellIntegral = FEEvaluation<dim, 
 *                                             degree, 
 *                                             n_points_1d, 
 *                                             dim + 2, 
 *                                             Number, 
 *                                             VectorizedArrayType>; 
 *         using FEFaceIntegral = FEFaceEvaluation<dim, 
 *                                                 degree, 
 *                                                 n_points_1d, 
 *                                                 dim + 2, 
 *                                                 Number, 
 *                                                 VectorizedArrayType>; 
 * 
 *         FECellIntegral phi(data); 
 *         FECellIntegral phi_temp(data); 
 *         FEFaceIntegral phi_m(data, true); 
 *         FEFaceIntegral phi_p(data, false); 
 * 
 *         Tensor<1, dim, VectorizedArrayType>     constant_body_force; 
 *         const Functions::ConstantFunction<dim> *constant_function = 
 *           dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get()); 
 * 
 *         if (constant_function) 
 *           constant_body_force = 
 *             evaluate_function<dim, VectorizedArrayType, dim>( 
 *               *constant_function, Point<dim, VectorizedArrayType>()); 
 * 
 *         const dealii::internal::EvaluatorTensorProduct< 
 *           dealii::internal::EvaluatorVariant::evaluate_evenodd, 
 *           dim, 
 *           n_points_1d, 
 *           n_points_1d, 
 *           VectorizedArrayType> 
 *           eval(AlignedVector<VectorizedArrayType>(), 
 *                data.get_shape_info().data[0].shape_gradients_collocation_eo, 
 *                AlignedVector<VectorizedArrayType>()); 
 * 
 *         AlignedVector<VectorizedArrayType> buffer(phi.static_n_q_points * 
 *                                                   phi.n_components); 
 * 
 * @endcode
 * 
 * 在所有单元格批次上循环。
 * 

 * 
 * 
 * @code
 *         for (unsigned int cell = cell_range.first; cell < cell_range.second; 
 *              ++cell) 
 *           { 
 *             phi.reinit(cell); 
 * 
 *             if (ai != Number()) 
 *               phi_temp.reinit(cell); 
 * 
 * @endcode
 * 
 * 从全局矢量中读取数值并计算正交点的数值。
 * 

 * 
 * 
 * @code
 *             if (ai != Number() && stage == 0) 
 *               { 
 *                 phi.read_dof_values(src); 
 * 
 *                 for (unsigned int i = 0; 
 *                      i < phi.static_dofs_per_component * (dim + 2); 
 *                      ++i) 
 *                   phi_temp.begin_dof_values()[i] = phi.begin_dof_values()[i]; 
 * 
 *                 phi.evaluate(EvaluationFlags::values); 
 *               } 
 *             else 
 *               { 
 *                 phi.gather_evaluate(src, EvaluationFlags::values); 
 *               } 
 * 
 * @endcode
 * 
 * 缓冲正交点的计算值，因为这些值在下一步被 FEEvaluation::submit_value() 所覆盖，但是，在后面的面积分中需要。
 * 

 * 
 * 
 * @code
 *             for (unsigned int i = 0; i < phi.static_n_q_points * (dim + 2); ++i) 
 *               buffer[i] = phi.begin_values()[i]; 
 * 
 * @endcode
 * 
 * 在单元格正交点上应用单元格积分。也可参见来自  step-67  的函数  <code>EulerOperator::local_apply_cell()</code>  。
 * 

 * 
 * 
 * @code
 *             for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *               { 
 *                 const auto w_q = phi.get_value(q); 
 *                 phi.submit_gradient(euler_flux<dim>(w_q), q); 
 *                 if (body_force.get() != nullptr) 
 *                   { 
 *                     const Tensor<1, dim, VectorizedArrayType> force = 
 *                       constant_function ? 
 *                         constant_body_force : 
 *                         evaluate_function<dim, VectorizedArrayType, dim>( 
 *                           *body_force, phi.quadrature_point(q)); 
 * 
 *                     Tensor<1, dim + 2, VectorizedArrayType> forcing; 
 *                     for (unsigned int d = 0; d < dim; ++d) 
 *                       forcing[d + 1] = w_q[0] * force[d]; 
 *                     for (unsigned int d = 0; d < dim; ++d) 
 *                       forcing[dim + 1] += force[d] * w_q[d + 1]; 
 * 
 *                     phi.submit_value(forcing, q); 
 *                   } 
 *               } 
 * 
 * @endcode
 * 
 * 用正交点中的测试函数的梯度进行测试。我们跳过插值回到元素的支持点，因为我们首先收集单元格正交点的所有贡献，只在最后一步进行插值。
 * 

 * 
 * 
 * @code
 *             { 
 *               auto *values_ptr   = phi.begin_values(); 
 *               auto *gradient_ptr = phi.begin_gradients(); 
 * 
 *               for (unsigned int c = 0; c < dim + 2; ++c) 
 *                 { 
 *                   if (dim >= 1 && body_force.get() == nullptr) 
 *                     eval.template gradients<0, false, false>( 
 *                       gradient_ptr + phi.static_n_q_points * 0, values_ptr); 
 *                   else if (dim >= 1) 
 *                     eval.template gradients<0, false, true>( 
 *                       gradient_ptr + phi.static_n_q_points * 0, values_ptr); 
 *                   if (dim >= 2) 
 *                     eval.template gradients<1, false, true>( 
 *                       gradient_ptr + phi.static_n_q_points * 1, values_ptr); 
 *                   if (dim >= 3) 
 *                     eval.template gradients<2, false, true>( 
 *                       gradient_ptr + phi.static_n_q_points * 2, values_ptr); 
 * 
 *                   values_ptr += phi.static_n_q_points; 
 *                   gradient_ptr += phi.static_n_q_points * dim; 
 *                 } 
 *             } 
 * 
 * @endcode
 * 
 * 在当前单元格的所有面上进行循环。
 * 

 * 
 * 
 * @code
 *             for (unsigned int face = 0; 
 *                  face < GeometryInfo<dim>::faces_per_cell; 
 *                  ++face) 
 *               { 
 * 
 * @endcode
 * 
 * 确定当前面的边界ID。由于我们在设置MatrixFree时，保证了所有填充的车道都有相同的边界ID，我们可以选择第一个车道的边界ID。
 * 

 * 
 * 
 * @code
 *                 const auto boundary_ids = 
 *                   data.get_faces_by_cells_boundary_id(cell, face); 
 * 
 *                 Assert(std::equal(boundary_ids.begin(), 
 *                                   boundary_ids.begin() + 
 *                                     data.n_active_entries_per_cell_batch(cell), 
 *                                   boundary_ids.begin()), 
 *                        ExcMessage("Boundary IDs of lanes differ.")); 
 * 
 *                 const auto boundary_id = boundary_ids[0]; 
 * 
 *                 phi_m.reinit(cell, face); 
 * 
 * @endcode
 * 
 * 通过简单的一维插值，将单元格正交点的值插到当前面的正交点上。
 * 

 * 
 * 
 * @code
 *                 internal::FEFaceNormalEvaluationImpl<dim, 
 *                                                      n_points_1d - 1, 
 *                                                      VectorizedArrayType>:: 
 *                   template interpolate_quadrature<true, false>( 
 *                     dim + 2, 
 *                     data.get_shape_info(), 
 *                     buffer.data(), 
 *                     phi_m.begin_values(), 
 *                     false, 
 *                     face); 
 * 
 * @endcode
 * 
 * 检查该面是内部面还是边界面，并根据这一信息选择不同的代码路径。
 * 

 * 
 * 
 * @code
 *                 if (boundary_id == numbers::internal_face_boundary_id) 
 *                   { 
 * 
 * @endcode
 * 
 * 处理和内部面。以下几行代码是对 step-67 中 <code>EulerDG::EulerOperator::local_apply_face</code> 函数的复制。
 * 

 * 
 * 
 * @code
 *                     phi_p.reinit(cell, face); 
 *                     phi_p.gather_evaluate(src, EvaluationFlags::values); 
 * 
 *                     for (unsigned int q = 0; q < phi_m.n_q_points; ++q) 
 *                       { 
 *                         const auto numerical_flux = 
 *                           euler_numerical_flux<dim>(phi_m.get_value(q), 
 *                                                     phi_p.get_value(q), 
 *                                                     phi_m.get_normal_vector(q)); 
 *                         phi_m.submit_value(-numerical_flux, q); 
 *                       } 
 *                   } 
 *                 else 
 *                   { 
 * 
 * @endcode
 * 
 * 处理一个边界面。下面这几行代码是对 step-67 中 <code>EulerDG::EulerOperator::local_apply_boundary_face</code> 函数的复制。
 * 

 * 
 * 
 * @code
 *                     for (unsigned int q = 0; q < phi_m.n_q_points; ++q) 
 *                       { 
 *                         const auto w_m    = phi_m.get_value(q); 
 *                         const auto normal = phi_m.get_normal_vector(q); 
 * 
 *                         auto rho_u_dot_n = w_m[1] * normal[0]; 
 *                         for (unsigned int d = 1; d < dim; ++d) 
 *                           rho_u_dot_n += w_m[1 + d] * normal[d]; 
 * 
 *                         bool at_outflow = false; 
 * 
 *                         Tensor<1, dim + 2, VectorizedArrayType> w_p; 
 * 
 *                         if (wall_boundaries.find(boundary_id) != 
 *                             wall_boundaries.end()) 
 *                           { 
 *                             w_p[0] = w_m[0]; 
 *                             for (unsigned int d = 0; d < dim; ++d) 
 *                               w_p[d + 1] = 
 *                                 w_m[d + 1] - 2. * rho_u_dot_n * normal[d]; 
 *                             w_p[dim + 1] = w_m[dim + 1]; 
 *                           } 
 *                         else if (inflow_boundaries.find(boundary_id) != 
 *                                  inflow_boundaries.end()) 
 *                           w_p = evaluate_function( 
 *                             *inflow_boundaries.find(boundary_id)->second, 
 *                             phi_m.quadrature_point(q)); 
 *                         else if (subsonic_outflow_boundaries.find( 
 *                                    boundary_id) != 
 *                                  subsonic_outflow_boundaries.end()) 
 *                           { 
 *                             w_p = w_m; 
 *                             w_p[dim + 1] = 
 *                               evaluate_function(*subsonic_outflow_boundaries 
 *                                                    .find(boundary_id) 
 *                                                    ->second, 
 *                                                 phi_m.quadrature_point(q), 
 *                                                 dim + 1); 
 *                             at_outflow = true; 
 *                           } 
 *                         else 
 *                           AssertThrow(false, 
 *                                       ExcMessage( 
 *                                         "Unknown boundary id, did " 
 *                                         "you set a boundary condition for " 
 *                                         "this part of the domain boundary?")); 
 * 
 *                         auto flux = euler_numerical_flux<dim>(w_m, w_p, normal); 
 * 
 *                         if (at_outflow) 
 *                           for (unsigned int v = 0; 
 *                                v < VectorizedArrayType::size(); 
 *                                ++v) 
 *                             { 
 *                               if (rho_u_dot_n[v] < -1e-12) 
 *                                 for (unsigned int d = 0; d < dim; ++d) 
 *                                   flux[d + 1][v] = 0.; 
 *                             } 
 * 
 *                         phi_m.submit_value(-flux, q); 
 *                       } 
 *                   } 
 * 
 * @endcode
 * 
 * 通过正交评估与单元相关的局部积分，并通过简单的一维插值加入到单元贡献中。
 * 

 * 
 * 
 * @code
 *                 internal::FEFaceNormalEvaluationImpl<dim, 
 *                                                      n_points_1d - 1, 
 *                                                      VectorizedArrayType>:: 
 *                   template interpolate_quadrature<false, true>( 
 *                     dim + 2, 
 *                     data.get_shape_info(), 
 *                     phi_m.begin_values(), 
 *                     phi.begin_values(), 
 *                     false, 
 *                     face); 
 *               } 
 * 
 * @endcode
 * 
 * 在单元格正交点中应用反质量矩阵。也请参见来自  <code>EulerDG::EulerOperator::local_apply_inverse_mass_matrix()</code>  的函数  step-67  。
 * 

 * 
 * 
 * @code
 *             for (unsigned int q = 0; q < phi.static_n_q_points; ++q) 
 *               { 
 *                 const auto factor = VectorizedArrayType(1.0) / phi.JxW(q); 
 *                 for (unsigned int c = 0; c < dim + 2; ++c) 
 *                   phi.begin_values()[c * phi.static_n_q_points + q] = 
 *                     phi.begin_values()[c * phi.static_n_q_points + q] * factor; 
 *               } 
 * 
 * @endcode
 * 
 * 将数值从配位空间转换到原始高斯-洛巴托空间。
 * 

 * 
 * 
 * @code
 *             internal::FEEvaluationImplBasisChange< 
 *               dealii::internal::EvaluatorVariant::evaluate_evenodd, 
 *               internal::EvaluatorQuantity::hessian, 
 *               dim, 
 *               degree + 1, 
 *               n_points_1d, 
 *               VectorizedArrayType, 
 *               VectorizedArrayType>::do_backward(dim + 2, 
 *                                                 data.get_shape_info() 
 *                                                   .data[0] 
 *                                                   .inverse_shape_values_eo, 
 *                                                 false, 
 *                                                 phi.begin_values(), 
 *                                                 phi.begin_dof_values()); 
 * 
 * @endcode
 * 
 * 执行Runge-Kutta更新并将结果写回全局向量。
 * 

 * 
 * 
 * @code
 *             if (ai == Number()) 
 *               { 
 *                 for (unsigned int q = 0; q < phi.static_dofs_per_cell; ++q) 
 *                   phi.begin_dof_values()[q] = bi * phi.begin_dof_values()[q]; 
 *                 phi.distribute_local_to_global(solution); 
 *               } 
 *             else 
 *               { 
 *                 if (stage != 0) 
 *                   phi_temp.read_dof_values(solution); 
 * 
 *                 for (unsigned int q = 0; q < phi.static_dofs_per_cell; ++q) 
 *                   { 
 *                     const auto K_i = phi.begin_dof_values()[q]; 
 * 
 *                     phi.begin_dof_values()[q] = 
 *                       phi_temp.begin_dof_values()[q] + (ai * K_i); 
 * 
 *                     phi_temp.begin_dof_values()[q] += bi * K_i; 
 *                   } 
 *                 phi.set_dof_values(dst); 
 *                 phi_temp.set_dof_values(solution); 
 *               } 
 *           } 
 *       }, 
 *       vec_ki, 
 *       current_ri, 
 *       true, 
 *       MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::values); 
 *   } 
 * 
 * @endcode
 * 
 * 从这里开始， step-67 的代码没有改变。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::initialize_vector( 
 *     LinearAlgebra::distributed::Vector<Number> &vector) const 
 *   { 
 *     data.initialize_dof_vector(vector); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary( 
 *     const types::boundary_id       boundary_id, 
 *     std::unique_ptr<Function<dim>> inflow_function) 
 *   { 
 *     AssertThrow(subsonic_outflow_boundaries.find(boundary_id) == 
 *                     subsonic_outflow_boundaries.end() && 
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(), 
 *                 ExcMessage("You already set the boundary with id " + 
 *                            std::to_string(static_cast<int>(boundary_id)) + 
 *                            " to another type of boundary before now setting " + 
 *                            "it as inflow")); 
 *     AssertThrow(inflow_function->n_components == dim + 2, 
 *                 ExcMessage("Expected function with dim+2 components")); 
 * 
 *     inflow_boundaries[boundary_id] = std::move(inflow_function); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary( 
 *     const types::boundary_id       boundary_id, 
 *     std::unique_ptr<Function<dim>> outflow_function) 
 *   { 
 *     AssertThrow(inflow_boundaries.find(boundary_id) == 
 *                     inflow_boundaries.end() && 
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(), 
 *                 ExcMessage("You already set the boundary with id " + 
 *                            std::to_string(static_cast<int>(boundary_id)) + 
 *                            " to another type of boundary before now setting " + 
 *                            "it as subsonic outflow")); 
 *     AssertThrow(outflow_function->n_components == dim + 2, 
 *                 ExcMessage("Expected function with dim+2 components")); 
 * 
 *     subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_wall_boundary( 
 *     const types::boundary_id boundary_id) 
 *   { 
 *     AssertThrow(inflow_boundaries.find(boundary_id) == 
 *                     inflow_boundaries.end() && 
 *                   subsonic_outflow_boundaries.find(boundary_id) == 
 *                     subsonic_outflow_boundaries.end(), 
 *                 ExcMessage("You already set the boundary with id " + 
 *                            std::to_string(static_cast<int>(boundary_id)) + 
 *                            " to another type of boundary before now setting " + 
 *                            "it as wall boundary")); 
 * 
 *     wall_boundaries.insert(boundary_id); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_body_force( 
 *     std::unique_ptr<Function<dim>> body_force) 
 *   { 
 *     AssertDimension(body_force->n_components, dim); 
 * 
 *     this->body_force = std::move(body_force); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::project( 
 *     const Function<dim> &                       function, 
 *     LinearAlgebra::distributed::Vector<Number> &solution) const 
 *   { 
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number, VectorizedArrayType> 
 *       phi(data, 0, 1); 
 *     MatrixFreeOperators::CellwiseInverseMassMatrix<dim, 
 *                                                    degree, 
 *                                                    dim + 2, 
 *                                                    Number, 
 *                                                    VectorizedArrayType> 
 *       inverse(phi); 
 *     solution.zero_out_ghost_values(); 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           phi.submit_dof_value(evaluate_function(function, 
 *                                                  phi.quadrature_point(q)), 
 *                                q); 
 *         inverse.transform_from_q_points_to_basis(dim + 2, 
 *                                                  phi.begin_dof_values(), 
 *                                                  phi.begin_dof_values()); 
 *         phi.set_dof_values(solution); 
 *       } 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   std::array<double, 3> EulerOperator<dim, degree, n_points_1d>::compute_errors( 
 *     const Function<dim> &                             function, 
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const 
 *   { 
 *     TimerOutput::Scope t(timer, "compute errors"); 
 *     double             errors_squared[3] = {}; 
 *     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number, VectorizedArrayType> 
 *       phi(data, 0, 0); 
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         phi.gather_evaluate(solution, EvaluationFlags::values); 
 *         VectorizedArrayType local_errors_squared[3] = {}; 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           { 
 *             const auto error = 
 *               evaluate_function(function, phi.quadrature_point(q)) - 
 *               phi.get_value(q); 
 *             const auto JxW = phi.JxW(q); 
 * 
 *             local_errors_squared[0] += error[0] * error[0] * JxW; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW; 
 *             local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW; 
 *           } 
 *         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell); 
 *              ++v) 
 *           for (unsigned int d = 0; d < 3; ++d) 
 *             errors_squared[d] += local_errors_squared[d][v]; 
 *       } 
 * 
 *     Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared); 
 * 
 *     std::array<double, 3> errors; 
 *     for (unsigned int d = 0; d < 3; ++d) 
 *       errors[d] = std::sqrt(errors_squared[d]); 
 * 
 *     return errors; 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   double EulerOperator<dim, degree, n_points_1d>::compute_cell_transport_speed( 
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const 
 *   { 
 *     TimerOutput::Scope t(timer, "compute transport speed"); 
 *     Number             max_transport = 0; 
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number, VectorizedArrayType> 
 *       phi(data, 0, 1); 
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         phi.gather_evaluate(solution, EvaluationFlags::values); 
 *         VectorizedArrayType local_max = 0.; 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           { 
 *             const auto solution = phi.get_value(q); 
 *             const auto velocity = euler_velocity<dim>(solution); 
 *             const auto pressure = euler_pressure<dim>(solution); 
 * 
 *             const auto          inverse_jacobian = phi.inverse_jacobian(q); 
 *             const auto          convective_speed = inverse_jacobian * velocity; 
 *             VectorizedArrayType convective_limit = 0.; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               convective_limit = 
 *                 std::max(convective_limit, std::abs(convective_speed[d])); 
 * 
 *             const auto speed_of_sound = 
 *               std::sqrt(gamma * pressure * (1. / solution[0])); 
 * 
 *             Tensor<1, dim, VectorizedArrayType> eigenvector; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               eigenvector[d] = 1.; 
 *             for (unsigned int i = 0; i < 5; ++i) 
 *               { 
 *                 eigenvector = transpose(inverse_jacobian) * 
 *                               (inverse_jacobian * eigenvector); 
 *                 VectorizedArrayType eigenvector_norm = 0.; 
 *                 for (unsigned int d = 0; d < dim; ++d) 
 *                   eigenvector_norm = 
 *                     std::max(eigenvector_norm, std::abs(eigenvector[d])); 
 *                 eigenvector /= eigenvector_norm; 
 *               } 
 *             const auto jac_times_ev   = inverse_jacobian * eigenvector; 
 *             const auto max_eigenvalue = std::sqrt( 
 *               (jac_times_ev * jac_times_ev) / (eigenvector * eigenvector)); 
 *             local_max = 
 *               std::max(local_max, 
 *                        max_eigenvalue * speed_of_sound + convective_limit); 
 *           } 
 * 
 *         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell); 
 *              ++v) 
 *           for (unsigned int d = 0; d < 3; ++d) 
 *             max_transport = std::max(max_transport, local_max[v]); 
 *       } 
 * 
 *     max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD); 
 * 
 *     return max_transport; 
 *   } 
 * 
 *   template <int dim>
 *   class EulerProblem 
 *   { 
 *   public: 
 *     EulerProblem(); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid_and_dofs(); 
 * 
 *     void output_results(const unsigned int result_number); 
 * 
 *     LinearAlgebra::distributed::Vector<Number> solution; 
 * 
 *     ConditionalOStream pcout;
 * 
 * 
 * #ifdef DEAL_II_WITH_P4EST 
 *     parallel::distributed::Triangulation<dim> triangulation; 
 * #else 
 *     Triangulation<dim> triangulation; 
 * #endif 
 * 
 *     FESystem<dim>        fe; 
 *     MappingQGeneric<dim> mapping; 
 *     DoFHandler<dim>      dof_handler; 
 * 
 *     TimerOutput timer; 
 * 
 *     EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator; 
 * 
 *     double time, time_step; 
 * 
 *     class Postprocessor : public DataPostprocessor<dim> 
 *     { 
 *     public: 
 *       Postprocessor(); 
 * 
 *       virtual void evaluate_vector_field( 
 *         const DataPostprocessorInputs::Vector<dim> &inputs, 
 *         std::vector<Vector<double>> &computed_quantities) const override; 
 * 
 *       virtual std::vector<std::string> get_names() const override; 
 * 
 *       virtual std::vector< 
 *         DataComponentInterpretation::DataComponentInterpretation> 
 *       get_data_component_interpretation() const override; 
 * 
 *       virtual UpdateFlags get_needed_update_flags() const override; 
 * 
 *     private: 
 *       const bool do_schlieren_plot; 
 *     }; 
 *   }; 
 * 
 *   template <int dim> 
 *   EulerProblem<dim>::Postprocessor::Postprocessor() 
 *     : do_schlieren_plot(dim == 2) 
 *   {} 
 * 
 *   template <int dim> 
 *   void EulerProblem<dim>::Postprocessor::evaluate_vector_field( 
 *     const DataPostprocessorInputs::Vector<dim> &inputs, 
 *     std::vector<Vector<double>> &               computed_quantities) const 
 *   { 
 *     const unsigned int n_evaluation_points = inputs.solution_values.size(); 
 * 
 *     if (do_schlieren_plot == true) 
 *       Assert(inputs.solution_gradients.size() == n_evaluation_points, 
 *              ExcInternalError()); 
 * 
 *     Assert(computed_quantities.size() == n_evaluation_points, 
 *            ExcInternalError()); 
 *     Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError()); 
 *     Assert(computed_quantities[0].size() == 
 *              dim + 2 + (do_schlieren_plot == true ? 1 : 0), 
 *            ExcInternalError()); 
 * 
 *     for (unsigned int q = 0; q < n_evaluation_points; ++q) 
 *       { 
 *         Tensor<1, dim + 2> solution; 
 *         for (unsigned int d = 0; d < dim + 2; ++d) 
 *           solution[d] = inputs.solution_values[q](d); 
 * 
 *         const double         density  = solution[0]; 
 *         const Tensor<1, dim> velocity = euler_velocity<dim>(solution); 
 *         const double         pressure = euler_pressure<dim>(solution); 
 * 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           computed_quantities[q](d) = velocity[d]; 
 *         computed_quantities[q](dim)     = pressure; 
 *         computed_quantities[q](dim + 1) = std::sqrt(gamma * pressure / density); 
 * 
 *         if (do_schlieren_plot == true) 
 *           computed_quantities[q](dim + 2) = 
 *             inputs.solution_gradients[q][0] * inputs.solution_gradients[q][0]; 
 *       } 
 *   } 
 * 
 *   template <int dim> 
 *   std::vector<std::string> EulerProblem<dim>::Postprocessor::get_names() const 
 *   { 
 *     std::vector<std::string> names; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       names.emplace_back("velocity"); 
 *     names.emplace_back("pressure"); 
 *     names.emplace_back("speed_of_sound"); 
 * 
 *     if (do_schlieren_plot == true) 
 *       names.emplace_back("schlieren_plot"); 
 * 
 *     return names; 
 *   } 
 * 
 *   template <int dim> 
 *   std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *   EulerProblem<dim>::Postprocessor::get_data_component_interpretation() const 
 *   { 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       interpretation; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       interpretation.push_back( 
 *         DataComponentInterpretation::component_is_part_of_vector); 
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
 * 
 *     if (do_schlieren_plot == true) 
 *       interpretation.push_back( 
 *         DataComponentInterpretation::component_is_scalar); 
 * 
 *     return interpretation; 
 *   } 
 * 
 *   template <int dim> 
 *   UpdateFlags EulerProblem<dim>::Postprocessor::get_needed_update_flags() const 
 *   { 
 *     if (do_schlieren_plot == true) 
 *       return update_values | update_gradients; 
 *     else 
 *       return update_values; 
 *   } 
 * 
 *   template <int dim> 
 *   EulerProblem<dim>::EulerProblem() 
 *     : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
 * #ifdef DEAL_II_WITH_P4EST 
 *     , triangulation(MPI_COMM_WORLD) 
 * #endif 
 *     , fe(FE_DGQ<dim>(fe_degree), dim + 2) 
 *     , mapping(fe_degree) 
 *     , dof_handler(triangulation) 
 *     , timer(pcout, TimerOutput::never, TimerOutput::wall_times) 
 *     , euler_operator(timer) 
 *     , time(0) 
 *     , time_step(0) 
 *   {} 
 * 
 *   template <int dim> 
 *   void EulerProblem<dim>::make_grid_and_dofs() 
 *   { 
 *     switch (testcase) 
 *       { 
 *         case 0: 
 *           { 
 *             Point<dim> lower_left; 
 *             for (unsigned int d = 1; d < dim; ++d) 
 *               lower_left[d] = -5; 
 * 
 *             Point<dim> upper_right; 
 *             upper_right[0] = 10; 
 *             for (unsigned int d = 1; d < dim; ++d) 
 *               upper_right[d] = 5; 
 * 
 *             GridGenerator::hyper_rectangle(triangulation, 
 *                                            lower_left, 
 *                                            upper_right); 
 *             triangulation.refine_global(2); 
 * 
 *             euler_operator.set_inflow_boundary( 
 *               0, std::make_unique<ExactSolution<dim>>(0)); 
 * 
 *             break; 
 *           } 
 * 
 *         case 1: 
 *           { 
 *             GridGenerator::channel_with_cylinder( 
 *               triangulation, 0.03, 1, 0, true); 
 * 
 *             euler_operator.set_inflow_boundary( 
 *               0, std::make_unique<ExactSolution<dim>>(0)); 
 *             euler_operator.set_subsonic_outflow_boundary( 
 *               1, std::make_unique<ExactSolution<dim>>(0)); 
 * 
 *             euler_operator.set_wall_boundary(2); 
 *             euler_operator.set_wall_boundary(3); 
 * 
 *             if (dim == 3) 
 *               euler_operator.set_body_force( 
 *                 std::make_unique<Functions::ConstantFunction<dim>>( 
 *                   std::vector<double>({0., 0., -0.2}))); 
 * 
 *             break; 
 *           } 
 * 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 *     triangulation.refine_global(n_global_refinements); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     euler_operator.reinit(mapping, dof_handler); 
 *     euler_operator.initialize_vector(solution); 
 * 
 *     std::locale s = pcout.get_stream().getloc(); 
 *     pcout.get_stream().imbue(std::locale("")); 
 *     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *           << " ( = " << (dim + 2) << " [vars] x " 
 *           << triangulation.n_global_active_cells() << " [cells] x " 
 *           << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )" 
 *           << std::endl; 
 *     pcout.get_stream().imbue(s); 
 *   } 
 * 
 *   template <int dim> 
 *   void EulerProblem<dim>::output_results(const unsigned int result_number) 
 *   { 
 *     const std::array<double, 3> errors = 
 *       euler_operator.compute_errors(ExactSolution<dim>(time), solution); 
 *     const std::string quantity_name = testcase == 0 ? "error" : "norm"; 
 * 
 *     pcout << "Time:" << std::setw(8) << std::setprecision(3) << time 
 *           << ", dt: " << std::setw(8) << std::setprecision(2) << time_step 
 *           << ", " << quantity_name << " rho: " << std::setprecision(4) 
 *           << std::setw(10) << errors[0] << ", rho * u: " << std::setprecision(4) 
 *           << std::setw(10) << errors[1] << ", energy:" << std::setprecision(4) 
 *           << std::setw(10) << errors[2] << std::endl; 
 * 
 *     { 
 *       TimerOutput::Scope t(timer, "output"); 
 * 
 *       Postprocessor postprocessor; 
 *       DataOut<dim>  data_out; 
 * 
 *       DataOutBase::VtkFlags flags; 
 *       flags.write_higher_order_cells = true; 
 *       data_out.set_flags(flags); 
 * 
 *       data_out.attach_dof_handler(dof_handler); 
 *       { 
 *         std::vector<std::string> names; 
 *         names.emplace_back("density"); 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           names.emplace_back("momentum"); 
 *         names.emplace_back("energy"); 
 * 
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *           interpretation; 
 *         interpretation.push_back( 
 *           DataComponentInterpretation::component_is_scalar); 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           interpretation.push_back( 
 *             DataComponentInterpretation::component_is_part_of_vector); 
 *         interpretation.push_back( 
 *           DataComponentInterpretation::component_is_scalar); 
 * 
 *         data_out.add_data_vector(dof_handler, solution, names, interpretation); 
 *       } 
 *       data_out.add_data_vector(solution, postprocessor); 
 * 
 *       LinearAlgebra::distributed::Vector<Number> reference; 
 *       if (testcase == 0 && dim == 2) 
 *         { 
 *           reference.reinit(solution); 
 *           euler_operator.project(ExactSolution<dim>(time), reference); 
 *           reference.sadd(-1., 1, solution); 
 *           std::vector<std::string> names; 
 *           names.emplace_back("error_density"); 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             names.emplace_back("error_momentum"); 
 *           names.emplace_back("error_energy"); 
 * 
 *           std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *             interpretation; 
 *           interpretation.push_back( 
 *             DataComponentInterpretation::component_is_scalar); 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             interpretation.push_back( 
 *               DataComponentInterpretation::component_is_part_of_vector); 
 *           interpretation.push_back( 
 *             DataComponentInterpretation::component_is_scalar); 
 * 
 *           data_out.add_data_vector(dof_handler, 
 *                                    reference, 
 *                                    names, 
 *                                    interpretation); 
 *         } 
 * 
 *       Vector<double> mpi_owner(triangulation.n_active_cells()); 
 *       mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD); 
 *       data_out.add_data_vector(mpi_owner, "owner"); 
 * 
 *       data_out.build_patches(mapping, 
 *                              fe.degree, 
 *                              DataOut<dim>::curved_inner_cells); 
 * 
 *       const std::string filename = 
 *         "solution_" + Utilities::int_to_string(result_number, 3) + ".vtu"; 
 *       data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD); 
 *     } 
 *   } 
 * 
 *   template <int dim> 
 *   void EulerProblem<dim>::run() 
 *   { 
 *     { 
 *       const unsigned int n_vect_number = VectorizedArrayType::size(); 
 *       const unsigned int n_vect_bits   = 8 * sizeof(Number) * n_vect_number; 
 * 
 *       pcout << "Running with " 
 *             << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) 
 *             << " MPI processes" << std::endl; 
 *       pcout << "Vectorization over " << n_vect_number << " " 
 *             << (std::is_same<Number, double>::value ? "doubles" : "floats") 
 *             << " = " << n_vect_bits << " bits (" 
 *             << Utilities::System::get_current_vectorization_level() << ")" 
 *             << std::endl; 
 *     } 
 * 
 *     make_grid_and_dofs(); 
 * 
 *     const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme); 
 * 
 *     LinearAlgebra::distributed::Vector<Number> rk_register_1; 
 *     LinearAlgebra::distributed::Vector<Number> rk_register_2; 
 *     rk_register_1.reinit(solution); 
 *     rk_register_2.reinit(solution); 
 * 
 *     euler_operator.project(ExactSolution<dim>(time), solution); 
 * 
 *     double min_vertex_distance = std::numeric_limits<double>::max(); 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         min_vertex_distance = 
 *           std::min(min_vertex_distance, cell->minimum_vertex_distance()); 
 *     min_vertex_distance = 
 *       Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD); 
 * 
 *     time_step = courant_number * integrator.n_stages() / 
 *                 euler_operator.compute_cell_transport_speed(solution); 
 *     pcout << "Time step size: " << time_step 
 *           << ", minimal h: " << min_vertex_distance 
 *           << ", initial transport scaling: " 
 *           << 1. / euler_operator.compute_cell_transport_speed(solution) 
 *           << std::endl 
 *           << std::endl; 
 * 
 *     output_results(0); 
 * 
 *     unsigned int timestep_number = 0; 
 * 
 *     while (time < final_time - 1e-12 && timestep_number < max_time_steps) 
 *       { 
 *         ++timestep_number; 
 *         if (timestep_number % 5 == 0) 
 *           time_step = 
 *             courant_number * integrator.n_stages() / 
 *             Utilities::truncate_to_n_digits( 
 *               euler_operator.compute_cell_transport_speed(solution), 3); 
 * 
 *         { 
 *           TimerOutput::Scope t(timer, "rk time stepping total"); 
 *           integrator.perform_time_step(euler_operator, 
 *                                        time, 
 *                                        time_step, 
 *                                        solution, 
 *                                        rk_register_1, 
 *                                        rk_register_2); 
 *         } 
 * 
 *         time += time_step; 
 * 
 *         if (static_cast<int>(time / output_tick) != 
 *               static_cast<int>((time - time_step) / output_tick) || 
 *             time >= final_time - 1e-12) 
 *           output_results( 
 *             static_cast<unsigned int>(std::round(time / output_tick))); 
 *       } 
 * 
 *     timer.print_wall_time_statistics(MPI_COMM_WORLD); 
 *     pcout << std::endl; 
 *   } 
 * 
 * } // namespace Euler_DG 
 * 
 * int main(int argc, char **argv) 
 * { 
 *   using namespace Euler_DG; 
 *   using namespace dealii; 
 * 
 *   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *   try 
 *     { 
 *       deallog.depth_console(0); 
 * 
 *       EulerProblem<dimension> euler_problem; 
 *       euler_problem.run(); 
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
 * 
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
 * 
 * @endcode
examples/step-76/doc/results.dox



<a name="Results"></a><h1>Results</h1>


在一台有40个进程的机器上以默认设置运行该程序，会产生以下输出。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 27.648.000 ( = 5 [vars] x 25.600 [cells] x 216 [dofs/cell/var] )
Time step size: 0.000295952, minimal h: 0.0075, initial transport scaling: 0.00441179
Time:       0, dt:   0.0003, norm rho:  5.385e-16, rho * u:  1.916e-16, energy: 1.547e-15
+--------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed         |     17.52s    10 |     17.52s |     17.52s    11 |
|                                      |                  |                               |
| Section                  | no. calls |   min time  rank |   avg time |   max time  rank |
+--------------------------------------+------------------+------------+------------------+
| compute errors           |         1 |  0.009594s    16 |  0.009705s |  0.009819s     8 |
| compute transport speed  |        22 |    0.1366s     0 |    0.1367s |    0.1368s    18 |
| output                   |         1 |     1.233s     0 |     1.233s |     1.233s    32 |
| rk time stepping total   |       100 |     8.746s    35 |     8.746s |     8.746s     0 |
| rk_stage - integrals L_h |       500 |     8.742s    36 |     8.742s |     8.743s     2 |
+--------------------------------------+------------------+------------+------------------+
@endcode



和以下视觉输出。

 <table align="center" class="doxtable" style="width:85%">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_010.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_025.png" alt="" width="100%">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_050.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_100.png" alt="" width="100%">
    </td>
  </tr>
</table> 

作为参考，使用FCL的步骤-67的结果是。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 27.648.000 ( = 5 [vars] x 25.600 [cells] x 216 [dofs/cell/var] )
Time step size: 0.000295952, minimal h: 0.0075, initial transport scaling: 0.00441179
Time:       0, dt:   0.0003, norm rho:  5.385e-16, rho * u:  1.916e-16, energy: 1.547e-15
+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     13.33s     0 |     13.34s |     13.35s    34 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |         1 |  0.007977s    10 |  0.008053s |  0.008161s    30 |
| compute transport speed       |        22 |    0.1228s    34 |    0.2227s |    0.3845s     0 |
| output                        |         1 |     1.255s     3 |     1.257s |     1.259s    27 |
| rk time stepping total        |       100 |     11.15s     0 |     11.32s |     11.42s    34 |
| rk_stage - integrals L_h      |       500 |     8.719s    10 |     8.932s |     9.196s     0 |
| rk_stage - inv mass + vec upd |       500 |     1.944s     0 |     2.377s |      2.55s    10 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



通过本教程中的修改，我们能够使Runge-Kutta阶段的速度提高27%。

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这些算法很容易扩展到更高的维度：高维的<a href="https://github.com/hyperdeal/hyperdeal/blob/a9e67b4e625ff1dde2fed93ad91cdfacfaa3acdf/include/hyper.deal/operators/advection/advection_operation.h#L219-L569">advection operator based on cell-centric loops</a>是hyper.deal库的一部分。以单元为中心的循环扩展到局部细化的网格则涉及更多。

<a name="ExtensiontothecompressibleNavierStokesequations"></a><h4>Extension to the compressible Navier-Stokes equations</h4>


本教程中提出的求解器也可以通过增加粘性项来扩展到可压缩的Navier-Stokes方程，这也是步骤67中的建议。为了尽量保持这里获得的性能，尽管有额外的椭圆项的成本，例如通过内部惩罚方法，该教程建议将基础从FE_DGQ切换到FE_DGQHermite，就像步骤59的教程程序一样。这种转换背后的原因是，在FE_DGQ的情况下，需要相邻单元的所有值（即 $k+1$ 层），而在FE_DGQHermite的情况下，只需要2层，这使得后者明显更适合于高度数。额外的层一方面要在通量计算过程中从主内存加载，另一方面要进行通信。利用本教程介绍的共享内存能力，第二点可以在单个计算节点上消除，或者在混合环境下减少其影响。

<a name="BlockGaussSeidellikepreconditioners"></a><h4>Block Gauss-Seidel-like preconditioners</h4>


以单元为中心的循环可用于创建块状高斯-赛德尔预处理，在一个过程中是乘法的，在整个过程中是加法的。这些类型的预处理器在通量计算过程中使用，与雅可比型预处理器相反，已经从相邻的单元中更新了数值。下面的伪代码直观地说明了这在原则上是如何实现的。

@code
// vector monitor if cells have been updated or not
Vector<Number> visit_flags(data.n_cell_batches () + data.n_ghost_cell_batches ());


// element centric loop with a modified kernel
data.template loop_cell_centric<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto cell_range) {


    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        // cell integral as usual (not shown)


        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          {
            const auto boundary_id = data.get_faces_by_cells_boundary_id(cell, face)[0];


            if (boundary_id == numbers::internal_face_boundary_id)
              {
                phi_p.reinit(cell, face);


                const auto flags = phi_p.read_cell_data(visit_flags);
                const auto all_neighbors_have_been_updated =
                  std::min(flags.begin(),
                           flags().begin() + data.n_active_entries_per_cell_batch(cell) == 1;


                if(all_neighbors_have_been_updated)
                  phi_p.gather_evaluate(dst, EvaluationFlags::values);
                else
                  phi_p.gather_evaluate(src, EvaluationFlags::values);


                // continue as usual (not shown)
              }
            else
              {
                // boundary integral as usual (not shown)
              }
          }


        // continue as above and apply your favorite algorithm to invert
        // the cell-local operator (not shown)


        // make cells as updated
        phi.set_cell_data(visit_flags, VectorizedArrayType(1.0));
      }
  },
  dst,
  src,
  true,
  MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::values);
@endcode



为此，我们可以利用MatrixFree的单元数据向量能力和VectorizedArray的基于范围的迭代能力。

请注意，在给定的例子中，我们处理 <code>VectorizedArrayType::size()</code> 个块，因为每个通道对应一个块。如果一个矢量寄存器处理的所有块都被更新了，我们就认为块被更新了。在笛卡尔网格的情况下，这是一个合理的方法，然而，对于一般的非结构化网格，这种保守的方法可能会导致预处理程序的效率下降。通过明确减少 <code>VectorizedArrayType</code> 使用的通道数量来减少并行处理的单元可能会提高预处理器的质量，但代价是每次迭代可能会更昂贵。这种两难境地把我们引向另一种 "扩展的可能性"：元素内的矢量化。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-76.cc"
*/
