

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2009 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE.md at 
 * the top level directory of deal.II. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Authors: Katharina Kormann, Martin Kronbichler, Uppsala University, 
 * 2009-2012, updated to MPI version with parallel vectors in 2016 
 */ 



// 首先包括deal.II库中的必要文件。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/timer.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/la_parallel_vector.h> 
#include <deal.II/lac/precondition.h> 

#include <deal.II/fe/fe_q.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer_matrix_free.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_matrix.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 

// 这包括有效实现无矩阵方法的数据结构，或者用MatrixFree类的更通用的有限元算子。

#include <deal.II/matrix_free/matrix_free.h> 
#include <deal.II/matrix_free/operators.h> 
#include <deal.II/matrix_free/fe_evaluation.h> 

#include <iostream> 
#include <fstream> 

namespace Step37 
{ 
  using namespace dealii; 

// 为了提高效率，在无矩阵实现中进行的操作需要在编译时了解循环长度，这些长度是由有限元的度数给出的。因此，我们收集了两个模板参数的值，可以在代码中的一个地方改变。当然，我们可以把有限元的度数作为一个运行时的参数，通过编译所有可能的度数（比如，1到6之间）的计算核，并在运行时选择合适的核。在这里，我们只是选择二阶 $Q_2$ 元素，并选择维度3作为标准。

  const unsigned int degree_finite_element = 2; 
  const unsigned int dimension             = 3; 
// @sect3{Equation data}  

// 我们为泊松问题定义了一个可变系数函数。它与 step-5 中的函数类似，但我们使用 $a(\mathbf x)=\frac{1}{0.05 + 2\|\bf x\|^2}$ 的形式，而不是不连续的形式。这只是为了证明这种实现的可能性，而不是在物理上有什么意义。我们定义系数的方式与早期教程程序中的函数相同。有一个新的函数，即有模板参数 @p value 的 @p number. 方法。
template <int dim> 
  class Coefficient : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    template <typename number> 
    number value(const Point<dim, number> &p, 
                 const unsigned int        component = 0) const; 
  }; 

// 这就是上面提到的新函数。评估抽象类型的系数  @p number.  它可能只是一个普通的双数，但也可能是一个有点复杂的类型，我们称之为VectorizedArray。这种数据类型本质上是一个短的双数数组，正如在介绍中所讨论的那样，它可以容纳几个单元格的数据。例如，我们在这里评估的系数不是像通常那样在一个简单的点上，而是交给一个Point<dim,VectorizedArray<double>>点，在AVX的情况下，它实际上是四个点的集合。不要把VectorizedArray中的条目与点的不同坐标混淆。事实上，数据的布局是这样的： <code>p[0]</code> 返回一个VectorizedArray，它又包含了第一个点和第二个点的x坐标。你可以使用例如  <code>p[0][j]</code>  单独访问坐标，j=0,1,2,3，但建议尽可能在一个VectorizedArray上定义操作，以便利用矢量操作。

// 在函数的实现中，我们假设数字类型重载了基本的算术运算，所以我们只需照常写代码。然后，基类函数 @p value 是由带有双倍类型的模板函数计算出来的，以避免重复代码。

  template <int dim> 
  template <typename number> 
  number Coefficient<dim>::value(const Point<dim, number> &p, 
                                 const unsigned int /*component*/) const 
  { 
    return 1. / (0.05 + 2. * p.square()); 
  } 

  template <int dim> 
  double Coefficient<dim>::value(const Point<dim> & p, 
                                 const unsigned int component) const 
  { 
    return value<double>(p, component); 
  } 
// @sect3{Matrix-free implementation}  

// 下面这个名为 <code>LaplaceOperator</code> 的类，实现了微分运算符。就所有的实用目的而言，它是一个矩阵，也就是说，你可以向它询问它的大小（成员函数  <code>m(), n()</code>  ），你可以将它应用于一个矢量（ <code>vmult()</code>  函数）。当然，与实数矩阵的区别在于，这个类实际上并不存储矩阵的<i>elements</i>，而只知道如何计算运算器应用于向量时的动作。

// 描述矩阵大小的基础结构，来自MatrixFree对象的初始化，以及通过vmult()和Tvmult()方法实现矩阵-向量乘积的各种接口，是由本类派生的 MatrixFreeOperator::Base 类提供的。这里定义的LaplaceOperator类只需要提供几个接口，即通过vmult()函数中使用的apply_add()方法来实现运算符的实际操作，以及计算底层矩阵对角线项的方法。我们需要对角线来定义多梯度平滑器。由于我们考虑的是一个具有可变系数的问题，我们进一步实现了一个可以填充系数值的方法。

// 注意文件 <code>include/deal.II/matrix_free/operators.h</code> 已经包含了通过类 MatrixFreeOperators::LaplaceOperator. 对拉普拉斯的实现。 出于教育目的，本教程程序中重新实现了该运算符，解释了其中的成分和概念。

// 这个程序利用了集成在deal.II中的有限元算子应用的数据缓存。这个数据缓存类被称为MatrixFree。它包含局部和全局自由度之间的映射信息（Jacobian）和索引关系。它还包含约束条件，如来自悬挂节点或迪里切特边界条件的约束。此外，它可以在所有单元上以%并行方式发出一个循环，确保只有不共享任何自由度的单元被处理（这使得循环在写入目标向量时是线程安全的）。与 @ref threads 模块中描述的WorkStream类相比，这是一个更先进的策略。当然，为了不破坏线程安全，我们在写进类全局结构时必须小心。

// 实现拉普拉斯算子的类有三个模板参数，一个是维度（正如许多deal.II类所携带的），一个是有限元的度数（我们需要通过FEEvaluation类来实现高效计算），还有一个是底层标量类型。我们希望对最终矩阵使用 <code>double</code> 数字（即双精度，64位浮点），但对多网格级矩阵使用浮点数（单精度，32位浮点数字）（因为那只是一个预处理程序，而浮点数的处理速度是两倍）。FEEvaluation类也需要一个模板参数，用于确定一维正交点的数量。在下面的代码中，我们把它硬编码为  <code>fe_degree+1</code>  。如果我们想独立于多项式程度来改变它，我们需要添加一个模板参数，就像在  MatrixFreeOperators::LaplaceOperator  类中做的那样。

// 顺便说一下，如果我们在同一个网格和自由度上实现了几个不同的操作（比如质量矩阵和拉普拉斯矩阵），我们将为每个操作者定义两个像现在这样的类（来自于 MatrixFreeOperators::Base 类），并让它们都引用一般问题类中的同一个MatrixFree数据缓存。通过 MatrixFreeOperators::Base 的接口要求我们只提供一组最小的函数。这个概念允许编写具有许多无矩阵操作的复杂应用代码。

//  @note  储存类型 <code>VectorizedArray<number></code> 的值需要注意。在这里，我们使用deal.II表类，它准备以正确的对齐方式保存数据。然而，存储例如一个 <code>std::vector<VectorizedArray<number> ></code> 是不可能用矢量化的。数据与内存地址的边界需要一定的对齐（基本上，在AVX的情况下，一个32字节的VectorizedArray需要从一个能被32整除的内存地址开始）。表类（以及它所基于的AlignedVector类）确保这种对齐方式得到尊重，而 std::vector 一般不这样做，这可能会导致一些系统在奇怪的地方出现分段故障，或者其他系统的性能不理想。

  template <int dim, int fe_degree, typename number> 
  class LaplaceOperator 
    : public MatrixFreeOperators:: 
        Base<dim, LinearAlgebra::distributed::Vector<number>> 
  { 
  public: 
    using value_type = number; 

    LaplaceOperator(); 

    void clear() override; 

    void evaluate_coefficient(const Coefficient<dim> &coefficient_function); 

    virtual void compute_diagonal() override; 

  private: 
    virtual void apply_add( 
      LinearAlgebra::distributed::Vector<number> &      dst, 
      const LinearAlgebra::distributed::Vector<number> &src) const override; 

    void 
    local_apply(const MatrixFree<dim, number> &                   data, 
                LinearAlgebra::distributed::Vector<number> &      dst, 
                const LinearAlgebra::distributed::Vector<number> &src, 
                const std::pair<unsigned int, unsigned int> &cell_range) const; 

    void local_compute_diagonal( 
      const MatrixFree<dim, number> &              data, 
      LinearAlgebra::distributed::Vector<number> & dst, 
      const unsigned int &                         dummy, 
      const std::pair<unsigned int, unsigned int> &cell_range) const; 

    Table<2, VectorizedArray<number>> coefficient; 
  }; 

// 这是 @p LaplaceOperator 类的构造函数。它所做的就是调用基类 MatrixFreeOperators::Base, 的默认构造函数，而基类又是基于Subscriptor类的，它断言这个类在超出范围后不会被访问，比如在一个预处理程序中。

  template <int dim, int fe_degree, typename number> 
  LaplaceOperator<dim, fe_degree, number>::LaplaceOperator() 
    : MatrixFreeOperators::Base<dim, 
                                LinearAlgebra::distributed::Vector<number>>() 
  {} 

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::clear() 
  { 
    coefficient.reinit(0, 0); 
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>:: 
      clear(); 
  } 

//  @sect4{Computation of coefficient}  

// 为了初始化系数，我们直接赋予它上面定义的系数类，然后选择带有矢量数的方法 <code>coefficient_function.value</code> （编译器可以从点数据类型中推导出来）。下面将解释FEEvaluation类（及其模板参数）的使用。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::evaluate_coefficient( 
    const Coefficient<dim> &coefficient_function) 
  { 
    const unsigned int n_cells = this->data->n_cell_batches(); 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*this->data); 

    coefficient.reinit(n_cells, phi.n_q_points); 
    for (unsigned int cell = 0; cell < n_cells; ++cell) 
      { 
        phi.reinit(cell); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          coefficient(cell, q) = 
            coefficient_function.value(phi.quadrature_point(q)); 
      } 
  } 

//  @sect4{Local evaluation of Laplace operator}  

// 这里是这个类的主要功能，矩阵-向量乘积的评估（或者，一般来说，有限元算子评估）。这是在一个函数中完成的，该函数需要四个参数，MatrixFree对象，目标和源向量，以及要处理的单元格范围。MatrixFree类中的方法 <code>cell_loop</code> 将在内部用一些单元格范围来调用这个函数，这些单元格范围是通过检查哪些单元格可以同时工作来获得的，这样写操作就不会引起任何竞赛条件。请注意，循环中使用的单元格范围并不是直接指当前网格中的（活动）单元格数量，而是一个单元格批次的集合。 换句话说，"单元 "可能是一个错误的开始，因为FEEvaluation将几个单元的数据分组在一起。这意味着在正交点的循环中，我们实际上是将几个单元的正交点作为一个块来看待。这样做是为了实现更高的矢量化程度。 这种 "单元 "或 "单元批 "的数量存储在MatrixFree中，可以通过 MatrixFree::n_cell_batches(). 查询。与deal.II单元迭代器相比，在这个类中，所有的单元都被布置在一个普通的数组中，不直接知道水平或相邻关系，这使得通过无符号整数索引单元成为可能。

// 拉普拉斯运算符的实现非常简单。首先，我们需要创建一个对象FEEvaluation，它包含计算核，并有数据字段来存储临时结果（例如，在几个单元格集合的所有正交点上评估的梯度）。请注意，临时结果不会使用大量的内存，而且由于我们用元素顺序指定模板参数，数据被存储在堆栈中（没有昂贵的内存分配）。通常，只需要设置两个模板参数，维度作为第一个参数，有限元的度数作为第二个参数（这等于每个维度的自由度数减去FE_Q元素的一个）。然而，在这里，我们也希望能够使用浮点数来计算多网格预处理，这是最后一个（第五个）模板参数。因此，我们不能依赖默认的模板参数，因此必须填写第三和第四个字段。第三个参数指定每个方向的正交点的数量，其默认值等于元素的度数加1。第四个参数设置分量的数量（在PDEs系统中也可以评估矢量值的函数，但默认是标量元素），最后一个参数设置数字类型。

// 接下来，我们在给定的单元格范围内循环，然后继续进行实际的实现。  <ol>  
// <li>  告诉FEEvaluation对象我们要处理的（宏）单元。   <li>  读入源向量的值（  @p read_dof_values),  包括约束的解析。这将存储 $u_\mathrm{cell}$ ，如介绍中所述。   <li>  计算单元格梯度（有限元函数的评价）。由于FEEvaluation可以结合值计算和梯度计算，它使用一个统一的接口来处理0到2阶之间的各种导数。我们只想要梯度，不想要值，也不想要二阶导数，所以我们在梯度槽（第二槽）中将函数参数设置为真，而在值槽（第一槽）中设置为假。还有一个用于Hessian的第三槽，默认为假，所以不需要给它。请注意，FEEvaluation类在内部以一种有效的方式评估形状函数，一次只处理一个维度（如介绍中提到的使用形状函数和正交点的张量积形式）。与FEValues中使用的在所有局部自由度和正交点上循环的天真方法相比，在 $d$ 维度上，这给出了等于 $\mathcal O(d^2 (p+1)^{d+1})$ 的多项式度数 $p$ 的复杂度，并花费了 $\mathcal O(d (p+1)^{2d})$  。   <li>  接下来是雅各布变换的应用，乘以变量系数和正交权重。FEEvaluation有一个访问函数 @p get_gradient ，可以应用Jacobian并返回实空间中的梯度。然后，我们只需要乘以（标量）系数，并让函数 @p submit_gradient 应用第二个雅各布式（用于测试函数）和正交权重及雅各布式行列式（JxW）。注意，提交的梯度存储在与 @p get_gradient. 中读取梯度的地方相同的数据字段中。因此，你需要确保在调用 @p submit_gradient 后不要再从同一正交点读取该特定正交点。一般来说，当 @p get_gradient 被多次使用时，复制其结果是个好主意。   <li>  接下来是对所有测试函数的正交点进行求和，对应于实际积分步骤。对于拉普拉斯算子，我们只是乘以梯度，所以我们用各自的参数集调用积分函数。如果你有一个方程，同时用测试函数的值和梯度进行测试，那么两个模板参数都需要设置为真。先调用积分函数的值，再单独调用梯度，会导致错误的结果，因为第二次调用会在内部覆盖第一次调用的结果。请注意，积分步骤的二次导数没有函数参数。   <li>  最终，介绍中提到的向量 $v_\mathrm{cell}$ 中的局部贡献需要被添加到结果向量中（并应用约束）。这是通过调用 @p distribute_local_to_global, 来完成的，该函数与AffineConstraints中的相应函数名称相同（只是我们现在将局部向量存储在FEEvaluation对象中，正如局部和全局自由度之间的指数一样）。   </ol>  
  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::local_apply( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     cell_range) const 
  { 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data); 

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        AssertDimension(coefficient.size(0), data.n_cell_batches()); 
        AssertDimension(coefficient.size(1), phi.n_q_points); 

        phi.reinit(cell); 
        phi.read_dof_values(src); 
        phi.evaluate(EvaluationFlags::gradients); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q), q); 
        phi.integrate(EvaluationFlags::gradients); 
        phi.distribute_local_to_global(dst); 
      } 
  } 

// 这个函数实现了对 Base::apply_add() 接口的所有单元的循环。这是用MatrixFree类的 @p cell_loop 来实现的，它接受这个类的operator()，参数为MatrixFree, OutVector, InVector, cell_range。当使用MPI并行化（但没有线程）时，如本教程程序中所做的，单元格循环对应于以下三行代码。

// 
// @code
//  src.update_ghost_values();
//  local_apply(*this->data, dst, src, std::make_pair(0U,
//                                                    data.n_cell_batches()));
//  dst.compress(VectorOperation::add);
//  @endcode

// 这里，两个调用update_ghost_values()和compress()为MPI执行处理器边界上的数据交换，一次用于源向量，我们需要从远程处理器拥有的条目中读取，一次用于目的向量，我们已经积累了部分残余，需要添加到所有者处理器的相应条目中。然而， MatrixFree::cell_loop 不仅抽象出这两个调用，而且还进行了一些额外的优化。一方面，它将把update_ghost_values()和compress()的调用拆开，以允许通信和计算的重叠。然后用三个代表从0到 MatrixFree::n_cell_batches(). 的单元格范围的分区来调用local_apply函数。另一方面，cell_loop也支持线程并行，在这种情况下，单元格范围被分割成更小的块，并以一种先进的方式安排，避免了几个线程对同一个向量条目的访问。这一特性在  step-48  中有解释。

// 注意，在单元格循环之后，受约束的自由度需要再次被触及，以实现合理的vmult()操作。由于装配循环会自动解决约束问题（就像 AffineConstraints::distribute_local_to_global() 的调用一样），它不会计算对受约束自由度的任何贡献，而是将各自的条目留为零。这将表示一个矩阵的受限自由度的行和列都是空的。然而，像CG这样的迭代求解器只对非星形矩阵有效。最简单的方法是将矩阵中对应于受限自由度的子块设置为同一矩阵，在这种情况下，矩阵的应用只是将右侧向量的元素复制到左侧。幸运的是，vmult()的实现 MatrixFreeOperators::Base 在apply_add()函数之外自动为我们做了这个，所以我们不需要在这里采取进一步的行动。

// 当使用MatrixFree和FEEvaluation的组合与MPI并行时，有一个方面需要注意&mdash; 用于访问向量的索引。出于性能的考虑，MatrixFree和FEEvaluation被设计为在MPI本地索引空间中访问向量，当与多个处理器一起工作时也是如此。在本地索引空间工作意味着除了不可避免的间接寻址外，在向量访问发生的地方不需要进行索引转换。然而，本地索引空间是模糊的：虽然标准的惯例是用0和本地大小之间的索引访问向量的本地拥有的范围，但对于重影项的编号并不那么明确，而且有些随意。对于矩阵-向量乘积，只有出现在本地拥有的单元格上的指数（加上那些通过悬挂节点约束引用的指数）是必要的。然而，在deal.II中，我们经常将重影元素上的所有自由度设置为重影向量条目，称为 @ref GlossLocallyRelevantDof "术语表中描述的本地相关DoF"。在这种情况下，尽管指的是同一个全局索引，但在两个可能的重影集中，重影向量条目的MPI本地索引一般会有所不同。为了避免问题，FEEvaluation通过一个名为 LinearAlgebra::distributed::Vector::partitioners_are_compatible. 的检查来检查用于矩阵-向量乘积的向量分区是否确实与MatrixFree中的索引分区相匹配。 为了方便， MatrixFreeOperators::Base 类包括一个机制来使鬼魂集适合正确的布局。这发生在向量的重影区域，所以请记住，在调用vmult()方法后，目标和源向量的重影区域都可能被修改。这是合法的，因为分布式deal.II向量的ghost区域是一个可变的部分，并按需填充。在矩阵-向量乘积中使用的向量在进入vmult()函数时不能被重影，所以没有信息丢失。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::apply_add( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src); 
  } 

// 下面的函数实现了算子对角线的计算。计算无矩阵算子评估的矩阵项，结果比评估算子更复杂。从根本上说，我们可以通过在<i>all</i>单位向量上应用算子来获得算子的矩阵表示。当然，这将是非常低效的，因为我们需要进行<i>n</i>运算符的评估来检索整个矩阵。此外，这种方法会完全忽视矩阵的稀疏性。然而，对于单个单元来说，这是一种方法，而且实际上效率并不低，因为单元内的所有自由度之间通常都存在着耦合。

// 我们首先将对角线向量初始化为正确的平行布局。这个向量被封装在基类 MatrixFreeOperators::Base. 中DiagonalMatrix类型的一个名为inverse_diagonal_entries的成员中，这个成员是一个共享指针，我们首先需要初始化它，然后获得代表矩阵中对角线条目的向量。至于实际的对角线计算，我们再次使用MatrixFree的cell_loop基础设施来调用一个名为local_compute_diagonal()的本地工作程序。由于我们只写进一个向量，而没有任何源向量，我们用一个<tt>unsigned int</tt>类型的假参数来代替源向量，以便与cell_loop接口确认。在循环之后，我们需要将受Dirichlet边界条件约束的向量条目设置为1（要么是MatrixFree内部AffineConstraints对象描述的边界上的条目，要么是自适应多网格中不同网格层次之间的索引）。这是通过函数 MatrixFreeOperators::Base::set_constrained_entries_to_one() 完成的，并与Base算子提供的矩阵-向量乘积中的设置相匹配。最后，我们需要反转对角线条目，这是基于Jacobi迭代的Chebyshev平滑器所要求的形式。在循环中，我们断言所有的条目都是非零的，因为它们应该从积分中获得正的贡献，或者被约束并被 @p set_constrained_entries_to_one() 以下的cell_loop处理。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::compute_diagonal() 
  { 
    this->inverse_diagonal_entries.reset( 
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>()); 
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal = 
      this->inverse_diagonal_entries->get_vector(); 
    this->data->initialize_dof_vector(inverse_diagonal); 
    unsigned int dummy = 0; 
    this->data->cell_loop(&LaplaceOperator::local_compute_diagonal, 
                          this, 
                          inverse_diagonal, 
                          dummy); 

    this->set_constrained_entries_to_one(inverse_diagonal); 

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i) 
      { 
        Assert(inverse_diagonal.local_element(i) > 0., 
               ExcMessage("No diagonal entry in a positive definite operator " 
                          "should be zero")); 
        inverse_diagonal.local_element(i) = 
          1. / inverse_diagonal.local_element(i); 
      } 
  } 

// 在本地计算循环中，我们通过循环本地矩阵中的所有列来计算对角线，并将条目1放在<i>i</i>槽中，将条目0放在所有其他槽中，也就是说，我们一次在一个单位向量上应用单元格微分运算。调用 FEEvaluation::evaluate, 的内部部分是对正交点的循环， FEEvalution::integrate, 则与local_apply函数完全相同。之后，我们挑出本地结果的第<i>i</i>个条目，并将其放入一个临时存储器（因为我们在下一次循环迭代时覆盖了 FEEvaluation::get_dof_value() 后面数组中的所有条目）。最后，临时存储被写到目标向量中。注意我们是如何使用 FEEvaluation::get_dof_value() 和 FEEvaluation::submit_dof_value() 来读取和写入FEEvaluation用于积分的数据字段，并在另一方面写入全局向量的。

// 鉴于我们只对矩阵的对角线感兴趣，我们简单地扔掉了沿途计算过的本地矩阵的所有其他条目。虽然计算完整的单元格矩阵，然后扔掉除对角线以外的所有东西看起来很浪费，但是整合的效率很高，所以计算并没有花费太多时间。请注意，对于多项式度数来说，每个元素的算子评估的复杂度是 $\mathcal O((p+1)^{d+1})$ ，所以计算整个矩阵要花费我们 $\mathcal O((p+1)^{2d+1})$ 次操作，与用FEValues计算对角线的复杂度 $\mathcal O((p+1)^{2d})$ 相差不大。由于FEEvaluation也由于矢量化和其他优化而大大加快了速度，所以用这个函数计算对角线实际上是最快的（简单的）变量。(有可能用 $\mathcal O((p+1)^{d+1})$ 操作中的和分解技术来计算对角线，这涉及到特别适应的内核&mdash;但是由于这种内核只在特定的环境下有用，而对角线计算通常不在关键路径上，所以它们没有在deal.II中实现。)

// 注意在向量上调用distribution_local_to_global来将对角线条目累积到全局矩阵的代码有一些限制。对于带有悬空节点约束的操作者来说，在distribution_local_to_global的调用中，将一个受约束的DoF的积分贡献分配给其他几个条目，这里使用的向量接口并不完全计算对角线条目，而是将一些位于本地矩阵对角线上的贡献，最终在全局矩阵的非对角线位置堆积到对角线上。如<a href="http:dx.doi.org/10.4208/cicp.101214.021015a">Kormann (2016), section 5.3</a>中所解释的，该结果在离散化精度上是正确的，但在数学上并不平等。在这个教程程序中，不会发生任何危害，因为对角线只用于没有悬空节点约束出现的多网格水平矩阵中。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal( 
    const MatrixFree<dim, number> &             data, 
    LinearAlgebra::distributed::Vector<number> &dst, 
    const unsigned int &, 
    const std::pair<unsigned int, unsigned int> &cell_range) const 
  { 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data); 

    AlignedVector<VectorizedArray<number>> diagonal(phi.dofs_per_cell); 

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        AssertDimension(coefficient.size(0), data.n_cell_batches()); 
        AssertDimension(coefficient.size(1), phi.n_q_points); 

        phi.reinit(cell); 
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i) 
          { 
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j) 
              phi.submit_dof_value(VectorizedArray<number>(), j); 
            phi.submit_dof_value(make_vectorized_array<number>(1.), i); 

            phi.evaluate(EvaluationFlags::gradients); 
            for (unsigned int q = 0; q < phi.n_q_points; ++q) 
              phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q), 
                                  q); 
            phi.integrate(EvaluationFlags::gradients); 
            diagonal[i] = phi.get_dof_value(i); 
          } 
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i) 
          phi.submit_dof_value(diagonal[i], i); 
        phi.distribute_local_to_global(dst); 
      } 
  } 

//  @sect3{LaplaceProblem class}  

// 这个类是基于  step-16  中的一个。然而，我们用我们的无矩阵实现取代了SparseMatrix<double>类，这意味着我们也可以跳过稀疏性模式。请注意，我们定义LaplaceOperator类时，将有限元的度数作为模板参数（该值在文件的顶部定义），我们使用浮点数来表示多网格级矩阵。

// 该类还有一个成员变量，用来记录在我们真正去解决这个问题之前设置整个数据链的所有详细时间。此外，还有一个输出流（默认情况下是禁用的），可以用来输出各个设置操作的细节，而不是默认情况下只打印出的摘要。

// 由于这个程序被设计成与MPI一起使用，我们也提供了通常的 @p pcout 输出流，只打印MPI等级为0的处理器的信息。这个程序使用的网格可以是基于p4est的分布式三角图（在deal.II被配置为使用p4est的情况下），否则它就是一个只在没有MPI的情况下运行的串行网格。

  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_rhs(); 
    void solve(); 
    void output_results(const unsigned int cycle) const; 

#ifdef DEAL_II_WITH_P4EST 
    parallel::distributed::Triangulation<dim> triangulation; 
#else 
    Triangulation<dim> triangulation; 
#endif 

    FE_Q<dim>       fe; 
    DoFHandler<dim> dof_handler; 

    MappingQ1<dim> mapping; 

    AffineConstraints<double> constraints;
    using SystemMatrixType = 
      LaplaceOperator<dim, degree_finite_element, double>; 
    SystemMatrixType system_matrix; 

    MGConstrainedDoFs mg_constrained_dofs; 
    using LevelMatrixType = LaplaceOperator<dim, degree_finite_element, float>; 
    MGLevelObject<LevelMatrixType> mg_matrices; 

    LinearAlgebra::distributed::Vector<double> solution; 
    LinearAlgebra::distributed::Vector<double> system_rhs; 

    double             setup_time; 
    ConditionalOStream pcout; 
    ConditionalOStream time_details; 
  }; 

// 当我们初始化有限元时，我们当然也要使用文件顶部指定的度数（否则，在某些时候会抛出一个异常，因为在模板化的LaplaceOperator类中定义的计算内核和MatrixFree读出的有限元信息将不匹配）。三角形的构造函数需要设置一个额外的标志，告诉网格要符合顶点上的2:1单元平衡，这对于几何多网格例程的收敛是必需的。对于分布式网格，我们还需要特别启用多网格的层次结构。

  template <int dim> 
  LaplaceProblem<dim>::LaplaceProblem() 
    : 
#ifdef DEAL_II_WITH_P4EST 
    triangulation( 
      MPI_COMM_WORLD, 
      Triangulation<dim>::limit_level_difference_at_vertices, 
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy) 
    , 
#else 
    triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
    , 
#endif 
    fe(degree_finite_element) 
    , dof_handler(triangulation) 
    , setup_time(0.) 
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
    , 

// LaplaceProblem类拥有一个额外的输出流，用于收集关于设置阶段的详细时间信息。这个流被称为time_details，默认情况下通过这里指定的 @p false 参数被禁用。对于详细的时间，去掉 @p false 参数可以打印出所有的细节。

    time_details(std::cout, 
                 false && Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
  {} 

//  @sect4{LaplaceProblem::setup_system}  

// 设置阶段与 step-16 类似，由于LaplaceOperator类的存在而有相关的变化。首先要做的是设置DoFHandler，包括多网格层次的自由度，以及初始化悬挂节点的约束和同质二列条件。由于我们打算用MPI的%并行方式使用这个程序，我们需要确保约束条件能知道本地相关的自由度，否则在使用超过几亿个自由度的时候，存储会爆炸，见  step-40  。

// 一旦我们创建了多网格dof_handler和约束条件，我们就可以为全局矩阵算子以及多网格方案的每一层调用reinit函数。主要的操作是为问题设置 <code> MatrixFree </code> 实例。 <code>LaplaceOperator</code> 类的基类， MatrixFreeOperators::Base, 被初始化为一个指向MatrixFree对象的共享指针。这样，我们可以在这里简单地创建它，然后将它分别传递给系统矩阵和水平矩阵。为了设置MatrixFree，我们需要激活MatrixFree的AdditionalData字段中的更新标志，使其能够存储实空间中的正交点坐标（默认情况下，它只缓存梯度（反转置的雅各布）和JxW值的数据）。请注意，如果我们调用 reinit 函数而不指定级别（即给出  <code>level = numbers::invalid_unsigned_int</code>  ），MatrixFree 将在活动单元上构建一个循环。在本教程中，除了MPI之外，我们不使用线程，这就是为什么我们通过将 MatrixFree::AdditionalData::tasks_parallel_scheme 设置为 MatrixFree::AdditionalData::none. 来明确地禁用它 最后，系数被评估，向量被初始化，如上所述。

  template <int dim> 
  void LaplaceProblem<dim>::setup_system() 
  { 
    Timer time; 
    setup_time = 0; 

    system_matrix.clear(); 
    mg_matrices.clear_elements(); 

    dof_handler.distribute_dofs(fe); 
    dof_handler.distribute_mg_dofs(); 

    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
          << std::endl; 

    IndexSet locally_relevant_dofs; 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 

    constraints.clear(); 
    constraints.reinit(locally_relevant_dofs); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    VectorTools::interpolate_boundary_values( 
      mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints); 
    constraints.close(); 
    setup_time += time.wall_time(); 
    time_details << "Distribute DoFs & B.C.     (CPU/wall) " << time.cpu_time() 
                 << "s/" << time.wall_time() << "s" << std::endl; 
    time.restart(); 

    { 
      typename MatrixFree<dim, double>::AdditionalData additional_data; 
      additional_data.tasks_parallel_scheme = 
        MatrixFree<dim, double>::AdditionalData::none; 
      additional_data.mapping_update_flags = 
        (update_gradients | update_JxW_values | update_quadrature_points); 
      std::shared_ptr<MatrixFree<dim, double>> system_mf_storage( 
        new MatrixFree<dim, double>()); 
      system_mf_storage->reinit(mapping, 
                                dof_handler, 
                                constraints, 
                                QGauss<1>(fe.degree + 1), 
                                additional_data); 
      system_matrix.initialize(system_mf_storage); 
    } 

    system_matrix.evaluate_coefficient(Coefficient<dim>()); 

    system_matrix.initialize_dof_vector(solution); 
    system_matrix.initialize_dof_vector(system_rhs); 

    setup_time += time.wall_time(); 
    time_details << "Setup matrix-free system   (CPU/wall) " << time.cpu_time() 
                 << "s/" << time.wall_time() << "s" << std::endl; 
    time.restart(); 

// 接下来，初始化所有层次上的多网格方法的矩阵。数据结构MGConstrainedDoFs保留了受边界条件约束的指数信息，以及不同细化层次之间的边缘指数，如 step-16 教程程序中所述。然后，我们穿过网格的各个层次，在每个层次上构建约束和矩阵。这与原始网格上的系统矩阵的构造密切相关，只是在访问层级信息而不是活动单元的信息时，在命名上略有不同。

    const unsigned int nlevels = triangulation.n_global_levels(); 
    mg_matrices.resize(0, nlevels - 1); 

    std::set<types::boundary_id> dirichlet_boundary; 
    dirichlet_boundary.insert(0); 
    mg_constrained_dofs.initialize(dof_handler); 
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, 
                                                       dirichlet_boundary); 

    for (unsigned int level = 0; level < nlevels; ++level) 
      { 
        IndexSet relevant_dofs; 
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
                                                      level, 
                                                      relevant_dofs); 
        AffineConstraints<double> level_constraints; 
        level_constraints.reinit(relevant_dofs); 
        level_constraints.add_lines( 
          mg_constrained_dofs.get_boundary_indices(level)); 
        level_constraints.close(); 

        typename MatrixFree<dim, float>::AdditionalData additional_data; 
        additional_data.tasks_parallel_scheme = 
          MatrixFree<dim, float>::AdditionalData::none; 
        additional_data.mapping_update_flags = 
          (update_gradients | update_JxW_values | update_quadrature_points); 
        additional_data.mg_level = level; 
        std::shared_ptr<MatrixFree<dim, float>> mg_mf_storage_level( 
          new MatrixFree<dim, float>()); 
        mg_mf_storage_level->reinit(mapping, 
                                    dof_handler, 
                                    level_constraints, 
                                    QGauss<1>(fe.degree + 1), 
                                    additional_data); 

        mg_matrices[level].initialize(mg_mf_storage_level, 
                                      mg_constrained_dofs, 
                                      level); 
        mg_matrices[level].evaluate_coefficient(Coefficient<dim>()); 
      } 
    setup_time += time.wall_time(); 
    time_details << "Setup matrix-free levels   (CPU/wall) " << time.cpu_time() 
                 << "s/" << time.wall_time() << "s" << std::endl; 
  } 

//  @sect4{LaplaceProblem::assemble_rhs}  

// 组装函数非常简单，因为我们所要做的就是组装右侧。多亏了FEEvaluation和所有缓存在MatrixFree类中的数据，我们从 MatrixFreeOperators::Base, 中查询，这可以在几行中完成。由于这个调用没有被包裹到 MatrixFree::cell_loop 中（这将是一个替代方案），我们一定不要忘记在装配结束时调用compress()，将右手边的所有贡献发送给各自自由度的所有者。

  template <int dim> 
  void LaplaceProblem<dim>::assemble_rhs() 
  { 
    Timer time; 

    system_rhs = 0; 
    FEEvaluation<dim, degree_finite_element> phi( 
      *system_matrix.get_matrix_free()); 
    for (unsigned int cell = 0; 
         cell < system_matrix.get_matrix_free()->n_cell_batches(); 
         ++cell) 
      { 
        phi.reinit(cell); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          phi.submit_value(make_vectorized_array<double>(1.0), q); 
        phi.integrate(EvaluationFlags::values); 
        phi.distribute_local_to_global(system_rhs); 
      } 
    system_rhs.compress(VectorOperation::add); 

    setup_time += time.wall_time(); 
    time_details << "Assemble right hand side   (CPU/wall) " << time.cpu_time() 
                 << "s/" << time.wall_time() << "s" << std::endl; 
  } 

//  @sect4{LaplaceProblem::solve}  

// 解决的过程与  step-16  中类似。我们先从转移的设置开始。对于 LinearAlgebra::distributed::Vector, 来说，有一个非常快速的转移类，叫做MGTransferMatrixFree，它用FEEvaluation中同样的快速和因子化核在网格层之间进行插值。

  template <int dim> 
  void LaplaceProblem<dim>::solve() 
  { 
    Timer                            time; 
    MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs); 
    mg_transfer.build(dof_handler); 
    setup_time += time.wall_time(); 
    time_details << "MG build transfer time     (CPU/wall) " << time.cpu_time() 
                 << "s/" << time.wall_time() << "s\n"; 
    time.restart(); 

// 作为一个平滑器，本教程程序使用切比雪夫迭代，而不是 step-16 中的SOR。（SOR将很难实现，因为我们没有明确的矩阵元素，而且很难使其在%并行中有效工作）。 平滑器是用我们的水平矩阵和切比雪夫平滑器的强制性附加数据初始化的。我们在这里使用一个相对较高的度数（5），因为矩阵-向量乘积是比较便宜的。我们选择在平滑器中平滑出 $[1.2 \hat{\lambda}_{\max}/15,1.2 \hat{\lambda}_{\max}]$ 的范围，其中 $\hat{\lambda}_{\max}$ 是对最大特征值的估计（系数1.2在PreconditionChebyshev中应用）。为了计算该特征值，Chebyshev初始化执行了几步没有预处理的CG算法。由于最高的特征值通常是最容易找到的，而且一个粗略的估计就足够了，我们选择10次迭代。最后，我们还设置了切比雪夫方法中的内部预处理类型，这是一个雅可比迭代。这由DiagonalMatrix类来表示，该类得到了由我们的LaplaceOperator类提供的反对角线条目。

// 在第0层，我们以不同的方式初始化平滑器，因为我们想使用切比雪夫迭代作为求解器。PreconditionChebyshev允许用户切换到求解器模式，其中迭代次数在内部选择为正确值。在附加数据对象中，通过将多项式的度数选择为 @p numbers::invalid_unsigned_int. 来激活这一设置，然后算法将攻击粗级矩阵中最小和最大之间的所有特征值。切比雪夫平滑器的步数是这样选择的：切比雪夫收敛估计值保证将残差减少到变量 @p  smoothing_range中指定的数字。注意，对于求解来说， @p smoothing_range 是一个相对的公差，并且选择小于1，在这种情况下，我们选择三个数量级，而当只对选定的特征值进行平滑时，它是一个大于1的数字。

// 从计算的角度来看，只要粗粒度适中，Chebyshev迭代是一个非常有吸引力的粗粒度求解器。这是因为Chebyshev方法只执行矩阵-向量乘积和向量更新，这通常比其他迭代方法中涉及的内积更好地并行到有几万个核心的最大集群规模。前者只涉及到（粗）网格中邻居之间的局部通信，而后者则需要在所有处理器上进行全局通信。

    using SmootherType = 
      PreconditionChebyshev<LevelMatrixType, 
                            LinearAlgebra::distributed::Vector<float>>; 
    mg::SmootherRelaxation<SmootherType, 
                           LinearAlgebra::distributed::Vector<float>> 
                                                         mg_smoother; 
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data; 
    smoother_data.resize(0, triangulation.n_global_levels() - 1); 
    for (unsigned int level = 0; level < triangulation.n_global_levels(); 
         ++level) 
      { 
        if (level > 0) 
          { 
            smoother_data[level].smoothing_range     = 15.; 
            smoother_data[level].degree              = 5; 
            smoother_data[level].eig_cg_n_iterations = 10; 
          } 
        else 
          { 
            smoother_data[0].smoothing_range = 1e-3; 
            smoother_data[0].degree          = numbers::invalid_unsigned_int; 
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m(); 
          } 
        mg_matrices[level].compute_diagonal(); 
        smoother_data[level].preconditioner = 
          mg_matrices[level].get_matrix_diagonal_inverse(); 
      } 
    mg_smoother.initialize(mg_matrices, smoother_data); 

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>> 
      mg_coarse; 
    mg_coarse.initialize(mg_smoother); 

// 下一步是设置悬挂节点情况下所需的接口矩阵。deal.II中的自适应多网格实现了一种叫做局部平滑的方法。这意味着最细级别的平滑只覆盖固定（最细）网格级别所定义的网格的局部部分，而忽略了计算域中终端单元比该级别更粗的部分。随着该方法向更粗的级别发展，越来越多的全局网格将被覆盖。在某个更粗的层次上，整个网格将被覆盖。由于多网格方法中的所有层次矩阵都覆盖了网格中的单一层次，所以在层次矩阵上不会出现悬空节点。在多网格层之间的界面上，在平滑的同时设置同质Dirichlet边界条件。然而，当残差被转移到下一个更粗的层次时，需要考虑到多网格界面的耦合。这是由所谓的界面（或边缘）矩阵来完成的，它计算了被具有同质Dirichlet条件的层次矩阵所遗漏的残差部分。我们参考 @ref mg_paper "Janssen和Kanschat的多网格论文 "以了解更多细节。

// 对于这些接口矩阵的实现，已经有一个预定义的类 MatrixFreeOperators::MGInterfaceOperator ，它将例程 MatrixFreeOperators::Base::vmult_interface_down() 和 MatrixFreeOperators::Base::vmult_interface_up() 包装在一个带有 @p vmult()和 @p Tvmult() 操作（最初是为矩阵编写的，因此期待这些名字）的新类中。请注意，vmult_interface_down是在多网格V周期的限制阶段使用的，而vmult_interface_up是在延长阶段使用的。

// 一旦接口矩阵被创建，我们完全按照 step-16 的方法设置剩余的多网格预处理基础设施，以获得一个可以应用于矩阵的 @p preconditioner 对象。

    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix( 
      mg_matrices); 

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>> 
      mg_interface_matrices; 
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1); 
    for (unsigned int level = 0; level < triangulation.n_global_levels(); 
         ++level) 
      mg_interface_matrices[level].initialize(mg_matrices[level]); 
    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface( 
      mg_interface_matrices); 

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg( 
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother); 
    mg.set_edge_matrices(mg_interface, mg_interface); 

    PreconditionMG<dim, 
                   LinearAlgebra::distributed::Vector<float>, 
                   MGTransferMatrixFree<dim, float>> 
      preconditioner(dof_handler, mg, mg_transfer); 

// 多网格程序的设置非常简单，与  step-16  相比，在求解过程中看不出有什么不同。所有的魔法都隐藏在  LaplaceOperator::vmult  操作的实现背后。请注意，我们通过标准输出打印出求解时间和累积的设置时间，也就是说，在任何情况下，而设置操作的详细时间只在构造函数中的detail_times标志被改变的情况下打印。

    SolverControl solver_control(100, 1e-12 * system_rhs.l2_norm()); 
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control); 
    setup_time += time.wall_time(); 
    time_details << "MG build smoother time     (CPU/wall) " << time.cpu_time() 
                 << "s/" << time.wall_time() << "s\n"; 
    pcout << "Total setup time               (wall) " << setup_time << "s\n"; 

    time.reset(); 
    time.start(); 
    constraints.set_zero(solution); 
    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    constraints.distribute(solution); 

    pcout << "Time solve (" << solver_control.last_step() << " iterations)" 
          << (solver_control.last_step() < 10 ? "  " : " ") << "(CPU/wall) " 
          << time.cpu_time() << "s/" << time.wall_time() << "s\n"; 
  } 

//  @sect4{LaplaceProblem::output_results}  

// 这里是数据输出，是  step-5  的简化版本。我们对细化过程中产生的每个网格使用标准的VTU（=压缩的VTK）输出。此外，我们还使用了一种针对速度而不是磁盘使用量进行优化的压缩算法。默认设置（针对磁盘使用进行优化）使得保存输出的时间是运行线性求解器的4倍，而将 DataOutBase::VtkFlags::compression_level 设置为 DataOutBase::VtkFlags::best_speed 则将其降低到只有线性求解的四分之一的时间。

// 当网格过大时，我们禁用输出。这个程序的一个变种已经在几十万个MPI行列上运行，网格单元多达1000亿个，经典的可视化工具无法直接访问。

  template <int dim> 
  void LaplaceProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    Timer time; 
    if (triangulation.n_global_active_cells() > 1000000) 
      return; 

    DataOut<dim> data_out; 

    solution.update_ghost_values(); 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(mapping); 

    DataOutBase::VtkFlags flags; 
    flags.compression_level = DataOutBase::VtkFlags::best_speed; 
    data_out.set_flags(flags); 
    data_out.write_vtu_with_pvtu_record( 
      "./", "solution", cycle, MPI_COMM_WORLD, 3); 

    time_details << "Time write output          (CPU/wall) " << time.cpu_time() 
                 << "s/" << time.wall_time() << "s\n"; 
  } 

//  @sect4{LaplaceProblem::run}  

// 运行该程序的函数与  step-16  中的函数非常相似。与2D相比，我们在3D中做了很少的细化步骤，但仅此而已。

// 在运行程序之前，我们先输出一些关于检测到的矢量化水平的信息，正如在介绍中所讨论的那样。

  template <int dim> 
  void LaplaceProblem<dim>::run() 
  { 
    { 
      const unsigned int n_vect_doubles = VectorizedArray<double>::size(); 
      const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles; 

      pcout << "Vectorization over " << n_vect_doubles 
            << " doubles = " << n_vect_bits << " bits (" 
            << Utilities::System::get_current_vectorization_level() << ")" 
            << std::endl; 
    } 

    for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle) 
      { 
        pcout << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_cube(triangulation, 0., 1.); 
            triangulation.refine_global(3 - dim); 
          } 
        triangulation.refine_global(1); 
        setup_system(); 
        assemble_rhs(); 
        solve(); 
        output_results(cycle); 
        pcout << std::endl; 
      }; 
  } 
} // namespace Step37 

//  @sect3{The <code>main</code> function}  

// 除了我们根据 step-40 设置了MPI框架外，主函数中没有任何意外。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace Step37; 

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1); 

      LaplaceProblem<dimension> laplace_problem; 
      laplace_problem.run(); 
    } 
  catch (std::exception &exc) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Exception on processing: " << std::endl 
                << exc.what() << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      return 1; 
    } 
  catch (...) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Unknown exception!" << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      return 1; 
    } 

  return 0; 
} 

