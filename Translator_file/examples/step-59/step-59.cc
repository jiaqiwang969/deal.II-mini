

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2018 - 2021 by the deal.II authors 
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
 * Authors: Katharina Kormann, Martin Kronbichler, 2018 
 */ 



// 包含的文件与  step-37  中的基本相同，只是用有限元类FE_DGQHermite代替了FE_Q。所有对面积分进行无矩阵计算的功能已经包含在`fe_evaluation.h`中。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/timer.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/la_parallel_vector.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/tensor_product_matrix.h> 

#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_tools.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 

#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer_matrix_free.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_matrix.h> 

#include <deal.II/numerics/vector_tools.h> 

#include <deal.II/matrix_free/matrix_free.h> 
#include <deal.II/matrix_free/fe_evaluation.h> 

#include <iostream> 
#include <fstream> 

namespace Step59 
{ 
  using namespace dealii; 

// 和 step-37 一样，为了简单起见，我们在程序顶部将维数和多项式程度收集为常数。与 step-37 不同的是，这次我们选择了一个真正的高阶方法，度数为8，任何不使用和因式分解的实现都会变得非常慢，而使用MatrixFree的实现则提供了与度数为2或3时基本相同的效率。此外，本教程程序中的所有类都是模板化的，因此，通过在`main()`函数中添加适当度数的实例，可以很容易地在运行时从输入文件或命令行参数中选择度数。

  const unsigned int degree_finite_element = 8; 
  const unsigned int dimension             = 3; 
// @sect3{Equation data}  

// 与 step-7 相类似，我们定义了一个分析解，我们试图用离散化重现这个分析解。由于本教程的目的是展示无矩阵方法，我们选择了一个最简单的可能性，即一个余弦函数，其导数对我们来说足够简单，可以通过分析计算。再往下看，我们在这里选择的波数2.4将与 $x$ -方向的域范围即2.5相匹配，这样我们在 $x = 2.5$ 得到一个周期性的解，包括 $6pi$ 或余弦的三个整波转。第一个函数定义了解和它的梯度，分别用于表达Dirichlet和Neumann边界条件的解析解。此外，一个代表解的负拉普拉斯的类被用来表示右手边（强制）函数，我们用它来匹配离散化版本中的给定分析解（制造解）。

  template <int dim> 
  class Solution : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int = 0) const override final 
    { 
      double val = 1.; 
      for (unsigned int d = 0; d < dim; ++d) 
        val *= std::cos(numbers::PI * 2.4 * p[d]); 
      return val; 
    } 

    virtual Tensor<1, dim> gradient(const Point<dim> &p, 
                                    const unsigned int = 0) const override final 
    { 
      const double   arg = numbers::PI * 2.4; 
      Tensor<1, dim> grad; 
      for (unsigned int d = 0; d < dim; ++d) 
        { 
          grad[d] = 1.; 
          for (unsigned int e = 0; e < dim; ++e) 
            if (d == e) 
              grad[d] *= -arg * std::sin(arg * p[e]); 
            else 
              grad[d] *= std::cos(arg * p[e]); 
        } 
      return grad; 
    } 
  }; 

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int = 0) const override final 
    { 
      const double arg = numbers::PI * 2.4; 
      double       val = 1.; 
      for (unsigned int d = 0; d < dim; ++d) 
        val *= std::cos(arg * p[d]); 
      return dim * arg * arg * val; 
    } 
  }; 

//  @sect3{Matrix-free implementation}  

// `LaplaceOperator`类与  step-37  中的相应类类似。一个重要的区别是，我们没有从  MatrixFreeOperators::Base  派生出这个类，因为我们想呈现  MatrixFree::loop()  的一些额外特性，这些特性在通用类  MatrixFreeOperators::Base.  中是不可用的。我们从Subscriptor类派生出这个类，以便能够在Chebyshev预处理程序中使用该操作符，因为该预处理程序通过SmartPointer存储基础矩阵。

// 鉴于我们手工实现了一个完整的矩阵接口，我们需要添加一个`initialize()`函数，一个`m()`函数，一个`vmult()`函数和一个`Tvmult()`函数，这些都是之前由  MatrixFreeOperators::Base.  我们的LaplaceOperator还包含一个成员函数`get_penalty_factor()`，根据  step-39  集中选择对称内部惩罚方法中的惩罚参数 。

  template <int dim, int fe_degree, typename number> 
  class LaplaceOperator : public Subscriptor 
  { 
  public: 
    using value_type = number; 

    LaplaceOperator() = default; 

    void initialize(std::shared_ptr<const MatrixFree<dim, number>> data); 

    void clear(); 

    types::global_dof_index m() const; 

    void initialize_dof_vector( 
      LinearAlgebra::distributed::Vector<number> &vec) const; 

    std::shared_ptr<const MatrixFree<dim, number>> get_matrix_free() const; 

    void vmult(LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src) const; 

    void Tvmult(LinearAlgebra::distributed::Vector<number> &      dst, 
                const LinearAlgebra::distributed::Vector<number> &src) const; 

    number get_penalty_factor() const 
    { 
      return 1.0 * fe_degree * (fe_degree + 1); 
    } 

  private: 
    void 
    apply_cell(const MatrixFree<dim, number> &                   data, 
               LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src, 
               const std::pair<unsigned int, unsigned int> &cell_range) const; 

    void 
    apply_face(const MatrixFree<dim, number> &                   data, 
               LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src, 
               const std::pair<unsigned int, unsigned int> &face_range) const; 

    void apply_boundary( 
      const MatrixFree<dim, number> &                   data, 
      LinearAlgebra::distributed::Vector<number> &      dst, 
      const LinearAlgebra::distributed::Vector<number> &src, 
      const std::pair<unsigned int, unsigned int> &     face_range) const; 

    std::shared_ptr<const MatrixFree<dim, number>> data; 
  }; 

// `%PreconditionBlockJacobi`类定义了我们对这个问题的自定义预处理程序。与基于矩阵对角线的 step-37 不同，我们在这里通过使用介绍中讨论的所谓快速对角线化方法来计算非连续Galerkin方法中对角线块的近似反演。

  template <int dim, int fe_degree, typename number> 
  class PreconditionBlockJacobi 
  { 
  public: 
    using value_type = number; 

    void clear() 
    { 
      cell_matrices.clear(); 
    } 

    void initialize(const LaplaceOperator<dim, fe_degree, number> &op); 

    void vmult(LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src) const; 

    void Tvmult(LinearAlgebra::distributed::Vector<number> &      dst, 
                const LinearAlgebra::distributed::Vector<number> &src) const 
    { 
      vmult(dst, src); 
    } 

  private: 
    std::shared_ptr<const MatrixFree<dim, number>> data; 
    std::vector<TensorProductMatrixSymmetricSum<dim, 
                                                VectorizedArray<number>, 
                                                fe_degree + 1>> 
      cell_matrices; 
  }; 

//这个独立的函数在`LaplaceOperator'和`%PreconditionBlockJacobi'类中都被用来调整鬼魂范围。这个函数是必要的，因为`vmult()`函数所提供的一些向量没有用包括正确的鬼魂条目布局的 `LaplaceOperator::initialize_dof_vector` 来正确初始化，而是来自MGTransferMatrixFree类，该类对无矩阵类的鬼魂选择没有概念。为了避免索引混乱，我们必须在对这些向量进行实际操作之前调整鬼域。由于向量在多网格平滑器和传输类中被保留下来，一个曾经被调整过重影范围的向量在对象的整个生命周期中都会保持这种状态，所以我们可以在函数的开始使用一个快捷方式来查看分布式向量的分区器对象（以共享指针的形式存储）是否与MatrixFree所期望的布局相同，它被存储在一个由 MatrixFree::get_dof_info(0),  访问的数据结构中 ]，其中的0表示从中提取的DoFHandler编号；我们在MatrixFree中只使用一个DoFHandler，所以这里唯一有效的编号是0。

  template <int dim, typename number> 
  void adjust_ghost_range_if_necessary( 
    const MatrixFree<dim, number> &                   data, 
    const LinearAlgebra::distributed::Vector<number> &vec) 
  { 
    if (vec.get_partitioner().get() == 
        data.get_dof_info(0).vector_partitioner.get()) 
      return; 

    LinearAlgebra::distributed::Vector<number> copy_vec(vec); 
    const_cast<LinearAlgebra::distributed::Vector<number> &>(vec).reinit( 
      data.get_dof_info(0).vector_partitioner); 
    const_cast<LinearAlgebra::distributed::Vector<number> &>(vec) 
      .copy_locally_owned_data_from(copy_vec); 
  } 

// 接下来的五个函数用于清除和初始化`LaplaceOperator`类，返回持有MatrixFree数据容器的共享指针，以及正确初始化向量和运算符大小，与 step-37 或者说 MatrixFreeOperators::Base. 的内容相同。
  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::clear() 
  { 
    data.reset(); 
  } 

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::initialize( 
    std::shared_ptr<const MatrixFree<dim, number>> data) 
  { 
    this->data = data; 
  } 

  template <int dim, int fe_degree, typename number> 
  std::shared_ptr<const MatrixFree<dim, number>> 
  LaplaceOperator<dim, fe_degree, number>::get_matrix_free() const 
  { 
    return data; 
  } 

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::initialize_dof_vector( 
    LinearAlgebra::distributed::Vector<number> &vec) const 
  { 
    data->initialize_dof_vector(vec); 
  } 

  template <int dim, int fe_degree, typename number> 
  types::global_dof_index LaplaceOperator<dim, fe_degree, number>::m() const 
  { 
    Assert(data.get() != nullptr, ExcNotInitialized()); 
    return data->get_dof_handler().n_dofs(); 
  } 

// 这个函数在向量`src`上实现了LaplaceOperator的动作，并将结果存储在向量`dst`中。与 step-37 相比，这个调用有四个新特性。

// 第一个新特性是上面提到的`adjust_ghost_range_if_necessary`函数，该函数需要使向量符合单元和面函数中FEEvaluation和FEFaceEvaluation所期望的布局。

// 第二个新特征是我们没有像 step-37 中那样实现`vmult_add()`函数（通过虚拟函数 MatrixFreeOperators::Base::vmult_add()), ，而是直接实现`vmult()`功能。由于单元和面的积分都将和到目的向量中，我们当然必须在某处将向量归零。对于DG元素，我们有两个选择&ndash；一个是使用 FEEvaluation::set_dof_values() 而不是下面`apply_cell`函数中的 FEEvaluation::distribute_local_to_global() 。这是因为MatrixFree中的循环布局是这样的：单元积分总是在面积分之前接触到一个给定的向量条目。然而，这实际上只适用于完全不连续的基数，其中每个单元都有自己的自由度，不与邻近的结果共享。另一种设置，即这里选择的设置，是让 MatrixFree::loop() 来处理向量的归零问题。这可以被认为是在代码的某个地方简单地调用`dst = 0;`。对于像 `LinearAlgebra::distributed::Vector`, 这样的支持性向量来说，实现起来就比较麻烦了，因为我们的目标是不要一次性将整个向量清零。在足够小的几千个向量项上进行归零操作的好处是，在 FEEvaluation::distribute_local_to_global() 和 FEFaceEvaluation::distribute_local_to_global(). 中再次访问之前，被归零的向量项会保留在缓存中，因为无矩阵运算符的评估真的很快，仅仅归零一个大的向量就会相当于运算符评估时间的25%，我们显然希望避免这种代价。对于 MatrixFree::cell_loop 和连续基数来说，也可以使用这种将向量归零的选项，尽管在 step-37 或 step-48 的教程程序中没有使用它。

// 第三个新特征是我们提供了在单元格、内面和边界面进行计算的函数方式。MatrixFree类有一个叫做`loop`的函数，它接收三个函数指针，用于三种情况，允许分开实现不同的东西。正如在 step-37 中所解释的，这些函数指针可以是 `std::function` 对象或类的成员函数。在这种情况下，我们使用指向成员函数的指针。

// 最后的新特征是可以给 MatrixFree::DataAccessOnFaces 类型的最后两个参数，这个类将面积分的数据访问类型传递给并行向量的MPI数据交换例程 LinearAlgebra::distributed::Vector::update_ghost_values() 和 LinearAlgebra::distributed::Vector::compress() 。其目的是不发送相邻元素的所有自由度，而是将数据量减少到手头计算真正需要的程度。数据交换是一个真正的瓶颈，特别是对于高自由度的DG方法来说，因此一个更严格的交换方式往往是有益的。枚举字段 MatrixFree::DataAccessOnFaces 可以取值`none`，这意味着根本不做面的积分，这类似于 MatrixFree::cell_loop(), 的值`values`，意味着只使用面的形状函数值（但不使用导数），而值`gradients`时，除了值之外还访问面的第一导数。值`unspecified`意味着所有的自由度将被交换给位于处理器边界的面，并指定在本地处理器上进行处理。

// 为了了解数据是如何被减少的，想想节点元素FE_DGQ的情况，节点点在元素表面，在一个单元的 $(k+1)^d$ 个自由度中，只有 $(k+1)^{d-1}$ 个自由度对 $d$ 个空间维度的多项式程度的面的值有贡献。类似的减少也可以用于内部惩罚方法，该方法评估面的值和一导数。当在一维中使用类Hermite基时，最多只有两个基函数对数值和导数有贡献。FE_DGQHermite类实现了这一概念的张量乘积，在介绍中已经讨论过。因此，每个面只需交换 $2(k+1)^{d-1}$ 个自由度，一旦 $k$ 个自由度大于4或5个，这显然是一种胜利。请注意，FE_DGQHermite的这种减少的交换在具有弯曲边界的网格上也是有效的，因为导数是在参考元素上取的，而几何体只在内部混合它们。因此，这与试图用连续的Hermite型形状函数获得 $C^1$ 的连续性不同，在这种情况下，非笛卡尔的情况会大大改变情况。显然，在非笛卡尔网格上，导数还包括超出法向导数的形状函数的切向导数，但这些也只需要元素表面的函数值。如果元素不提供任何压缩，循环会自动交换受影响单元的所有条目。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::vmult( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    adjust_ghost_range_if_necessary(*data, dst); 
    adjust_ghost_range_if_necessary(*data, src); 
    data->loop(&LaplaceOperator::apply_cell, 
               &LaplaceOperator::apply_face, 
               &LaplaceOperator::apply_boundary, 
               this, 
               dst, 
               src, 
               /* zero_dst =  */ true,
               MatrixFree<dim, number>::DataAccessOnFaces::gradients,
               MatrixFree<dim, number>::DataAccessOnFaces::gradients);
  } 

// 由于拉普拉斯是对称的，`Tvmult()`（多网格平滑界面需要）操作被简单地转发给`vmult()`的情况。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::Tvmult( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    vmult(dst, src); 
  } 

//单元格的操作与 step-37 非常相似。不过我们在这里没有使用系数。第二个区别是，我们用一个单一的函数调用 FEEvaluation::gather_evaluate() 代替了 FEEvaluation::read_dof_values() 后面的 FEEvaluation::evaluate() 这两个步骤，在内部调用这两个单独方法的序列。同样， FEEvaluation::integrate_scatter() 实现了 FEEvaluation::integrate() 之后的 FEEvaluation::distribute_local_to_global(). 的序列。 在这种情况下，这些新函数只是节省了两行代码。然而，我们用它们来与FEFaceEvaluation进行类比，在那里它们更重要，如下所述。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::apply_cell( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     cell_range) const 
  { 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data); 
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        phi.reinit(cell); 
        phi.gather_evaluate(src, EvaluationFlags::gradients); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          phi.submit_gradient(phi.get_gradient(q), q); 
        phi.integrate_scatter(EvaluationFlags::gradients, dst); 
      } 
  } 

// 面部操作实现了与 step-39 类似的内部惩罚方法的条款，正如介绍中所解释的。我们需要两个评估器对象来完成这个任务，一个用于处理来自内部面的两边之一的单元格的解，另一个用于处理来自另一边的解。面积分的评价器被称为FEFaceEvaluation，并在构造函数的第二个槽中接受一个布尔参数，以指示评价器应属于两边中的哪一边。在FEFaceEvaluation和MatrixFree中，我们称两边中的一边为 "内部"，另一边为 "外部"。`外部'这个名字是指两边的评价器将返回相同的法向量。对于 "内部 "一侧，法向量指向外部，而另一侧则指向内部，并且与该单元的外部法向量相对应。除了新的类名之外，我们再次得到了一系列的项目，与 step-37 中讨论的类似，但在这种情况下是针对内部面的。请注意，MatrixFree的数据结构形成了面的批次，类似于单元积分的单元批次。一批中的所有面涉及不同的单元格编号，但在参考单元格中具有相同的面编号，具有相同的细化配置（无细化或相同的子面）和相同的方向，以保持SIMD操作的简单和高效。

// 注意，除了法线方向的逻辑决定外，内部与外部没有任何隐含的意义，这在内部是相当随机的。我们绝对不能依赖分配内部与外部标志的某种模式，因为这个决定是为了MatrixFree设置例程中的访问规则性和统一性而做出的。由于大多数正常的DG方法都是保守的，也就是说，通量在接口的两边看起来都是一样的，所以如果内部/外部标志被调换，法线向量得到相反的符号，那么数学就不会有任何改变。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::apply_face( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     face_range) const 
  { 
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_inner(data, 
                                                                         true); 
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_outer(data, 
                                                                         false); 
    for (unsigned int face = face_range.first; face < face_range.second; ++face) 
      { 

// 在给定的一批面孔上，我们首先更新指向当前面孔的指针，然后访问矢量。如上所述，我们把访问向量和评估结合起来。在面积分的情况下，对于FE_DGQHermite基础的特殊情况，可以减少对向量的数据访问，正如上面解释的数据交换。由于 $2(k+1)^{d-1}$ 个单元自由度中只有 $(k+1)^d$ 个单元自由度被非零值或形状函数的导数所乘，这种结构可以被用于评估，大大减少了数据访问。减少数据访问不仅是有益的，因为它减少了飞行中的数据，从而有助于缓存，而且当从单元格索引列表中相距较远的单元格收集数值时，对面的数据访问往往比单元格积分更不规则。

        phi_inner.reinit(face); 
        phi_inner.gather_evaluate(src, 
                                  EvaluationFlags::values | 
                                    EvaluationFlags::gradients); 
        phi_outer.reinit(face); 
        phi_outer.gather_evaluate(src, 
                                  EvaluationFlags::values | 
                                    EvaluationFlags::gradients); 

// 接下来的两个语句是计算内部惩罚法的惩罚参数。正如在介绍中所解释的，我们希望有一个像 $\frac{1}{h_\text{i}}$ 这样的长度 $h_\text{i}$ 法线到面的缩放比例。对于一般的非笛卡尔网格，这个长度必须由反雅各布系数乘以实坐标的法向量的乘积来计算。从这个 "dim "分量的向量中，我们必须最终挑选出与参考单元的法线方向一致的分量。在MatrixFree中存储的几何数据中，雅各布式中的分量被应用，使得后一个方向总是最后一个分量`dim-1`（这很有利，因为参考单元的导数排序可以与面的方向无关）。这意味着我们可以简单地访问最后一个分量`dim-1`，而不必在`data.get_face_info(face).internal_face_no`和`data.get_face_info(face).exterior_face_no`中查找局部面的编号。最后，我们还必须取这些因素的绝对值，因为法线可能指向正或负的方向。

        const VectorizedArray<number> inverse_length_normal_to_face = 
          0.5 * (std::abs((phi_inner.get_normal_vector(0) * 
                           phi_inner.inverse_jacobian(0))[dim - 1]) + 
                 std::abs((phi_outer.get_normal_vector(0) * 
                           phi_outer.inverse_jacobian(0))[dim - 1])); 
        const VectorizedArray<number> sigma = 
          inverse_length_normal_to_face * get_penalty_factor(); 

// 在正交点的循环中，我们最终计算了对内部惩罚方案的所有贡献。根据介绍中的公式，测试函数的值被乘以解决方案中的跳跃乘以惩罚参数和实空间中的法向导数的平均值的差值。由于内侧和外侧的两个评估器由于跳跃而得到不同的符号，我们在这里用不同的符号传递结果。测试函数的正态导数会被内侧和外侧的解决方案中的负跳跃所乘。这个术语，被称为邻接一致性术语，根据其与原始一致性术语的关系，在代码中还必须包括 $\frac{1}{2}$ 的系数，由于测试函数槽中的平均数，它得到了二分之一的系数。

        for (unsigned int q = 0; q < phi_inner.n_q_points; ++q) 
          { 
            const VectorizedArray<number> solution_jump = 
              (phi_inner.get_value(q) - phi_outer.get_value(q)); 
            const VectorizedArray<number> average_normal_derivative = 
              (phi_inner.get_normal_derivative(q) + 
               phi_outer.get_normal_derivative(q)) * 
              number(0.5); 
            const VectorizedArray<number> test_by_value = 
              solution_jump * sigma - average_normal_derivative; 

            phi_inner.submit_value(test_by_value, q); 
            phi_outer.submit_value(-test_by_value, q); 

            phi_inner.submit_normal_derivative(-solution_jump * number(0.5), q); 
            phi_outer.submit_normal_derivative(-solution_jump * number(0.5), q); 
          } 

// 一旦我们完成了正交点的循环，我们就可以对面的积分循环进行和因子化操作，并将结果加到结果向量中，使用`integrate_scatter`函数。`scatter'这个名字反映了使用与`gather_evaluate'相同的模式将矢量数据分布到矢量中的分散位置。像以前一样，整合+写操作的组合允许我们减少数据访问。

        phi_inner.integrate_scatter(EvaluationFlags::values | 
                                      EvaluationFlags::gradients, 
                                    dst); 
        phi_outer.integrate_scatter(EvaluationFlags::values | 
                                      EvaluationFlags::gradients, 
                                    dst); 
      } 
  } 

// 边界面函数大体上沿用了内部面函数。唯一的区别是，我们没有一个单独的FEFaceEvaluation对象为我们提供外部值  $u^+$  ，但我们必须从边界条件和内部值  $u^-$  来定义它们。正如介绍中所解释的，我们在Dirichlet边界上使用 $u^+ = -u^- + 2 g_\text{D}$ 和 $\mathbf{n}^-\cdot \nabla u^+ = \mathbf{n}^-\cdot \nabla u^-$ ，在Neumann边界上使用 $u^+=u^-$ 和 $\mathbf{n}^-\cdot \nabla u^+ = -\mathbf{n}^-\cdot \nabla u^- + 2 g_\text{N}$  。由于这个操作实现了同质部分，即矩阵-向量乘积，我们必须在这里忽略边界函数 $g_\text{D}$ 和 $g_\text{N}$ ，并在 `LaplaceProblem::compute_rhs()`. 中把它们加到右侧。 ] 注意，由于通过 $u^+$ 将解 $u^-$ 扩展到外部，我们可以保持所有因子 $0.5$ 与内面函数相同，也可参见 step-39 中的讨论。

// 在这一点上有一个问题。下面的实现使用一个布尔变量`is_dirichlet`来切换Dirichlet和Neumann情况。然而，我们解决了一个问题，我们还想在一些边界上施加周期性的边界条件，即沿着 $x$ 方向的边界。人们可能会问，这里应该如何处理这些条件。答案是MatrixFree会自动将周期性边界视为技术上的边界，即两个相邻单元的解值相遇的内面，必须用适当的数值通量来处理。因此，周期性边界上的所有面将出现在`apply_face()`函数中，而不是这个函数中。

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::apply_boundary( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     face_range) const 
  { 
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_inner(data, 
                                                                         true); 
    for (unsigned int face = face_range.first; face < face_range.second; ++face) 
      { 
        phi_inner.reinit(face); 
        phi_inner.gather_evaluate(src, 
                                  EvaluationFlags::values | 
                                    EvaluationFlags::gradients); 

        const VectorizedArray<number> inverse_length_normal_to_face = 
          std::abs((phi_inner.get_normal_vector(0) * 
                    phi_inner.inverse_jacobian(0))[dim - 1]); 
        const VectorizedArray<number> sigma = 
          inverse_length_normal_to_face * get_penalty_factor(); 

        const bool is_dirichlet = (data.get_boundary_id(face) == 0); 

        for (unsigned int q = 0; q < phi_inner.n_q_points; ++q) 
          { 
            const VectorizedArray<number> u_inner = phi_inner.get_value(q); 
            const VectorizedArray<number> u_outer = 
              is_dirichlet ? -u_inner : u_inner; 
            const VectorizedArray<number> normal_derivative_inner = 
              phi_inner.get_normal_derivative(q); 
            const VectorizedArray<number> normal_derivative_outer = 
              is_dirichlet ? normal_derivative_inner : -normal_derivative_inner; 
            const VectorizedArray<number> solution_jump = (u_inner - u_outer); 
            const VectorizedArray<number> average_normal_derivative = 
              (normal_derivative_inner + normal_derivative_outer) * number(0.5); 
            const VectorizedArray<number> test_by_value = 
              solution_jump * sigma - average_normal_derivative; 
            phi_inner.submit_normal_derivative(-solution_jump * number(0.5), q); 
            phi_inner.submit_value(test_by_value, q); 
          } 
        phi_inner.integrate_scatter(EvaluationFlags::values | 
                                      EvaluationFlags::gradients, 
                                    dst); 
      } 
  } 

// 接下来我们来看看预处理程序的初始化。正如介绍中所解释的，我们想从一维质量和拉普拉斯矩阵的乘积中构造一个（近似的）单元矩阵的逆。我们的首要任务是计算一维矩阵，我们通过首先创建一个一维有限元来实现。在这里，我们没有预见到FE_DGQHermite<1>，而是从DoFHandler获得有限元的名称，用1替换 @p dim 参数（2或3）来创建一个一维名称，并通过使用FETools来构造一维元素。

  template <int dim, int fe_degree, typename number> 
  void PreconditionBlockJacobi<dim, fe_degree, number>::initialize( 
    const LaplaceOperator<dim, fe_degree, number> &op) 
  { 
    data = op.get_matrix_free(); 

    std::string name = data->get_dof_handler().get_fe().get_name(); 
    name.replace(name.find('<') + 1, 1, "1"); 
    std::unique_ptr<FiniteElement<1>> fe_1d = FETools::get_fe_by_name<1>(name); 

// 至于在单位元素上计算一维矩阵，我们简单地写下在矩阵的行和列以及正交点上的典型装配程序会做什么。我们一劳永逸地选择相同的拉普拉斯矩阵，对内部面使用系数0.5（但可能由于网格的原因，在不同方向上的缩放比例不同）。因此，我们在Dirichlet边界（正确的系数应该是导数项为1，惩罚项为2，见 step-39 ）或在Neumann边界（系数应该为0）犯了一个小错误。由于我们只在多网格方案中使用这个类作为平滑器，这个错误不会有任何重大影响，只是影响平滑质量。

    const unsigned int                                 N = fe_degree + 1; 
    FullMatrix<double>                                 laplace_unscaled(N, N); 
    std::array<Table<2, VectorizedArray<number>>, dim> mass_matrices; 
    std::array<Table<2, VectorizedArray<number>>, dim> laplace_matrices; 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        mass_matrices[d].reinit(N, N); 
        laplace_matrices[d].reinit(N, N); 
      } 

    QGauss<1> quadrature(N); 
    for (unsigned int i = 0; i < N; ++i) 
      for (unsigned int j = 0; j < N; ++j) 
        { 
          double sum_mass = 0, sum_laplace = 0; 
          for (unsigned int q = 0; q < quadrature.size(); ++q) 
            { 
              sum_mass += (fe_1d->shape_value(i, quadrature.point(q)) * 
                           fe_1d->shape_value(j, quadrature.point(q))) * 
                          quadrature.weight(q); 
              sum_laplace += (fe_1d->shape_grad(i, quadrature.point(q))[0] * 
                              fe_1d->shape_grad(j, quadrature.point(q))[0]) * 
                             quadrature.weight(q); 
            } 
          for (unsigned int d = 0; d < dim; ++d) 
            mass_matrices[d](i, j) = sum_mass; 

// 接下来两个语句组装的左右边界项似乎有一些任意的符号，但这些都是正确的，可以通过查看 step-39 并插入1D情况下的法向量的值-1和1来验证。

          sum_laplace += 
            (1. * fe_1d->shape_value(i, Point<1>()) * 
               fe_1d->shape_value(j, Point<1>()) * op.get_penalty_factor() + 
             0.5 * fe_1d->shape_grad(i, Point<1>())[0] * 
               fe_1d->shape_value(j, Point<1>()) + 
             0.5 * fe_1d->shape_grad(j, Point<1>())[0] * 
               fe_1d->shape_value(i, Point<1>())); 

          sum_laplace += 
            (1. * fe_1d->shape_value(i, Point<1>(1.0)) * 
               fe_1d->shape_value(j, Point<1>(1.0)) * op.get_penalty_factor() - 
             0.5 * fe_1d->shape_grad(i, Point<1>(1.0))[0] * 
               fe_1d->shape_value(j, Point<1>(1.0)) - 
             0.5 * fe_1d->shape_grad(j, Point<1>(1.0))[0] * 
               fe_1d->shape_value(i, Point<1>(1.0))); 

          laplace_unscaled(i, j) = sum_laplace; 
        } 

// 接下来，我们通过单元格，将缩放后的矩阵传递给TensorProductMatrixSymmetricSum，以实际计算代表逆的广义特征值问题。由于矩阵的近似构造为 $A\otimes M + M\otimes A$ ，并且每个元素的权重是恒定的，我们可以在拉普拉斯矩阵上应用所有的权重，并简单地保持质量矩阵不被缩放。在单元格的循环中，我们要利用MatrixFree类提供的几何体压缩，并检查当前的几何体是否与上一批单元格上的几何体相同，在这种情况下就没有什么可做的。一旦调用了`reinit()`，就可以通过 FEEvaluation::get_mapping_data_index_offset() 访问这种压缩。

// 一旦我们通过FEEvaluation访问函数访问了反雅各布系数（我们取第4个正交点的，因为它们在笛卡尔单元的所有正交点上都应该是一样的），我们检查它是对角线的，然后提取原始雅各布系数的行列式，即反雅各布系数的行列式，并根据质量矩阵的一维拉普拉斯乘以 $d-1$ 份，设置权重为 $\text{det}(J) / h_d^2$  。

    cell_matrices.clear(); 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data); 
    unsigned int old_mapping_data_index = numbers::invalid_unsigned_int; 
    for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 

        if (phi.get_mapping_data_index_offset() == old_mapping_data_index) 
          continue; 

        Tensor<2, dim, VectorizedArray<number>> inverse_jacobian = 
          phi.inverse_jacobian(0); 

        for (unsigned int d = 0; d < dim; ++d) 
          for (unsigned int e = 0; e < dim; ++e) 
            if (d != e) 
              for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v) 
                AssertThrow(inverse_jacobian[d][e][v] == 0., 
                            ExcNotImplemented()); 

        VectorizedArray<number> jacobian_determinant = inverse_jacobian[0][0]; 
        for (unsigned int e = 1; e < dim; ++e) 
          jacobian_determinant *= inverse_jacobian[e][e]; 
        jacobian_determinant = 1. / jacobian_determinant; 

        for (unsigned int d = 0; d < dim; ++d) 
          { 
            const VectorizedArray<number> scaling_factor = 
              inverse_jacobian[d][d] * inverse_jacobian[d][d] * 
              jacobian_determinant; 

// 一旦我们知道了拉普拉斯矩阵的比例系数，我们就将这个权重应用于未被缩放的DG拉普拉斯矩阵，并将数组发送到TensorProductMatrixSymmetricSum类，以计算介绍中提到的广义特征值问题。

            for (unsigned int i = 0; i < N; ++i) 
              for (unsigned int j = 0; j < N; ++j) 
                laplace_matrices[d](i, j) = 
                  scaling_factor * laplace_unscaled(i, j); 
          } 
        if (cell_matrices.size() <= phi.get_mapping_data_index_offset()) 
          cell_matrices.resize(phi.get_mapping_data_index_offset() + 1); 
        cell_matrices[phi.get_mapping_data_index_offset()].reinit( 
          mass_matrices, laplace_matrices); 
      } 
  } 

// 在DG背景下，用于近似块状Jacobi预处理的vmult函数非常简单。我们只需要读取当前单元格批次的值，对张量积矩阵阵列中的给定条目进行逆运算，并将结果写回。在这个循环中，我们覆盖了`dst`中的内容，而不是首先将条目设置为零。这对于DG方法来说是合法的，因为每个单元都有独立的自由度。此外，我们手动写出所有单元批的循环，而不是通过 MatrixFree::cell_loop(). 我们这样做是因为我们知道我们在这里不需要通过MPI网络进行数据交换，因为所有的计算都是在每个处理器上的本地单元上完成的。

  template <int dim, int fe_degree, typename number> 
  void PreconditionBlockJacobi<dim, fe_degree, number>::vmult( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    adjust_ghost_range_if_necessary(*data, dst); 
    adjust_ghost_range_if_necessary(*data, src); 

    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data); 
    for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 
        phi.read_dof_values(src); 
        cell_matrices[phi.get_mapping_data_index_offset()].apply_inverse( 
          ArrayView<VectorizedArray<number>>(phi.begin_dof_values(), 
                                             phi.dofs_per_cell), 
          ArrayView<const VectorizedArray<number>>(phi.begin_dof_values(), 
                                                   phi.dofs_per_cell)); 
        phi.set_dof_values(dst); 
      } 
  } 

// LaplaceProblem类的定义与  step-37  非常相似。一个区别是我们将元素度作为模板参数添加到类中，这将允许我们通过在`main()`函数中创建不同的实例，更容易在同一个程序中包含多个度。第二个区别是选择了FE_DGQHermite这个元素，它是专门用于这种方程的。

  template <int dim, int fe_degree> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void compute_rhs(); 
    void solve(); 
    void analyze_results() const; 

#ifdef DEAL_II_WITH_P4EST 
    parallel::distributed::Triangulation<dim> triangulation; 
#else 
    Triangulation<dim> triangulation; 
#endif 

    FE_DGQHermite<dim> fe; 
    DoFHandler<dim>    dof_handler; 

    MappingQ1<dim> mapping; 

    using SystemMatrixType = LaplaceOperator<dim, fe_degree, double>; 
    SystemMatrixType system_matrix; 

    using LevelMatrixType = LaplaceOperator<dim, fe_degree, float>; 
    MGLevelObject<LevelMatrixType> mg_matrices; 

    LinearAlgebra::distributed::Vector<double> solution; 
    LinearAlgebra::distributed::Vector<double> system_rhs; 

    double             setup_time; 
    ConditionalOStream pcout; 
    ConditionalOStream time_details; 
  }; 

  template <int dim, int fe_degree> 
  LaplaceProblem<dim, fe_degree>::LaplaceProblem() 
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
    fe(fe_degree) 
    , dof_handler(triangulation) 
    , setup_time(0.) 
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
    , time_details(std::cout, 
                   false && 
                     Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
  {} 

// 设置函数在两个方面与  step-37  不同。首先是我们不需要为不连续的Ansatz空间插值任何约束，而只是将一个假的AffineConstraints对象传入 Matrixfree::reinit().  第二个变化是因为我们需要告诉MatrixFree同时初始化面的数据结构。我们通过为内部面和边界面分别设置更新标志来做到这一点。在边界面，我们需要函数值、它们的梯度、JxW值（用于积分）、法向量和正交点（用于边界条件的评估），而对于内部面，我们只需要形状函数值、梯度、JxW值和法向量。只要`mapping_update_flags_inner_faces`或`mapping_update_flags_boundary_faces`中的一个与UpdateFlags的默认值`update_default`不同，MatrixFree中的面数据结构总是被建立的。

  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::setup_system() 
  { 
    Timer time; 
    setup_time = 0; 

    system_matrix.clear(); 
    mg_matrices.clear_elements(); 

    dof_handler.distribute_dofs(fe); 
    dof_handler.distribute_mg_dofs(); 

    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
          << std::endl; 

    setup_time += time.wall_time(); 
    time_details << "Distribute DoFs               " << time.wall_time() << " s" 
                 << std::endl; 
    time.restart(); 

    AffineConstraints<double> dummy; 
    dummy.close(); 

    { 
      typename MatrixFree<dim, double>::AdditionalData additional_data; 
      additional_data.tasks_parallel_scheme = 
        MatrixFree<dim, double>::AdditionalData::none; 
      additional_data.mapping_update_flags = 
        (update_gradients | update_JxW_values | update_quadrature_points); 
      additional_data.mapping_update_flags_inner_faces = 
        (update_gradients | update_JxW_values | update_normal_vectors); 
      additional_data.mapping_update_flags_boundary_faces = 
        (update_gradients | update_JxW_values | update_normal_vectors | 
         update_quadrature_points); 
      const auto system_mf_storage = 
        std::make_shared<MatrixFree<dim, double>>(); 
      system_mf_storage->reinit( 
        mapping, dof_handler, dummy, QGauss<1>(fe.degree + 1), additional_data); 
      system_matrix.initialize(system_mf_storage); 
    } 

    system_matrix.initialize_dof_vector(solution); 
    system_matrix.initialize_dof_vector(system_rhs); 

    setup_time += time.wall_time(); 
    time_details << "Setup matrix-free system      " << time.wall_time() << " s" 
                 << std::endl; 
    time.restart(); 

    const unsigned int nlevels = triangulation.n_global_levels(); 
    mg_matrices.resize(0, nlevels - 1); 

    for (unsigned int level = 0; level < nlevels; ++level) 
      { 
        typename MatrixFree<dim, float>::AdditionalData additional_data; 
        additional_data.tasks_parallel_scheme = 
          MatrixFree<dim, float>::AdditionalData::none; 
        additional_data.mapping_update_flags = 
          (update_gradients | update_JxW_values); 
        additional_data.mapping_update_flags_inner_faces = 
          (update_gradients | update_JxW_values); 
        additional_data.mapping_update_flags_boundary_faces = 
          (update_gradients | update_JxW_values); 
        additional_data.mg_level = level; 
        const auto mg_mf_storage_level = 
          std::make_shared<MatrixFree<dim, float>>(); 
        mg_mf_storage_level->reinit(mapping, 
                                    dof_handler, 
                                    dummy, 
                                    QGauss<1>(fe.degree + 1), 
                                    additional_data); 

        mg_matrices[level].initialize(mg_mf_storage_level); 
      } 
    setup_time += time.wall_time(); 
    time_details << "Setup matrix-free levels      " << time.wall_time() << " s" 
                 << std::endl; 
  } 

// 右手边的计算比  step-37  中的计算要复杂一些。现在的单元项包括分析解的负拉普拉斯，`RightHandSide'，为此我们需要首先将VectorizedArray字段的Point，即一批点，通过分别评估VectorizedArray中的所有通道，拆成一个点。请记住，通道的数量取决于硬件；对于不提供矢量化的系统（或deal.II没有本征），它可能是1，但在最近的Intel架构的AVX-512上也可能是8或16。

  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::compute_rhs() 
  { 
    Timer time; 
    system_rhs                          = 0; 
    const MatrixFree<dim, double> &data = *system_matrix.get_matrix_free(); 
    FEEvaluation<dim, fe_degree>   phi(data); 
    RightHandSide<dim>             rhs_func; 
    Solution<dim>                  exact_solution; 
    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            VectorizedArray<double> rhs_val = VectorizedArray<double>(); 
            Point<dim, VectorizedArray<double>> point_batch = 
              phi.quadrature_point(q); 
            for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v) 
              { 
                Point<dim> single_point; 
                for (unsigned int d = 0; d < dim; ++d) 
                  single_point[d] = point_batch[d][v]; 
                rhs_val[v] = rhs_func.value(single_point); 
              } 
            phi.submit_value(rhs_val, q); 
          } 
        phi.integrate_scatter(EvaluationFlags::values, system_rhs); 
      } 

// 其次，我们还需要应用Dirichlet和Neumann边界条件。一旦Dirichlet边界上的外部求解值 $u^+ = -u^- + 2 g_\text{D}$ 和 $\mathbf{n}^-\cdot \nabla u^+ = \mathbf{n}^-\cdot \nabla u^-$ 以及Neumann边界上的 $u^+=u^-$ 和 $\mathbf{n}^-\cdot \nabla u^+ = -\mathbf{n}^-\cdot \nabla u^- + 2 g_\text{N}$ 被插入并以边界函数 $g_\text{D}$ 和 $g_\text{N}$ 展开，这个函数就是到函数 `LaplaceOperator::apply_boundary()` 所缺的部分。需要记住的一点是，我们把边界条件移到右手边，所以符号与我们强加在解的部分相反。

// 我们可以通过 MatrixFree::loop 部分发出单元格和边界部分，但我们选择手动写出所有面的完整循环，以了解MatrixFree中面的索引布局是如何设置的：内部面和边界面都共享索引范围，所有批次的内部面的数字都比批次的边界单元格低。两种变体的单一索引使我们可以很容易地在两种情况下使用相同的数据结构FEFaceEvaluation，它附着在同一个数据域上，只是位置不同。内层面的批次（其中一个批次是由于将几个面合并成一个面进行矢量化）的数量由 MatrixFree::n_inner_face_batches(), 给出，而边界面的批次数量由 MatrixFree::n_boundary_face_batches(). 给出。
    FEFaceEvaluation<dim, fe_degree> phi_face(data, true); 
    for (unsigned int face = data.n_inner_face_batches(); 
         face < data.n_inner_face_batches() + data.n_boundary_face_batches(); 
         ++face) 
      { 
        phi_face.reinit(face); 

        const VectorizedArray<double> inverse_length_normal_to_face = 
          std::abs((phi_face.get_normal_vector(0) * 
                    phi_face.inverse_jacobian(0))[dim - 1]); 
        const VectorizedArray<double> sigma = 
          inverse_length_normal_to_face * system_matrix.get_penalty_factor(); 

        for (unsigned int q = 0; q < phi_face.n_q_points; ++q) 
          { 
            VectorizedArray<double> test_value = VectorizedArray<double>(), 
                                    test_normal_derivative = 
                                      VectorizedArray<double>(); 
            Point<dim, VectorizedArray<double>> point_batch = 
              phi_face.quadrature_point(q); 

            for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v) 
              { 
                Point<dim> single_point; 
                for (unsigned int d = 0; d < dim; ++d) 
                  single_point[d] = point_batch[d][v]; 

// MatrixFree类让我们查询当前面批的边界_id。请记住，MatrixFree为矢量化设置了批次，使一个批次中的所有面孔都有相同的属性，其中包括它们的`边界_id'。因此，我们可以在这里为当前面的索引`face`查询该id，并在Dirichlet情况下（我们在函数值上添加一些东西）或Neumann情况下（我们在法线导数上添加一些东西）施加。

                if (data.get_boundary_id(face) == 0) 
                  test_value[v] = 2.0 * exact_solution.value(single_point); 
                else 
                  { 
                    Tensor<1, dim> normal; 
                    for (unsigned int d = 0; d < dim; ++d) 
                      normal[d] = phi_face.get_normal_vector(q)[d][v]; 
                    test_normal_derivative[v] = 
                      -normal * exact_solution.gradient(single_point); 
                  } 
              } 
            phi_face.submit_value(test_value * sigma - test_normal_derivative, 
                                  q); 
            phi_face.submit_normal_derivative(-0.5 * test_value, q); 
          } 
        phi_face.integrate_scatter(EvaluationFlags::values | 
                                     EvaluationFlags::gradients, 
                                   system_rhs); 
      } 

// 由于我们手动运行了单元格的循环，而不是使用 MatrixFree::loop(), ，我们不能忘记与MPI进行数据交换。

// 或者说，对于DG元素来说，我们不需要这样做，因为每个单元都有自己的自由度，单元和边界积分只对本地拥有的单元进行评估。与相邻子域的耦合只通过内表面积分来实现，我们在这里没有做这个。也就是说，在这里调用这个函数并没有什么坏处，所以我们这样做是为了提醒大家在 MatrixFree::loop(). 里面发生了什么。
    system_rhs.compress(VectorOperation::add); 
    setup_time += time.wall_time(); 
    time_details << "Compute right hand side       " << time.wall_time() 
                 << " s\n"; 
  } 

// `solve()`函数几乎逐字复制自  step-37  。我们设置了相同的多网格成分，即水平转移、平滑器和粗略的网格求解器。唯一不同的是，我们没有使用拉普拉斯的对角线作为用于平滑的切比雪夫迭代的预处理，而是使用我们新解决的类`%PreconditionBlockJacobi`。不过，机制是一样的。

  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::solve() 
  { 
    Timer                            time; 
    MGTransferMatrixFree<dim, float> mg_transfer; 
    mg_transfer.build(dof_handler); 
    setup_time += time.wall_time(); 
    time_details << "MG build transfer time        " << time.wall_time() 
                 << " s\n"; 
    time.restart(); 

    using SmootherType = 
      PreconditionChebyshev<LevelMatrixType, 
                            LinearAlgebra::distributed::Vector<float>, 
                            PreconditionBlockJacobi<dim, fe_degree, float>>; 
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
            smoother_data[level].degree              = 3; 
            smoother_data[level].eig_cg_n_iterations = 10; 
          } 
        else 
          { 
            smoother_data[0].smoothing_range = 2e-2; 
            smoother_data[0].degree          = numbers::invalid_unsigned_int; 
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m(); 
          } 
        smoother_data[level].preconditioner = 
          std::make_shared<PreconditionBlockJacobi<dim, fe_degree, float>>(); 
        smoother_data[level].preconditioner->initialize(mg_matrices[level]); 
      } 
    mg_smoother.initialize(mg_matrices, smoother_data); 

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>> 
      mg_coarse; 
    mg_coarse.initialize(mg_smoother); 

    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix( 
      mg_matrices); 

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg( 
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother); 

    PreconditionMG<dim, 
                   LinearAlgebra::distributed::Vector<float>, 
                   MGTransferMatrixFree<dim, float>> 
      preconditioner(dof_handler, mg, mg_transfer); 

    SolverControl solver_control(10000, 1e-12 * system_rhs.l2_norm()); 
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control); 
    setup_time += time.wall_time(); 
    time_details << "MG build smoother time        " << time.wall_time() 
                 << "s\n"; 
    pcout << "Total setup time              " << setup_time << " s\n"; 

    time.reset(); 
    time.start(); 
    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    pcout << "Time solve (" << solver_control.last_step() << " iterations)    " 
          << time.wall_time() << " s" << std::endl; 
  } 

// 由于我们已经用分析法解决了一个问题，我们想通过计算数值结果与分析法的L2误差来验证我们实现的正确性。

  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::analyze_results() const 
  { 
    Vector<float> error_per_cell(triangulation.n_active_cells()); 
    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      solution, 
                                      Solution<dim>(), 
                                      error_per_cell, 
                                      QGauss<dim>(fe.degree + 2), 
                                      VectorTools::L2_norm); 
    pcout << "Verification via L2 error:    " 
          << std::sqrt( 
               Utilities::MPI::sum(error_per_cell.norm_sqr(), MPI_COMM_WORLD)) 
          << std::endl; 
  } 

// `run()`函数设置了初始网格，然后以常规方式运行多网格程序。作为一个域，我们选择一个矩形，在 $x$ -方向上有周期性的边界条件，在 $y$ 方向上的正面（即索引号为2的面，边界id等于0）有迪里希特条件，在背面以及 $z$ 方向上的两个面为三维情况（边界id等于1）有纽曼条件。与 $y$ 和 $z$ 方向相比， $x$ 方向的域的范围有些不同（考虑到 "解决方案 "的定义，我们希望在这里实现周期性的解决方案）。

  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::run() 
  { 
    const unsigned int n_ranks = 
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); 
    pcout << "Running with " << n_ranks << " MPI process" 
          << (n_ranks > 1 ? "es" : "") << ", element " << fe.get_name() 
          << std::endl 
          << std::endl; 
    for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle) 
      { 
        pcout << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            Point<dim> upper_right; 
            upper_right[0] = 2.5; 
            for (unsigned int d = 1; d < dim; ++d) 
              upper_right[d] = 2.8; 
            GridGenerator::hyper_rectangle(triangulation, 
                                           Point<dim>(), 
                                           upper_right); 
            triangulation.begin_active()->face(0)->set_boundary_id(10); 
            triangulation.begin_active()->face(1)->set_boundary_id(11); 
            triangulation.begin_active()->face(2)->set_boundary_id(0); 
            for (unsigned int f = 3; 
                 f < triangulation.begin_active()->n_faces(); 
                 ++f) 
              triangulation.begin_active()->face(f)->set_boundary_id(1); 

            std::vector<GridTools::PeriodicFacePair< 
              typename Triangulation<dim>::cell_iterator>> 
              periodic_faces; 
            GridTools::collect_periodic_faces( 
              triangulation, 10, 11, 0, periodic_faces); 
            triangulation.add_periodicity(periodic_faces); 

            triangulation.refine_global(6 - 2 * dim); 
          } 
        triangulation.refine_global(1); 
        setup_system(); 
        compute_rhs(); 
        solve(); 
        analyze_results(); 
        pcout << std::endl; 
      }; 
  } 
} // namespace Step59 

// `main()`函数中没有任何意外。我们通过`MPI_Init()`类调用`MPI_InitFinalize`，传入文件顶部设置的关于维度和度的两个参数，然后运行拉普拉斯问题。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace Step59; 

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1); 

      LaplaceProblem<dimension, degree_finite_element> laplace_problem; 
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


