//include/deal.II-translator/matrix_free/operators_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_operators_h
#define dealii_matrix_free_operators_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>


DEAL_II_NAMESPACE_OPEN


namespace MatrixFreeOperators
{
  namespace BlockHelper
  {
    // workaroud for unifying non-block vector and block vector implementations
    // a non-block vector has one block and the only subblock is the vector
    // itself
    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value,
                            unsigned int>::type
    n_blocks(const VectorType &vector)
    {
      return vector.n_blocks();
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value,
                            unsigned int>::type
    n_blocks(const VectorType &)
    {
      return 1;
    }

    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value,
                            typename VectorType::BlockType &>::type
    subblock(VectorType &vector, unsigned int block_no)
    {
      return vector.block(block_no);
    }

    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value,
                            const typename VectorType::BlockType &>::type
    subblock(const VectorType &vector, unsigned int block_no)
    {
      AssertIndexRange(block_no, vector.n_blocks());
      return vector.block(block_no);
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value,
                            VectorType &>::type
    subblock(VectorType &vector, unsigned int)
    {
      return vector;
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value,
                            const VectorType &>::type
    subblock(const VectorType &vector, unsigned int)
    {
      return vector;
    }

    template <typename VectorType>
    typename std::enable_if<IsBlockVector<VectorType>::value, void>::type
    collect_sizes(VectorType &vector)
    {
      vector.collect_sizes();
    }

    template <typename VectorType>
    typename std::enable_if<!IsBlockVector<VectorType>::value, void>::type
    collect_sizes(const VectorType &)
    {}
  } // namespace BlockHelper

  /**
   * 无矩阵运算符的抽象基类，可以用于最细的网格或几何多网格中的某个层次。
   * 派生类必须实现apply_add()方法以及compute_diagonal()来初始化受保护成员inverse_diagonal_entries和/或diagonal_entries。在非对称运算符的情况下，应该额外实现Tapply_add()。
   * 目前，唯一支持的向量是 LinearAlgebra::distributed::Vector 和
   * LinearAlgebra::distributed::BlockVector.  <h4>Selective use of blocks in
   * MatrixFree</h4>
   * MatrixFree允许使用几个DoFHandler/AffineConstraints组合，方法是在
   * MatrixFree::reinit
   * 函数中传递一个带有指向各自对象的指针的 std::vector
   * 。该类支持在底层MatrixFree对象中只选择一些块，通过可选的整数列表来指定所选的块。
   * 只在选定的块上构造无矩阵运算符的一个应用是 step-32
   * 教程程序的设置。这个问题有三个<i>blocks</i>，一个用于速度，一个用于压力，一个用于温度。用于时间演化的时滞方案将温度方程从速度和压力的斯托克斯系统中分割出来。然而，还有一些交叉项，如进入温度平流-扩散方程的速度或进入速度右侧的温度。为了确保MatrixFree对不同的块使用相同的整数索引，人们需要将所有三个块放入同一个MatrixFree对象。然而，当求解一个线性系统时，所涉及的运算符要么在斯托克斯求解器中针对前两个，要么在温度求解器中针对最后一个。在前一种情况下，一个有两个分量的BlockVector将被选择，在MatrixFree中选择块{0,
   * 1}，而在后一种情况下，将使用一个选择块{2}的非块向量。
   * 选择的第二个应用是在具有牛顿式迭代的问题或具有不均匀边界条件的问题中。在这种情况下，我们必须处理两组不同的约束。一套约束条件适用于解向量，其中可能包括悬挂节点约束或周期性约束，但没有对不均匀的Dirichlet边界的约束。在非线性迭代之前，边界值被设置为向量中的预期值，代表初始猜测。在牛顿方法的每次迭代中，一个受零Dirichlet边界条件约束的线性系统被解决，然后被加入到初始猜测中。这个设置可以通过使用一个指向同一DoFHandler对象的两个指针的向量和一个指向两个AffineConstraints对象的两个指针的向量来实现。如果第一个AffineConstraints对象是包括零Dirichlet约束的对象，我们将给initialize()函数一个
   * std::vector<unsigned  int>(1,
   * 0)，也就是说，一个长度为1的向量，正好选择索引为0的第一个AffineConstraints对象。
   * 对PDEs系统来说，MatrixFree的不同块与方程的不同物理成分相关，仅仅为了边界条件的目的而添加另一个具有不同AffineConstraints参数的块可能导致繁琐的索引处理。相反，可以用不同的AffineConstraints对象设置第二个MatrixFree实例，但对块的解释相同，并使用它来插值不均匀的边界条件（也可参见
   * step-37 教程程序中结果部分的讨论）。
   * @code
   * matrix_free_inhomogeneous.reinit(dof_handler, constraints_no_dirichlet,
   *                                quadrature, additional_data);
   * operator_inhomogeneous.initialize(matrix_free_inhomogeneous,
   *                                 selected_blocks);
   * LinearAlgebra::distributed::Vector<double> inhomogeneity;
   * matrix_free_inhomogeneous.initialize_dof_vector(inhomogeneity);
   * constraints_with_dirichlet.distribute(inhomogeneity);
   * operator_inhomogeneous.vmult(system_rhs, inhomogeneity);
   * system_rhs=
   *
   * -1.;
   * // proceed with other terms from right hand side...
   * @endcode
   *
   *
   */
  template <int dim,
            typename VectorType = LinearAlgebra::distributed::Vector<double>,
            typename VectorizedArrayType =
              VectorizedArray<typename VectorType::value_type>>
  class Base : public Subscriptor
  {
  public:
    /**
     * 数字别名。
     *
     */
    using value_type = typename VectorType::value_type;

    /**
     * 先决条件类所需的size_type。
     *
     */
    using size_type = typename VectorType::size_type;

    /**
     * 默认构造函数。
     *
     */
    Base();

    /**
     * 虚拟解构器。
     *
     */
    virtual ~Base() override = default;

    /**
     * 释放所有内存并返回到与调用默认构造函数后相同的状态。
     *
     */
    virtual void
    clear();

    /**
     * 初始化操作者在精细的规模上。
     * 可选的选择向量允许从底层的MatrixFree对象中只选择一些组件，例如，只选择一个。矢量中的
     * @p selected_row_blocks[i]
     * 项选择了DoFHandler和AffineConstraints对象，该对象被作为
     * @p selected_row_blocks[i]-th 参数给到 MatrixFree::reinit()
     * 调用。
     * 行和列的不同参数也使得选择非对角线块或矩形块成为可能。如果行向量为空，则选择所有组件，否则其大小必须小于或等于
     * MatrixFree::n_components()
     * ，并且所有指数需要是唯一的，并且在0和
     * MatrixFree::n_components(). 的范围内
     * 如果列选择向量为空，其取值与行选择相同，定义一个对角块。
     *
     */
    void
    initialize(std::shared_ptr<
                 const MatrixFree<dim, value_type, VectorizedArrayType>> data,
               const std::vector<unsigned int> &selected_row_blocks =
                 std::vector<unsigned int>(),
               const std::vector<unsigned int> &selected_column_blocks =
                 std::vector<unsigned int>());

    /**
     * 在一个层次 @p level
     * 上对单个有限元素进行初始化操作。
     * 可选的选择向量允许从底层的MatrixFree对象中只选择一些组件，例如只选择一个。矢量中的
     * @p selected_row_blocks[i]
     * 项选择DoFHandler和AffineConstraints对象，该对象被作为 @p
     * selected_row_blocks[i]-th 参数给到 MatrixFree::reinit() 调用。
     * 由于多网格运算符总是与反转矩阵相关，因此代表对角线块，相对于非水平初始化函数，行和列的向量是相同的。如果是空的，所有的组件都被选中。
     *
     */
    void
    initialize(std::shared_ptr<
                 const MatrixFree<dim, value_type, VectorizedArrayType>> data,
               const MGConstrainedDoFs &        mg_constrained_dofs,
               const unsigned int               level,
               const std::vector<unsigned int> &selected_row_blocks =
                 std::vector<unsigned int>());

    /**
     * 对多个FiniteElement对象进行水平 @p level 初始化操作。
     * 可选的选择向量允许只从底层的MatrixFree对象中选择一些组件，例如，只选择一个。矢量中的
     * @p selected_row_blocks[i]
     * 项选择DoFHandler和AffineConstraints对象，该对象被作为 @p
     * selected_row_blocks[i]-th 参数给到 MatrixFree::reinit() 调用。
     * 由于多网格运算符总是与反转矩阵相关，因此代表对角线块，相对于非水平初始化函数，行和列的向量是相同的。如果是空的，所有的组件都被选中。
     *
     */
    void
    initialize(std::shared_ptr<
                 const MatrixFree<dim, value_type, VectorizedArrayType>> data_,
               const std::vector<MGConstrainedDoFs> &mg_constrained_dofs,
               const unsigned int                    level,
               const std::vector<unsigned int> &     selected_row_blocks =
                 std::vector<unsigned int>());

    /**
     * 返回共域（或范围）空间的维度。
     *
     */
    size_type
    m() const;

    /**
     * 返回域空间的维数。
     *
     */
    size_type
    n() const;

    /**
     * 接口的vmult运算符。
     *
     */
    void
    vmult_interface_down(VectorType &dst, const VectorType &src) const;

    /**
     * 用于接口的vmult运算符。
     *
     */
    void
    vmult_interface_up(VectorType &dst, const VectorType &src) const;

    /**
     * 矩阵-向量乘法。
     *
     */
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * 转置矩阵-向量乘法。
     *
     */
    void
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * 加法 矩阵-向量乘法。
     *
     */
    void
    vmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * 加法转置矩阵-向量乘法。
     *
     */
    void
    Tvmult_add(VectorType &dst, const VectorType &src) const;

    /**
     * 返回矩阵条目(row,col)的值。在无矩阵情况下，当对角线被初始化时，这个函数只对row==col有效。
     *
     */
    value_type
    el(const unsigned int row, const unsigned int col) const;

    /**
     * 确定此对象的内存消耗（以字节为单位）的估计值。
     *
     */
    virtual std::size_t
    memory_consumption() const;

    /**
     * MatrixFree对象的initialize_dof_vector()的一个封装器。
     *
     */
    void
    initialize_dof_vector(VectorType &vec) const;

    /**
     * 计算此运算符的对角线。
     * 派生类需要实现这个函数并相应地调整和填充保护成员inverse_diagonal_entries和/或diagonal_entries的大小。
     *
     */
    virtual void
    compute_diagonal() = 0;

    /**
     * 获得对用该操作符存储的MatrixFree对象的读取权限。
     *
     */
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
    get_matrix_free() const;

    /**
     * 获得对该运算符的反对角线的读取权限。
     *
     */
    const std::shared_ptr<DiagonalMatrix<VectorType>> &
    get_matrix_diagonal_inverse() const;

    /**
     * 获取对这个运算符的对角线的读取权限。
     *
     */
    const std::shared_ptr<DiagonalMatrix<VectorType>> &
    get_matrix_diagonal() const;


    /**
     * 应用雅可比预处理程序，该程序将<tt>src</tt>向量的每个元素乘以各自对角线元素的逆值，并将结果乘以松弛因子<tt>omega</tt>。
     *
     */
    void
    precondition_Jacobi(VectorType &      dst,
                        const VectorType &src,
                        const value_type  omega) const;

  protected:
    /**
     * 在调用mult_add()内的apply_add()或Tapply_add()之前，执行与约束有关的必要操作。
     *
     */
    void
    preprocess_constraints(VectorType &dst, const VectorType &src) const;

    /**
     * 在调用mult_add()内的apply_add()或Tapply_add()后，执行与约束有关的必要操作。
     *
     */
    void
    postprocess_constraints(VectorType &dst, const VectorType &src) const;

    /**
     * 将 @p dst
     * 的约束条目（包括来自悬挂节点和边缘约束）设为1。
     *
     */
    void
    set_constrained_entries_to_one(VectorType &dst) const;

    /**
     * 对 @p src 应用运算符，并在 @p dst. 中添加结果。
     *
     */
    virtual void
    apply_add(VectorType &dst, const VectorType &src) const = 0;

    /**
     * 对 @p src 应用转置运算符，并在 @p dst. 中添加结果
     * 默认实现是调用apply_add()。
     *
     */
    virtual void
    Tapply_add(VectorType &dst, const VectorType &src) const;

    /**
     * MatrixFree对象将与该操作符一起使用。
     *
     */
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
      data;

    /**
     * 一个指向对角线矩阵的共享指针，将对角线元素作为一个向量存储。
     *
     */
    std::shared_ptr<DiagonalMatrix<VectorType>> diagonal_entries;

    /**
     * 一个指向对角线矩阵的共享指针，该矩阵将对角线元素的倒数作为一个向量存储。
     *
     */
    std::shared_ptr<DiagonalMatrix<VectorType>> inverse_diagonal_entries;

    /**
     * 一个向量，它定义了MatrixFree的子组件对矩阵表示的行的选择。
     *
     */
    std::vector<unsigned int> selected_rows;

    /**
     * 一个向量，它定义了矩阵表示的列的MatrixFree的子组件的选择。
     *
     */
    std::vector<unsigned int> selected_columns;

  private:
    /**
     * 如果操作者在GMG背景下使用，边缘上的DoF的索引。
     *
     */
    std::vector<std::vector<unsigned int>> edge_constrained_indices;

    /**
     * 辅助向量。
     *
     */
    mutable std::vector<std::vector<std::pair<value_type, value_type>>>
      edge_constrained_values;

    /**
     * 一个标志，决定该运算符在GMG背景下是否有接口矩阵。
     *
     */
    bool have_interface_matrices;

    /**
     * 实现vmult_add (  @p transpose  = false) 和Tvmult_add (  @p transpose
     * = true) 的%函数。
     *
     */
    void
    mult_add(VectorType &      dst,
             const VectorType &src,
             const bool        transpose) const;

    /**
     * 根据底层MatrixFree类的存储要求，调整向量的ghost范围。这在mult_add()以及vmult_interface_up()和vmult_interface_down()方法中使用，以确保单元格循环能够以正确的本地索引访问ghost索引。
     *
     */
    void
    adjust_ghost_range_if_necessary(const VectorType &vec,
                                    const bool        is_row) const;
  };



  /**
   * 辅助类，提供自适应几何多网格中需要的接口vmult/Tvmult方法。  @p OperatorType 类应派生于 MatrixFreeOperators::Base 类。    deal.II中的自适应多网格实现了一种叫做局部平滑的方法。这意味着最细级别的平滑只覆盖固定（最细）网格级别所定义的网格的局部部分，而忽略了计算域中终端单元比该级别更粗的部分。随着该方法向更粗的级别发展，越来越多的全局网格将被覆盖。在某个更粗的层次上，整个网格将被覆盖。由于多网格方法中的所有层次矩阵都覆盖了网格中的单一层次，所以在层次矩阵上不会出现悬空节点。在多网格层之间的界面上，在平滑的同时设置同质Dirichlet边界条件。然而，当残差被转移到下一个更粗的层次时，需要考虑到多网格界面的耦合。  这是由所谓的界面（或边缘）矩阵来完成的，它计算了被具有同质Dirichlet条件的层次矩阵所遗漏的残差部分。我们参考 @ref mg_paper "Janssen和Kanschat的多网格论文 "
   * 以了解更多细节。
   * 对于这些接口矩阵的实现，大部分基础设施已经到位，由
   * MatrixFreeOperators::Base
   * 通过两个乘法例程vmult_interface_down()和vmult_interface_up()提供。MGInterfaceOperator所做的唯一事情是包装这些操作，并使它们可以通过
   * @p vmult() 和 @p Tvmult
   * 接口访问，正如多网格例程（最初是为矩阵编写的，因此期待这些名称）所期望的那样。
   * 请注意，vmult_interface_down在多网格V周期的限制阶段使用，而vmult_interface_up在延长阶段使用。
   *
   */
  template <typename OperatorType>
  class MGInterfaceOperator : public Subscriptor
  {
  public:
    /**
     * 编号别名。
     *
     */
    using value_type = typename OperatorType::value_type;

    /**
     * 尺寸类型。
     *
     */
    using size_type = typename OperatorType::size_type;

    /**
     * 默认构造函数。
     *
     */
    MGInterfaceOperator();

    /**
     * 清除指向OperatorType对象的指针。
     *
     */
    void
    clear();

    /**
     * 用一个操作符来初始化这个类  @p operator_in.
     *
     */
    void
    initialize(const OperatorType &operator_in);

    /**
     * vmult运算符，更多信息见类描述。
     *
     */
    template <typename VectorType>
    void
    vmult(VectorType &dst, const VectorType &src) const;

    /**
     * Tvmult运算符，更多信息见类的描述。
     *
     */
    template <typename VectorType>
    void
    Tvmult(VectorType &dst, const VectorType &src) const;

    /**
     * OperatorType对象的initialize_dof_vector()的封装器。
     *
     */
    template <typename VectorType>
    void
    initialize_dof_vector(VectorType &vec) const;


  private:
    /**
     * 指向操作员类的常量指针。
     *
     */
    SmartPointer<const OperatorType> mf_base_operator;
  };



  /**
   * 该类实现了质量矩阵的反作用于元素的操作，适用于评估对象的正交点与单元自由度一样多的特殊情况。它使用FEEvaluation的算法，为DGQ元素产生精确的质量矩阵。该算法使用正交点上的反一维形状矩阵的张量乘积，因此反操作与在每个单元上应用正向操作一样昂贵。当然，对于连续有限元来说，这种运算并不能产生质量运算的逆运算，因为元素之间的耦合不能被这种运算所考虑。
   * 方程可能包含可变系数，所以用户需要提供一个数组，用于局部系数的逆运算（该类提供了一个辅助方法'fill_inverse_JxW_values'来获得常数系数运算的逆运算）。
   *
   */
  template <int dim,
            int fe_degree,
            int n_components             = 1,
            typename Number              = double,
            typename VectorizedArrayType = VectorizedArray<Number>>
  class CellwiseInverseMassMatrix
  {
    static_assert(
      std::is_same<Number, typename VectorizedArrayType::value_type>::value,
      "Type of Number and of VectorizedArrayType do not match.");

  public:
    /**
     * 构造函数。从FEEval类中的ShapeInfo字段初始化形状信息。
     *
     */
    CellwiseInverseMassMatrix(
      const FEEvaluationBase<dim,
                             n_components,
                             Number,
                             false,
                             VectorizedArrayType> &fe_eval);

    /**
     * 在一个输入数组上应用反质量矩阵操作。假设传递的输入和输出数组的大小正确，即
     * FEEvaluation::dofs_per_cell
     * 长。本地系数的逆值（也包含逆JxW值）必须作为第一个参数传递。允许在系数中传递一个以上的分量。
     *
     */
    void
    apply(const AlignedVector<VectorizedArrayType> &inverse_coefficient,
          const unsigned int                        n_actual_components,
          const VectorizedArrayType *               in_array,
          VectorizedArrayType *                     out_array) const;

    /**
     * 在一个输入数组上应用反质量矩阵操作，使用由传递给该类构造函数的`fe_eval`参数提供的JxW值的逆值。注意，用户代码必须在底层评估器上调用
     * FEEvaluation::reinit() ，以使 FEEvaluationBase::JxW()
     * 方法返回正确单元的信息。假设输入和输出数组的指针在长度
     * FEEvaluation::dofs_per_cell,
     * 上是有效的，这个长度是这个函数处理的条目数。`in_array`和`out_array`参数可以指向相同的内存位置。
     *
     */
    void
    apply(const VectorizedArrayType *in_array,
          VectorizedArrayType *      out_array) const;

    /**
     * 这个操作执行从正交点中给出的数据到这个对象的实际基础的投影。这个投影也可以解释为从正交点的拉格朗日插值多项式到当前`fe_eval`对象的基础的变化。
     * 在一个数组上调用这个函数为
     * @code
     * inverse_mass.transform_from_q_points_to_basis(1, array,
     *                                             phi.begin_dof_values());
     * @endcode
     * 相当于
     * @code
     * for (unsigned int q=0; q<phi.n_q_points; ++q)
     * phi.submit_value(array[q], q);
     * phi.integrate(EvaluationFlags::values);
     * inverse_mass.apply(coefficients, 1, phi.begin_dof_values(),
     *                  phi.begin_dof_values());
     * @endcode
     * 只要 @p coefficients
     * 持有正交权重的逆值，没有额外的系数。这种设置突出了基本的投影，测试了一个右手边，并应用了一个反质量矩阵。这个函数既适用于例子中描述的标量情况，也适用于逐个组件排列的多个组件。
     * 与更繁琐的替代方法相比，给定的程序要快得多，因为它可以绕过
     * @p integrate()
     * 步骤和正交点转换的前半部分，将张量积调用的数量从3*dim*n_components减少到dim*n_components。
     *
     */
    void
    transform_from_q_points_to_basis(const unsigned int n_actual_components,
                                     const VectorizedArrayType *in_array,
                                     VectorizedArrayType *out_array) const;

    /**
     * 用JxW值的逆值填充给定数组，即系数为1的质量矩阵。非单位系数必须与该数组相乘（以反形式）。
     *
     */
    void
    fill_inverse_JxW_values(
      AlignedVector<VectorizedArrayType> &inverse_jxw) const;

  private:
    /**
     * 对FEEvaluation对象的引用，用于获取JxW_值。
     *
     */
    const FEEvaluationBase<dim,
                           n_components,
                           Number,
                           false,
                           VectorizedArrayType> &fe_eval;
  };



  /**
   * 这个类实现了质量矩阵的动作操作。
   * 请注意，这个类只支持Base操作的非阻塞矢量变体，因为在apply函数中只使用一个FEEvaluation对象。
   *
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d   = fe_degree + 1,
            int n_components    = 1,
            typename VectorType = LinearAlgebra::distributed::Vector<double>,
            typename VectorizedArrayType =
              VectorizedArray<typename VectorType::value_type>>
  class MassOperator : public Base<dim, VectorType, VectorizedArrayType>
  {
  public:
    /**
     * 数字别名。
     *
     */
    using value_type =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;

    /**
     * 先决条件类需要的size_type。
     *
     */
    using size_type =
      typename Base<dim, VectorType, VectorizedArrayType>::size_type;

    /**
     * 构造函数。
     *
     */
    MassOperator();

    /**
     * 对于预处理，我们在对角线条目处存储一个块状的质量矩阵。
     *
     */
    virtual void
    compute_diagonal() override;

  private:
    /**
     * 在一个输入向量上应用质量矩阵运算。假设通过的输入和输出向量使用initialize_dof_vector()进行了正确的初始化。
     *
     */
    virtual void
    apply_add(VectorType &dst, const VectorType &src) const override;

    /**
     * 对于这个运算符，只有一个单元格的贡献。
     *
     */
    void
    local_apply_cell(
      const MatrixFree<dim, value_type, VectorizedArrayType> &data,
      VectorType &                                            dst,
      const VectorType &                                      src,
      const std::pair<unsigned int, unsigned int> &           cell_range) const;
  };



  /**
   * 这个类实现了拉普拉斯矩阵的作用的操作，即 $ L_{ij} =
   * \int_\Omega c(\mathbf x) \mathbf \nabla N_i(\mathbf x) \cdot \mathbf
   * \nabla N_j(\mathbf x)\,d \mathbf x$ ，其中 $c(\mathbf x)$
   * 是标量异质性系数。
   * 请注意，这个类只支持Base操作的非阻塞矢量变体，因为在apply函数中只使用一个FEEvaluation对象。
   *
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d   = fe_degree + 1,
            int n_components    = 1,
            typename VectorType = LinearAlgebra::distributed::Vector<double>,
            typename VectorizedArrayType =
              VectorizedArray<typename VectorType::value_type>>
  class LaplaceOperator : public Base<dim, VectorType, VectorizedArrayType>
  {
  public:
    /**
     * 数字别名。
     *
     */
    using value_type =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;

    /**
     * 先决条件类需要的size_type。
     *
     */
    using size_type =
      typename Base<dim, VectorType, VectorizedArrayType>::size_type;

    /**
     * 构造函数。
     *
     */
    LaplaceOperator();

    /**
     * 对角线是通过计算每个元素的局部对角线矩阵并将其分配给全局对角线来近似计算的。这将导致在有悬空节点的元素上出现错误的结果，但仍是一个可接受的近似值，可用于预处理程序中。
     *
     */
    virtual void
    compute_diagonal() override;

    /**
     * 设置异质标量系数 @p scalar_coefficient
     * ，在正交点使用。表的行数需要与底层MatrixFree对象中的单元格批次一样多，
     * MatrixFree::n_cell_batches().
     * 批次的数量与无矩阵运算符不在单个单元格上工作有关，而是由于矢量化而一次对单元格的批次工作。表可以采取两种不同的列数。
     * 一种情况是选择它等于`dim`维度上的正交点总数，这是`n_q_points_1d`模板参数的`dim`次方。这里，`(*scalar_coefficient)(cell,q)`对应于单元批`cell`和正交点指数`q`上的系数值。第二个支持的变体是一个单列的表，在这种情况下，一个单元格的所有正交点都使用相同的可变系数值。
     * 这种表可以通过以下方式初始化
     * @code
     * std::shared_ptr<Table<2, VectorizedArray<double> > > coefficient;
     * coefficient = std::make_shared<Table<2, VectorizedArray<double> > >();
     * {
     * FEEvaluation<dim,fe_degree,n_q_points_1d,1,double> fe_eval(mf_data);
     * const unsigned int n_cells = mf_data.n_cell_batches();
     * const unsigned int n_q_points = fe_eval.n_q_points;
     * coefficient->reinit(n_cells, n_q_points);
     * for (unsigned int cell=0; cell<n_cells; ++cell)
     *   {
     *     fe_eval.reinit(cell);
     *     for (unsigned int q=0; q<n_q_points; ++q)
     *       (*coefficient)(cell,q) =
     *         function.value(fe_eval.quadrature_point(q));
     *   }
     * }
     * @endcode
     * 其中 <code>mf_data</code> 是一个MatrixFree对象，
     * <code>function</code>
     * 是一个函数，它提供了以下方法<code>VectorizedArray<double>
     * value(const Point<dim, VectorizedArray<double> > &p_vec)</code>。
     * 如果这个函数没有被调用，那么系数就被认为是统一的。
     * 这个函数的参数是一个指向这样一个表的共享指针。该类存储了这个表的共享指针，而不是一个深度拷贝，并使用它来形成拉普拉斯矩阵。因此，你可以更新该表并重新使用当前对象，以获得具有该更新系数的拉普拉斯矩阵的动作。另外，如果表的值只需要填充一次，原来的共享指针也可以在用户代码中超出范围，clear()命令或这个类的析构器将删除表。
     *
     */
    void
    set_coefficient(
      const std::shared_ptr<Table<2, VectorizedArrayType>> &scalar_coefficient);

    /**
     * 将所有数据结构重新设置为与新构建的对象相同的状态。
     *
     */
    virtual void
    clear() override;

    /**
     * 读/写对拉普拉斯运算器中使用的系数的访问。
     * 如果之前没有通过set_coefficient()函数设置系数，该函数将抛出一个错误。
     *
     */
    std::shared_ptr<Table<2, VectorizedArrayType>>
    get_coefficient();

  private:
    /**
     * 在一个输入矢量上应用拉普拉斯矩阵运算。假设传递的输入和输出向量已经用initialize_dof_vector()函数正确初始化。
     *
     */
    virtual void
    apply_add(VectorType &dst, const VectorType &src) const override;

    /**
     * 在一个单元格上应用拉普拉斯算子。
     *
     */
    void
    local_apply_cell(
      const MatrixFree<dim, value_type, VectorizedArrayType> &data,
      VectorType &                                            dst,
      const VectorType &                                      src,
      const std::pair<unsigned int, unsigned int> &           cell_range) const;

    /**
     * 在一个单元格上应用拉普拉斯算子的对角线部分。
     *
     */
    void
    local_diagonal_cell(
      const MatrixFree<dim, value_type, VectorizedArrayType> &data,
      VectorType &                                            dst,
      const VectorType &,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    /**
     * 在单元格上应用拉普拉斯算子 @p cell. 。
     *
     */
    void
    do_operation_on_cell(
      FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, value_type>
        &                phi,
      const unsigned int cell) const;

    /**
     * 用户提供的异质性系数。
     *
     */
    std::shared_ptr<Table<2, VectorizedArrayType>> scalar_coefficient;
  };



  // ------------------------------------ inline functions ---------------------

  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline CellwiseInverseMassMatrix<dim,
                                   fe_degree,
                                   n_components,
                                   Number,
                                   VectorizedArrayType>::
    CellwiseInverseMassMatrix(
      const FEEvaluationBase<dim,
                             n_components,
                             Number,
                             false,
                             VectorizedArrayType> &fe_eval)
    : fe_eval(fe_eval)
  {}



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<dim,
                            fe_degree,
                            n_components,
                            Number,
                            VectorizedArrayType>::
    fill_inverse_JxW_values(
      AlignedVector<VectorizedArrayType> &inverse_jxw) const
  {
    constexpr unsigned int dofs_per_component_on_cell =
      Utilities::pow(fe_degree + 1, dim);
    Assert(inverse_jxw.size() > 0 &&
             inverse_jxw.size() % dofs_per_component_on_cell == 0,
           ExcMessage(
             "Expected diagonal to be a multiple of scalar dof per cells"));

    // compute values for the first component
    for (unsigned int q = 0; q < dofs_per_component_on_cell; ++q)
      inverse_jxw[q] = 1. / fe_eval.JxW(q);
    // copy values to rest of vector
    for (unsigned int q = dofs_per_component_on_cell; q < inverse_jxw.size();)
      for (unsigned int i = 0; i < dofs_per_component_on_cell; ++i, ++q)
        inverse_jxw[q] = inverse_jxw[i];
  }



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<
    dim,
    fe_degree,
    n_components,
    Number,
    VectorizedArrayType>::apply(const VectorizedArrayType *in_array,
                                VectorizedArrayType *      out_array) const
  {
    internal::CellwiseInverseMassMatrixImplBasic<dim, VectorizedArrayType>::
      template run<fe_degree>(n_components, fe_eval, in_array, out_array);
  }



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<dim,
                            fe_degree,
                            n_components,
                            Number,
                            VectorizedArrayType>::
    apply(const AlignedVector<VectorizedArrayType> &inverse_coefficients,
          const unsigned int                        n_actual_components,
          const VectorizedArrayType *               in_array,
          VectorizedArrayType *                     out_array) const
  {
    internal::CellwiseInverseMassMatrixImplFlexible<dim, VectorizedArrayType>::
      template run<fe_degree>(
        n_actual_components,
        fe_eval.get_shape_info().data.front().inverse_shape_values_eo,
        inverse_coefficients,
        in_array,
        out_array);
  }



  template <int dim,
            int fe_degree,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  inline void
  CellwiseInverseMassMatrix<dim,
                            fe_degree,
                            n_components,
                            Number,
                            VectorizedArrayType>::
    transform_from_q_points_to_basis(const unsigned int n_actual_components,
                                     const VectorizedArrayType *in_array,
                                     VectorizedArrayType *      out_array) const
  {
    internal::CellwiseInverseMassMatrixImplTransformFromQPoints<
      dim,
      VectorizedArrayType>::template run<fe_degree>(n_actual_components,
                                                    fe_eval,
                                                    in_array,
                                                    out_array);
  }



  //----------------- Base operator -----------------------------
  template <int dim, typename VectorType, typename VectorizedArrayType>
  Base<dim, VectorType, VectorizedArrayType>::Base()
    : Subscriptor()
    , have_interface_matrices(false)
  {}



  template <int dim, typename VectorType, typename VectorizedArrayType>
  typename Base<dim, VectorType, VectorizedArrayType>::size_type
  Base<dim, VectorType, VectorizedArrayType>::m() const
  {
    Assert(data.get() != nullptr, ExcNotInitialized());
    typename Base<dim, VectorType, VectorizedArrayType>::size_type total_size =
      0;
    for (const unsigned int selected_row : selected_rows)
      total_size += data->get_vector_partitioner(selected_row)->size();
    return total_size;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  typename Base<dim, VectorType, VectorizedArrayType>::size_type
  Base<dim, VectorType, VectorizedArrayType>::n() const
  {
    Assert(data.get() != nullptr, ExcNotInitialized());
    typename Base<dim, VectorType, VectorizedArrayType>::size_type total_size =
      0;
    for (const unsigned int selected_column : selected_columns)
      total_size += data->get_vector_partitioner(selected_column)->size();
    return total_size;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::clear()
  {
    data.reset();
    inverse_diagonal_entries.reset();
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  typename Base<dim, VectorType, VectorizedArrayType>::value_type
  Base<dim, VectorType, VectorizedArrayType>::el(const unsigned int row,
                                                 const unsigned int col) const
  {
    (void)col;
    Assert(row == col, ExcNotImplemented());
    Assert(inverse_diagonal_entries.get() != nullptr &&
             inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    return 1.0 / (*inverse_diagonal_entries)(row, row);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize_dof_vector(
    VectorType &vec) const
  {
    Assert(data.get() != nullptr, ExcNotInitialized());
    AssertDimension(BlockHelper::n_blocks(vec), selected_rows.size());
    for (unsigned int i = 0; i < BlockHelper::n_blocks(vec); ++i)
      {
        const unsigned int index = selected_rows[i];
        if (!BlockHelper::subblock(vec, index)
               .partitioners_are_compatible(
                 *data->get_dof_info(index).vector_partitioner))
          data->initialize_dof_vector(BlockHelper::subblock(vec, index), index);

        Assert(BlockHelper::subblock(vec, index)
                 .partitioners_are_globally_compatible(
                   *data->get_dof_info(index).vector_partitioner),
               ExcInternalError());
      }
    BlockHelper::collect_sizes(vec);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize(
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
                                     data_,
    const std::vector<unsigned int> &given_row_selection,
    const std::vector<unsigned int> &given_column_selection)
  {
    data = data_;

    selected_rows.clear();
    selected_columns.clear();
    if (given_row_selection.empty())
      for (unsigned int i = 0; i < data_->n_components(); ++i)
        selected_rows.push_back(i);
    else
      {
        for (unsigned int i = 0; i < given_row_selection.size(); ++i)
          {
            AssertIndexRange(given_row_selection[i], data_->n_components());
            for (unsigned int j = 0; j < given_row_selection.size(); ++j)
              if (j != i)
                Assert(given_row_selection[j] != given_row_selection[i],
                       ExcMessage("Given row indices must be unique"));

            selected_rows.push_back(given_row_selection[i]);
          }
      }
    if (given_column_selection.size() == 0)
      selected_columns = selected_rows;
    else
      {
        for (unsigned int i = 0; i < given_column_selection.size(); ++i)
          {
            AssertIndexRange(given_column_selection[i], data_->n_components());
            for (unsigned int j = 0; j < given_column_selection.size(); ++j)
              if (j != i)
                Assert(given_column_selection[j] != given_column_selection[i],
                       ExcMessage("Given column indices must be unique"));

            selected_columns.push_back(given_column_selection[i]);
          }
      }

    edge_constrained_indices.clear();
    edge_constrained_indices.resize(selected_rows.size());
    edge_constrained_values.clear();
    edge_constrained_values.resize(selected_rows.size());
    have_interface_matrices = false;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize(
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
                                     data_,
    const MGConstrainedDoFs &        mg_constrained_dofs,
    const unsigned int               level,
    const std::vector<unsigned int> &given_row_selection)
  {
    std::vector<MGConstrainedDoFs> mg_constrained_dofs_vector(
      1, mg_constrained_dofs);
    initialize(data_, mg_constrained_dofs_vector, level, given_row_selection);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::initialize(
    std::shared_ptr<const MatrixFree<dim, value_type, VectorizedArrayType>>
                                          data_,
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs,
    const unsigned int                    level,
    const std::vector<unsigned int> &     given_row_selection)
  {
    AssertThrow(level != numbers::invalid_unsigned_int,
                ExcMessage("level is not set"));

    selected_rows.clear();
    selected_columns.clear();
    if (given_row_selection.empty())
      for (unsigned int i = 0; i < data_->n_components(); ++i)
        selected_rows.push_back(i);
    else
      {
        for (unsigned int i = 0; i < given_row_selection.size(); ++i)
          {
            AssertIndexRange(given_row_selection[i], data_->n_components());
            for (unsigned int j = 0; j < given_row_selection.size(); ++j)
              if (j != i)
                Assert(given_row_selection[j] != given_row_selection[i],
                       ExcMessage("Given row indices must be unique"));

            selected_rows.push_back(given_row_selection[i]);
          }
      }
    selected_columns = selected_rows;

    AssertDimension(mg_constrained_dofs.size(), selected_rows.size());
    edge_constrained_indices.clear();
    edge_constrained_indices.resize(selected_rows.size());
    edge_constrained_values.clear();
    edge_constrained_values.resize(selected_rows.size());

    data = data_;

    for (unsigned int j = 0; j < selected_rows.size(); ++j)
      {
        if (data_->n_cell_batches() > 0)
          {
            AssertDimension(level, data_->get_cell_iterator(0, 0, j)->level());
          }

        // setup edge_constrained indices
        std::vector<types::global_dof_index> interface_indices;
        mg_constrained_dofs[j]
          .get_refinement_edge_indices(level)
          .fill_index_vector(interface_indices);
        edge_constrained_indices[j].clear();
        edge_constrained_indices[j].reserve(interface_indices.size());
        edge_constrained_values[j].resize(interface_indices.size());
        const IndexSet &locally_owned =
          data->get_dof_handler(selected_rows[j]).locally_owned_mg_dofs(level);
        for (const auto interface_index : interface_indices)
          if (locally_owned.is_element(interface_index))
            edge_constrained_indices[j].push_back(
              locally_owned.index_within_set(interface_index));
        have_interface_matrices |=
          Utilities::MPI::max(
            static_cast<unsigned int>(edge_constrained_indices[j].size()),
            data->get_vector_partitioner()->get_mpi_communicator()) > 0;
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::set_constrained_entries_to_one(
    VectorType &dst) const
  {
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        const std::vector<unsigned int> &constrained_dofs =
          data->get_constrained_dofs(selected_rows[j]);
        for (const auto constrained_dof : constrained_dofs)
          BlockHelper::subblock(dst, j).local_element(constrained_dof) = 1.;
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          BlockHelper::subblock(dst, j).local_element(
            edge_constrained_indices[j][i]) = 1.;
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult(VectorType &      dst,
                                                    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    dst = Number(0.);
    vmult_add(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult_add(
    VectorType &      dst,
    const VectorType &src) const
  {
    mult_add(dst, src, false);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::Tvmult_add(
    VectorType &      dst,
    const VectorType &src) const
  {
    mult_add(dst, src, true);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::adjust_ghost_range_if_necessary(
    const VectorType &src,
    const bool        is_row) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    for (unsigned int i = 0; i < BlockHelper::n_blocks(src); ++i)
      {
        const unsigned int mf_component =
          is_row ? selected_rows[i] : selected_columns[i];
        // If both vectors use the same partitioner -> done
        if (BlockHelper::subblock(src, i).get_partitioner().get() ==
            data->get_dof_info(mf_component).vector_partitioner.get())
          continue;

        // If not, assert that the local ranges are the same and reset to the
        // current partitioner
        Assert(BlockHelper::subblock(src, i)
                   .get_partitioner()
                   ->locally_owned_size() ==
                 data->get_dof_info(mf_component)
                   .vector_partitioner->locally_owned_size(),
               ExcMessage(
                 "The vector passed to the vmult() function does not have "
                 "the correct size for compatibility with MatrixFree."));

        // copy the vector content to a temporary vector so that it does not get
        // lost
        LinearAlgebra::distributed::Vector<Number> copy_vec(
          BlockHelper::subblock(src, i));
        this->data->initialize_dof_vector(
          BlockHelper::subblock(const_cast<VectorType &>(src), i),
          mf_component);
        BlockHelper::subblock(const_cast<VectorType &>(src), i)
          .copy_locally_owned_data_from(copy_vec);
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::preprocess_constraints(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    adjust_ghost_range_if_necessary(src, false);
    adjust_ghost_range_if_necessary(dst, true);

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            edge_constrained_values[j][i] = std::pair<Number, Number>(
              BlockHelper::subblock(src, j).local_element(
                edge_constrained_indices[j][i]),
              BlockHelper::subblock(dst, j).local_element(
                edge_constrained_indices[j][i]));
            BlockHelper::subblock(const_cast<VectorType &>(src), j)
              .local_element(edge_constrained_indices[j][i]) = 0.;
          }
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::mult_add(
    VectorType &      dst,
    const VectorType &src,
    const bool        transpose) const
  {
    AssertDimension(dst.size(), src.size());
    AssertDimension(BlockHelper::n_blocks(dst), BlockHelper::n_blocks(src));
    AssertDimension(BlockHelper::n_blocks(dst), selected_rows.size());
    preprocess_constraints(dst, src);
    if (transpose)
      Tapply_add(dst, src);
    else
      apply_add(dst, src);
    postprocess_constraints(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::postprocess_constraints(
    VectorType &      dst,
    const VectorType &src) const
  {
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        const std::vector<unsigned int> &constrained_dofs =
          data->get_constrained_dofs(selected_rows[j]);
        for (const auto constrained_dof : constrained_dofs)
          BlockHelper::subblock(dst, j).local_element(constrained_dof) +=
            BlockHelper::subblock(src, j).local_element(constrained_dof);
      }

    // reset edge constrained values, multiply by unit matrix and add into
    // destination
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            BlockHelper::subblock(const_cast<VectorType &>(src), j)
              .local_element(edge_constrained_indices[j][i]) =
              edge_constrained_values[j][i].first;
            BlockHelper::subblock(dst, j).local_element(
              edge_constrained_indices[j][i]) =
              edge_constrained_values[j][i].second +
              edge_constrained_values[j][i].first;
          }
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult_interface_down(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    AssertDimension(dst.size(), src.size());
    adjust_ghost_range_if_necessary(src, false);
    adjust_ghost_range_if_necessary(dst, true);

    dst = Number(0.);

    if (!have_interface_matrices)
      return;

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
        {
          edge_constrained_values[j][i] = std::pair<Number, Number>(
            BlockHelper::subblock(src, j).local_element(
              edge_constrained_indices[j][i]),
            BlockHelper::subblock(dst, j).local_element(
              edge_constrained_indices[j][i]));
          BlockHelper::subblock(const_cast<VectorType &>(src), j)
            .local_element(edge_constrained_indices[j][i]) = 0.;
        }

    apply_add(dst, src);

    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        unsigned int c = 0;
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            for (; c < edge_constrained_indices[j][i]; ++c)
              BlockHelper::subblock(dst, j).local_element(c) = 0.;
            ++c;

            // reset the src values
            BlockHelper::subblock(const_cast<VectorType &>(src), j)
              .local_element(edge_constrained_indices[j][i]) =
              edge_constrained_values[j][i].first;
          }
        for (; c < BlockHelper::subblock(dst, j).locally_owned_size(); ++c)
          BlockHelper::subblock(dst, j).local_element(c) = 0.;
      }
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::vmult_interface_up(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    AssertDimension(dst.size(), src.size());
    adjust_ghost_range_if_necessary(src, false);
    adjust_ghost_range_if_necessary(dst, true);

    dst = Number(0.);

    if (!have_interface_matrices)
      return;

    VectorType src_cpy(src);
    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      {
        unsigned int c = 0;
        for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
          {
            for (; c < edge_constrained_indices[j][i]; ++c)
              BlockHelper::subblock(src_cpy, j).local_element(c) = 0.;
            ++c;
          }
        for (; c < BlockHelper::subblock(src_cpy, j).locally_owned_size(); ++c)
          BlockHelper::subblock(src_cpy, j).local_element(c) = 0.;
      }

    apply_add(dst, src_cpy);

    for (unsigned int j = 0; j < BlockHelper::n_blocks(dst); ++j)
      for (unsigned int i = 0; i < edge_constrained_indices[j].size(); ++i)
        BlockHelper::subblock(dst, j).local_element(
          edge_constrained_indices[j][i]) = 0.;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::Tvmult(
    VectorType &      dst,
    const VectorType &src) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    dst = Number(0.);
    Tvmult_add(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  std::size_t
  Base<dim, VectorType, VectorizedArrayType>::memory_consumption() const
  {
    return inverse_diagonal_entries.get() != nullptr ?
             inverse_diagonal_entries->memory_consumption() :
             sizeof(*this);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  std::shared_ptr<const MatrixFree<
    dim,
    typename Base<dim, VectorType, VectorizedArrayType>::value_type,
    VectorizedArrayType>>
  Base<dim, VectorType, VectorizedArrayType>::get_matrix_free() const
  {
    return data;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  const std::shared_ptr<DiagonalMatrix<VectorType>> &
  Base<dim, VectorType, VectorizedArrayType>::get_matrix_diagonal_inverse()
    const
  {
    Assert(inverse_diagonal_entries.get() != nullptr &&
             inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    return inverse_diagonal_entries;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  const std::shared_ptr<DiagonalMatrix<VectorType>> &
  Base<dim, VectorType, VectorizedArrayType>::get_matrix_diagonal() const
  {
    Assert(diagonal_entries.get() != nullptr && diagonal_entries->m() > 0,
           ExcNotInitialized());
    return diagonal_entries;
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::Tapply_add(
    VectorType &      dst,
    const VectorType &src) const
  {
    apply_add(dst, src);
  }



  template <int dim, typename VectorType, typename VectorizedArrayType>
  void
  Base<dim, VectorType, VectorizedArrayType>::precondition_Jacobi(
    VectorType &                                                          dst,
    const VectorType &                                                    src,
    const typename Base<dim, VectorType, VectorizedArrayType>::value_type omega)
    const
  {
    Assert(inverse_diagonal_entries.get() && inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    inverse_diagonal_entries->vmult(dst, src);
    dst *= omega;
  }



  //------------------------- MGInterfaceOperator ------------------------------

  template <typename OperatorType>
  MGInterfaceOperator<OperatorType>::MGInterfaceOperator()
    : Subscriptor()
    , mf_base_operator(nullptr)
  {}



  template <typename OperatorType>
  void
  MGInterfaceOperator<OperatorType>::clear()
  {
    mf_base_operator = nullptr;
  }



  template <typename OperatorType>
  void
  MGInterfaceOperator<OperatorType>::initialize(const OperatorType &operator_in)
  {
    mf_base_operator = &operator_in;
  }



  template <typename OperatorType>
  template <typename VectorType>
  void
  MGInterfaceOperator<OperatorType>::vmult(VectorType &      dst,
                                           const VectorType &src) const
  {
#ifndef DEAL_II_MSVC
    static_assert(
      std::is_same<typename VectorType::value_type, value_type>::value,
      "The vector type must be based on the same value type as this "
      "operator");
#endif

    Assert(mf_base_operator != nullptr, ExcNotInitialized());

    mf_base_operator->vmult_interface_down(dst, src);
  }



  template <typename OperatorType>
  template <typename VectorType>
  void
  MGInterfaceOperator<OperatorType>::Tvmult(VectorType &      dst,
                                            const VectorType &src) const
  {
#ifndef DEAL_II_MSVC
    static_assert(
      std::is_same<typename VectorType::value_type, value_type>::value,
      "The vector type must be based on the same value type as this "
      "operator");
#endif

    Assert(mf_base_operator != nullptr, ExcNotInitialized());

    mf_base_operator->vmult_interface_up(dst, src);
  }



  template <typename OperatorType>
  template <typename VectorType>
  void
  MGInterfaceOperator<OperatorType>::initialize_dof_vector(
    VectorType &vec) const
  {
    Assert(mf_base_operator != nullptr, ExcNotInitialized());

    mf_base_operator->initialize_dof_vector(vec);
  }



  //-----------------------------MassOperator----------------------------------

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::MassOperator()
    : Base<dim, VectorType, VectorizedArrayType>()
  {}



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::compute_diagonal()
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    Assert((Base<dim, VectorType, VectorizedArrayType>::data.get() != nullptr),
           ExcNotInitialized());

    this->inverse_diagonal_entries =
      std::make_shared<DiagonalMatrix<VectorType>>();
    this->diagonal_entries = std::make_shared<DiagonalMatrix<VectorType>>();
    VectorType &inverse_diagonal_vector =
      this->inverse_diagonal_entries->get_vector();
    VectorType &diagonal_vector = this->diagonal_entries->get_vector();
    this->initialize_dof_vector(inverse_diagonal_vector);
    this->initialize_dof_vector(diagonal_vector);
    inverse_diagonal_vector = Number(1.);
    apply_add(diagonal_vector, inverse_diagonal_vector);

    this->set_constrained_entries_to_one(diagonal_vector);
    inverse_diagonal_vector = diagonal_vector;

    const unsigned int locally_owned_size =
      inverse_diagonal_vector.locally_owned_size();
    for (unsigned int i = 0; i < locally_owned_size; ++i)
      inverse_diagonal_vector.local_element(i) =
        Number(1.) / inverse_diagonal_vector.local_element(i);

    inverse_diagonal_vector.update_ghost_values();
    diagonal_vector.update_ghost_values();
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::apply_add(VectorType &      dst,
                                               const VectorType &src) const
  {
    Base<dim, VectorType, VectorizedArrayType>::data->cell_loop(
      &MassOperator::local_apply_cell, this, dst, src);
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  MassOperator<dim,
               fe_degree,
               n_q_points_1d,
               n_components,
               VectorType,
               VectorizedArrayType>::
    local_apply_cell(
      const MatrixFree<
        dim,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type,
        VectorizedArrayType> &                     data,
      VectorType &                                 dst,
      const VectorType &                           src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    FEEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>
      phi(data, this->selected_rows[0]);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::values);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_value(phi.get_value(q), q);
        phi.integrate(EvaluationFlags::values);
        phi.distribute_local_to_global(dst);
      }
  }


  //-----------------------------LaplaceOperator----------------------------------

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::LaplaceOperator()
    : Base<dim, VectorType, VectorizedArrayType>()
  {}



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::clear()
  {
    Base<dim, VectorType, VectorizedArrayType>::clear();
    scalar_coefficient.reset();
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    set_coefficient(
      const std::shared_ptr<Table<2, VectorizedArrayType>> &scalar_coefficient_)
  {
    scalar_coefficient = scalar_coefficient_;
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  std::shared_ptr<Table<2, VectorizedArrayType>>
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::get_coefficient()
  {
    Assert(scalar_coefficient.get(), ExcNotInitialized());
    return scalar_coefficient;
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::compute_diagonal()
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    Assert((Base<dim, VectorType, VectorizedArrayType>::data.get() != nullptr),
           ExcNotInitialized());

    this->inverse_diagonal_entries =
      std::make_shared<DiagonalMatrix<VectorType>>();
    this->diagonal_entries = std::make_shared<DiagonalMatrix<VectorType>>();
    VectorType &inverse_diagonal_vector =
      this->inverse_diagonal_entries->get_vector();
    VectorType &diagonal_vector = this->diagonal_entries->get_vector();
    this->initialize_dof_vector(inverse_diagonal_vector);
    this->initialize_dof_vector(diagonal_vector);

    this->data->cell_loop(&LaplaceOperator::local_diagonal_cell,
                          this,
                          diagonal_vector,
                           /*unused*/  diagonal_vector);
    this->set_constrained_entries_to_one(diagonal_vector);

    inverse_diagonal_vector = diagonal_vector;

    for (unsigned int i = 0; i < inverse_diagonal_vector.locally_owned_size();
         ++i)
      if (std::abs(inverse_diagonal_vector.local_element(i)) >
          std::sqrt(std::numeric_limits<Number>::epsilon()))
        inverse_diagonal_vector.local_element(i) =
          1. / inverse_diagonal_vector.local_element(i);
      else
        inverse_diagonal_vector.local_element(i) = 1.;

    inverse_diagonal_vector.update_ghost_values();
    diagonal_vector.update_ghost_values();
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::apply_add(VectorType &      dst,
                                                  const VectorType &src) const
  {
    Base<dim, VectorType, VectorizedArrayType>::data->cell_loop(
      &LaplaceOperator::local_apply_cell, this, dst, src);
  }

  namespace Implementation
  {
    template <typename VectorizedArrayType>
    bool
    non_negative(const VectorizedArrayType &n)
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        if (n[v] < 0.)
          return false;

      return true;
    }
  } // namespace Implementation



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    do_operation_on_cell(
      FEEvaluation<
        dim,
        fe_degree,
        n_q_points_1d,
        n_components,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type> &phi,
      const unsigned int cell) const
  {
    phi.evaluate(EvaluationFlags::gradients);
    if (scalar_coefficient.get())
      {
        Assert(scalar_coefficient->size(1) == 1 ||
                 scalar_coefficient->size(1) == phi.n_q_points,
               ExcMessage("The number of columns in the coefficient table must "
                          "be either 1 or the number of quadrature points " +
                          std::to_string(phi.n_q_points) +
                          ", but the given value was " +
                          std::to_string(scalar_coefficient->size(1))));
        if (scalar_coefficient->size(1) == phi.n_q_points)
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            {
              Assert(Implementation::non_negative(
                       (*scalar_coefficient)(cell, q)),
                     ExcMessage("Coefficient must be non-negative"));
              phi.submit_gradient((*scalar_coefficient)(cell, q) *
                                    phi.get_gradient(q),
                                  q);
            }
        else
          {
            Assert(Implementation::non_negative((*scalar_coefficient)(cell, 0)),
                   ExcMessage("Coefficient must be non-negative"));
            const VectorizedArrayType coefficient =
              (*scalar_coefficient)(cell, 0);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(coefficient * phi.get_gradient(q), q);
          }
      }
    else
      {
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            phi.submit_gradient(phi.get_gradient(q), q);
          }
      }
    phi.integrate(EvaluationFlags::gradients);
  }



  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    local_apply_cell(
      const MatrixFree<
        dim,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type,
        VectorizedArrayType> &                     data,
      VectorType &                                 dst,
      const VectorType &                           src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;
    FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> phi(
      data, this->selected_rows[0]);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        do_operation_on_cell(phi, cell);
        phi.distribute_local_to_global(dst);
      }
  }


  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename VectorType,
            typename VectorizedArrayType>
  void
  LaplaceOperator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorType,
                  VectorizedArrayType>::
    local_diagonal_cell(
      const MatrixFree<
        dim,
        typename Base<dim, VectorType, VectorizedArrayType>::value_type,
        VectorizedArrayType> &data,
      VectorType &            dst,
      const VectorType &,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    using Number =
      typename Base<dim, VectorType, VectorizedArrayType>::value_type;

    FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> phi(
      data, this->selected_rows[0]);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        VectorizedArrayType local_diagonal_vector[phi.static_dofs_per_cell];
        for (unsigned int i = 0; i < phi.dofs_per_component; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_component; ++j)
              phi.begin_dof_values()[j] = VectorizedArrayType();
            phi.begin_dof_values()[i] = 1.;
            do_operation_on_cell(phi, cell);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_component; ++i)
          for (unsigned int c = 0; c < phi.n_components; ++c)
            phi.begin_dof_values()[i + c * phi.dofs_per_component] =
              local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);
      }
  }


} // end of namespace MatrixFreeOperators


DEAL_II_NAMESPACE_CLOSE

#endif


