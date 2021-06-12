//include/deal.II-translator/A-headers/concepts_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


/**
 *    @defgroup Concepts Concepts, or expectations on template parameters
 * 有时对一个对象的类型施加约束而不要求它属于一个特定的继承层次是有用的。这些通常在C++社区中被称为
 * <em> 概念 </em>
 * 。本模块列出了deal.II中常用的概念，并对其意图进行了简要描述。在deal.II中列出类型的约束条件的惯例是在模板中提供概念的名称作为
 * <code>typename</code>
 * ：例如，矢量的类型取决于底层字段的类型，因此它被定义为模板。
 * @code
 * template <typename Number>
 * class Vector;
 * @endcode
 * 这里的重点是，你正在创建一个可以存储 @p Number.
 * 类型元素的向量，但这上面有一些基本假设。例如，deal.II
 * Vector类并不打算仅仅作为一个集合使用（与
 * <code>std::vector</code>
 * 不同），而是定义了向量空间的操作，如向量的加法，或向量的规范。因此，用户可以为
 * @p Number
 * 指定的数据类型必须满足某些条件（即，它必须符合或
 * "模拟 "一个
 * "概念"）。具体来说，该类型必须表示代表数学上称为
 * "场
 * "的元素的对象（你可以认为是，嗯，"数字"：我们可以加、乘、除、取绝对值的东西，等等）。概念的意义在于描述
 * <em> 一个类型必须满足哪些条件 </em>
 * 才能在特定的环境中成为有效的模板参数。
 * 本页描述了在整个deal.II中使用的一些概念的这些条件。具体来说，在上面的例子中，下面讨论的 @ref ConceptNumber "数字概念 "
 * 描述了可以作为向量类的参数的类型。
 * 概念作为C++的语言扩展已经被提议了很久了。它们将允许我们描述一个类或函数具有某些属性，以便成为一个合格的模板参数。例如，它允许我们在C++代码中表达，例如，
 * GridTools::find_closest_vertex(),
 * 的第一个参数必须有一个代表实际网格的类型。
 *
 * - 目前我们只能用文字来描述，见下文。使用C++的概念将允许我们在代码中描述这一点，并且试图用一个对象作为第一参数来调用这样一个函数，但实际上它不是一个网格，这将产生一个编译器错误，使不匹配的情况变得清晰。
 * 不幸的是，这些对C++的建议从未进入任何官方的C++标准；然而，它们被建议用于C++20。一旦我们的绝大多数用户拥有支持这一标准的编译器，我们可能会开始使用它们。
 * 关于这个主题的更多信息可以在<a
 * href="https://en.wikipedia.org/wiki/Concepts_(C%2B%2B)">this wikipedia
 * page</a>找到。
 *
 *  <dl> <dt class="concepts">  @anchor  ConceptDoFHandlerType
 * <b>DoFHandlerType</b></dt>。 <dd>  deal.II包括DoFHandler和
 * hp::DoFHandler
 * 这两个管理网格上自由度的对象。虽然两者没有任何继承关系，但它们足够相似，许多函数只需要类似于DoFHandler的东西就能正常工作。
 * </dd> <dt class="concepts">  @anchor  ConceptMatrixType
 * <b>MatrixType</b></dt>。
 * <dd>  deal.II中的许多函数和类需要一个知道如何计算矩阵-向量积（成员函数  <code>vmult</code>  ）、转置矩阵-向量积（成员函数  <code>Tvmult</code>  ）以及 "乘法和加法 "
 * 等效物的对象  <code>vmult_add</code> and <code>Tvmult_add</code>
 * 。有些函数只需要 <code>vmult</code> and <code>Tvmult</code>
 * ，但如果模板需要MatrixType参数的话，一个对象应该实现所有四个成员函数。编写满足这些条件的类是很常见的，所以编写LinearOperator类是为了让事情变得更简单；更多信息请参见
 * @ref LAOperators  。 对 <code>MatrixType</code>
 * 的一种看法是，假装它是一个具有如下签名的基类（这几乎就是SparseMatrix提供的接口）。
 *
 * @code
 * class MatrixType
 * {
 * public:
 * template <typename VectorType>
 * virtual void vmult(VectorType &u, const VectorType &v) const =0;
 *
 * template <typename VectorType>
 * virtual void Tvmult(VectorType &u, const VectorType &v) const =0;
 *
 * template <typename VectorType>
 * virtual void vmult_add(VectorType &u, const VectorType &v) const =0;
 *
 * template <typename VectorType>
 * virtual void Tvmult_add(VectorType &u, const VectorType &v) const =0;
 * };
 * @endcode
 *
 * C++中的模板函数不能是虚拟的（这是deal.II中不使用这种方法的主要原因），所以用继承的方式实现这个接口是行不通的，但这仍然是思考这个模板概念的一个好方法。人们可以使用LinearOperator类来实现
 * <code>vmult_add</code> and <code>Tvmult_add</code>
 * ，而不是手动实现它们。   </dd> <dt class="concepts">  @anchor
 * ConceptMeshType <b>MeshType</b></dt>。 <dd>
 * 网格可以被认为是顶点和连接点的数组，但一个更有成效的观点是将它们视为<i>collections
 * of
 * cells</i>。在C++中，集合通常被称为<i>containers</i>（典型的容器是
 * std::vector,  std::list,
 * 等），它们的特点是能够对集合中的元素进行迭代。<tt>MeshType</tt>概念是指任何定义了适当方法（如
 * DoFHandler::begin_active()) 和<tt>typedefs</tt>（如
 * DoFHandler::active_cell_iterator)
 * ）的容器，用于管理单元的集合。 Triangulation、DoFHandler和
 * hp::DoFHandler
 * 的实例都可以被视为单元格的容器。事实上，这些类的公共接口中最重要的部分仅仅由获得其元素的迭代器的能力组成。由于接口的这些部分是通用的，也就是说，这些函数在所有的类中都有相同的名字，所以可以编写一些操作，这些操作实际上并不关心它们是在三角化还是在DoF处理对象上工作。例如，在GridTools命名空间中就有大量的例子，强调了网格和DoF处理程序都可以简单地被视为单元格的集合（容器）这一抽象的力量。
 * 另一方面，网格是不同于 std::vector 或 std::list
 * 的非标准容器，因为它们可以以多种方式被切分。例如，我们可以在活动单元的子集上迭代，也可以在所有单元上迭代；同样，单元被组织成层次，我们可以只为一个层次上的单元获得迭代器范围。然而，一般来说，所有实现单元格容器概念的类都使用相同的函数名称来提供相同的功能。
 * 可以用任何一个类来调用的函数通过接受一个模板参数来表示，例如
 * @code
 * template <template <int, int> class MeshType>
 * @endcode
 * 或
 * @code
 * template <typename MeshType>
 * @endcode
 * 满足这个概念的类被统称为  <em>  网格类  </em>  。<tt>MeshType</tt>的确切定义在很大程度上依赖于库的内部结构，但它可以被概括为具有以下属性的任何类。   <ol>   <li>  一个名为<tt>typedef</tt>的<tt>active_cell_iterator</tt>。     </li>   <li>  一个<tt>get_triangulation()</tt>方法，返回单元格集合的基础几何描述（Triangulation类之一）的引用。如果网格恰好是一个Triangulation，那么网格只是返回一个对其本身的引用。     </li>   <li>  一个方法<tt>begin_active()</tt>，返回一个指向第一个活动单元的迭代器。     </li>   <li>  一个静态成员值<tt>dimension</tt>，包含对象所处的维度。     </li>   <li>  一个静态成员值<tt>space_dimension</tt>，包含对象的维度（例如，一个2D表面在3D环境中会有<tt>space_dimension = 2</tt>）。     </li>   </ol>   </dd> 。
 * <dt class="concepts">  @anchor  ConceptNumber <b>Number</b></dt>。 <dd>
 * 这个概念描述了作为向量或矩阵条目有意义的标量，这通常是场元素的一些有限精度的近似。典型的例子是
 * <code>double</code> 和 <code>float</code>, but deal.II supports
 * <code>std::complex&lt;T&gt;</code> 的浮点类型 <code>T</code>
 * 在很多地方也是如此。   </dd> <dt class="concepts">  @anchor
 * ConceptPolynomialType <b>PolynomialType</b></dt>。
 * <dd>
 * 更多信息见 @ref Polynomials
 * 中的描述。在某些情况下，任何满足类似于接口的
 * @code
 * template <int dim>
 * class PolynomialType
 * {
 * virtual void compute (const Point<dim>            &unit_point,
 *                       std::vector<Tensor<1,dim> > &values,
 *                       std::vector<Tensor<2,dim> > &grads,
 *                       std::vector<Tensor<3,dim> > &grad_grads) const =0;
 * }
 * @endcode
 *
 * 为了实现有限元的目的，可将其视为多项式。   </dd> <dt
 * class="concepts">  @anchor  ConceptPreconditionerType
 * <b>PreconditionerType</b></dt>。 <dd>  这基本上是
 * <code>MatrixType</code> 的同义词，但通常只要求定义
 * <code>vmult()</code> and <code>Tvmult()</code> 。大多数时候，定义
 * <code>Tvmult()</code> 是不必要的。人们应该认为
 * <code>vmult()</code>
 * 是将线性算子的逆运算应用于向量的一些近似值，而不是线性算子对向量的作用，用于预处理器类。
 * </dd> <dt class="concepts">  @anchor  ConceptRelaxationType
 * <b>RelaxationType</b></dt>。
 * <dd>  这是一个能够对多网格方法进行放松的对象。我们可以认为满足这个约束的对象具有以下接口以及 @ref ConceptMatrixType  "MatrixType "
 * 所要求的约束。
 * @code
 * class RelaxationType
 * {
 * public:
 * template <typename VectorType>
 * virtual void step(VectorType &u, const VectorType &v) const =0;
 *
 * template <typename VectorType>
 * virtual void Tstep(VectorType &u, const VectorType &v) const =0;
 * };
 * @endcode
 * 其中这两个成员函数执行平滑方案的一个步骤（或这种步骤的转置）。换句话说，这些函数执行的操作是
 * $u = u
 *
 * - P^{-1} (A u
 *
 * - v)$  和  $u = u
 *
 * - P^{-T} (A u
 *
 * - v)$  。   </dd> <dt class="concepts">  @anchor
 * ConceptSparsityPatternType <b>SparsityPatternType</b></dt>。
 * <dd>
 * 几乎所有的函数（除了 SparsityTools::distribute_sparsity_pattern)
 * 这个明显的例外）都以稀疏性模式为参数，可以采用常规的SparsityPattern或DynamicSparsityPattern，甚至是块状稀疏性模式之一。更多信息见
 * @ref Sparsity  。   </dd> <dt class="concepts">  @anchor
 * ConceptStreamType <b>StreamType</b></dt>。 <dd>
 * 在C++中派生新的流类，众所周知是很困难的。为了解决这个问题，一些函数接受一个定义了
 * <code>operator&lt;&lt;</code>
 * 的参数，这样就可以轻松地输出到任何一种输出流。
 * </dd> <dt class="concepts">  @anchor  ConceptVectorType
 * <b>VectorType</b></dt>。 <dd>
 * deal.II支持许多不同的向量类，包括与其他库的向量的绑定。这些与标准库中的向量类似（即它们定义了
 * <code>begin()</code>, <code>end()</code>  ,  <code>operator[]</code>, and
 * <code>size()</code>  ），但也定义了像  <code>add()</code>
 * 这样的数字操作。VectorType的一些例子包括Vector,
 * TrilinosWrappers::MPI::Vector,  和 BlockVector。   </dd> </dl>
 *
 */


