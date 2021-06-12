//include/deal.II-translator/A-headers/update_flags_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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
 *    @defgroup UpdateFlags The interplay of UpdateFlags, Mapping, and
 * FiniteElement in FEValues <h2>Introduction</h2>
 * 为了计算单个单元对全局矩阵和右手边的局部贡献，我们通常采用两种技术。
 *
 *
 * - 首先，积分从实际单元 $K$ 转换到单位/参考单元 $\hat K$  。例如，对于拉普拉斯方程，我们将@f[
 * A^K_{ij} = \int_K \nabla \varphi_i(\bf x) \cdot \nabla \varphi_j(\bf x) \;
 * dx
 * @f]转换为@f[
 * A^K_{ij} = \int_{\hat K} \left[ J^{-1}(\hat{\bf x}) \hat \nabla
 * \varphi_i(\hat{\bf x}) \right] \cdot \left[ J^{-1}(\hat{\bf x}) \hat \nabla
 * \varphi_j(\hat{\bf x}) \right] \; |\textrm{det}\; J(\hat{\bf x})| \;\;
 * d\hat x, @f]，其中帽子表示参考坐标，而 $J(\hat{\bf x}_q)$
 * 是映射 $\bf x = \bf F_K(\hat{\bf x})$ 的雅各布 $\frac{\partial \bf
 * F_K(\hat{\bf x})}{\partial\bf \hat x}$  。
 *
 *
 * - 第二，这个积分然后通过正交法进行近似。这就产生了公式@f[
 * A^K_{ij} = \sum_{q}J^{-1}\left[(\hat{\bf x}_q) \hat \nabla
 * \varphi_i(\hat{\bf x}_q)\right] \cdot \left[J^{-1}(\hat{\bf x}_q) \hat
 * \nabla \varphi_j(\hat{\bf x}_q)\right]\ |\textrm{det}\ J(\hat{\bf x}_q)|
 * w_q, @f]，其中 $q$ 表示正交点的索引， $\hat{\bf x}_q$
 * 其在参考单元上的位置，以及 $w_q$ 其重量。
 * 为了在应用程序代码中评估这样的表达式，我们必须访问三种不同的对象：描述参考单元上正交点的位置
 * $\hat{\bf x}_q$ 和权重 $w_q$
 * 的正交对象；描述单元上形状函数梯度 $\hat\nabla
 * \varphi_i(\hat{\bf x}_q)$
 * 的有限元对象；以及提供雅各布系数和其行列式的映射对象。处理所有这些对象会很麻烦而且容易出错。
 * 另一方面，这三种对象几乎总是一起出现，事实上，除了一起使用正交、有限元或映射对象，deal.II的应用代码很少对它们做任何事情。由于这个原因，deal.II使用FEValues抽象，结合了形状函数、实际网格单元的几何信息和参考单元的正交规则。在构建时，它在上述三个类别中各取一个对象。之后，它可以为一个具体的网格单元
 * "重新初始化"，然后提供映射的正交点和权重，映射的形状函数值和导数，以及从参考单元到实际网格单元转换的一些属性。
 * 由于任何这些值的计算都是潜在的昂贵的（例如，当使用高阶元素的高阶映射时），FEValues类只计算它被明确要求的值。为此，它在构建时需要一个UpdateFlags类型的标志列表，指定每次访问一个单元时应该更新哪些数量。在上面的例子中，你想要实数单元上形状函数的梯度，这由标志
 * <code>update_gradients</code>
 * 编码，以及雅各布系数的行列式乘以正交权重的乘积，这用术语
 * <code>JxW</code> and encoded in the flag <code>update_JxW_values</code>
 * 记忆性地编码了。因为这些标志是由整数中的单比特表示的，产生一个<i>set
 * of flags</i>相当于在一个整数中设置多个比特，这是用操作
 * <code>update_gradients | update_JxW_values</code>
 * 来促进的（换句话说，也许有点混乱，操作 @"this
 * 操作<i>and</i>，操作 @" 是由表达
 * @"single-bit-in-an-integer-for-this-operation <i>binary-or</i>
 * 单位-整数-该操作  @").
 * 为了使操作更便宜，FEValues和它所依赖的映射和有限元对象实际上只计算你在更新标志中指定的那些信息（加上一些计算指定内容所需的信息，见下文），而不是可能在一个单元上计算的所有信息。这种优化使得在单元格上进行迭代装配的成本大大降低，但这也意味着我们应该注意提供尽可能少的标志集。
 * 此外，一旦你传递了一组你想要的标志，填充FEValues数据字段的函数就能够区分必须在每个单元上重新计算的值（例如映射梯度）和单元间不发生变化的量（例如不同单元上同一正交点的通常
 * $Q_p$
 * 有限元的形状函数值；但是这一特性对于Raviart-Thomas元素的形状函数不成立，它必须随本地单元旋转）。这允许进一步优化底层装配的计算。
 *
 *  <h2> Tracking dependencies </h2>
 * 假设你想计算如上所示的拉普拉斯矩阵。在这种情况下，你需要指定
 * <code>update_gradients</code> 标志（用于 $\nabla\varphi_i(\bf x_q)$
 * ）和 <code>update_JxW_values</code> 标志（用于计算
 * $|\textrm{det}\; J(\bf x_q)|w_q$
 * ）。然而，在内部，有限元要求计算完整的雅各布矩阵的逆，
 * $J^{-1}(\bf x_q)$
 * （而不仅仅是矩阵的行列式），为了计算雅各布矩阵的逆，还需要先计算雅各布矩阵。
 * 由于这些是对用户不重要的要求，所以没有必要在用户代码中指定。相反，给定一组更新标志，FEValues对象首先询问有限元对象需要计算哪些信息，以满足用户在更新标志中提供的要求。因此，有限元对象可以向更新标志添加其他标志（例如，在上面的例子中，FE_Q对象将向列表中添加
 * <code>update_covariant_transformation</code> ，因为这是从
 * $\hat\nabla\hat\varphi_i(\hat{\bf x}_q)$ 到 $\nabla\varphi_i(\bf x_q)$
 * 的必要转换）。有了这些更新的标志，FEValues就会通过调用
 * Mapping::requires_update_flags().
 * 来询问映射是否也要在列表中添加更多的标志，以满足用户和有限元对象的需要（这种先询问有限元，再询问映射的程序不需要迭代，因为映射从来不需要有限元类计算的信息，而有限元类通常需要映射计算的信息）。使用这个最终的列表，FEValues对象然后要求有限元对象和映射对象都创建临时结构，将一些可以一次性计算的临时信息存储到其中，这些标志将在以后我们访问的每个单元上重新计算数据时使用。
 *
 *  <h2>Update once or each</h2>
 * 如上所述，我们现在已经确定了满足用户所需信息的最后一套东西，这些信息是由他们提供的更新标志所传达的。然后，这些信息通常会在随后的整合循环中对用户代码访问的每个单元进行查询。
 * 鉴于许多映射或有限元类的计算都是潜在的昂贵的，FEValues采用了一个系统，鼓励映射和有限元对象预先计算那些无需参考具体单元就能计算的信息，并在要求访问网格的特定单元时利用这些信息。一个例子是，普通FE_Q元素的形状函数的值是在参考单元上定义的，而实际单元上的值正好是参考单元上的值
 *
 * 因此，没有必要对每个单元格的形状函数进行评估，只需在开始时进行一次评估，将数值存储在某个地方，当访问一个具体的单元格时，只需将这些数值从临时位置复制到输出结构中即可。(但是请注意，这是FE_Q元素所特有的：如果我们使用FE_RaviartThomas元素就不是这样了，因为在那里，计算一个单元上的形状函数值需要知道映射的Jacobian，这取决于我们访问的单元的几何形状；因此，对于这个元素，简单地复制预先计算的信息并不足以评估特定单元上的形状函数值。)
 * 为了适应这种结构，映射和有限元类都可以在内部将更新标志分成两组，通常被称为
 * <code>update_once</code> and <code>update_each</code>
 * （尽管这些名称没有出现在任何公共接口中）。前者包含所有那些在FEValues对象开始与映射或有限元交互时可以预先计算一次的信息，而后者则包含那些对应于需要在每个单元上计算的标志。例如，如果
 * <code>update_flags=update_values</code> ，那么FE_Q类将设置
 * <code>update_once=update_values</code> 和 <code>update_each=0</code>
 * ，而Raviart-Thomas元素将以相反的方式进行。
 * 这些标志集的目的是相互排斥的。另一方面，没有任何东西可以为映射或有限元类之外的东西提供这种分解。
 *
 * - 它是一个纯粹的内部分解。
 *
 *  <h2>Generation of the actual data</h2>
 * 如上所述，数据在两个不同的时间被计算：一次是在开始时在参考单元格上，另一次是每当我们移动到一个实际单元格时。接下来将讨论每个步骤所涉及的函数。
 *
 *  <h3>Initialization</h3>
 * 在我们还没有访问第一个真实单元之前，计算参考单元的数据是一个两步的过程。首先，FEValues、FEFaceValues和FESubfaceValues的构造函数分别需要让Mapping和FiniteElement对象建立内部数据结构。这些结构在以下意义上是内部的：FEValues对象要求有限元和映射对象各创建一个 FiniteElement::InternalDataBase 和 Mapping::InternalDataBase 类型的对象；如果实际的有限元和映射类希望存储一些超出这些基类已经提供的数据，实际上可以创建衍生类型的对象。这其中涉及的函数有  <ul>   <li>Mapping::get_data()   <li>Mapping::get_face_data()   <li>Mapping::get_subface_data()   <li>FiniteElement::get_data()   <li>FiniteElement::get_face_data()   <li>FiniteElement::get_subface_data()   </ul>  。
 * 然后，FEValues对象接管了这些对象的所有权，并将在FEValues对象的生命周期结束时销毁它们。之后，FEValues对象要求FiniteElement和Mapping对象向这些InternalDataBase对象中填充与参考单元上可以和需要计算的内容有关的数据。这是在这些函数中完成的。   <ul>   <li>FEValues::initialize()   <li>FEFaceValues::initialize()   <li>FESubfaceValues::initialize()   </ul> 。
 *
 *  <h3>Reinitialization for a mesh cell</h3>
 * 一旦初始化结束，我们调用 FEValues::reinit,  FEFaceValues::reinit 或 FESubfaceValues::reinit 来移动到一个具体的单元或面，我们需要计算 "update_each "的各种数据。这是在以下函数中完成的。   <ul>   <li>FEValues::reinit()  调用 Mapping::fill_fe_values(),  然后 FiniteElement::fill_fe_values()   <li>FEFaceValues::reinit()  调用 Mapping::fill_fe_face_values(),  然后 FiniteElement::fill_fe_face_values()   <li>FESubfaceValues::reinit()  调用 Mapping::fill_fe_subface_values(),  然后 FiniteElement::fill_fe_subface_values()   </ul>  。
 * 这是计算存储在 internal::FEValues::MappingRelatedData 和
 * internal::FEValues::FiniteElementRelatedData
 * 对象中的FEValues实际数据字段的地方。这些函数首先调用Mapping中的函数，这样，有限元所需的所有映射数据都可以得到。然后，调用FiniteElement函数。
 *
 * @ingroup feall
 *
 */


