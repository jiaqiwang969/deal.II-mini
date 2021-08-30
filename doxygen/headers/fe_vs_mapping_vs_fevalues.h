//include/deal.II-translator/A-headers/fe_vs_mapping_vs_fevalues_0.txt
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
 *    @defgroup FE_vs_Mapping_vs_FEValues How Mapping, FiniteElement, and
 * FEValues work together <h2>Introduction</h2>
 * 大多数人只创建一次有限元（以及可能的映射）对象，但实际上从未调用过任何成员函数
 *
 * - 他们只是通过FEValues接口使用它们进行装配。大多数人唯一的其他互动是通过读取 FiniteElementData::dofs_per_cell 变量，但这也只是在构造时设置的。换句话说，人们从来没有观察到FiniteElement或Mapping对象实际上是<i>do</i>的东西。
 *
 * - 而这完全是设计好的。
 * 因此，本文档是为那些对编写有限元或映射类感兴趣，并想了解FEValues如何工作以及与FiniteElement和Mapping类互动的人准备的。在下文中，我们将不对FEValues（作用于单元）、FEFaceValues（作用于面）和FESubfaceValues（作用于单元的面的子女）进行区分，因为它们在概念上都是一样的。因此，在下面的文字中，"FEValues
 * "这个术语将被普遍用于所有这三个类。
 *
 *  <h2>Who is responsible for what?</h2>
 * 在详细介绍数据和控制流之前，让我们定义哪个类负责提供什么样的信息。
 * <h3>%FEValues objects</h3>
 * FEValues是一个抽象的概念，它来自于这样的观察：人们在有限元代码中所做的几乎所有事情都只需要在正交点评估有限元的形状函数。例如，可以用正交@f[
 *   A^K_{ij} = \sum_q \nabla \varphi_i(\bf x_q) \cdot \nabla \varphi_j(\bf x_q) \;
 *   |\text{det}\; J(\bf x_q)| w_q,
 * @f]对@f[
 * A^K_{ij} = \int_K \nabla \varphi_i(\bf x) \cdot \nabla \varphi_j(\bf x) \;
 * dx
 * @f]形式的积分进行逼近，但在想要生成图形输出时也同样有效：在那里我们只需要知道网格顶点的有限元场的值，这也可以写成在正交点评估一切
 *
 * - 这些正交点就是单元的顶点（例如由QTrapez提供）。
 * FEValues的作用是为用户提供形状函数的值，以及它们的梯度，等等，在正交点。一些几何信息也是如此，例如，正交点的法向量。为此，它在FEValuesBase基类中提供了大量的成员函数，允许用户查询有关形状函数和几何信息的所有信息，但仅限于FEValues对象被初始化的正交点。
 * FEValues本身并不实际计算这些信息。它实际上只是提供了一个存储信息的地方，然后协调映射和有限元类之间的互动，让它们计算所要求的信息并将结果存储在FEValues提供的位置。
 * 最后，请记住，FEValues可以提供一系列令人难以置信的信息，但几乎所有的信息在任何情况下都是不必要的。例如，为了计算上述积分，不需要知道形状函数的二阶导数，也不需要知道正交点的法向量。为此，FEValues在与Mapping和FiniteElement类的交互中使用UpdateFlags来确定实际需要计算的内容。这在
 * @ref UpdateFlags  中有稍微详细的讨论。
 *
 *  <h3>Mappings</h3>
 * 映射（即从映射基类派生的类）负责与从参考（单位）单元
 * $[0,1]^\text{dim}$  到每个实际单元  $K\subset{\mathbb
 * R}^\text{spacedim}$  的映射有关的一切。这是由一个映射函数
 * $\mathbf F_K:[0,1]^\text{dim} \mapsto K$
 * 促成的。因此，映射类实现了一些接口，允许评估
 * $\mathbf F_K$  从参考单元向前映射点  $\hat{\mathbf x}$  到  $K$
 * ，并使用  $\mathbf F_K^{-1}$
 * 从实际单元向后映射到参考单元。映射提供的其他常见操作是将向量（你可以认为是连接到参考单元上的点
 * $\hat{\mathbf x}$
 * 并指向某些方向的向量）映射到真实单元上的等效向量。例如，这就是我们需要对形状函数的梯度所做的工作：这些是定义在参考单元上的向量，我们需要将这些梯度映射到实数单元上
 * $K$
 * 。类似的操作也可以为矩阵（等级为2的张量，而不是等级为1的向量）和高阶张量定义。
 * 许多这样的映射不仅需要映射 $\mathbf F_K$
 * 本身，还需要这个映射的梯度，通常被称为雅各布
 * $J_K=\hat\nabla \mathbf F_K$ ，以及高阶导数。
 * 由于FEValues只需要在正交点评估这些东西，所以映射一般不需要提供在<i>arbitrary</i>点评估的能力。相反，正如我们将在下面看到的，它们将被初始化为使用在参考单元上定义的一组正交点，然后为一个特定的单元
 * "重新初始化"，然后所有进一步的操作将只需要在真实单元上的这些正交点评估
 * $\mathbf F_K$ 。
 * 映射类具有双重作用：(i)计算几何信息（如法向量、雅各布定理等），并将其放入数据结构中，FEValues可以将其提供给用户；(ii)提供有限元所需的支持，将形状函数及其导数从参考单元映射到实际单元。
 *
 *  <h3>Finite elements</h3>
 * 有限元类（即从FiniteElement派生的类）负责在参考单元上定义其形状函数、导数和许多其他方面，但也负责在实际单元上计算映射值和导数（显然是在映射对象的帮助下）。在目前的讨论中，只有后一个角色是重要的。
 * 与映射一样，这里对我们来说重要的是有限元类可以在给定的正交点上提供这些信息，并且它们可以将计算的信息放入FEValues提供的结构中，然后FEValues的成员函数可以通过FEValuesBase中的成员函数将其传递给用户。
 *
 *  <h2>What to compute?</h2>
 * 假设用户想要计算形状函数的梯度，比如说计算上面的积分。然后他们会通过给出update_gradients标志来初始化一个FEValues对象（从
 * step-3
 * 开始，基本上每个教程程序都会这样做）。这表明用户希望FEValues对象能够提供真实单元上形状函数的梯度，但没有表示希望得到任何其他信息。
 * 然后，FEValues将首先找出映射对象和有限元对象之间的实际需求，以实现这一目标。这在运行FEValues构造函数时已经发生了。因为映射不依赖于有限元（尽管后者依赖于前者），FEValues首先通过
 * FiniteElement::requires_update_flags()
 * 询问有限元需要哪些<i>other</i>的信息来实现用户的请求。例如，如果有限元是FE_Q类型，那么它将确定为了计算实单元
 * $K$
 * 上形状函数的梯度，它需要计算参考单元上形状函数的梯度（这是它自己可以做到的，不需要任何外部帮助），但是这些参考梯度必须在每个正交点上乘以映射的雅各布系数的逆值
 * $J^{-1}_K$ 。这个乘法通常被称为<i>covariant
 * transformation</i>，因此FE_Q的 FiniteElement::requires_update_flags()
 * 函数的实现（在中间类FE_Poly中提供）将同时返回原始的update_gradients标志以及update_covariant_transformation。
 * 在第二步中，FEValues对象将调用映射中的相应函数，
 * Mapping::requires_update_flags()
 * 以确定提供update_gradients和update_covariant_transformation的要求。前者不在映射的范围内，所以被忽略了。后者通常需要先计算雅各布矩阵
 * $J_K$ ，一个典型的映射类将通过在列表中添加
 * update_contravariant_transformation 来表示。
 *
 *  <h2>Pre-computing things</h2>
 * 此时，FEValues对象已经找出了一套完整的标志，表明大家要计算什么来满足用户的要求。下一步，仍然是在构建FEValues对象的过程中，源于这样的认识：许多东西可以预先计算一次，然后在我们每次移动到一个真正的单元时重复使用。一个例子是，为了计算真实单元上形状函数的梯度，我们需要知道参考单元上形状函数的梯度（在参考单元上的正交点），而且这些梯度总是相同的：每次我们访问一个新单元时，这些值都会保持不变，所以每次都重新计算它们是低效的。对于一些映射类所计算的一些信息也可以提出类似的论点。
 * 因此，FEValues对象同时初始化了映射和它所指向的有限元对象，使用正交对象和上一节所述的最后一组更新标志来计算。这个初始化包括预先计算这些类在给定更新标志集后可以预先计算的内容，然后存储这些信息供以后使用。
 * 然后问题来了：在哪里存储这些信息。在实践中，我们不希望将这些信息存储在映射或有限元对象本身，因为这意味着（i）一次只能有一个FEValues对象使用任何给定的映射或有限元对象，以及（ii）这些对象不能在多线程环境下使用。
 * 相反，该方法是这样的。
 *
 * - FEValues调用 Mapping::get_data() （以及FEFaceValues调用 Mapping::get_face_data(), 和FESubfaceValues调用 Mapping::get_subface_data()) ，带有正交对象和最后一组更新标志。从Mapping派生出来的类中的这些函数的实现将分配一个从 Mapping::InternalDataBase 派生出来的类型的对象，在那里他们基本上可以存储他们认为有用的任何东西，以便以后重新使用。   Mapping::InternalDataBase 本身实际上并没有提供任何重要的成员变量，但真正留给派生类的是他们认为在这个时候他们可以有效地预先计算和存储的东西。如果一个映射没有什么需要预先计算的（或者映射类的作者很懒，不想考虑什么可能被预先计算），那么这样的类将简单地从 Mapping::InternalDataBase 中派生出它自己的InternalData对象，而没有实际添加任何成员变量。
 * 这样产生的对象就会被返回到FEValues中的调用站点，并由FEValues对象存储。以后每当FEValues对象想要从映射中获得任何信息时，它都会被交还，从而为映射对象提供了读取其先前存储的数据的能力。
 *
 *
 * - 其次，FEValues也调用 FiniteElement::get_data() （FEFaceValues调用 Mapping::get_face_data(), ，FESubfaceValues调用 Mapping::get_subface_data()), ，再次调用正交对象和最后一组更新标志。这些函数的作用与它们在映射中的对应函数基本相同，而且这样初始化的对象，这次是源自 FiniteElement::InternalDataBase, 的类型，每当FEValues对象在以后的时间里想要从有限元对象那里得到一些东西时，总是会被反馈给有限元。
 * 这种方法允许我们同时从多个FEValues对象中使用有限元和映射对象，也可能同时从多个线程中使用。重点是，有限元或映射对象的每个用户都会持有他们自己的、唯一的、从
 * <code>get_data()</code>
 * 函数返回的对象，而且所有发生的事情都发生在这些对象上，而不是映射或有限元对象本身的成员变量上。
 *
 *  <h2>Computing on a given cell</h2>
 * 之前的所有步骤都发生在创建FEValues对象的时候。到此为止，我们所做的都是设置数据结构，但从用户的角度来看，到目前为止还没有计算出任何有用的东西。这只发生在
 * FEValues::reinit()  在具体单元格上被调用时  $K$  。
 * 然后FEValues所做的事情是，按照这个顺序。
 *
 *
 * - FEValues计算出该单元是否是调用了 FEValues::reinit() 的前一个单元的平移或其他类似的简单转换。这个结果存储在 CellSimilarly::Similarity 对象中，然后将被传递给映射和有限元，以潜在地简化一些计算。例如，如果当前单元只是前一个单元的平移，那么就不需要重新计算映射的雅各布矩阵 $J_K$ （或其逆），因为它将与前一个单元相同。
 *
 *
 * - 接下来， FEValues::reinit() 调用 Mapping::fill_fe_values() （显然，FEFaceValues调用 Mapping::fill_fe_face_values() ，FESSubfaceValues调用 Mapping::fill_fe_subface_values()).  这个函数的参数包括我们被要求访问的单元格（或面，或子面），以及上面的单元格相似性参数，对我们之前从 Mapping::get_data(), 获得的对象的引用，以及对 internal::FEValues::MappingRelatedData 类型的对象的引用，映射应该把它的结果写入其中。特别是，它需要计算之前由更新标志指定的所有映射相关信息，然后将它们写进输出对象。   在输出对象中，映射需要填充的字段的例子是JxW值的计算、雅各布矩阵及其反值的计算，以及单元格（如果dim小于spacedim）和面的法向量。
 *
 *
 * - 最后， FEValues::reinit() 调用 FiniteElement::fill_fe_values() （显然，FEFaceValues调用 FiniteElement::fill_fe_face_values() ，FESSubfaceValues调用 FiniteElement::fill_fe_subface_values()).  这个函数的参数包括我们被要求访问的单元格（或面，或子面），以及上面的单元格相似性参数，对我们之前从 FiniteElement::get_data(), 获得的对象的引用，以及对 internal::FEValues::MappingRelatedData 类型的对象的引用，映射应该将其结果写入其中。
 * 除了这些， FiniteElement::fill_fe_values()
 * 函数还接收对正在使用的映射对象的引用，以及我们之前从
 * Mapping::get_data(). 中收到的 Mapping::InternalDataBase
 * 对象。原因是，通常，有限元希望将形状函数的值或梯度从参考单元映射到实际单元，而这些映射由各种
 * Mapping::transform() 函数提供便利
 *
 * - 这都需要对FEValues对象先前从映射中获得的内部对象的引用。这可能最好是通过查看实际代码来理解，在 FE_Poly::fill_fe_values(), 中可以找到一个简单而有启发性的例子，该函数适用于一般标量、多项式有限元基。
 * 与映射一样， FiniteElement::fill_fe_values()
 * 函数然后使用它们之前在构建FEValues对象时计算的任何信息（即当它调用
 * FiniteElement::get_data()),
 * 时，使用这个和映射中的函数来计算更新标志所指定的任何请求。
 * 这一切完成后，我们终于可以向FEValues的所有者提供对最初通过更新标志要求的字段的访问。
 *
 * @ingroup feall
 *
 */


