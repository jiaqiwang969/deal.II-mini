//include/deal.II-translator/A-headers/constraints_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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
 *    @defgroup constraints Constraints on degrees of freedom
 * @ingroup dofs
 * 本模块处理自由度的约束问题。处理约束的中心类是AffineConstraints类。
 * 约束通常来自几个方面，例如。
 *
 *
 * - 如果你有迪里希特型边界条件， $u|_{\partial\Omega}=g$ ，通常通过要求边界上的自由度有特定的值来执行，例如 $x_{12}=42$ ，如果边界条件 $g(\mathbf x)$ 要求自由度12位置的有限元解 $u(\mathbf x)$ 有42值。这样的约束是由那些接受AffineConstraints参数的 VectorTools::interpolate_boundary_values 函数版本产生的（尽管也有其他处理Dirichlet条件的方法，使用 MatrixTools::apply_boundary_values, 见例如 step-3 和 step-4  ）。
 *
 *
 * - 如果你有边界条件，设定了解决方案的某一部分数值，例如没有法向通量， $\mathbf n \cdot
 * \mathbf u=0$ （如发生在流动问题中，由
 * VectorTools::compute_no_normal_flux_constraints
 * 函数处理）或规定的切向分量， $\mathbf{n}\times\mathbf{u}=
 * \mathbf{n}\times\mathbf{f}$ （如发生在电磁问题中，由
 * VectorTools::project_boundary_values_curl_conforming
 * 函数处理）。对于前一种情况，例如，设想我们在法线矢量具有
 * $\frac 1{\sqrt{14}} (1,2,3)^T$
 * 形式的顶点，在这个顶点的流场的 $x$ -、 $y$ -和 $z$
 * 分量与自由度12、28和40相关。那么无正态流条件意味着我们需要有条件
 * $\frac 1{\sqrt{14}} (x_{12}+2x_{28}+3x_{40})=0$  。
 * 规定的切向分量会导致类似的约束，尽管右手边经常有一些东西。
 *
 *
 * - 如果你有悬挂节点约束，例如在这样的网格中。          @image html hanging_nodes.png ""  我们假设右下角的两个红色自由度之一是 $x_{12}$ ，其左右两边的黄色邻居是 $x_{28}$ 和 $x_{40}$  。那么，要求有限元函数是连续的，就相当于要求 $x_{12}=
   \frac 12 (x_{28}+x_{40})$  。类似的情况发生在hp自适应有限元方法的背景下。   例如，当在网格的两个标记单元上使用Q1和Q2元素（即使用FE_Q(1)和FE_Q(2)）时  @image html hp-refinement-simple.png  有三个约束：首先  $x_2=\frac 12 x_0 + \frac 12 x_1$  ，然后  $x_4=\frac 14 x_0 + \frac 34 x_1$  ，最后是身份  $x_3=x_1$  。即使所有的单元格都使用相同的有限元，也会出现类似的约束条件作为悬挂节点。在所有这些情况下，你将使用 DoFTools::make_hanging_node_constraints 函数来计算这种约束。
 *
 *
 * - 其他线性约束，例如，当你试图为一个问题施加某个平均值时，否则就没有唯一的解决方案。在 step-11 的教程程序中给出了这样一个例子。
 * 在所有这些例子中，对自由度的约束是线性的，而且可能是不均匀的。换句话说，它们总是具有
 * $x_{i_1} = \sum_{j=2}^M a_{i_j} x_{i_j} + b_i$
 * 的形式。处理存储和使用这些约束的deal.II类是AffineConstraints。
 *
 *  <h3>Eliminating constraints</h3>
 * 在建立全局系统矩阵和右手边时，可以不考虑约束条件，即简单地在单元上循环，将局部贡献加入全局矩阵和右手边对象。为了进行实际计算，你必须对线性系统进行
 * "浓缩"：消除受约束的自由度并将适当的值分配给无约束的自由度。这改变了有限元计算中使用的稀疏矩阵的稀疏模式，因此是一个相当昂贵的操作。事情的一般方案是，你建立你的系统，使用
 * AffineConstraints::condense()
 * 函数消除（浓缩）约束节点，然后你解决剩余的系统，最后你使用
 * AffineConstraints::distribute()
 * 函数从无约束节点的值中计算出约束节点的值。请注意，
 * AffineConstraints::condense()
 * 函数适用于线性系统的矩阵和右手边，而
 * AffineConstraints::distribute() 函数则适用于解向量。
 * 这种先建立线性系统，再消除约束自由度的方案效率很低，如果约束条件多，矩阵满，即特别是对3d和/或高阶或hp-finite元素，则是一个瓶颈。此外，它不可能在一个进程可能无法接触到矩阵元素的情况下实现%的并行计算。因此，我们提供了建立线性系统的第二种方法，使用下面讨论的
 * AffineConstraints::add_entries_local_to_global() 和
 * AffineConstraints::distribute_local_to_global()
 * 函数。得到的线性系统与调用 AffineConstraints::condense()
 * 函数后得到的线性系统是等价的。
 *
 *
 * @note
 * 这两种应用约束的方式都是将矩阵对角线的值设置为与矩阵中其他项相同大小的<i>positive</i>项。因此，你需要设置你的问题，使描述主要矩阵贡献的弱形式不是<i>negative
 * definite</i>。否则，像CG这样的迭代求解器会崩溃，或者像GMRES那样慢得多。
 *
 *
 * @note
 * 虽然这两种方式是<i>equivalent</i>，即通过任何一种方式计算的线性系统的解是相同的，但线性系统本身不一定具有相同的矩阵和右侧向量条目。具体来说，由于我们计算的方式不同，对应于受限自由度的矩阵对角线和右手边条目可能不同；但是，它们总是以这样的方式选择，即线性系统的解是相同的。
 * <h4>Condensing matrices and sparsity patterns</h4>
 * 如上所述，使用约束条件的第一种方式是在不考虑约束条件的情况下建立线性系统，然后将其
 * "浓缩 "掉。浓缩一个矩阵分四个步骤进行。
 *
 *
 * - 首先是建立稀疏模式（例如，使用 DoFTools::make_sparsity_pattern()); ）。
 *
 * - 那么浓缩矩阵的稀疏模式是由原始稀疏模式和约束条件组成的。
 *
 *
 * - 第三，全局矩阵的组装。
 *
 *
 * - 第四，矩阵最终被浓缩。
 * 在浓缩过程中，我们实际上没有改变稀疏模式、矩阵和向量的行数或列数。相反，凝结函数将非零条目添加到矩阵的稀疏模式中（其中有受限节点），矩阵的凝结过程将产生额外的非零元素。在浓缩过程本身中，受约束的行和列被分配到无约束节点的行和列中。受约束的自由度保持原位。为了不干扰求解过程，这些行和列用零和主对角线上的一个适当的正值来填充（我们选择其他对角线元素大小的平均值，以确保新的对角线条目具有与其他条目相同的大小顺序；这保留了矩阵的缩放特性）。右手边的相应数值被设置为零。这样一来，受约束的节点在方程组求解时将始终得到零值，并且不会再与其他节点耦合。
 * 与创建一个新的、更小的矩阵相比，保留矩阵中的条目有一个好处，即只需要一个矩阵和稀疏模式，因此需要的内存更少。此外，浓缩过程的成本较低，因为不是所有的而是只有矩阵中的受限值必须被复制。另一方面，求解过程将花费更长的时间，因为矩阵向量的乘法会在受约束的行中产生乘以零的结果。此外，矢量的大小更大，对于那些使用较大数量的辅助矢量的迭代求解方法（例如使用显式正交程序的方法）来说，会导致更多的内存消耗。尽管如此，这个过程由于其较低的内存消耗而更有效率。
 * 浓缩函数存在于不同的参数类型中。SparsityPattern,
 * SparseMatrix 和 BlockSparseMatrix。请注意，对于
 * PETScWrappers::SparseMatrix()
 * 类型的参数或其他PETSc或Trilinos矩阵封装类，没有任何版本。这是因为相对来说，要得到PETSc矩阵的稀疏结构的表示，并有效地修改它们是很困难的；这一点尤其适用于矩阵实际分布在一个计算机集群中的情况。如果你想使用PETSc/Trilinos矩阵，你可以复制一个已经浓缩的deal.II矩阵，或者以已经浓缩的形式组装PETSc/Trilinos矩阵，见下面的讨论。
 *
 *  <h4>Condensing vectors</h4>
 * 浓缩向量的工作原理与上面描述的矩阵的工作原理完全相同。请注意，缩合是一个等价的操作，也就是说，对一个向量或矩阵做一次以上的缩合操作与只做一次的结果相同：一旦一个对象被缩合，进一步的缩合操作就不会再改变它了。
 * 与矩阵凝结函数相反，矢量凝结函数存在于PETSc和Trilinos矢量的变体中。然而，使用它们通常很昂贵，应该避免。你应该使用与上述相同的技术来避免使用它们。
 *
 *  <h4>Avoiding explicit condensation</h4>
 * 有时，人们希望在一个线性系统建立之后，根本就避免对它进行显式凝结。想这样做有两个主要原因。
 * <ul>   <li>  缩合是一个昂贵的操作，特别是当有许多约束条件和/或矩阵有许多非零项时。对于三维或高多项式程度的计算，以及hp-finite element方法来说，这两种情况都很典型，例如见 @ref hp_paper  "hp-paper"。这是hp教程程序中讨论的情况， @ref  step_27 "  step-27  "，以及 step-22  和  @ref step_31  "  step-31  "。
 * <li>  你使用的矩阵可能没有 AffineConstraints::condense() 函数（例如，PETSc和Trilinos封装类就是这种情况，我们无法访问矩阵的底层表示，因此无法有效地实现 AffineConstraints::condense() 操作）。这种情况在  step-17  、  step-18  、  step-31  和  step-32  中讨论。   </ul>
 * 在这种情况下，一种可能性是在将局部条目转移到全局矩阵和向量的时刻就将其分配到最终目的地，同样在最初设置的时候就在浓缩的形式中建立一个稀疏的模式。
 * AffineConstraints类也为这些操作提供了支持。例如，
 * AffineConstraints::add_entries_local_to_global()
 * 函数将非零条目添加到一个稀疏模式对象中。它不仅添加了一个给定的条目，而且还添加了所有的条目，如果当前的条目对应于以后要消除的受限自由度，我们就必须写到这些条目。类似地，在将局部贡献复制到全局矩阵或向量时，可以使用
 * AffineConstraints::distribute_local_to_global()
 * 函数直接分配向量和矩阵中的条目。这些调用使得后续调用
 * AffineConstraints::condense()
 * 变得没有必要。关于它们的使用例子，请看上面提到的教程程序。
 * 注意，尽管它们的名字描述了函数的真正作用，
 * AffineConstraints::distribute_local_to_global()
 * 函数必须应用于矩阵和右手边的向量，而下面讨论的
 * AffineConstraints::distribute()
 * 函数则应用于求解线性系统后的向量。
 *
 *  <h3>Distributing constraints</h3>
 * 在求解浓缩方程组后，解向量必须被 "分配"：通过调用
 * AffineConstraints::condense()
 * 对原始线性系统的修改，导致一个线性系统对所有无约束的自由度都能正确求解，但对有约束的自由度的值却没有定义。为了得到这些自由度的正确值，你需要将无约束的值也
 * "分配 "给它们的有约束的同事。这是由
 * AffineConstraints::distribute()
 * 函数完成的。分布的操作在某种意义上撤销了凝结过程，但应该注意的是，它不是逆向操作。基本上，分布将受约束的节点的值设置为从约束中计算出来的值，给定的是无约束的节点的值加上可能的不均匀性。
 *
 *  <h3>Treatment of inhomogeneous constraints</h3>
 * 如果一些约束线有不均匀性（如果约束来自于不均匀边界条件的实现，这就是典型的情况），情况就会比仅仅由于悬挂节点的约束更复杂一些。这是因为消除矩阵中的非对角线值会在向量中消除的行中产生贡献。这意味着，不均匀性只能用同时作用于矩阵和向量的函数来处理。这意味着，如果在没有任何矩阵的情况下调用相应的凝结函数（或者如果矩阵之前已经被凝结过），所有的不均匀性都会被忽略。
 * 使用AffineConstraints类来实现Dirichlet边界条件在 step-22
 * 教程程序中讨论。另一个利用AffineConstraints的例子是
 * step-41
 * 。这里的情况要复杂一些，因为我们有一些不在边界上的约束。在创建AffineConstraints对象后，有两种方法来应用不均匀约束。
 * 第一种方法。
 *
 *
 * - 将 AffineConstraints::distribute_local_to_global() 函数应用于系统矩阵和右侧，参数use_inhomogeneities_for_rhs = false（即默认）。
 *
 *
 * - 使用 AffineConstraints::set_zero() 函数将不均匀约束部分的解设为零（或者从等于零的解矢量开始）。
 *
 *
 * - 解决()线性系统
 *
 *
 * - 将 AffineConstraints::distribute() 应用于解决方案中
 * 第二种方法。
 *
 *
 * - 使用参数use_inhomogeneities_for_rhs = true的 AffineConstraints::distribute_local_to_global() 函数，并将其应用于系统矩阵和右手方
 *
 *
 * - 将解的有关分量设置为不均匀的约束值（例如使用 AffineConstraints::distribute()) ）。
 *
 *
 * - 解决()线性系统
 *
 *
 * - 根据求解器现在你必须对解应用 AffineConstraints::distribute() 函数，因为求解器可以改变解中的约束值。对于一个基于Krylov的求解器来说，这应该不是严格意义上的需要，但是仍然有可能在不均匀值和解的值之间存在机器精度的差异，而且如果你有额外的约束，例如来自悬挂节点的约束，你可能无论如何都想调用 AffineConstraints::distribute() 。
 * 当然，这两种方法都导致了相同的最终答案，但方式不同。使用第一种方法（即在
 * <code>use_inhomogeneities_for_rhs = false</code> 中使用
 * AffineConstraints::distribute_local_to_global()),
 * 时，我们建立的线性系统在所有那些自由度受到约束的地方，右手边的条目都是零，在这些线的矩阵对角线上有一些正值。因此，线性系统的解向量对于不均匀约束的自由度会有一个零值，我们需要调用
 * AffineConstraints::distribute()
 * 来给这些自由度以正确的非零值。
 * 另一方面，在第二种方法中，对于不均匀约束自由度的矩阵对角线元素和相应的右手边条目，使线性系统的解已经具有正确的值（例如，如果约束条件是
 * $x_{13}=42$ ，那么如果矩阵除了对角线条目外是空的，则行
 * $13$ ，而 $b_{13}/A_{13,13}=42$ ，这样 $Ax=b$ 的解必须如愿满足
 * $x_{13}=42$ ）。因此，我们不需要在求解后调用
 * AffineConstraints::distribute()
 * 来修复解的不均匀约束成分，尽管这样做也无妨。
 * 还有一个问题，即采取哪种方法，以及为什么我们需要将第一种方法中的解向量的值设为零。这两个问题的答案都与迭代求解器解决线性系统的方式有关。为此，考虑到我们通常在残差下降到右手边法线的某个分数以下时停止迭代，或者说，在初始残差的某个分数以下时停止迭代。现在考虑这个问题。
 *
 *
 * - 在第一种方法中，受限自由度的右手边条目为零，也就是说，右手边的准则实际上只包括我们关心的那些部分。另一方面，如果我们从一个解向量开始，而这个解向量在受约束的条目中不为零，那么初始残差就会非常大，因为目前解向量中的值与线性系统的解（在这些部分中为零）不匹配。   因此，如果我们一旦将初始残差减少了某一系数就停止迭代，那么我们可能在一次迭代后就达到了阈值，因为受限自由度被迭代求解器在一次迭代中就解决了。如果初始残差是由这些自由度主导的，那么我们在第一步就看到了急剧的减少，尽管在这短短的一次迭代中我们在线性系统的其余部分并没有真正取得什么进展。我们可以通过以下方式来避免这个问题：一旦残差的规范达到<i>norm of the right hand side</i>的某个分数就停止迭代，或者我们可以将解的成分设置为零（从而减少初始残差），然后迭代直到达到<i>norm of the initial
 * residual</i>的某个分数。
 *
 *
 * - 在第二种方法中，如果迭代中的起始向量为零，我们会遇到同样的问题，因为此时残差可能被受限自由度所支配，其值与我们在解中希望的值不一致。我们可以通过调用 AffineConstraints::distribute() <i>before</i>求解线性系统（必要时在求解后再进行第二次），将解向量的相应元素设置为正确的值，从而再次规避这个问题。
 * 除了这些考虑，考虑我们有 $x_{3}=\tfrac 12 x_1 + \tfrac 12$
 * 这种不均匀约束的情况，例如，从 $x_{3}=\tfrac 12 (x_1 + x_2)$
 * 形式的悬挂节点约束，其中 $x_2$ 本身被边界值约束到
 * $x_2=1$  。在这种情况下，AffineConstraints容器当然不能找出
 * $x_3$
 * 的最终值，因此，不能正确设置解向量的第三分量。因此，第二种方法将不起作用，你应该采取第一种方法。
 *
 *  <h3>Dealing with conflicting constraints</h3>
*有些情况下，自由度受到不止一种方式的约束，有时是相互冲突的方式。例如，考虑下面的情况。       @image html conflicting_constraints.png ""  这里，蓝色标记的自由度 $x_0$ 是一个悬挂节点。如果我们使用三线有限元，即FE_Q(1)，那么它将带有约束条件  $x_0=\frac 12 (x_{1}+x_{2})$  。另一方面，它在边界上，如果我们施加了边界条件 $u|_{\partial\Omega}=g$ ，那么我们将有约束 $x_0=g_0$ ，其中 $g_0$ 是这个自由度位置上的边界函数 $g(\mathbf x)$ 的值。
 * 那么，哪一个会赢？或者说：哪一个<i>should</i>赢？这个问题没有好的答案。
 *
 *
 * - 如果悬挂的节点约束是最终执行的约束，那么对于一般的边界函数，所得到的解不再满足边界条件  $g$  。
 *
 *
 * - 如果反其道而行之，在这一点上，解决方案将不满足悬挂节点的约束，因此将不满足所选元素的规则性属性（例如，尽管使用 $Q_1$ 元素，但将不连续）。
 *
 *
 * -如果你考虑弯曲的边界，情况就会变得完全没有希望，因为那时边的中点（即悬挂的节点）一般不在母边上。因此，无论两个竞争约束的优先级如何，解决方案都不会是
 * $H^1$
 * 符合要求的。如果悬空节点约束获胜，那么解决方案将既不符合要求，也没有正确的边界值。换句话说，"正确
 * "的解决方案是什么并不完全清楚。在大多数情况下，这并不重要：无论是哪种情况，由不符合性或不正确的边界值引入的误差最差也会与离散化的整体误差处于同一等级。
 * 也就是说，如果你知道你想要的是什么，你应该怎么做。
 *
 *
 * - 如果你想让悬挂的节点约束获胜，那么首先通过 DoFTools::make_hanging_node_constraints() 函数建立这些约束。   然后用 VectorTools::interpolate_boundary_values() 将边界值插值到同一个AffineConstraints对象中。如果后一个函数遇到一个已经被约束的边界节点，它将简单地忽略这个节点的边界值，不触及约束。
 *
 * - 如果你想让边界值约束获胜，就像上面那样建立悬空节点约束，并使用这些约束用 AffineConstraints::distribute_local_to_global() 函数来组装矩阵（或者，另一种方法是组装矩阵，然后对其使用 AffineConstraints::condense() ）。在第二步，使用 VectorTools::interpolate_boundary_values() 函数返回 std::map ，并将其作为 MatrixTools::apply_boundary_values() 的输入，将边界节点设置为正确的值。
 * 两种行为也可以通过建立两个独立的AffineConstraints对象，并以特定的第二个参数调用
 * AffineConstraints::merge() 函数来实现。
 *
 *  <h3>Applying constraints indirectly with a LinearOperator</h3>
 * 有时候，直接压缩或消除线性方程组中的约束是不可取的，也是不可能的。特别是如果没有底层矩阵对象可以被压缩（或在装配过程中照顾到约束）。如果系统是由LinearOperator描述的，通常就是这种情况。
 * 在这种情况下，我们可以用修改后的系统@f[
 * (C^T A C + Id_c) \tilde x = C^T (b
 *
 * - A\,k)
 * @f]代替[1]（M. S. Shephard.作为直接刚度装配过程的一部分，通过变换应用的线性多点约束。<i>International Journal for Numerical Methods in Engineering</i> 20(11):2107-2112, 1985).
 * 这里， $A$ 是一个给定的（无约束的）系统矩阵，对于它，我们只假设可以应用于一个向量，但不一定可以访问单个矩阵条目。   $b$ 是线性方程组 $A\,x=b$ 的相应右手边。矩阵 $C$ 描述了存储在AffineConstraints对象中的线性约束的同质部分，向量 $k$ 是相应的不均匀性的向量。更确切地说，应用于向量 $x$ 的 AffineConstraints::distribute() 操作是@f[
 *  x \leftarrow C\,x+k.
 * @f]的操作。最后， $Id_c$ 表示约束自由度子空间上的身份。
 * 然后通过分配约束条件来恢复服从这些约束的 $A\,x=b$
 * 的相应解。   $x=C\tilde x+k$  .
 * 整个系统可以通过以下代码片段来设置和解决。
 * @code
 * #include <deal.II/lac/constrained_linear_operator.h>
 *
 * // ...
 *
 * // system_matrix
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * - unconstrained and assembled system matrix
 * // right_hand_side
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * - unconstrained and assembled right hand side
 * // affine_constraints
 *
 * - an AffineConstraints object
 * // solver
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * - an appropriate, iterative solver
 * // preconditioner
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * - a preconditioner
 *
 * const auto op_a = linear_operator(system_matrix);
 * const auto op_amod = constrained_linear_operator(affine_constraints, op_a);
 * Vector<double> rhs_mod = constrained_right_hand_side(affine_constraints,
 *                                                    op_a,
 *                                                    right_hand_side);
 *
 * solver.solve(op_amod, solution, rhs_mod, preconditioner);
 * affine_constraints.distribute(solution);
 * @endcode
 *
 */


