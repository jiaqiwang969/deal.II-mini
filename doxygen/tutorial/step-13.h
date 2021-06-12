/**
@page step_13 The step-13 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Backgroundandpurpose">Background and purpose</a>
        <li><a href="#Whattheprogramdoes">What the program does</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Evaluationofthesolution">Evaluation of the solution</a>
      <ul>
        <li><a href="#Pointevaluation">%Point evaluation</a>
        <li><a href="#Generatingoutput">Generating output</a>
        <li><a href="#Otherevaluations">Other evaluations</a>
      </ul>
        <li><a href="#TheLaplacesolverclasses">The Laplace solver classes</a>
      <ul>
        <li><a href="#Anabstractbaseclass">An abstract base class</a>
        <li><a href="#Ageneralsolverclass">A general solver class</a>
        <li><a href="#Aprimalsolver">A primal solver</a>
        <li><a href="#Globalrefinement">Global refinement</a>
        <li><a href="#LocalrefinementbytheKellyerrorindicator">Local refinement by the Kelly error indicator</a>
      </ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Thedriverroutines">The driver routines</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-13/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Backgroundandpurpose"></a><h3>Background and purpose</h3>



在这个例子程序中，我们将不太关注描述如何使用deal.II及其设施的新方法，而是介绍编写模块化和可扩展有限元程序的方法。其主要原因是现代研究软件的规模和复杂性：实现现代误差估计概念和自适应求解方法的应用程序往往变得相当大。例如，当这个程序是在2002年写的时候，deal.II的主要作者的三个最大的应用程序，在写这个例子程序的时候。<ol>  <li>  一个用非连续加尔金有限元法求解守恒双曲方程的程序。33,775行代码；  <li>  一个参数估计程序。28,980行代码；  <li>  一个波浪方程求解器：21,020行代码。   </ol> 

这个库本身--不包括示例程序和测试套件--在2002年春天有略多于15万行的代码。当然，现在已经大了好几倍了）。)这些应用程序的规模是一个人，甚至是一个有经验的程序员，所能处理的边缘。




上面的数字相当清楚地说明了一件事：没有被分解成较小的、大部分独立的片段的单体程序是没有办法生存的，因为即使是作者也会很快失去对程序不同部分之间各种依赖关系的概述。只有数据封装，例如使用面向对象的编程方法，以及通过定义小而固定的接口来实现模块化，才能帮助构建数据流和相互依赖。如果一个以上的人在开发一个程序，这也是一个绝对的先决条件，因为否则的话，混乱会很快出现，因为如果一个开发人员需要知道另一个人是否改变了不同模块的内部结构，如果它们没有被干净地分开的话。




在前面的例子中，你已经看到了库本身是如何被分成几个复合体的，每个复合体都建立在底层复合体之上，但相对独立于其他复合体。<ol>  <li>  三角形类复合物，以及相关的迭代器类； <li>  有限元类； <li>  DoFHandler类复合物，以及相关的迭代器，建立在三角形和有限元类之上； <li>  实现单元和实数单元间映射的类； <li>  FEValues类复合物，建立在有限元和映射之上。   </ol>  除了这些，还有大量的小类，当然还有以下的 "工具 "模块。<ol>  <li>  以各种图形格式输出；  <li>  线性代数类。   </ol>  这些复数也可以在deal.II手册网站的首页上以流程图的形式找到。




这个程序的目标是给出一个例子，说明一个相对简单的有限元程序的结构，使我们最终得到一组尽可能相互独立的模块。这使得我们可以在一端改变程序，而不必担心另一端的程序会被破坏，只要我们不触及两端交流的接口。当然，C++中的接口是抽象基类的声明。




在这里，我们将（再次）实现一个拉普拉斯求解器，尽管与之前的例子程序相比有一些不同。<ol>  <li>  实现数值求解方程过程的类不再负责驱动 "求解-估计误差-再求解 "的过程，而是将其委托给外部函数。这首先允许在更大的范围内将其作为一个构件，在这里拉普拉斯方程的解可能只是其中的一部分（例如，在非线性问题中，拉普拉斯方程可能要在每个非线性步骤中解决）。它还允许围绕该类建立一个框架，允许使用其他方程的求解器（但具有相同的外部接口），以备对不同类型的偏微分方程评估一些技术。   <li>  它将评估计算出的解的过程分割成一组单独的类。原因是，人们通常对偏微分方程的解本身不感兴趣，而是对它的某些方面感兴趣。例如，人们可能希望在弹性计算中计算某一边界的牵引力，或者在某一位置的接收器上计算地震波的信号。有时，人们可能对这些方面中的几个方面感兴趣。由于解的评估是通常不影响解的过程，我们把它拆成一个单独的模块，以便独立于解算器类的开发来开发这种评估过滤器。   <li>  将实现网格细化的类与计算解的类分开。   <li>  将我们要介绍的测试案例的描述与程序的其他部分分开。   <li>  使用WorkStream设施对线性系统的装配进行并行化。这是在 @ref threads "多处理器访问共享内存的并行计算 "文档模块中可以找到的广泛描述。该实现基本上遵循步骤9中已经描述过的内容。   </ol> 




该程序所做的事情并不新鲜。事实上，这更像是以前程序的混合体，从以前的例子中拆解了各种部分和功能。读者应该关注的是它们在这个程序中的安排方式，即程序中使用的软件设计技术，以达到实现所需数学方法的目的。然而，我们必须强调，软件设计在某种程度上也是一个主观的问题：不同的人有不同的编程背景，对编程的 "正确 "风格有不同的看法；因此，这个程序只表达了作者认为有用的做法，如果你对所选择的方式感到不舒服，不一定要采用这种风格来编写成功的数值软件。不过，它应该作为一个案例研究，用启发读者的想法来达到理想的目的。




一旦你完成了这个程序，你会注意到它的结构已经有些复杂了。然而，它只有大约850行代码，没有注释。在真正的应用程序中，当然会有注释和类文件，这将使其达到1200行。然而，与上面列出的应用程序相比，这仍然是很小的，因为它们的规模是它们的20到25倍。对于这么大的程序，从一开始就进行适当的设计是不可缺少的。否则，一旦它变得过于庞大而无法管理，就必须在其生命中的某一时刻重新设计它。




尽管如此，上面列出的三个程序都经历了重大的修改，甚至是重写。例如，波浪程序，在它还明显较小的时候，曾经被完全撕成了碎片，只是为了以更模块化的形式再次组装。那时，已经不可能在不影响代码的旧部分的情况下增加功能了（代码的主要问题是数据流：在时间依赖的应用中，主要的问题是什么时候把数据存储到磁盘，什么时候再重新加载；如果这不是以一种有组织的方式进行的，那么你最终会发现数据释放得太早，加载得太晚，或者根本没有释放）。尽管本例程序因此吸取了几年的经验，但它的设计肯定不是没有缺陷的，特别是可能不适合目标不同的应用。它应该作为一个灵感，让你以模块化的方式编写自己的应用程序，以避免过于紧密耦合的代码的陷阱。




<a name="Whattheprogramdoes"></a><h3>What the program does</h3>



程序实际做什么甚至不是这个程序的重点，程序的结构更重要。然而，用几句话来描述就是：求解给定右手边的拉普拉斯方程，使其解为函数  $u(x,t)=\exp(x+\sin(10y+5x^2))$  。计算的目标是得到解在点 $x_0=(0.5,0.5)$ 处的值，并比较我们在两种细化标准下解决这个值的准确性，即全局细化和通过Kelly等人的误差指标细化，我们已经在以前的例子中使用过。




像往常一样，这些结果将在本文件的相应部分进行讨论。在这样做的过程中，我们将发现一个关于两个细化标准的相对性能的略微令人恼火的观察。在以后的例子程序中，在这个程序的基础上，我们将设计一个不同的方法，希望它能比这里讨论的技术表现更好。




现在，所有的理论和传闻背景都说了这么多。了解一个项目的最好方法是看它，所以它就在这里。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 像所有的程序一样，我们从库中的include文件列表开始，像往常一样，它们的标准顺序是 <code>base</code>  --  <code>lac</code> -- <code>grid</code> -- <code>dofs</code>  --  <code>fe</code> -- <code>numerics</code> （因为每一类大致都是建立在前面的基础上），然后是C++标准头文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/table_handler.h> 
 * #include <deal.II/base/thread_management.h> 
 * #include <deal.II/base/work_stream.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 现在是C++标准头文件。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * #include <list> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step13 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Evaluationofthesolution"></a> 
 * <h3>Evaluation of the solution</h3>
 * 

 * 
 * 至于程序本身，我们首先定义了评估拉普拉斯方程解的类。事实上，它们可以评估每一种解，只要它是由一个 <code>DoFHandler</code> 对象和一个解向量描述的。我们首先在这里定义它们，甚至在实际生成要评估的解的类之前，因为我们需要声明一个抽象的基类，以便解算器类可以引用。
 * 

 * 
 * 从抽象的角度来看，我们声明一个纯粹的基类，它提供了一个评估算子（），它将对解进行评估（无论派生类如何考虑 <code>evaluation</code>  ）。由于这是该基类唯一真正的功能（除了一些簿记机器），我们通常把这样一个只有 <code>operator()</code> 的类称为C++术语中的 <code>functor</code> ，因为它的使用就像一个函数对象。
 * 

 * 
 * 这种函数类型的对象随后将被传递给求解器对象，后者将其应用于刚刚计算出的解决方案。然后，评估对象可以从解决方案中提取他们喜欢的任何数量。将这些评估函数放入一个单独的类的层次结构的好处是，在设计上它们不能使用求解器对象的内部结构，因此独立于求解器工作方式的变化。此外，在不修改求解器类的情况下编写另一个评价类是很容易的，这就加快了编程速度（不能使用另一个类的内部结构也意味着你不必担心它们--对评价器的编程通常是一个相当快的任务），以及编译速度（如果求解器和评价类被放在不同的文件中：求解器只需要看到抽象基类的声明，因此在增加一个新的评价类或修改一个旧的类时不需要被重新编译）。 与此相关的是，你可以在其他项目中重复使用这些评估类，解决不同的方程。
 * 

 * 
 * 为了提高代码在不同模块中的分离度，我们把评估类放到了一个自己的命名空间中。这使得在同一个程序中实际解决不同的方程更加容易，通过现有的构件进行组装。这样做的原因是，用于类似目的的类往往具有相同的名称，尽管它们是在不同的背景下开发的。为了能够在一个程序中一起使用它们，有必要将它们放在不同的命名空间中。我们在这里就是这样做的。
 * 

 * 
 * 
 * @code
 *   namespace Evaluation 
 *   { 
 * 
 * @endcode
 * 
 * 现在是评估类的抽象基类：它的主要目的是声明一个纯虚函数 <code>operator()</code> ，接收一个 <code>DoFHandler</code> 对象和解向量。为了能够只使用指向这个基类的指针，它还必须声明一个虚拟的析构器，但这个析构器什么也不做。除此之外，它只提供了一点簿记功能：由于我们通常想在后续的细化水平上评估解决方案，我们存储了当前细化周期的编号，并提供了一个函数来改变这个编号。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class EvaluationBase 
 *     { 
 *     public: 
 *       virtual ~EvaluationBase() = default; 
 * 
 *       void set_refinement_cycle(const unsigned int refinement_cycle); 
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler, 
 *                               const Vector<double> & solution) const = 0; 
 * 
 *     protected: 
 *       unsigned int refinement_cycle; 
 *     }; 
 * 
 *     template <int dim> 
 *     void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step) 
 *     { 
 *       refinement_cycle = step; 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Pointevaluation"></a> 
 * <h4>%Point evaluation</h4>
 * 

 * 
 * 下一件事是实现实际的评估类。正如介绍中指出的，我们想从解决方案中提取一个点值，所以第一个类在它的 <code>operator()</code> 中做这个。实际的点是通过构造函数给这个类的，还有一个表格对象，它将把它的发现放入其中。
 * 

 * 
 * 如果我们不能依靠知道实际使用的有限元，那么找出任意点的有限元域的值是相当困难的，因为这样我们就不能，例如，在节点之间进行插值。因此，为了简单起见，我们在这里假设我们要评估场的点实际上是一个节点。如果在求解的过程中，我们发现我们在所有顶点上循环时没有遇到这个点，那么我们就必须抛出一个异常，以便向调用的函数发出信号，说明出了问题，而不是默默地忽略这个错误。
 * 

 * 
 * 在  step-9  示例程序中，我们已经看到如何使用  <code>DeclExceptionN</code>  宏来声明这样一个异常类。我们在这里再次使用这种机制。
 * 

 * 
 * 由此可见，这个类的实际声明应该是很明显的。请注意，即使我们没有明确地列出一个析构器，编译器也会生成一个隐含的析构器，而且它和基类的析构器一样是虚拟的。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class PointValueEvaluation : public EvaluationBase<dim> 
 *     { 
 *     public: 
 *       PointValueEvaluation(const Point<dim> &evaluation_point, 
 *                            TableHandler &    results_table); 
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler, 
 *                               const Vector<double> & solution) const override; 
 * 
 *       DeclException1( 
 *         ExcEvaluationPointNotFound, 
 *         Point<dim>, 
 *         << "The evaluation point " << arg1 
 *         << " was not found among the vertices of the present grid."); 
 * 
 *  
 *       const Point<dim> evaluation_point; 
 *       TableHandler &   results_table; 
 *     }; 
 * 
 * @endcode
 * 
 * 至于定义，构造函数是微不足道的，只是接收数据并将其存储在对象本地的。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     PointValueEvaluation<dim>::PointValueEvaluation( 
 *       const Point<dim> &evaluation_point, 
 *       TableHandler &    results_table) 
 *       : evaluation_point(evaluation_point) 
 *       , results_table(results_table) 
 *     {} 
 * 
 * @endcode
 * 
 * 现在是本类中主要感兴趣的函数，即点值的计算。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void PointValueEvaluation<dim>:: 
 *          operator()(const DoFHandler<dim> &dof_handler, 
 *                const Vector<double> & solution) const 
 *     { 
 * 
 * @endcode
 * 
 * 首先分配一个变量，用来保存点值。用一个明显是假的值来初始化它，这样如果我们不能把它设置成一个合理的值，我们就会马上注意到。这在像本函数这样小的函数中可能没有必要，因为我们在这里可以很容易地看到所有可能的执行路径，但事实证明它对更复杂的情况是有帮助的，所以我们在这里也采用了这个策略。
 * 

 * 
 * 
 * @code
 *       double point_value = 1e20; 
 * @endcode
 * 
 * 然后
 * 循环所有单元格及其所有顶点，并检查顶点是否与评估点匹配。如果是这样，就提取点的值，设置一个标志，表示我们已经找到了感兴趣的点，然后退出循环。
 * 

 * 
 * 
 * @code
 *       bool evaluation_point_found = false; 
 *       for (const auto &cell : dof_handler.active_cell_iterators()) 
 *         if (!evaluation_point_found) 
 *           for (const auto vertex : cell->vertex_indices()) 
 *             if (cell->vertex(vertex) == evaluation_point) 
 *               { 
 * 
 * @endcode
 * 
 * 为了从全局解决方案矢量中提取点值，挑选属于感兴趣的顶点的那个分量，如果解决方案是矢量值的，则取其第一个分量。
 * 

 * 
 * 
 * @code
 *                 point_value = solution(cell->vertex_dof_index(vertex, 0)); 
 * 
 * @endcode
 * 
 * 请注意，我们在这里做了一个假设，这个假设并不总是有效的，如果这是实际应用的代码，而不是一个教程程序，就应该在类的声明中记录下来：我们假设用于我们试图评估的解决方案的有限元实际上有与顶点相关的自由度。例如，这对不连续元素来说是不成立的，因为形状函数的支持点恰好位于顶点，但不与顶点相关，而是与单元内部相关，因为与顶点相关意味着那里的连续性。这对于面向边缘的元素等也是不成立的。            理想情况下，我们会在函数开始时检查这一点，例如通过一个类似<code>Assert (dof_handler.get_fe().dofs_per_vertex  @>  0, ExcNotImplemented())</code>的语句，这应该可以在异常触发时很清楚地说明问题所在。在这种情况下，我们省略了它（这的确是不好的风格），但是知道这一点在这里并没有什么坏处，因为如果我们要求语句 <code>cell-@>vertex_dof_index(vertex,0)</code> 给我们顶点的DoF索引，如果没有的话，语句就会失败。            我们再次强调，这种对允许的有限元的限制应该在类的文档中说明。
 * 

 * 
 * 由于我们找到了正确的点，我们现在设置相应的标志并退出最里面的循环。由于设置了标志，外循环也将被终止。
 * 

 * 
 * 
 * @code
 *                 evaluation_point_found = true; 
 *                 break; 
 *               }; 
 * 
 * @endcode
 * 
 * 最后，我们要确定我们确实已经找到了评估点，因为如果不是这样，我们就不能在那里给出一个合理的解的值，反正剩下的计算也是无用的。所以通过 <code>AssertThrow</code> 程序中已经使用的 step-9 宏，确保我们确实找到了这个点。如果不是这样，这个宏就会抛出一个作为第二个参数给它的类型的异常，但与直接的 <code>throw</code> 语句相比，它在异常对象中填充了一组额外的信息，例如，产生异常的源文件和行号，以及失败的条件。如果你在你的主函数里有一个 <code>catch</code> 子句（就像这个程序一样），你会捕捉到所有没有在中间某个地方捕捉到的、因而已经处理过的异常，这些额外的信息会帮助你找出发生了什么以及哪里出了问题。
 * 

 * 
 * 
 * @code
 *       AssertThrow(evaluation_point_found, 
 *                   ExcEvaluationPointNotFound(evaluation_point)); 
 * 
 * @endcode
 * 
 * 注意，我们在其他示例程序中也使用了 <code>Assert</code> 宏。它与这里使用的 <code>AssertThrow</code> 宏不同的是，它只是中止程序，而不是抛出一个异常，而且它只在调试模式下这样做。它是用来检查作为参数传递给函数的向量大小的正确宏，以及类似的。
 * 

 * 
 * 然而，这里的情况是不同的：我们是否找到评估点可能会在不同的细化过程中发生变化（例如，如果点周围的四个单元被粗化掉了，那么在细化和粗化之后，点可能会消失）。这是在调试模式下无法预测的事情，但应该经常检查，在生产运行中也是如此。因此，这里使用了 <code>AssertThrow</code> 宏。
 * 

 * 
 * 现在，如果我们确信我们已经找到了评估点，我们可以把结果加入到结果表中。
 * 

 * 
 * 
 * @code
 *       results_table.add_value("DoFs", dof_handler.n_dofs()); 
 *       results_table.add_value("u(x_0)", point_value); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="Generatingoutput"></a> 
 * <h4>Generating output</h4>
 * 

 * 
 * 一种不同的，也许略显奇怪的 <code>evaluation</code> 的解决方案是将其以图形格式输出到一个文件中。因为在评估函数中，我们得到了一个 <code>DoFHandler</code> 对象和解决方案的向量，我们已经有了做这件事所需要的一切，所以我们可以在评估类中做这件事。实际上这样做而不是把它放到计算解决方案的类中的原因是，这样我们有更多的灵活性：如果我们选择只输出它的某些方面，或者根本不输出它。在任何情况下，我们都不需要修改求解器类，我们只需要修改其中的一个模块，就可以构建这个程序了。如上所述，这种形式的封装可以帮助我们保持程序的每个部分相当简单，因为接口保持简单，不可能访问隐藏的数据。
 * 

 * 
 * 由于这个生成输出的类是从普通的 <code>EvaluationBase</code> 基类派生出来的，它的主要接口是 <code>operator()</code> 函数。此外，它有一个构造函数，接收一个字符串，该字符串将被用作文件名的基本部分，输出将被发送到该文件名中（我们将用一个数字来增加它，表示细化周期的数量--基类手头有这个信息--以及一个后缀），构造函数还接收一个值，表示要求的格式，即我们将为哪个图形程序生成输出（然后我们也将从这个值中生成我们写入的文件名后缀）。
 * 

 * 
 * 关于输出格式，DataOutBase命名空间提供了一个枚举字段 DataOutBase::OutputFormat ，列出了所有支持的输出格式的名称。在编写本程序时，支持的图形格式由枚举值 <code>ucd</code> 、 <code>gnuplot</code>, <code>povray</code>, <code>eps</code> 、 <code>gmv</code>, <code>tecplot</code>, <code>tecplot_binary</code> 、 <code>dx</code>, <code>vtk</code> 等表示，但这个列表肯定会随着时间而增加。现在，在该基类的各种函数中，你可以使用这种类型的值来获得关于这些图形格式的信息（例如每种格式的文件所使用的默认后缀），你可以调用一个通用的 <code>write</code> 函数，然后根据给它的第二个参数的值表示所需的输出格式，将其分支到我们在以前的例子中已经使用的 <code>write_gnuplot</code>, <code>write_ucd</code> 等函数。这种机制使得编写一个可扩展的程序变得很简单，它可以在运行时决定使用哪种输出格式，同时也使得编写程序的方式变得相当简单，它可以利用新实现的输出格式，而不需要改变应用程序。
 * 

 * 
 * 在这两个字段中，即基本名称和输出格式描述符，构造函数取值并存储它们，以便以后由实际的评估函数使用。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class SolutionOutput : public EvaluationBase<dim> 
 *     { 
 *     public: 
 *       SolutionOutput(const std::string &             output_name_base, 
 *                      const DataOutBase::OutputFormat output_format); 
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler, 
 *                               const Vector<double> & solution) const override; 
 * 
 *     private: 
 *       const std::string               output_name_base; 
 *       const DataOutBase::OutputFormat output_format; 
 *     }; 
 * 
 *     template <int dim> 
 *     SolutionOutput<dim>::SolutionOutput( 
 *       const std::string &             output_name_base, 
 *       const DataOutBase::OutputFormat output_format) 
 *       : output_name_base(output_name_base) 
 *       , output_format(output_format) 
 *     {} 
 * 
 * @endcode
 * 
 * 按照上面的描述，生成实际输出的函数现在相对简单了。与以前的例子程序相比，唯一特别有趣的特征是使用了 DataOutBase::default_suffix 函数，返回给定格式文件的通常后缀（例如，".eps "用于封装的postscript文件，".gnuplot "用于Gnuplot文件），以及带有第二个参数的通用 DataOut::write() 函数，该函数根据作为第二个参数的格式描述符的值，在内部分支到不同图形格式的实际输出函数。
 * 

 * 
 * 还要注意，我们必须在 <code>this-@></code> 前加上前缀，以访问依赖模板的基类的成员变量。这里的原因，以及在程序中更进一步的原因，与 step-7 示例程序中描述的相同（在那里寻找 <code>two-stage name lookup</code> ）。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void SolutionOutput<dim>::operator()(const DoFHandler<dim> &dof_handler, 
 *                                          const Vector<double> & solution) const 
 *     { 
 *       DataOut<dim> data_out; 
 *       data_out.attach_dof_handler(dof_handler); 
 *       data_out.add_data_vector(solution, "solution"); 
 *       data_out.build_patches(); 
 * 
 *       std::ofstream out(output_name_base + "-" + 
 *                         std::to_string(this->refinement_cycle) + 
 *                         data_out.default_suffix(output_format)); 
 * 
 *       data_out.write(out, output_format); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="Otherevaluations"></a> 
 * <h4>Other evaluations</h4>
 * 

 * 
 * 在实际应用中，人们会在这里添加一个其他可能的评价类的列表，代表人们可能感兴趣的数量。对于这个例子，这些就足够了，所以我们关闭命名空间。
 * 

 * 
 * 
 * @code
 *   } // namespace Evaluation 
 * @endcode
 * 
 * 
 * <a name="TheLaplacesolverclasses"></a> 
 * <h3>The Laplace solver classes</h3>
 * 

 * 
 * 在定义了我们想知道的解决方案之后，我们现在应该关心如何去获得它。我们将把所有我们需要的东西都打包到一个自己的命名空间中，原因和上面的评估差不多。
 * 

 * 
 * 由于我们在前面的例子中已经相当详细地讨论了拉普拉斯求解器，所以下面就没有什么新东西了。相反，我们在很大程度上拆解了以前的例子，并以稍微不同的形式把它们放到这个例子程序中。因此，我们将主要讨论与以前的例子的不同之处。
 * 

 * 
 * 基本上，正如在介绍中已经说过的，这个例子中缺乏新的东西是故意的，因为它更多地是为了展示软件设计的实践，而不是数学。因此，下面解释的重点将更多地放在实际的实现上。
 * 

 * 
 * 
 * @code
 *   namespace LaplaceSolver 
 *   { 
 * @endcode
 * 
 * 
 * <a name="Anabstractbaseclass"></a> 
 * <h4>An abstract base class</h4>
 * 

 * 
 * 在定义拉普拉斯求解器时，我们首先声明一个抽象的基类，它本身没有任何功能，只是接受和存储一个指向三角形的指针，以便以后使用。
 * 

 * 
 * 这个基类是非常通用的，也可以用于任何其他静止问题。它提供了一些函数的声明，这些函数将在派生类中分别解决一个问题，用评估对象的列表对解决方案进行后处理，以及细化网格。在基类中，这些函数本身都没有做什么。
 * 

 * 
 * 由于缺乏实际功能，声明非常抽象的基类的编程风格类似于Smalltalk或Java程序中使用的风格，所有的类都是从完全抽象的类派生出来的 <code>Object</code>  ，甚至是数字表示。作者承认，他并不特别喜欢在C++中使用这种风格，因为它将风格置于理性之上。此外，它提倡对一切事物使用虚拟函数（例如，在Java中，所有的函数本身就是虚拟的），然而，这在许多应用中被证明是相当低效的，在这些应用中，函数往往只是访问数据，而不是进行计算，因此很快就会返回；这样，虚拟函数的开销就会很大。笔者的观点是，只要至少有一部分实际实现的代码可以被共享，从而被分离到基类中，就应该有抽象的基类。
 * 

 * 
 * 除了这些理论上的问题，我们在这里还有一个很好的理由，这个理由在下面会让读者更清楚。基本上，我们希望能够有一个不同的拉普拉斯求解器家族，这些求解器的差别很大，以至于无法找到更大的共同功能子集。因此，我们只是声明了这样一个抽象的基类，在构造函数中获取一个指向三角形的指针，并从此存储它。由于这个三角剖分将在所有的计算中使用，我们必须确保这个三角剖分在最后使用之前是有效的。我们通过保留一个 <code>SmartPointer</code> 到这个三角剖分来做到这一点，正如 step-7 中所解释的。
 * 

 * 
 * 注意，虽然指针本身被声明为常数（即在这个对象的整个生命周期中，指针指向同一个对象），但它没有被声明为指向一个常数三角的指针。事实上，通过这种方式，我们允许派生类在 <code>refine_grid</code> 函数中细化或粗化三角结构。
 * 

 * 
 * 最后，我们有一个函数 <code>n_dofs</code> 只是驱动函数的一个工具，用来决定我们是否要继续进行网格细化。它返回当前模拟的自由度数量。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class Base 
 *     { 
 *     public: 
 *       Base(Triangulation<dim> &coarse_grid); 
 *       virtual ~Base() = default; 
 * 
 *       virtual void solve_problem() = 0; 
 *       virtual void postprocess( 
 *         const Evaluation::EvaluationBase<dim> &postprocessor) const = 0; 
 *       virtual void         refine_grid()                            = 0; 
 *       virtual unsigned int n_dofs() const                           = 0; 
 * 
 *     protected: 
 *       const SmartPointer<Triangulation<dim>> triangulation; 
 *     }; 
 * 
 * @endcode
 * 
 * 仅有的两个非抽象函数的实现就相当无聊了。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     Base<dim>::Base(Triangulation<dim> &coarse_grid) 
 *       : triangulation(&coarse_grid) 
 *     {} 
 * @endcode
 * 
 * 
 * <a name="Ageneralsolverclass"></a> 
 * <h4>A general solver class</h4>
 * 

 * 
 * 下面是主类，它实现了组装线性系统的矩阵，解决它，并在解决方案上调用后处理器对象。它实现了基类中声明的  <code>solve_problem</code>  和  <code>postprocess</code>  函数。然而，它并没有实现 <code>refine_grid</code> 方法，因为网格细化将在一些派生类中实现。
 * 

 * 
 * 它还声明了一个新的抽象虚函数， <code>assemble_rhs</code>  ，需要在子类中重载。原因是我们将实现两个不同的类，它们将实现不同的方法来组装右手边的向量。这个函数在以下情况下可能也很有趣：右手边不仅仅取决于一个连续函数，还取决于其他东西，例如另一个离散问题的解，等等。后者经常发生在非线性问题中。
 * 

 * 
 * 正如我们之前提到的，这门课的实际内容并不是新的，而是以前的例子中已经使用过的各种技术的混合。因此，我们将不对它们进行详细讨论，而是让读者参考这些程序。
 * 

 * 
 * 基本上，用几句话来说，这个类的构造函数接收指向一个三角形、一个有限元和一个代表边界值的函数对象的指针。这些东西或者被传递给基类的构造函数，或者被存储起来并在以后用来生成一个 <code>DoFHandler</code> 对象。由于有限元和正交公式应该是匹配的，所以它也被传递给一个正交对象。
 * 

 * 
 * <code>solve_problem</code> 为实际求解设置数据结构，调用函数来组装线性系统，并求解它。
 * 

 * 
 * <code>postprocess</code> 函数最后接收一个评估对象并将其应用于计算出的解决方案。
 * 

 * 
 * <code>n_dofs</code> 函数最后实现了基类的纯虚拟函数。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class Solver : public virtual Base<dim> 
 *     { 
 *     public: 
 *       Solver(Triangulation<dim> &      triangulation, 
 *              const FiniteElement<dim> &fe, 
 *              const Quadrature<dim> &   quadrature, 
 *              const Function<dim> &     boundary_values); 
 *       virtual ~Solver() override; 
 * 
 *       virtual void solve_problem() override; 
 * 
 *       virtual void postprocess(
 *         const Evaluation::EvaluationBase<dim> &postprocessor) const override; 
 * 
 *       virtual unsigned int n_dofs() const override; 
 * 
 * @endcode
 * 
 * 在这个类的保护部分，我们首先有一些成员变量，其用途在前面的例子中应该很清楚。
 * 

 * 
 * 
 * @code
 *     protected: 
 *       const SmartPointer<const FiniteElement<dim>> fe; 
 *       const SmartPointer<const Quadrature<dim>>    quadrature; 
 *       DoFHandler<dim>                              dof_handler; 
 *       Vector<double>                               solution; 
 *       const SmartPointer<const Function<dim>>      boundary_values; 
 * 
 * @endcode
 * 
 * 然后我们声明一个抽象函数，该函数将用于组装右手边的内容。如上所述，在各种情况下，这个动作的必要性有很大的不同，所以我们将其推迟到派生类中。
 * 

 * 
 * 
 * @code
 *       virtual void assemble_rhs(Vector<double> &rhs) const = 0; 
 * 
 * @endcode
 * 
 * 接下来，在私有部分，我们有一个小类，它代表了整个线性系统，即一个矩阵、一个右手边和一个解向量，以及应用于它的约束，如那些由于悬挂节点而产生的约束。它的构造函数初始化了各种子对象，还有一个函数实现了共轭梯度法作为求解器。
 * 

 * 
 * 
 * @code
 *     private: 
 *       struct LinearSystem 
 *       { 
 *         LinearSystem(const DoFHandler<dim> &dof_handler); 
 * 
 *         void solve(Vector<double> &solution) const; 
 * 
 *         AffineConstraints<double> hanging_node_constraints; 
 *         SparsityPattern           sparsity_pattern; 
 *         SparseMatrix<double>      matrix; 
 *         Vector<double>            rhs; 
 *       }; 
 * 
 * @endcode
 * 
 * 最后，有一组函数将被用来组装实际的系统矩阵。这一组的主函数 <code>assemble_linear_system()</code> 使用以下两个辅助函数，在多核系统上并行计算矩阵。这样做的机制与  step-9  示例程序相同，并遵循  @ref threads  中概述的 WorkStream 概念。主函数还调用了组装右手边的虚拟函数。
 * 

 * 
 * 
 * @code
 *       struct AssemblyScratchData 
 *       { 
 *         AssemblyScratchData(const FiniteElement<dim> &fe, 
 *                             const Quadrature<dim> &   quadrature); 
 *         AssemblyScratchData(const AssemblyScratchData &scratch_data); 
 * 
 *         FEValues<dim> fe_values;
 *       };
 * 
 *  
 *  
 * 
 *       struct AssemblyCopyData 
 *       { 
 *         FullMatrix<double>                   cell_matrix; 
 *         std::vector<types::global_dof_index> local_dof_indices; 
 *       }; 
 * 
 *       void assemble_linear_system(LinearSystem &linear_system); 
 * 
 *       void local_assemble_matrix( 
 *         const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *         AssemblyScratchData &                                 scratch_data, 
 *         AssemblyCopyData &                                    copy_data) const; 
 * 
 *       void copy_local_to_global(const AssemblyCopyData &copy_data, 
 *                                 LinearSystem &          linear_system) const; 
 *     }; 
 * 
 * @endcode
 * 
 * 现在是该类的构造函数。它没有做什么，只是存储了给定对象的指针，并生成了 <code>DoFHandler</code> 对象，初始化了给定的三角形的指针。这使得DoF处理程序存储该指针，但并没有生成有限元编号（我们只在 <code>solve_problem</code> 函数中要求这样做）。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     Solver<dim>::Solver(Triangulation<dim> &      triangulation, 
 *                         const FiniteElement<dim> &fe, 
 *                         const Quadrature<dim> &   quadrature, 
 *                         const Function<dim> &     boundary_values) 
 *       : Base<dim>(triangulation) 
 *       , fe(&fe) 
 *       , quadrature(&quadrature) 
 *       , dof_handler(triangulation) 
 *       , boundary_values(&boundary_values) 
 *     {} 
 * 
 * @endcode
 * 
 * 解构器很简单，它只是清除存储在DoF处理程序对象中的信息以释放内存。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     Solver<dim>::~Solver() 
 *     { 
 *       dof_handler.clear(); 
 *     } 
 * 
 * @endcode
 * 
 * 下一个函数是解决这个问题的主要工作：它用给这个对象的构造函数的有限元来设置DoF处理程序对象，创建一个表示线性系统的对象（即矩阵、右手向量和解向量），调用函数来组装它，最后解决它。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void Solver<dim>::solve_problem() 
 *     { 
 *       dof_handler.distribute_dofs(*fe); 
 *       solution.reinit(dof_handler.n_dofs()); 
 * 
 *       LinearSystem linear_system(dof_handler); 
 *       assemble_linear_system(linear_system); 
 *       linear_system.solve(solution); 
 *     } 
 * 
 * @endcode
 * 
 * 如上所述， <code>postprocess</code> 函数接收一个评估对象，并将其应用于计算的解决方案。这个函数可以被多次调用，对用户要求的每一个解的评估都要调用一次。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void Solver<dim>::postprocess( 
 *       const Evaluation::EvaluationBase<dim> &postprocessor) const 
 *     { 
 *       postprocessor(dof_handler, solution); 
 *     } 
 * 
 * @endcode
 * 
 * <code>n_dofs</code> 函数应该是不言自明的。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     unsigned int Solver<dim>::n_dofs() const 
 *     { 
 *       return dof_handler.n_dofs(); 
 *     } 
 * 
 * @endcode
 * 
 * 下面的函数在每一步中组装矩阵和要解决的线性系统的右手边。我们将在几个层面上并行地做事情。首先，请注意，我们需要组装矩阵和右手边。这些都是独立的操作，我们应该并行地进行这些操作。为此，我们使用 @ref threads 文档模块中讨论的 "任务 "概念。本质上，我们想说的是 "这里有一些需要处理的事情，只要有CPU核可用就去做"，然后再做其他事情，当我们需要第一个操作的结果时，就等待它的完成。在第二层，我们想使用与我们在 step-9 中已经使用过的完全相同的策略来组装矩阵，即WorkStream概念。
 * 

 * 
 * 虽然我们可以考虑在做另一件事的时候在后台组装右侧或组装矩阵，但我们将选择前一种方法，只是因为调用 <code>Solver::assemble_rhs</code> 比调用 WorkStream::run 及其许多参数要简单得多。在任何情况下，代码看起来像这样，以组装整个线性系统。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void Solver<dim>::assemble_linear_system(LinearSystem &linear_system) 
 *     { 
 *       Threads::Task<void> rhs_task = 
 *         Threads::new_task(&Solver<dim>::assemble_rhs, *this, linear_system.rhs); 
 * 
 *       auto worker = 
 *         [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *                AssemblyScratchData &scratch_data, 
 *                AssemblyCopyData &   copy_data) { 
 *           this->local_assemble_matrix(cell, scratch_data, copy_data); 
 *         }; 
 * 
 *       auto copier = [this, &linear_system](const AssemblyCopyData &copy_data) { 
 *         this->copy_local_to_global(copy_data, linear_system); 
 *       }; 
 * 
 *       WorkStream::run(dof_handler.begin_active(), 
 *                       dof_handler.end(), 
 *                       worker, 
 *                       copier, 
 *                       AssemblyScratchData(*fe, *quadrature), 
 *                       AssemblyCopyData()); 
 *       linear_system.hanging_node_constraints.condense(linear_system.matrix); 
 * 
 * @endcode
 * 
 * 上面的语法需要一些解释。 WorkStream::run 有多个版本，期待不同的参数。在  step-9  中，我们使用了一个版本，它需要一对迭代器、一对指向具有非常具体的参数列表的成员函数的指针、一个指向这些成员函数必须工作的对象的指针或引用，以及一个抓取和复制数据对象。这有点限制性，因为这样调用的成员函数的参数列表必须与 WorkStream::run 所期望的完全一致：本地装配函数需要接收一个迭代器、一个抓取对象和一个复制对象；而复制-本地-全局函数需要接收的正是一个复制对象。但是，如果我们想要的东西稍微更通用一些呢？例如，在目前的程序中，copy-local-to-global函数需要知道将本地贡献写入哪个线性系统对象中，也就是说，它还必须接受一个 <code>LinearSystem</code> 参数。这在使用成员函数指针的方法中是行不通的。
 * 

 * 
 * 幸运的是，C++提供了一条出路。这些被称为函数对象。本质上， WorkStream::run 想要做的不是调用一个成员函数。它想调用一些函数，这些函数在第一种情况下需要一个迭代器、一个抓取对象和一个拷贝对象，而在第二种情况下需要一个拷贝对象。不管这些是成员函数、全局函数，还是其他什么，对WorkStream来说，真的不是很关心。因此，有第二个版本的函数只接收函数对象--具有  <code>operator()</code>  的对象，因此可以像函数一样被调用，不管它们真正代表什么。产生这种函数对象的典型方法是使用一个<a href="http:en.wikipedia.org/wiki/Anonymous_function">lambda function</a>，它用固定的值来包装函数调用，包括各个参数。所有属于外层函数签名的参数在lambda函数中被指定为常规的函数参数。固定值使用捕获列表（`[...]`）传递到lambda函数中。可以使用捕获默认值，也可以明确列出所有要绑定到lambda的变量。为了清楚起见，我们决定在这里省略捕获默认值，但是捕获列表同样可以是`[&]`，这意味着所有使用的变量都通过引用复制到lambda中。
 * 

 * 
 * 在这一点上，我们已经组装好了矩阵，并将其浓缩。右手边可能已经完全组装好了，也可能还没有，但是我们接下来想浓缩右手边的向量。我们只有在这个向量的组装完成后才能这样做，所以我们必须等待任务的完成；在计算机科学中，等待任务通常被称为 "加入 "任务，解释了我们下面调用的函数的名称。
 * 

 * 
 * 既然这个任务可能已经完成，也可能没有完成，既然我们可能要等它完成，我们不妨试着把其他需要完成的事情装进这个空隙。因此，我们首先插值边界值，然后再等待右手边的工作。当然，另一种可能性是在一个单独的任务中也插值边界值，因为这样做与我们到目前为止在这个函数中所做的其他事情无关。请自由地找到正确的语法，为这个插值创建一个任务，并在这个函数的顶部启动它，同时装配右手边。(你会发现这稍微有点复杂，因为 VectorTools::interpolate_boundary_values(), 有多个版本，所以简单地取地址 <code>&VectorTools::interpolate_boundary_values</code> 会产生一组重载函数，不能马上传递给 Threads::new_task() --你必须通过将地址表达式转换为函数指针类型，选择你想要的这个重载集合中的哪个元素，这是你想在任务中调用的特定版本的函数。)
 * 

 * 
 * 
 * @code
 *       std::map<types::global_dof_index, double> boundary_value_map; 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                0, 
 *                                                *boundary_values, 
 *                                                boundary_value_map); 
 * 
 *       rhs_task.join(); 
 *       linear_system.hanging_node_constraints.condense(linear_system.rhs); 
 * 
 * @endcode
 * 
 * 现在我们有了完整的线性系统，我们也可以处理边界值，需要从矩阵和右手边消除。
 * 

 * 
 * 
 * @code
 *       MatrixTools::apply_boundary_values(boundary_value_map, 
 *                                          linear_system.matrix, 
 *                                          solution, 
 *                                          linear_system.rhs); 
 *     } 
 * 
 * @endcode
 * 
 * 这组函数的后半部分是处理每个单元上的局部装配，并将局部贡献复制到全局矩阵对象中。这与  step-9  中描述的工作方式完全相同。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     Solver<dim>::AssemblyScratchData::AssemblyScratchData( 
 *       const FiniteElement<dim> &fe, 
 *       const Quadrature<dim> &   quadrature) 
 *       : fe_values(fe, quadrature, update_gradients | update_JxW_values) 
 *     {} 
 * 
 *     template <int dim> 
 *     Solver<dim>::AssemblyScratchData::AssemblyScratchData( 
 *       const AssemblyScratchData &scratch_data) 
 *       : fe_values(scratch_data.fe_values.get_fe(), 
 *                   scratch_data.fe_values.get_quadrature(), 
 *                   update_gradients | update_JxW_values) 
 *     {} 
 * 
 *     template <int dim> 
 *     void Solver<dim>::local_assemble_matrix( 
 *       const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *       AssemblyScratchData &                                 scratch_data, 
 *       AssemblyCopyData &                                    copy_data) const 
 *     { 
 *       const unsigned int dofs_per_cell = fe->n_dofs_per_cell(); 
 *       const unsigned int n_q_points    = quadrature->size(); 
 * 
 *       copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
 * 
 *       copy_data.local_dof_indices.resize(dofs_per_cell); 
 * 
 *       scratch_data.fe_values.reinit(cell); 
 * 
 *       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             copy_data.cell_matrix(i, j) += 
 *               (scratch_data.fe_values.shape_grad(i, q_point) * 
 *                scratch_data.fe_values.shape_grad(j, q_point) * 
 *                scratch_data.fe_values.JxW(q_point)); 
 * 
 *       cell->get_dof_indices(copy_data.local_dof_indices); 
 *     } 
 * 
 *     template <int dim> 
 *     void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data, 
 *                                            LinearSystem &linear_system) const 
 *     { 
 *       for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i) 
 *         for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j) 
 *           linear_system.matrix.add(copy_data.local_dof_indices[i], 
 *                                    copy_data.local_dof_indices[j], 
 *                                    copy_data.cell_matrix(i, j)); 
 *     } 
 * 
 * @endcode
 * 
 * 现在是实现线性系统类中动作的函数。首先，构造函数将所有数据元素初始化为正确的大小，并设置了一些额外的数据结构，例如由于悬挂节点而产生的约束。由于设置悬空节点和找出矩阵的非零元素是独立的，所以我们以并行方式进行（如果库被配置为使用并发，至少是这样；否则，这些动作是按顺序执行的）。注意，我们只启动一个线程，并在主线程中做第二个动作。由于只生成了一个任务，我们在这里不使用 <code>Threads::TaskGroup</code> 类，而是直接使用创建的一个任务对象来等待这个特定任务的退出。
 * 

 * 
 * 注意，占用 <code>DoFTools::make_hanging_node_constraints</code> 函数的地址有点麻烦，因为它实际上有三个，每个支持的空间维度都有一个。在C++中，获取重载函数的地址有些复杂，因为在这种情况下，操作符 <code>&</code> 返回的更像是一组值（所有具有该名称的函数的地址），然后选择正确的函数是下一步的工作。如果上下文决定采取哪一个（例如通过分配给一个已知类型的函数指针），那么编译器可以自己做，但如果这组指针应作为一个采取模板的函数的参数，编译器可以选择所有的，而不偏向于一个。因此，我们必须向编译器说明我们想要哪一个；为此，我们可以使用cast，但为了更清楚，我们把它分配给一个具有正确类型的临时 <code>mhnc_p</code> （简称<code>pointer to make_hanging_node_constraints</code>），并使用这个指针代替。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     Solver<dim>::LinearSystem::LinearSystem(const DoFHandler<dim> &dof_handler) 
 *     { 
 *       hanging_node_constraints.clear(); 
 * 
 *       void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) = 
 *         &DoFTools::make_hanging_node_constraints; 
 * 
 * @endcode
 * 
 * 启动一个辅助任务，然后在主线程上继续进行
 * 

 * 
 * 
 * @code
 *       Threads::Task<void> side_task = 
 *         Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints); 
 * 
 *       DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *       DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 * 
 * @endcode
 * 
 * 等到边上的任务完成后再继续前进
 * 

 * 
 * 
 * @code
 *       side_task.join(); 
 * 
 *       hanging_node_constraints.close(); 
 *       hanging_node_constraints.condense(dsp); 
 *       sparsity_pattern.copy_from(dsp); 
 * 
 * @endcode
 * 
 * 最后初始化矩阵和右手边的向量
 * 

 * 
 * 
 * @code
 *       matrix.reinit(sparsity_pattern); 
 *       rhs.reinit(dof_handler.n_dofs()); 
 *     } 
 * 
 * @endcode
 * 
 * 该类的第二个函数只是通过预处理的共轭梯度法来解决线性系统。这一点之前已经被广泛讨论过了，所以我们不再赘述了。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const 
 *     { 
 *       SolverControl            solver_control(1000, 1e-12); 
 *       SolverCG<Vector<double>> cg(solver_control); 
 * 
 *       PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *       preconditioner.initialize(matrix, 1.2); 
 * 
 *       cg.solve(matrix, solution, rhs, preconditioner); 
 * 
 *       hanging_node_constraints.distribute(solution); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="Aprimalsolver"></a> 
 * <h4>A primal solver</h4>
 * 

 * 
 * 在上一节中，我们实现了一个拉普拉斯求解器的基类，该基类缺乏组装右手边向量的功能，但是，由于其中的原因，我们已经解释了。现在我们实现了一个相应的类，它可以在问题的右边以函数对象的形式给出的情况下完成这一工作。
 * 

 * 
 * 这个类的动作和你在以前的例子中已经看到的差不多，所以简单解释一下就够了：构造函数和底层类的数据相同（它把所有的信息传递给底层类），除了一个表示问题右侧的函数对象。这个对象的指针被存储起来（同样作为一个 <code>SmartPointer</code> ，以确保这个函数对象只要还被这个类使用就不会被删除）。
 * 

 * 
 * 这个类的唯一功能部分是 <code>assemble_rhs</code> 方法，它的作用和它的名字一样。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class PrimalSolver : public Solver<dim> 
 *     { 
 *     public: 
 *       PrimalSolver(Triangulation<dim> &      triangulation, 
 *                    const FiniteElement<dim> &fe, 
 *                    const Quadrature<dim> &   quadrature, 
 *                    const Function<dim> &     rhs_function, 
 *                    const Function<dim> &     boundary_values); 
 * 
 *     protected: 
 *       const SmartPointer<const Function<dim>> rhs_function; 
 *       virtual void assemble_rhs(Vector<double> &rhs) const override; 
 *     }; 
 * 
 * @endcode
 * 
 * 这个类的构造函数基本上做了上面宣布的事情......
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &      triangulation, 
 *                                     const FiniteElement<dim> &fe, 
 *                                     const Quadrature<dim> &   quadrature, 
 *                                     const Function<dim> &     rhs_function, 
 *                                     const Function<dim> &     boundary_values) 
 *       : Base<dim>(triangulation) 
 *       , Solver<dim>(triangulation, fe, quadrature, boundary_values) 
 *       , rhs_function(&rhs_function) 
 *     {} 
 * 
 * @endcode
 * 
 * ... 和 <code>assemble_rhs</code> 函数一样。因为在前面的几个例子程序中已经解释过了，所以我们就不多说了。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const 
 *     { 
 *       FEValues<dim> fe_values(*this->fe, 
 *                               *this->quadrature, 
 *                               update_values | update_quadrature_points | 
 *                                 update_JxW_values); 
 * 
 *       const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell(); 
 *       const unsigned int n_q_points    = this->quadrature->size(); 
 * 
 *       Vector<double>                       cell_rhs(dofs_per_cell);
 *       std::vector<double>                  rhs_values(n_q_points); 
 *       std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *       for (const auto &cell : this->dof_handler.active_cell_iterators()) 
 *         { 
 *           cell_rhs = 0; 
 *           fe_values.reinit(cell); 
 *           rhs_function->value_list(fe_values.get_quadrature_points(), 
 *                                    rhs_values); 
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               cell_rhs(i) += fe_values.shape_value(i, q_point) * // 
 *                              rhs_values[q_point] *               // 
 *                              fe_values.JxW(q_point); 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             rhs(local_dof_indices[i]) += cell_rhs(i); 
 *         }; 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Globalrefinement"></a> 
 * <h4>Global refinement</h4>
 * 

 * 
 * 至此，除了 <code>refine_grid</code> 函数外，抽象基类的所有函数都已实现。现在我们将有两个类为 <code>PrimalSolver</code> 类实现这个函数，一个做全局细化，一个做局部细化的形式。
 * 

 * 
 * 第一个做全局细化的类相当简单：它的主函数只是调用 <code>triangulation-@>refine_global (1);</code> ，它做所有的工作。
 * 

 * 
 * 注意，由于 <code>Base</code> 类的基类是虚拟的，我们必须声明一个构造函数来初始化直接的基类和抽象的虚拟类。
 * 

 * 
 * 除了这个技术上的复杂性之外，这个类可能很简单，可以不做进一步的评论。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class RefinementGlobal : public PrimalSolver<dim> 
 *     { 
 *     public: 
 *       RefinementGlobal(Triangulation<dim> &      coarse_grid, 
 *                        const FiniteElement<dim> &fe, 
 *                        const Quadrature<dim> &   quadrature, 
 *                        const Function<dim> &     rhs_function, 
 *                        const Function<dim> &     boundary_values); 
 * 
 *       virtual void refine_grid() override; 
 *     }; 
 * 
 *     template <int dim> 
 *     RefinementGlobal<dim>::RefinementGlobal( 
 *       Triangulation<dim> &      coarse_grid, 
 *       const FiniteElement<dim> &fe, 
 *       const Quadrature<dim> &   quadrature, 
 *       const Function<dim> &     rhs_function, 
 *       const Function<dim> &     boundary_values) 
 *       : Base<dim>(coarse_grid) 
 *       , PrimalSolver<dim>(coarse_grid, 
 *                           fe, 
 *                           quadrature, 
 *                           rhs_function, 
 *                           boundary_values) 
 *     {} 
 * 
 *     template <int dim> 
 *     void RefinementGlobal<dim>::refine_grid() 
 *     { 
 *       this->triangulation->refine_global(1); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="LocalrefinementbytheKellyerrorindicator"></a> 
 * <h4>Local refinement by the Kelly error indicator</h4>
 * 

 * 
 * 第二个实现细化策略的类使用了之前各种示例程序中使用的凯利细化指标。由于这个指标已经在deal.II库中用自己的类实现了，所以这里没有太多的事情要做，只是调用计算指标的函数，然后用它来选择一些单元进行细化和粗化，并相应地对网格进行细化。
 * 

 * 
 * 同样，现在应该足够标准了，可以省去更多的注释。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class RefinementKelly : public PrimalSolver<dim> 
 *     { 
 *     public: 
 *       RefinementKelly(Triangulation<dim> &      coarse_grid, 
 *                       const FiniteElement<dim> &fe, 
 *                       const Quadrature<dim> &   quadrature, 
 *                       const Function<dim> &     rhs_function, 
 *                       const Function<dim> &     boundary_values); 
 * 
 *       virtual void refine_grid() override; 
 *     }; 
 * 
 *     template <int dim> 
 *     RefinementKelly<dim>::RefinementKelly(Triangulation<dim> &      coarse_grid, 
 *                                           const FiniteElement<dim> &fe, 
 *                                           const Quadrature<dim> &   quadrature, 
 *                                           const Function<dim> &rhs_function, 
 *                                           const Function<dim> &boundary_values) 
 *       : Base<dim>(coarse_grid) 
 *       , PrimalSolver<dim>(coarse_grid, 
 *                           fe, 
 *                           quadrature, 
 *                           rhs_function, 
 *                           boundary_values) 
 *     {} 
 * 
 *     template <int dim> 
 *     void RefinementKelly<dim>::refine_grid() 
 *     { 
 *       Vector<float> estimated_error_per_cell( 
 *         this->triangulation->n_active_cells()); 
 *       KellyErrorEstimator<dim>::estimate( 
 *         this->dof_handler, 
 *         QGauss<dim - 1>(this->fe->degree + 1), 
 *         std::map<types::boundary_id, const Function<dim> *>(), 
 *         this->solution, 
 *         estimated_error_per_cell); 
 *       GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation, 
 *                                                       estimated_error_per_cell, 
 *                                                       0.3, 
 *                                                       0.03); 
 *       this->triangulation->execute_coarsening_and_refinement(); 
 *     } 
 * 
 *   } // namespace LaplaceSolver 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 由于这又是一个学术性的例子，我们想对精确解和计算解进行相互比较。为此，我们需要声明代表精确解的函数类（用于比较和Dirichlet边界值），以及一个表示方程右边的类（这只是应用于我们想恢复的精确解的拉普拉斯算子）。
 * 

 * 
 * 在这个例子中，让我们选择函数 $u(x,y)=exp(x+sin(10y+5x^2))$ 作为精确解。在超过两个维度的情况下，只需用 <code>y</code> replaced by <code>z</code> 重复正弦系数，以此类推。鉴于此，以下两类可能是直接从以前的例子中得出的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Solution : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double Solution<dim>::value(const Point<dim> & p, 
 *                               const unsigned int component) const 
 *   { 
 *     (void)component; 
 *     AssertIndexRange(component, 1); 
 *     double q = p(0); 
 *     for (unsigned int i = 1; i < dim; ++i) 
 *       q += std::sin(10 * p(i) + 5 * p(0) * p(0)); 
 *     const double exponential = std::exp(q); 
 *     return exponential; 
 *   } 
 * 
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double RightHandSide<dim>::value(const Point<dim> & p, 
 *                                    const unsigned int component) const 
 *   { 
 *     (void)component; 
 *     AssertIndexRange(component, 1); 
 *     double q = p(0); 
 *     for (unsigned int i = 1; i < dim; ++i) 
 *       q += std::sin(10 * p(i) + 5 * p(0) * p(0)); 
 *     const double u  = std::exp(q); 
 *     double       t1 = 1, t2 = 0, t3 = 0; 
 *     for (unsigned int i = 1; i < dim; ++i) 
 *       { 
 *         t1 += std::cos(10 * p(i) + 5 * p(0) * p(0)) * 10 * p(0); 
 *         t2 += 10 * std::cos(10 * p(i) + 5 * p(0) * p(0)) - 
 *               100 * std::sin(10 * p(i) + 5 * p(0) * p(0)) * p(0) * p(0); 
 *         t3 += 100 * std::cos(10 * p(i) + 5 * p(0) * p(0)) * 
 *                 std::cos(10 * p(i) + 5 * p(0) * p(0)) - 
 *               100 * std::sin(10 * p(i) + 5 * p(0) * p(0)); 
 *       }; 
 *     t1 = t1 * t1; 
 * 
 *     return -u * (t1 + t2 + t3); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Thedriverroutines"></a> 
 * <h3>The driver routines</h3>
 * 

 * 
 * 现在缺少的只是实际选择各种选项的函数，以及在连续的更细的网格上运行模拟，监测网格细化的进展。
 * 

 * 
 * 我们在下面的函数中做到了这一点：它接收一个求解器对象和一个后处理（评估）对象的列表，并在间歇性的网格细化中运行它们。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void run_simulation( 
 *     LaplaceSolver::Base<dim> &                          solver, 
 *     const std::list<Evaluation::EvaluationBase<dim> *> &postprocessor_list) 
 *   { 
 * 
 * @endcode
 * 
 * 我们将给出一个我们目前正在计算的步骤的指示器，以便让用户知道一些事情仍在发生，并且程序没有处于无尽的循环中。这就是这个状态行的标题。
 * 

 * 
 * 
 * @code
 *     std::cout << "Refinement cycle: "; 
 * 
 * @endcode
 * 
 * 然后开始一个循环，只有当自由度数大于20000时才会结束（当然你可以改变这个限制，如果你需要更多--或者更少--你的程序的准确性）。
 * 

 * 
 * 
 * @code
 *     for (unsigned int step = 0; true; ++step) 
 *       { 
 * 
 * @endcode
 * 
 * 然后给这个迭代的 <code>alive</code> 指示。注意， <code>std::flush</code> 是需要的，以使文本真正出现在屏幕上，而不是只出现在某个缓冲区中，而这个缓冲区只有在我们下一次发出结束线时才会被刷新。
 * 

 * 
 * 
 * @code
 *         std::cout << step << " " << std::flush; 
 * 
 * @endcode
 * 
 * 现在在现在的网格上解决问题，并在其上运行评估器。迭代器进入列表的长类型名称有点烦人，但如果需要的话，可以用别名来缩短。
 * 

 * 
 * 
 * @code
 *         solver.solve_problem(); 
 * 
 *         for (const auto &postprocessor : postprocessor_list) 
 *           { 
 *             postprocessor->set_refinement_cycle(step); 
 *             solver.postprocess(*postprocessor); 
 *           }; 
 * 
 * @endcode
 * 
 * 现在检查是否需要更多的迭代，或者是否应该结束循环。
 * 

 * 
 * 
 * @code
 *         if (solver.n_dofs() < 20000) 
 *           solver.refine_grid(); 
 *         else 
 *           break; 
 *       }; 
 * 
 * @endcode
 * 
 * 最后结束我们显示状态报告的那一行。
 * 

 * 
 * 
 * @code
 *     std::cout << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 最后一个函数是接受一个求解器的名字（目前允许使用 "kelly "和 "global"），用一个粗网格（这里是无处不在的单位方格）和一个有限元对象（这里也是无处不在的双线性对象）创建一个求解器对象，并使用该求解器来要求在一连串的细化网格上解决问题。
 * 

 * 
 * 该函数还设置了两个评估函数，一个是在(0.5,0.5)点评估解决方案，另一个是将解决方案写入一个文件。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void solve_problem(const std::string &solver_name) 
 *   { 
 * 
 * @endcode
 * 
 * 第一个小任务：告诉用户将发生什么。因此，写一个标题行，并在下面写上与第一个标题相同长度的所有'-'字符的行。
 * 

 * 
 * 
 * @code
 *     const std::string header = 
 *       "Running tests with \"" + solver_name + "\" refinement criterion:"; 
 *     std::cout << header << std::endl 
 *               << std::string(header.size(), '-') << std::endl; 
 * 
 * @endcode
 * 
 * 然后设置三角法、有限元等。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim> triangulation; 
 *     GridGenerator::hyper_cube(triangulation, -1, 1); 
 *     triangulation.refine_global(2); 
 *     const FE_Q<dim>    fe(1); 
 *     const QGauss<dim>  quadrature(4); 
 *     RightHandSide<dim> rhs_function; 
 *     Solution<dim>      boundary_values; 
 * 
 * @endcode
 * 
 * 创建一个由该函数的参数指示的解算器对象。如果该名称不被识别，则抛出一个异常! 各自的求解器对象被存储在一个 `std::unique_ptr` 中，以避免使用后不得不删除指针。
 * 

 * 
 * 
 * @code
 *     std::unique_ptr<LaplaceSolver::Base<dim>> solver; 
 *     if (solver_name == "global") 
 *       solver = std::make_unique<LaplaceSolver::RefinementGlobal<dim>>( 
 *         triangulation, fe, quadrature, rhs_function, boundary_values); 
 *     else if (solver_name == "kelly") 
 *       solver = std::make_unique<LaplaceSolver::RefinementKelly<dim>>( 
 *         triangulation, fe, quadrature, rhs_function, boundary_values); 
 *     else 
 *       AssertThrow(false, ExcNotImplemented()); 
 * 
 * @endcode
 * 
 * 接下来创建一个表对象，其中将存储点（0.5,0.5）的数值解的值，并创建一个相应的评估对象。
 * 

 * 
 * 
 * @code
 *     TableHandler                          results_table; 
 *     Evaluation::PointValueEvaluation<dim> postprocessor1(Point<dim>(0.5, 0.5), 
 *                                                          results_table); 
 * 
 * @endcode
 * 
 * 还会生成一个评估器，将解决方案写出来。
 * 

 * 
 * 
 * @code
 *     Evaluation::SolutionOutput<dim> postprocessor2(std::string("solution-") + 
 *                                                      solver_name, 
 *                                                    DataOutBase::gnuplot); 
 * 
 * @endcode
 * 
 * 把这两个评价对象放在一个列表中...
 * 

 * 
 * 
 * @code
 *     std::list<Evaluation::EvaluationBase<dim> *> postprocessor_list; 
 *     postprocessor_list.push_back(&postprocessor1); 
 *     postprocessor_list.push_back(&postprocessor2); 
 * 
 * @endcode
 * 
 * 然后，我们可以将其传递给在连续细化的网格上实际运行模拟的函数。
 * 

 * 
 * 
 * @code
 *     run_simulation(*solver, postprocessor_list); 
 * 
 * @endcode
 * 
 * 当这一切完成后，写出点评估的结果。
 * 

 * 
 * 
 * @code
 *     results_table.write_text(std::cout); 
 * 
 * @endcode
 * 
 * 在所有结果之后再写上一行空白。
 * 

 * 
 * 
 * @code
 *     std::cout << std::endl; 
 *   } 
 * } // namespace Step13 
 * 
 * @endcode
 * 
 * 关于主函数没有什么可说的。它沿用了之前所有例子中的模式，试图捕捉被抛出的异常，并在我们得到一些信息时尽可能多地显示出来。剩下的就不言自明了。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       Step13::solve_problem<2>("global"); 
 *       Step13::solve_problem<2>("kelly"); 
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
 *     }; 
 * 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-13/doc/results.dox



<a name="Results"></a><h1>Results</h1>





这个程序的结果并不那么有趣--毕竟它的目的不是为了演示一些新的数学思想，也不是为了演示如何用deal.II编程，而是为了利用我们在前面的例子中所开发的材料，形成一些演示以模块化和可扩展的方式建立现代有限元软件的方法。




尽管如此，我们当然会展示程序的结果。最感兴趣的是点值计算，为此我们实现了相应的评估类。该程序的结果（即输出）看起来如下。

@code
Running tests with "global" refinement criterion:


-------------------------------------------------
Refinement cycle: 0 1 2 3 4 5 6
DoFs  u(x_0)
   25 1.2868
   81 1.6945
  289 1.4658
 1089 1.5679
 4225 1.5882
16641 1.5932
66049 1.5945


Running tests with "kelly" refinement criterion:


------------------------------------------------
Refinement cycle: 0 1 2 3 4 5 6 7 8 9 10 11
DoFs  u(x_0)
   25 1.2868
   47 0.8775
   89 1.5365
  165 1.2974
  316 1.6442
  589 1.5221
 1093 1.5724
 2042 1.5627
 3766 1.5916
 7124 1.5876
13111 1.5942
24838 1.5932
@endcode




这里令人惊讶的是，精确的数值是1.59491554...，而且计算该解显然出奇的复杂，甚至只达到百分之一的精度，尽管该解是平滑的（事实上是无限常可微）。这种平滑性显示在程序生成的图形输出中，这里是粗网格和凯利细化指标的前9个细化步骤。


 <table width="80%" align="center">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-0.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-1.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-3.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-4.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-5.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-6.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-7.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-8.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-9.png" alt="">
    </td>
  </tr>
</table> 


当我们已经在观看图片时，这是第八个网格，从顶部看。


 <img src="https://www.dealii.org/images/steps/developer/step-13.grid-kelly-8.png" alt=""> 


然而，我们还没有完成对点值计算的评估。事实上，将两个细化标准的误差 $e=|u(x_0)-u_h(x_0)|$ 绘制成图，可以得到下面的图片。


 <img src="https://www.dealii.org/images/steps/developer/step-13.error.png" alt=""> 





这幅图 <em> 令人不安的是，自适应网格细化不仅没有像人们通常期望的那样比全局细化好，甚至明显更差，因为它的收敛是不规则的，在使用后续网格的值时，阻止了所有的外推技术!另一方面，全局细化提供了一个完美的 $1/N$ 或 $h^{-2}$ 收敛历史，并提供了各种机会，甚至可以通过外推法来改善点值。因此，在这个例子中，全局网格细化必须被认为是优越的!这更令人惊讶，因为评估点不是在左边的某个地方，那里的网格是粗糙的，而是在右边，自适应细化也应该细化评估点周围的网格。




因此，我们以一个问题来结束对这个例子程序的讨论。

<p align="center"> <strong>  <em>  如果适应性不比全局细化好，那么它有什么问题？ </em>  </strong>





 <em>  在这个例子的结尾处进行练习。 </em>  有一个简单的原因导致适应性网格解决方案的不良和不规则行为。通过观察每个步骤中评估点周围的网格，可以很简单地找到这个原因--这个数据在程序的输出文件中。因此，一个练习是修改网格细化程序，使问题（一旦你注意到它）得以避免。第二个练习是检查结果是否比全局细化要好，如果是的话，是否能达到更好的收敛顺序（就自由度数而言），或者只达到一个更好的常数。




(  <em>  对于没有耐心的人来说，非常简短的回答。 </em> 在误差较大的步骤中，网格在评估点上是不规则的，即一些相邻的单元有悬空的节点；这破坏了一些超级近似的效果，而全局细化的网格可以从中受益。答案2：这个快速黑客

@code
  bool refinement_indicated = false;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (const auto v : cell->vertex_indices())
	  if (cell->vertex(v) == Point<dim>(.5,.5))
	    {
	      cell->clear_coarsen_flag();
	      refinement_indicated |= cell->refine_flag_set();
	    }
  if (refinement_indicated)
    for (const auto &cell : triangulation.active_cell_iterators())
      for (const auto v : cell->vertex_indices())
	    if (cell->vertex(v) == Point<dim>(.5,.5))
	      cell->set_refine_flag ();
@endcode

在执行细化之前，在Kelly细化类的细化函数中，将改善结果（练习：代码是怎么做的？不过，行为仍然是不规则的，所以不可能有关于收敛顺序的结果）。)


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-13.cc"
*/
