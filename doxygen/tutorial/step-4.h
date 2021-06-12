/**
@page step_4 The step-4 tutorial program
This tutorial depends on step-3.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeStep4codeclasstemplate">The <code>Step4</code> class template</a>
        <li><a href="#Righthandsideandboundaryvalues">Right hand side and boundary values</a>
        <li><a href="#ImplementationofthecodeStep4codeclass">Implementation of the <code>Step4</code> class</a>
      <ul>
        <li><a href="#Step4Step4">Step4::Step4</a>
        <li><a href="#Step4make_grid">Step4::make_grid</a>
        <li><a href="#Step4setup_system">Step4::setup_system</a>
        <li><a href="#Step4assemble_system">Step4::assemble_system</a>
        <li><a href="#Step4solve">Step4::solve</a>
        <li><a href="#Step4output_results">Step4::output_results</a>
        <li><a href="#Step4run">Step4::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-4/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{12,13} 

deal.II有一个独特的功能，我们称之为 "无维度编程"。你可能已经注意到，在前面的例子中，许多类的后缀都是角括号中的数字。这是为了表明，例如，二维和三维空间的三角形是不同的，但是相关的数据%类型。我们完全可以把它们称为 <code>Triangulation2d</code> and <code>Triangulation3d</code> 而不是 <code>Triangulation@<2@></code> 和 <code>Triangulation@<3@></code> 来命名这两个类，但这有一个重要的缺点：假设你有一个功能完全相同的函数，但在2D或3D三角形上，取决于我们目前想在哪个维度上解方程（如果你不相信一个函数在所有维度上都做同样的事情是常见的情况，看看下面的代码就知道了，2D和3D之间几乎没有区别！）。我们将不得不把同一个函数写两次，一次在 <code>Triangulation2d</code> 上工作，一次在 <code>Triangulation3d</code> 上工作。这在编程中是一个不必要的障碍，并且导致了保持两个函数同步的麻烦（最好是），或者在两个版本不同步时难以发现错误（最坏的情况是；这可能是更常见的情况）。





这种障碍可以通过使用C++语言提供的一些模板魔法来规避：模板化的类和函数并不是真正的类或函数，而只是取决于一个尚未定义的数据类型参数或在定义时也未知的数值的一种模式。然而，如果你向它提供了所需的信息，编译器可以从这些模板中建立适当的类或函数。当然，模板的部分内容可以依赖于模板参数，它们将在编译时被解析为特定的模板参数。例如，考虑下面这段代码。

@code
  template <int dim>
  void make_grid (Triangulation<dim> &triangulation)
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
  };
@endcode






在编译器看到这个函数的时候，它对 <code>dim</code> 的实际值并不了解。编译器唯一拥有的是一个模板，即蓝图，如果给定 <code>make_grid</code> 的特定值有一个未知的值，编译器暂时没有可以生成的代码。




然而，如果以后下来，编译器会遇到一些代码，例如，看起来像这样。

@code
  Triangulation<2> triangulation;
  make_grid (triangulation);
@endcode

那么编译器将推断出请求将函数 <code>make_grid</code> 替换为 <code>dim==2</code> ，并将上述模板编译为一个到处都用2替换了dim的函数，也就是说，它将编译该函数，就好像它被定义为

@code
  void make_grid (Triangulation<2> &triangulation)
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
  };
@endcode






然而，值得注意的是，函数 <code>GridGenerator::hyper_cube</code> 也取决于维度，所以在这种情况下，编译器将调用函数 <code>GridGenerator::hyper_cube@<2@></code> ，而如果dim是3，它将调用 <code>GridGenerator::hyper_cube@<3@></code> ，这可能是（实际上是）一个完全无关的函数。




对成员变量也可以这样做。考虑一下下面的函数，它可能反过来调用上面的函数。

@code
  template <int dim>
  void make_grid_and_dofs (Triangulation<dim> &triangulation)
  {
    make_grid (triangulation);


    DoFHandler<dim> dof_handler(triangulation);
    ...
  };
@endcode

这个函数有一个类型为  <code>DoFHandler@<dim@></code>  的成员变量。同样，编译器在知道哪个维度之前不能编译这个函数。如果你像上面那样为一个特定的维度调用这个函数，编译器将使用模板，用调用的维度替换所有出现的dim，并编译它。如果你为不同的维度多次调用该函数，它将多次编译，每次都调用正确的 <code>make_grid</code> 函数，并为成员变量保留适当的内存量；注意， <code>DoFHandler</code> 的大小可能，事实上也确实取决于空间维度。




deal.II库是围绕这个独立于维度的编程概念建立的，因此允许你以一种不需要区分空间维度的方式来编程。应该注意的是，只有在极少数的地方才有必要使用 <code>if</code>s or <code>switch</code> es来实际比较尺寸。然而，由于编译器必须为每个维度单独编译每个函数，即使在那里，它在编译时也知道 <code>dim</code> 的值，因此将能够优化掉 <code>if</code> 语句和未使用的分支。




在这个例子程序中，我们将展示如何独立编程维度（事实上，这比你必须照顾到维度还要简单），我们将把上一个例子的拉普拉斯问题扩展到一个同时在两个和三个空间维度运行的程序。其他的扩展是使用非恒定的右手边函数和非零边界值。




 @note  在使用模板时，C++强加了各种语法限制，有时让人有点难以理解为什么到底要这样写。一个典型的例子是，在很多地方都需要使用关键字 <code>typename</code> 。如果你已经不完全熟悉，那么在<a
href="http://www.dealii.org/">deal.II homepage</a>中链接的deal.II常见问题解答（FAQ）中解释了其中的几个困难。

<！--我们需要一个空行来正确结束上述块。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * 前面几个（很多）include文件已经在前面的例子中使用过了，所以我们在这里不再解释它们的含义。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 这是新的，但是：在前面的例子中，我们从线性求解器得到了一些不需要的输出。如果我们想抑制它，我们必须包括这个文件，并在程序的某个地方添加一行字（见下面的main()函数）。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/logstream.h> 
 * 
 * @endcode
 * 
 * 最后一步，和以前的程序一样，是将所有deal.II的类和函数名导入全局命名空间中。
 * 

 * 
 * 
 * @code
 * using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep4codeclasstemplate"></a> 
 * <h3>The <code>Step4</code> class template</h3>
 * 

 * 
 * 这又是前面例子中的 <code>Step4</code> 类。唯一不同的是，我们现在把它声明为一个带有模板参数的类，而模板参数当然是我们要解决拉普拉斯方程的空间维度。当然，几个成员变量也取决于这个维度，特别是Triangulation类，它必须分别表示四边形或六面体。除此以外，一切都和以前一样。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * class Step4 
 * { 
 * public: 
 *   Step4(); 
 *   void run(); 
 * 
 * private: 
 *   void make_grid(); 
 *   void setup_system(); 
 *   void assemble_system(); 
 *   void solve(); 
 *   void output_results() const; 
 * 
 *   Triangulation<dim> triangulation; 
 *   FE_Q<dim>          fe; 
 *   DoFHandler<dim>    dof_handler; 
 * 
 *   SparsityPattern      sparsity_pattern; 
 *   SparseMatrix<double> system_matrix; 
 * 
 *   Vector<double> solution; 
 *   Vector<double> system_rhs; 
 * }; 
 * @endcode
 * 
 * 
 * <a name="Righthandsideandboundaryvalues"></a> 
 * <h3>Right hand side and boundary values</h3>
 * 

 * 
 * 在下文中，我们又声明了两个类，表示右手边和非均质的Dirichlet边界值。两者都是一个二维空间变量的函数，所以我们也将它们声明为模板。
 * 

 * 
 * 这些类中的每一个都是从一个共同的、抽象的基类Function派生出来的，它声明了所有函数都必须遵循的共同接口。特别是，具体的类必须重载 <code>value</code> 函数，该函数接收二维空间中的一个点作为参数，并将该点的值作为 <code>double</code> 变量返回。
 * 

 * 
 * <code>value</code> 函数需要第二个参数，我们在这里将其命名为 <code>component</code>  : 这只适用于矢量值函数，你可能想访问点 <code>p</code> 处的矢量的某个分量。然而，我们的函数是标量的，所以我们不需要担心这个参数，在函数的实现中也不会使用它。在库的头文件中，Function基类对 <code>value</code> 函数的声明中，分量的默认值为0，所以我们在访问右侧的 <code>value</code> 函数时，只需要一个参数，即我们要评估函数的点。然后，对于标量函数，可以简单地省略分量的值。
 * 

 * 
 * 函数对象在库中很多地方都有使用（例如，在 step-3 中我们使用了一个 Functions::ZeroFunction 实例作为 VectorTools::interpolate_boundary_values) 的参数，这是我们定义一个继承自Function的新类的第一个教程。由于我们只调用 Function::value(), ，我们可以只用一个普通的函数（这就是 step-5 中的做法），但由于这是一个教程，为了举例说明，我们继承了Function。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * class RightHandSide : public Function<dim> 
 * { 
 * public: 
 *   virtual double value(const Point<dim> & p, 
 *                        const unsigned int component = 0) const override; 
 * }; 
 * 
 * template <int dim> 
 * class BoundaryValues : public Function<dim> 
 * { 
 * public: 
 *   virtual double value(const Point<dim> & p, 
 *                        const unsigned int component = 0) const override; 
 * }; 
 * 
 * @endcode
 * 
 * 如果你不熟悉上述函数声明中的关键字 "virtual "和 "override "是什么意思，你可能会想看看你最喜欢的C++书籍或在线教程，如http:www.cplusplus.com/doc/tutorial/polymorphism/ 。从本质上讲，这里发生的事情是Function<dim>是一个 "抽象 "基类，它声明了某种 "接口"--一组可以在这类对象上调用的函数。但它实际上并没有*实现*这些函数：它只是说 "Function对象是这样的"，但它实际上是什么样的函数，则留给实现了`value()`函数的派生类。
 * 

 * 
 * 从另一个类中派生出一个类，通常称为 "is-a "关系函数。在这里，`RightHandSide`类 "是一个 "函数类，因为它实现了Function基类所描述的接口。("value() "函数的实际实现在下面的代码块中)。那么`virtual`关键字意味着 "是的，这里的函数可以被派生类覆盖"，而`override`关键字意味着 "是的，这实际上是一个我们知道已经被声明为基类一部分的函数"。覆盖 "关键字不是严格必要的，但它是防止打字错误的一个保险。如果我们把函数的名字或一个参数的类型弄错了，编译器会警告我们说："你说这个函数覆盖了基类中的一个函数，但实际上我不知道有任何这样的函数有这个名字和这些参数。"
 * 

 * 
 * 但回到这里的具体案例。在本教程中，我们选择2D中的函数 $4(x^4+y^4)$ ，或者3D中的 $4(x^4+y^4+z^4)$ 作为右手边。我们可以用空间维度上的if语句来写这个区别，但这里有一个简单的方法，通过使用一个短循环，也允许我们在一维（或四维，如果你想这样做）中使用相同的函数。 幸运的是，编译器在编译时就知道循环的大小（记住，在你定义模板时，编译器不知道 <code>dim</code> 的值，但当它后来遇到语句或声明 <code>RightHandSide@<2@></code> 时，它将采取模板，用2替换所有出现的dim，并编译出结果函数）。 换句话说，在编译这个函数的时候，主体将被执行的次数是已知的，编译器可以将循环所需的开销降到最低；结果将和我们马上使用上面的公式一样快。
 * 

 * 
 * 最后要注意的是， <code>Point@<dim@></code> 表示二维空间中的一个点，它的各个组成部分（即 $x$ 、 $y$ 、...坐标）可以像C和C++中一样用（）运算符访问（事实上，[]运算符也同样有效），索引从0开始。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * double RightHandSide<dim>::value(const Point<dim> &p, 
 *                                  const unsigned int /*component*/) const 
 * { 
 *   double return_value = 0.0; 
 *   for (unsigned int i = 0; i < dim; ++i) 
 *     return_value += 4.0 * std::pow(p(i), 4.0); 
 * 
 *   return return_value; 
 * } 
 * 
 * @endcode
 * 
 * 作为边界值，我们选择二维的 $x^2+y^2$ ，三维的 $x^2+y^2+z^2$ 。这恰好等于从原点到我们想评估函数的点的矢量的平方，而不考虑维度。所以这就是我们的返回值。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * double BoundaryValues<dim>::value(const Point<dim> &p, 
 *                                   const unsigned int /*component*/) const 
 * { 
 *   return p.square(); 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeStep4codeclass"></a> 
 * <h3>Implementation of the <code>Step4</code> class</h3>
 * 

 * 
 * 接下来是利用上述函数的类模板的实现。和以前一样，我们将把所有东西写成模板，这些模板有一个形式参数 <code>dim</code> ，在我们定义模板函数时，我们假设这个参数是未知的。只有在以后，编译器才会发现 <code>Step4@<2@></code> (in the <code>main</code> 函数的声明，实际上），并在编译整个类时将 <code>dim</code> 替换成2，这个过程被称为 "模板的实例化"。这样做的时候，它也会用 <code>RightHandSide@<dim@></code> 的实例替换 <code>RightHandSide@<2@></code> ，并从类模板中实例化后一个类。
 * 

 * 
 * 事实上，编译器也会在 <code>main()</code> 中找到一个 <code>Step4@<3@></code> 声明。这将导致它再次回到一般的 <code>Step4@<dim@></code> 模板，替换所有出现的 <code>dim</code> ，这次是3，并第二次编译这个类。注意这两个实例  <code>Step4@<2@></code>  和  <code>Step4@<3@></code>  是完全独立的类；它们唯一的共同特征是它们都是从同一个通用模板中实例化出来的，但是它们不能相互转换，例如，它们没有共享代码（两个实例都是完全独立编译的）。
 * 

 * 
 * 
 * <a name="Step4Step4"></a> 
 * <h4>Step4::Step4</h4>
 * 

 * 
 * 在这个介绍之后，这里是  <code>Step4</code>  类的构造函数。它指定了所需的有限元素的多项式程度，并将DoFHandler与三角形关联起来，就像在前面的例子程序中一样，  step-3  。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * Step4<dim>::Step4() 
 *   : fe(1) 
 *   , dof_handler(triangulation) 
 * {} 
 * @endcode
 * 
 * 
 * <a name="Step4make_grid"></a> 
 * <h4>Step4::make_grid</h4>
 * 

 * 
 * 网格的创建在本质上是与维度有关的东西。然而，只要领域在二维或三维中足够相似，库就可以为你抽象。在我们的例子中，我们想再次在二维的正方形 $[-1,1]\times [-1,1]$ 上求解，或者在三维的立方体 $[-1,1] \times [-1,1] \times [-1,1]$ 上求解；两者都可以被称为 GridGenerator::hyper_cube(), ，因此我们可以在任何维度上使用同一个函数。当然，在二维和三维中创建超立方体的函数有很大的不同，但这是你不需要关心的事情。让库来处理这些困难的事情。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step4<dim>::make_grid() 
 * { 
 *   GridGenerator::hyper_cube(triangulation, -1, 1); 
 *   triangulation.refine_global(4); 
 * 
 *   std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
 *             << std::endl 
 *             << "   Total number of cells: " << triangulation.n_cells() 
 *             << std::endl; 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step4setup_system"></a> 
 * <h4>Step4::setup_system</h4>
 * 

 * 
 * 这个函数看起来和前面的例子完全一样，尽管它执行的动作在细节上有很大的不同，如果 <code>dim</code> 刚好是3。从用户的角度来看，唯一显著的区别是所产生的单元格数量，在三个空间维度中比两个空间维度中要高得多
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step4<dim>::setup_system() 
 * { 
 *   dof_handler.distribute_dofs(fe); 
 * 
 *   std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *             << std::endl; 
 * 
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *   sparsity_pattern.copy_from(dsp); 
 * 
 *   system_matrix.reinit(sparsity_pattern); 
 * 
 *   solution.reinit(dof_handler.n_dofs()); 
 *   system_rhs.reinit(dof_handler.n_dofs()); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step4assemble_system"></a> 
 * <h4>Step4::assemble_system</h4>
 * 

 * 
 * 与前面的例子不同，我们现在想使用一个非恒定的右侧函数和非零边界值。这两个任务都是很容易实现的，只需在矩阵和右手边的组合中增加几行代码即可。
 * 

 * 
 * 更有趣的是，我们将矩阵和右手边的向量维度独立组装起来的方式：与二维的情况根本没有区别。由于这个函数中使用的重要对象（正交公式、FEValues）也通过模板参数的方式依赖于维度，它们可以为这个函数所编译的维度正确设置一切。通过使用模板参数声明所有可能依赖于维度的类，库可以为你完成几乎所有的工作，你不需要关心大多数事情。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step4<dim>::assemble_system() 
 * { 
 *   QGauss<dim> quadrature_formula(fe.degree + 1); 
 * 
 * @endcode
 * 
 * 我们希望有一个非恒定的右手，所以我们使用上面声明的类的一个对象来生成必要的数据。由于这个右侧对象只在本函数中局部使用，所以我们在这里把它声明为一个局部变量。
 * 

 * 
 * 
 * @code
 *   RightHandSide<dim> right_hand_side; 
 * 
 * @endcode
 * 
 * 与之前的例子相比，为了评估非恒定右手函数，我们现在还需要我们目前所在单元上的正交点（之前，我们只需要FEValues对象中的形状函数的值和梯度，以及正交权重， FEValues::JxW() ）。我们可以通过给FEValues对象添加#update_quadrature_points标志来让它为我们做事。
 * 

 * 
 * 
 * @code
 *   FEValues<dim> fe_values(fe, 
 *                           quadrature_formula, 
 *                           update_values | update_gradients | 
 *                             update_quadrature_points | update_JxW_values); 
 * 
 * @endcode
 * 
 * 然后我们再次定义与前面程序中相同的缩写。这个变量的值当然取决于我们现在使用的维度，但是FiniteElement类为你做了所有必要的工作，你不需要关心与维度有关的部分。
 * 

 * 
 * 
 * @code
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *   Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 接下来，我们又要在所有的单元格上进行循环，并汇集局部贡献。 请注意，一个单元在两个空间维度上是一个四边形，但在三维上是一个六面体。事实上， <code>active_cell_iterator</code> 的数据类型是不同的，这取决于我们所处的维度，但对外界来说，它们看起来是一样的，你可能永远不会看到区别。在任何情况下，真正的类型是通过使用`auto`来隐藏的。
 * 

 * 
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators()) 
 *     { 
 *       fe_values.reinit(cell); 
 *       cell_matrix = 0; 
 *       cell_rhs    = 0; 
 * 
 * @endcode
 * 
 * 现在我们要把本地矩阵和右手边组合起来。这个过程和前面的例子完全一样，但是现在我们重新调整循环的顺序（我们可以安全地这样做，因为它们是相互独立的），并尽可能地合并本地矩阵和本地向量的循环，使事情变得更快。
 * 

 * 
 * 组装右手边与我们在 step-3 中的做法有唯一的区别：我们没有使用值为1的常数右手边，而是使用代表右手边的对象并在正交点对其进行评估。
 * 

 * 
 * 
 * @code
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
 *         for (const unsigned int i : fe_values.dof_indices()) 
 *           { 
 *             for (const unsigned int j : fe_values.dof_indices()) 
 *               cell_matrix(i, j) += 
 *                 (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
 *                  fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
 *                  fe_values.JxW(q_index));           // dx 
 * 
 *             const auto &x_q = fe_values.quadrature_point(q_index); 
 *             cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
 *                             right_hand_side.value(x_q) *        // f(x_q) 
 *                             fe_values.JxW(q_index));            // dx 
 *           } 
 * 
 * @endcode
 * 
 * 作为对这些循环的最后说明：当我们将局部贡献集合到 <code>cell_matrix(i,j)</code> 时，我们必须将形状函数 $i$ 和 $j$ 在点号q_index的梯度相乘并与标量权重JxW相乘。这就是实际发生的情况。  <code>fe_values.shape_grad(i,q_index)</code> 返回一个 <code>dim</code> 维向量，由 <code>Tensor@<1,dim@></code> 对象表示，将其与 <code>fe_values.shape_grad(j,q_index)</code> 的结果相乘的运算器*确保两个向量的 <code>dim</code> 分量被适当收缩，结果是一个标量浮点数，然后与权重相乘。在内部，这个操作符*确保对向量的所有 <code>dim</code> 分量都能正确发生，无论 <code>dim</code> 是2、3还是其他空间维度；从用户的角度来看，这并不值得费心，然而，如果想独立编写代码维度，事情就会简单很多。
 * 

 * 
 * 随着本地系统的组装，转移到全局矩阵和右手边的工作与之前完全一样，但在这里我们再次合并了一些循环以提高效率。
 * 

 * 
 * 
 * @code
 *       cell->get_dof_indices(local_dof_indices); 
 *       for (const unsigned int i : fe_values.dof_indices()) 
 *         { 
 *           for (const unsigned int j : fe_values.dof_indices()) 
 *             system_matrix.add(local_dof_indices[i], 
 *                               local_dof_indices[j], 
 *                               cell_matrix(i, j)); 
 * 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i); 
 *         } 
 *     } 
 * 
 * @endcode
 * 
 * 作为这个函数的最后一步，我们希望在这个例子中拥有非均质的边界值，与之前的例子不同。这是一个简单的任务，我们只需要用一个描述我们想使用的边界值的类的对象（即上面声明的 <code>BoundaryValues</code> 类）来替换那里使用的 Functions::ZeroFunction 。
 * 

 * 
 * 函数 VectorTools::interpolate_boundary_values() 只对标有边界指标0的面起作用（因为我们在下面的第二个参数中说该函数应该对其起作用）。如果有的面的边界指标不是0，那么函数interpolate_boundary_values将对这些面不起作用。对于拉普拉斯方程来说，什么都不做相当于假设在边界的这些部分，零诺伊曼边界条件成立。
 * 

 * 
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values; 
 *   VectorTools::interpolate_boundary_values(dof_handler, 
 *                                            0, 
 *                                            BoundaryValues<dim>(), 
 *                                            boundary_values); 
 *   MatrixTools::apply_boundary_values(boundary_values, 
 *                                      system_matrix, 
 *                                      solution, 
 *                                      system_rhs); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step4solve"></a> 
 * <h4>Step4::solve</h4>
 * 

 * 
 * 解决线性方程组是在大多数程序中看起来几乎相同的事情。特别是，它与维度无关，所以这个函数是从前面的例子中逐字复制的。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step4<dim>::solve() 
 * { 
 *   SolverControl            solver_control(1000, 1e-12); 
 *   SolverCG<Vector<double>> solver(solver_control); 
 *   solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 
 * 
 * @endcode
 * 
 * 不过我们做了一个补充：由于我们抑制了线性求解器的输出，我们必须手工打印迭代次数。
 * 

 * 
 * 
 * @code
 *   std::cout << "   " << solver_control.last_step() 
 *             << " CG iterations needed to obtain convergence." << std::endl; 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step4output_results"></a> 
 * <h4>Step4::output_results</h4>
 * 

 * 
 * 这个函数也做了  step-3  中各自的工作。这里也没有改变维度的独立性。
 * 

 * 
 * 由于程序将同时运行拉普拉斯求解器的2D和3D版本，我们使用文件名中的维度为每次运行生成不同的文件名（在一个更好的程序中，我们将检查 <code>dim</code> 是否可以有2或3以外的其他值，但为了简洁起见，我们在这里忽略了这一点）。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step4<dim>::output_results() const 
 * { 
 *   DataOut<dim> data_out; 
 * 
 *   data_out.attach_dof_handler(dof_handler); 
 *   data_out.add_data_vector(solution, "solution"); 
 * 
 *   data_out.build_patches(); 
 * 
 *   std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk"); 
 *   data_out.write_vtk(output); 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4run"></a> 
 * <h4>Step4::run</h4>
 * 

 * 
 * 这是一个对所有事情都有最高级别控制的函数。除了一行额外的输出外，它与前面的例子相同。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step4<dim>::run() 
 * { 
 *   std::cout << "Solving problem in " << dim << " space dimensions." 
 *             << std::endl; 
 * 
 *   make_grid(); 
 *   setup_system(); 
 *   assemble_system(); 
 *   solve(); 
 *   output_results(); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 这是主函数。它看起来也大多像 step-3 中的内容，但如果你看下面的代码，注意我们是如何首先创建一个 <code>Step4@<2@></code> 类型的变量（迫使编译器用 <code>dim</code> replaced by <code>2</code> 编译类模板）并运行一个2d模拟，然后我们用3d做整个事情。
 * 

 * 
 * 在实践中，这可能不是你经常做的事情（你可能要么想解决一个2D的问题，要么想解决一个3D的问题，但不会同时解决这两个问题）。然而，它展示了一种机制，我们可以在一个地方简单地改变我们想要的维度，从而迫使编译器为我们要求的维度重新编译独立的类模板。这里的重点在于，我们只需要改变一个地方。这使得在计算速度较快的2D环境下调试程序变得非常简单，然后将一个地方切换到3，在3D环境下运行计算量大得多的程序，进行 "真实 "的计算。
 * 

 * 
 * 这两个区块中的每一个都用大括号括起来，以确保 <code>laplace_problem_2d</code> 这个变量在我们继续为3D情况分配内存之前就已经超出了范围（并释放了它所持有的内存）。如果没有额外的大括号， <code>laplace_problem_2d</code> 变量只会在函数结束时被销毁，也就是在运行完3d问题后被销毁，而且会在3d运行时不必要地占用内存，而实际使用它。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   { 
 *     Step4<2> laplace_problem_2d; 
 *     laplace_problem_2d.run(); 
 *   } 
 * 
 *   { 
 *     Step4<3> laplace_problem_3d; 
 *     laplace_problem_3d.run(); 
 *   } 
 * 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-4/doc/results.dox



<a name="Results"></a><h1>Results</h1>



程序的输出看起来如下（迭代次数可能会有一到两次的变化，这取决于你的计算机，因为这通常取决于浮点运算的舍入精度，而这在不同的处理器之间是不同的）。

@code
Solving problem in 2 space dimensions.
   Number of active cells: 256
   Total number of cells: 341
   Number of degrees of freedom: 289
   26 CG iterations needed to obtain convergence.
Solving problem in 3 space dimensions.
   Number of active cells: 4096
   Total number of cells: 4681
   Number of degrees of freedom: 4913
   30 CG iterations needed to obtain convergence.
@endcode

很明显，在三个空间维度中，单元格的数量，因此也是自由度的数量要高得多。这里看不到的是，除了矩阵中更多的行和列之外，在三个空间维度中，矩阵的每一行也有明显更多的条目。这就导致了解方程组时需要付出更多的数值努力，当你实际运行程序时，你可以从两个求解步骤的运行时间中感受到这一点。




该程序产生两个文件。   <code>solution-2d.vtk</code> 和 <code>solution-3d.vtk</code> ，可以用VisIt或Paraview程序查看（如果你没有这些程序，你可以很容易地在程序中改变输出格式，使你更容易查看）。解决方案的可视化是一门艺术，但它也可以很有趣，所以你应该玩一玩你最喜欢的可视化工具，熟悉它的功能。下面是我想出的2D解决方案。

<p align="center">  <img src="https://www.dealii.org/images/steps/developer/step-4.solution-2d.png" alt="">   </p>  。

(  @dealiiVideoLectureSeeAlso{11,32})  图片显示了所考虑的问题的解决方案，是一个三维图。可以看出，该解在域的内部几乎是平的，而在边界附近有较高的曲率。当然，这是因为对于拉普拉斯方程来说，解的曲率等于右手边，而右手边被选为四次多项式，在内部几乎为零，只有在接近域的边界时才急剧上升；右手边函数的最大值在域的角落，在那里解的移动也最迅速。很高兴看到解沿着域的边界遵循理想的二次边界值。将计算出的解与分析出的解进行验证也是很有用的。关于这一技术的解释，请参见步骤7。

另一方面，尽管图片中没有明确显示网格线，但你可以看到它们在解决方案中的小疙瘩。这清楚地表明，解决方案还没有被计算到非常高的精度，为了得到更好的解决方案，我们可能必须在更细的网格上进行计算。

在三个空间维度上，可视化就比较困难了。左图显示了解决方案和它在域的表面上计算出来的网格。这很好，但它的缺点是完全掩盖了内部的情况。右图是通过显示解的恒定值的表面（如左上角的图例所示），试图将内部的情况也可视化。如果我们把各个表面弄得稍微透明一些，这样就有可能透过它们看到后面的东西，那么等值面图片看起来就最好了。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-4.solution-3d.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-4.contours-3d.png" alt="">
    </td>
  </tr>
</table> 

 @note  关于可视化的最后一句话：可视化的想法是给人以洞察力，这与显示信息是不同的。特别是，在一张图片上很容易显示过多的信息，但在显示更多的信息的同时，也使人们更难收集到洞察力。举个例子，我用来生成这些图片的程序，VisIt，默认情况下在每个轴上都有刻度线，在 $x$ 轴上贴上一个大胖标签 "X轴"，其他轴也是如此，在左上方显示提取数据的文件名，在右下方显示用户的名字以及时间和日期。这些在这里都不重要：轴同样容易辨认，因为左下方的三脚架仍然可见，而且我们从程序中知道域是 $[-1,1]^3$ ，所以不需要刻度线。因此，我关掉了图片中所有不相干的东西：可视化的艺术在于把图片缩减到那些对看清自己想看的东西很重要的部分，而不是其他。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>



基本上，玩这个程序的可能性与前一个程序相同，只是它们现在也适用于3D情况。为了获得灵感，请阅读<a href="step_3.html#extensions"
target="body">possible extensions in the documentation of step 3</a>。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-4.cc"
*/
