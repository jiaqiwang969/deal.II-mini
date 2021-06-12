/**
@page step_53 The step-53 tutorial program
This tutorial depends on step-49.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Wheregeometryandmeshesintersect">Where geometry and meshes intersect</a>
        <li><a href="#Theexamplecase">The example case</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#DescribingtopographyAfricaTopography">Describing topography: AfricaTopography</a>
        <li><a href="#DescribingthegeometryAfricaGeometry">Describing the geometry: AfricaGeometry</a>
        <li><a href="#Creatingthemesh">Creating the mesh</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Issueswithadaptivelyrefinedmeshesgeneratedthisway">Issues with adaptively refined meshes generated this way</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-53/doc/intro.dox

 <br> 

<i>This program was contributed by Wolfgang Bangerth and Luca Heltai, using
data provided by D. Sarah Stamps.</i>

 @note  这个程序阐述了几何学的概念和实现它的类。这些类被归入 @ref
manifold "三角形的流形描述 "的文档模块。其他信息见那里。

 @note  本教程也可作为Jupyter Python笔记本，使用deal.II python接口。该笔记本与原始的C++程序在同一目录下可用。渲染的笔记本也可以在<a
href="https://github.com/dealii/dealii/blob/master/example/step-53/step-53.ipynb">github</a>上查看。


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


现实问题的偏微分方程往往是在具有复杂几何形状的域上提出的。为了提供几个例子，请考虑这些情况。

- 在有限元方法的两个可以说是最重要的工业应用中，空气动力学和更普遍的流体动力学是其中之一。今天，计算机模拟被用于每架飞机、汽车、火车和船舶的设计。在这些情况下，提出偏微分方程的领域是飞机周围的空气，包括机翼、襟翼和发动机；汽车周围的空气，包括车轮、轮井、后视镜，在赛车的情况下，还有各种空气动力设备；火车周围的空气，包括车轮和车厢之间的空隙。在船舶的情况下，该领域是指有船舵和螺旋桨的船舶周围的水。

- 有限元方法的两大应用中的另一个是结构工程，其领域是桥梁、飞机机舱和机翼，以及其他通常形状复杂的固体物体。

- 有限元建模也经常被用来描述地震波的产生和传播。在这些情况下，人们需要准确地表示地壳中断层的几何形状。由于断层相交，有一定的倾角，而且往往不完全是直线，所以域经常是非常复杂的。我们可以举出更多复杂几何形状的例子，在这些例子中，我们要提出和解决一个偏微分方程。这表明，"真实 "世界比我们在这之前的几乎所有教程中所展示的要复杂得多。

因此，这个程序致力于展示如何用具体的应用来处理复杂的几何图形。特别是，它所展示的是我们如何使网格符合我们想要解决的领域。另一方面，该程序没有展示的是如何为一个领域创建一个粗的。获得粗网格的过程被称为 "网格生成"，有许多高质量的程序在这方面做得比我们做得更好。然而，deal.II确实有能力读取由网格生成器生成的多种格式的网格，然后使其适合给定的形状，可以通过变形网格或多次细化网格直至适合。从http://www.dealii.org/ 引用的deal.II常见问题页面提供了网格生成器的资源。




<a name="Wheregeometryandmeshesintersect"></a><h3>Where geometry and meshes intersect</h3>


让我们假设你有一个复杂的领域，并且你已经有一个粗略的网格，在某种程度上代表了这个领域的一般特征。那么在两种情况下，有必要向deal.II程序描述你的几何形状的细节。

- 网格细化。每当一个单元被细化时，有必要在三角网中引入新的顶点。在最简单的情况下，我们假设构成三角网的对象是直线段、双线性表面或三线性体。然后，下一个顶点被简单地放在旧顶点的中间。然而，对于弯曲的边界，或者如果我们想解决一个嵌入高维空间的弯曲的低维流形上的PDE，这是不充分的，因为它将不尊重实际的几何。因此，我们将不得不告诉Triangulation在哪里放置新的点。

- 积分。当使用高阶有限元方法时，经常需要使用边界的曲线近似来计算积分，即把单元的每个边缘或面描述为曲线，而不是直线段或双线性补丁。当然，在积分边界条款时也是如此（例如，不均匀的诺伊曼边界条件）。为了整合的目的，各种Mapping类提供了从参考单元到实际单元的转换。

在这两种情况下，我们需要一种方法来提供关于单个单元、其面和边的域的几何信息。这就是Manifold类开始发挥作用的地方。Manifold是一个抽象的基类，它只定义了一个接口，Triangulation和Mapping类可以通过这个接口查询领域的几何信息。从概念上讲，Manifold看待世界的方式与数学分支学科几何学看待世界的方式并无二致：域本质上只是一个点的集合，以某种方式配备了点之间的距离概念，这样我们就可以在其他一些点的 "中间 "获得一个点。

deal.II提供了一些实现Manifold所提供的接口的类，用于各种常见的几何形状。另一方面，在这个程序中，我们将只考虑一种非常常见的、简单得多的情况，即我们要解决的域的（一部分）可以通过转换一个简单得多的域（我们将称之为 "参考域"）来描述的情况。在数学语言中，这意味着该（部分）域是一个<a
href="http://en.wikipedia.org/wiki/Chart_%28topology%29">chart</a>。图表由一个平滑函数描述，该函数从较简单的域映射到图表（"向前推 "函数）及其逆向（"向后拉 "函数）。如果域作为一个整体不是一个图表（例如，球体的表面），那么它通常可以被描述为一个图表的集合（例如，北半球和南半球各自是一个图表），然后域可以被描述为一个<a
href="http://en.wikipedia.org/wiki/Atlas_%28topology%29">atlas</a>。

如果一个域可以被分解成一个图集，那么我们需要做的就是为每个图集提供回拉和推送函数。在deal.II中，这意味着提供一个从ChartManifold派生的类，而这正是我们在这个程序中要做的。




<a name="Theexamplecase"></a><h3>The example case</h3>


为了说明如何在deal.II中使用图表描述几何形状，我们将考虑一个源于<a
href="https://aspect.geodynamics.org">ASPECT mantle convection code</a>的应用的案例，使用D. Sarah Stamps提供的数据集。在具体应用中，我们对描述<a
href="http://en.wikipedia.org/wiki/East_African_rift">East African Rift</a>下的地幔流动感兴趣，这是一个两个大陆板块漂移的区域。不拐弯抹角，我们想要描述的几何形状看起来是这样的。

 <img src="https://www.dealii.org/images/steps/developer/step-53.topo.png" alt=""> 

特别是，虽然你在这里看不到，但顶部的表面不仅仅是由高程着色，实际上，它是按照正确的地形变形的。虽然实际的应用在这里并不重要，但几何学是相关的。我们感兴趣的领域是地球的一部分，范围从表面到500公里的深度，从格林威治子午线以东26度到35度，从赤道以北5度到南10度。

这种对几何学的描述建议从一个盒子 $\hat U=[26,35]\times[-10,5]\times[-500000,0]$ 开始（以度、度数和米为单位），并提供一个地图 $\varphi$ ，以便 $\varphi^{-1}(\hat U)=\Omega$ 其中 $\Omega$ 是我们寻求的领域。  然后， $(\Omega,\varphi)$ 是一个图表， $\varphi$ 是回拉运算符， $\varphi^{-1}$ 是前推运算符。如果我们需要一个点 $q$ 是其他点 $q_i\in\Omega$ 的 "平均值"，那么ChartManifold类首先应用回拉得到 $\hat q_i=\varphi(q_i)$ ，将其平均到一个点 $\hat p$ ，然后计算出 $p=\varphi^{-1}(\hat p)$ 。

因此，我们这里的目标是实现一个描述 $\varphi$ 和 $\varphi^{-1}$ 的类。如果地球是一个球体，那么这并不困难：如果我们用 $(\hat \phi,\hat \theta,\hat d)$ 表示 $\hat U$ 的点（即经度向东计算，纬度向北计算，海拔相对于零深度），那么

@f[
  \mathbf x = \varphi^{-1}(\hat \phi,\hat \theta,\hat d)
  = (R+\hat d) (\cos\hat \phi\cos\hat \theta, \sin\hat \phi\cos\hat \theta, \sin\hat \theta)^T


@f]

提供直角坐标系中的坐标，其中 $R$ 是球体的半径。然而，地球不是一个球体。

<ol>  <li> 它在两极是扁平的，在赤道是较大的：半主轴比半副轴长约22公里。我们将使用<a href="http://en.wikipedia.org/wiki/WGS84">WGS 84</a>的地球形状参考标准来说明这一点。在WGS 84中，用于从经度、纬度和海拔获得直角坐标位置的公式是

@f[
  \mathbf x = \varphi_\text{WGS84}^{-1}(\phi,\theta,d)
  = \left(
    \begin{array}{c}
     (\bar R(\theta)+d) \cos\phi\cos\theta, \\
     (\bar R(\theta)+d) \sin\phi\cos\theta, \\
     ((1-e^2)\bar R(\theta)+d) \sin\theta
    \end{array}
    \right),


@f]

  其中  $\bar R(\theta)=\frac{R}{\sqrt{1-(e \sin\theta)^2}}$  ，而半径和椭圆度由  $R=6378137\text{m}, e=0.081819190842622$  给出。在这个公式中，我们假设正弦和余弦的参数是以度数而不是弧度来计算的（尽管我们将不得不在代码中改变这个假设）。

 <li>  它的地形是山脉和山谷的形式。我们将使用真实的地形数据来说明这一点（见下文对这些数据来源的描述）。使用这个数据集，我们可以在地球表面的经纬度网格上查找高程。从方框 $\hat U=[26,35]\times[-10,5]\times[-500000,0]$ 开始，我们将首先在垂直方向上拉伸它，然后再把它交给WGS 84函数：如果 $h(\hat\phi,\hat\theta)$ 是经度 $\hat\phi$ 和纬度 $\hat\theta$ 的高度，那么我们定义

@f[
  (\phi,\theta,d) =
  \varphi_\text{topo}^{-1}(\hat\phi,\hat\theta,\hat d)
  = \left(
      \hat\phi,
      \hat\theta,
      \hat d + \frac{\hat d+500000}{500000}h(\hat\phi,\hat\theta)
    \right).


@f]

  使用这个函数，盒子 $\hat U$ 的顶面被移到正确的地形上，底面保持原来的位置，中间的点被线性内插。   </ol> 

利用这两个函数，我们就可以将整个推送函数 $\varphi^{-1}: \hat U \rightarrow \Omega$ 定义为

@f[
  \mathbf x
  =
  \varphi^{-1}(\hat\phi,\hat\theta,\hat d)
  =
  \varphi_\text{WGS84}^{-1}(\varphi_\text{topo}^{-1}(\hat\phi,\hat\theta,\hat d)).


@f]

此外，我们将不得不定义这个函数的逆运算，即回拉运算，我们可以将其写为

@f[
  (\hat\phi,\hat\theta,\hat d)
  =
  \varphi(\mathbf x)
  =
  \varphi_\text{topo}(\varphi_\text{WGS84}(\mathbf x)).


@f]

我们可以通过倒置上面的公式得到这个函数的一个分量。

@f[
  (\hat\phi,\hat\theta,\hat d) =
  \varphi_\text{topo}(\phi,\theta,d)
  = \left(
      \phi,
      \theta,
      500000\frac{d-h(\phi,\theta)}{500000+h(\phi,\theta)}
    \right).


@f]

计算 $\varphi_\text{WGS84}(\mathbf x)$ 也是可能的，不过要笨拙得多。我们不会在这里展示这个公式，而是只提供程序中的实现。




<a name="Implementation"></a><h3>Implementation</h3>


在这个程序中，我们需要解决一些问题。在最大范围内，我们需要编写一个实现ChartManifold接口的类。这涉及到一个函数 <code>push_forward()</code> ，该函数在参考域 $\hat U$ 中取一个点，并使用上面概述的函数 $\varphi^{-1}$ 将其转换为实空间，以及实现 <code>pull_back()</code> 的反函数 $\varphi$  。我们将在下面的 <code>AfricaGeometry</code> 类中这样做，该类本质上看起来像这样。

@code
  class AfricaGeometry : public ChartManifold<3,3>
  {
  public:
    virtual
    Point<3>
    pull_back(const Point<3> &space_point) const;


    virtual
    Point<3>
    push_forward(const Point<3> &chart_point) const;


  private:
    ... some member variables and other member functions...;
  };
@endcode



上述转换有两个部分：WGS 84转换和地形转换。因此， <code>AfricaGeometry</code> 类将有额外的（非虚拟）成员函数 <code>AfricaGeometry::push_forward_wgs84()</code> 和 <code>AfricaGeometry::push_forward_topo()</code> 来实现这两部分，以及相应的回拉函数。

WGS 84的转换函数并不特别有趣（尽管它们实现的公式令人印象深刻）。更有趣的部分是地形变换。回顾一下，为此我们需要评估高程函数 $h(\hat\phi,\hat\theta)$  。当然，这没有公式。地球就是这样，人们能做的最好的事情就是从一些表格中查找海拔高度。事实上，这就是我们要做的。

我们使用的数据最初是由<a
href="http://en.wikipedia.org/wiki/Shuttle_Radar_Topography_Mission">Shuttle
Radar Topography Mission</a>创建的，是从美国地质调查局（USGS）下载的，并由D. Sarah Stamps处理，他还编写了WGS 84转换函数的初始版本。这样处理过的地形数据被储存在一个文件 <code>topography.txt.gz</code> 中，解压后看起来是这样的。

@code
6.983333 25.000000 700
6.983333 25.016667 692
6.983333 25.033333 701
6.983333 25.050000 695
6.983333 25.066667 710
6.983333 25.083333 702
...


-11.983333 35.950000 707


-11.983333 35.966667 687


-11.983333 35.983333 659
@endcode

数据格式为 <code>latitude longitude elevation</code> ，其中前两栏以赤道以北的度数和格林威治子午线以东的度数提供。最后一列是以WGS84零点以上的米数为单位。

在转换函数中，我们需要对给定的经度 $\hat\phi$ 和纬度 $\hat\theta$ 评估 $h(\hat\phi,\hat\theta)$  。一般来说，这个数据点是不可用的，我们将不得不在相邻的数据点之间进行内插。编写这样一个插值程序并不特别困难，但它有点乏味，而且容易出错。幸运的是，我们可以以某种方式将这个数据集塞进一个现有的类中。   Functions::InterpolatedUniformGridData  .不幸的是，这个类并不完全适合，所以我们需要绕过它。问题来自于我们初始化这个类的方式：在其最简单的形式下，它需要一个数值流，它假设在 $x-y$ 平面（或者，这里是 $\phi-\theta$ 平面）形成一个等距的网格。这就是它们在这里所做的，某种程度上：它们的顺序是纬度第一，经度第二；更尴尬的是，第一列从最大的数值开始往下数，而不是通常的其他方式。

现在，虽然教程程序是为了说明如何用deal.II编码，但它们不一定要满足与生产代码相同的质量标准。在生产代码中，我们会写一个函数来读取数据，并（i）自动确定第一列和第二列的外延，（ii）自动确定每个方向的数据点的数量，（iii）无论数据的排列顺序如何，都要进行插值，如果有必要的话，在读取和呈现给 Functions::InterpolatedUniformGridData 类之间切换顺序。

另一方面，辅导课程最好是短小精悍，展示关键点，而不是纠缠于不重要的方面，从而掩盖了我们真正想要展示的东西。因此，我们将允许自己有一点回旋余地。

- 由于这个程序只针对东非裂谷地区周围的特定几何形状，并且由于这正是数据文件所描述的区域，我们将在程序中硬编码有 $1139\times 660$ 个数据。

- 我们将硬编码数据的边界  $[-11.98333^\circ,6.983333^\circ]\times[25^\circ,35.98333^\circ]$  。

- 我们将对 Functions::InterpolatedUniformGridData 类撒谎：该类将只看到这个数据文件最后一列的数据，我们将假装数据的排列方式是：在第一个坐标方向上有1139个数据点，这些数据点按<i>ascending</i>的顺序排列，但在一个区间 $[-6.983333^\circ,11.98333^\circ]$ （不是否定的边界）。然后，当我们需要查询某个纬度 $\hat\theta$ 的东西时，我们可以向内插表类索取 $-\hat\theta$ 的数值。有了这个小技巧，我们就可以避免在从文件中读取数据的时候切换顺序。

所有这些都需要一个本质上看起来像这样的类。

@code
  class AfricaTopography
  {
  public:
    AfricaTopography ()
      :
      topography_data (...initialize somehow...)
    {}


    double value (const double lon, const double lat) const
    {
      return topography_data.value (Point<2>(-lat * 180/numbers::PI,
                                             lon * 180/numbers::PI));
    }


  private:
    const Functions::InterpolatedUniformGridData<2> topography_data;
  };
@endcode



注意 <code>value()</code> 函数如何否定了纬度。它还将我们在其他地方使用的格式 $\phi,\theta$ 转换为表格中使用的纬度-经度格式。最后，它的参数以弧度为单位，因为我们在程序中的其他地方也是这样做的，但随后将它们转换为用于表格查询的基于度的系统。正如你在下面的实现中所看到的，该函数还有几个（静态）成员函数，我们将在初始化 <code>topography_data</code> 成员变量时调用：该变量的类类型有一个构造函数，允许我们在构造时正确设置一切，而不是在以后填充数据，但这个构造函数需要一些不能就地构造的对象（至少在C++98中不能）。因此，我们要在初始化中传递的每个对象的构造都发生在一些静态成员函数中。

在讨论了我们要实施的事情的大体轮廓之后，让我们去看程序，并展示它在实践中是如何做的。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 让我们从这里需要的包含文件开始。显然，我们需要描述三角形的文件（  <code>tria.h</code>  ），以及允许我们创建和输出三角形的文件（  <code>grid_generator.h</code>  和  <code>grid_out.h</code>  ）。此外，我们需要声明Manifold和ChartManifold类的头文件，我们将需要这些类来描述几何体（ <code>manifold.h</code> ）。然后我们还需要以下头文件中的 GridTools::transform() 函数；这个函数的用途将在我们使用它时讨论。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_out.h> 
 * #include <deal.II/grid/manifold.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * 
 * @endcode
 * 
 * 其余的包含文件与读取地形数据有关。正如介绍中所解释的，我们将从一个文件中读取它，然后使用下面头文件中第一个声明的 Functions::InterpolatedUniformGridData  类。因为数据很大，所以我们读取的文件是以gzip压缩数据的形式存储的，我们利用BOOST提供的一些功能来直接读取gzipped数据。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/function_lib.h> 
 * 
 * #include <boost/iostreams/filtering_stream.hpp> 
 * #include <boost/iostreams/filter/gzip.hpp> 
 * #include <boost/iostreams/device/file.hpp> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * #include <memory> 
 * 
 * @endcode
 * 
 * 上事的最后部分是打开一个命名空间，把所有东西都放进去，然后把dealii命名空间导入其中。
 * 

 * 
 * 
 * @code
 * namespace Step53 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="DescribingtopographyAfricaTopography"></a> 
 * <h3>Describing topography: AfricaTopography</h3>
 * 

 * 
 * 这个程序的第一个重要部分是描述地形 $h(\hat phi,\hat \theta)$ 作为经度和纬度的函数的类。正如在介绍中所讨论的那样，我们在这里将使我们的生活更容易一些，不以最普遍的方式来写这个类，而是只为我们在这里感兴趣的特定目的来写：插值从一个非常具体的数据文件中获得的数据，该文件包含了关于世界上一个特定地区的信息，我们知道该地区的范围。
 * 

 * 
 * 该类的总体布局已经在上面讨论过了。下面是它的声明，包括我们在初始化 <code>topography_data</code> 成员变量时需要的三个静态成员函数。
 * 

 * 
 * 
 * @code
 *   class AfricaTopography 
 *   { 
 *   public: 
 *     AfricaTopography(); 
 * 
 *     double value(const double lon, const double lat) const; 
 * 
 *   private: 
 *     const Functions::InterpolatedUniformGridData<2> topography_data; 
 * 
 *     static std::vector<double> get_data(); 
 *   }; 
 * 
 * @endcode
 * 
 * 让我们来看看这个类的实现。该类的有趣部分是构造函数和 <code>value()</code> 函数。前者初始化了 Functions::InterpolatedUniformGridData 成员变量，我们将使用这个构造函数，它要求我们传入我们要插值的二维数据集的端点（这里由区间 $[-6.983333, 11.98333]$ 给出）。 ]，使用介绍中讨论的切换端点的技巧，和 $[25, 35.983333]$ ，都是以度数给出的），数据被分割成的区间数（纬度方向379，经度方向219，总共 $380\times 220$ 个数据点），和一个包含数据的表对象。然后，数据的大小当然是 $380\times 220$ ，我们通过提供一个迭代器给下面 std::vector 函数返回的 <code>get_data()</code> 对象的83,600个元素中的第一个来初始化它。注意，我们在这里调用的所有成员函数都是静态的，因为(i)它们不访问类的任何成员变量，(ii)因为它们是在对象没有完全初始化的时候调用的。
 * 

 * 
 * 
 * @code
 *   AfricaTopography::AfricaTopography() 
 *     : topography_data({{std::make_pair(-6.983333, 11.966667), 
 *                         std::make_pair(25, 35.95)}}, 
 *                       {{379, 219}}, 
 *                       Table<2, double>(380, 220, get_data().begin())) 
 *   {} 
 * 
 *   double AfricaTopography::value(const double lon, const double lat) const 
 *   { 
 *     return topography_data.value( 
 *       Point<2>(-lat * 180 / numbers::PI, lon * 180 / numbers::PI)); 
 *   } 
 * 
 * @endcode
 * 
 * 唯一一个更有意义的函数是 <code>get_data()</code> 函数。它返回一个临时向量，其中包含描述高度的所有83600个数据点，并从文件 <code>topography.txt.gz</code> 中读取。因为文件被gzip压缩了，所以我们不能直接通过类型为 std::ifstream, 的对象来读取它，但在BOOST库中有一些方便的方法（见http:www.boost.org），允许我们从压缩的文件中读取，而不用先在磁盘上解压缩。其结果是，基本上，只是另一个输入流，就所有的实际目的而言，看起来就像我们一直使用的那些输入流。
 * 

 * 
 * 当读取数据时，我们读取三列数据，但忽略了前两列。最后一列的数据被附加到一个数组中，我们返回的数组将被复制到 <code>topography_data</code> 的表中，并被初始化。由于BOOST.iostreams库在输入文件不存在、不可读或不包含正确的数据行数时没有提供非常有用的异常，我们捕捉它可能产生的所有异常并创建我们自己的异常。为此，在 <code>catch</code> 子句中，我们让程序运行到一个 <code>AssertThrow(false, ...)</code> 语句中。由于条件总是假的，这总是会触发一个异常。换句话说，这相当于写了 <code>throw ExcMessage("...")</code> ，但它也填补了异常对象中的某些字段，这些字段以后会被打印在屏幕上，识别出发生异常的函数、文件和行。
 * 

 * 
 * 
 * @code
 *   std::vector<double> AfricaTopography::get_data() 
 *   { 
 *     std::vector<double> data; 
 * 
 * @endcode
 * 
 * 创建一个流，我们从gzipped数据中读取
 * 

 * 
 * 
 * @code
 *     boost::iostreams::filtering_istream in; 
 *     in.push(boost::iostreams::basic_gzip_decompressor<>()); 
 *     in.push(boost::iostreams::file_source("topography.txt.gz")); 
 * 
 *     for (unsigned int line = 0; line < 83600; ++line) 
 *       { 
 *         try 
 *           { 
 *             double lat, lon, elevation; 
 *             in >> lat >> lon >> elevation; 
 * 
 *             data.push_back(elevation); 
 *           } 
 *         catch (...) 
 *           { 
 *             AssertThrow(false, 
 *                         ExcMessage("Could not read all 83,600 data points " 
 *                                    "from the file <topography.txt.gz>!")); 
 *           } 
 *       } 
 * 
 *     return data; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="DescribingthegeometryAfricaGeometry"></a> 
 * <h3>Describing the geometry: AfricaGeometry</h3>
 * 

 * 
 * 下面的类是本程序的主类。它的结构已经在介绍中详细描述过了，不需要再多做介绍。
 * 

 * 
 * 
 * @code
 *   class AfricaGeometry : public ChartManifold<3, 3> 
 *   { 
 *   public: 
 *     virtual Point<3> pull_back(const Point<3> &space_point) const override; 
 * 
 *     virtual Point<3> push_forward(const Point<3> &chart_point) const override; 
 * 
 *     virtual std::unique_ptr<Manifold<3, 3>> clone() const override; 
 * 
 *   private: 
 *     static const double R; 
 *     static const double ellipticity; 
 * 
 *     const AfricaTopography topography; 
 * 
 *     Point<3> push_forward_wgs84(const Point<3> &phi_theta_d) const; 
 *     Point<3> pull_back_wgs84(const Point<3> &x) const; 
 * 
 *     Point<3> push_forward_topo(const Point<3> &phi_theta_d_hat) const; 
 *     Point<3> pull_back_topo(const Point<3> &phi_theta_d) const; 
 *   }; 
 * 
 *   const double AfricaGeometry::R           = 6378137; 
 *   const double AfricaGeometry::ellipticity = 8.1819190842622e-2; 
 * 
 * @endcode
 * 
 * 如果你读过介绍，实现起来也是非常简单的。特别是，回拉和前推函数都只是WGS 84和地形图映射各自函数的串联。
 * 

 * 
 * 
 * @code
 *   Point<3> AfricaGeometry::pull_back(const Point<3> &space_point) const 
 *   { 
 *     return pull_back_topo(pull_back_wgs84(space_point)); 
 *   } 
 * 
 *   Point<3> AfricaGeometry::push_forward(const Point<3> &chart_point) const 
 *   { 
 *     return push_forward_wgs84(push_forward_topo(chart_point)); 
 *   } 
 * 
 * @endcode
 * 
 * 下一个函数是Manifold基类的接口所要求的，它允许克隆AfricaGeometry类。注意，虽然该函数返回一个  `std::unique_ptr<Manifold<3,3>>`,  我们在内部创建了一个 `unique_ptr<AfricaGeometry>`。换句话说，这个库需要一个指向基类的指针，我们通过创建一个指向派生类的指针来提供这个指针。
 * 

 * 
 * 
 * @code
 *   std::unique_ptr<Manifold<3, 3>> AfricaGeometry::clone() const 
 *   { 
 *     return std::make_unique<AfricaGeometry>(); 
 *   } 
 * 
 * @endcode
 * 
 * 下面的两个函数就定义了对应于地球WGS84参考形状的正向和反向变换。正向变换遵循介绍中所示的公式。反变换要复杂得多，至少不是直观的。它还存在一个问题，即它返回一个角度，在函数结束时，如果它应该从那里逃出来，我们需要将其夹回区间 $[0,2\pi]$ 。
 * 

 * 
 * 
 * @code
 *   Point<3> AfricaGeometry::push_forward_wgs84(const Point<3> &phi_theta_d) const 
 *   { 
 *     const double phi   = phi_theta_d[0]; 
 *     const double theta = phi_theta_d[1]; 
 *     const double d     = phi_theta_d[2]; 
 * 
 *     const double R_bar = R / std::sqrt(1 - (ellipticity * ellipticity * 
 *                                             std::sin(theta) * std::sin(theta))); 
 * 
 *     return {(R_bar + d) * std::cos(phi) * std::cos(theta), 
 *             (R_bar + d) * std::sin(phi) * std::cos(theta), 
 *             ((1 - ellipticity * ellipticity) * R_bar + d) * std::sin(theta)}; 
 *   } 
 * 
 *   Point<3> AfricaGeometry::pull_back_wgs84(const Point<3> &x) const 
 *   { 
 *     const double b   = std::sqrt(R * R * (1 - ellipticity * ellipticity)); 
 *     const double ep  = std::sqrt((R * R - b * b) / (b * b)); 
 *     const double p   = std::sqrt(x(0) * x(0) + x(1) * x(1)); 
 *     const double th  = std::atan2(R * x(2), b * p); 
 *     const double phi = std::atan2(x(1), x(0)); 
 *     const double theta = 
 *       std::atan2(x(2) + ep * ep * b * std::pow(std::sin(th), 3), 
 *                  (p - 
 *                   (ellipticity * ellipticity * R * std::pow(std::cos(th), 3)))); 
 *     const double R_bar = 
 *       R / (std::sqrt(1 - ellipticity * ellipticity * std::sin(theta) * 
 *                            std::sin(theta))); 
 *     const double R_plus_d = p / std::cos(theta); 
 * 
 *     Point<3> phi_theta_d; 
 *     if (phi < 0) 
 *       phi_theta_d[0] = phi + 2 * numbers::PI; 
 *     else if (phi > 2 * numbers::PI) 
 *       phi_theta_d[0] = phi - 2 * numbers::PI; 
 *     else 
 *       phi_theta_d[0] = phi; 
 *     phi_theta_d[1] = theta; 
 *     phi_theta_d[2] = R_plus_d - R_bar; 
 *     return phi_theta_d; 
 *   } 
 * 
 * @endcode
 * 
 * 与此相反，地形变换完全按照介绍中的描述进行。因此，没有什么可以补充的。
 * 

 * 
 * 
 * @code
 *   Point<3> 
 *   AfricaGeometry::push_forward_topo(const Point<3> &phi_theta_d_hat) const 
 *   { 
 *     const double d_hat = phi_theta_d_hat[2]; 
 *     const double h = topography.value(phi_theta_d_hat[0], phi_theta_d_hat[1]); 
 *     const double d = d_hat + (d_hat + 500000) / 500000 * h; 
 *     return {phi_theta_d_hat[0], phi_theta_d_hat[1], d}; 
 *   } 
 * 
 *   Point<3> AfricaGeometry::pull_back_topo(const Point<3> &phi_theta_d) const 
 *   { 
 *     const double d     = phi_theta_d[2]; 
 *     const double h     = topography.value(phi_theta_d[0], phi_theta_d[1]); 
 *     const double d_hat = 500000 * (d - h) / (500000 + h); 
 *     return {phi_theta_d[0], phi_theta_d[1], d_hat}; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Creatingthemesh"></a> 
 * <h3>Creating the mesh</h3>
 * 

 * 
 * 在描述了几何体的属性之后，现在是处理用于离散它的网格的时候了。为此，我们为几何体和三角形创建对象，然后继续创建一个与参考域 $1\times 2\times 1$ 相对应的 $\hat U=[26,35]\times[-10,5]\times[-500000,0]$ 矩形网格。我们选择这个数目的细分，因为它导致了单元格大致上像立方体，而不是在某个方向上被拉伸。
 * 

 * 
 * 当然，我们实际上对参考域的网格划分不感兴趣。我们感兴趣的是对真实域的网格划分。因此，我们将使用 GridTools::transform() 函数，它只是根据一个给定的变换来移动三角形的每个点。它想要的变换函数是一个将参考域中的一个点作为其单一参数的函数，并返回我们想要映射到的域中的相应位置。当然，这正是我们使用的几何学的前推函数。我们用一个lambda函数来包装它，以获得转换所需的那种函数对象。
 * 

 * 
 * 
 * @code
 *   void run() 
 *   { 
 *     AfricaGeometry   geometry; 
 *     Triangulation<3> triangulation; 
 * 
 *     { 
 *       const Point<3> corner_points[2] = { 
 *         Point<3>(26 * numbers::PI / 180, -10 * numbers::PI / 180, -500000), 
 *         Point<3>(35 * numbers::PI / 180, 5 * numbers::PI / 180, 0)}; 
 *       std::vector<unsigned int> subdivisions(3); 
 *       subdivisions[0] = 1; 
 *       subdivisions[1] = 2; 
 *       subdivisions[2] = 1; 
 *       GridGenerator::subdivided_hyper_rectangle( 
 *         triangulation, subdivisions, corner_points[0], corner_points[1], true); 
 * 
 *       GridTools::transform( 
 *         [&geometry](const Point<3> &chart_point) { 
 *           return geometry.push_forward(chart_point); 
 *         }, 
 *         triangulation); 
 *     } 
 * 
 * @endcode
 * 
 * 下一步是向三角计算说明，在细化网格时，每当需要一个新的点时，都要使用我们的几何对象。我们通过告诉三角计算对所有流形指示器为零的物体使用我们的几何体，然后继续用流形指示器为零标记所有单元及其边界面和边。这确保了三角计算在每次需要新的顶点时都会参考我们的几何对象。由于流形指标是由母体继承给子体的，这也会在几个递归细化步骤之后发生。
 * 

 * 
 * 
 * @code
 *     triangulation.set_manifold(0, geometry); 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       cell->set_all_manifold_ids(0); 
 * 
 * @endcode
 * 
 * 最后一步是在最初的 $1\times 2\times 1$ 粗略网格之外细化该网格。我们可以在全局范围内细化若干次，但由于本教程程序的目的，我们实际上只对靠近表面的情况感兴趣，所以我们只是对所有在边界上有一个指标为5的面的单元进行6次细化。在我们上面使用的 GridGenerator::subdivided_hyper_rectangle() 函数的文档中查找，发现边界指标5对应于域的顶面（这就是上面调用 GridGenerator::subdivided_hyper_rectangle() 的最后一个 <code>true</code> 参数的含义：通过给每个边界分配一个独特的边界指标来给边界 "着色"）。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < 6; ++i) 
 *       { 
 *         for (const auto &cell : triangulation.active_cell_iterators()) 
 *           for (const auto &face : cell->face_iterators()) 
 *             if (face->boundary_id() == 5) 
 *               { 
 *                 cell->set_refine_flag(); 
 *                 break; 
 *               } 
 *         triangulation.execute_coarsening_and_refinement(); 
 * 
 *         std::cout << "Refinement step " << i + 1 << ": " 
 *                   << triangulation.n_active_cells() << " cells, " 
 *                   << GridTools::minimal_cell_diameter(triangulation) / 1000 
 *                   << "km minimal cell diameter" << std::endl; 
 *       } 
 * 
 * @endcode
 * 
 * 做完这一切，我们现在可以将网格输出到一个自己的文件中。
 * 

 * 
 * 
 * @code
 *     const std::string filename = "mesh.vtu"; 
 *     std::ofstream     out(filename); 
 *     GridOut           grid_out; 
 *     grid_out.write_vtu(triangulation, out); 
 *   } 
 * } // namespace Step53 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * 最后是主函数，它采用了从  step-6  开始的所有教程程序中使用的相同方案。这里没有什么可做的，只需要调用单一的  <code>run()</code>  函数。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       Step53::run(); 
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
 * } 
 * 
 * 
 * @endcode
examples/step-53/doc/results.dox



<a name="Results"></a><h1>Results</h1>


运行程序会产生一个网格文件 <code>mesh.vtu</code> ，我们可以用任何可以读取VTU文件格式的常规可视化程序来进行可视化。如果只看网格本身，实际上很难看到任何不只是看起来像一块完全圆形的球体的东西（尽管如果修改程序，使其确实产生一个球体，并同时看它们，整体球体和WGS 84形状之间的差异是相当明显的）。很明显，地球实际上是一个相当平坦的地方。当然，我们已经从卫星图片中知道了这一点。然而，我们可以通过对细胞的体积进行着色来找出更多的东西。这既产生了沿顶面的轻微色调变化，也为可视化程序提供了应用其着色算法的东西（因为单元的顶面现在不再只是与球体的切线，而是倾斜的）。

 <img src="https://www.dealii.org/images/steps/developer/step-53.mesh.png" alt=""> 

然而，至少就视觉化而言，这仍然不是太令人印象深刻。相反，让我们以一种可视化的方式，使我们显示出沿顶面的实际海拔。换句话说，我们想要一张这样的图片，有难以置信的细节。

 <img src="https://www.dealii.org/images/steps/developer/step-53.topo.png" alt=""> 

这张照片的放大显示了相当清楚的垂直位移（这里，从西-西北方向看裂谷上空，<a href="http://en.wikipedia.org/wiki/Mount_Stanley">Mount Stanley</a>、<a href="http://en.wikipedia.org/wiki/Mount_Speke">Mount Speke</a>和<a href="http://en.wikipedia.org/wiki/Mount_Baker_%28Uganda%29">Mount Baker</a>的三座山峰在<a href="http://en.wikipedia.org/wiki/Rwenzori_Mountains">Rwenzori Range</a>、<a href="http://en.wikipedia.org/wiki/Lake_George_%28Uganda%29">Lake
George</a>和向<a href="http://en.wikipedia.org/wiki/Lake_Victoria">Lake Victoria</a>的巨大平坦处）。

 <img src="https://www.dealii.org/images/steps/developer/step-53.topozoom.png" alt=""> 


这些图片是经过三个小的修改后产生的。<ol>  <li>  在这两张图片中，第一张图片的顶面增加了第七个网格细化，第二张图片总共增加了九个。在第二张图片中，水平方向的网格大小约为1.5km，垂直方向的网格大小略低于1km。(这张图也是用一个更有分辨率的数据集制作的；但是，它太大了，不能作为教程的一部分分发)。

    <li>  增加以下函数，在给定一个点 <code>x</code> 时，通过将该点转换为参考WGS 84坐标并只保留深度变量来计算海拔高度（因此，该函数是 <code>AfricaGeometry::pull_back_wgs84()</code> 函数的简化版本）。

@code
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>



double get_elevation (const Point<3> &x)
  {
    const double R           = 6378137;
    const double ellipticity = 8.1819190842622e-2;


    const double b     = std::sqrt(R * R * (1 - ellipticity * ellipticity));
    const double ep    = std::sqrt((R * R - b * b) / (b * b));
    const double p     = std::sqrt(x(0) * x(0) + x(1) * x(1));
    const double th    = std::atan2(R * x(2), b * p);
    const double theta = std::atan2((x(2) + ep * ep * b * std::sin(th) * std::sin(th) * std::sin(th)),
                                      (p - (ellipticity * ellipticity * R  * (std::cos(th) * std::cos(th) * std::cos(th)))));
    const double R_bar = R / (std::sqrt(1 - ellipticity * ellipticity * std::sin(theta) * std::sin(theta)));
    const double R_plus_d = p / std::cos(theta);


    return R_plus_d - R_bar;
  }
@endcode



    <li>  在 <code>run()</code> 函数的底部添加以下一块。

@code
      FE_Q<3>       fe(1);
      DoFHandler<3> dof_handler (triangulation);
      dof_handler.distribute_dofs(fe);


      Vector<double> elevation (dof_handler.n_dofs());
      {
        std::map<unsigned int,double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 5,
                                                 ScalarFunctionFromFunctionObject<3>(get_elevation),
                                                 boundary_values);
        for (std::map<unsigned int,double>::const_iterator p = boundary_values.begin();
             p!=boundary_values.end(); ++p)
          elevation[p->first] = p->second;
      }


      DataOut<3>    data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector (elevation, "elevation");
      data_out.build_patches();


      std::ofstream out ("data.vtu");
      data_out.write_vtu (out);
@endcode

 </ol>  这最后一段代码首先在网格上创建一个 $Q_1$ 有限元空间。然后使用 VectorTools::interpolate_boundary_values() 对顶部边界的每个节点（边界指标为5的节点）进行高程函数的评估。我们在这里用ScalarFunctionFromFunctionObject类来包装对 <code>get_elevation()</code> 的调用，使一个普通的C++函数看起来像一个派生自Function类的对象，我们想在 VectorTools::interpolate_boundary_values(). 中使用。然后像往常一样用DataOut输出这个向量，并可以如上图所示进行可视化。




<a name="Issueswithadaptivelyrefinedmeshesgeneratedthisway"></a><h3>Issues with adaptively refined meshes generated this way</h3>


如果你放大上图所示的网格并仔细观察，你会发现在悬空节点处，连接到悬空节点的两条小边与相邻单元的大边的位置不完全相同。这可以通过使用不同的表面描述来更清楚地显示出来，在这种描述中，我们放大了垂直地形以增强效果（由Alexander Grayver提供）。

 <img src="https://www.dealii.org/images/steps/developer/step-53.smooth-geometry.png" alt=""> 

那么这里发生了什么？部分原因是，这只是视觉化的结果，但也有一个潜在的真正原因。

 <ul>   <li>  当你使用任何一个常见的可视化程序对网格进行可视化时，它们真正显示的只是一组在三维空间中被绘制成直线的边缘。这是因为几乎所有用于可视化的数据文件格式都只将六面体单元描述为三维空间中八个顶点的集合，而不允许任何更复杂的描述。这就是为什么 DataOut::build_patches() 需要一个可以设置为大于1的参数的主要原因）。这些线性边缘可能是你进行实际计算的单元格的边缘，也可能不是，这取决于你在使用FEValues进行积分时使用何种映射。当然，在默认情况下，FEValues使用的是线性映射（即MappingQ1类的对象），在这种情况下，一个3D单元确实完全由其8个顶点描述，它所填充的体积是这些点之间的三线插值，从而产生了线性边缘。但是，你也可以使用三次方、三次立方、甚至更高阶的映射，在这些情况下，每个单元的体积将由二次方、三次方或高阶多项式曲线来限定。然而，你只能在可视化程序中看到这些带有线性边缘的曲线，因为如前所述，文件格式不允许描述细胞的真实几何形状。

    <li>  也就是说，为了简单起见，让我们假设你确实在使用三线性映射，那么上面显示的图像就是你形成积分的单元的忠实代表。在这种情况下，一般来说，悬空节点上的小单元并不与大单元紧密贴合，而是留有间隙，或者可能与大单元相交。这是为什么呢？   因为当三角测量需要在它想要细化的边缘上增加一个新的顶点时，它会询问流形描述这个新的顶点应该在哪里，流形描述通过（在从ChartManifold派生的几何体的情况下）将线的相邻点拉回到参考域，平均它们的位置，并将这个新的位置推到真实域，从而适当地返回这样一个点。但是这个新的位置通常不是沿着相邻顶点之间的直线（在实空间），因此，形成精炼边缘的两条小直线并不完全位于形成悬挂节点的未精炼边的一条大直线上。   </ul> 

如果你使用MappingQ类的高阶映射，情况会稍微复杂一些，但没有根本的不同。我们暂且采用二次映射（高阶映射也没有什么根本性的变化）。那么你需要把你所整合的单元格的每条边想象成一条二次曲线，尽管你实际上不会看到可视化程序以这种方式绘制。但请先想象一下。那么，MappingQ采用哪条二次曲线呢？它是经过边缘末端的两个顶点以及中间的一个点的二次曲线，它从流形中查询。在未精炼一侧的长边的情况下，这当然正是悬空节点的位置，所以描述长边的二次曲线确实经过了悬空节点，这与线性映射的情况不同。但是两条小边也是二次曲线；例如，左边的小边会经过长边的左顶点和悬挂节点，再加上它从流形中查询到的一个点。因为如前所述，流形沿左小边中途返回的点很少完全在描述长边的二次曲线上，二次短边通常不会与二次长边的左半边重合，右短边也是如此。换句话说，同样的，大单元格和它在悬挂节点上的小邻居的几何形状并不相依相偎。

这一切引出了两个问题：第一，这是否重要，第二，这是否可以解决。让我们在下文中讨论这些问题。

 <ul>   <li>  这重要吗？几乎可以肯定的是，这取决于你所求解的方程。例如，众所周知，在复杂几何体上求解气体动力学的欧拉方程需要高度精确的边界描述，以确保衡量靠近边界的流动的量的收敛。另一方面，具有椭圆成分的方程（如拉普拉斯方程或斯托克斯方程）通常对这些问题比较宽容：无论如何都要做正交来近似积分，鉴于每个悬空节点的重叠或间隙的体积只有 ${\cal O}(h^d)$ ，进一步近似几何可能不会像人们担心的那样造成伤害。]，即使是线性映射，对于度数为 $p$ 的映射，也只有 ${\cal
  O}(h^{d+p-1})$ 。 (你可以通过考虑到在2D中，间隙/重叠是一个三角形，底 $h$ ，高 ${\cal
  O}(h)$ ；在3D中，它是一个类似金字塔的结构，底面积 $h^2$ ，高 ${\cal O}(h)$  。类似的考虑也适用于高阶映射，其中空隙/重叠的高度为  ${\cal O}(h^p)$  .)换句话说，如果你使用线性元素的线性映射，你所积分的体积的误差已经与使用通常的高斯正交的积分误差处于同一水平。当然，对于更高阶的元素，人们将不得不选择匹配的映射对象。

  关于为什么不值得过分担心这个问题的另一个观点是，在数值分析界肯定没有人说这些问题是使用复杂几何体时需要注意的一个主要问题。如果从业者之间似乎并不经常讨论这个问题，如果有的话，那么它至少不是人们认定的一个普遍问题。

  这个问题与在弯曲的边界上有悬挂的节点没有什么不同，在这种情况下，边界的几何描述通常会将悬挂的节点拉到边界上，而大的边缘仍然是直的，使得相邻的小单元和大单元不能相互匹配。虽然这种行为从一开始就存在于deal.II中，在流形描述出现之前的15年，但在邮件列表的讨论或与同事的交谈中，它从未出现过。

    <li>  能否修复？原则上是的，但这是一个复杂的问题。我们暂且假设我们只会使用MappingQ1类，也就是线性映射。在这种情况下，每当三角化类需要沿着一条将成为悬空节点的边建立一个新的顶点时，它就会直接取相邻顶点的平均值<i>in real
  space</i>，也就是说，不问流形的描述。这样一来，该点就位于长直边上，两条短直边将与一条长直边匹配。只有当所有相邻的单元都被细化，并且该点不再是一个悬挂的节点时，我们才会用通过流形得到的坐标来替换它的坐标。这在实现上可能比较麻烦，但肯定是可行的。

  更复杂的问题出现了，因为人们可能想使用高阶的MappingQ对象。在这种情况下，Triangulation类可以自由选择悬挂节点的位置（因为长边的二次曲线可以选择通过悬挂节点的方式），但是MappingQ类在确定中边点的位置时，必须确保如果该边是相邻的较粗单元的长边的一半，那么中点不能从流形中获得，而必须沿着长二次边选择。对于立方体（和所有其他奇数）映射，这个问题又有点复杂了，因为人们通常安排立方体边沿着边的1/3和2/3点走，因此必然要通过悬挂的节点，但这可能是可以解决的。在任何情况下，即使如此，这也有两个问题。

  - 在细化三角函数时，三角函数类不可能知道将使用什么映射。事实上，在同一个程序中，一个三角函数在不同的情况下被使用的情况并不罕见。如果使用的映射决定了我们是否可以自由选择一个点，那么，三角剖分应该如何定位新的顶点？

  - 映射是纯粹的局部构造：它们只对孤立的单元起作用，而这正是有限元方法的重要特征之一。要问一条边的一个顶点是否是悬空节点，需要查询一个单元的邻域；此外，这样的查询不只是涉及到三维中一个单元的6个面的邻域，而是可能需要遍历与一条边相连的大量的其他单元。即使可以做到这一点，人们仍然需要根据邻域的样子做不同的事情，产生的代码可能非常复杂，难以维护，而且可能很慢。

  因此，至少在目前，这些想法都没有被实施。这导致了不连续的几何形状的不良后果，但是，正如上面所讨论的，这种影响在实际操作中似乎并不构成问题。

 </ul> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-53.cc"
*/
