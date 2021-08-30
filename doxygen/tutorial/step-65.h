/**
@page step_65 The step-65 tutorial program
This tutorial depends on step-49.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Workingwithmanifolds">Working with manifolds</a>
      <ul>
        <li><a href="#Whatwewant">What we want</a>
      </ul>
        <li><a href="#TheclassTransfiniteInterpolationManifold">The class TransfiniteInterpolationManifold</a>
      <ul>
        <li><a href="#Howitworks">How it works</a>
      </ul>
        <li><a href="#Transfiniteinterpolationisexpensiveandhowtodealwithit">Transfinite interpolation is expensive and how to deal with it</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Analyticalsolutionandcoefficient">Analytical solution and coefficient</a>
        <li><a href="#ThePoissonProblemclass">The PoissonProblem class</a>
        <li><a href="#Gridcreationandinitializationofthemanifolds">Grid creation and initialization of the manifolds</a>
        <li><a href="#Setupofdatastructures">Setup of data structures</a>
        <li><a href="#Assemblyofthesystemmatrixandrighthandside">Assembly of the system matrix and right hand side</a>
        <li><a href="#Solutionofthelinearsystem">Solution of the linear system</a>
        <li><a href="#Outputofthesolutionandcomputationoferrors">Output of the solution and computation of errors</a>
        <li><a href="#ThePoissonProblemrunfunction">The PoissonProblem::run() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-65/doc/intro.dox



 <br> 

<i>
This program was contributed by Martin Kronbichler.
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序介绍了一个高级流形类--TransfiniteInterpolationManifold，以及如何解决其主要缺点--相对较高的成本。

<a name="Workingwithmanifolds"></a><h3>Working with manifolds</h3>


<a name="Whatwewant"></a><h4>What we want</h4>


在许多应用中，有限元网格必须能够表示一个相对复杂的几何体。在step-1、step-49和step-53教程程序中，已经介绍了deal.II库中可用的一些生成网格的技术。给定一个基础网格，deal.II能够通过将单元格均匀地或仅在计算域的选定部分细分为子单元来创建一个更细的网格。除了GridGenerator命名空间中的基本网格划分功能外，deal.II还提供了一些接口，可以使用命名空间GridIn的功能读入由（仅四边形和六边形）网格生成器生成的网格，例如在步骤5中演示。外部生成的网格的一个基本限制是，网格中生成的单元所提供的信息只包括顶点的位置和它们的连接性，而没有最初创建这个网格的网格生成器所提供的底层几何体的背景。一旦网格在deal.II中被细化，需要放置更多的点，这就成了问题。第54步教程程序显示了如何通过使用OpenCASCADE库方面的CAD曲面来克服这个限制，第53步通过在源代码内以编程方式提供同类信息。

在deal.II中，网格细化过程中新点的放置或高阶映射的定义由流形对象控制，详见 @ref manifold  "流形模块"。举个例子，考虑以下二维环形的情况（图片取自流形模块）。如果我们从10个单元的初始网格开始，在不附加任何流形的情况下全局细化网格三次，我们将得到以下网格。

 @image html hypershell-nothing.png "" 

图片看起来是这样的，因为在默认情况下，deal.II只知道通过平均父单元的顶点位置来放置子单元的顶点。这就产生了一个多边形域，其面是原始（粗略的网格）单元的边缘。很明显，我们必须对三角形的边界面附加一个曲线描述，以便在网格细化时重现圆形，就像下图一样。

 @image html hypershell-boundary-only.png "" 

这样就好多了。如果我们继续细化网格，至少现在的内外边界已经接近真实的圆了。然而，这幅图中的网格对于环形来说仍然不是最佳的，因为从一个单元到下一个单元的<i>interior</i>线在某些顶点上有扭结，人们宁愿使用下面的网格。

 @image html hypershell-all.png "" 

在这最后一种（最佳）情况下，也是由 GridGenerator::hyper_shell(), 产生的默认情况，曲线流形描述（在这种情况下是极地流形描述）不仅适用于边界面，而且适用于整个域。每当三角化要求一个新的点，例如，当它将一个单元细化为四个子单元时，边缘或单元的中点，它将沿着极坐标系统中各自的中点来放置它们。相比之下，上面那种只有边界受制于极地流形的情况，只有沿边界的中点会沿着曲线描述放置，而内部的中点会通过直角坐标系中周围点的合适平均值来计算（更多细节见 @ref manifold  "流形模块"）。

在这一点上，人们可能会认为曲线形的体积描述是一种方式。这一般来说是不会错的，尽管有时并不那么容易描述这到底应该如何工作。这里有几个例子。

- 想象一下，上面的网格实际上是一个圆盘，而不只是一个环。   在这种情况下，极地流形会在原点退化，不会产生合理的新点。事实上，为那些应该 "看起来是圆的 "但可能在原点或接近原点的东西定义流形描述是令人惊讶的非常困难的。

- 当人们试图将球形流形附加到整个体积上时，类似的事情也会发生在三维球的原点；在这种情况下，新流形点的计算会因异常而中止。

- CAD几何体通常只描述域的边界，就像我们在上面第二张图中只把流形附在边界上一样。同样地，Step-54只使用CAD几何体来生成表面网格（也许是因为这是解决相关问题所需要的），但是如果想解决那里描述的水或船周围空气中的问题，我们就需要有一个体积网格。那么问题来了，我们应该如何准确地描述域的内部应该发生什么。

这些简单的例子清楚地表明，对于许多有趣的情况，我们必须从对全体积的分析性曲率描述的愿望中退一步。将需要有<i>some</i>种导致曲率也在内部的信息，但必须有可能做到这一点，而不需要实际写下描述那种几何的明确公式。

那么，如果我们不在内部做任何事情，只把表面描述为流形，会发生什么呢？有时，如上图所示的环，结果并不可怕。但有时却很糟糕。考虑一个环状体的情况（例如用 GridGenerator::torus()) 生成的环状体，只在表面附加一个TorusManifold对象，内部单元和面没有附加流形，细化前有六个单元在环状方向。如果对网格进行一次细化，我们将得到如下的网格，图中显示的是网格的上半部分被剪切掉了。

 @image html torus_no_inner_manifold.png "" 

这显然是次优的。事实上，如果我们从少于上面显示的环形方向的六个单元开始，映射在某些区域实际上是倒置的，因为沿着内部单元放置的新点与边界相交，因为它们不是沿着环形方向的圆形。环形的简单情况仍然可以被固定下来，因为我们知道环形方向是沿着圆柱坐标系的，所以在曲面上附加一个TorusManifold，结合CylindricalManifold，在环形方向上有适当的周期性，应用于所有的内部实体，将产生一个高质量的网格，如下所示，现在有两个顶部的单元被隐藏起来。

 @image html torus_cylindrical_inner_manifold.png "" 

这个网格是相当不错的，但显然它与对体积的良好描述有关，而我们在其他情况下缺乏这种描述。实际上，在这种情况下也有一个不完美之处，因为我们可以看到域的内部的两个相邻单元的一些不自然的扭结，这些扭结被顶部的两个边界单元所隐藏，与下面的设置相反（由 GridGenerator::torus() 应用的默认流形并使用TransfiniteInterpolationManifold）。

 @image html torus_transfinite_manifold.png "" 

<a name="TheclassTransfiniteInterpolationManifold"></a><h3>The class TransfiniteInterpolationManifold</h3>


为了找到一个更好的策略，让我们再看看二维圆盘（这也是沿环形方向旋转的环形实体的基础）。正如我们上面所学到的，我们只能将弯曲的极坐标描述应用于边界（或离原点足够远的单元的边缘），但最终必须过渡到对圆盘中心的直线描述。如果我们在单元的内部使用平坦的流形（也就是说，新的顶点是通过相邻的现有顶点的平均化而产生的），而极地流形只用于圆盘的边界，那么在进行四次全局细化之后，我们会得到以下网格。

 @image html circular_mesh_only_boundary_manifold.png "" 

这并不是一个可怕的网格。同时，如果你知道最初的粗略网格是由中间的一个正方形和周围的四个盖子组成的，那么就不难看出这个网格的每一个细化步骤都是为了得到上面的图片。

虽然deal.II的三角测量类在创建新点时试图将信息从边界传播到内部，但这种算法的影响范围是有限的。

 @image html circular_mesh_boundary_cells.png "" 

上图突出了盘面上那些触及边界的单元，在当时只看一个单元的情况下，原则上可以考虑边界信息。显然，随着网格的细化，可以考虑到一些曲率的区域变得更加有限，因此在网格中形成了看似不规则的点。当计算最左边图片中任何一个边界单元的中心时，理想的位置是外圈和中间的单元之间的中点。这正是三角剖分类中第一个细化步骤所使用的内容。然而，对于第二次细化，所有的内部边缘以及内部单元格层只能根据平面流形描述来加点。

在这一点上，我们意识到需要什么来创建一个更好的网格。对于<i>all</i>子单元中的<i>any</i>新点，即在最左边的红色阴影层内创建的单元，我们要计算相对于各自粗单元覆盖区域内的曲率的插值。这可以通过在上图最左边面板的粗略网格的高亮单元中添加TransfiniteInterpolationManifold类来实现。该类遵守一般流形接口，即在其定义域内给定任何一组点，它可以计算符合流形的加权平均数（使用的公式将在一分钟后给出）。这些加权平均数在细化网格的时候，或者在给定单元上评估符合该流形的高阶映射（如MappingQGeneric或MappingC1）的时候都会用到。在圆盘粗网格的阴影单元上使用该流形（即不仅仅是最外层的单元），在经过四个全局步骤的细化后产生以下网格。

 @image html circular_mesh_transfinite_interpolation.png "" 

这个网格的线条仍有一些扭结，但它们仅限于粗大的网格单元之间的面，而网格的其他部分则如人们所希望的那样平滑。事实上，给定一个直边的中心单元，这个表示法是最好的，因为所有的网格单元都遵循一个平滑的过渡，从内部的方形块的直边到边界的圆形。(我们可以做得更好一些，在中心的方形块中也允许一些曲率，随着中心的接近，这些曲率最终会消失)。




<a name="Howitworks"></a><h4>How it works</h4>


在一个有一条曲线和三条直线的圆盘的简单情况下，我们可以明确写下如何实现形状的混合。为此，将物理单元（如顶部的单元）映射回参考坐标系 $(\xi,\eta)\in (0,1)^2$ 是很有用的，在那里我们计算某些点之间的平均数。如果我们使用一个由四个顶点 $(x_0,y_0), (x_1,y_1), (x_2,y_2), (x_3, y_3)$ 跨越的简单双线性地图，一个点 $(\xi, \eta)\in (0,1)^2$ 的图像将是

@f{align*}{
(x,y) = (1-\xi)(1-\eta) (x_0,y_0) + \xi(1-\eta) (x_1,y_1) +
       (1-\xi)\eta  (x_2,y_2) + \xi\eta  (x_3,y_3).


@f}



对于曲面的情况，我们要修改这个公式。对于圆盘粗略网格的顶部单元，我们可以假设点 $(x_0,y_0)$ 和 $(x_1,y_1)$ 位于下端的直线上，点 $(x_2,y_2)$ 和 $(x_3,y_3)$ 沿顶部由一个四分之一圆连接。然后我们将点 $(\xi, \eta)$ 映射为

@f{align*}{
(x,y) = (1-\eta) \big[(1-\xi) (x_0,y_0) + \xi (x_1,y_1)\big] +
      \eta \mathbf{c}_3(\xi),


@f}

其中 $\mathbf{c}_3(\xi)$ 是一条曲线，用一个arclength参数 $\xi\in (0,1)$ 来描述四分之一圆的 $(x,y)$ 坐标。这代表了单元格的直线下边缘和曲线上边缘之间的线性插值，也是上图的基础。

这个公式很容易被推广到所有四条边都由曲线而不是直线描述的情况。我们称这四个函数为 $\xi$ 或 $\eta$ ，在水平和垂直方向上，分别为四边形的左、右、下和上边缘的参数。然后，内插法为

@f{align*}{
(x,y) =& (1-\xi)\mathbf{c}_0(\eta) + \xi \mathbf{c}_1(\eta)
        +(1-\eta)\mathbf{c}_2(\xi) + \eta \mathbf{c}_3(\xi)\\
       &-\big[(1-\xi)(1-\eta) (x_0,y_0) + \xi(1-\eta) (x_1,y_1) +
        (1-\xi)\eta  (x_2,y_2) + \xi\eta  (x_3,y_3)\big].


@f}



这个公式假设边界曲线与顶点 $(x_0,y_0), (x_1,y_1), (x_2,y_2), (x_3, y_3)$ ，例如 $\mathbf{c}_0(0)
= (x_0,y_0)$ 或 $\mathbf{c}_0(1) = (x_2,y_2)$ 匹配和重合。公式第二行的双线性插值的减法确保了在边界上完全遵循规定的曲线：沿着四条边中的每一条，我们需要减去在角上评估的两条相邻边的贡献，这时的贡献只是一个顶点位置。很容易检查出，如果四条曲线中的三条 $\mathbf{c}_i$ 是直的，从而与双线性插值重合，那么上面的圆的公式就再现了。

这个公式被称为转折内插，由<a
href="https://doi.org/10.1002%2Fnme.1620070405">Gordon and Hall</a>在1973年提出。尽管转义插值本质上只表示边界曲线的线性混合，但插值完全遵循每个实数 $\xi\in (0,1)$ 或 $\eta\in (0,1)$ 的边界曲线，也就是说，它插值的点数是无限的，这也是Gordon和Hall将这种插值的变体称为转义插值的最初动机。另一种解释是，无限插值从左右和上下线性插值，我们需要从中减去双线性插值，以确保在域的内部有一个单位重量。

无限插值很容易被推广到三个空间维度。在这种情况下，插值允许为一个三维单元的任何一个四边形混合6种不同的表面描述，为一个单元的线混合12种边缘描述。同样，为了确保地图的一致性，有必要减去边缘的贡献，再加上顶点的贡献，使曲线遵循规定的表面或边缘描述。在三维的情况下，也可以使用从曲线边缘到相邻面和相邻单元的转折插值。

在处理.II中的转置插值是通用的，因为它可以处理任意的曲线。它将以 $d$ -维空间的原始坐标来评估曲线，但有一个（或两个，在三维的边缘情况下）坐标固定在 $0$ 或 $1$ ，以确保任何其他流形类，包括CAD文件（如果需要），都可以应用于开箱。无限插值是网格生成器的标准成分，因此在deal.II库中集成这一功能的主要优势是在自适应细化和粗化网格时启用，以及用于创建更高程度的映射，使用流形来插入网格顶点以外的额外点。

关于无限插值的最后一句话，我们提到在没有体积流形描述的情况下，deal.II中的网格细化策略也是基于无限插值的权重，在这个意义上是最优的。不同的是，默认算法一次只看到一个单元，因此将只对那些接触曲面流形的单元应用最优算法。相反，在整个<i>patches</i>个单元上使用转义映射（源于一个较粗的单元），可以用转义插值法，将信息从边界传播到远处的单元。




<a name="Transfiniteinterpolationisexpensiveandhowtodealwithit"></a><h3>Transfinite interpolation is expensive and how to deal with it</h3>


一个具有无限流形描述的网格通常分两步建立。第一步是创建一个粗略的网格（或者从文件中读入），并在一些网格实体上附加一个曲面流形。对于上面的圆盘例子，我们将极坐标流形附加到沿外圆的面（这是由 GridGenerator::hyper_ball()). 自动完成的）在我们开始细化网格之前，我们再为网格的所有内部单元和边缘分配一个TransfiniteInterpolationManifold，当然，这需要基于我们分配给这些实体的一些流形ID（除了边界上的圆以外的一切）。我们是否也给圆盘的内部正方形分配一个TransfiniteInterpolationManifold并不重要，因为对具有直边（或3D中的平坦面）的粗单元进行Transfinite插值，只是产生具有直边（平坦面）的细分子。

之后，当网格被细化或基于此网格建立高阶映射时，单元将查询底层流形对象的新点。这个过程需要一组周围的点，例如一个二维单元的四个顶点，以及一组对这些点的权重，用于定义一个新的点。对于一个单元格的中间点，四个顶点中的每一个都将得到0.25的权重。对于无限插值流形，建立加权和的过程需要一些严肃的工作。根据结构，我们要根据周围点的参考坐标 $\xi$ 和 $\eta$ （或三维的 $\xi, \eta, \zeta$ ）来组合这些点。然而，deal.II中流形类的接口并没有得到周围点的参考坐标（因为它们不是全局存储的），而只是物理坐标。因此，无限插值流形必须做的第一步是反转映射，并在无限插值的一个粗网格单元内找到参考坐标（例如，上述磁盘网格的四个阴影粗网格单元之一）。这种反演是通过牛顿迭代（或者说，基于有限差分的牛顿方案与布罗伊登方法相结合）完成的，并根据上述公式多次查询无限期插值。这些查询中的每一次都可能依次调用一个昂贵的流形，例如球的球形描述，并且本身就很昂贵。由于deal.II的Manifold接口类只提供了一组点，所以转置插值最初甚至不知道周围的点集属于哪个粗大的网格单元，需要根据一些启发式方法在几个单元中搜索。就<a
href="https://en.wikipedia.org/wiki/Atlas_(topology)#Charts">charts</a>而言，我们可以把无限插值的实现描述为基于<a
href="https://en.wikipedia.org/wiki/Atlas_(topology)">atlas</a>的实现。三角形初始粗网格的每个单元都代表一个有自己参考空间的图表，而周围的流形提供了一种从图表空间（即参考单元）到物理空间的转换方式。粗网格单元的图表的集合是一个图集，像往常一样，在图集中查找东西时，人们做的第一件事就是找到正确的图表。

一旦找到周围点的参考坐标，参考坐标系中的新点将通过简单的加权和计算出来。最后，参考点被插入到无限插值的公式中，这就得到了所需的新点。

在许多情况下，曲面流形不仅在网格细化过程中使用，而且还用于确保计算域单元内边界的曲面表示。这对于保证复杂几何上的高阶多项式的高阶收敛性是必须的，但有时对于线性形状函数也需要一个精确的几何。这通常是通过对单元的多项式描述来实现的，如果表示曲面网格元素的多项式程度与数值解的多项式程度相同，则称为等参数概念。如果几何体的度数高于或低于解的度数，则分别称为超参数或次参数的几何体表示。在deal.II中，多项式表示的标准类是MappingQGeneric。例如，如果在三维中使用这个类的多项式度数 $4$ ，总共需要125个（即 $(4+1)^3$ ）点来进行内插。在这些点中，8个是单元的顶点，已经可以从网格中获得，但其他117个点需要由流形提供。如果使用无限插值流形，我们可以想象，通过回拉到一些尚未确定的粗略单元的参考坐标，然后在117个点中的每个点上进行后续的前推，是一个大量的工作，可能非常耗时。

更糟糕的是，许多程序的结构是这样的：对于同一个单元，映射被独立地查询了几次。它的主要用途是组装线性系统，即通过FEValues对象的`mapping`参数计算系统矩阵和右手边。然而，边界值的插值、数值误差的计算、输出的写入以及误差估计器的评估也必须涉及相同的映射，以确保对解向量的一致解释。因此，即使是解决一次的线性静止问题，也会多次评估映射的点。对于上面提到的三维案例，这意味着通过昂贵的算法多次计算每个单元的117个点。对于非线性或随时间变化的问题，情况更加紧迫，这些操作要反复进行。

由于通过无限插值进行流形描述的成本比对平流形的类似查询要高几百倍，所以只计算一次额外的点并在所有后续调用中使用这些点是有意义的。deal.II库提供的MappingQCache类正是为了这个目的。与系统矩阵所消耗的内存相比，这个缓存通常不会太大，在看这个教程程序的结果时就会明白。MappingQCache的用法很简单。一旦网格被设置好（或在细化过程中被改变），我们就调用 MappingQCache::initialize() ，将所需的三角形以及所需的映射作为参数。然后，初始化会遍历网格的所有单元，并查询给定映射的附加点。这些点会被保存为单元的标识符，以便以后当映射计算一些与单元相关的量（如参考坐标和物理坐标之间的映射的Jacobian）时，它们可以被返回。

最后，我们提到，TransfiniteInterpolationManifold也使得网格的细化变得更加昂贵。在这种情况下，MappingQCache无济于事，因为它将计算随后不能再使用的点；目前在deal.II中不存在一个更有效的机制。然而，网格细化也包含许多其他昂贵的步骤，所以与其他计算相比，它并不是一个大问题。它在每个时间步长或非线性迭代中也最多只发生一次。

<a name="Thetestcase"></a><h3>The test case</h3>


在这个教程程序中，TransfiniteInterpolationManifold与MappingQCache的结合被举例说明。这个测试案例相对简单，占用了许多典型程序中涉及的解决阶段，例如步骤6的教程程序。作为一个几何体，我们选择了TransfiniteInterpolationManifold的一个使用原型，即涉及一个球形的设置，而这个球形又被一个立方体所包围。例如，这样的设置将用于嵌入背景介质中的球形包容物，如果该包容物具有不同的材料特性，需要通过元素界面来跟踪两种材料之间的界面。这里给出了一个可视化的网格。

 <img src="https://www.dealii.org/images/steps/developer/step-65-mesh.png" alt=""> 

对于这种情况，我们要在域内的表面附加一个球形描述，并使用转折插值来平滑地切换到外部立方体的直线和球中心的立方体。

在该程序中，我们将遵循有限元程序的典型流程，从DoFHandler和稀疏模式的设置开始，组装一个线性系统来解决带有跳跃系数的泊松方程，用一个简单的迭代方法解决，用 VectorTools::integrate_difference() 计算一些数值误差以及误差估计器。我们记录了每个部分的时间，并运行了两次代码。在第一次运行中，我们将一个MappingQGeneric对象分别交给程序的每个阶段，在这些阶段中，点被反复计算。在第二次运行中，我们用MappingQCache代替。


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
 * 本教程的包含文件与  step-6  中的基本相同。重要的是，我们将使用的TransfiniteInterpolationManifold类是由`deal.II/grid/manifold_lib.h`提供。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/timer.h> 
 * 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/vector.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/mapping_q_generic.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * #include <fstream> 
 * 
 * @endcode
 * 
 * 唯一的新include文件是MappingQCache类的文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/mapping_q_cache.h> 
 * 
 * namespace Step65 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Analyticalsolutionandcoefficient"></a> 
 * <h3>Analytical solution and coefficient</h3>
 * 

 * 
 * 在这个教程程序中，我们要解决泊松方程，其系数沿半径为0.5的球体跳跃，并使用一个恒定的右手边值 $f(\mathbf{x}) = -3$  。（这个设置与 step-5 和 step-6 相似，但系数和右手边的具体数值不同）。由于系数的跳跃，分析解必须有一个结点，即系数从一个值切换到另一个值。为了保持简单，我们选择了一个在所有分量中都是二次的分析解，即在半径为0.5的球中为 $u(x,y,z) = x^2 + y^2 + z^2$ ，在域的外部为 $u(x,y,z) = 0.1(x^2 + y^2 + z^2) + 0.25-0.025$ 。这个分析解在内球的系数为0.5，外球的系数为5的情况下与右手边兼容。它也是沿着半径为0.5的圆连续的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ExactSolution : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> &p, 
 *                          const unsigned int /*component*/ = 0) const override 
 *     { 
 *       if (p.norm_square() < 0.25) 
 *         return p.norm_square(); 
 *       else 
 *         return 0.1 * p.norm_square() + (0.25 - 0.025); 
 *     } 
 * 
 *     virtual Tensor<1, dim> 
 *     gradient(const Point<dim> &p, 
 *              const unsigned int /*component*/ = 0) const override 
 *     { 
 *       if (p.norm_square() < 0.25) 
 *         return 2. * p; 
 *       else 
 *         return 0.2 * p; 
 *     } 
 *   }; 
 * 
 *   template <int dim> 
 *   double coefficient(const Point<dim> &p) 
 *   { 
 *     if (p.norm_square() < 0.25) 
 *       return 0.5; 
 *     else 
 *       return 5.0; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePoissonProblemclass"></a> 
 * <h3>The PoissonProblem class</h3>
 * 

 * 
 * 泊松问题的实现与我们在  step-5  教程中使用的非常相似。两个主要的区别是，我们向程序中的各个步骤传递了一个映射对象，以便在两种映射表示法之间进行切换，正如介绍中所解释的那样，还有一个`计时器'对象（TimerOutput类型），将用于测量各种情况下的运行时间。(映射对象的概念在 step-10 和 step-11 中首次提出，如果你想查一下这些类的用途的话)。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class PoissonProblem 
 *   { 
 *   public: 
 *     PoissonProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void create_grid(); 
 *     void setup_system(const Mapping<dim> &mapping); 
 *     void assemble_system(const Mapping<dim> &mapping); 
 *     void solve(); 
 *     void postprocess(const Mapping<dim> &mapping); 
 * 
 *     Triangulation<dim> triangulation; 
 *     FE_Q<dim>          fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     AffineConstraints<double> constraints; 
 *     SparsityPattern           sparsity_pattern; 
 *     SparseMatrix<double>      system_matrix; 
 *     Vector<double>            solution; 
 *     Vector<double>            system_rhs; 
 * 
 *     TimerOutput timer; 
 *   }; 
 * 
 * @endcode
 * 
 * 在构造函数中，我们设置了定时器对象来记录墙的时间，但在正常执行过程中是安静的。我们将在 `PoissonProblem::run()` 函数中查询它的计时细节。此外，我们为正在使用的有限元选择了一个相对较高的多项式三度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   PoissonProblem<dim>::PoissonProblem() 
 *     : fe(3) 
 *     , dof_handler(triangulation) 
 *     , timer(std::cout, TimerOutput::never, TimerOutput::wall_times) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="Gridcreationandinitializationofthemanifolds"></a> 
 * <h3>Grid creation and initialization of the manifolds</h3>
 * 

 * 
 * 接下来的函数介绍了TransfiniteInterpolationManifold的典型用法。第一步是创建所需的网格，这可以通过GridGenerator的两个网格的组合来完成。内球网格是很简单的。我们以原点为中心运行 GridGenerator::hyper_cube() ，半径为0.5（第三个函数参数）。第二个网格更有趣，构建方法如下。我们希望有一个在内部是球形的，但在外表面是平的网格。此外，内球的网格拓扑结构应该与外球的网格兼容，即它们的顶点重合，这样才能使两个网格合并起来。从 GridGenerator::hyper_shell 出来的网格满足了内侧的要求，如果它是用 $2d$ 的粗大单元创建的（在3D中我们将使用6个粗大单元）&ndash；这与球的边界面的单元数量相同。对于外表面，我们利用这样一个事实：没有流形附着的壳表面的6个面将退化为立方体的表面。我们仍然缺少的是外壳边界的半径。由于我们想要一个范围为 $[-1, 1]$ 的立方体，而6单元壳将其8个外顶点放在8条对角线上，我们必须将点 $(\pm 1, \pm 1, \pm 1)$ 转化为半径。显然，在 $d$ 维度上，半径必须是 $\sqrt{d}$ ，也就是说，对于我们要考虑的三维情况，半径是 $\sqrt{3}$ 。
 * 

 * 
 * 这样，我们就有了一个计划。在创建了球的内部三角形和外壳的三角形之后，我们将这两个网格合并，但是将GridGenerator中的函数可能从产生的三角形中设置的所有流形移除，以确保我们对流形有充分的控制。特别是，我们希望在细化过程中在边界上添加的额外点能够遵循平坦的流形描述。为了开始添加更合适的流形ID的过程，我们给所有的网格实体（单元、面、线）分配流形ID 0，这些实体以后将与TransfiniteInterpolationManifold相关联。然后，我们必须识别沿着半径为0.5的球体的面和线，并给它们标记一个不同的流形ID，以便随后给这些面和线分配一个SphericalManifold。由于我们在调用 GridGenerator::hyper_ball(), 后丢弃了所有预先存在的流形，我们手动检查了网格的单元格和所有的面。如果四个顶点的半径都是0.5，我们就在球体上找到了一个面，或者像我们在程序中写的那样，有  $r^2-0.25 \approx 0$  。注意，我们调用`cell->face(f)->set_all_manifold_ids(1)`来设置面和周围线上的流形id。此外，我们希望通过一个材料ID来区分球内和球外的单元，以便于可视化，对应于介绍中的图片。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PoissonProblem<dim>::create_grid() 
 *   { 
 *     Triangulation<dim> tria_inner; 
 *     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5); 
 * 
 *     Triangulation<dim> tria_outer; 
 *     GridGenerator::hyper_shell( 
 *       tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim); 
 * 
 *     GridGenerator::merge_triangulations(tria_inner, tria_outer, triangulation); 
 * 
 *     triangulation.reset_all_manifolds(); 
 *     triangulation.set_all_manifold_ids(0); 
 * 
 *     for (const auto &cell : triangulation.cell_iterators()) 
 *       { 
 *         for (const auto &face : cell->face_iterators()) 
 *           { 
 *             bool face_at_sphere_boundary = true; 
 *             for (const auto v : face->vertex_indices()) 
 *               { 
 *                 if (std::abs(face->vertex(v).norm_square() - 0.25) > 1e-12) 
 *                   face_at_sphere_boundary = false; 
 *               } 
 *             if (face_at_sphere_boundary) 
 *               face->set_all_manifold_ids(1); 
 *           } 
 *         if (cell->center().norm_square() < 0.25) 
 *           cell->set_material_id(1); 
 *         else 
 *           cell->set_material_id(0); 
 *       } 
 * 
 * @endcode
 * 
 * 有了所有单元格、面和线的适当标记，我们可以将流形对象附加到这些数字上。流形ID为1的实体将得到一个球形流形，而流形ID为0的其他实体将被分配到TransfiniteInterpolationManifold。正如介绍中提到的，我们必须通过调用 TransfiniteInterpolationManifold::initialize() 显式初始化当前网格的流形，以获取粗略的网格单元和连接到这些单元边界的流形。我们还注意到，我们在这个函数中本地创建的流形对象是允许超出范围的（就像它们在函数范围结束时那样），因为Triangulation对象在内部复制它们。
 * 

 * 
 * 在连接了所有的流形之后，我们最后将去细化网格几次，以创建一个足够大的测试案例。
 * 

 * 
 * 
 * @code
 *     triangulation.set_manifold(1, SphericalManifold<dim>()); 
 * 
 *     TransfiniteInterpolationManifold<dim> transfinite_manifold; 
 *     transfinite_manifold.initialize(triangulation); 
 *     triangulation.set_manifold(0, transfinite_manifold); 
 * 
 *     triangulation.refine_global(9 - 2 * dim); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Setupofdatastructures"></a> 
 * <h3>Setup of data structures</h3>
 * 

 * 
 * 下面的函数在其他教程中是众所周知的，它枚举了自由度，创建了一个约束对象并为线性系统设置了一个稀疏矩阵。唯一值得一提的是，该函数接收了一个映射对象的引用，然后我们将其传递给 VectorTools::interpolate_boundary_values() 函数，以确保我们的边界值在用于装配的高阶网格上被评估。在本例中，这并不重要，因为外表面是平的，但对于弯曲的外单元，这将导致边界值的更精确的近似。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PoissonProblem<dim>::setup_system(const Mapping<dim> &mapping) 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 *     std::cout << "   Number of active cells:       " 
 *               << triangulation.n_global_active_cells() << std::endl; 
 *     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 * 
 *     { 
 *       TimerOutput::Scope scope(timer, "Compute constraints"); 
 * 
 *       constraints.clear(); 
 * 
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *       VectorTools::interpolate_boundary_values( 
 *         mapping, dof_handler, 0, ExactSolution<dim>(), constraints); 
 * 
 *       constraints.close(); 
 *     } 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
 * 
 *     sparsity_pattern.copy_from(dsp); 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Assemblyofthesystemmatrixandrighthandside"></a> 
 * <h3>Assembly of the system matrix and right hand side</h3>
 * 

 * 
 * 组装线性系统的函数在前面的教程程序中也是众所周知的。有一点需要注意的是，我们将正交点的数量设置为多项式的度数加2，而不是像其他大多数教程中的度数加1。这是因为我们期望有一些额外的精度，因为映射也涉及到比解的多项式多一度的程度。
 * 

 * 
 * 汇编中唯一有点不寻常的代码是我们计算单元格矩阵的方式。我们没有使用正交点索引、行和矩阵列的三个嵌套循环，而是首先收集形状函数的导数，乘以系数和积分因子`JxW`的乘积的平方根，放在一个单独的矩阵`partial_matrix`中。为了计算单元矩阵，我们在 "partial_matrix.mTmult(cell_matrix, partial_matrix); "一行中执行 "cell_matrix = partial_matrix * transpose(partial_matrix)"。为了理解这一点，我们要知道矩阵与矩阵的乘法是对`partial_matrix`的各列进行求和。如果我们用 
 * $a(\mathbf{x}_q)$ 表示系数，临时矩阵的条目是 $\sqrt{\text{det}(J) w_q a(x)} \frac{\partial \varphi_i(\boldsymbol
 * \xi_q)}{\partial x_k}$ 。如果我们将该矩阵的第<i>i</i>行与第<i>j</i>列相乘，我们计算出一个涉及 $\sum_q \sum_{k=1}^d \sqrt{\text{det}(J) w_q a(x)} \frac{\partial
 * \varphi_i(\boldsymbol \xi_q)}{\partial x_k} \sqrt{\text{det}(J) w_q a(x)}
 * \frac{\partial \varphi_j(\boldsymbol \xi_q)}{\partial x_k} = \sum_q
 * \sum_{k=1}^d\text{det}(J) w_q a(x)\frac{\partial \varphi_i(\boldsymbol
 * \xi_q)}{\partial x_k} \frac{\partial \varphi_j(\boldsymbol
 * \xi_q)}{\partial x_k}$ 的嵌套和，这正是拉普拉斯方程的双线性形式所需的条款。
 * 

 * 
 * 选择这种有点不寻常的方案的原因是由于计算三维中相对较高的多项式程度的单元矩阵所涉及的繁重工作。由于我们想在这个教程程序中强调映射的成本，我们最好以优化的方式进行装配，以便不追逐已经被社区解决的瓶颈。矩阵-矩阵乘法是HPC背景下最好的优化内核之一， FullMatrix::mTmult() 函数将调用到那些优化的BLAS函数。如果用户在配置deal.II时提供了一个好的BLAS库（如OpenBLAS或英特尔的MKL），那么单元矩阵的计算将执行到接近处理器的峰值算术性能。顺便提一下，尽管有优化的矩阵-矩阵乘法，但目前的策略在复杂性方面是次优的，因为要做的工作与 $(p+1)^9$ 度 $p$ 的运算成正比（这也适用于用FEValues的通常评估）。我们可以通过利用形状函数的张量乘积结构，用 $\mathcal O((p+1)^7)$ 的操作来计算单元格矩阵，就像交易二中的无矩阵框架那样。我们参考 step-37 和张量积感知评估器FEEvaluation的文档，以了解如何实现更有效的单元矩阵计算的细节。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PoissonProblem<dim>::assemble_system(const Mapping<dim> &mapping) 
 *   { 
 *     TimerOutput::Scope scope(timer, "Assemble linear system"); 
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 2); 
 *     FEValues<dim>     fe_values(mapping, 
 *                             fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 *     FullMatrix<double> partial_matrix(dofs_per_cell, dim * n_q_points); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_rhs = 0.; 
 *         fe_values.reinit(cell); 
 * 
 *         for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) 
 *           { 
 *             const double current_coefficient = 
 *               coefficient(fe_values.quadrature_point(q_index)); 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               { 
 *                 for (unsigned int d = 0; d < dim; ++d) 
 *                   partial_matrix(i, q_index * dim + d) = 
 *                     std::sqrt(fe_values.JxW(q_index) * current_coefficient) * 
 *                     fe_values.shape_grad(i, q_index)[d]; 
 *                 cell_rhs(i) += 
 *                   (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
 *                    (-dim) *                            // f(x_q) 
 *                    fe_values.JxW(q_index));            // dx 
 *               } 
 *           } 
 * 
 *         partial_matrix.mTmult(cell_matrix, partial_matrix); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         constraints.distribute_local_to_global( 
 *           cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Solutionofthelinearsystem"></a> 
 * <h3>Solution of the linear system</h3>
 * 

 * 
 * 对于线性系统的求解，我们选择一个简单的雅可比条件共轭梯度求解器，类似于早期教程中的设置。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PoissonProblem<dim>::solve() 
 *   { 
 *     TimerOutput::Scope scope(timer, "Solve linear system"); 
 * 
 *     SolverControl            solver_control(1000, 1e-12); 
 *     SolverCG<Vector<double>> solver(solver_control); 
 * 
 *     PreconditionJacobi<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix); 
 * 
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner); 
 *     constraints.distribute(solution); 
 * 
 *     std::cout << "   Number of solver iterations:  " 
 *               << solver_control.last_step() << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Outputofthesolutionandcomputationoferrors"></a> 
 * <h3>Output of the solution and computation of errors</h3>
 * 

 * 
 * 在下一个函数中，我们对解决方案做了各种后处理步骤，所有这些步骤都以这种或那种方式涉及映射。
 * 

 * 
 * 我们做的第一个操作是把解决方案以及材料ID写到VTU文件中。这与其他许多教程程序中的做法类似。这个教程程序中提出的新内容是，我们要确保写到文件中用于可视化的数据实际上是deal.II内部使用的数据的忠实代表。这是因为大多数可视化数据格式只用顶点坐标表示单元，但没有办法表示deal.II中使用高阶映射时的曲线边界--换句话说，你在可视化工具中看到的东西实际上不是你正在计算的东西。顺带一提，在使用高阶形状函数时也是如此。大多数可视化工具只呈现双线性/三线性的表示。这在 DataOut::build_patches().) 中有详细的讨论。
 * 

 * 
 * 所以我们需要确保高阶表示被写入文件中。我们需要考虑两个特别的话题。首先，我们通过 DataOutBase::VtkFlags 告诉DataOut对象，我们打算把元素的细分解释为高阶拉格朗日多项式，而不是双线性补丁的集合。最近的可视化程序，如ParaView 5.5版或更新的程序，然后可以呈现高阶解决方案（更多细节见<a
 * href="https:github.com/dealii/dealii/wiki/Notes-on-visualizing-high-order-output">wiki
 * page</a>）。其次，我们需要确保映射被传递给 DataOut::build_patches() 方法。最后，DataOut类默认只打印<i>boundary</i>单元的曲面，所以我们需要确保通过映射将内部单元也打印成曲面。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PoissonProblem<dim>::postprocess(const Mapping<dim> &mapping) 
 *   { 
 *     { 
 *       TimerOutput::Scope scope(timer, "Write output"); 
 * 
 *       DataOut<dim> data_out; 
 * 
 *       DataOutBase::VtkFlags flags; 
 *       flags.write_higher_order_cells = true; 
 *       data_out.set_flags(flags); 
 * 
 *       data_out.attach_dof_handler(dof_handler); 
 *       data_out.add_data_vector(solution, "solution"); 
 * 
 *       Vector<double> material_ids(triangulation.n_active_cells()); 
 *       for (const auto &cell : triangulation.active_cell_iterators()) 
 *         material_ids[cell->active_cell_index()] = cell->material_id(); 
 *       data_out.add_data_vector(material_ids, "material_ids"); 
 * 
 *       data_out.build_patches(mapping, 
 *                              fe.degree, 
 *                              DataOut<dim>::curved_inner_cells); 
 * 
 *       std::ofstream file( 
 *         ("solution-" + 
 *          std::to_string(triangulation.n_global_levels() - 10 + 2 * dim) + 
 *          ".vtu") 
 *           .c_str()); 
 * 
 *       data_out.write_vtu(file); 
 *     } 
 * 
 * @endcode
 * 
 * 后处理函数的下一个操作是对照分析解计算 $L_2$ 和 $H^1$ 误差。由于分析解是一个二次多项式，我们期望在这一点上得到一个非常准确的结果。如果我们是在一个具有平面面的简单网格上求解，并且系数的跳动与单元间的面对齐，那么我们会期望数值结果与分析解相吻合，直至舍去精度。然而，由于我们使用的是跟随球体的变形单元，这些单元只能由4度的多项式跟踪（比有限元的度数多一个），我们会发现在 $10^{-7}$ 附近有一个误差。我们可以通过增加多项式的度数或细化网格来获得更多的精度。
 * 

 * 
 * 
 * @code
 *     { 
 *       TimerOutput::Scope scope(timer, "Compute error norms"); 
 * 
 *       Vector<double> norm_per_cell_p(triangulation.n_active_cells()); 
 * 
 *       VectorTools::integrate_difference(mapping, 
 *                                         dof_handler, 
 *                                         solution, 
 *                                         ExactSolution<dim>(), 
 *                                         norm_per_cell_p, 
 *                                         QGauss<dim>(fe.degree + 2), 
 *                                         VectorTools::L2_norm); 
 *       std::cout << "   L2 error vs exact solution:   " 
 *                 << norm_per_cell_p.l2_norm() << std::endl; 
 * 
 *       VectorTools::integrate_difference(mapping, 
 *                                         dof_handler, 
 *                                         solution, 
 *                                         ExactSolution<dim>(), 
 *                                         norm_per_cell_p, 
 *                                         QGauss<dim>(fe.degree + 2), 
 *                                         VectorTools::H1_norm); 
 *       std::cout << "   H1 error vs exact solution:   " 
 *                 << norm_per_cell_p.l2_norm() << std::endl; 
 *     } 
 * 
 * @endcode
 * 
 * 我们在这里做的最后一个后处理操作是用KellyErrorEstimator计算出一个误差估计。我们使用了与 step-6 教程程序中完全相同的设置，只是我们还交出了映射，以确保误差是沿着曲线元素评估的，与程序的其余部分一致。然而，我们并没有真正使用这里的结果来驱动网格适应步骤（会沿着球体细化材料界面周围的网格），因为这里的重点是这个操作的成本。
 * 

 * 
 * 
 * @code
 *     { 
 *       TimerOutput::Scope scope(timer, "Compute error estimator"); 
 * 
 *       Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 *       KellyErrorEstimator<dim>::estimate( 
 *         mapping, 
 *         dof_handler, 
 *         QGauss<dim - 1>(fe.degree + 1), 
 *         std::map<types::boundary_id, const Function<dim> *>(), 
 *         solution, 
 *         estimated_error_per_cell); 
 *       std::cout << "   Max cell-wise error estimate: " 
 *                 << estimated_error_per_cell.linfty_norm() << std::endl; 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePoissonProblemrunfunction"></a> 
 * <h3>The PoissonProblem::run() function</h3>
 * 

 * 
 * 最后，我们定义了`run()`函数，控制我们如何执行这个程序（由main()函数以常规方式调用）。我们首先调用`create_grid()`函数，用适当的流形设置我们的几何体。然后我们运行两个求解器链的实例，从方程的设置开始，组装线性系统，用一个简单的迭代求解器求解，以及上面讨论的后处理。这两个实例在使用映射的方式上有所不同。第一个使用传统的MappingQGeneric映射对象，我们将其初始化为比有限元多一级的程度；毕竟，我们期望几何表示是瓶颈，因为分析解只是二次多项式。实际上，事情在相当程度上是相互关联的，因为实坐标中多项式的评估涉及到高阶多项式的映射，而高阶多项式代表一些光滑的有理函数。因此，高阶多项式还是有回报的，所以进一步增加映射的度数是没有意义的)。一旦第一遍完成，我们就让定时器打印出各个阶段的计算时间的摘要。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PoissonProblem<dim>::run() 
 *   { 
 *     create_grid(); 
 * 
 *     { 
 *       std::cout << std::endl 
 *                 << "====== Running with the basic MappingQGeneric class ====== " 
 *                 << std::endl 
 *                 << std::endl; 
 * 
 *       MappingQGeneric<dim> mapping(fe.degree + 1); 
 *       setup_system(mapping); 
 *       assemble_system(mapping); 
 *       solve(); 
 *       postprocess(mapping); 
 * 
 *       timer.print_summary(); 
 *       timer.reset(); 
 *     } 
 * 
 * @endcode
 * 
 * 对于第二个实例，我们转而设置了MappingQCache类。它的使用非常简单。在构建好它之后（考虑到我们希望它在其他情况下显示正确的度数功能，所以用度数），我们通过 MappingQCache::initialize() 函数填充缓存。在这个阶段，我们为缓存指定我们想要使用的映射（很明显，与之前的MappingQGeneric相同，以便重复相同的计算），然后再次运行相同的函数，现在交出修改后的映射。最后，我们再次打印重置后的累计壁挂时间，看看这些时间与原来的设置相比如何。
 * 

 * 
 * 
 * @code
 *     { 
 *       std::cout 
 *         << "====== Running with the optimized MappingQCache class ====== " 
 *         << std::endl 
 *         << std::endl; 
 * 
 *       MappingQCache<dim> mapping(fe.degree + 1); 
 *       { 
 *         TimerOutput::Scope scope(timer, "Initialize mapping cache"); 
 *         mapping.initialize(MappingQGeneric<dim>(fe.degree + 1), triangulation); 
 *       } 
 *       std::cout << "   Memory consumption cache:     " 
 *                 << 1e-6 * mapping.memory_consumption() << " MB" << std::endl; 
 * 
 *       setup_system(mapping); 
 *       assemble_system(mapping); 
 *       solve(); 
 *       postprocess(mapping); 
 * 
 *       timer.print_summary(); 
 *     } 
 *   } 
 * } // namespace Step65 
 * 
 * int main() 
 * { 
 *   Step65::PoissonProblem<3> test_program; 
 *   test_program.run(); 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-65/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Programoutput"></a><h3>Program output</h3>


如果我们用三度的多项式运行这个程序的三维版本，我们会得到以下程序输出。

@code
> make run
Scanning dependencies of target \step-65
[ 33%] Building CXX object CMakeFiles/\step-65.dir/\step-65.cc.o
[ 66%] Linking CXX executable \step-65
[ 66%] Built target \step-65
[100%] Run \step-65 with Release configuration


====== Running with the basic MappingQGeneric class ======


   Number of active cells:       6656
   Number of degrees of freedom: 181609
   Number of solver iterations:  285
   L2 error vs exact solution:   8.99339e-08
   H1 error vs exact solution:   6.45341e-06
   Max cell-wise error estimate: 0.00743406



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      49.4s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble linear system          |         1 |       5.8s |        12% |
| Compute constraints             |         1 |     0.109s |      0.22% |
| Compute error estimator         |         1 |      16.5s |        33% |
| Compute error norms             |         1 |      9.11s |        18% |
| Solve linear system             |         1 |      9.92s |        20% |
| Write output                    |         1 |      4.85s |       9.8% |
+---------------------------------+-----------+------------+------------+


====== Running with the optimized MappingQCache class ======


   Memory consumption cache:     22.9981 MB
   Number of active cells:       6656
   Number of degrees of freedom: 181609
   Number of solver iterations:  285
   L2 error vs exact solution:   8.99339e-08
   H1 error vs exact solution:   6.45341e-06
   Max cell-wise error estimate: 0.00743406



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      18.4s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble linear system          |         1 |      1.44s |       7.8% |
| Compute constraints             |         1 |   0.00336s |         0% |
| Compute error estimator         |         1 |     0.476s |       2.6% |
| Compute error norms             |         1 |     0.505s |       2.7% |
| Initialize mapping cache        |         1 |      4.96s |        27% |
| Solve linear system             |         1 |      9.95s |        54% |
| Write output                    |         1 |     0.875s |       4.8% |
+---------------------------------+-----------+------------+------------+


[100%] Built target run
@endcode



在讨论时序之前，我们看一下MappingQCache对象的内存消耗。我们的程序打印出它占用了23MB的内存。如果我们把这个数字与单个（解或右侧）向量的内存消耗联系起来，即1.5MB（即181,609个元素乘以每条8字节的双精度），或者与系统矩阵和稀疏模式所消耗的内存联系起来（274MB），我们就会发现，考虑到它的好处，这并不是一个过于沉重的数据结构。

在计时器方面，我们看到程序的整体运行时间明显改善了2.7倍。如果我们不考虑迭代求解器，它在两种情况下都是一样的（考虑到我们使用的简单预处理程序，以及稀疏矩阵-向量乘积对立方多项式的操作浪费的事实，它并不是最佳的），优势是几乎5倍。这对于一个线性静止问题来说是相当令人印象深刻的，而对于时间依赖性和非线性问题来说，成本的节省确实要突出得多，因为在这些问题上，会多次调用装配。如果我们研究一下各个组件，我们就会更清楚地了解发生了什么，以及为什么缓存会如此高效。在MappingQGeneric的情况下，基本上每一个涉及到映射的操作都需要至少5秒的时间来运行。规范的计算运行了两个 VectorTools::integrate_difference() 函数，每个都需要5秒。(约束条件的计算比较便宜，因为它只在边界的单元格中评估映射，用于边界条件的插值)。如果我们将这5秒与填充MappingQCache的时间相比较，即5.2秒（对于所有单元，而不仅仅是活动单元），很明显，在MappingQGeneric的情况下，映射支持点的计算比其他一切都重要。也许最引人注目的结果是误差估计器的时间，标记为 "计算误差估计器"，其中MappingQGeneric的实现需要17.3秒，而MappingQCache的变体不到0.5秒。前者如此昂贵的原因（比如说，比装配贵三倍）是误差估计涉及到对面的量的评估，其中网格中的每个面都要求映射的额外点，而这些点又要经过非常昂贵的TransfiniteInterpolationManifold类。由于每个单元有六个面，这种情况比在装配中发生得更频繁。同样，MappingQCache很好地消除了重复评估，将所有涉及流形的昂贵步骤汇总到一个被重复使用的初始化调用中。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-65.cc"
*/
