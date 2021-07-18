//include/deal.II-translator/A-headers/manifold_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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
 *    @defgroup manifold Manifold description for triangulations
 * <h3>Overview</h3>
 * 本模块中的类涉及到流形的描述，在流形中，一个三角形所描述的域生活在其中。这种流形描述在一些情况下是必要的。
 * <ul>
 * <li>
 * 网格细化。每当一个单元被细化时，有必要在三角网中引入新的顶点。在最简单的情况下，我们假设构成三角网的对象是直线段、双线性曲面或三线性体。然后，下一个顶点被简单地放到旧顶点的中间（这里的
 * "中间
 * "是指先前存在的顶点位置的一个合适的平均值）。这是Triangulation类的默认行为，并由FlatManifold类描述。
 * 另一方面，如果处理弯曲的几何体，或者需要在某些方向上进行更密集细化的几何体，这并不是合适的做法。因此，从Manifold基类派生的类描述了一个域的几何形状。然后，我们可以使用
 * Triangulation::set_manifold()
 * 函数将一个从该基类派生的类对象附加到三角形对象上，并将其与manifold_id相关联（见
 * types::manifold_id), 使用 TriaAccessor::set_manifold_id()
 * 函数将此manifold_id用于三角形的单元、面或边上，这些单元、面或边应由该manifold描述，然后三角形将询问manifold对象在网格细化时应将新顶点定位在哪里。已经有几个类支持最常见的几何形状，例如，CylindricalManifold或PolarManifold，它们分别代表以圆柱坐标或极坐标描述空间时获得的几何形状。默认情况下，所有使用GridGenerator命名空间中的函数生成的弯曲几何体都会将正确的Manifold对象附加到域的弯曲部分。
 * <li>
 * 集成。当使用高阶有限元方法时，经常需要使用边界的曲线近似，而不是直线近似来计算单元项（如单元对线性系统的矩阵和右手边的贡献）。这种曲线元素的实际实现发生在Mapping类中（见
 * @ref mapping
 * 模块），然而它从这里描述的类中获得关于域的边界信息。当然，在整合边界项时也是如此（例如，不均匀的诺伊曼边界条件）。
 * <li>
 * 非零维的域。在Triangulation被嵌入高维空间的情况下，即只要Triangulation类的第二个模板参数被明确指定且大于第一个模板参数（例子见
 * step-34
 * ），流形描述对象不仅可以作为描述域的边界的几何形状的工具，而且可以作为描述域本身的工具，以防域是一个事实上是弯曲的流形。在这种情况下，人们可以使用
 * Triangulation::set_manifold()
 * 函数来指示在细化曲线时，或在使用高阶映射计算积分时，应使用何种流形描述。
 * </ul>  许多其他的例子，以及在deal.II中实现的许多理论基础，在 @ref geometry_paper  "几何学论文 "中提供。
 * 在deal.II中，流形被看作是一个点的集合，同时还有一个点与点之间距离的概念（在流形上）。新的点通常是通过在流形上提供一个局部坐标系来获得的，识别局部坐标系中的现有点（使用局部地图将其拉回，以获得其局部坐标），通过现有点的加权和找到局部坐标系中的新点，并将该点在实空间中转换回来（使用局部地图将其向前推）。实现这一机制的主要类是ChartManifold类，这也是用户可能为复杂的几何形状而重载的类。
 * 虽然这个过程在大多数感兴趣的情况下是非琐碎的，但对于大多数琐碎的几何体，如圆柱体、球体或壳体，deal.II提供了合理的实现。更复杂的例子可以用
 * step-53 和 step-54 中的技术来描述。
 * 从总体上看，这个模块的类与库中的其他各种部分相互作用。
 * @dot digraph G { graph[rankdir="TB",bgcolor="transparent"];
 *
 * node [fontname="FreeSans",fontsize=15, shape=box,height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10];
 *
 * tria       [label="Triangulation",    URL="\ref grid"]; fe
 * [label="Finite elements",    URL="\ref feall"]; mapping
 * [label="Mapping",          URL="\ref mapping"]; quadrature
 * [label="Quadrature",       URL="\ref Quadrature"]; dh
 * [label="DoFHandler",       URL="\ref dofs"]; fevalues   [label="FEValues",
 * URL="\ref feaccess"]; systems    [label="Linear systems",   URL="\ref
 * LAC"]; solvers    [label="Linear solvers",   URL="\ref Solvers"]; output
 * [label="Graphical output", URL="\ref output"]; manifold
 * [label="Manifold",         URL="\ref manifold", fillcolor="deepskyblue"];
 *
 * tria
 *
 * -> dh              [color="black",style="solid"]; fe
 *
 * -> dh                [color="black",style="solid"]; fe
 *
 * -> fevalues          [color="black",style="solid"]; mapping
 *
 * -> fevalues     [color="black",style="solid"]; quadrature
 *
 * -> fevalues  [color="black",style="solid"]; dh
 *
 * -> systems           [color="black",style="solid"]; fevalues
 *
 * -> systems     [color="black",style="solid"]; systems
 *
 * -> solvers      [color="black",style="solid"]; solvers
 *
 * -> output       [color="black",style="solid"]; manifold
 *
 * -> tria        [color="black",style="solid"]; manifold
 *
 * -> mapping     [color="black",style="solid"];
 *
 * { rank=same mapping
 *
 * -> quadrature [dir="none", color="transparent"]; quadrature
 *
 * -> fe      [dir="none", color="transparent"]; fe
 *
 * -> tria            [dir="none", color="transparent"]; }
 *
 * node [fontname="FreeSans",fontsize=12, shape=record,height=0.2,width=0.4,
 * color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
 * edge [color="gray55", weight=1];
 *
 * opencascade [label="OpenCASCADE"]; opencascade
 *
 * -> manifold [dir="none"];
 *
 *
 * node [fontname="FreeSans",fontsize=12, shape=ellipse,height=0.2,width=0.4,
 * color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
 * edge [color="gray55", weight=1];
 *
 * gmsh        [label="gmsh", URL="\ref Gmsh"]; gmsh
 *
 * -> tria       [dir="none"]; } @enddot   <h3>An example</h3> step-1
 * 已经提供了一个处理曲面几何体的简单例子，尽管那里没有详细说明。默认情况下，GridGenerator中的函数会在需要时将流形附加到网格上。在下面的每个代码片段中，我们都会调用
 * Triangulation::reset_all_manifolds()
 * 来移除这些流形，并在例子本身中处理所有的流形附件，以使流形选择的影响清晰。
 * 考虑一下那里显示的 <code>second_grid()</code>
 * 函数的这个小变化，我们只是简单地将<i>every</i>单元格细化了几次。
 * @code
 * const Point<2> center (1,0);
 * const double inner_radius = 0.5,
 *             outer_radius = 1.0;
 * Triangulation<2> triangulation;
 * GridGenerator::hyper_shell (triangulation,
 *                            center, inner_radius, outer_radius,
 *                            10);
 * // as noted above: disable all non-Cartesian manifolds
 * // for demonstration purposes:
 * triangulation.reset_all_manifolds();
 *
 * triangulation.refine_global (3);
 * @endcode
 * 这段代码导致了一个看起来像这样的网格。
*  @image html hypershell-nothing.png ""
 * 我们的意图是要得到一个类似于环形的网格。然而，由于我们没有对三角形进行描述，所发生的情况是，我们从我们告诉
 * GridGenerator::hyper_shell()
 * 要创建的圆周方向的10个粗单元开始，然后每个单元被全局细化3次。每次细化都需要一个新的顶点，它被放在现有顶点的中间，而不考虑我们可能的意图（但在代码中忽略了描述）。
 * 这很容易补救。考虑一下这段代码。
 * @code
 * const Point<2> center (1,0);
 * const SphericalManifold<2> manifold(center);
 * const double inner_radius = 0.5,
 *             outer_radius = 1.0;
 * Triangulation<2> triangulation;
 * GridGenerator::hyper_shell (triangulation,
 *                            center, inner_radius, outer_radius,
 *                            10);
 * // again disable all manifolds for demonstration purposes
 * triangulation.reset_all_manifolds();
 * triangulation.set_all_manifold_ids_on_boundary(0);
 * triangulation.set_manifold (0, manifold);
 *
 * triangulation.refine_global (3);
 * @endcode
 * 这个代码更好，产生了以下的网格。
*  @image html hypershell-boundary-only.png ""
 * 这个网格看起来更好，它忠实地再现了域的圆形内边界和外边界。然而，仍然可以在切线上发现20个结点。它们是由于每次细化单元时，内部线的新顶点只是被放置在现有线的中间（边界线的处理方式不同，因为我们附加了一个流形对象）。在有10个单元的第一次细化中，我们得到了改进的点，因为根据下面关于混合不同流形的描述，两个外部边界都提供了一个弯曲的描述。换句话说，第一次细化后的新点最终出现在可能处于直线的几何中间的地方，但不在围绕中心的圆上。
 * 这一点可以通过分配流形描述来弥补，不仅是沿边界的线，还有径向的线和单元（反过来，这些线会继承到网格细化后产生的新线）。这正是
 * GridGenerator::hyper_shell()
 * 的默认做法。为了演示，我们禁用默认的Manifold行为，然后手动复制它。
 * @code
 * const Point<2> center (1,0);
 * const SphericalManifold<2> manifold(center);
 * const double inner_radius = 0.5,
 *             outer_radius = 1.0;
 * Triangulation<2> triangulation;
 * GridGenerator::hyper_shell (triangulation,
 *                            center, inner_radius, outer_radius,
 *                            10);
 * // again disable all manifolds for demonstration purposes
 * triangulation.reset_all_manifolds();
 * // reenable the manifold:
 * triangulation.set_all_manifold_ids(0);
 * triangulation.set_manifold (0, manifold);
 * triangulation.refine_global (3);
 * @endcode
 * 这导致了以下的网状结构。
*  @image html hypershell-all.png ""
 * 那么，这有什么关系呢？毕竟，最后两个网格描述的是完全相同的领域，而且我们知道，无论选择什么样的单元，只要最大的单元的直径为零，在网格细化后，我们就能得到正确的解。
 * 这个问题有两个答案。首先，求解偏微分方程到一定精度的数值努力通常取决于单元格的<i>quality</i>，因为
 * $\|u-u_h\|_{H^1} \le Ch^p \|u\|_{H^{p+1}}$
 * 形式的误差估计中的常数 $C$
 * 取决于所有单元格中最小周长与最大内含圆的半径的最大比率等因素（对于三角形；或者对于其他类型的单元格的适当概括）。因此，创建具有尽可能好的单元的网格是值得的。可以说，对于上面所示的网格来说，这并不是一个问题，但有时却是一个问题。例如，考虑一下下面的代码和网格。
 * @code
 * const Point<2> center (1,0);
 * const SphericalManifold<2> manifold(center);
 * const double inner_radius = 0.5,
 *             outer_radius = 1.0;
 * Triangulation<2> triangulation;
 * GridGenerator::hyper_shell (triangulation,
 *                            center, inner_radius, outer_radius,
 *                            3);    // three circumferential cells
 * triangulation.reset_all_manifolds();
 * triangulation.set_all_manifold_ids_on_boundary(0);
 * triangulation.set_manifold (0, manifold);
 *
 * triangulation.refine_global (3);
 * @endcode
 *
*  @image html hypershell-boundary-only-3.png ""
 * 这里，我们在开始时只创建了三个圆周单元，并对它们进行细化，得到了所示的网格。很明显，尽管第一次细化将新的点放在中间，但我们有长宽比不好的单元。
 * 如果我们进一步推动，从半径0.8和1.0之间的更薄的圆周的粗略网格开始，并且只有三个单元（这在这里可能是不合适的，因为我们知道这是不够的，但对于在网格生成器中生成的复杂几何形状来说，可能也是不可能避免的），我们观察到以下情况。
 *
 * @code
 * const Point<2> center (1,0);
 * const SphericalManifold<2> manifold(center);
 * const double inner_radius = 0.8,
 *             outer_radius = 1.0;
 * Triangulation<2> triangulation;
 * GridGenerator::hyper_shell (triangulation,
 *                            center, inner_radius, outer_radius,
 *                            3);    // three circumferential cells
 * triangulation.reset_all_manifolds();
 * triangulation.set_all_manifold_ids_on_boundary(0);
 * triangulation.set_manifold (0, manifold);
 *
 * triangulation.refine_global (3);
 * @endcode
 *
*  @image html hypershell-boundary-thin-3.png ""
 * 这个网格在细化后既没有正确的几何形状，也没有所有单元的正面积，这对于有限元方法的工作是必要的。然而，即使从这样一个不合适的网格开始，我们也可以通过使用上述相同的代码，不仅给边界，而且给内部单元和边缘附加一个合适的几何描述，使事情顺利进行。
 * @code
 * const Point<2> center (1,0);
 * const double inner_radius = 0.8,
 *             outer_radius = 1.0;
 * Triangulation<2> triangulation;
 * GridGenerator::hyper_shell (triangulation,
 *                            center, inner_radius, outer_radius,
 *                            3);    // three circumferential cells
 *
 * triangulation.refine_global (3);
 * @endcode
 *
*  @image html hypershell-all-3.png ""
 * 在这最后一个例子中，我们终于让GridGenerator完成了它的工作，我们保留了默认的流形配置，即每个单元和面都是SphericalManifold。
 * 在这里，即使从一个初始的、选择不当的网格开始，也保留了我们的能力，可以充分地将网格细化为一个对我们有利的网格。这个例子可能是在这里制造的，但它是相关的，例如在 GridGenerator::hyper_shell() 产生的3D背景下（见这个函数的文档）。它也与 @ref GlossDistorted "关于扭曲细胞的词汇条目 "
 * 中讨论的案例有关。
 * @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
 * 。 <h3>Computing the weights for combining different manifold
 * descriptions</h3>
 * 在现实应用中，经常发生需要将不同的流形描述结合起来的情况。最简单的情况是，一个曲面描述只适用于边界，而不适用于计算域的内部。一个球的流形描述也属于这种情况，因为它需要将圆形部分的球形流形与域中心的直角描述结合起来，而球形流形在这里是无效的。
 * 一般来说，在deal.II中混合不同的流形描述的过程是通过所谓的transfinite插值实现的。它在二维中的公式，例如，在<a
 * href="https://en.wikipedia.org/wiki/Transfinite_interpolation">
 * Wikipedia</a>上有描述。给定图表上的一个点 $(u,v)$
 * ，该点在实空间中的图像由以下公式给出
 * @f{align*}{
 * \mathbf S(u,v) &= (1-v)\mathbf c_0(u)+v \mathbf c_1(u) + (1-u)\mathbf c_2(v) + u \mathbf c_3(v) \\
 * &\quad
 *
 * - \left[(1-u)(1-v) \mathbf x_0 + u(1-v) \mathbf x_1 + (1-u)v \mathbf x_2 + uv \mathbf x_3 \right]
 * @f}
 * 其中 $\bf x_0, \bf x_1, \bf x_2, \bf x_3$
 * 表示限定图像空间的四个顶点， $\bf c_0, \bf c_1, \bf c_2, \bf
 * c_3$ 是描述单元格线条的四条曲线。
 * 如果我们想根据流形找到单元格的中心（这也是在细化网格时使用的），图表是单元格
 * $(0,1)^2$ ，我们想在点 $(u,v) = (0.5, 0.5)$
 * 上评估这个公式。在这种情况下， $\mathbf c_0(0.5)$
 * 是下层面的中点位置（在deal.II的排序中以2为索引），它是由自己的流形衍生出来的，
 * $\mathbf c_1(0.5)$
 * 是上层面的中点位置（在deal.II中以3为索引）， $\mathbf
 * c_2(0.5)$ 是左边面的中点（以0为索引），而 $\mathbf c_3(0.5)$
 * 是右边面的中点。在这个公式中，面中的四个中点的权重相当于
 * $\frac{\displaystyle 1}{\displaystyle 2}$ ，四个顶点的权重相当于
 * $-\frac{\displaystyle 1}{\displaystyle 4}$
 * 。这些权重乍一看很奇怪，因为顶点的权重是负的，但这个机制是我们想要的。如果一个单元格在两个相对的面上有弯曲的描述，但在另外两个面上有直线，顶点中
 * $-\frac{\displaystyle 1}{\displaystyle 4}$
 * 的负权重与径向的两条直线的中心平衡，得到权重
 * $\frac{\displaystyle 1}{\displaystyle 2}$
 * 。因此，在曲线方向的两个中心点上取平均值，正好把新点放在中间。
 * 在三个空间维度上，面的中点的权重为 $+\frac{\displaystyle
 * 1}{\displaystyle 2}$ ，线的中点的权重为 $-\frac{\displaystyle
 * 1}{\displaystyle 4}$ ，顶点的权重为 $\frac{\displaystyle
 * 1}{\displaystyle 8}$
 * ，再次平衡了不同的实体。如果一个单元格的所有周围都是直的，那么这个公式就简化为八个顶点中每个顶点的明显权重
 * $\frac{\displaystyle 1}{\displaystyle 8}$ 。
 * 在MappingQGeneric类中，通过评估各自Gauss-Lobatto点的边界曲线
 * $(u_i,v_i)$
 * 并将其与上述公式相结合，实现了这一概念对曲线单元的多项式表示的支持点的泛化，即Gauss-Lobatto正交的节点。这些权重已经被验证为产生最佳收敛率
 * $\mathcal O(h^{k+1})$  ，也适用于非常高的多项式度数，例如
 * $k=10$  。
 * 在文献中，也使用了其他的边界描述。在9.0版本之前，deal.II使用了一种叫做拉普拉斯平滑的东西，其中应用于圆周上的节点以获得内部节点的位置的权重是通过解决单位元素上的拉普拉斯方程来确定的。然而，这导致边界层接近于曲线描述，即从单元到实心单元的映射的高导数的奇异性。
 * 如果从弯曲的边界描述到内部的直线描述的过渡做错了，通常不可能达到高阶收敛率。例如，单个单元内的拉普拉斯平滑会导致从参考单元到实数单元的映射的第四导数出现奇点，将边界处单元的收敛率限制在3（如果在二维测量全局L2误差，则为3.5）。其他更粗糙的策略，比如完全忽略两个不同流形的存在，只是简单地计算直坐标系中高阶映射的附加点，会导致更差的收敛率。另一方面，目前在deal.II中的实现，在这方面已经得到了广泛的验证，应该表现得很好。
 * 将曲线边界表示与平面内部表示相混合的不良策略显然也反映了网格质量。例如，在上述只有3个圆周单元的情况下，导致以下的网格是用拉普拉斯流形平滑，而不是像deal.II中实现的从边界插值。
*  @image html hypershell-boundary-only-3-old.png ""
 * 为了使用一个更实际的例子，考虑对一个球的细化，在球面上附加一个SphericalManifold。拉普拉斯类型的平滑给出了以下相当差的网格。
*  @image html hyperball-mesh-smoothing-laplace.png ""
 * 如果我们改用从转折插值得到的权重，情况就会有很大的改善。
*  @image html hyperball-mesh-smoothing-interpolate.png ""
 * 当然，我们可以通过将TransfiniteInterpolationManifold应用于整个域（除了连接SphericalManifold的边界）来得到更好的网格，如该类中的数字所示，但原则上，在deal.II中实现的网格平滑已经和单从边界描述得到的一样好。
 *
 * @ingroup grid   @author  Luca Heltai, 2013, Martin Kronbichler, 2017
 *
 */


