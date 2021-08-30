/**
@page step_47 The step-47 tutorial program
This tutorial depends on step-12.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Whatstheissue"> What's the issue? </a>
        <li><a href="#Whattodoinstead"> What to do instead? </a>
        <li><a href="#DerivationoftheC0IPmethod"> Derivation of the C0IP method </a>
      <ul>
        <li><a href="#ConvergenceRates">Convergence Rates </a>
      </ul>
        <li><a href="#OtherBoundaryConditions">Other Boundary Conditions</a>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Themainclass">The main class</a>
      <ul>
        <li><a href="#Assemblingthelinearsystem">Assembling the linear system</a>
        <li><a href="#Solvingthelinearsystemandpostprocessing">Solving the linear system and postprocessing</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#TestresultsoniQsub2subiiQsub2subiwithigammapp1iigammapp1i">Test results on <i>Q<sub>2</sub></i><i>Q<sub>2</sub></i> with <i>&gamma; = p(p+1)</i><i>&gamma; = p(p+1)</i> </a>
        <li><a href="#TestresultsoniQsub3subiiQsub3subiwithigammapp1iigammapp1i">Test results on <i>Q<sub>3</sub></i><i>Q<sub>3</sub></i> with <i>&gamma; = p(p+1)</i><i>&gamma; = p(p+1)</i> </a>
        <li><a href="#TestresultsoniQsub4subiiQsub4subiwithigammapp1iigammapp1i">Test results on <i>Q<sub>4</sub></i><i>Q<sub>4</sub></i> with <i>&gamma; = p(p+1)</i><i>&gamma; = p(p+1)</i> </a>
        <li><a href="#TestresultsoniQsub2subiiQsub2subiwithigamma1iigamma1i">Test results on <i>Q<sub>2</sub></i><i>Q<sub>2</sub></i> with <i>&gamma; = 1</i><i>&gamma; = 1</i> </a>
        <li><a href="#TestresultsoniQsub2subiiQsub2subiwithigamma2iigamma2i">Test results on <i>Q<sub>2</sub></i><i>Q<sub>2</sub></i> with <i>&gamma; = 2</i><i>&gamma; = 2</i> </a>
        <li><a href="#Conclusionsforthechoiceofthepenaltyparameter"> Conclusions for the choice of the penalty parameter </a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#Derivationforthesimplysupportedplates"> Derivation for the simply supported plates </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-47/doc/intro.dox

 <br> 

<i>
This program was contributed by Natasha Sharma, Guido Kanschat, Timo
Heister, Wolfgang Bangerth, and Zhuoran Wang.


The first author would like to acknowledge the support of NSF Grant
No. DMS-1520862.
Timo Heister and Wolfgang Bangerth acknowledge support through NSF
awards DMS-1821210, EAR-1550901, and OAC-1835673.
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个程序处理的是<a
href="https://en.wikipedia.org/wiki/Biharmonic_equation">biharmonic
equation</a>。

@f{align*}{
  \Delta^2 u(\mathbf x) &= f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega.


@f}

这个方程出现在薄结构的建模中，如体育场的屋顶。当然，这些物体在现实中是三维的，其横向范围与垂直厚度的长宽比很大，但人们通常可以通过对内力在垂直方向上的变化作出假设，将这些结构非常准确地建模为二维的。这些假设导致了上面的方程式。

该模型通常有两种不同的类型，取决于施加什么样的边界条件。第一种情况。

@f{align*}{
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \Delta u(\mathbf x) &= h(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega,


@f}

对应于连接到高度为 $g(\mathbf x)$ 的墙顶的薄结构的边缘，这样作用在结构上的弯曲力为 $h(\mathbf x)$ ；在大多数物理情况下，我们会有 $h=0$ ，对应于结构只是坐在墙顶。

在边界值的第二种可能情况下，我们将有

@f{align*}{
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \frac{\partial u(\mathbf x)}{\partial \mathbf n} &= j(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega.


@f}

这相当于一个 "钳制 "的结构，对于这个结构来说，非零的 $j(\mathbf x)$ 意味着对水平面有一定的角度。

就像拉普拉斯方程的迪里希特和诺依曼边界条件一样，当然有可能在边界的一部分有一种边界条件，而在其余部分有另一种。




<a name="Whatstheissue"></a><h3> What's the issue? </h3>


该方程的基本问题是它的解有四个导数。在我们在步骤3、步骤4和其他几个教程中处理的拉普拉斯方程的情况下，人们乘以一个测试函数，进行积分，通过部分积分，最后在测试函数和试验函数上只得到一个导数--对于全局连续的函数，人们可以做到这一点，但在单元之间的界面上可能有结点。导数可能不在界面上定义，但那是在一个较低维的流形上（所以不在积分值中显示出来）。

但是对于双调子方程，如果按照同样的程序，在整个领域（即所有单元的联盟）上使用积分，最终会在测试函数和试验函数上各产生两个导数。如果使用通常的片状多项式函数，并在单元格界面上有结点，那么第一个导数将产生一个不连续的梯度，第二个导数在界面上有delta函数--但由于测试函数和试验函数的第二个导数都产生一个delta函数，我们将尝试对两个delta函数的乘积进行积分。例如，在1d中， $\varphi_i$ 是通常的片状线性 "帽子函数"，我们会得到这样的积分

@f{align*}{
  \int_0^L (\Delta \varphi_i) (\Delta \varphi_j)
  =
  \int_0^L
  \frac 1h \left[\delta(x-x_{i-1}) - 2\delta(x-x_i) + \delta(x-x_{i+1})\right]
  \frac 1h \left[\delta(x-x_{j-1}) - 2\delta(x-x_j) + \delta(x-x_{j+1})\right]


@f}

其中 $x_i$ 是定义形状函数 $\varphi_i$ 的节点位置， $h$ 是网格大小（假设为均匀）。问题是，积分中的delta函数是用以下关系定义的

@f{align*}{
  \int_0^L \delta(x-\hat x) f(x) \; dx
  =
  f(\hat x).


@f}

但这只有在以下情况下才行得通：(i)  $f(\cdot)$ 实际上在 $\hat x$ 处定义良好，(ii) 它是有限的。另一方面，以下形式的积分

@f{align*}{
\int_0^L \delta(x-x_i) \delta (x-x_i)


@f}

是没有意义的。类似的推理也可以适用于2D和3D的情况。

换句话说。这种试图在整个领域内进行整合，然后再按部分进行整合的方法不可能成功。

历史上，数值分析家试图通过发明 "C<sup>1</sup>连续 "的有限元来解决这个问题，也就是说，使用的形状函数不仅是连续的，而且还有连续的一导数。这是诸如阿吉里斯元素、克拉夫-托歇尔元素和其他元素的领域，这些元素都是在20世纪60年代末开发的。从二十一世纪的角度来看，它们的结构只能说是怪异的。如果想使用一般的网格，它们的实现也是非常麻烦的。因此，它们在很大程度上已经失去了作用，deal.II目前不包含这些形状函数的实现。




<a name="Whattodoinstead"></a><h3> What to do instead? </h3>


那么，如何解决此类问题呢？这在一定程度上取决于边界条件。如果有第一组边界条件，即，如果方程是

@f{align*}{
  \Delta^2 u(\mathbf x) &= f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \Delta u(\mathbf x) &= h(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega,


@f}

那么下面的诀窍就能起作用（至少如果域是凸的，见下文）。正如我们通过引入第二个变量从常规拉普拉斯方程中得到步骤20的混合拉普拉斯方程一样，我们可以在这里引入一个变量 $v=\Delta u$ ，然后可以用下面的 "混合 "系统代替上面的方程。

@f{align*}{


  -\Delta u(\mathbf x) +v(\mathbf x) &= 0
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\


  -\Delta v(\mathbf x) &= -f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  v(\mathbf x) &= h(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega.


@f}

换句话说，我们最终得到的实质上是 $u,v$ 的两个耦合的拉普拉斯方程组，每个方程组都有迪里希勒型边界条件。我们知道如何解决这样的问题，使用第20步或第22步的技术为这个系统构造良好的求解器和预处理器应该不是很困难。所以这种情况很容易处理。

 @note 值得指出的是，这只适用于边界有角的域，如果该域也是凸的--换句话说，如果没有重入角。   这听起来是一个相当随意的条件，但考虑到以下两个事实，它是有意义的。原始双调方程的解必须满足  $u\in H^2(\Omega)$  。另一方面，上面的混合系统重述表明， $u$ 和 $v$ 都满足 $u,v\in H^1(\Omega)$ ，因为这两个变量只解决一个泊松方程。换句话说，如果我们想确保混合问题的解 $u$ 也是原来的偏谐方程的解，那么我们需要能够以某种方式保证 $-\Delta u=v$ 的解实际上比只是 $H^1(\Omega)$ 更光滑。这一点可以作如下论证。对于凸域，<a href="https://en.wikipedia.org/wiki/Elliptic_operator#Elliptic_regularity_theorem">"elliptic
  regularity"</a>意味着，如果右侧 $v\in H^s$ ，那么 $u\in H^{s+2}$ 如果域是凸的并且边界足够光滑。(如果域的边界足够光滑，也可以保证这一点--但边界没有角的域在现实生活中不是很实用。)   我们知道 $v\in H^1$ ，因为它解决了方程 $-\Delta v=f$ ，但我们仍然留下了边界凸性的条件；我们可以证明多边形的凸域足以保证在这种情况下 $u\in H^2$ （光滑的有界凸域将导致 $u\in H^3$ ，但我们不需要这么多规则）。另一方面，如果域不是凸的，我们就不能保证混合系统的解在 $H^2$ 中，因此可能会得到一个不能等于原始偏谐方程的解。

更复杂的情况是，如果我们有 "钳制 "的边界条件，也就是说，如果方程看起来像这样。

@f{align*}{
  \Delta^2 u(\mathbf x) &= f(\mathbf x)
  \qquad \qquad &&\forall \mathbf x \in \Omega, \\
  u(\mathbf x) &= g(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega, \\
  \frac{\partial u(\mathbf x)}{\partial \mathbf n} &= j(\mathbf x) \qquad \qquad
  &&\forall \mathbf x \in \partial\Omega.


@f}

混合系统的相同技巧在这里不起作用，因为我们最终会对 $u$ 采用<i>both</i>的迪里切特和诺伊曼边界条件，但对 $v$ 则没有。


在20世纪90年代和21世纪初，这一难题的解决方案随着非连续Galerkin方法的出现而到来。与使用<i>discontinuous</i>形状函数处理拉普拉斯方程，通过惩罚不连续的大小来获得每个形状函数上有一个导数的方程的方案一样，我们可以使用一个使用<i>continuous</i>（但不是 $C^1$ 连续）形状函数的方案，惩罚导数的跳跃来获得每个形状函数上有两个导数的方案。与拉普拉斯方程的内部惩罚（IP）方法相类似，这种用于双调方程的方案通常被称为 $C^0$ IP（或C0IP）方法，因为它使用 $C^0$ （连续但不连续可微）形状函数与内部惩罚公式。




<a name="DerivationoftheC0IPmethod"></a><h3> Derivation of the C0IP method </h3>


我们以Susanne Brenner和Li-Yeng Sung在 "C  $^0$  多边形域上线性四阶边界值问题的内部惩罚方法" @cite Brenner2005 中提出的 $C^0$ IP方法为基础，该方法是针对具有 "钳制 "边界条件的双谐波方程而提出的。

如前所述，这种方法依赖于使用 $C^0$ 拉格朗日有限元，其中 $C^1$ 的连续性要求被放宽，并被内部惩罚技术所取代。为了推导这个方法，我们考虑一个 $C^0$ 形状函数 $v_h$ ，它在 $\partial\Omega$ 上消失。我们引入符号 $ \mathbb{F} $ 作为 $\mathbb{T}$ 的所有面的集合， $ \mathbb{F}^b $ 作为边界面的集合， $ \mathbb{F}^i $ 作为内部面的集合，供下面进一步使用。由于 $v_h$ 的高阶导数在每个界面 $e\in \mathbb{F}$ 上有两个值（由两个单元 $K_{+},K_{-} \in \mathbb{T}$ 共享），我们通过在 $e$ 上定义以下单值函数来应对这一不连续性。

@f{align*}{
  \jump{\frac{\partial^k v_h}{\partial \mathbf n^k}}
  &=
  \frac{\partial^k v_h|_{K_+}}{\partial \mathbf n^k} \bigg |_e


  - \frac{\partial^k v_h|_{K_-}}{\partial \mathbf n^k} \bigg |_e,
  \\
  \average{\frac{\partial^k v_h}{\partial \mathbf n^k}}
  &=
  \frac{1}{2}
  \bigg( \frac{\partial^k v_h|_{K_+}}{\partial \mathbf n^k} \bigg |_e
  + \frac{\partial^k v_h|_{K_-}}{\partial \mathbf n^k} \bigg |_e \bigg )


@f}

为 $k =1,2$ （即为梯度和二阶导数矩阵），其中 $\mathbf n$ 表示从 $K_+$ 指向 $K_-$ 的一个单位向量法线。在文献中，这些函数分别被称为 "跳跃 "和 "平均 "操作。

为了得到 $C^0$ IP的近似值 $u_h$ ，我们将双调方程乘以 $v_h$ ，然后对 $\Omega$ 进行积分。如上所述，我们不能用这些形状函数对 $\Omega$ 的所有部分进行积分，但我们可以对每个单元进行积分，因为这些形状函数只是每个单元的多项式。因此，我们首先在每个网格单元 $K \in {\mathbb{T}}$ 上使用下面的分项积分公式。

@f{align*}{
  \int_K v_h (\Delta^2 w_h)
  &= \int_K v_h (\nabla\cdot\nabla) (\Delta w_h)
  \\
  &= -\int_K \nabla v_h \cdot (\nabla \Delta w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n).


@f}

在这一点上，我们有两个选择。我们可以再一次整合域项的 $\nabla\Delta w_h$ ，得到

@f{align*}{
  \int_K v_h (\Delta^2 w_h)
  &= \int_K (\Delta v_h) (\Delta w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n)


     -\int_{\partial K} (\nabla v_h \cdot \mathbf n) \Delta w_h.


@f}

由于各种原因，这被证明是一个对我们的目的没有用处的变体。

相反，我们要做的是认识到 $\nabla\Delta w_h = \text{grad}\,(\text{div}\,\text{grad}\, w_h)$  ，我们可以将这些操作重新排序为 $\nabla\Delta w_h = \text{div}\,(\text{grad}\,\text{grad}\, w_h)$ ，其中我们通常写成 $\text{grad}\,\text{grad}\, w_h = D^2 w_h$ ，表示这是第二导数的 "黑森 "矩阵。通过这样的重新排序，我们现在可以整合发散，而不是梯度算子，我们得到以下结果。

@f{align*}{
  \int_K v_h (\Delta^2 w_h)
  &= \int_K (\nabla \nabla v_h) : (\nabla \nabla w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n)


     -\int_{\partial K} (\nabla v_h \otimes \mathbf n) : (\nabla\nabla w_h)
  \\
  &= \int_K (D^2 v_h) : (D^2 w_h)
     +\int_{\partial K} v_h (\nabla \Delta w_h \cdot \mathbf n)


     -\int_{\partial K} (\nabla v_h) \cdot (D^2 w_h \mathbf n).


@f}

这里，冒号表示对其左边和右边的矩阵的指数进行双缩，即两个张量之间的标量乘积。两个向量 $a \otimes b$ 的外积可以得到矩阵 $(a \otimes b)_{ij} = a_i b_j$  。

然后，我们对所有单元格 $K \in  \mathbb{T}$ 进行求和，并考虑到这意味着每个内部面在求和中出现两次。因此，如果我们把所有的东西分成细胞内部的积分之和和细胞界面的单独之和，我们就可以使用上面定义的跳跃和平均运算符。还有两个步骤。首先，由于我们的形状函数是连续的，形状函数的梯度可能是不连续的，但连续性保证了实际上只有梯度的法向分量是不连续的，而切向分量是连续的。其次，当网格大小为零时，产生的离散公式并不稳定，为了得到一个稳定的公式，收敛到正确的解，我们需要增加以下条款。

@f{align*}{


-\sum_{e \in \mathbb{F}} \int_{e}
  \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
  \jump{\frac{\partial u_h}{\partial \mathbf n}}
+ \sum_{e \in \mathbb{F}}
  \frac{\gamma}{h_e}\int_e
  \jump{\frac{\partial v_h}{\partial \mathbf n}}
  \jump{\frac{\partial u_h}{\partial \mathbf n}}.


@f}

然后，在进行出现的取消后，我们得出以下双调子方程的C0IP表述：找到 $u_h$ ，使 $u_h =
g$ 对 $\partial \Omega$ 和

@f{align*}{
\mathcal{A}(v_h,u_h)&=\mathcal{F}(v_h) \quad \text{holds for all test functions } v_h,


@f}

其中

@f{align*}{
\mathcal{A}(v_h,u_h):=&\sum_{K \in \mathbb{T}}\int_K D^2v_h:D^2u_h \ dx
\\
&


 -\sum_{e \in \mathbb{F}} \int_{e}
  \jump{\frac{\partial v_h}{\partial \mathbf n}}
  \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds


 -\sum_{e \in \mathbb{F}} \int_{e}
 \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
 \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
\\
&+ \sum_{e \in \mathbb{F}}
 \frac{\gamma}{h_e}
 \int_e
 \jump{\frac{\partial v_h}{\partial \mathbf n}}
 \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds,


@f}

和

@f{align*}{
\mathcal{F}(v_h)&:=\sum_{K \in \mathbb{T}}\int_{K} v_h f \ dx


-
\sum_{e \in \mathbb{F}, e\subset\partial\Omega}
\int_e \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}} j \ ds
+
\sum_{e \in \mathbb{F}, e\subset\partial\Omega}
\frac{\gamma}{h_e}
\int_e \jump{\frac{\partial v_h}{\partial \mathbf n}} j \ ds.


@f}

这里， $\gamma$ 是惩罚参数，它既弱化了边界条件的执行

@f{align*}{
\frac{\partial u(\mathbf x)}{\partial \mathbf n} = j(\mathbf x)


@f}

在边界界面 $e \in \mathbb{F}^b$ 上，也确保在极限 $h\rightarrow 0$ 中， $u_h$ 收敛为 $C^1$ 连续函数。   $\gamma$ 被选择为足够大以保证方法的稳定性。我们将在下面的程序中讨论我们的选择。




<a name="ConvergenceRates"></a><h4>Convergence Rates </h4> 在多边形域上，双调方程的弱解 $u$ 存在于 $H^{2 +\alpha}(\Omega)$ 中，其中 $\alpha \in(1/2, 2]$ 是由 $\Omega$ 的角的内角决定。例如，只要 $\Omega$ 是凸的， $\alpha=1$ ； $\alpha$ 可能小于1，如果域有重心角，但如果所有内角之一接近 $\pi$ ， $\alpha$ 就接近于 $1$  。


现在假设 $C^0$ IP解 $u_h$ 被 $C^0$ 形状函数近似，其程度为多项式 $p \ge 2$  。然后，上面概述的离散化产生了下面讨论的收敛率。


<b>Convergence in the $C^0$ IP-norm</b>

理想情况下，我们希望测量 "能量准则" $\|D^2(u-u_h)\|$ 中的收敛性。然而，这并不可行，因为同样地，离散解 $u_h$ 并没有两个（弱）导数。相反，我们可以定义一个离散的（ $C^0$  IP）半规范，"等同于 "能量规范，如下所示。

@f{align*}{
 |u_h|_{h}^2 :=
 \sum\limits_{K \in \mathbb{T}} \big|u_h\big|_{H^2(K)}^2
 +
 \sum\limits_{e \in \mathbb{F} }
 \frac{\gamma }{h_e} \left\|
 \jump{\frac{\partial u_h}{\partial \mathbf n}} \right\|_{L^2(e)}^2.


@f}



在这个半规范中，上面提到的论文中的理论得出，我们可以期望

@f{align*}{
 |u-u_h|_{h}^2 = {\cal O}(h^{p-1}),


@f}

我们知道，拉普拉斯方程的通常离散化的收敛率是真实的，这与我们所期望的差不多。

当然，只有在精确解足够平滑的情况下，这才是真的。事实上，如果 $f \in H^m(\Omega)$ 有 $m \ge 0$ ， $u \in H^{2+\alpha}(\Omega)$ 其中 $ 2 < 2+\alpha  \le m+4$ ，那么 $C^0$ IP方法的收敛率是 $\mathcal{O}(h^{\min\{p-1, \alpha\}})$ 。换句话说，只有当解平滑到 $\alpha\ge p-1$ 时，才能期待最佳收敛率；这只有在以下情况下才会发生：(i) 域是凸的，边界足够平滑，(ii)  $m\ge p-3$  。当然，在实践中，解决方案是什么就是什么（与我们选择的多项式程度无关），那么最后一个条件可以等同于说，如果 $p$ 不大，那么选择 $m$ 大肯定没有意义。换句话说， $p$ 唯一合理的选择是 $p\le
m+3$ ，因为更大的多项式度数不会导致更高的收敛阶数。

就本程序而言，我们有点懒得实际实现这个等价的语义规则--尽管这并不难，而且会成为一个很好的练习。相反，我们将简单地在程序中检查 "坏的" $H^2$ 语义规则是什么？

@f{align*}{
 \left(|u_h|^\circ_{H^2}\right)^2
 :=
 \sum\limits_{K \in \mathbb{T}} \big|u_h\big|_{H^2(K)}^2
 =
 \sum\limits_{K \in \mathbb{T}} \big|D^2 u_h\big|_{L_2}^2


@f}

产量。从理论的角度来看，这个准则的收敛率当然不可能比<i>worse</i>的收敛率高，因为它只包含了必要条件的一个子集，但至少可以想象到它会更好。还有一种情况是，即使程序中存在一个错误，我们也能得到最佳收敛率，而这个错误只会在 $|\cdot|_h$ 中出现的额外条款中显示出次优的收敛率。但是，人们可能希望，如果我们在破碎规范和下面讨论的规范中得到最优速率，那么这个程序确实是正确的。结果部分将证明，我们在所有显示的规范中都得到了最优率。


<b>Convergence in the $L_2$-norm</b>

在  $L_2$  -norm中的最佳收敛率是  $\mathcal{O}(h^{p+1})$  提供的  $p \ge 3$  。更多细节可以在  @cite Engel2002  的定理4.6中找到。

下面的程序中默认的是选择 $p=2$ 。在这种情况下，该定理并不适用，事实上，人们只能得到 $\mathcal{O}(h^2)$ 而不是 $\mathcal{O}(h^3)$ ，我们将在结果部分展示。


<b>Convergence in the $H^1$-seminorm</b>

鉴于我们在最好的情况下对相当于 $H^2$ 半规范的规范期待 $\mathcal{O}(h^{p-1})$ ，对 $L_2$ 规范期待 $\mathcal{O}(h^{p+1})$ ，人们可能会问在 $H^1$ 半规范中会发生什么，它介于其他两种规范之间。一个合理的猜测是，我们应该期待 $\mathcal{O}(h^{p})$  。可能在某个地方有一篇论文证明了这一点，但我们在下面也验证了这个猜想在实验上是真实的。




<a name="OtherBoundaryConditions"></a><h3>Other Boundary Conditions</h3>


我们注意到，对于具有其他边界条件的偏谐方程--例如，对于第一组边界条件即 $u(\mathbf x) =
g(\mathbf x)$ 和 $\Delta u(\mathbf x)= h(\mathbf x)$ 上的 $\partial\Omega$ --的IP方法的推导，可以通过对 $\mathcal{A}(\cdot,\cdot)$ 和 $\mathcal{F}(\cdot)$ 的适当修改得到，这些修改在书中的章节 @cite Brenner2011 中描述。




<a name="Thetestcase"></a><h3>The testcase</h3>


剩下要描述的最后一步是这个程序的求解内容。像往常一样，三角函数是一个既好又坏的选择，因为它不在我们可能寻求解决方案的任何多项式空间中，同时又比实际解决方案通常更平滑（这里，它在 $C^\infty$ 中，而实际解决方案在凸多边形域上通常只在 $H^3$ 左右，如果域不是凸的，则在 $H^2$ 和 $H^3$ 之间）。但是，由于我们没有办法用相对简单的公式来描述现实问题的解决方案，所以我们就用下面的方法，在单位平方的域上 $\Omega$  。

@f{align*}{
  u = \sin(\pi x) \sin(\pi y).


@f}

因此，我们需要选择以下条件作为边界条件。

@f{align*}{
  g &= u|_{\partial\Omega} = \sin(\pi x) \sin(\pi y)|_{\partial\Omega},
  \\
  j &= \frac{\partial u}{\partial\mathbf n}|_{\partial\Omega}
  \\
    &= \left.\begin{pmatrix}
                \pi\cos(\pi x) \sin(\pi y) \\
                \pi\sin(\pi x) \cos(\pi y)
             \end{pmatrix}\right|_{\partial\Omega} \cdot \mathbf n.


@f}

右手边很容易计算为

@f{align*}{
  f = \Delta^2 u = 4 \pi^4 \sin(\pi x) \sin(\pi y).


@f}

该程序有 `ExactSolution::Solution` 和 `ExactSolution::RightHandSide` 类来编码这一信息。


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
 * 前面的几个include文件已经在前面的例子中使用过了，所以我们在这里不再解释它们的含义。该程序的主要结构与例如 step-4 的结构非常相似，因此我们包含了许多相同的头文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/sparse_direct.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/mapping_q.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * @endcode
 * 
 * 最有趣的两个头文件将是这两个。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_interface_values.h> 
 * #include <deal.II/meshworker/mesh_loop.h> 
 * 
 * @endcode
 * 
 * 其中第一个文件负责提供FEInterfaceValues类，该类可用于评估单元间界面的形状函数（或其梯度）的跳跃或平均值等数量。这个类在评估C0IP公式中出现的惩罚项时将相当有用。
 * 

 * 
 * 
 * @code
 * #include <fstream> 
 * #include <iostream> 
 * #include <cmath> 
 * 
 * namespace Step47 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 在下面的命名空间中，让我们定义精确解，我们将与数值计算的解进行比较。它的形式是 $u(x,y) = \sin(\pi x) \sin(\pi y)$ （只实现了2d的情况），该命名空间还包含一个对应于产生该解的右手边的类。
 * 

 * 
 * 
 * @code
 *   namespace ExactSolution 
 *   { 
 *     using numbers::PI; 
 * 
 *     template <int dim> 
 *     class Solution : public Function<dim> 
 *     { 
 *     public: 
 *       static_assert(dim == 2, "Only dim==2 is implemented."); 
 * 
 *       virtual double value(const Point<dim> &p, 
 *                            const unsigned int /*component*/ = 0) const override 
 *       { 
 *         return std::sin(PI * p[0]) * std::sin(PI * p[1]); 
 *       } 
 * 
 *       virtual Tensor<1, dim> 
 *       gradient(const Point<dim> &p, 
 *                const unsigned int /*component*/ = 0) const override 
 *       { 
 *         Tensor<1, dim> r; 
 *         r[0] = PI * std::cos(PI * p[0]) * std::sin(PI * p[1]); 
 *         r[1] = PI * std::cos(PI * p[1]) * std::sin(PI * p[0]); 
 *         return r; 
 *       } 
 * 
 *       virtual void 
 *       hessian_list(const std::vector<Point<dim>> &       points, 
 *                    std::vector<SymmetricTensor<2, dim>> &hessians, 
 *                    const unsigned int /*component*/ = 0) const override 
 *       { 
 *         for (unsigned i = 0; i < points.size(); ++i) 
 *           { 
 *             const double x = points[i][0]; 
 *             const double y = points[i][1]; 
 * 
 *             hessians[i][0][0] = -PI * PI * std::sin(PI * x) * std::sin(PI * y); 
 *             hessians[i][0][1] = PI * PI * std::cos(PI * x) * std::cos(PI * y); 
 *             hessians[i][1][1] = -PI * PI * std::sin(PI * x) * std::sin(PI * y); 
 *           } 
 *       } 
 *     }; 
 * 
 *     template <int dim> 
 *     class RightHandSide : public Function<dim> 
 *     { 
 *     public: 
 *       static_assert(dim == 2, "Only dim==2 is implemented"); 
 * 
 *       virtual double value(const Point<dim> &p, 
 *                            const unsigned int /*component*/ = 0) const override 
 * 
 *       { 
 *         return 4 * std::pow(PI, 4.0) * std::sin(PI * p[0]) * 
 *                std::sin(PI * p[1]); 
 *       } 
 *     }; 
 *   } // namespace ExactSolution 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainclass"></a> 
 * <h3>The main class</h3>
 * 

 * 
 * 以下是本教程程序的主类。它具有许多其他教程程序的结构，其内容和后面的构造函数应该没有什么特别令人惊讶的地方。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BiharmonicProblem 
 *   { 
 *   public: 
 *     BiharmonicProblem(const unsigned int fe_degree); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid(); 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void compute_errors(); 
 *     void output_results(const unsigned int iteration) const; 
 * 
 *     Triangulation<dim> triangulation; 
 * 
 *     MappingQ<dim> mapping; 
 * 
 *     FE_Q<dim>                 fe; 
 *     DoFHandler<dim>           dof_handler; 
 *     AffineConstraints<double> constraints; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 *   }; 
 * 
 *   template <int dim> 
 *   BiharmonicProblem<dim>::BiharmonicProblem(const unsigned int fe_degree) 
 *     : mapping(1) 
 *     , fe(fe_degree) 
 *     , dof_handler(triangulation) 
 *   {} 
 * 
 * @endcode
 * 
 * 接下来是创建初始网格（一次精炼的单元格）和设置每个网格的约束、向量和矩阵的函数。同样，这两个函数与之前的许多教程程序基本没有变化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BiharmonicProblem<dim>::make_grid() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, 0., 1.); 
 *     triangulation.refine_global(1); 
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "Total number of cells: " << triangulation.n_cells() 
 *               << std::endl; 
 *   } 
 * 
 *   template <int dim> 
 *   void BiharmonicProblem<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 * 
 *     constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              ExactSolution::Solution<dim>(), 
 *                                              constraints); 
 *     constraints.close(); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, true); 
 *     sparsity_pattern.copy_from(dsp); 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Assemblingthelinearsystem"></a> 
 * <h4>Assembling the linear system</h4>
 * 

 * 
 * 下面的几段代码更有意思。它们都与线性系统的组装有关。虽然组装单元格内部项的难度不大--这在本质上就像组装拉普拉斯方程的相应项一样，你已经在 step-4 或 step-6 中看到了这是如何工作的，例如，困难在于公式中的惩罚项。这需要在单元格的界面上对形状函数的梯度进行评估。因此，至少需要使用两个FEFaceValues对象，但如果其中一个面是自适应细化的，那么实际上需要一个FEFaceValues和一个FESubfaceValues对象；我们还需要跟踪哪些形状函数在哪里，最后我们需要确保每个面只被访问一次。所有这些对于我们真正想要实现的逻辑（即双线性形式中的惩罚项）来说都是一笔不小的开销。因此，我们将使用FEInterfaceValues类--这是deal.II中的一个辅助类，它允许我们抽象出两个FEFaceValues或FESubfaceValues对象，直接访问我们真正关心的东西：跳跃、平均等。
 * 

 * 
 * 但这还没有解决我们的问题，即当我们在所有单元格和它们的所有面中循环时，必须跟踪我们已经访问过哪些面。为了使这个过程更简单，我们使用了 MeshWorker::mesh_loop() 函数，它为这个任务提供了一个简单的接口：基于WorkStream命名空间文档中概述的想法， MeshWorker::mesh_loop() 需要三个函数对单元、内部面和边界面进行工作。这些函数在抓取对象上工作，以获得中间结果，然后将其计算结果复制到复制数据对象中，由一个复制器函数将其复制到全局矩阵和右侧对象中。
 * 

 * 
 * 然后，下面的结构提供了这种方法所需的从头开始和复制对象。你可以查阅WorkStream命名空间以及 @ref threads "多处理器并行计算 "模块，了解更多关于它们通常如何工作的信息。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct ScratchData 
 *   { 
 *     ScratchData(const Mapping<dim> &      mapping, 
 *                 const FiniteElement<dim> &fe, 
 *                 const unsigned int        quadrature_degree, 
 *                 const UpdateFlags         update_flags, 
 *                 const UpdateFlags         interface_update_flags) 
 *       : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags) 
 *       , fe_interface_values(mapping, 
 *                             fe, 
 *                             QGauss<dim - 1>(quadrature_degree), 
 *                             interface_update_flags) 
 *     {} 
 * 
 *     ScratchData(const ScratchData<dim> &scratch_data) 
 *       : fe_values(scratch_data.fe_values.get_mapping(), 
 *                   scratch_data.fe_values.get_fe(), 
 *                   scratch_data.fe_values.get_quadrature(), 
 *                   scratch_data.fe_values.get_update_flags()) 
 *       , fe_interface_values(scratch_data.fe_values.get_mapping(), 
 *                             scratch_data.fe_values.get_fe(), 
 *                             scratch_data.fe_interface_values.get_quadrature(), 
 *                             scratch_data.fe_interface_values.get_update_flags()) 
 *     {} 
 * 
 *     FEValues<dim>          fe_values; 
 *     FEInterfaceValues<dim> fe_interface_values; 
 *   }; 
 * 
 *   struct CopyData 
 *   { 
 *     CopyData(const unsigned int dofs_per_cell) 
 *       : cell_matrix(dofs_per_cell, dofs_per_cell) 
 *       , cell_rhs(dofs_per_cell) 
 *       , local_dof_indices(dofs_per_cell) 
 *     {} 
 * 
 *     CopyData(const CopyData &) = default; 
 * 
 *     CopyData(CopyData &&) = default; 
 * 
 *     ~CopyData() = default; 
 * 
 *     CopyData &operator=(const CopyData &) = default; 
 * 
 *     CopyData &operator=(CopyData &&) = default; 
 * 
 *     struct FaceData 
 *     { 
 *       FullMatrix<double>                   cell_matrix; 
 *       std::vector<types::global_dof_index> joint_dof_indices; 
 *     }; 
 * 
 *     FullMatrix<double>                   cell_matrix; 
 *     Vector<double>                       cell_rhs; 
 *     std::vector<types::global_dof_index> local_dof_indices; 
 *     std::vector<FaceData>                face_data; 
 *   }; 
 * 
 * @endcode
 * 
 * 更有趣的部分是我们实际组装线性系统的地方。从根本上说，这个函数有五个部分。
 * 

 * 
 * - `cell_worker`λ函数的定义，这是一个定义在`assemble_system()`函数中的小函数，它将负责计算单个单元上的局部积分。它将在`ScratchData`类的副本上工作，并将其结果放入相应的`CopyData`对象。
 * 

 * 
 * - `face_worker` lambda函数的定义，它将对单元格之间的界面上的所有项进行积分。
 * 

 * 
 * - 定义了`boundary_worker`函数，对位于域的边界上的单元面做同样的工作。
 * 

 * 
 * - `copier`函数的定义，该函数负责将前面三个函数中的所有数据复制到单个单元的复制对象中，并复制到全局矩阵和右侧。
 * 

 * 
 * 第五部分是我们把所有这些都集中在一起。
 * 

 * 
 * 让我们轮流浏览一下这些组装所需的每一块。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BiharmonicProblem<dim>::assemble_system() 
 *   { 
 *     using Iterator = typename DoFHandler<dim>::active_cell_iterator; 
 * 
 * @endcode
 * 
 * 第一部分是`cell_worker'，它在细胞内部进行组装。它是一个（lambda）函数，以一个单元格（输入）、一个抓取对象和一个复制对象（输出）为参数。它看起来像许多其他教程程序的装配函数，或者至少是所有单元格的循环主体。
 * 

 * 
 * 我们在这里整合的条款是单元格对全局矩阵的贡献
 * @f{align*}{
 * A^K_{ij} = \int_K \nabla^2\varphi_i(x) : \nabla^2\varphi_j(x) dx
 * @f} ，
 * 以及对右侧向量的贡献
 * @f{align*}{
 * f^K_i = \int_K \varphi_i(x) f(x) dx
 * @f}
 * 

 * 
 * 我们使用与组装 step-22 相同的技术来加速该函数。我们不在最里面的循环中调用`fe_values.shape_hessian(i, qpoint)`，而是创建一个变量`hessian_i`，在循环中对`i`进行一次评估，在循环中对`j`重新使用如此评估的值。为了对称，我们对变量`hessian_j`也做了同样的处理，尽管它确实只用了一次，而且我们可以在计算两个项之间标量乘积的指令中留下对`fe_values.shape_hessian(j,qpoint)`的调用。
 * 

 * 
 * 
 * @code
 *     auto cell_worker = [&](const Iterator &  cell, 
 *                            ScratchData<dim> &scratch_data, 
 *                            CopyData &        copy_data) { 
 *       copy_data.cell_matrix = 0; 
 *       copy_data.cell_rhs    = 0; 
 * 
 *       FEValues<dim> &fe_values = scratch_data.fe_values; 
 *       fe_values.reinit(cell); 
 * 
 *       cell->get_dof_indices(copy_data.local_dof_indices); 
 * 
 *       const ExactSolution::RightHandSide<dim> right_hand_side; 
 * 
 *       const unsigned int dofs_per_cell = 
 *         scratch_data.fe_values.get_fe().n_dofs_per_cell(); 
 * 
 *       for (unsigned int qpoint = 0; qpoint < fe_values.n_quadrature_points; 
 *            ++qpoint) 
 *         { 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               const Tensor<2, dim> &hessian_i = 
 *                 fe_values.shape_hessian(i, qpoint); 
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   const Tensor<2, dim> &hessian_j = 
 *                     fe_values.shape_hessian(j, qpoint); 
 * 
 *                   copy_data.cell_matrix(i, j) += 
 *                     scalar_product(hessian_i,   // nabla^2 phi_i(x) 
 *                                    hessian_j) * // nabla^2 phi_j(x) 
 *                     fe_values.JxW(qpoint);      // dx 
 *                 } 
 * 
 *               copy_data.cell_rhs(i) += 
 *                 fe_values.shape_value(i, qpoint) * // phi_i(x) 
 *                 right_hand_side.value( 
 *                   fe_values.quadrature_point(qpoint)) * // f(x) 
 *                 fe_values.JxW(qpoint);                  // dx 
 *             } 
 *         } 
 *     }; 
 * 
 * @endcode
 * 
 * 下一个构建模块是在网格的每个内部面组装惩罚项。正如 MeshWorker::mesh_loop(), 文档中所描述的，这个函数接收到的参数表示一个单元和它的相邻单元，以及（对于这两个单元中的每一个）我们必须整合的面（以及潜在的子面）。同样地，我们也得到了一个从头开始的对象，以及一个用于放置结果的拷贝对象。
 * 

 * 
 * 这个函数本身有三个部分。在顶部，我们初始化FEInterfaceValues对象，并创建一个新的 `CopyData::FaceData` 对象来存储我们的输入。这将被推到`copy_data.face_data`变量的末尾。我们需要这样做，因为我们对一个给定单元进行积分的面（或子面）的数量因单元而异，而且这些矩阵的大小也不同，取决于面或子面相邻的自由度。正如 MeshWorker::mesh_loop(), 文档中所讨论的，每次访问一个新的单元时，复制对象都会被重置，所以我们推到`copy_data.face_data()`末尾的内容实际上就是后来的`copier`函数在复制每个单元的贡献到全局矩阵和右侧对象时所能看到的。
 * 

 * 
 * 
 * @code
 *     auto face_worker = [&](const Iterator &    cell, 
 *                            const unsigned int &f, 
 *                            const unsigned int &sf, 
 *                            const Iterator &    ncell, 
 *                            const unsigned int &nf, 
 *                            const unsigned int &nsf, 
 *                            ScratchData<dim> &  scratch_data, 
 *                            CopyData &          copy_data) { 
 *       FEInterfaceValues<dim> &fe_interface_values = 
 *         scratch_data.fe_interface_values; 
 *       fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf); 
 * 
 *       copy_data.face_data.emplace_back(); 
 *       CopyData::FaceData &copy_data_face = copy_data.face_data.back(); 
 * 
 *       copy_data_face.joint_dof_indices = 
 *         fe_interface_values.get_interface_dof_indices(); 
 * 
 *       const unsigned int n_interface_dofs = 
 *         fe_interface_values.n_current_interface_dofs(); 
 *       copy_data_face.cell_matrix.reinit(n_interface_dofs, n_interface_dofs); 
 * 
 * @endcode
 * 
 * 第二部分涉及到确定惩罚参数应该是什么。通过观察双线性形式中各种项的单位，很明显，惩罚必须具有 $\frac{\gamma}{h_K}$ 的形式（即，超过长度尺度的一个），但如何选择无维数 $\gamma$ 并不是先验的。从拉普拉斯方程的不连续Galerkin理论来看，人们可能猜想正确的选择是 $\gamma=p(p+1)$ 是正确的选择，其中 $p$ 是所用有限元的多项式程度。我们将在本程序的结果部分更详细地讨论这个选择。
 * 

 * 
 * 在上面的公式中， $h_K$  是单元格  $K$  的大小。但这也不是很简单的事情。如果使用高度拉伸的单元格，那么一个更复杂的理论说， $h$ 应该被单元格 $K$ 的直径取代，该直径是有关边缘方向的法线。 事实证明，在deal.II中有一个函数用于此。其次，当从一个面的两个不同侧面看时， $h_K$ 可能是不同的。
 * 

 * 
 * 为了安全起见，我们取这两个值的最大值。我们将注意到，如果使用自适应网格细化所产生的悬空节点，有可能需要进一步调整这一计算方法。
 * 

 * 
 * 
 * @code
 *       const unsigned int p = fe.degree; 
 *       const double       gamma_over_h = 
 *         std::max((1.0 * p * (p + 1) / 
 *                   cell->extent_in_direction( 
 *                     GeometryInfo<dim>::unit_normal_direction[f])), 
 *                  (1.0 * p * (p + 1) / 
 *                   ncell->extent_in_direction( 
 *                     GeometryInfo<dim>::unit_normal_direction[nf]))); 
 * 
 * @endcode
 * 
 * 最后，像往常一样，我们在正交点和指数`i`和`j`上循环，把这个面或子面的贡献加起来。然后将这些数据存储在上面创建的`copy_data.face_data`对象中。至于单元格工作者，如果可能的话，我们将平均数和跳跃的评估从循环中拉出来，引入局部变量来存储这些结果。然后组件只需要在最里面的循环中使用这些局部变量。关于这段代码实现的具体公式，回顾一下，双线性形式的接口项如下。
 * @f{align*}{
 * -\sum_{e \in \mathbb{F}} \int_{e}
 * \jump{ \frac{\partial v_h}{\partial \mathbf n}}
 * \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds
 * -\sum_{e \in \mathbb{F}} \int_{e}
 * \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
 * \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
 * + \sum_{e \in \mathbb{F}}
 * \frac{\gamma}{h_e}
 * \int_e
 * \jump{\frac{\partial v_h}{\partial \mathbf n}}
 * \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds.
 * @f}
 * 

 * 
 * 
 * @code
 *       for (unsigned int qpoint = 0; 
 *            qpoint < fe_interface_values.n_quadrature_points; 
 *            ++qpoint) 
 *         { 
 *           const auto &n = fe_interface_values.normal(qpoint); 
 * 
 *           for (unsigned int i = 0; i < n_interface_dofs; ++i) 
 *             { 
 *               const double av_hessian_i_dot_n_dot_n = 
 *                 (fe_interface_values.average_hessian(i, qpoint) * n * n); 
 *               const double jump_grad_i_dot_n = 
 *                 (fe_interface_values.jump_gradient(i, qpoint) * n); 
 * 
 *               for (unsigned int j = 0; j < n_interface_dofs; ++j) 
 *                 { 
 *                   const double av_hessian_j_dot_n_dot_n = 
 *                     (fe_interface_values.average_hessian(j, qpoint) * n * n); 
 *                   const double jump_grad_j_dot_n = 
 *                     (fe_interface_values.jump_gradient(j, qpoint) * n); 
 * 
 *                   copy_data_face.cell_matrix(i, j) += 
 *                     (-av_hessian_i_dot_n_dot_n       // - {grad^2 v n n } 
 *                        * jump_grad_j_dot_n           // [grad u n] 
 *                      - av_hessian_j_dot_n_dot_n      // - {grad^2 u n n } 
 *                          * jump_grad_i_dot_n         // [grad v n] 
 *                      +                               // + 
 *                      gamma_over_h *                  // gamma/h 
 *                        jump_grad_i_dot_n *           // [grad v n] 
 *                        jump_grad_j_dot_n) *          // [grad u n] 
 *                     fe_interface_values.JxW(qpoint); // dx 
 *                 } 
 *             } 
 *         } 
 *     }; 
 * 
 * @endcode
 * 
 * 第三块是对处于边界的面做同样的装配。当然，想法和上面一样，唯一不同的是，现在有惩罚条款也进入了右手边。
 * 

 * 
 * 和以前一样，这个函数的第一部分只是设置了一些辅助对象。
 * 

 * 
 * 
 * @code
 *     auto boundary_worker = [&](const Iterator &    cell, 
 *                                const unsigned int &face_no, 
 *                                ScratchData<dim> &  scratch_data, 
 *                                CopyData &          copy_data) { 
 *       FEInterfaceValues<dim> &fe_interface_values = 
 *         scratch_data.fe_interface_values; 
 *       fe_interface_values.reinit(cell, face_no); 
 *       const auto &q_points = fe_interface_values.get_quadrature_points(); 
 * 
 *       copy_data.face_data.emplace_back(); 
 *       CopyData::FaceData &copy_data_face = copy_data.face_data.back(); 
 * 
 *       const unsigned int n_dofs = 
 *         fe_interface_values.n_current_interface_dofs(); 
 *       copy_data_face.joint_dof_indices = 
 *         fe_interface_values.get_interface_dof_indices(); 
 * 
 *       copy_data_face.cell_matrix.reinit(n_dofs, n_dofs); 
 * 
 *       const std::vector<double> &JxW = fe_interface_values.get_JxW_values(); 
 *       const std::vector<Tensor<1, dim>> &normals = 
 *         fe_interface_values.get_normal_vectors(); 
 * 
 *       const ExactSolution::Solution<dim> exact_solution; 
 *       std::vector<Tensor<1, dim>>        exact_gradients(q_points.size()); 
 *       exact_solution.gradient_list(q_points, exact_gradients); 
 * 
 * @endcode
 * 
 * 从正面看，由于我们现在只处理与面相邻的一个单元（因为我们在边界上），惩罚因子 $\gamma$ 的计算大大简化了。
 * 

 * 
 * 
 * @code
 *       const unsigned int p = fe.degree; 
 *       const double       gamma_over_h = 
 *         (1.0 * p * (p + 1) / 
 *          cell->extent_in_direction( 
 *            GeometryInfo<dim>::unit_normal_direction[face_no])); 
 * 
 * @endcode
 * 
 * 第三块是术语的组合。由于这些条款包含了矩阵的条款和右手边的条款，所以现在稍微有些麻烦。前者与上面所说的内部面完全相同，如果我们只是适当地定义了跳跃和平均（这就是FEInterfaceValues类所做的）。后者需要我们评估边界条件 $j(\mathbf x)$ ，在当前情况下（我们知道确切的解决方案），我们从 $j(\mathbf x) = \frac{\partial u(\mathbf x)}{\partial {\mathbf n}}$ 中计算出来。然后，要添加到右侧向量的项是  $\frac{\gamma}{h_e}\int_e \jump{\frac{\partial v_h}{\partial \mathbf n}} j \ ds$  。
 * 

 * 
 * 
 * @code
 *       for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint) 
 *         { 
 *           const auto &n = normals[qpoint]; 
 * 
 *           for (unsigned int i = 0; i < n_dofs; ++i) 
 *             { 
 *               const double av_hessian_i_dot_n_dot_n = 
 *                 (fe_interface_values.average_hessian(i, qpoint) * n * n); 
 *               const double jump_grad_i_dot_n = 
 *                 (fe_interface_values.jump_gradient(i, qpoint) * n); 
 * 
 *               for (unsigned int j = 0; j < n_dofs; ++j) 
 *                 { 
 *                   const double av_hessian_j_dot_n_dot_n = 
 *                     (fe_interface_values.average_hessian(j, qpoint) * n * n); 
 *                   const double jump_grad_j_dot_n = 
 *                     (fe_interface_values.jump_gradient(j, qpoint) * n); 
 * 
 *                   copy_data_face.cell_matrix(i, j) += 
 *                     (-av_hessian_i_dot_n_dot_n  // - {grad^2 v n n} 
 *                        * jump_grad_j_dot_n      //   [grad u n] 
 * 
 * @endcode
 * 
 * 
 * 
 * 

 * 
 * 
 * @code
 *                      - av_hessian_j_dot_n_dot_n // - {grad^2 u n n} 
 *                          * jump_grad_i_dot_n    //   [grad v n] 
 * 
 * @endcode
 * 
 * 
 * 
 * 

 * 
 * 
 * @code
 *                      + gamma_over_h             //  gamma/h 
 *                          * jump_grad_i_dot_n    // [grad v n] 
 *                          * jump_grad_j_dot_n    // [grad u n] 
 *                      ) * 
 *                     JxW[qpoint]; // dx 
 *                 } 
 * 
 *               copy_data.cell_rhs(i) += 
 *                 (-av_hessian_i_dot_n_dot_n *       // - {grad^2 v n n } 
 *                    (exact_gradients[qpoint] * n)   //   (grad u_exact . n) 
 *                  +                                 // + 
 *                  gamma_over_h                      //  gamma/h 
 *                    * jump_grad_i_dot_n             // [grad v n] 
 *                    * (exact_gradients[qpoint] * n) // (grad u_exact . n) 
 *                  ) * 
 *                 JxW[qpoint]; // dx 
 *             } 
 *         } 
 *     }; 
 * 
 * @endcode
 * 
 * 第四部分是一个小函数，它将上面的单元格、内部和边界面装配程序产生的数据复制到全局矩阵和右手向量中。这里真的没有什么可做的。我们分配单元格矩阵和右侧贡献，就像我们在其他几乎所有的教程程序中使用约束对象那样。然后，我们还必须对面矩阵的贡献做同样的处理，这些贡献已经获得了面（内部和边界）的内容，并且`面_工作`和`边界_工作`已经添加到`copy_data.face_data`阵列中。
 * 

 * 
 * 
 * @code
 *     auto copier = [&](const CopyData &copy_data) { 
 *       constraints.distribute_local_to_global(copy_data.cell_matrix, 
 *                                              copy_data.cell_rhs, 
 *                                              copy_data.local_dof_indices, 
 *                                              system_matrix, 
 *                                              system_rhs); 
 * 
 *       for (auto &cdf : copy_data.face_data) 
 *         { 
 *           constraints.distribute_local_to_global(cdf.cell_matrix, 
 *                                                  cdf.joint_dof_indices, 
 *                                                  system_matrix); 
 *         } 
 *     }; 
 * 
 * @endcode
 * 
 * 在设置了所有这些之后，剩下的就是创建一个从头开始和复制数据的对象，并调用 MeshWorker::mesh_loop() 函数，然后遍历所有的单元格和面，调用它们各自的工作器，然后是复制器函数，将东西放入全局矩阵和右侧。作为一个额外的好处， MeshWorker::mesh_loop() 以并行方式完成所有这些工作，使用你的机器恰好有多少个处理器核心。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1; 
 *     ScratchData<dim>   scratch_data(mapping, 
 *                                   fe, 
 *                                   n_gauss_points, 
 *                                   update_values | update_gradients | 
 *                                     update_hessians | update_quadrature_points | 
 *                                     update_JxW_values, 
 *                                   update_values | update_gradients | 
 *                                     update_hessians | update_quadrature_points | 
 *                                     update_JxW_values | update_normal_vectors); 
 *     CopyData           copy_data(dof_handler.get_fe().n_dofs_per_cell()); 
 *     MeshWorker::mesh_loop(dof_handler.begin_active(), 
 *                           dof_handler.end(), 
 *                           cell_worker, 
 *                           copier, 
 *                           scratch_data, 
 *                           copy_data, 
 *                           MeshWorker::assemble_own_cells | 
 *                             MeshWorker::assemble_boundary_faces | 
 *                             MeshWorker::assemble_own_interior_faces_once, 
 *                           boundary_worker, 
 *                           face_worker); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Solvingthelinearsystemandpostprocessing"></a> 
 * <h4>Solving the linear system and postprocessing</h4>
 * 

 * 
 * 到此为止，节目基本上结束了。其余的函数并不太有趣或新颖。第一个函数只是用一个直接求解器来求解线性系统（也见 step-29  ）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BiharmonicProblem<dim>::solve() 
 *   { 
 *     std::cout << "   Solving system..." << std::endl; 
 * 
 *     SparseDirectUMFPACK A_direct; 
 *     A_direct.initialize(system_matrix); 
 *     A_direct.vmult(solution, system_rhs); 
 * 
 *     constraints.distribute(solution); 
 *   } 
 * 
 * @endcode
 * 
 * 下一个函数评估了计算出的解和精确解之间的误差（在这里是已知的，因为我们选择了右手边和边界值的方式，所以我们知道相应的解）。在下面的前两个代码块中，我们计算了 $L_2$ 准则和 $H^1$ 半准则下的误差。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BiharmonicProblem<dim>::compute_errors() 
 *   { 
 *     { 
 *       Vector<float> norm_per_cell(triangulation.n_active_cells()); 
 *       VectorTools::integrate_difference(mapping, 
 *                                         dof_handler, 
 *                                         solution, 
 *                                         ExactSolution::Solution<dim>(), 
 *                                         norm_per_cell, 
 *                                         QGauss<dim>(fe.degree + 2), 
 *                                         VectorTools::L2_norm); 
 *       const double error_norm = 
 *         VectorTools::compute_global_error(triangulation, 
 *                                           norm_per_cell, 
 *                                           VectorTools::L2_norm); 
 *       std::cout << "   Error in the L2 norm           :     " << error_norm 
 *                 << std::endl; 
 *     } 
 * 
 *     { 
 *       Vector<float> norm_per_cell(triangulation.n_active_cells()); 
 *       VectorTools::integrate_difference(mapping, 
 *                                         dof_handler, 
 *                                         solution, 
 *                                         ExactSolution::Solution<dim>(), 
 *                                         norm_per_cell, 
 *                                         QGauss<dim>(fe.degree + 2), 
 *                                         VectorTools::H1_seminorm); 
 *       const double error_norm = 
 *         VectorTools::compute_global_error(triangulation, 
 *                                           norm_per_cell, 
 *                                           VectorTools::H1_seminorm); 
 *       std::cout << "   Error in the H1 seminorm       : " << error_norm 
 *                 << std::endl; 
 *     } 
 * 
 * @endcode
 * 
 * 现在也计算一下 $H^2$ 半正态误差的近似值。实际的 $H^2$ 半规范要求我们对解决方案 $u_h$ 的二阶导数进行积分，但是考虑到我们使用的拉格朗日形状函数， $u_h$ 当然在单元间的界面上有结点，因此二阶导数在界面是奇异的。因此，我们实际上只对单元的内部进行积分，而忽略了界面的贡献。这不是*等同于问题的能量准则，但是仍然可以让我们了解误差收敛的速度。
 * 

 * 
 * 我们注意到，我们可以通过定义一个等同于能量准则的准则来解决这个问题。这将涉及到不仅要像我们下面做的那样将细胞内部的积分相加，而且还要为 $u_h$ 的导数在界面上的跳跃添加惩罚项，并对这两种项进行适当的缩放。我们将把这个问题留给以后的工作。
 * 

 * 
 * 
 * @code
 *     { 
 *       const QGauss<dim>            quadrature_formula(fe.degree + 2); 
 *       ExactSolution::Solution<dim> exact_solution; 
 *       Vector<double> error_per_cell(triangulation.n_active_cells()); 
 * 
 *       FEValues<dim> fe_values(mapping, 
 *                               fe, 
 *                               quadrature_formula, 
 *                               update_values | update_hessians | 
 *                                 update_quadrature_points | update_JxW_values); 
 * 
 *       FEValuesExtractors::Scalar scalar(0); 
 *       const unsigned int         n_q_points = quadrature_formula.size(); 
 * 
 *       std::vector<SymmetricTensor<2, dim>> exact_hessians(n_q_points); 
 *       std::vector<Tensor<2, dim>>          hessians(n_q_points); 
 *       for (auto &cell : dof_handler.active_cell_iterators()) 
 *         { 
 *           fe_values.reinit(cell); 
 *           fe_values[scalar].get_function_hessians(solution, hessians); 
 *           exact_solution.hessian_list(fe_values.get_quadrature_points(), 
 *                                       exact_hessians); 
 * 
 *           double local_error = 0; 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *             { 
 *               local_error += 
 *                 ((exact_hessians[q_point] - hessians[q_point]).norm_square() * 
 *                  fe_values.JxW(q_point)); 
 *             } 
 *           error_per_cell[cell->active_cell_index()] = std::sqrt(local_error); 
 *         } 
 * 
 *       const double error_norm = error_per_cell.l2_norm(); 
 *       std::cout << "   Error in the broken H2 seminorm: " << error_norm 
 *                 << std::endl; 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 同样无趣的是生成图形输出的函数。它看起来和  step-6  中的一模一样，比如说。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   BiharmonicProblem<dim>::output_results(const unsigned int iteration) const 
 *   { 
 *     std::cout << "   Writing graphical output..." << std::endl; 
 * 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "solution"); 
 *     data_out.build_patches(); 
 * 
 *     const std::string filename = 
 *       ("output_" + Utilities::int_to_string(iteration, 6) + ".vtu"); 
 *     std::ofstream output_vtu(filename); 
 *     data_out.write_vtu(output_vtu); 
 *   } 
 * 
 * @endcode
 * 
 * `run()`函数的情况也是如此。就像在以前的程序中一样。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BiharmonicProblem<dim>::run() 
 *   { 
 *     make_grid(); 
 * 
 *     const unsigned int n_cycles = 4; 
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
 *       { 
 *         std::cout << "Cycle " << cycle << " of " << n_cycles << std::endl; 
 * 
 *         triangulation.refine_global(1); 
 *         setup_system(); 
 * 
 *         assemble_system(); 
 *         solve(); 
 * 
 *         output_results(cycle); 
 * 
 *         compute_errors(); 
 *         std::cout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step47 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * 最后是 "main() "函数。同样，这里没有什么可看的。它看起来和以前的教程程序中的一样。有一个变量，可以选择我们要用来解方程的元素的多项式程度。因为我们使用的C0IP公式要求元素的度数至少为2，所以我们用一个断言来检查，无论为多项式度数设置什么都是有意义的。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step47; 
 * 
 *       const unsigned int fe_degree = 2; 
 *       Assert(fe_degree >= 2, 
 *              ExcMessage("The C0IP formulation for the biharmonic problem " 
 *                         "only works if one uses elements of polynomial " 
 *                         "degree at least 2.")); 
 * 
 *       BiharmonicProblem<2> biharmonic_problem(fe_degree); 
 *       biharmonic_problem.run(); 
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
 * 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-47/doc/results.dox



<a name="Results"></a><h1>Results</h1>


我们用介绍中讨论的右手边和边界值运行程序。这些将产生域 $\Omega = (0,1)^2$ 上的解 $u = \sin(\pi x) \sin(\pi y)$ 。我们用 $Q_2$ 、 $Q_3$ 和 $Q_4$ 元素来测试这个设置，我们可以通过`main()`函数中的`fe_degree`变量来改变。通过网格细化， $L_2$ 的收敛率、 $H^1$ 的近似值率和 $H^2$ 的近似值收敛率对于 $u$ 应该分别为2、2、1（如介绍中所述， $L_2$ 的规范为次优）；对于 $Q_3$ 为4、3、2；而对于 $Q_4$ 为5、4、3。

从文献来看，并不立即清楚惩罚参数 $\gamma$ 应该是什么。例如， @cite Brenner2009 指出它需要大于1，并选择 $\gamma=5$  。FEniCS/Dolphin教程选择它为 $\gamma=8$  ，见https://fenicsproject.org/docs/dolfin/1.6.0/python/demo/documented/biharmonic/python/documentation.html 。   @cite Wells2007 使用的 $\gamma$ 值大于Kirchhoff板的元素所属的边数（见他们的第4.2节）。这表明也许 $\gamma = 1$ ,  $2$ , 太小了；另一方面， $p(p+1)$ 的值也是合理的，其中 $p$ 是多项式的度数。通过与拉普拉斯方程的不连续Galerkin公式相比较，人们期望最后一个选择是可行的（例如，见步骤39和步骤74的讨论），而且它在这里也将证明是可行的。但是我们应该检查 $\gamma$ 的哪个值是正确的，下面我们将这样做；改变 $\gamma$ 在`assemble_system()`中定义的两个`face_worker`和`boundary_worker`函数中很容易。




<a name="TestresultsoniQsub2subiiQsub2subiwithigammapp1iigammapp1i"></a><h3>Test results on <i>Q<sub>2</sub></i><i>Q<sub>2</sub></i> with <i>&gamma; = p(p+1)</i><i>&gamma; = p(p+1)</i> </h3>


我们用不同的细化网格运行代码，得到以下收敛率。

 <table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>   8.780e-03 </td><td>       </td><td>  7.095e-02   </td><td>           </td><td>  1.645 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   3.515e-03   </td><td>  1.32 </td><td> 2.174e-02  </td><td>     1.70     </td><td> 8.121e-01  </td><td>  1.018  </td>
  </tr>
  <tr>
   <td>   4                  </td><td>   1.103e-03   </td><td>  1.67   </td><td> 6.106e-03    </td><td>  1.83        </td><td>   4.015e-01 </td><td> 1.016  </td>
  </tr>
  <tr>
   <td>   5                  </td><td>  3.084e-04  </td><td>  1.83   </td><td>  1.622e-03   </td><td>    1.91        </td><td> 1.993e-01 </td><td>  1.010   </td>
  </tr>
</table>  我们可以看到， $L_2$ 的收敛率在2左右， $H^1$  -seminorm收敛率在2左右， $H^2$  -seminorm收敛率在1左右。后两者与理论上的预期收敛率相符；对于前者，我们没有定理，但鉴于介绍中的评论，对于它是次优的也不奇怪。




<a name="TestresultsoniQsub3subiiQsub3subiwithigammapp1iigammapp1i"></a><h3>Test results on <i>Q<sub>3</sub></i><i>Q<sub>3</sub></i> with <i>&gamma; = p(p+1)</i><i>&gamma; = p(p+1)</i> </h3>



 <table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>    2.045e-04 </td><td>       </td><td>   4.402e-03   </td><td>           </td><td> 1.641e-01 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   1.312e-05   </td><td> 3.96  </td><td>  5.537e-04  </td><td>   2.99     </td><td> 4.096e-02 </td><td>  2.00  </td>
  </tr>
  <tr>
   <td>   4                  </td><td>   8.239e-07 </td><td>  3.99  </td><td> 6.904e-05   </td><td> 3.00     </td><td> 1.023e-02 </td><td> 2.00 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>   5.158e-08  </td><td>  3.99 </td><td> 8.621e-06 </td><td>  3.00      </td><td> 2.558e-03  </td><td>  2.00  </td>
  </tr>
</table>  我们可以看到， $L_2$  收敛率在4左右， $H^1$  -seminorm 收敛率在3左右， $H^2$  -seminorm 收敛率在2左右。当然，这符合我们的理论预期。




<a name="TestresultsoniQsub4subiiQsub4subiwithigammapp1iigammapp1i"></a><h3>Test results on <i>Q<sub>4</sub></i><i>Q<sub>4</sub></i> with <i>&gamma; = p(p+1)</i><i>&gamma; = p(p+1)</i> </h3>


 <table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>    6.510e-06 </td><td>       </td><td> 2.215e-04   </td><td>           </td><td>  1.275e-02 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   2.679e-07  </td><td>  4.60  </td><td> 1.569e-05  </td><td>   3.81    </td><td> 1.496e-03 </td><td>  3.09  </td>
  </tr>
  <tr>
   <td>   4                  </td><td>   9.404e-09  </td><td> 4.83   </td><td> 1.040e-06    </td><td> 3.91       </td><td> 1.774e-04 </td><td> 3.07 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>   7.943e-10 </td><td>  3.56  </td><td>   6.693e-08 </td><td> 3.95     </td><td> 2.150e-05  </td><td> 3.04    </td>
  </tr>
</table>  我们可以看到， $L_2$  norm收敛率在5左右， $H^1$  -seminorm收敛率在4左右，而 $H^2$  -seminorm收敛率在3左右。在最细的网格上， $L_2$ 规范收敛率比我们的理论预期小得多，因为线性求解器由于舍入而成为限制因素。当然在这种情况下， $L_2$ 误差也已经非常小了。




<a name="TestresultsoniQsub2subiiQsub2subiwithigamma1iigamma1i"></a><h3>Test results on <i>Q<sub>2</sub></i><i>Q<sub>2</sub></i> with <i>&gamma; = 1</i><i>&gamma; = 1</i> </h3>


为了与上述结果进行比较，现在让我们也考虑一下我们简单地选择 $\gamma=1$ 的情况。

 <table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>   7.350e-02 </td><td>       </td><td>   7.323e-01   </td><td>           </td><td> 10.343 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>   6.798e-03   </td><td> 3.43  </td><td> 1.716e-01   </td><td>   2.09    </td><td>4.836 </td><td>  1.09 </td>
  </tr>
  <tr>
   <td>   4                  </td><td>  9.669e-04   </td><td> 2.81   </td><td> 6.436e-02    </td><td> 1.41      </td><td>  3.590 </td><td> 0.430 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>   1.755e-04 </td><td> 2.46 </td><td>  2.831e-02  </td><td>    1.18      </td><td>3.144  </td><td>  0.19  </td>
  </tr>
</table>  虽然 $L_2$ 规范的收敛率 $u$ 或多或少符合理论预期，但 $H^1$  -seminorm和 $H^2$  -seminorm似乎并没有像预期那样收敛。比较 $\gamma = 1$ 和 $\gamma = p(p+1)$ 的结果，很明显， $\gamma = p(p+1)$ 是一个更好的惩罚。鉴于 $\gamma=1$ 对于 $Q_2$ 元素来说已经太小了，如果用 $Q_3$ 元素重复实验，结果更加令人失望，这可能就不奇怪了。我们再次只得到2、1、0的收敛率--也就是说，不比 $Q_2$ 元素好（尽管误差的大小更小）。然而，也许令人惊讶的是，当使用 $Q_4$ 元素时，人们获得了或多或少的预期收敛顺序。无论如何，这种不确定性表明， $\gamma=1$ 充其量是一个有风险的选择，最坏的情况是一个不可靠的选择，我们应该选择 $\gamma$ 更大。




<a name="TestresultsoniQsub2subiiQsub2subiwithigamma2iigamma2i"></a><h3>Test results on <i>Q<sub>2</sub></i><i>Q<sub>2</sub></i> with <i>&gamma; = 2</i><i>&gamma; = 2</i> </h3>


由于 $\gamma=1$ 显然太小了，人们可能猜想 $\gamma=2$ 实际上可能效果更好。下面是在这种情况下得到的结果。

 <table align="center" class="doxtable">
  <tr>
   <th>Number of refinements </th><th>  $\|u-u_h^\circ\|_{L_2}$ </th><th>  Conv. rates  </th><th>  $|u-u_h|_{H^1}$ </th><th> Conv. rates </th><th> $|u-u_h|_{H^2}$ </th><th> Conv. rates </th>
  </tr>
  <tr>
   <td>   2                  </td><td>   4.133e-02 </td><td>       </td><td>  2.517e-01   </td><td>           </td><td> 3.056 </td><td>   </td>
  </tr>
  <tr>
   <td>   3                  </td><td>  6.500e-03   </td><td>2.66  </td><td> 5.916e-02  </td><td>  2.08    </td><td>1.444 </td><td>  1.08 </td>
  </tr>
  <tr>
   <td>   4                  </td><td> 6.780e-04   </td><td> 3.26  </td><td> 1.203e-02    </td><td> 2.296      </td><td> 6.151e-01 </td><td> 1.231 </td>
  </tr>
  <tr>
   <td>   5                  </td><td> 1.622e-04 </td><td> 2.06 </td><td>  2.448e-03  </td><td>   2.297     </td><td> 2.618e-01  </td><td> 1.232  </td>
  </tr>
</table>  在这种情况下，收敛率或多或少符合理论预期，但与 $\gamma =
p(p+1)$ 的结果相比，变化更大。同样，我们可以对 $Q_3$ 和 $Q_4$ 元素重复这种实验。在这两种情况下，我们会发现我们获得了大致的预期收敛率。那么，更感兴趣的可能是比较误差的绝对大小。在上表中，对于 $Q_2$ 情况，最细网格上的误差在 $\gamma=p(p+1)$ 和 $\gamma=2$ 情况下是相当的，而对于 $Q_3$ ， $\gamma=2$ 的误差要比 $\gamma=p(p+1)$ 的大很多。对于 $Q_4$ 的情况也是如此。




<a name="Conclusionsforthechoiceofthepenaltyparameter"></a><h3> Conclusions for the choice of the penalty parameter </h3>


关于应该使用哪种 "合理 "的惩罚参数的结论是， $\gamma=p(p+1)$ 产生了预期的结果。因此，它是目前编写的代码所使用的。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


这个方案有一些明显的扩展，会有意义。

- 该程序使用一个正方形域和一个均匀的网格。真正的问题不是这样的，我们应该验证在其他形状的域上的收敛性，特别是在弯曲的边界上。人们也可能对使用自适应网格细化来解决规则性较差的区域感兴趣。

- 从更多的理论角度来看，上面的收敛结果只使用了 "破损的" $H^2$ 半规范 $|\cdot|^\circ_{H^2}$ ，而不是 "等同的 "规范 $|\cdot|_h$  。这足以让我们相信，这个程序并没有从根本上被破坏。然而，测量我们有理论结果的实际规范的误差可能是有趣的。例如，使用FEInterfaceValues类与 MeshWorker::mesh_loop() 结合，实现这一补充应该不会太困难，其精神与我们用于组装线性系统的精神相同。


<a name="Derivationforthesimplysupportedplates"></a>  <h4> Derivation for the simply supported plates </h4>


  类似于实施中所涉及的 "夹持 "边界条件，我们将推导出 $C^0$ IP有限元方案，用于简单支撑板。   @f{align*}{
    \Delta^2 u(\mathbf x) &= f(\mathbf x)
    \qquad \qquad &&\forall \mathbf x \in \Omega,
    u(\mathbf x) &= g(\mathbf x) \qquad \qquad
    &&\forall \mathbf x \in \partial\Omega, \\
    \Delta u(\mathbf x) &= h(\mathbf x) \qquad \qquad
    &&\forall \mathbf x \in \partial\Omega.
  @f}

  我们用测试函数 $v_h$ 乘以双调方程，并对 $ K $ 进行积分，得到。   @f{align*}{
    \int_K v_h (\Delta^2 u_h)
     &= \int_K (D^2 v_h) : (D^2 u_h)
       + \int_{\partial K} v_h \frac{\partial (\Delta u_h)}{\partial \mathbf{n}}


       -\int_{\partial K} (\nabla v_h) \cdot (\frac{\partial \nabla u_h}{\partial \mathbf{n}}).
  @f}



  将所有单元格 $K \in  \mathbb{T}$ 相加，因为 $\Delta u_h$ 的法线方向在两条单元格共享的每条内边上指向相反的方向， $v_h = 0$ 在 $\partial \Omega$ 上，@f{align*}{
  \sum_{K \in \mathbb{T}} \int_{\partial K} v_h \frac{\partial (\Delta u_h)}{\partial \mathbf{n}} = 0,
  @f}

  并通过细胞界面上的跳跃定义，@f{align*}{


  -\sum_{K \in \mathbb{T}} \int_{\partial K} (\nabla v_h) \cdot (\frac{\partial \nabla u_h}{\partial \mathbf{n}}) = -\sum_{e \in \mathbb{F}} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} (\frac{\partial^2 u_h}{\partial \mathbf{n^2}}).
  @f} 。

  我们将域的内部面和边界面分开，@f{align*}{


  -\sum_{K \in \mathbb{T}} \int_{\partial K} (\nabla v_h) \cdot (\frac{\partial \nabla u_h}{\partial \mathbf{n}}) = -\sum_{e \in \mathbb{F}^i} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} (\frac{\partial^2 u_h}{\partial \mathbf{n^2}})


  - \sum_{e \in \partial \Omega} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} h,
  @f} 。

  其中 $\mathbb{F}^i$ 是内部面的集合。   这使我们得出@f{align*}{
  \sum_{K \in \mathbb{T}} \int_K (D^2 v_h) : (D^2 u_h) \ dx - \sum_{e \in \mathbb{F}^i} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} (\frac{\partial^2 u_h}{\partial \mathbf{n^2}}) \ ds
  = \sum_{K \in \mathbb{T}}\int_{K} v_h f  \ dx + \sum_{e\subset\partial\Omega} \int_{e} \jump{\frac{\partial v_h}{\partial \mathbf{n}}} h \ ds.
  @f}。



  为了使离散问题对称化和稳定化，我们加入了对称化和稳定化项。   我们最终得到双调方程的 $C^0$ IP有限元方案：找到 $u_h$ ，使 $u_h =g$ 上的 $\partial \Omega$ 和@f{align*}{
  \mathcal{A}(v_h,u_h)&=\mathcal{F}(v_h) \quad \text{holds for all test functions } v_h,
  @f}

  其中@f{align*}{
  \mathcal{A}(v_h,u_h):=&\sum_{K \in \mathbb{T}}\int_K D^2v_h:D^2u_h \ dx
  \\
  &


   -\sum_{e \in \mathbb{F}^i} \int_{e}
    \jump{\frac{\partial v_h}{\partial \mathbf n}}
    \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds


   -\sum_{e \in \mathbb{F}^i} \int_{e}
   \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
   \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
  \\
  &+ \sum_{e \in \mathbb{F}^i}
   \frac{\gamma}{h_e}
   \int_e
   \jump{\frac{\partial v_h}{\partial \mathbf n}}
   \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds,
  @f}

  和@f{align*}{
  \mathcal{F}(v_h)&:=\sum_{K \in \mathbb{T}}\int_{K} v_h f \ dx
  +
  \sum_{e\subset\partial\Omega}
  \int_e \jump{\frac{\partial v_h}{\partial \mathbf n}} h \ ds.
  @f}

  这个边界案例的实现与 "钳制 "版本类似，只是在系统组装时不再需要`边界_工人'，并且根据配方改变右手边。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-47.cc"
*/
