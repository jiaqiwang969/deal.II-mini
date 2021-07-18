//include/deal.II-translator/numerics/error_estimator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_error_estimator_h
#define dealii_error_estimator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <deal.II/fe/component_mask.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int, int>
class DoFHandler;
template <int, int>
class Mapping;
template <int>
class Quadrature;

namespace hp
{
  template <int>
  class QCollection;
}
#endif


/**
 * 实现Kelly, De S. R. Gago, Zienkiewicz和Babuska的误差指标（见
 * @cite KGZB83
 * ）及其对hp-FEM的修改。这个误差指标试图通过整合解的梯度沿每个单元的面的跳跃来近似每个单元的误差。
 * 它可以被理解为梯度恢复估计器；完整的讨论见Ainsworth和Oden的调查报告《有限元分析中的后验误差估计》（Wiley，2000）。
 * 在最初的Kelly误差估计器中，每个面对单元误差的贡献是以单元对角线为尺度的。然而，在修改后的版本中，我们采用了一个缩放系数，它取决于相邻元素的面的对角线和多项式度。两者之间的选择是通过类中定义的枚举器来完成的。
 *
 *
 * @note
 * 尽管名字是这样，Kelly估计器并不是真正的后验误差估计器，即使只应用于泊松问题。它为网格细化提供了很好的提示，但估计值是不值得信任的。对于高阶试验空间，这里计算的积分往往比误差本身更快归零，因此排除了作为误差估计的数值。然而，可以利用下面讨论的修改版本，通过增加残差（体积）部分来获得可靠的误差估计器。
 * 误差估计器实际上只对广义泊松方程 $-\nabla\cdot a(x) \nabla
 * u = f$
 * 的误差进行估计，该方程具有Dirichlet边界条件或涉及正交导数的广义诺伊曼边界条件
 * $a\frac{du}{dn} = g$  。
 * 误差估计器返回每个单元的估计误差向量，可用于输入
 * GridRefinement::refine_and_coarsen_fixed_fraction(),
 * GridRefinement::refine_and_coarsen_fixed_number(),
 * 和类似函数。这个向量包含数据类型为 @p float, 而不是 @p
 * double, 的元素，因为准确性在当前情况下并不重要。
 *
 *  <h3>Implementation</h3>
 * 原则上，误差估计的实现很简单：让@f[
 * \eta_K^2
 * =
 * \sum_{F\in\partial K}
 *   c_F \int_{\partial K_F} \jump{a \frac{\partial u_h}{\partial n}}^2
 * @f]成为单元 $K$ 的误差估计器。  $\jump{\cdot}$ 表示面的方括号中的函数的跳转， $c_F$ 是下面讨论的一个因子。这是Kelly等人在上面提到的论文中得出的误差估计器的界面项的一般形式。然后，整体误差估计被计算为@f[
 * \eta^2 = \sum_K \eta_K^2 @f]，这样 $\eta \approx \|\nabla (u-u_h)\|$
 * 为拉普拉斯方程。这一类的函数计算出一个对应于 $\eta_K$
 * （即上述数量的平方根）的值向量。 在Ainsworth的论文中 $
 * c_F=\frac {h_K}{24} $
 * ，但这个因素有点深奥，源于插值估计和稳定常数，这对泊松问题可能成立，但对更普遍的情况可能不成立。另外，我们考虑当
 * $c_F=\frac {h_F}{2p_F}$ ，其中 $h_F$ 是面的直径，
 * $p_F=max(p^+,p^-)$ 是相邻元素的最大多项式程度；或 $c_F=h_K$
 * 。这些因素之间的选择是通过枚举器完成的，在所有函数中作为最后一个参数提供。
 * 为了进行积分，我们使用了FEFaceValues和FESubfaceValues类。整合是通过在所有单元上循环进行的，并在尚未处理的面上进行整合。这样我们就避免了对面的整合两次，每次访问相邻的一个单元时都要进行一次。在所有单元的第二个循环中，我们将每个单元的面的贡献（是跳跃的积分平方乘以某个因子）相加，并取其平方根。
 * 积分是使用本类声明的 estimate()
 * 函数的调用者提供的面的正交公式进行的。对于线性试验函数（FE_Q(1)），有两个点的QGauss甚至QMidpoint规则实际上可能就足够了。对于高阶元素，有必要利用带有`fe.degree+1`高斯点的高阶正交公式。
 * 我们将每个面的贡献存储在C++标准库提供的 @p map,
 * 中，指向该面的迭代器是进入地图的关键。当对所有单元格进行第二次循环时，我们必须将各面的贡献相加并取其平方根。对于凯利估计器，与
 * $\frac {h_K}{24}$
 * 的乘法是在第二次循环中完成的。通过这样做，我们避免了决定与哪个
 * $h_K$
 * 相乘的问题，即与面的一侧的单元相乘还是与面的另一侧的单元相乘。而对于hp-estimator来说，
 * @p map 存储的是乘以 $\frac {h_F}{2p_F}$
 * 的积分，然后在第二个循环中进行求和。 $h_K$  (  $h_F$  )
 * 被认为是单元（面）对角线的最大长度。对于没有变形角的或多或少的均匀单元（面），这与单元（面）的直径相吻合。
 *
 *  <h3>Vector-valued functions</h3>
 * 如果要估计误差的有限元场是矢量值的，即有限元有一个以上的分量，那么上述所有内容可以同时应用于所有或只有一些分量。该类的主函数接收一个标志列表，表示误差估计器将应用于哪些分量；默认情况下，它是一个只有
 * @p trues, 的列表，意味着所有变量都将被处理。
 * 如果一个场的不同分量有不同的物理意义（例如斯托克斯方程中的速度和压力），对所有分量使用相同的系数是毫无意义的。在这种情况下，你可以通过一个函数，其分量与有限元场中的分量一样多，然后误差估计器的每个分量将由这个系数函数中的各自分量来加权。在另一种情况下，当所有分量具有相同的意义时（例如Lam&eacute;的弹性方程中的位移），你可以指定一个标量系数，然后将用于所有分量。
 *
 *  <h3>Boundary values</h3>
 * 如果面处于边界，即没有可以计算梯度跳跃的相邻单元，则有两种可能性。  <ul>   <li>  该面属于一个Dirichlet边界。那么这个面就不被考虑，这可以从对偶问题技术的角度来证明，如果边界能被所用的有限元准确地逼近，那么这个面就应该是准确的（即对于线性有限元是线性边界，对于等参量二次元是二次边界，等等）。对于不能精确近似的边界，应该考虑面的差值 $z-z_h$ ， $z$ 是对偶问题的解，在真实边界为零， $z_h$ 是近似值，在大多数情况下，在数值边界为零。由于在数值边界上 $z$ 一般不会为零，我们会在这里得到另一个项，但由于实际原因，这个项被忽略了，希望这里的误差会比我们想要估计的能量误差更快地趋于零。
 * 尽管没有必要进行积分，但在面的贡献列表中，我们为这个面存储了一个零，这使得将不同面对单元的贡献相加更容易。
 * <li>  该面属于一个诺伊曼边界。 在这种情况下，面
 * $F\in\partial K$ 的贡献看起来像\f[ n_F\int_F \left|g-a\frac{\partial
 * u_h}{\partial n}\right|^2 ds \f]，其中 $g$ 是Neumann边界函数，
 * $n_F=\frac {h_K}{24}$ 和 $n_F=\frac {h_F}{p}$
 * 分别是Kelly和hp-estimator的。如果有限元是矢量值的，那么显然表示诺伊曼边界条件的函数也需要是矢量值的。
 * <li>  没有考虑其他边界条件。  </ul>
 * 在实践中，如果你有罗宾边界条件或者懒得准确描述诺伊曼值，那么这很少是一个问题：如果你在地图中没有说任何关于边界的特定部分，那么凯利指标将简单地假定解决方案在边界的那一部分是正确的，而不去碰它。当然，如果你有一个有诺伊曼或罗宾的边界，这并不完全正确，数字解的法向导数和诺伊曼值之间会有差异，这些法向导数应该等于。因此，如果我们简单地忽略边界的这些部分，我们会低估误差。在实践中，这很少出现问题
 *
 * -你可能这次没有细化单元，但你可能会在下一个细化步骤中细化它，一切又都好了。毕竟，除了拉普拉斯方程之外，对于所有的问题，凯利指标只是一个指标，而不是一个估计器，因此它所计算的值无论如何都不是精确的误差表示。
 *
 *  <h3>Handling of hanging nodes</h3>
 * 沿着有悬挂节点的面进行积分是相当棘手的，因为其中一个元素必须向上或向下移动一级。关于这个话题的更多技术问题，请参见FESubFaceValues类的文档。
 * 在praxi中，由于我们只对每个面进行一次积分，所以当我们在与子面相邻的两个单元中较粗的一个单元中进行积分（子面被定义为一个面的子；从粗的单元看，它是一个子面，而从精的单元看，它是它的一个面）。原因是这样的话，寻找相邻关系的信息就比较容易了，但这都是实际的推理，没有什么根本性的。
 * 由于我们从面的粗面进行积分，所以我们可以随时得到母面，并将母面的积分结果（沿子面的积分之和）也存储在上述的积分图中。这样做会消耗一些超过需要的内存，但会使单元格的面贡献的求和变得更容易，因为那时我们手头有来自所有单元格的所有面的信息，不需要考虑明确地确定一个面是否被提炼。这同样适用于边界面，见上文。
 *
 *  <h3>Multiple solution vectors</h3>
 * 在某些情况下，例如在时间依赖性问题中，人们希望一次计算同一网格上的几个解向量的误差估计，具有相同的系数、边界条件对象等，例如在几个连续的时间步骤上的解。然后，人们可以为每个解多次调用这一类的函数。然而，计算误差估计值的最大因素（就计算时间而言）是FEFaceValues和FESubFaceValues对象的初始化，以及对所有面和子面进行迭代。如果求解向量生活在同一个网格上，通过同时处理所有的求解向量，每个单元只初始化一次FEFaceValues对象，同时对所有的求解向量只循环一次，可以大大减少这种工作量。出于这个原因，除了这个类中的
 * @p estimate
 * 函数接收一个输入向量并返回一个输出向量外，还有一个函数可以同时接受几个输入和输出向量。
 *
 *
 * @ingroup numerics
 *
 *
 */
template <int dim, int spacedim = dim>
class KellyErrorEstimator
{
public:
  /**
   * 给予类函数的枚举类型，用于决定面积分的缩放系数。
   *
   */
  enum Strategy
  {
    //! Kelly error estimator with the factor $\frac {h_K}{24}$.
    cell_diameter_over_24 = 0,
    //! the boundary residual estimator with the factor $\frac {h_F}{2
    //! max(p^+,p^-)}$.
    face_diameter_over_twice_max_degree,
    //! Kelly error estimator with the factor $h_K$.
    cell_diameter
  };

  /**
   * 上面描述的误差估计器的实现。你可以给出一个系数，但有一个默认值，表示数值为1的常数系数。系数函数可以是一个标量函数，在这种情况下，它被用于有限元的所有分量，或者是一个矢量值的函数，其分量与有限元中的分量一样多；在后一种情况下，每个分量由系数中的各自分量加权计算。
   * 如果DoFHandler对象使用的有限元是矢量值的，你可以给出一个你想评估的分量列表。然后你必须将位向量
   * @p component_mask 中的那些条目设置为真（见 @ref
   * GlossComponentMask
   * ），在误差估计器中使用相应的分量。默认情况下是使用所有的组件，这可以通过提供一个包含所有设置项的位向量或者一个空的位向量来实现。
   * @p subdomain_id
   * 参数表示我们是否要计算所有单元的指标（如果它的值是默认的，
   * <tt>numbers::invalid_unsigned_int</tt>),
   * ），或者只计算属于某个子域的单元的给定指标。后一种情况用于并行计算，即所有处理器节点都存储有全局网格，并且可以自己计算所有单元的所有指标，但是让每个进程只计算它所拥有的单元的指标，并让它们在事后交换结果信息，这样做更有效率。这尤其适用于这样的情况：网格非常大，计算
   * @em
   * 每个单元的指标过于昂贵，而只计算局部单元的指标是可以接受的。请注意，如果你只要求计算某个子域的指标，你必须确保这个函数能够访问
   * @em
   * 所有自由度的正确节点值。这是因为该函数也需要访问相邻单元的自由度值，即使它们属于不同的子域。
   * @p material_id
   * 参数也有类似的含义：如果不设置为其默认值（即
   * numbers::invalid_material_id),
   * ，则意味着只对具有该特定材料ID的单元计算指标。
   * @p n_threads
   * 参数用来指示用于计算误差估计器的线程数。现在这个参数被忽略了，线程数会自动确定。保留该参数是为了与旧版本的库兼容。
   * @p strategy 参数用于选择单元格面上的积分的比例因子。
   * @note  如果作为该函数参数的DoFHandler对象建立在
   * parallel::distributed::Triangulation,
   * 的基础上，该函数将跳过所有不属于本地的单元的计算。在这种情况下，subdomain_id参数的唯一有效值（除了无效值）是与当前处理器相关的子域id，如
   * parallel::distributed::Triangulation::locally_owned_subdomain().
   * 所报告的那样。
   * 即使没有对我们不属于本地的单元进行计算，错误指示向量的长度仍然必须等于
   * parallel::distributed::Triangulation::n_locally_owned_active_cells().
   * 所报告的网格中活动单元的数量。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficients   = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);

  /**
   * 调用 @p estimate 函数，见上文，<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficients   = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);

  /**
   * 与上述函数相同，但接受一个以上的解向量，并为每个解向量返回一个错误向量。关于这个函数存在的原因，请看这个类的一般文档。
   * 由于我们不想强迫这个函数的用户复制他们的解向量，解向量的向量采取指向解的指针，而不是向量的向量。这使得解向量在内存中的某个地方变得更简单，而不是在某个特殊的地方收集它们。(注意，不可能构建引用的向量，所以我们不得不使用指针的向量)。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);

  /**
   * 调用 @p estimate 函数，见上文，<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数需要一个正交集，用于hp-finite
   * element dof处理程序。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim - 1> & quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficients   = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数为hp-finite element
   * dof处理程序取一个正交集合。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim - 1> & quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficients   = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数为hp-finite element
   * dof处理程序取一个正交集合。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim - 1> & quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数为hp-finite element
   * dof处理程序取一个正交集合。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim - 1> & quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcInvalidComponentMask,
                   "You provided a ComponentMask argument that is invalid. "
                   "Component masks need to be either default constructed "
                   "(in which case they indicate that every component is "
                   "selected) or need to have a length equal to the number "
                   "of vector components of the finite element in use "
                   "by the DoFHandler object. In the latter case, at "
                   "least one component needs to be selected.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcInvalidCoefficient,
    "If you do specify the argument for a (possibly "
    "spatially variable) coefficient function for this function, "
    "then it needs to refer to a coefficient that is either "
    "scalar (has one vector component) or has as many vector "
    "components as there are in the finite element used by "
    "the DoFHandler argument.");
  /**
   * 异常情况
   *
   */
  DeclException3(ExcInvalidBoundaryFunction,
                 types::boundary_id,
                 int,
                 int,
                 << "You provided a function map that for boundary indicator "
                 << arg1 << " specifies a function with " << arg2
                 << " vector components. However, the finite "
                    "element in use has "
                 << arg3
                 << " components, and these two numbers need to match.");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcIncompatibleNumberOfElements,
                 int,
                 int,
                 << "The number of input vectors, " << arg1
                 << " needs to be equal to the number of output vectors, "
                 << arg2
                 << ". This is not the case in your call of this function.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcNoSolutions,
                   "You need to specify at least one solution vector as "
                   "input.");
};



/**
 * 这是对1d的一般模板的特殊化。对于1d来说，实现方式的不同足以证明这种特殊化的合理性。1d和其他空间维度的基本区别在于，在1d中，没有单元格的面，只有线段之间的顶点，所以我们必须以不同的方式计算跳跃项。然而，这个类提供了与一般模板完全相同的公共函数，因此，用户不会看到任何区别。
 *
 *
 */
template <int spacedim>
class KellyErrorEstimator<1, spacedim>
{
public:
  /**
   * 给予类函数的枚举类型，用于决定面积分的缩放系数。
   *
   */
  enum Strategy
  {
    //! Kelly error estimator with the factor $\frac {h_K}{24}$.
    cell_diameter_over_24 = 0,
    //! the boundary residual estimator with the factor $\frac {h_F}{2
    //! max(p^+,p^-)}$.
    face_diameter_over_twice_max_degree,
    //! Kelly error estimator with the factor $h_K$.
    cell_diameter
  };

  /**
   * 上面描述的误差估计器的实现。你可以给出一个系数，但有一个默认值，表示数值为1的常数系数。系数函数可以是一个标量函数，在这种情况下，它被用于有限元的所有分量，或者是一个矢量值的函数，其分量与有限元中的分量一样多；在后一种情况下，每个分量都由系数中的各自分量加权。
   * 如果DoFHandler对象使用的有限元是矢量值的，你可以给出一个你想评估的分量列表。然后，你必须将位向量
   * @p component_mask
   * 中的那些条目设置为真，在误差估计器中使用各自的分量。默认情况下是使用所有的分量，这可以通过提供一个具有所有设置项的位向量或者一个空的位向量来实现。所有其他参数与用于2d及以上版本的一般情况相同。
   * 参数n_threads不再使用，将被忽略。
   * 多线程目前没有在1d中实现，但我们保留了相应的参数，以便与一般情况下的函数签名兼容。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<1, spacedim> &   mapping,
    const DoFHandler<1, spacedim> &dof,
    const Quadrature<0> &          quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficient    = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);

  /**
   * 调用 @p estimate
   * 函数，见上文，使用<tt>mapping=MappingQGeneric1<1>()</tt>。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<1, spacedim> &dof,
    const Quadrature<0> &          quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficients   = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);

  /**
   * 与上述函数相同，但接受一个以上的解向量，并为每个解向量返回一个误差向量。关于这个函数存在的原因，请看这个类的一般文档。
   * 由于我们不想强迫这个函数的用户复制他们的解向量，解向量的向量采取指向解的指针，而不是向量的向量。这使得解向量在内存中的某个地方变得更简单，而不是在某个特殊的地方收集它们。(注意，不可能构建引用的向量，所以我们不得不使用指针的向量)。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<1, spacedim> &   mapping,
    const DoFHandler<1, spacedim> &dof,
    const Quadrature<0> &          quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);

  /**
   * 调用 @p estimate
   * 函数，见上文，使用<tt>mapping=MappingQGeneric1<1>()</tt>。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<1, spacedim> &dof,
    const Quadrature<0> &          quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数需要一个正交集，用于hp-finite
   * element dof处理程序。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<1, spacedim> &   mapping,
    const DoFHandler<1, spacedim> &dof,
    const hp::QCollection<0> &     quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficients   = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数为hp-finite element
   * dof处理程序取一个正交集合。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<1, spacedim> &dof,
    const hp::QCollection<0> &     quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                       neumann_bc,
    const InputVector &       solution,
    Vector<float> &           error,
    const ComponentMask &     component_mask = ComponentMask(),
    const Function<spacedim> *coefficients   = nullptr,
    const unsigned int        n_threads      = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id   = numbers::invalid_subdomain_id,
    const types::material_id  material_id    = numbers::invalid_material_id,
    const Strategy            strategy       = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数为hp-finite element
   * dof处理程序取一个正交集合。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const Mapping<1, spacedim> &   mapping,
    const DoFHandler<1, spacedim> &dof,
    const hp::QCollection<0> &     quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);


  /**
   * 相当于上面的函数集，只是这个函数为hp-finite element
   * dof处理程序取一个正交集合。
   *
   */
  template <typename InputVector>
  static void
  estimate(
    const DoFHandler<1, spacedim> &dof,
    const hp::QCollection<0> &     quadrature,
    const std::map<types::boundary_id,
                   const Function<spacedim, typename InputVector::value_type> *>
      &                                     neumann_bc,
    const std::vector<const InputVector *> &solutions,
    std::vector<Vector<float> *> &          errors,
    const ComponentMask &                   component_mask = ComponentMask(),
    const Function<spacedim> *              coefficients   = nullptr,
    const unsigned int        n_threads    = numbers::invalid_unsigned_int,
    const types::subdomain_id subdomain_id = numbers::invalid_subdomain_id,
    const types::material_id  material_id  = numbers::invalid_material_id,
    const Strategy            strategy     = cell_diameter_over_24);

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcInvalidComponentMask,
                   "You provided a ComponentMask argument that is invalid. "
                   "Component masks need to be either default constructed "
                   "(in which case they indicate that every component is "
                   "selected) or need to have a length equal to the number "
                   "of vector components of the finite element in use "
                   "by the DoFHandler object. In the latter case, at "
                   "least one component needs to be selected.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcInvalidCoefficient,
    "If you do specify the argument for a (possibly "
    "spatially variable) coefficient function for this function, "
    "then it needs to refer to a coefficient that is either "
    "scalar (has one vector component) or has as many vector "
    "components as there are in the finite element used by "
    "the DoFHandler argument.");
  /**
   * 异常情况
   *
   */
  DeclException3(ExcInvalidBoundaryFunction,
                 types::boundary_id,
                 int,
                 int,
                 << "You provided a function map that for boundary indicator "
                 << arg1 << " specifies a function with " << arg2
                 << " vector components. However, the finite "
                    "element in use has "
                 << arg3
                 << " components, and these two numbers need to match.");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcIncompatibleNumberOfElements,
                 int,
                 int,
                 << "The number of input vectors, " << arg1
                 << " needs to be equal to the number of output vectors, "
                 << arg2
                 << ". This is not the case in your call of this function.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcNoSolutions,
                   "You need to specify at least one solution vector as "
                   "input.");
};



DEAL_II_NAMESPACE_CLOSE

#endif


