

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Jean-Paul Pelteret, 2021
 */



// 我们首先包括所有必要的deal.II头文件和一些C++相关的文件。这第一个头文件将使我们能够访问一个数据结构，使我们能够在其中存储任意的数据。

#include <deal.II/algorithms/general_data_storage.h>

// 接下来是一些核心类，包括一个提供时间步进的实现。

#include <deal.II/base/discrete_time.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// 然后是一些标题，定义了一些有用的坐标变换和运动学关系，这些关系在非线性弹性中经常出现。

#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/transformations.h>

// 下面两个标头提供了我们进行自动微分所需的所有功能，并使用deal.II可以利用的符号计算机代数系统。所有自动微分和符号微分封装类的头文件，以及任何需要的辅助数据结构，都被收集在这些统一的头文件中。

#include <deal.II/differentiation/ad.h>
#include <deal.II/differentiation/sd.h>

// 包括这个头文件使我们有能力将输出写入文件流中。

#include <fstream>

// 按照惯例，整个教程程序被定义在它自己独特的命名空间中。

namespace Step71
{
  using namespace dealii;
  // @sect3{An introductory example: The fundamentals of automatic and symbolic
  // differentiation}

  // 自动和象征性的区分有一些神奇和神秘的特质。尽管在一个项目中使用它们会因多种原因而受益，但了解如何使用这些框架或如何利用它们的障碍可能会超过试图将它们（可靠地）整合到工作中的开发者的耐心。

  // 尽管作者希望能够成功地说明这些工具是如何被整合到有限元建模的工作流程中的，但最好还是先退一步，从基础开始。因此，一开始，我们先看看如何使用这两个框架来区分一个
  // "简单
  // "的数学函数，这样就可以牢固地建立和理解基本的操作（包括它们的顺序和功能），并使其复杂程度降到最低。在本教程的第二部分，我们将把这些基本原理付诸实践，并在此基础上进一步发展。

  // 伴随着对使用框架的算法步骤的描述，我们将对它们*可能在后台做的事情有一个简化的看法。这种描述在很大程度上是为了帮助理解，我们鼓励读者查看 @ref auto_symb_diff 模块文档，以获得对这些工具实际工作的更正式描述。

  //  @sect4{An analytical function}
  namespace SimpleExample
  {
    // 为了让读者相信这些工具在实践中确实有用，让我们选择一个函数，用手计算分析导数并不难。只是它的复杂程度足以让你考虑是否真的要去做这个练习，也可能让你怀疑你是否完全确定你对其导数的计算和实现是正确的。当然，问题的关键在于，函数的微分在某种意义上是相对公式化的，应该是计算机擅长的事情--如果我们能在现有的软件基础上理解这些规则，我们就不必费力地自己做了。

    // 我们为此选择了双变量三角函数 $f(x,y) = \cos\left(\frac{y}{x}\right)$
    // 。注意，这个函数是以数字类型为模板的。这样做是因为我们经常（但不总是）可以使用特殊的自动微分和符号类型来替代实数或复数类型，然后这些类型将执行一些基本的计算，例如评估一个函数值及其导数。我们将利用这一特性，确保我们只需要定义一次我们的函数，然后就可以在我们希望对其进行微分操作的任何情况下重新使用。

    template <typename NumberType>
    NumberType
    f(const NumberType &x, const NumberType &y)
    {
      return std::cos(y / x);
    }

    // 我们没有立即揭示这个函数的导数，而是向前声明返回导数的函数，并将它们的定义推迟到以后。正如函数名称所暗示的，它们分别返回导数
    // $\frac{df(x,y)}{dx}$  。

    double
    df_dx(const double x, const double y);
    // $\frac{df(x,y)}{dy}$  :

    double
    df_dy(const double x, const double y);
    // $\frac{d^{2}f(x,y)}{dx^{2}}$  :

    double
    d2f_dx_dx(const double x, const double y);
    // $\frac{d^{2}f(x,y)}{dx dy}$  :

    double
    d2f_dx_dy(const double x, const double y);
    // $\frac{d^{2}f(x,y)}{dy dx}$  :

    double
    d2f_dy_dx(const double x, const double y);

    // 最后是  $\frac{d^{2}f(x,y)}{dy^{2}}$  。

    double
    d2f_dy_dy(const double x, const double y);
    // @sect4{Computing derivatives using automatic differentiation}

    // 首先，我们将使用AD作为工具，为我们自动计算导数。我们将用参数`x`和`y`来评估函数，并期望得到的值和所有的导数都能在给定的公差范围内匹配。

    void
    run_and_verify_ad(const double x, const double y, const double tol = 1e-12)
    {
      // 我们的函数 $f(x,y)$
      // 是一个标量值函数，其参数代表代数计算或张量计算中遇到的典型输入变量。由于这个原因，
      // Differentiation::AD::ScalarFunction
      // 类是合适的包装类，可以用来做我们需要的计算。(作为比较，如果函数参数代表有限元单元的自由度，我们会希望以不同的方式处理它们)。问题的空间维度是不相关的，因为我们没有矢量或张量值的参数需要容纳，所以`dim`模板参数被任意分配为1的值。
      // 第二个模板参数规定了将使用哪个AD框架（deal.II支持几个外部AD框架），以及这个框架提供的基础数字类型将被使用。这个数字类型影响了微分运算的最大顺序，以及用于计算它们的基础算法。鉴于其模板性质，这个选择是一个编译时的决定，因为许多（但不是全部）AD库利用编译时的元编程，以有效的方式实现这些特殊的数字类型。第三个模板参数说明了结果类型是什么；在我们的例子中，我们要处理的是
      // "双数"。

      constexpr unsigned int                     dim = 1;
      constexpr Differentiation::AD::NumberTypes ADTypeCode =
        Differentiation::AD::NumberTypes::sacado_dfad_dfad;
      using ADHelper =
        Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>;

      // 我们有必要在我们的 @p ADHelper 类中预先登记函数 $f(x,y)$ 有多少个参数（我们将称之为 "独立变量"）。这些参数是`x`和`y`，所以显然有两个。

      constexpr unsigned int n_independent_variables = 2;

      // 我们现在有足够的信息来创建和初始化一个辅助类的实例。我们还可以得到具体的数字类型，它将在所有后续计算中使用。这很有用，因为我们可以从这里开始通过引用这个类型来编写一切，如果我们想改变使用的框架或数字类型（例如，如果我们需要更多的微分运算），那么我们只需要调整`ADTypeCode`模板参数。

      ADHelper ad_helper(n_independent_variables);
      using ADNumberType = typename ADHelper::ad_type;

      // 下一步是在辅助类中注册自变量的数值。这样做是因为函数和它的导数将正好针对这些参数进行评估。由于我们按照`{x,y}`的顺序注册它们，变量`x`将被分配到分量号`0`，而`y`将是分量号`1`--这个细节将在接下来的几行中使用。

      ad_helper.register_independent_variables({x, y});

      // 我们现在要求辅助类向我们提供自变量及其自动区分的表示。这些被称为
      // "敏感变量"，因为从现在开始，我们对组件`独立变量_ad`所做的任何操作都会被AD框架跟踪和记录，并且在我们要求计算它们的导数时，会被考虑。帮助器返回的是一个可自动微分的
      // "向量"，但是我们可以确定，第2个元素代表 "x"，第1个元素代表
      // "y"。为了完全确保这些变量的数字类型没有任何歧义，我们给所有的自动微分变量加上`ad'的后缀。

      const std::vector<ADNumberType> independent_variables_ad =
        ad_helper.get_sensitive_variables();
      const ADNumberType &x_ad = independent_variables_ad[0];
      const ADNumberType &y_ad = independent_variables_ad[1];

      // 我们可以立即将自变量的敏感表示法传递给我们的模板函数，计算出  $f(x,y)$
      // 。这也会返回一个可自动微分的数字。

      const ADNumberType f_ad = f(x_ad, y_ad);

      // 所以现在要问的自然是，我们把这些特殊的`x_ad`和`y_ad`变量传递给函数`f`，而不是原来的`double`变量`x`和`y`，实际上计算了什么？换句话说，这一切与我们想要确定的导数的计算有什么关系？或者，更简洁地说。这个返回的`ADNumberType`对象有什么特别之处，使它有能力神奇地返回导数？

      // 从本质上讲，这*可以*做的是以下几点。这个特殊的数字可以被看作是一个数据结构，它存储了函数值，以及规定的导数数量。对于一个期望有两个参数的一次可导数，它可能看起来像这样。

      //
      // @code
      //  struct ADNumberType
      //  {
      //    double value; The value of the object
      //    double derivatives[2]; Array of derivatives of the object with
      //                  respect to x and y
      //  };
      //  @endcode

      // 对于我们的自变量`x_ad`，`x_ad.value`的起始值只是它的赋值（即这个变量代表的实值）。导数`x_ad.derivatives[0]`将被初始化为`1'，因为`x'是第2个独立变量和
      // $\frac{d(x)}{dx} = 1$
      // 。导数`x.derivatives[1]`将被初始化为零，因为第一个自变量是`y`和
      // $\frac{d(x)}{dy} = 0$  。

      // 为了使函数导数有意义，我们必须假设这个函数不仅在分析意义上是可微的，而且在评估点`x,y`也是可微的。我们可以利用这两个假设：当我们在数学运算中使用这种数字类型时，AD框架可以**的
      //重载操作（例如，`%operator+()`, `%operator*()`以及`%sin()`, `%exp()`,
      //等等），使返回的结果具有预期值。同时，它将通过对被重载的确切函数的了解和对连锁规则的严格应用来计算导数。因此，`%sin()`函数（其参数`a`本身是自变量`x`和`y`的一个函数
      //*可能*被定义如下。

      //
      // @code
      //  ADNumberType sin(const ADNumberType &a)
      //  {
      //    ADNumberType output;



      // //  For the input argument "a", "a.value" is simply its value.
      // //    output.value = sin(a.value);

      // //  We know that the derivative of sin(a) is cos(a), but we need
      // //  to also consider the chain rule and that the input argument
      // //  `a` is also differentiable with respect to the original
      // //  independent variables `x` and `y`. So `a.derivatives[0]`
      // //  and `a.derivatives[1]` respectively represent the partial
      // //  derivatives of `a` with respect to its inputs `x` and `y`.
      //    output.derivatives[0] = cos(a.value)*a.derivatives[0];
      //    output.derivatives[1] = cos(a.value)*a.derivatives[1];



      //    return output;
      //  }
      //  @endcode

      // 当然，所有这些也可以用于二阶甚至高阶导数。

      // 所以现在很明显，通过上述表示，`ADNumberType`携带了一些额外的数据，这些数据代表了可微调函数相对于原始（敏感）自变量的各种导数。因此应该注意到，使用它们会产生计算开销（因为我们在做导数计算时要计算额外的函数），以及存储这些结果的内存开销。因此，规定的微分运算的级数最好保持在最低水平，以限制计算成本。例如，我们可以自己计算第一级导数，然后使用
      // Differentiation::AD::VectorFunction
      // 辅助类来确定依赖函数集合的梯度，这将是原始标量函数的第二级导数。

      // 还值得注意的是，由于链式规则是无差别应用的，我们只看到计算的起点和终点`{x,y}`
      // $\rightarrow$
      // `f(x,y)`，我们永远只能查询到`f`的总导数；部分导数（上例中的`a.导数[0]`和`a.导数[1]`）是中间值，对我们是隐藏的。

      // 好的，既然我们现在至少知道了`f_ad'代表什么，以及它的编码是什么，让我们把所有的东西用于实际的用途。为了获得那些隐藏的派生结果，我们将最终结果注册到帮助类中。在这之后，我们不能再改变`f_ad`的值，也不能让这些变化反映在帮助者类返回的结果中。

      ad_helper.register_dependent_variable(f_ad);

      // 下一步是提取导数（特别是函数梯度和Hessian）。为此，我们首先创建一些临时数据结构（结果类型为`double`）来存储导数（注意，所有的导数都是一次性返回的，而不是单独返回）...

      Vector<double>     Df(ad_helper.n_dependent_variables());
      FullMatrix<double> D2f(ad_helper.n_dependent_variables(),
                             ad_helper.n_independent_variables());

      // ...
      // 然后我们要求助手类计算这些导数，以及函数值本身。就这样了。我们得到了我们想得到的一切。

      const double computed_f = ad_helper.compute_value();
      ad_helper.compute_gradient(Df);
      ad_helper.compute_hessian(D2f);

      // 我们可以通过与分析方案的比较来说服自己，AD框架是正确的。(或者，如果你像作者一样，你会做相反的事情，宁愿验证你对分析方案的实现是正确的！)

      AssertThrow(std::abs(f(x, y) - computed_f) < tol,
                  ExcMessage(std::string("Incorrect value computed for f. ") +
                             std::string("Hand-calculated value: ") +
                             Utilities::to_string(f(x, y)) +
                             std::string(" ; ") +
                             std::string("Value computed by AD: ") +
                             Utilities::to_string(computed_f)));

      // 因为我们知道自变量的排序，所以我们知道梯度的哪个部分与哪个导数有关......。

      const double computed_df_dx = Df[0];
      const double computed_df_dy = Df[1];

      AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for df/dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(df_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_df_dx)));
      AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for df/dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(df_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_df_dy)));

      // .......对于Hessian也是如此。

      const double computed_d2f_dx_dx = D2f[0][0];
      const double computed_d2f_dx_dy = D2f[0][1];
      const double computed_d2f_dy_dx = D2f[1][0];
      const double computed_d2f_dy_dy = D2f[1][1];

      AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dx_dx)));
      AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dx_dy)));
      AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dy_dx)));
      AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_d2f_dy_dy)));
    }

    // 这很不错。在计算这个三角函数的二阶导数时并没有太多的工作。

    //  @sect4{Hand-calculated derivatives of the analytical solution}

    // 既然我们现在知道了让AD框架为我们计算这些导数需要多少
    // "执行工作"，让我们把它与手工计算并在几个独立的函数中实现的同样的导数进行比较。

    // 这里是 $f(x,y) = \cos\left(\frac{y}{x}\right)$ 的两个一阶导数。

    //  $\frac{df(x,y)}{dx} = \frac{y}{x^2} \sin\left(\frac{y}{x}\right)$
    double
    df_dx(const double x, const double y)
    {
      Assert(x != 0.0, ExcDivideByZero());
      return y * std::sin(y / x) / (x * x);
    }
    // $\frac{df(x,y)}{dx} = -\frac{1}{x} \sin\left(\frac{y}{x}\right)$
    double
    df_dy(const double x, const double y)
    {
      return -std::sin(y / x) / x;
    }

    // 这里是 $f(x,y)$ 的四个二次导数。

    //  $\frac{d^{2}f(x,y)}{dx^{2}} = -\frac{y}{x^4} (2x
    //  \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$
    double
    d2f_dx_dx(const double x, const double y)
    {
      return -y * (2 * x * std::sin(y / x) + y * std::cos(y / x)) /
             (x * x * x * x);
    }
    // $\frac{d^{2}f(x,y)}{dx dy} = \frac{1}{x^3} (x
    // \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$
    double
    d2f_dx_dy(const double x, const double y)
    {
      return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x);
    }
    // $\frac{d^{2}f(x,y)}{dy dx} = \frac{1}{x^3} (x
    // \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$
    // （正如预期的那样，根据[施瓦茨定理]（https:en.wikipedia.org/wiki/Symmetry_of_second_derivatives））。

    double
    d2f_dy_dx(const double x, const double y)
    {
      return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x);
    }
    // $\frac{d^{2}f(x,y)}{dy^{2}} = -\frac{1}{x^2}
    // \cos\left(\frac{y}{x}\right)$
    double
    d2f_dy_dy(const double x, const double y)
    {
      return -(std::cos(y / x)) / (x * x);
    }

    // 嗯......上面有很多地方我们可以引入错误，特别是在应用链式规则的时候。虽然它们不是银弹，但至少这些AD框架可以作为一个验证工具，确保我们没有犯任何错误（无论是计算还是执行），从而对我们的结果产生负面影响。

    // 当然，这个例子的重点是，我们可能选择了一个相对简单的函数 $f(x,y)$
    // ，我们可以手工验证AD框架计算的导数是否正确。但是AD框架并不关心这个函数是否简单。它可能是一个复杂得多的表达式，或者取决于两个以上的变量，它仍然能够计算出导数--唯一的区别是，*我们*不能再想出导数来验证AD框架的正确性。

    //  @sect4{Computing derivatives using symbolic differentiation}

    // 我们现在要用符号微分法重复同样的练习。术语 "符号微分
    // "有点误导，因为微分只是计算机代数系统（CAS）（即符号框架）提供的一个工具。然而，在有限元建模和应用的背景下，它是CAS最常见的用途，因此将是我们关注的重点。再一次，我们将提供参数值`x`和`y`来评估我们的函数
    // $f(x,y) = \cos\left(\frac{y}{x}\right)$
    // 和它的导数，并提供一个公差来测试返回结果的正确性。

    void
    run_and_verify_sd(const double x, const double y, const double tol = 1e-12)
    {
      // 我们需要做的第一步是形成符号变量，代表我们希望对其进行微分的函数参数。同样，这些将是我们问题的独立变量，因此在某种意义上是原始变量，与其他变量没有任何关系。我们通过初始化一个符号类型
      // Differentiation::SD::Expression,
      // 来创建这些类型的（独立）变量，这个符号类型是对符号框架所使用的一组类的包装，有一个唯一的标识。在这种情况下，这个标识符，一个
      // `std::string`, 对于 $x$ 的参数来说，是简单的 "x"，同样，对于依赖函数的
      // $y$ 参数来说，也是
      // "y"。像以前一样，我们将用`sd`作为符号变量名称的后缀，这样我们就可以清楚地看到哪些变量是符号性的（而不是数字性的）。

      const Differentiation::SD::Expression x_sd("x");
      const Differentiation::SD::Expression y_sd("y");

      // 使用计算 $f(x,y)$
      // 的模板化函数，我们可以将这些独立变量作为参数传递给该函数。返回的结果将是另一个符号类型，代表用于计算
      // $\cos\left(\frac{y}{x}\right)$  的操作序列。

      const Differentiation::SD::Expression f_sd = f(x_sd, y_sd);

      // 在这一点上，打印出表达式`f_sd`是合法的，如果我们这样做的话
      // @code
      //  std::cout << "f(x,y) = " << f_sd << std::endl;
      //  @endcode
      //我们会看到`f(x,y) = cos(y/x)`打印到控制台。

      // 你可能会注意到，我们在构建我们的符号函数`f_sd`时，没有说明我们可能要如何使用它。与上面显示的AD方法相比，我们从调用`f(x_sd,
      // y_sd)`返回的不是函数`f`在某个特定点的评价，而实际上是在一个通用的、尚未确定的点的评价的符号表示。这是使符号框架（CAS）不同于自动区分框架的关键点之一。每个变量`x_sd`和`y_sd`，甚至复合依赖函数`f_sd`，在某种意义上分别是数值的
      // "占位符
      // "和操作的组合。事实上，用于组成函数的各个组件也是占位符。操作序列被编码成一个树状的数据结构（概念上类似于[抽象语法树](https:en.wikipedia.org/wiki/Abstract_syntax_tree)）。

      // 一旦我们形成了这些数据结构，我们就可以把我们可能想对它们进行的任何操作推迟到以后的某个时间。这些占位符中的每一个都代表了一些东西，但我们有机会在任何方便的时间点上定义或重新定义它们所代表的东西。因此，对于这个特定的问题，我们想把
      // "x "和 "y
      // "与*一些*数值（类型尚未确定）联系起来是有道理的，但我们可以在概念上（如果有意义的话）给
      // "y/x "这个比率赋值，而不是单独给 "x "和 "y "这些变量赋值。我们还可以将
      // "x "或 "y
      // "与其他一些符号函数`g(a,b)`联系起来。这些操作中的任何一个都涉及到对所记录的操作树的操作，以及用其他东西替换树上的突出节点（以及该节点的子树）。这里的关键词是
      // "替换"，事实上，在 Differentiation::SD
      // 命名空间中有许多函数的名称中都有这个词。

      // 这种能力使框架完全通用。在有限元模拟的背景下，我们通常会对我们的符号类型进行的操作类型是函数组合、微分、替换（部分或完全）和评估（即符号类型向其数字对应物的转换）。但如果你需要，一个CAS的能力往往不止这些。它可以形成函数的反导数（积分），对形成函数的表达式进行简化（例如，用
      // $1$ 替换 $(\sin a)^2 + (\cos a)^2$
      // ；或者，更简单：如果函数做了像`1+2`这样的运算，CAS可以用`3`替换它），等等。变量所代表的*表达式是从函数
      // $f$ 的实现方式中得到的，但CAS可以对其进行任何功能的操作。

      // 具体来说，为了计算因果函数相对于各个自变量的一阶导数的符号表示，我们使用
      // Differentiation::SD::Expression::differentiate()
      // 函数，自变量作为其参数。每次调用都会导致CAS通过组成`f_sd`的运算树，并对表达式树的每个节点进行相对于给定符号参数的微分。

      const Differentiation::SD::Expression df_dx_sd = f_sd.differentiate(x_sd);
      const Differentiation::SD::Expression df_dy_sd = f_sd.differentiate(y_sd);

      // 为了计算二阶导数的符号表示，我们只需对自变量的一阶导数进行微分。所以要计算高阶导数，我们首先需要计算低阶导数。由于调用
      // "differentiate()
      // "的返回类型是一个表达式，我们原则上可以通过将两个调用连在一起，直接从标量上执行双倍微分。但是在这个特殊的情况下，这是不需要的，因为我们手头有中间结果）。)

      const Differentiation::SD::Expression d2f_dx_dx_sd =
        df_dx_sd.differentiate(x_sd);
      const Differentiation::SD::Expression d2f_dx_dy_sd =
        df_dx_sd.differentiate(y_sd);
      const Differentiation::SD::Expression d2f_dy_dx_sd =
        df_dy_sd.differentiate(x_sd);
      const Differentiation::SD::Expression d2f_dy_dy_sd =
        df_dy_sd.differentiate(y_sd);

      // 使用语句
      // @code
      //  std::cout << "df_dx_sd: " << df_dx_sd << std::endl;
      //  std::cout << "df_dy_sd: " << df_dy_sd << std::endl;
      //  std::cout << "d2f_dx_dx_sd: " << d2f_dx_dx_sd << std::endl;
      //  std::cout << "d2f_dx_dy_sd: " << d2f_dx_dy_sd << std::endl;
      //  std::cout << "d2f_dy_dx_sd: " << d2f_dy_dx_sd << std::endl;
      //  std::cout << "d2f_dy_dy_sd: " << d2f_dy_dy_sd << std::endl;
      //  @endcode
      //  打印由CAS计算的第一和第二导数的表达式，得到以下输出。
      //  @code{.sh}
      //  df_dx_sd: y*sin(y/x)/x**2
      //  df_dy_sd: -sin(y/x)/x
      //  d2f_dx_dx_sd: -y**2*cos(y/x)/x**4 - 2*y*sin(y/x)/x**3
      //  d2f_dx_dy_sd: sin(y/x)/x**2 + y*cos(y/x)/x**3
      //  d2f_dy_dx_sd: sin(y/x)/x**2 + y*cos(y/x)/x**3
      //  d2f_dy_dy_sd: -cos(y/x)/x**2
      //  @endcode
      //  这与前面介绍的这些导数的分析表达式相比，效果很好。

      // 现在我们已经形成了函数及其导数的符号表达式，我们想对函数的主要参数`x`和`y`的数字值进行评估。为了达到这个目的，我们构造了一个*替代图*，它将符号值映射到它们的数字对应值。

      const Differentiation::SD::types::substitution_map substitution_map =
        Differentiation::SD::make_substitution_map(
          std::pair<Differentiation::SD::Expression, double>{x_sd, x},
          std::pair<Differentiation::SD::Expression, double>{y_sd, y});

      // 这个过程的最后一步是将所有的符号变量和操作转换成数值，并产生这个操作的数值结果。为了做到这一点，我们在上面已经提到的步骤中，将替换图与符号变量结合起来。"替换"。

      // 一旦我们把这个替换图传递给CAS，它就会把符号变量的每个实例（或者更一般的，子表达式）替换成它的数字对应物，然后把这些结果在操作树上传播，如果可能的话，简化树上的每个节点。如果运算树被简化为一个单一的值（也就是说，我们已经将所有的独立变量替换成了它们的数字对应值），那么评估就完成了。

      // 由于C++的强类型特性，我们需要指示CAS将其对结果的表示转换为内在的数据类型（本例中为`double'）。这就是
      // "评估
      // "步骤，通过模板类型我们定义了这个过程的返回类型。方便的是，如果我们确定我们已经进行了完整的替换，这两个步骤可以一次完成。

      const double computed_f =
        f_sd.substitute_and_evaluate<double>(substitution_map);

      AssertThrow(std::abs(f(x, y) - computed_f) < tol,
                  ExcMessage(std::string("Incorrect value computed for f. ") +
                             std::string("Hand-calculated value: ") +
                             Utilities::to_string(f(x, y)) +
                             std::string(" ; ") +
                             std::string("Value computed by AD: ") +
                             Utilities::to_string(computed_f)));

      // 我们可以对第一个导数做同样的处理......

      const double computed_df_dx =
        df_dx_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_df_dy =
        df_dy_sd.substitute_and_evaluate<double>(substitution_map);

      AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for df/dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(df_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by AD: ") +
                    Utilities::to_string(computed_df_dx)));
      AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol, 
 
 
 
 
 
                    Utilities::to_string(computed_df_dy)));

      // ...以及二阶导数。请注意，我们可以在这些操作中重复使用相同的替换图，因为我们希望针对相同的`x`和`y`值评估所有这些函数。修改置换图中的值，就可以得到相同的符号表达式的评估结果，同时给自变量分配不同的值。我们也可以很高兴地让每个变量在一次中代表一个实值，在下一次中代表一个复值。

      const double computed_d2f_dx_dx =
        d2f_dx_dx_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_d2f_dx_dy =
        d2f_dx_dy_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_d2f_dy_dx =
        d2f_dy_dx_sd.substitute_and_evaluate<double>(substitution_map);
      const double computed_d2f_dy_dy =
        d2f_dy_dy_sd.substitute_and_evaluate<double>(substitution_map);

      AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dx_dx)));
      AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dx_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dx_dy)));
      AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dx. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dy_dx)));
      AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol,
                  ExcMessage(
                    std::string("Incorrect value computed for d2f/dy_dy. ") +
                    std::string("Hand-calculated value: ") +
                    Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") +
                    std::string("Value computed by SD: ") +
                    Utilities::to_string(computed_d2f_dy_dy)));
    }
    // @sect4{The SimpleExample::run() function}

    // 用来驱动这些初始例子的函数是直接的。我们将任意选择一些值来评估该函数（尽管知道`x
    // = 0`是不允许的），然后将这些值传递给使用AD和SD框架的函数。

    void
    run()
    {
      const double x = 1.23;
      const double y = 0.91;

      std::cout << "Simple example using automatic differentiation..."
                << std::endl;
      run_and_verify_ad(x, y);
      std::cout << "... all calculations are correct!" << std::endl;

      std::cout << "Simple example using symbolic differentiation."
                << std::endl;
      run_and_verify_sd(x, y);
      std::cout << "... all calculations are correct!" << std::endl;
    }

  } // namespace SimpleExample
  // @sect3{A more complex example: Using automatic and symbolic differentiation
  // to compute derivatives at continuum points}

  // 现在我们已经介绍了自动分化和符号分化背后的原理，我们将通过制定两个耦合的磁力学构成法将其付诸实施：一个是与速率无关的，另一个则表现为与速率有关的行为。

  // 正如你在介绍中记得的那样，我们将考虑的材料构成法则要比上面的简单例子复杂得多。这不仅仅是因为我们将考虑的函数
  // $\psi_{0}$ 的形式，而且特别是因为 $\psi_{0}$
  // 不仅仅取决于两个标量变量，而是取决于一大堆*张量，每个张量都有几个组成部分。在某些情况下，这些是*对称*张量，对于这些张量来说，只有一个分量子集实际上是独立的，我们必须考虑计算
  // $\frac{\partial\psi_{0}}{\partial \mathbf{C}}$ 这样的导数的实际意义，其中
  // $\mathbf C$
  // 是一个对称张量。希望这一切将在下面变得清晰。我们也将清楚地看到，用手来做这件事，在最好的情况下，将是非常*繁琐，而在最坏的情况下，充满了难以发现的错误。

  namespace CoupledConstitutiveLaws
  {
    // @sect4{Constitutive parameters}

    // 我们先描述一下能量函数描述中出现的各种材料参数  $\psi_{0}$  。

    // ConstitutiveParameters类被用来保存这些数值。所有参数的值（包括构成参数和流变参数）都来自于
    // @cite Pelteret2018a
    // ，并给出了能够产生大致代表真实的、实验室制造的磁活性聚合物的构成响应的值，当然，这里使用的具体数值对本程序的目的没有影响。

    // 前四个构成参数分别代表

    // 弹性剪切模量 $\mu_{e}$  。

    // --磁饱和时的弹性剪切模量  $\mu_{e}^{\infty}$  。

    // - 弹性剪切模量的饱和磁场强度  $h_{e}^{\text{sat}}$  ，以及

    // 泊松比  $\nu$  。

    class ConstitutiveParameters : public ParameterAcceptor
    {
    public:
      ConstitutiveParameters();

      double mu_e       = 30.0e3;
      double mu_e_inf   = 250.0e3;
      double mu_e_h_sat = 212.2e3;
      double nu_e       = 0.49;

      // 接下来的四个，只与速率相关的材料有关，是以下的参数

      // - 粘弹性剪切模量  $\mu_{v}$  。

      // - 磁饱和时的粘弹性剪切模量  $\mu_{v}^{\infty}$  。

      // 粘弹性剪切模量的饱和磁场强度 $h_{v}^{\text{sat}}$  ，以及

      // --特征松弛时间  $\tau$  。

      double mu_v       = 20.0e3;
      double mu_v_inf   = 35.0e3;
      double mu_v_h_sat = 92.84e3;
      double tau_v      = 0.6;

      // 最后一个参数是相对磁导率  $\mu_{r}$  。

      double mu_r = 6.0;

      bool initialized = false;
    };

    // 参数是通过ParameterAcceptor框架初始化的，该框架在  step-60 中有详细讨论。

    ConstitutiveParameters::ConstitutiveParameters()
      : ParameterAcceptor("/Coupled Constitutive Laws/Constitutive Parameters/")
    {
      add_parameter("Elastic shear modulus", mu_e);
      add_parameter("Elastic shear modulus at magnetic saturation", mu_e_inf);
      add_parameter(
        "Saturation magnetic field strength for elastic shear modulus",
        mu_e_h_sat);
      add_parameter("Poisson ratio", nu_e);

      add_parameter("Viscoelastic shear modulus", mu_v);
      add_parameter("Viscoelastic shear modulus at magnetic saturation",
                    mu_v_inf);
      add_parameter(
        "Saturation magnetic field strength for viscoelastic shear modulus",
        mu_v_h_sat);
      add_parameter("Characteristic relaxation time", tau_v);

      add_parameter("Relative magnetic permeability", mu_r);

      parse_parameters_call_back.connect([&]() { initialized = true; });
    }
    // @sect4{Constitutive laws: Base class}

    // 由于我们将为同一类材料制定两种构成法，因此定义一个基类以确保它们有统一的接口是有意义的。

    // 类的声明从构造函数开始，它将接受一组构成参数，这些参数与材料定律本身一起决定了材料的响应。

    template <int dim>
    class Coupled_Magnetomechanical_Constitutive_Law_Base
    {
    public:
      Coupled_Magnetomechanical_Constitutive_Law_Base(
        const ConstitutiveParameters &constitutive_parameters);

      // 我们将在一个方法中计算和存储这些值，而不是随意计算和返回动力学变量或其线性化。然后这些缓存的结果将在请求时返回。我们将把为什么要这样做的精确解释推迟到以后的阶段。现在重要的是看到这个函数接受所有的场变量，即磁场矢量 $\boldsymbol{\mathbb{H}}$ 和右Cauchy-Green变形张量 $\mathbf{C}$ ，以及时间离散器。除了 @p constitutive_parameters, 之外，这些都是计算材料响应所需的基本量。

      virtual void
      update_internal_data(const SymmetricTensor<2, dim> &C,
                           const Tensor<1, dim> &         H,
                           const DiscreteTime &           time) = 0;

      // 接下来的几个函数提供了探测材料响应的接口，这些响应是由于施加的变形和磁荷载引起的。

      // 由于该类材料可以用自由能 $\psi_{0}$ 来表示，我们可以计算出......

      virtual double
      get_psi() const = 0;

      // ... 以及两个动力学量。

      // 磁感应矢量  $\boldsymbol{\mathbb{B}}$  ，和

      // --皮奥拉-基尔霍夫总应力张量 $\mathbf{S}^{\text{tot}}$  。
      virtual Tensor<1, dim>
      get_B() const = 0;

      virtual SymmetricTensor<2, dim>
      get_S() const = 0;

      // .......以及动力学量的线性化，它们是。

      // --磁静力学正切张量  $\mathbb{D}$  。

      // - 总的参考性磁弹性耦合张量 $\mathfrak{P}^{\text{tot}}$  ，以及

      // --总的参考弹性正切张量 $\mathcal{H}^{\text{tot}}$  。

      virtual SymmetricTensor<2, dim>
      get_DD() const = 0;

      virtual Tensor<3, dim>
      get_PP() const = 0;

      virtual SymmetricTensor<4, dim>
      get_HH() const = 0;

      // 我们还将定义一个方法，为这个类实例提供一个机制，在进入下一个时间段之前做任何额外的任务。同样，这样做的原因将在稍后变得清晰。

      virtual void
      update_end_of_timestep()
      {}

      // 在该类的`保护'部分，我们存储了一个对支配材料响应的构成参数实例的引用。为了方便起见，我们还定义了一些函数来返回各种构成参数（包括明确定义的，以及计算的）。

      与材料的弹性响应有关的参数依次是： //。

        // - 弹性剪切模量。

        // - 饱和磁场下的弹性剪切模量。

        // - 弹性剪切模量的饱和磁场强度。

        // - 泊松比。

        // 泊松比、Lam&eacute;参数，以及

        // 体积模量。

        protected : const ConstitutiveParameters &constitutive_parameters;

      double
      get_mu_e() const;

      double
      get_mu_e_inf() const;

      double
      get_mu_e_h_sat() const;

      double
      get_nu_e() const;

      double
      get_lambda_e() const;

      double
      get_kappa_e() const;

      // 与材料的弹性响应有关的参数依次是

      // - 粘弹性剪切模量。

      // -- 磁饱和时的粘弹性剪切模量。

      // - 粘弹性剪切模量的饱和磁场强度，以及

      粘弹性剪切模量的饱和磁场强度，和 //--特征松弛时间。

        double
        get_mu_v() const;

      double
      get_mu_v_inf() const;

      double
      get_mu_v_h_sat() const;

      double
      get_tau_v() const;

      // 与材料的磁响应有关的参数依次是：。

      // 相对磁导率，以及

      // - 磁导率常数 $\mu_{0}$
      // （其实不是一个材料常数，而是一个普遍的常数，为了简单起见，我们在这里分组）。

      // 我们还将实现一个函数，从时间离散性中返回时间步长。

      double
      get_mu_r() const;

      constexpr double
      get_mu_0() const;
      double
      get_delta_t(const DiscreteTime &time) const;
    };

    // 在下文中，让我们从实现刚才定义的类的几个相对琐碎的成员函数开始。

    template <int dim>
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::
      Coupled_Magnetomechanical_Constitutive_Law_Base(
        const ConstitutiveParameters &constitutive_parameters)
      : constitutive_parameters(constitutive_parameters)
    {
      Assert(get_kappa_e() > 0, ExcInternalError());
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e() const
    {
      return constitutive_parameters.mu_e;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_inf() const
    {
      return constitutive_parameters.mu_e_inf;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_h_sat() const
    {
      return constitutive_parameters.mu_e_h_sat;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_nu_e() const
    {
      return constitutive_parameters.nu_e;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_lambda_e() const
    {
      return 2.0 * get_mu_e() * get_nu_e() / (1.0 - 2.0 * get_nu_e());
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_kappa_e() const
    {
      return (2.0 * get_mu_e() * (1.0 + get_nu_e())) /
             (3.0 * (1.0 - 2.0 * get_nu_e()));
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v() const
    {
      return constitutive_parameters.mu_v;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_inf() const
    {
      return constitutive_parameters.mu_v_inf;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_h_sat() const
    {
      return constitutive_parameters.mu_v_h_sat;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_tau_v() const
    {
      return constitutive_parameters.tau_v;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_r() const
    {
      return constitutive_parameters.mu_r;
    }

    template <int dim>
    constexpr double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_0() const
    {
      return 4.0 * numbers::PI * 1e-7;
    }

    template <int dim>
    double
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_delta_t(
      const DiscreteTime &time) const
    {
      return time.get_previous_step_size();
    }
    // @sect4{Magnetoelastic constitutive law (using automatic differentiation)}

    // 我们将首先考虑一种非耗散性材料，即受磁超弹性构成法则支配的材料，在浸入磁场时表现出僵硬。正如介绍中所述，这种材料的储能密度函数可能由
    // @f[
    //    \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  = \frac{1}{2} \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //      \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //      \right]
    //  + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    //  - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //      \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //      \boldsymbol{\mathbb{H}} \right]
    //  @f]
    //  和
    //  @f[
    //   f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //  = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //      \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //      \boldsymbol{\mathbb{H}}}
    //        {\left(h_{e}^{\text{sat}}\right)^{2}} \right) .
    //  @f]
    //  给出。

    // 现在来看看实现这种行为的类。由于我们希望这个类能完全描述一种材料，所以我们将它标记为
    // "final"，这样继承树就在这里终止了。在类的顶部，我们定义了辅助类型，我们将在标量能量密度函数的AD计算中使用它。请注意，我们希望它能返回
    // "double
    // "类型的值。我们还必须指定空间维度的数量，`dim'，以便建立矢量、张量和对称张量场与它们所含分量数量之间的联系。用于ADHelper类的具体的`ADTypeCode`将在实际使用该类的时候作为模板参数提供。

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    class Magnetoelastic_Constitutive_Law_AD final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
      using ADHelper =
        Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>;
      using ADNumberType = typename ADHelper::ad_type;

    public:
      Magnetoelastic_Constitutive_Law_AD(
        const ConstitutiveParameters &constitutive_parameters);

      // 由于基类的公共接口是纯 "虚拟
      // "的，这里我们将声明这个类将覆盖所有这些基类方法。

      virtual void
      update_internal_data(const SymmetricTensor<2, dim> &C,
                           const Tensor<1, dim> &         H,
                           const DiscreteTime &) override;

      virtual double
      get_psi() const override;

      virtual Tensor<1, dim>
      get_B() const override;

      virtual SymmetricTensor<2, dim>
      get_S() const override;

      virtual SymmetricTensor<2, dim>
      get_DD() const override;

      virtual Tensor<3, dim>
      get_PP() const override;

      virtual SymmetricTensor<4, dim>
      get_HH() const override;

      // 在这个类的`private`部分，我们需要定义一些提取器，这些提取器将帮助我们设置自变量，随后得到与因变量相关的计算值。如果这个类是在有限元问题的背景下使用，那么这些提取器中的每一个都（很可能）与解场的一个分量（在本例中，位移和磁标势）的梯度有关。正如你现在可能推断的那样，这里
      // "C "表示右Cauchy-Green张量，"H "表示磁场向量。

    private:
      const FEValuesExtractors::Vector             H_components;
      const FEValuesExtractors::SymmetricTensor<2> C_components;

      // 这是一个自动微分助手的实例，我们将设置它来完成与构成法则有关的所有微分计算......

      ADHelper ad_helper;

      // ... 以下三个成员变量将存储来自 @p ad_helper. 的输出。  @p ad_helper 一次性返回关于所有场变量的导数，因此我们将保留完整的梯度向量和Hessian矩阵。我们将从中提取我们真正感兴趣的单个条目。

      double             psi;
      Vector<double>     Dpsi;
      FullMatrix<double> D2psi;
    };

    // 在设置字段组件提取器时，对于它们的顺序是完全任意的。但重要的是，这些提取器没有重叠的索引。这些提取器的组件总数定义了 @p ad_helper 需要跟踪的独立变量的数量，并且我们将对其进行导数。由此产生的数据结构 @p Dpsi 和 @p D2psi 也必须有相应的大小。一旦 @p ad_helper 被配置好（它的输入参数是 $\mathbf{C}$ 和 $\boldsymbol{\mathbb{H}}$ 的组件总数），我们就可以直接询问它使用多少个独立变量。

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::
      Magnetoelastic_Constitutive_Law_AD(
        const ConstitutiveParameters &constitutive_parameters)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , H_components(0)
      , C_components(Tensor<1, dim>::n_independent_components)
      , ad_helper(Tensor<1, dim>::n_independent_components +
                  SymmetricTensor<2, dim>::n_independent_components)
      , psi(0.0)
      , Dpsi(ad_helper.n_independent_variables())
      , D2psi(ad_helper.n_independent_variables(),
              ad_helper.n_independent_variables())
    {}

    // 如前所述，由于自动微分库的工作方式， @p ad_helper 将总是同时返回能量密度函数相对于所有场变量的导数。由于这个原因，在函数`get_B()`、`get_S()`等中计算导数是没有意义的，因为我们会做很多额外的计算，然后直接丢弃。因此，处理这个问题的最好方法是用一个单一的函数调用来完成所有的前期计算，然后我们在需要时提取存储的数据。这就是我们在 "update_internal_data() "方法中要做的。由于材料是与速率无关的，我们可以忽略DiscreteTime参数。

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    void
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &)
    {
      Assert(determinant(C) > 0, ExcInternalError());

      // 由于我们在每个时间步骤中都会重复使用 @p ad_helper 数据结构，所以我们需要在使用前清除它的所有陈旧信息。

      ad_helper.reset();

      // 下一步是设置所有字段组件的值。这些定义了 "点"，我们将围绕这个点计算函数梯度及其线性化。我们之前创建的提取器提供了字段和 @p ad_helper 中的注册表之间的关联 -- 它们将被反复使用，以确保我们对哪个变量对应于`H`或`C`的哪个分量有正确的解释。

      ad_helper.register_independent_variable(H, H_components);
      ad_helper.register_independent_variable(C, C_components);

      // 现在我们已经完成了最初的设置，我们可以检索我们字段的AD对应关系。这些是真正的能量函数的独立变量，并且对用它们进行的计算是
      // "敏感的"。请注意，AD数被视为一种特殊的数字类型，可以在许多模板化的类中使用（在这个例子中，作为Tensor和SymmetricTensor类的标量类型）。

      const Tensor<1, dim, ADNumberType> H_ad =
        ad_helper.get_sensitive_variables(H_components);
      const SymmetricTensor<2, dim, ADNumberType> C_ad =
        ad_helper.get_sensitive_variables(C_components);

      // 我们还可以在许多以标量类型为模板的函数中使用它们。因此，对于我们需要的这些中间值，我们可以进行张量运算和一些数学函数。由此产生的类型也将是一个自动可分的数字，它对这些函数中的操作进行编码。

      const ADNumberType det_F_ad = std::sqrt(determinant(C_ad));
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      AssertThrow(det_F_ad > ADNumberType(0.0),
                  ExcMessage("Volumetric Jacobian must be positive."));

      // 接下来我们将计算出在磁场影响下导致剪切模量变化（增加）的比例函数......

      const ADNumberType f_mu_e_ad =
        1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
                std::tanh((2.0 * H_ad * H_ad) /
                          (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()));

      // ...然后我们就可以定义材料的储能密度函数。我们将在后面看到，这个例子足够复杂，值得使用AD，至少可以验证一个无辅助的实现。

      const ADNumberType psi_ad =
        0.5 * this->get_mu_e() * f_mu_e_ad *
          (trace(C_ad) - dim - 2.0 * std::log(det_F_ad))                 //
        + this->get_lambda_e() * std::log(det_F_ad) * std::log(det_F_ad) //
        - 0.5 * this->get_mu_0() * this->get_mu_r() * det_F_ad *
            (H_ad * C_inv_ad * H_ad); //

      // 储存的能量密度函数实际上是这个问题的因变量，所以作为 "配置 "阶段的最后一步，我们用 @p ad_helper. 注册其定义。
      ad_helper.register_dependent_variable(psi_ad);

      // 最后，我们可以检索存储的能量密度函数的结果值，以及它相对于输入字段的梯度和Hessian，并将它们缓存起来。

      psi = ad_helper.compute_value();
      ad_helper.compute_gradient(Dpsi);
      ad_helper.compute_hessian(D2psi);
    }

    // 下面的几个函数可以查询 $\psi_{0}$
    // 的存储值，并提取梯度向量和Hessian矩阵的所需成分。我们再次利用提取器来表达我们希望检索的总梯度向量和Hessian矩阵的哪些部分。它们只返回能量函数的导数，所以对于我们的动能变量的定义和它们的线性化，还需要进行一些操作来形成所需的结果。

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    double
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_psi() const
    {
      return psi;
    }

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    Tensor<1, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_B() const
    {
      const Tensor<1, dim> dpsi_dH =
        ad_helper.extract_gradient_component(Dpsi, H_components);
      return -dpsi_dH;
    }

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    SymmetricTensor<2, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_S() const
    {
      const SymmetricTensor<2, dim> dpsi_dC =
        ad_helper.extract_gradient_component(Dpsi, C_components);
      return 2.0 * dpsi_dC;
    }

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    SymmetricTensor<2, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_DD() const
    {
      const Tensor<2, dim> dpsi_dH_dH =
        ad_helper.extract_hessian_component(D2psi, H_components, H_components);
      return -symmetrize(dpsi_dH_dH);
    }

    // 请注意，对于耦合项来说，提取器参数的顺序特别重要，因为它决定了定向导数的提取顺序。因此，如果我们在调用`extract_hessian_component()`时颠倒了提取器的顺序，那么我们实际上是在检索
    // $\left[ \mathfrak{P}^{\text{tot}} \right]^{T}$  的一部分。

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    Tensor<3, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_PP() const
    {
      const Tensor<3, dim> dpsi_dC_dH =
        ad_helper.extract_hessian_component(D2psi, C_components, H_components);
      return -2.0 * dpsi_dC_dH;
    }

    template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
    SymmetricTensor<4, dim>
    Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_HH() const
    {
      const SymmetricTensor<4, dim> dpsi_dC_dC =
        ad_helper.extract_hessian_component(D2psi, C_components, C_components);
      return 4.0 * dpsi_dC_dC;
    }
    // @sect4{Magneto-viscoelastic constitutive law (using symbolic algebra and
    // differentiation)}

    // 我们要考虑的第二个材料定律将是一个代表具有单一耗散机制的磁涡弹材料。我们将考虑这种材料的自由能密度函数定义为
    // @f{align*}{
    //    \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
    //    \right)
    //  &= \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  + \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)
    //  \\ \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  &= \frac{1}{2} \mu_{e} f_{\mu_{e}^{ME}} \left( \boldsymbol{\mathbb{H}}
    //  \right)
    //      \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //      \right]
    //  + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    //  - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //      \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //      \boldsymbol{\mathbb{H}} \right]
    //  \\ \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)
    //  &= \frac{1}{2} \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //  \right)
    //      \left[ \mathbf{C}_{v} : \left[
    //        \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C} \right] - d - \ln\left(
    //        \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    //  @f}
    //  ，其中
    //  @f[
    //    f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)
    //  = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //      \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //      \boldsymbol{\mathbb{H}}}
    //        {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //  @f]

    // @f[
    //    f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)
    //  = 1 + \left[ \frac{\mu_{v}^{\infty}}{\mu_{v}} - 1 \right]
    //      \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //      \boldsymbol{\mathbb{H}}}
    //        {\left(h_{v}^{\text{sat}}\right)^{2}} \right),
    //  @f]
    //  与内部粘性变量
    //  @f[
    //  \mathbf{C}_{v}^{(t)}
    //  = \frac{1}{1 + \frac{\Delta t}{\tau_{v}}} \left[
    //      \mathbf{C}_{v}^{(t-1)}
    //    + \frac{\Delta t}{\tau_{v}}
    //      \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //      \mathbf{C} \right]^{-1}
    //    \right]
    //  @f]
    //  的演化规律相结合，该演化规律采用一阶后向差分近似法进行离散。

    // 再一次，让我们看看在一个具体的类中是如何实现的。我们现在将利用SD方法，而不是之前类中使用的AD框架。为了支持这一点，这个类的构造函数不仅接受 @p constitutive_parameters, ，而且还接受两个额外的变量，这些变量将被用来初始化一个 Differentiation::SD::BatchOptimizer. 我们将在后面给出更多的背景。

    template <int dim>
    class Magnetoviscoelastic_Constitutive_Law_SD final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
    public:
      Magnetoviscoelastic_Constitutive_Law_SD(
        const ConstitutiveParameters &               constitutive_parameters,
        const Differentiation::SD::OptimizerType     optimizer_type,
        const Differentiation::SD::OptimizationFlags optimization_flags);

      // 和自动区分助手一样， Differentiation::SD::BatchOptimizer
      // 将一次性返回一个结果集合。因此，为了只做一次，我们将利用与之前类似的方法，在`update_internal_data()`函数中做所有昂贵的计算，并将结果缓存起来，以便分层提取。

      virtual void
      update_internal_data(const SymmetricTensor<2, dim> &C,
                           const Tensor<1, dim> &         H,
                           const DiscreteTime &           time) override;

      virtual double
      get_psi() const override;

      virtual Tensor<1, dim>
      get_B() const override;

      virtual SymmetricTensor<2, dim>
      get_S() const override;

      virtual SymmetricTensor<2, dim>
      get_DD() const override;

      virtual Tensor<3, dim>
      get_PP() const override;

      virtual SymmetricTensor<4, dim>
      get_HH() const override;

      // 因为我们要处理的是一个与速率有关的材料，所以我们必须在适当的时候更新历史变量。这将是这个函数的目的。

      virtual void
      update_end_of_timestep() override;

      // 在该类的`private`部分，我们将希望跟踪内部的粘性变形，所以下面两个（实值的，非符号的）成员变量分别持有

      // - 内部变量时间步长（如果嵌入非线性求解器框架，则为牛顿步长）的值，以及

      // - 内部变量在前一个时间步长的值。

      // （我们将这些变量标记为
      // "Q"，以便于识别；在计算的海洋中，不一定容易将`Cv`或`C_v`与`C`区分开来）。

    private:
      SymmetricTensor<2, dim> Q_t;
      SymmetricTensor<2, dim> Q_t1;

      // 由于我们将使用符号类型，我们需要定义一些符号变量，以便与框架一起使用。(它们都以
      // "SD "为后缀，以方便区分符号类型或表达式与实值类型或标量。)
      // 这可以在前面做一次（甚至有可能作为 "静态
      // "变量），以尽量减少与创建这些变量相关的开销。为了实现通用编程的终极目标，我们甚至可以用符号来描述构成参数，*有可能*允许一个类的实例在这些值的不同输入下被重复使用。

      // 这些是代表弹性、粘性和磁性材料参数的符号标量（定义的顺序与它们在 @p ConstitutiveParameters 类中出现的顺序基本相同）。我们还存储了一个符号表达式， @p delta_t_sd, ，表示时间步长）。)

      const Differentiation::SD::Expression mu_e_sd;
      const Differentiation::SD::Expression mu_e_inf_sd;
      const Differentiation::SD::Expression mu_e_h_sat_sd;
      const Differentiation::SD::Expression lambda_e_sd;
      const Differentiation::SD::Expression mu_v_sd;
      const Differentiation::SD::Expression mu_v_inf_sd;
      const Differentiation::SD::Expression mu_v_h_sat_sd;
      const Differentiation::SD::Expression tau_v_sd;
      const Differentiation::SD::Expression delta_t_sd;
      const Differentiation::SD::Expression mu_r_sd;

      // 接下来我们定义一些代表独立场变量的张量符号变量，在此基础上，能量密度函数被参数化。

      const Tensor<1, dim, Differentiation::SD::Expression>          H_sd;
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_sd;

      // 同样，我们也有内部粘性变量的符号表示（包括它的当前值和它在前一个时间段的值）。

      const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t_sd;
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t1_sd;

      // 我们还应该存储从属表达式的定义。虽然我们只计算一次，但我们需要它们从下面声明的 @p optimizer 中检索数据。此外，当序列化一个像这样的材料类时（不是作为本教程的一部分），我们要么需要把这些表达式也序列化，要么需要在重新加载时重建它们。

      Differentiation::SD::Expression                          psi_sd;
      Tensor<1, dim, Differentiation::SD::Expression>          B_sd;
      SymmetricTensor<2, dim, Differentiation::SD::Expression> S_sd;
      SymmetricTensor<2, dim, Differentiation::SD::Expression> BB_sd;
      Tensor<3, dim, Differentiation::SD::Expression>          PP_sd;
      SymmetricTensor<4, dim, Differentiation::SD::Expression> HH_sd;

      // 然后，下一个变量是用于评估从属函数的优化器。更具体地说，它提供了加速评估符号依赖表达式的可能性。这是一个重要的工具，因为对冗长表达式的本地评估（不使用加速方法，而是直接对符号表达式进行评估）会非常慢。
      // Differentiation::SD::BatchOptimizer
      // 类提供了一种机制，可以将符号表达式树转化为另一种代码路径，例如，在各种从属表达式之间共享中间结果（意味着这些中间值每次评估只计算一次）和/或使用即时编译器编译代码（从而检索评估步骤的接近原生性能）。

      // 执行这种代码转换在计算上是非常昂贵的，所以我们存储了优化器，使其在每个类实例中只做一次。这也进一步促使我们决定将构成参数本身变成符号化。这样我们就可以在几种材料（当然是相同的能量函数）和潜在的多个连续体点（如果嵌入到有限元模拟中）中重复使用这个 @p optimizer 的单一实例。

      // 正如模板参数所指定的，数值结果将是<tt>double</tt>类型。

      Differentiation::SD::BatchOptimizer<double> optimizer;

      // 在评估阶段，我们必须将符号变量映射到它们的实值对应物。下一个方法将提供这个功能。

      // 这个类的最后一个方法将配置  @p optimizer.  。
      Differentiation::SD::types::substitution_map
      make_substitution_map(const SymmetricTensor<2, dim> &C,
                            const Tensor<1, dim> &         H,
                            const double                   delta_t) const;

      void
      initialize_optimizer();
    };

    // 由于静止变形状态是材料被认为是完全松弛的状态，内部粘性变量被初始化为同一张量，即
    // $\mathbf{C}_{v} = \mathbf{I}$
    // 。代表构成参数、时间步长、场和内部变量的各种符号变量都有一个唯一的标识符。优化器被传递给两个参数，这两个参数声明了应该应用哪种优化（加速）技术，以及CAS应该采取哪些额外步骤来帮助提高评估期间的性能。

    template <int dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::
      Magnetoviscoelastic_Constitutive_Law_SD(
        const ConstitutiveParameters &               constitutive_parameters,
        const Differentiation::SD::OptimizerType     optimizer_type,
        const Differentiation::SD::OptimizationFlags optimization_flags)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , Q_t(Physics::Elasticity::StandardTensors<dim>::I)
      , Q_t1(Physics::Elasticity::StandardTensors<dim>::I)
      , mu_e_sd("mu_e")
      , mu_e_inf_sd("mu_e_inf")
      , mu_e_h_sat_sd("mu_e_h_sat")
      , lambda_e_sd("lambda_e")
      , mu_v_sd("mu_v")
      , mu_v_inf_sd("mu_v_inf")
      , mu_v_h_sat_sd("mu_v_h_sat")
      , tau_v_sd("tau_v")
      , delta_t_sd("delta_t")
      , mu_r_sd("mu_r")
      , H_sd(Differentiation::SD::make_vector_of_symbols<dim>("H"))
      , C_sd(Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("C"))
      , Q_t_sd(
          Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t"))
      , Q_t1_sd(
          Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t1"))
      , optimizer(optimizer_type, optimization_flags)
    {
      initialize_optimizer();
    }

    // 替换图只是将以下所有数据配对在一起。

    // - 构成参数（从基类中获取的值）。

    // - 时间步长（从时间离散器中获取其值）。

    // 场值（其值由调用该 @p Magnetoviscoelastic_Constitutive_Law_SD 实例的外部函数规定），以及

    // 当前和之前的内部粘性变形（其值存储在这个类实例中）。

    template <int dim>
    Differentiation::SD::types::substitution_map
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::make_substitution_map(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const double                   delta_t) const
    {
      return Differentiation::SD::make_substitution_map(
        std::make_pair(mu_e_sd, this->get_mu_e()),
        std::make_pair(mu_e_inf_sd, this->get_mu_e_inf()),
        std::make_pair(mu_e_h_sat_sd, this->get_mu_e_h_sat()),
        std::make_pair(lambda_e_sd, this->get_lambda_e()),
        std::make_pair(mu_v_sd, this->get_mu_v()),
        std::make_pair(mu_v_inf_sd, this->get_mu_v_inf()),
        std::make_pair(mu_v_h_sat_sd, this->get_mu_v_h_sat()),
        std::make_pair(tau_v_sd, this->get_tau_v()),
        std::make_pair(delta_t_sd, delta_t),
        std::make_pair(mu_r_sd, this->get_mu_r()),
        std::make_pair(H_sd, H),
        std::make_pair(C_sd, C),
        std::make_pair(Q_t_sd, Q_t),
        std::make_pair(Q_t1_sd, Q_t1));
    }

    // 由于符号表达式的 "自然 "使用，配置 @p optimizer 的大部分过程看起来与构建自动区分帮助器的过程非常相似。尽管如此，我们还是要再次详细说明这些步骤，以强调这两个框架的不同之处。

    // 该函数从符号化编码变形梯度行列式的表达式开始（用右Cauchy-Green变形张量表示，即我们的主要场变量），以及
    // $\mathbf{C}$ 本身的逆。

    template <int dim>
    void
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::initialize_optimizer()
    {
      const Differentiation::SD::Expression det_F_sd =
        std::sqrt(determinant(C_sd));
      const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_inv_sd =
        invert(C_sd);

      // 接下来是自由能密度函数的弹性部分的饱和函数的符号表示，然后是自由能密度函数的磁弹性贡献。这一切都与我们之前看到的结构相同。

      const Differentiation::SD::Expression f_mu_e_sd =
        1.0 +
        (mu_e_inf_sd / mu_e_sd - 1.0) *
          std::tanh((2.0 * H_sd * H_sd) / (mu_e_h_sat_sd * mu_e_h_sat_sd));

      const Differentiation::SD::Expression psi_ME_sd =
        0.5 * mu_e_sd * f_mu_e_sd *
          (trace(C_sd) - dim - 2.0 * std::log(det_F_sd)) +
        lambda_e_sd * std::log(det_F_sd) * std::log(det_F_sd) -
        0.5 * this->get_mu_0() * mu_r_sd * det_F_sd * (H_sd * C_inv_sd * H_sd);

      // 此外，我们还定义了自由能密度函数的磁-粘弹性贡献。实现这一点所需的第一个组件是一个缩放函数，它将使粘性剪切模量在磁场影响下发生变化（增加）（见
      // @cite Pelteret2018a
      // ，公式29）。此后，我们可以计算能量密度函数的耗散分量；其表达式见 @cite
      // Pelteret2018a （公式28），这是对 @cite Linder2011a
      // （公式46）中提出的能量密度函数的直接扩展。

      const Differentiation::SD::Expression f_mu_v_sd =
        1.0 +
        (mu_v_inf_sd / mu_v_sd - 1.0) *
          std::tanh((2.0 * H_sd * H_sd) / (mu_v_h_sat_sd * mu_v_h_sat_sd));

      const Differentiation::SD::Expression psi_MVE_sd =
        0.5 * mu_v_sd * f_mu_v_sd *
        (Q_t_sd * (std::pow(det_F_sd, -2.0 / dim) * C_sd) - dim -
         std::log(determinant(Q_t_sd)));

      // 从这些构件中，我们可以定义材料的总自由能密度函数。

      psi_sd = psi_ME_sd + psi_MVE_sd;

      // 目前，对中科院来说，变量 @p Q_t_sd 似乎是独立于 @p C_sd. 的，我们的张量符号表达式 @p Q_t_sd 只是有一个与之相关的标识符，没有任何东西将其与另一个张量符号表达式 @p C_sd. 联系起来。因此，相对于 @p C_sd 的任何导数将忽略这种内在的依赖关系，正如我们从进化规律可以看到的，实际上是 $\mathbf{C}_{v} = \mathbf{C}_{v} \left( \mathbf{C}, t \right)$  。这意味着，相对于 $\mathbf{C}$ 推导任何函数 $f = f(\mathbf{C}, \mathbf{Q})$ 将返回部分导数 $\frac{\partial f(\mathbf{C}, \mathbf{Q})}{\partial \mathbf{C}}
      //  \Big\vert_{\mathbf{Q}}$ ，而不是总导数 $\frac{d f(\mathbf{C},
      //  \mathbf{Q}(\mathbf{C}))}{d \mathbf{C}} = \frac{\partial f(\mathbf{C},
      //  \mathbf{Q}(\mathbf{C}))}{\partial \mathbf{C}} \Big\vert_{\mathbf{Q}} +
      //  \frac{\partial f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))}{\partial
      //  \mathbf{Q}} \Big\vert_{\mathbf{C}} : \frac{d
      //  \mathbf{Q}(\mathbf{C}))}{d \mathbf{C}}$  。

      // 相比之下，在当前的AD库中，总导数将总是被返回。这意味着对于这类材料模型来说，计算出的动能变量是不正确的，这使得AD成为从能量密度函数中推导出（连续点水平）这种耗散性材料的构成法的不正确工具。

      // 正是这种特定的控制水平描述了SD和AD框架之间的一个决定性差异。在几行中，我们将对内部变量 @p Q_t_sd 的表达式进行操作，使其产生正确的线性化。
      //但是，
      //首先，我们将计算动能变量的符号表达式，即磁感应向量和Piola-Kirchhoff应力张量。执行微分的代码相当接近于模仿理论中所述的定义。

      B_sd = -Differentiation::SD::differentiate(psi_sd, H_sd);
      S_sd = 2.0 * Differentiation::SD::differentiate(psi_sd, C_sd);

      // 因为下一步是对上述内容进行线性化，所以现在是告知CAS  @p Q_t_sd  对  @p C_sd,  的明确依赖性的适当时机，即说明  $\mathbf{C}_{v} = \mathbf{C}_{v} \left( \mathbf{C}, t\right)$  。这意味着未来所有关于 @p C_sd 的微分运算将考虑到这种依赖关系（即计算总导数）。换句话说，我们将转换一些表达式，使其内在参数化从 $f(\mathbf{C}, \mathbf{Q})$ 变为 $f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))$  .

      // 为了做到这一点，我们考虑时间离散的演化规律。由此，我们有了内部变量在其历史上的明确表达，以及主要场变量。这就是它在这个表达式中描述的内容。

      const SymmetricTensor<2, dim, Differentiation::SD::Expression>
        Q_t_sd_explicit =
          (1.0 / (1.0 + delta_t_sd / tau_v_sd)) *
          (Q_t1_sd +
           (delta_t_sd / tau_v_sd * std::pow(det_F_sd, 2.0 / dim) * C_inv_sd));

      // 接下来我们产生一个中间替换图，它将在一个表达式中找到 @p Q_t_sd （我们的标识符）的每个实例，并用 @p Q_t_sd_explicit. 中的完整表达式来替换它。
      const Differentiation::SD::types::substitution_map
        substitution_map_explicit = Differentiation::SD::make_substitution_map(
          std::make_pair(Q_t_sd, Q_t_sd_explicit));

      // 我们可以在两个动力学变量上进行这种替换，并立即将替换后的结果与场变量进行区分。(如果你愿意，这可以分成两步进行，中间的结果储存在一个临时变量中)。同样，如果你忽略了代换所产生的
      // "复杂性"，这些将运动变量线性化并产生三个切向张量的调用与理论中所述的非常相似。

      BB_sd = symmetrize(Differentiation::SD::differentiate(
        Differentiation::SD::substitute(B_sd, substitution_map_explicit),
        H_sd));
      PP_sd = -Differentiation::SD::differentiate(
        Differentiation::SD::substitute(S_sd, substitution_map_explicit), H_sd);
      HH_sd =
        2.0 *
        Differentiation::SD::differentiate(
          Differentiation::SD::substitute(S_sd, substitution_map_explicit),
          C_sd);

      // 现在我们需要告诉 @p optimizer 我们需要提供哪些条目的数值，以便它能成功地进行计算。这些基本上充当了 @p optimizer 必须评估的所有从属函数的输入参数。它们统称为问题的自变量、历史变量、时间步长和构成参数（因为我们没有在能量密度函数中硬编码它们）。

      // 因此，我们真正想要的是为它提供一个符号集合，我们可以这样完成。
      // @code
      //  optimizer.register_symbols(Differentiation::SD::make_symbol_map(
      //    mu_e_sd, mu_e_inf_sd, mu_e_h_sat_sd, lambda_e_sd,
      //    mu_v_sd, mu_v_inf_sd, mu_v_h_sat_sd, tau_v_sd,
      //    delta_t_sd, mu_r_sd,
      //    H_sd, C_sd,
      //    Q_t_sd, Q_t1_sd));
      //  @endcode
      //  但这实际上都已经被编码为替换图的键。这样做还意味着我们需要在两个地方（这里和构建替换图时）管理这些符号，这很烦人，而且如果这个材料类被修改或扩展，可能会出现错误。由于我们此时对数值不感兴趣，所以如果替换图中与每个键项相关的数值被填入无效的数据也没有关系。所以我们将简单地创建一个假的替换图，并从中提取符号。请注意，任何传递给 @p optimizer 的替换图都必须至少包含这些符号的条目。

      optimizer.register_symbols(
        Differentiation::SD::Utilities::extract_symbols(
          make_substitution_map({}, {}, 0)));

      // 然后我们通知优化器我们想要计算哪些数值，在我们的情况下，这包括所有的因变量（即能量密度函数及其各种导数）。

      optimizer.register_functions(psi_sd, B_sd, S_sd, BB_sd, PP_sd, HH_sd);

      // 最后一步是最终确定优化器。通过这个调用，它将确定一个等价的代码路径，一次性评估所有的从属函数，但计算成本比直接评估符号表达式时要低。注意：这是一个昂贵的调用，所以我们希望尽可能少地执行它。我们在类的构造函数中完成了这一过程，实现了每个类实例只被调用一次的目标。

      optimizer.optimize();
    }

    // 由于 @p optimizer 的配置是在前面完成的，所以每次我们想计算动能变量或它们的线性化（导数）时，要做的事情就很少了。

    template <int dim>
    void
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &           time)
    {
      // 为了更新内部历史变量，我们首先需要计算一些基本量，这一点我们之前已经看到了。我们还可以向时间离散器询问用于从上一个时间步长迭代到当前时间步长的时间步长。

      const double delta_t = this->get_delta_t(time);

      const double                  det_F = std::sqrt(determinant(C));
      const SymmetricTensor<2, dim> C_inv = invert(C);
      AssertThrow(det_F > 0.0,
                  ExcMessage("Volumetric Jacobian must be positive."));

      // 现在，我们可以按照演化规律给出的定义，结合所选择的时间离散化方案，更新（实值）内部粘性变形张量。

      Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) *
            (Q_t1 + (delta_t / this->get_tau_v()) * std::pow(det_F, 2.0 / dim) *
                      C_inv);

      // 接下来我们向优化器传递我们希望自变量、时间步长和（本调用隐含的）构成参数所代表的数值。

      const auto substitution_map = make_substitution_map(C, H, delta_t);

      // 在进行下一次调用时，用于（数值）评估从属函数的调用路径要比字典替换更快。

      optimizer.substitute(substitution_map);
    }

    // 在调用了`update_internal_data()`之后，从优化器中提取数据就有效了。在进行评估时，我们需要从优化器中提取数据的确切符号表达式。这意味着我们需要在优化器的生命周期内存储所有因变量的符号表达式（自然，对输入变量也有同样的暗示）。

    template <int dim>
    double
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_psi() const
    {
      return optimizer.evaluate(psi_sd);
    }

    template <int dim>
    Tensor<1, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_B() const
    {
      return optimizer.evaluate(B_sd);
    }

    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_S() const
    {
      return optimizer.evaluate(S_sd);
    }

    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_DD() const
    {
      return optimizer.evaluate(BB_sd);
    }

    template <int dim>
    Tensor<3, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_PP() const
    {
      return optimizer.evaluate(PP_sd);
    }

    template <int dim>
    SymmetricTensor<4, dim>
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_HH() const
    {
      return optimizer.evaluate(HH_sd);
    }

    // 当在时间上向前移动时，内部变量的 "当前 "状态瞬间定义了 "前
    // "时间段的状态。因此，我们记录历史变量的值，作为下一个时间步长的 "过去值
    // "使用。

    template <int dim>
    void
    Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_end_of_timestep()
    {
      Q_t1 = Q_t;
    }
    // @sect3{A more complex example (continued): Parameters and hand-derived
    // material classes}

    // 现在我们已经看到了AD和SD框架如何在定义这些构成法则方面做了大量的工作，为了验证，我们将手工实现相应的类，并对框架与本地实现做一些初步的基准测试。

    // 为了保证作者的理智，下面记录的（希望是准确的）是动能变量和它们的切线的完整定义，以及一些中间计算过程。由于构成法则类的结构和设计已经在前面概述过了，我们将略过它，只是在
    // "update_internal_data()
    // "方法的定义中对各阶段的计算进行划分。将导数计算（及其适度表达的变量名）与出现在类描述中的文档定义联系起来应该是很容易的。然而，我们将借此机会介绍两种实现构成法类的不同范式。第二种将比第一种提供更多的灵活性（从而使其更容易扩展，在作者看来），但要牺牲一些性能。

    //  @sect4{Magnetoelastic constitutive law (hand-derived)}

    // 从前面提到的储存能量中，对于这种磁弹性材料，定义为
    // @f[
    //    \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  = \frac{1}{2} \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //      \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //      \right]
    //  + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    //  - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //      \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //      \boldsymbol{\mathbb{H}} \right]
    //  @f]
    //  与
    //  @f[
    //   f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //  = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //      \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //      \boldsymbol{\mathbb{H}}}
    //        {\left(h_{e}^{\text{sat}}\right)^{2}} \right) ,
    //  \\ \text{det}(\mathbf{F}) = \sqrt{\text{det}(\mathbf{C})}
    //  @f]
    //  ，对应于磁感应向量和总Piola-Kirchhoff应力张量的第一导数是
    //  @f[
    //   \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //   \right)
    //  \dealcoloneq - \frac{d \psi_{0}}{d \boldsymbol{\mathbb{H}}}
    //  = - \frac{1}{2} \mu_{e} \left[ \text{tr}(\mathbf{C}) - d - 2 \ln
    //  (\text{det}(\mathbf{F}))
    //        \right] \frac{d f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}}
    //        \right)}{d \boldsymbol{\mathbb{H}}}
    //  + \mu_{0} \mu_{r} \text{det}(\mathbf{F}) \left[ \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}}
    //      \right]
    //  @f]

    // @f{align}
    //   \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //   \right)
    //  \dealcoloneq 2 \frac{d \psi_{0} \left( \mathbf{C},
    //  \boldsymbol{\mathbb{H}} \right)}{d \mathbf{C}}
    //  &= \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //      \left[ \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
    //      - 2 \frac{1}{\text{det}(\mathbf{F})}
    //      \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \right]
    //  + 4 \lambda_{e} \ln \left(\text{det}(\mathbf{F}) \right)
    //      \frac{1}{\text{det}(\mathbf{F})} \frac{d\,\text{det}(\mathbf{F})}{d
    //      \mathbf{C}}
    //  - \mu_{0} \mu_{r} \left[
    //      \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //      \boldsymbol{\mathbb{H}} \right] \frac{d\,\text{det}(\mathbf{F})}{d
    //      \mathbf{C}} + \text{det}(\mathbf{F}) \frac{d \left[
    //      \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //      \boldsymbol{\mathbb{H}}
    //        \right]}{d \mathbf{C}} \right]
    //  \\ &= \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //      \left[ \mathbf{I} - \mathbf{C}^{-1} \right]
    //  + 2 \lambda_{e} \ln \left(\text{det}(\mathbf{F}) \right) \mathbf{C}^{-1}
    //  - \mu_{0} \mu_{r} \left[
    //      \frac{1}{2}  \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //      \cdot \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F})
    //      \mathbf{C}^{-1}
    //  - \text{det}(\mathbf{F})
    //      \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //        \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //        \right]
    //  @f}
    //  与
    //  @f[
    //    \frac{d f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)}{d
    //    \boldsymbol{\mathbb{H}}}
    //  = \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //    \text{sech}^{2} \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //    \boldsymbol{\mathbb{H}}}
    //      {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //    \left[ \frac{4} {\left(h_{e}^{\text{sat}}\right)^{2}}
    //    \boldsymbol{\mathbb{H}} \right]
    //  @f] 。

    // @f[
    //    \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
    //  = \mathbf{I}
    //  \quad \text{(the second-order identity tensor)}
    //  @f]

    // @f[
    //    \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}}
    //  = \frac{1}{2} \text{det}(\mathbf{F}) \mathbf{C}^{-1}
    //  @f]

    // @f[
    //  \frac{d C^{-1}_{ab}}{d C_{cd}}
    //  = - \text{sym} \left( C^{-1}_{ac} C^{-1}_{bd} \right)
    //  = -\frac{1}{2} \left[ C^{-1}_{ac} C^{-1}_{bd} + C^{-1}_{ad} C^{-1}_{bc}
    //  \right]
    //  @f]

    // @f[
    //    \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //    \boldsymbol{\mathbb{H}} \right]}{d \mathbf{C}}
    //  = - \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //    \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //  @f]
    //  在上面的一个推导中使用对称算子 $\text{sym} \left( \bullet \right)$
    //  有助于确保所产生的秩-4张量，由于 $\mathbf{C}$
    //  的对称性而持有小的对称性，仍然将秩-2对称张量映射为秩-2对称张量。参见SymmetricTensor类文档和
    //  step-44 的介绍，并进一步解释在四阶张量的背景下对称性的含义。

    //每个运动学变量相对于其参数的线性化是
    // @f[
    //  \mathbb{D} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  = \frac{d \boldsymbol{\mathbb{B}}}{d \boldsymbol{\mathbb{H}}}
    //  = - \frac{1}{2} \mu_{e} \left[ \text{tr}(\mathbf{C}) - d - 2 \ln
    //  (\text{det}(\mathbf{F}))
    //      \right] \frac{d^{2} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}}
    //      \right)}{d \boldsymbol{\mathbb{H}} \otimes d
    //      \boldsymbol{\mathbb{H}}}
    //  + \mu_{0} \mu_{r} \text{det}(\mathbf{F}) \mathbf{C}^{-1}
    //  @f]

    // @f{align}
    //  \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right) = - \frac{d \mathbf{S}^{\text{tot}}}{d \boldsymbol{\mathbb{H}}}
    //  &= - \mu_{e}
    //      \left[ \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
    //      - 2 \frac{1}{\text{det}(\mathbf{F})}
    //      \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \right]
    //        \otimes \frac{d f_{\mu_{e} \left( \boldsymbol{\mathbb{H}}
    //        \right)}}{d \boldsymbol{\mathbb{H}}}
    //  + \mu_{0} \mu_{r} \left[
    //      \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \otimes
    //        \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //        \boldsymbol{\mathbb{H}}
    //          \right]}{d \boldsymbol{\mathbb{H}}} \right]
    //  + \text{det}(\mathbf{F})
    //      \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //      \cdot \boldsymbol{\mathbb{H}}
    //        \right]}{d \mathbf{C} \otimes d \boldsymbol{\mathbb{H}}}
    //  \\ &= - \mu_{e}
    //      \left[ \mathbf{I} - \mathbf{C}^{-1} \right] \otimes
    //        \frac{d f_{\mu_{e} \left( \boldsymbol{\mathbb{H}} \right)}}{d
    //        \boldsymbol{\mathbb{H}}}
    //  + \mu_{0} \mu_{r} \left[
    //      \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes
    //        \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //        \right]
    //  + \text{det}(\mathbf{F})
    //      \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //      \cdot \boldsymbol{\mathbb{H}}
    //        \right]}{d \mathbf{C} \otimes \mathbf{C} \boldsymbol{\mathbb{H}}}
    //  @f}


    // @f{align}
    //  \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right) = 2 \frac{d \mathbf{S}^{\text{tot}}}{d \mathbf{C}}
    //  &= 2 \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //      \left[ - \frac{d \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    //    + 4 \lambda_{e} \left[ \mathbf{C}^{-1} \otimes \left[
    //    \frac{1}{\text{det}(\mathbf{F})} \frac{d \, \text{det}(\mathbf{F})}{d
    //    \mathbf{C}} \right] + \ln \left(\text{det}(\mathbf{F}) \right) \frac{d
    //    \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    //  \\ &- \mu_{0} \mu_{r}  \left[
    //   \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes \frac{d \left[
    //   \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //   \boldsymbol{\mathbb{H}} \right]}{d \mathbf{C}}
    //  + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}} \right] \mathbf{C}^{-1} \otimes \frac{d \,
    //  \text{det}(\mathbf{F})}{d \mathbf{C}}
    //  + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F}) \frac{d
    //  \mathbf{C}^{-1}}{d \mathbf{C}}
    //  \right]
    //  \\ &+ 2 \mu_{0} \mu_{r} \left[ \left[
    //      \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //        \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //        \right] \otimes \frac{d \, \text{det}(\mathbf{F})}{d \mathbf{C}}
    //      - \text{det}(\mathbf{F})
    //      \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //      \cdot \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
    //      \mathbf{C}}
    //  \right]
    //  \\ &= 2 \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    //      \left[ - \frac{d \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    //   + 4 \lambda_{e} \left[ \frac{1}{2} \mathbf{C}^{-1} \otimes
    //   \mathbf{C}^{-1} + \ln \left(\text{det}(\mathbf{F}) \right) \frac{d
    //   \mathbf{C}^{-1}}{d \mathbf{C}} \right]
    //  \\ &- \mu_{0} \mu_{r}  \left[
    //   - \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes \left[ \left[
    //   \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //    \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \right]
    //  + \frac{1}{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F})  \mathbf{C}^{-1}
    //  \otimes \mathbf{C}^{-1}
    //  + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F}) \frac{d
    //  \mathbf{C}^{-1}}{d \mathbf{C}}
    //  \right]
    //  \\ &+ 2 \mu_{0} \mu_{r} \left[ \frac{1}{2} \text{det}(\mathbf{F}) \left[
    //      \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
    //        \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
    //        \right] \otimes \mathbf{C}^{-1}
    //      - \text{det}(\mathbf{F})
    //      \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
    //      \cdot \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
    //      \mathbf{C}}
    //  \right]
    //  @f}
    //  与
    //  @f[
    //   \frac{d^{2} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)}{d
    //   \boldsymbol{\mathbb{H}} \otimes d \boldsymbol{\mathbb{H}}}
    //  = -2 \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //    \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //    \boldsymbol{\mathbb{H}}}
    //      {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //    \text{sech}^{2} \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //    \boldsymbol{\mathbb{H}}}
    //      {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //    \left[ \frac{4} {\left(h_{e}^{\text{sat}}\right)^{2}} \mathbf{I}
    //    \right]
    //  @f]

    // @f[
    //  \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}}
    //          \right]}{d \boldsymbol{\mathbb{H}}}
    //  = 2 \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
    //  @f]

    // @f[
    //  \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
    //  \boldsymbol{\mathbb{H}}} \Rightarrow \frac{d^{2} \left[ \mathbb{H}_{e}
    //  C^{-1}_{ef} \mathbb{H}_{f}
    //        \right]}{d C_{ab} d \mathbb{H}_{c}}
    //  = - C^{-1}_{ac} C^{-1}_{be} \mathbb{H}_{e} - C^{-1}_{ae} \mathbb{H}_{e}
    //  C^{-1}_{bc}
    //  @f]

    // @f{align}
    //  \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //  \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d \mathbf{C}}
    //  &= -\frac{d \left[\left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
    //  \right] \otimes
    //        \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
    //        \right]\right]}{d \mathbf{C}}
    //  \\ \Rightarrow
    //  \frac{d^{2} \left[ \mathbb{H}_{e} C^{-1}_{ef} \mathbb{H}_{f}
    //        \right]}{d C_{ab} d C_{cd}}
    //  &= \text{sym} \left( C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{cf}
    //  \mathbb{H}_{f} C^{-1}_{bd}
    //            + C^{-1}_{ce} \mathbb{H}_{e} C^{-1}_{bf} \mathbb{H}_{f}
    //            C^{-1}_{ad} \right)
    //  \\ &= \frac{1}{2} \left[
    //       C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{cf} \mathbb{H}_{f} C^{-1}_{bd}
    //     + C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{df} \mathbb{H}_{f} C^{-1}_{bc}
    //     + C^{-1}_{ce} \mathbb{H}_{e} C^{-1}_{bf} \mathbb{H}_{f} C^{-1}_{ad}
    //     + C^{-1}_{be} \mathbb{H}_{e} C^{-1}_{df} \mathbb{H}_{f} C^{-1}_{ac}
    //    \right]
    //  @f}

    // 好吧，很快就升级了--尽管 $\psi_{0}$ 和 $f_{\mu_e}$
    // 的定义可能已经给出了一些提示，说明计算动能场和它们的线性化需要一些努力，但最终的定义可能比最初想象的要复杂一些。了解了我们现在所做的，也许可以说我们真的不想计算这些函数相对于其参数的一、二次导数--不管我们在微积分课上做得如何，或者我们可能是多么好的程序员。

    // 在最终实现这些的类方法定义中，我们以稍微不同的方式组成这些计算。一些中间步骤也被保留下来，以便从另一个角度说明如何系统地计算导数。此外，一些计算被分解得更少或更进一步，以重用一些中间值，并希望能帮助读者跟随导数的操作。

    template <int dim>
    class Magnetoelastic_Constitutive_Law final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
    public:
      Magnetoelastic_Constitutive_Law(
        const ConstitutiveParameters &constitutive_parameters);

      virtual void
      update_internal_data(const SymmetricTensor<2, dim> &C,
                           const Tensor<1, dim> &         H,
                           const DiscreteTime &) override;

      virtual double
      get_psi() const override;

      virtual Tensor<1, dim>
      get_B() const override;

      virtual SymmetricTensor<2, dim>
      get_S() const override;

      virtual SymmetricTensor<2, dim>
      get_DD() const override;

      virtual Tensor<3, dim>
      get_PP() const override;

      virtual SymmetricTensor<4, dim>
      get_HH() const override;

    private:
      double                  psi;
      Tensor<1, dim>          B;
      SymmetricTensor<2, dim> S;
      SymmetricTensor<2, dim> BB;
      Tensor<3, dim>          PP;
      SymmetricTensor<4, dim> HH;
    };

    template <int dim>
    Magnetoelastic_Constitutive_Law<dim>::Magnetoelastic_Constitutive_Law(
      const ConstitutiveParameters &constitutive_parameters)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , psi(0.0)
    {}

    // 对于这个类的更新方法，我们将简单地预先计算一个中间值的集合（用于函数求值、导数计算等），并
    // "手动
    // "安排它们的顺序，以使其重复使用最大化。这意味着我们必须自己管理，并决定哪些值必须在其他值之前计算，同时保持代码本身的某种秩序或结构的模样。这很有效，但也许有点乏味。它对类的未来扩展也没有太大的帮助，因为所有这些值都是这个单一方法的局部。

    // 有趣的是，这种预先计算在多个地方使用的中间表达式的基本技术有一个名字：[共同子表达式消除（CSE）]（https：en.wikipedia.org/wiki/Common_subexpression_elimination）。它是计算机代数系统在承担评估类似表达式的任务时用来减少计算费用的一种策略。

    template <int dim>
    void
    Magnetoelastic_Constitutive_Law<dim>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &)
    {
      const double                  det_F = std::sqrt(determinant(C));
      const SymmetricTensor<2, dim> C_inv = invert(C);
      AssertThrow(det_F > 0.0,
                  ExcMessage("Volumetric Jacobian must be positive."));

      // 磁弹性能的饱和函数。

      const double two_h_dot_h_div_h_sat_squ =
        (2.0 * H * H) / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat());
      const double tanh_two_h_dot_h_div_h_sat_squ =
        std::tanh(two_h_dot_h_div_h_sat_squ);

      const double f_mu_e =
        1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
                tanh_two_h_dot_h_div_h_sat_squ;

      // 饱和函数的一阶导数，注意到  $\frac{d \tanh(x)}{dx} =
      // \text{sech}^{2}(x)$  。

      const double dtanh_two_h_dot_h_div_h_sat_squ =
        std::pow(1.0 / std::cosh(two_h_dot_h_div_h_sat_squ), 2.0);
      const Tensor<1, dim> dtwo_h_dot_h_div_h_sat_squ_dH =
        2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) * H;

      const Tensor<1, dim> df_mu_e_dH =
        (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
        (dtanh_two_h_dot_h_div_h_sat_squ * dtwo_h_dot_h_div_h_sat_squ_dH);

      // 饱和度函数的二阶导数，注意  $\frac{d \text{sech}^{2}(x)}{dx} = -2
      // \tanh(x) \text{sech}^{2}(x)$  。

      const double d2tanh_two_h_dot_h_div_h_sat_squ =
        -2.0 * tanh_two_h_dot_h_div_h_sat_squ * dtanh_two_h_dot_h_div_h_sat_squ;
      const SymmetricTensor<2, dim> d2two_h_dot_h_div_h_sat_squ_dH_dH =
        2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) *
        Physics::Elasticity::StandardTensors<dim>::I;

      const SymmetricTensor<2, dim> d2f_mu_e_dH_dH =
        (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
        (d2tanh_two_h_dot_h_div_h_sat_squ *
           symmetrize(outer_product(dtwo_h_dot_h_div_h_sat_squ_dH,
                                    dtwo_h_dot_h_div_h_sat_squ_dH)) +
         dtanh_two_h_dot_h_div_h_sat_squ * d2two_h_dot_h_div_h_sat_squ_dH_dH);

      // 从场/运动学变量中直接获得的一些中间量。

      const double         log_det_F         = std::log(det_F);
      const double         tr_C              = trace(C);
      const Tensor<1, dim> C_inv_dot_H       = C_inv * H;
      const double         H_dot_C_inv_dot_H = H * C_inv_dot_H;

      // 中间量的一阶导数。

      const SymmetricTensor<2, dim> d_tr_C_dC =
        Physics::Elasticity::StandardTensors<dim>::I;
      const SymmetricTensor<2, dim> ddet_F_dC     = 0.5 * det_F * C_inv;
      const SymmetricTensor<2, dim> dlog_det_F_dC = 0.5 * C_inv;

      const Tensor<1, dim> dH_dot_C_inv_dot_H_dH = 2.0 * C_inv_dot_H;

      SymmetricTensor<4, dim> dC_inv_dC;
      for (unsigned int A = 0; A < dim; ++A)
        for (unsigned int B = A; B < dim; ++B)
          for (unsigned int C = 0; C < dim; ++C)
            for (unsigned int D = C; D < dim; ++D)
              dC_inv_dC[A][B][C][D] -=               //
                0.5 * (C_inv[A][C] * C_inv[B][D]     //
                       + C_inv[A][D] * C_inv[B][C]); //

      const SymmetricTensor<2, dim> dH_dot_C_inv_dot_H_dC =
        -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H));

      // 中间量的二阶导数。

      const SymmetricTensor<4, dim> d2log_det_F_dC_dC = 0.5 * dC_inv_dC;

      const SymmetricTensor<4, dim> d2det_F_dC_dC =
        0.5 * (outer_product(C_inv, ddet_F_dC) + det_F * dC_inv_dC);

      const SymmetricTensor<2, dim> d2H_dot_C_inv_dot_H_dH_dH = 2.0 * C_inv;

      Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH;
      for (unsigned int A = 0; A < dim; ++A)
        for (unsigned int B = 0; B < dim; ++B)
          for (unsigned int C = 0; C < dim; ++C)
            d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -=
              C_inv[A][C] * C_inv_dot_H[B] + //
              C_inv_dot_H[A] * C_inv[B][C];  //

      SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC;
      for (unsigned int A = 0; A < dim; ++A)
        for (unsigned int B = A; B < dim; ++B)
          for (unsigned int C = 0; C < dim; ++C)
            for (unsigned int D = C; D < dim; ++D)
              d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] +=
                0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] +
                       C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] +
                       C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] +
                       C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]);

      // 储存的能量密度函数。

      psi =
        (0.5 * this->get_mu_e() * f_mu_e) *
          (tr_C - dim - 2.0 * std::log(det_F)) +
        this->get_lambda_e() * (std::log(det_F) * std::log(det_F)) -
        (0.5 * this->get_mu_0() * this->get_mu_r()) * det_F * (H * C_inv * H);

      // 动能量。

      B = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) *
            df_mu_e_dH //
          + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
              dH_dot_C_inv_dot_H_dH; //

      S = 2.0 * (0.5 * this->get_mu_e() * f_mu_e) *                        //
            (d_tr_C_dC - 2.0 * dlog_det_F_dC)                              //
          + 2.0 * this->get_lambda_e() * (2.0 * log_det_F * dlog_det_F_dC) //
          - 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *            //
              (H_dot_C_inv_dot_H * ddet_F_dC                               //
               + det_F * dH_dot_C_inv_dot_H_dC);                           //

      // 动能量的线性化。

      BB = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * //
             d2f_mu_e_dH_dH                                             //
           + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
               d2H_dot_C_inv_dot_H_dH_dH; //

      PP = -2.0 * (0.5 * this->get_mu_e()) *                                  //
             outer_product(Tensor<2, dim>(d_tr_C_dC - 2.0 * dlog_det_F_dC),   //
                           df_mu_e_dH)                                        //
           +                                                                  //
           2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                //
             (outer_product(Tensor<2, dim>(ddet_F_dC), dH_dot_C_inv_dot_H_dH) //
              + det_F * d2H_dot_C_inv_dot_H_dC_dH);                           //

      HH =
        4.0 * (0.5 * this->get_mu_e() * f_mu_e) * (-2.0 * d2log_det_F_dC_dC) //
        + 4.0 * this->get_lambda_e() *                                       //
            (2.0 * outer_product(dlog_det_F_dC, dlog_det_F_dC)               //
             + 2.0 * log_det_F * d2log_det_F_dC_dC)                          //
        - 4.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                //
            (H_dot_C_inv_dot_H * d2det_F_dC_dC                               //
             + outer_product(ddet_F_dC, dH_dot_C_inv_dot_H_dC)               //
             + outer_product(dH_dot_C_inv_dot_H_dC, ddet_F_dC)               //
             + det_F * d2H_dot_C_inv_dot_H_dC_dC);                           //
    }

    template <int dim>
    double
    Magnetoelastic_Constitutive_Law<dim>::get_psi() const
    {
      return psi;
    }

    template <int dim>
    Tensor<1, dim>
    Magnetoelastic_Constitutive_Law<dim>::get_B() const
    {
      return B;
    }

    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoelastic_Constitutive_Law<dim>::get_S() const
    {
      return S;
    }

    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoelastic_Constitutive_Law<dim>::get_DD() const
    {
      return BB;
    }

    template <int dim>
    Tensor<3, dim>
    Magnetoelastic_Constitutive_Law<dim>::get_PP() const
    {
      return PP;
    }

    template <int dim>
    SymmetricTensor<4, dim>
    Magnetoelastic_Constitutive_Law<dim>::get_HH() const
    {
      return HH;
    }
    // @sect4{Magneto-viscoelastic constitutive law (hand-derived)}

    // 如前所述，我们将考虑的具有一种耗散机制的磁涡流材料的自由能密度函数定义为
    // @f[
    //    \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
    //    \right)
    //  = \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  + \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)
    //  @f] 。

    // @f[
    //    \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  = \frac{1}{2} \mu_{e} f_{\mu_{e}^{ME}} \left( \boldsymbol{\mathbb{H}}
    //  \right)
    //      \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    //      \right]
    //  + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
    //  - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    //      \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    //      \boldsymbol{\mathbb{H}} \right]
    //  @f]

    // @f[
    //    \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //    \boldsymbol{\mathbb{H}} \right)
    //  = \frac{1}{2} \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //  \right)
    //      \left[ \mathbf{C}_{v} : \left[
    //        \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C} \right] - d - \ln\left(
    //        \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    //  @f]
    //  与
    //  @f[
    //    f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)
    //  = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    //      \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //      \boldsymbol{\mathbb{H}}}
    //        {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
    //  @f]
    //
    // @f[
    //    f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)
    //  = 1 + \left[ \frac{\mu_{v}^{\infty}}{\mu_{v}} - 1 \right]
    //      \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    //      \boldsymbol{\mathbb{H}}}
    //        {\left(h_{v}^{\text{sat}}\right)^{2}} \right)
    //  @f]
    //  和演变规律
    //  @f[
    //   \dot{\mathbf{C}}_{v} \left( \mathbf{C} \right)
    //  = \frac{1}{\tau} \left[
    //        \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //          \mathbf{C}\right]^{-1}
    //      - \mathbf{C}_{v} \right]
    //  @f] ，
    //  其本身是以 $\mathbf{C}$ 为参数的。根据设计，能量 $\psi_{0}^{ME} \left(
    //  \mathbf{C}, \boldsymbol{\mathbb{H}} \right)$
    //  的磁弹性部分与前面介绍的磁弹性材料的磁弹性部分是相同的。因此，对于源于这部分能量的各种贡献的导数，请参考前面的章节。我们将继续强调来自这些条款的具体贡献，用
    //  $ME$ 对突出的条款进行上标，而来自磁弹性部分的贡献则用 $MVE$
    //  上标。此外，阻尼项的磁饱和函数 $f_{\mu_{v}}^{MVE} \left(
    //  \boldsymbol{\mathbb{H}} \right)$ 具有与弹性项相同的形式（即
    //  $f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)$
    //  ），因此其导数的结构与之前看到的相同；唯一的变化是三个构成参数，现在与粘性剪切模量
    //  $\mu_{v}$ 而非弹性剪切模量 $\mu_{e}$ 相关。

    // 对于这种磁-粘弹性材料，对应于磁感应矢量和Piola-Kirchhoff总应力张量的第一导数是
    // @f[
    //   \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \mathbf{C}_{v},
    //   \boldsymbol{\mathbb{H}} \right)
    //  \dealcoloneq - \frac{\partial \psi_{0} \left( \mathbf{C},
    //  \mathbf{C}_{v}, \boldsymbol{\mathbb{H}} \right)}{\partial
    //  \boldsymbol{\mathbb{H}}} \Big\vert_{\mathbf{C}, \mathbf{C}_{v}} \equiv
    //  \boldsymbol{\mathbb{B}}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right)
    //  + \boldsymbol{\mathbb{B}}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) =  - \frac{d \psi_{0}^{ME} \left(
    //  \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}}}
    //     - \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //     \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
    //  @f]

    // @f[
    //   \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
    //   \boldsymbol{\mathbb{H}} \right)
    //  \dealcoloneq 2 \frac{\partial \psi_{0} \left( \mathbf{C},
    //  \mathbf{C}_{v}, \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
    //  \Big\vert_{\mathbf{C}_{v}, \boldsymbol{\mathbb{H}}} \equiv
    //  \mathbf{S}^{\text{tot}, ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right)
    //  + \mathbf{S}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}}
    //      \right)
    //  =  2 \frac{d \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right)}{d \mathbf{C}}
    //   + 2 \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //   \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
    //  @f]
    //  ，其中粘性贡献为
    //  @f[
    //    \boldsymbol{\mathbb{B}}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //    \boldsymbol{\mathbb{H}} \right)
    //  = - \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
    //  \Big\vert_{\mathbf{C}, \mathbf{C}_{v}} = - \frac{1}{2} \mu_{v}
    //      \left[ \mathbf{C}_{v} : \left[
    //        \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C} \right] - d - \ln\left(
    //        \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    //        \frac{\partial f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //        \right)}{\partial \boldsymbol{\mathbb{H}}}
    //  @f]

    // @f[
    //    \mathbf{S}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //    \boldsymbol{\mathbb{H}}
    //      \right)
    //  = 2 \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
    //  \Big\vert_{\mathbf{C}_{v}, \boldsymbol{\mathbb{H}}} = \mu_{v}
    //  f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    //         \left[  \left[ \mathbf{C}_{v} : \mathbf{C} \right] \left[ -
    //         \frac{1}{d}
    //         \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //         \mathbf{C}^{-1} \right]
    //         + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //         \mathbf{C}_{v}
    //   \right]
    //  @f]
    //  和
    //  @f[
    //  \frac{\partial f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //  \right)}{\partial \boldsymbol{\mathbb{H}}} \equiv \frac{d
    //  f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)}{d
    //  \boldsymbol{\mathbb{H}}} .
    //  @f]
    //  时间微缩的演化规律，
    //  @f[
    //  \mathbf{C}_{v}^{(t)} \left( \mathbf{C} \right)
    //  = \frac{1}{1 + \frac{\Delta t}{\tau_{v}}} \left[
    //      \mathbf{C}_{v}^{(t-1)}
    //    + \frac{\Delta t}{\tau_{v}}
    //      \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //      \mathbf{C} \right]^{-1}
    //    \right]
    //  @f]
    //  也将决定内部变量相对于场变量的线性化是如何构成的。

    // 注意，为了获得这种耗散材料的磁感应矢量和总Piola-Kirchhoff应力张量的*正确表达式，我们必须严格遵守应用Coleman-Noll程序的结果：我们必须取*部分导数*。
    //自由能密度函数与场变量的关系。(对于我们的非耗散性磁弹性材料，取部分导数或全部导数都会有同样的结果，所以之前没有必要提请大家注意这一点)。操作的关键部分是冻结内部变量
    //$\mathbf{C}_{v}^{(t)} \left( \mathbf{C} \right)$ ，同时计算
    //$\psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v} \left( \mathbf{C}
    //\right), \boldsymbol{\mathbb{H}} \right)$ 相对于 $\mathbf{C}$ 的导数--
    //$\mathbf{C}_{v}^{(t)}$ 对 $\mathbf{C}$
    //的依赖性不被考虑。当决定是使用AD还是SD来执行这个任务时，选择是很清楚的--只有符号框架提供了一个机制来完成这个任务；如前所述，AD只能返回总导数，所以它不适合这个任务。

    // 为了对事情进行总结，我们将介绍这种速度依赖性耦合材料的材料切线。两个动能变量相对于其参数的线性化是
    // @f[
    //  \mathbb{D} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
    //  \right) = \frac{d \boldsymbol{\mathbb{B}}}{d \boldsymbol{\mathbb{H}}}
    //  \equiv \mathbb{D}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
    //  \right)
    //  + \mathbb{D}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) = \frac{d
    //  \boldsymbol{\mathbb{B}}^{ME}}{d \boldsymbol{\mathbb{H}}}
    //  + \frac{d \boldsymbol{\mathbb{B}}^{MVE}}{d \boldsymbol{\mathbb{H}}}
    //  @f]
    //
    // @f[
    //  \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) = - \frac{d \mathbf{S}^{\text{tot}}}{d
    //  \boldsymbol{\mathbb{H}}} \equiv \mathfrak{P}^{\text{tot}, ME} \left(
    //  \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
    //  + \mathfrak{P}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) = - \frac{d \mathbf{S}^{\text{tot},
    //  ME}}{d \boldsymbol{\mathbb{H}}}
    //  - \frac{d \mathbf{S}^{\text{tot}, MVE}}{d \boldsymbol{\mathbb{H}}}
    //  @f]
    //
    // @f[
    //  \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) = 2 \frac{d \mathbf{S}^{\text{tot}}}{d
    //  \mathbf{C}} \equiv \mathcal{H}^{\text{tot}, ME} \left( \mathbf{C},
    //  \boldsymbol{\mathbb{H}} \right)
    //  + \mathcal{H}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) = 2 \frac{d \mathbf{S}^{\text{tot},
    //  ME}}{d \mathbf{C}}
    //  + 2 \frac{d \mathbf{S}^{\text{tot}, MVE}}{d \mathbf{C}}
    //  @f]
    //  其中粘性贡献的切线为
    //  @f[
    //  \mathbb{D}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) = - \frac{1}{2} \mu_{v}
    //      \left[ \mathbf{C}_{v} : \left[
    //        \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //        \mathbf{C} \right] - d - \ln\left(
    //        \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
    //        \frac{\partial^{2} f_{\mu_{v}^{MVE}} \left(
    //        \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}
    //        \otimes d \boldsymbol{\mathbb{H}}}
    //  @f]
    //
    // @f[
    //  \mathfrak{P}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right) = - \mu_{v}
    //         \left[  \left[ \mathbf{C}_{v} : \mathbf{C} \right] \left[ -
    //         \frac{1}{d}
    //         \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //         \mathbf{C}^{-1} \right]
    //         + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //         \mathbf{C}_{v}
    //   \right] \otimes \frac{d f_{\mu_{v}^{MVE}} \left(
    //   \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}}}
    //  @f]
    //
    // @f{align}
    //  \mathcal{H}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
    //  \boldsymbol{\mathbb{H}} \right)
    //  &= 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    //    \left[ - \frac{1}{d}
    //    \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //    \mathbf{C}^{-1} \right] \otimes
    //    \left[ \mathbf{C}_{v} + \mathbf{C} : \frac{d \mathbf{C}_{v}}{d
    //    \mathbf{C}} \right]
    //  \\ &+ 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    //  \left[ \mathbf{C}_{v} : \mathbf{C} \right]
    //    \left[
    //      \frac{1}{d^{2}}
    //      \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //      \mathbf{C}^{-1} \otimes \mathbf{C}^{-1}
    //      - \frac{1}{d}
    //      \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //      \frac{d \mathbf{C}^{-1}}{d \mathbf{C}}
    //    \right]
    //  \\ &+ 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
    //    \left[
    //      -\frac{1}{d}
    //      \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //      \mathbf{C}_{v} \otimes \mathbf{C}^{-1}
    //      + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    //      \frac{d \mathbf{C}_{v}}{d \mathbf{C}}
    //    \right]
    //  @f}
    //  与
    //  @f[
    //  \frac{\partial^{2} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
    //  \right)}{\partial \boldsymbol{\mathbb{H}} \otimes
    //  d \boldsymbol{\mathbb{H}}} \equiv \frac{d^{2} f_{\mu_{v}^{MVE}} \left(
    //  \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}} \otimes d
    //  \boldsymbol{\mathbb{H}}}
    //  @f]
    //  ，从演化定律来看，
    // //  @f[
    //  \frac{d \mathbf{C}_{v}}{d \mathbf{C}}
    //  \equiv \frac{d \mathbf{C}_{v}^{(t)}}{d \mathbf{C}}
    //   = \frac{\frac{\Delta t}{\tau_{v}} }{1 + \frac{\Delta t}{\tau_{v}}}
    //   \left[
    //      \frac{1}{d}
    //      \left[\text{det}\left(\mathbf{F}\right)\right]^{\frac{2}{d}}
    //      \mathbf{C}^{-1} \otimes \mathbf{C}^{-1}
    //     + \left[\text{det}\left(\mathbf{F}\right)\right]^{\frac{2}{d}}
    //     \frac{d \mathbf{C}^{-1}}{d \mathbf{C}}
    //    \right] .
    //  @f]
    //  注意，只是 $\mathcal{H}^{\text{tot}, MVE}$
    //  的最后一项包含内部变量的切线。这个特殊演化规律的线性化是线性的。关于非线性演化定律的例子，这种线性化必须以迭代的方式求解，见
    //  @cite Koprowski  -Theiss2011a。

    template <int dim>
    class Magnetoviscoelastic_Constitutive_Law final
      : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
    {
    public:
      Magnetoviscoelastic_Constitutive_Law(
        const ConstitutiveParameters &constitutive_parameters);

      virtual void
      update_internal_data(const SymmetricTensor<2, dim> &C,
                           const Tensor<1, dim> &         H,
                           const DiscreteTime &           time) override;

      virtual double
      get_psi() const override;

      virtual Tensor<1, dim>
      get_B() const override;

      virtual SymmetricTensor<2, dim>
      get_S() const override;

      virtual SymmetricTensor<2, dim>
      get_DD() const override;

      virtual Tensor<3, dim>
      get_PP() const override;

      virtual SymmetricTensor<4, dim>
      get_HH() const override;

      virtual void
      update_end_of_timestep() override;

    private:
      SymmetricTensor<2, dim> Q_t;
      SymmetricTensor<2, dim> Q_t1;

      double                  psi;
      Tensor<1, dim>          B;
      SymmetricTensor<2, dim> S;
      SymmetricTensor<2, dim> BB;
      Tensor<3, dim>          PP;
      SymmetricTensor<4, dim> HH;

      // 一个用于存储所有中间计算的数据结构。我们很快就会准确地看到如何利用这一点来使我们实际进行计算的那部分代码变得干净和容易（好吧，至少是更容易）遵循和维护。但是现在，我们可以说，它将允许我们把计算中间量的导数的那部分代码从使用它们的地方移开。

      mutable GeneralDataStorage cache;

      // 接下来的两个函数是用来更新场和内部变量的状态的，在我们进行任何详细的计算之前会被调用。

      void
      set_primary_variables(const SymmetricTensor<2, dim> &C,
                            const Tensor<1, dim> &         H) const;

      void
      update_internal_variable(const DiscreteTime &time);

      // 该类接口的其余部分专门用于计算自由能密度函数及其所有导数所需的组件的方法。

      // 运动学变量，或称场变量。

      const Tensor<1, dim> &
      get_H() const;

      const SymmetricTensor<2, dim> &
      get_C() const;

      // 饱和度函数的一般化表述，所需的构成参数作为参数传递给每个函数。

      double
      get_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      double
      get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      double
      get_f_mu(const double mu,
               const double mu_inf,
               const double mu_h_sat) const;

      // 饱和度函数一阶导数的一般化表述，所需的构成参数作为参数传递给每个函数。

      double
      get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      Tensor<1, dim>
      get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const;

      Tensor<1, dim>
      get_df_mu_dH(const double mu,
                   const double mu_inf,
                   const double mu_h_sat) const;

      // 饱和度函数二阶导数的广义公式，所需的构成参数作为参数传递给每个函数。

      double
      get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;

      SymmetricTensor<2, dim>
      get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const;

      SymmetricTensor<2, dim>
      get_d2f_mu_dH_dH(const double mu,
                       const double mu_inf,
                       const double mu_h_sat) const;

      // 从场/运动学变量中直接获得的中间量。

      const double &
      get_det_F() const;

      const SymmetricTensor<2, dim> &
      get_C_inv() const;

      const double &
      get_log_det_F() const;



      const Tensor<1, dim> &
      get_C_inv_dot_H() const;



      // 中间量的一阶导数。

      const SymmetricTensor<4, dim> &
      get_dC_inv_dC() const;

      const SymmetricTensor<2, dim> &
      get_d_tr_C_dC() const;

      const SymmetricTensor<2, dim> &
      get_ddet_F_dC() const;

      const SymmetricTensor<2, dim> &
      get_dlog_det_F_dC() const;



      // 内部变量相对于场变量的导数。注意，我们只需要内部变量的这个导数，因为这个变量只是作为动力学变量线性化的一部分而被微分。

      const SymmetricTensor<4, dim> &
      get_dQ_t_dC(const DiscreteTime &time) const;

      // 中间量的二阶导数。

      const SymmetricTensor<4, dim> &
      get_d2log_det_F_dC_dC() const;

      const SymmetricTensor<4, dim> &
      get_d2det_F_dC_dC() const;



      const Tensor<3, dim> &
      get_d2H_dot_C_inv_dot_H_dC_dH() const;

      const SymmetricTensor<4, dim> &
      get_d2H_dot_C_inv_dot_H_dC_dC() const;
    };

    template <int dim>
    Magnetoviscoelastic_Constitutive_Law<
      dim>::Magnetoviscoelastic_Constitutive_Law(const ConstitutiveParameters
                                                   &constitutive_parameters)
      : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
          constitutive_parameters)
      , Q_t(Physics::Elasticity::StandardTensors<dim>::I)
      , Q_t1(Physics::Elasticity::StandardTensors<dim>::I)
      , psi(0.0)
    {}

    template <int dim>
    void
    Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_data(
      const SymmetricTensor<2, dim> &C,
      const Tensor<1, dim> &         H,
      const DiscreteTime &           time)
    {
      // 记录应用的变形状态以及磁载荷。此后，根据新的变形状态更新内部（粘性）变量。

      set_primary_variables(C, H);
      update_internal_variable(time);

      // 根据当前磁场获取弹性和粘性饱和函数的值...


                                     this->get_mu_e_inf(), 
                                     this->get_mu_e_h_sat());

                                     const double f_mu_v =
                                       get_f_mu(this->get_mu_v(),
                                                this->get_mu_v_inf(),
                                                this->get_mu_v_h_sat());

                                     // ... 以及它们的一阶导数...

                                     const Tensor<1, dim> df_mu_e_dH =
                                       get_df_mu_dH(this->get_mu_e(),
                                                    this->get_mu_e_inf(),
                                                    this->get_mu_e_h_sat());

                                     const Tensor<1, dim> df_mu_v_dH =
                                       get_df_mu_dH(this->get_mu_v(),
                                                    this->get_mu_v_inf(),
                                                    this->get_mu_v_h_sat());

                                     // ...以及它们的二阶导数。

                                     const SymmetricTensor<2, dim>
                                       d2f_mu_e_dH_dH = get_d2f_mu_dH_dH(
                                         this->get_mu_e(),
                                         this->get_mu_e_inf(),
                                         this->get_mu_e_h_sat());

                                     const SymmetricTensor<2, dim>
                                       d2f_mu_v_dH_dH = get_d2f_mu_dH_dH(
                                         this->get_mu_v(),
                                         this->get_mu_v_inf(),
                                         this->get_mu_v_h_sat());

                                     // 中间量。请注意，由于我们是从一个缓存中获取这些值，而这个缓存的寿命比这个函数调用的寿命长，所以我们可以对结果进行别名，而不是从缓存中复制这个值。

                                     const double &det_F = get_det_F();
                                     const SymmetricTensor<2, dim> &C_inv =
                                       get_C_inv();

                                     const double &log_det_F = get_log_det_F();
                                     const double &tr_C      = get_trace_C();
                                     const double &H_dot_C_inv_dot_H =
                                       get_H_dot_C_inv_dot_H();

                                     // 中间值的第一导数，以及内部变量相对于右Cauchy-Green变形张量的那部分。

                                     const SymmetricTensor<2, dim> &d_tr_C_dC =
                                       get_d_tr_C_dC();
                                     const SymmetricTensor<2, dim> &ddet_F_dC =
                                       get_ddet_F_dC();
                                     const SymmetricTensor<2, dim>
                                       &dlog_det_F_dC = get_dlog_det_F_dC();

                                     const SymmetricTensor<4, dim> &dQ_t_dC =
                                       get_dQ_t_dC(time);

                                     const Tensor<1, dim>
                                       &dH_dot_C_inv_dot_H_dH =
                                         get_dH_dot_C_inv_dot_H_dH();

                                     const SymmetricTensor<2, dim>
                                       &dH_dot_C_inv_dot_H_dC =
                                         get_dH_dot_C_inv_dot_H_dC();

                                     // 中间值的二阶导数。

                                     const SymmetricTensor<4, dim>
                                       &d2log_det_F_dC_dC =
                                         get_d2log_det_F_dC_dC();

                                     const SymmetricTensor<4, dim>
                                       &d2det_F_dC_dC = get_d2det_F_dC_dC();

                                     const SymmetricTensor<2, dim>
                                       &d2H_dot_C_inv_dot_H_dH_dH =
                                         get_d2H_dot_C_inv_dot_H_dH_dH();

                                     const Tensor<3, dim>
                                       &d2H_dot_C_inv_dot_H_dC_dH =
                                         get_d2H_dot_C_inv_dot_H_dC_dH();

                                     const SymmetricTensor<4, dim>
                                       &d2H_dot_C_inv_dot_H_dC_dC =
                                         get_d2H_dot_C_inv_dot_H_dC_dC();

                                     // 由于线性化的定义变得特别冗长，我们将把自由能密度函数分解成三个相加的部分。

                                     // --类似 "新胡克 "的项。

                                     // -- 与速度有关的项，以及

                                     // --类似于储存在磁场中的能量的项。

                                     // 为了保持一致，这些贡献中的每一个都将被单独加入到我们想要计算的变量中，其顺序也是如此。

                                     // 所以，首先这是能量密度函数本身。

                                     psi =
                                       (0.5 * this->get_mu_e() * f_mu_e) *
                                         (tr_C - dim - 2.0 * std::log(det_F)) +
                                       this->get_lambda_e() *
                                         (std::log(det_F) * std::log(det_F));
                                     psi +=
                                       (0.5 * this->get_mu_v() * f_mu_v) *
                                       (Q_t *
                                          (std::pow(det_F, -2.0 / dim) * C) -
                                        dim - std::log(determinant(Q_t)));
                                     psi -= (0.5 * this->get_mu_0() *
                                             this->get_mu_r()) *
                                            det_F * (H * C_inv * H);

                                     // ...然后是磁感应强度和Piola-Kirchhoff应力。

                                     B = -(0.5 * this->get_mu_e() *
                                           (tr_C - dim - 2.0 * log_det_F)) *
                                         df_mu_e_dH;
                                     B -= (0.5 * this->get_mu_v()) *
                                          (Q_t *
                                             (std::pow(det_F, -2.0 / dim) * C) -
                                           dim - std::log(determinant(Q_t))) *
                                          df_mu_v_dH;
                                     B += 0.5 * this->get_mu_0() *
                                          this->get_mu_r() * det_F *
                                          dH_dot_C_inv_dot_H_dH;

                                     S =
                                       2.0 *
                                         (0.5 * this->get_mu_e() * f_mu_e) * //
                                         (d_tr_C_dC - 2.0 * dlog_det_F_dC)   //
                                       + 2.0 * this->get_lambda_e() *
                                           (2.0 * log_det_F * dlog_det_F_dC); //
                                     S +=
                                       2.0 * (0.5 * this->get_mu_v() * f_mu_v) *
                                       ((Q_t * C) *
                                          ((-2.0 / dim) *
                                           std::pow(det_F, -2.0 / dim - 1.0) *
                                           ddet_F_dC) +
                                        std::pow(det_F, -2.0 / dim) *
                                          Q_t); // dC/dC = II
                                     S -= 2.0 *
                                          (0.5 * this->get_mu_0() *
                                           this->get_mu_r()) *               //
                                          (H_dot_C_inv_dot_H * ddet_F_dC     //
                                           + det_F * dH_dot_C_inv_dot_H_dC); //

                                     // ......
                                     // 最后是由于动能变量的线性化而产生的切线。

                                     BB = -(0.5 * this->get_mu_e() *
                                            (tr_C - dim - 2.0 * log_det_F)) *
                                          d2f_mu_e_dH_dH;
                                     BB -= (0.5 * this->get_mu_v()) *
                                           (Q_t * (std::pow(det_F, -2.0 / dim) *
                                                   C) -
                                            dim - std::log(determinant(Q_t))) *
                                           d2f_mu_v_dH_dH;
                                     BB += 0.5 * this->get_mu_0() *
                                           this->get_mu_r() * det_F *
                                           d2H_dot_C_inv_dot_H_dH_dH;

                                     PP = -2.0 * (0.5 * this->get_mu_e()) *
                                          outer_product(Tensor<2, dim>(
                                                          d_tr_C_dC -
                                                          2.0 * dlog_det_F_dC),
                                                        df_mu_e_dH);
                                     PP -=
                                       2.0 * (0.5 * this->get_mu_v()) *
                                       outer_product(
                                         Tensor<2, dim>(
                                           (Q_t * C) *
                                             ((-2.0 / dim) *
                                              std::pow(det_F,
                                                       -2.0 / dim - 1.0) *
                                              ddet_F_dC) +
                                           std::pow(det_F, -2.0 / dim) * Q_t),
                                         df_mu_v_dH);
                                     PP +=
                                       2.0 *
                                       (0.5 * this->get_mu_0() *
                                        this->get_mu_r()) *
                                       (outer_product(Tensor<2, dim>(ddet_F_dC),
                                                      dH_dot_C_inv_dot_H_dH) +
                                        det_F * d2H_dot_C_inv_dot_H_dC_dH);

                                     HH =
                                       4.0 * (0.5 * this->get_mu_e() * f_mu_e) *
                                         (-2.0 * d2log_det_F_dC_dC) //
                                       +
                                       4.0 * this->get_lambda_e() * //
                                         (2.0 * outer_product(dlog_det_F_dC,
                                                              dlog_det_F_dC) //
                                          + 2.0 * log_det_F *
                                              d2log_det_F_dC_dC); //
                                     HH +=
                                       4.0 * (0.5 * this->get_mu_v() * f_mu_v) *
                                       (outer_product(
                                          (-2.0 / dim) *
                                            std::pow(det_F, -2.0 / dim - 1.0) *
                                            ddet_F_dC,
                                          C * dQ_t_dC + Q_t) +
                                        (Q_t * C) *
                                          (outer_product(
                                             ddet_F_dC,
                                             (-2.0 / dim) * (-2.0 / dim - 1.0) *
                                               std::pow(det_F,
                                                        -2.0 / dim - 2.0) *
                                               ddet_F_dC) +
                                           ((-2.0 / dim) *
                                            std::pow(det_F, -2.0 / dim - 1.0) *
                                            d2det_F_dC_dC)) +
                                        outer_product(
                                          Q_t,
                                          (-2.0 / dim) *
                                            std::pow(det_F, -2.0 / dim - 1.0) *
                                            ddet_F_dC) +
                                        std::pow(det_F, -2.0 / dim) * dQ_t_dC);
                                     HH -=
                                       4.0 *
                                       (0.5 * this->get_mu_0() *
                                        this->get_mu_r()) *               //
                                       (H_dot_C_inv_dot_H * d2det_F_dC_dC //
                                        +
                                        outer_product(ddet_F_dC,
                                                      dH_dot_C_inv_dot_H_dC) //
                                        + outer_product(dH_dot_C_inv_dot_H_dC,
                                                        ddet_F_dC)            //
                                        + det_F * d2H_dot_C_inv_dot_H_dC_dC); //

                                     // 现在我们已经用完了存储在缓存中的所有临时变量，我们可以把它清除掉，以释放一些内存。

                                     cache.reset();
    }

    template <int dim>
    double
    Magnetoviscoelastic_Constitutive_Law<dim>::get_psi() const
    {
      return psi;
    }

    template <int dim>
    Tensor<1, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_B() const
    {
      return B;
    }

    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_S() const
    {
      return S;
    }

    template <int dim>
    SymmetricTensor<2, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_DD() const
    {
      return BB;
    }

    template <int dim>
    Tensor<3, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_PP() const
    {
      return PP;
    }

    template <int dim>
    SymmetricTensor<4, dim>
    Magnetoviscoelastic_Constitutive_Law<dim>::get_HH() const
    {
      return HH;
    }

    template <int dim>
    void
    Magnetoviscoelastic_Constitutive_Law<dim>::update_end_of_timestep()
    {
      Q_t1 = Q_t;
    }

    template <int dim>
    void
    Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_variable(
      const DiscreteTime &time)
    {
      const double delta_t = this->get_delta_t(time);

      Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) *
            (Q_t1 + (delta_t / this->get_tau_v()) *
                      std::pow(get_det_F(), 2.0 / dim) * get_C_inv());
    }

    // 接下来的几个函数实现了饱和度函数的广义表述，以及它的各种导数。

    template <int dim>
    double
    Magnetoviscoelastic_Constitutive_Law<dim>::get_two_h_dot_h_div_h_sat_squ(
      const double mu_h_sat) const
    {
      const Tensor<1, dim> &H = get_H();
      return (2.0 * H * H) / (mu_h_sat * mu_h_sat);
    }

    template <int dim>
    double
    Magnetoviscoelastic_Constitutive_Law<
      dim>::get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
    {
      return std::tanh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat));
    }

    // 一个比例函数，它将使剪切模量在磁场的影响下发生变化（增加）。


    double
    Magnetoviscoelastic_Constitutive_Law<dim>::get_f_mu(
      const double mu,
      const double mu_inf,
      const double mu_h_sat) const
    {
      return 1.0 +
             (mu_inf / mu - 1.0) * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat);
    }

    // 缩放函数的一阶导数

    template <int dim>
    double
    Magnetoviscoelastic_Constitutive_Law<
      dim>::get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
    {
      return std::pow(1.0 / std::cosh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat)),
                      2.0);
    }

    template <int dim>
    Tensor<1, dim>
    Magnetoviscoelastic_Constitutive_Law<
      dim>::get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const
    {
      return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) * get_H();
    }

    template <int dim>
      Tensor<1, dim>
      Magnetoviscoelastic_Constitutive_Law<dim>::get_df_mu_dH(
        const double mu,
        const double mu_inf,
        const double mu_h_sat) const template <int dim>
      return (mu_inf / mu - 1.0) *
      (get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
       get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat));
  }

  template <int dim>
  double
  Magnetoviscoelastic_Constitutive_Law<
    dim>::get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
  {
    return -2.0 * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
           get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat);
  }

  template <int dim>
  SymmetricTensor<2, dim>
  Magnetoviscoelastic_Constitutive_Law<
    dim>::get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const
  {
    return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) *
           Physics::Elasticity::StandardTensors<dim>::I;
  }

  template <int dim>
  SymmetricTensor<2, dim>
  Magnetoviscoelastic_Constitutive_Law<dim>::get_d2f_mu_dH_dH(
    const double mu,
    const double mu_inf,
    const double mu_h_sat) const
  {
    return (mu_inf / mu - 1.0) *
           (get_d2tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
              symmetrize(
                outer_product(get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat),
                              get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat))) +
            get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
              get_d2two_h_dot_h_div_h_sat_squ_dH_dH(mu_h_sat));
  }

  // 对于我们为这个材料类采用的缓存计算方法，所有计算的根基是场变量，以及不可改变的辅助数据，如构成参数和时间步长。因此，我们需要以与其他变量不同的方式将它们输入缓存，因为它们是由类本身之外规定的输入。这个函数只是将它们从输入参数中直接添加到缓存中，同时检查那里是否有等效的数据（我们希望每个时间步长或牛顿迭代只调用一次`update_internal_data()`方法）。

  template <int dim>
  void
  Magnetoviscoelastic_Constitutive_Law<dim>::set_primary_variables(
    const SymmetricTensor<2, dim> &C,
    const Tensor<1, dim> &         H) const
  {
    // 设置  $\boldsymbol{\mathbb{H}}$  的值。

    const std::string name_H("H");
    Assert(!cache.stores_object_with_name(name_H),
           ExcMessage(
             "The primary variable has already been added to the cache."));
    cache.add_unique_copy(name_H, H);

    // 设置  $\mathbf{C}$  的值。

    const std::string name_C("C");
    Assert(!cache.stores_object_with_name(name_C),
           ExcMessage(
             "The primary variable has already been added to the cache."));
    cache.add_unique_copy(name_C, C);
  }

  // 此后，我们可以在任何时间点从缓存中获取它们。

  template <int dim>
  const Tensor<1, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_H() const
  {
    const std::string name("H");
    Assert(cache.stores_object_with_name(name),
           ExcMessage("Primary variables must be added to the cache."));
    return cache.template get_object_with_name<Tensor<1, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<2, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_C() const
  {
    const std::string name("C");
    Assert(cache.stores_object_with_name(name),
           ExcMessage("Primary variables must be added to the cache."));
    return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
  }

  // 当我们需要主要变量时，保证它们在缓存中，我们不能从它们中计算出所有的中间值（无论是直接，还是间接）。

  // 如果缓存中还没有存储我们要找的值，那么我们就快速计算，把它存储在缓存中，然后返回刚刚存储在缓存中的值。这样我们就可以把它作为一个引用返回，避免复制对象。同样的道理也适用于复合函数可能依赖的任何值。换句话说，如果在我们目前感兴趣的计算之前有一个依赖链，那么在我们继续使用这些值之前，我们可以保证解决这些依赖关系。尽管从缓存中获取数据是有成本的，但
  // "已解决的依赖关系
  // "的概念可能足够方便，使其值得看一下这个额外的成本。如果这些材料定律被嵌入到有限元框架中，那么额外的成本甚至可能不会被注意到。

  template <int dim>
  const double &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_det_F() const
  {
    const std::string name("det_F");
    if (cache.stores_object_with_name(name) == false)
      {
        const double det_F = std::sqrt(determinant(get_C()));
        AssertThrow(det_F > 0.0,
                    ExcMessage("Volumetric Jacobian must be positive."));
        cache.add_unique_copy(name, det_F);
      }

    return cache.template get_object_with_name<double>(name);
  }

  template <int dim>
  const SymmetricTensor<2, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv() const
  {
    const std::string name("C_inv");
    if (cache.stores_object_with_name(name) == false)
      {
        cache.add_unique_copy(name, invert(get_C()));
      }

    return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
  }

  template <int dim>
  const double &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_log_det_F() const
  {
    const std::string name("log(det_F)");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, std::log(get_det_F()));

    return cache.template get_object_with_name<double>(name);
  }

  template <int dim>
  const double &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_trace_C() const
  {
    const std::string name("trace(C)");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, trace(get_C()));

    return cache.template get_object_with_name<double>(name);
  }

  template <int dim>
  const Tensor<1, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv_dot_H() const
  {
    const std::string name("C_inv_dot_H");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, get_C_inv() * get_H());

    return cache.template get_object_with_name<Tensor<1, dim>>(name);
  }

  template <int dim>
  const double &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_H_dot_C_inv_dot_H() const
  {
    const std::string name("H_dot_C_inv_dot_H");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, get_H() * get_C_inv_dot_H());

    return cache.template get_object_with_name<double>(name);
  }

  template <int dim>
  const SymmetricTensor<4, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_dQ_t_dC(
    const DiscreteTime &time) const
  {
    const std::string name("dQ_t_dC");
    if (cache.stores_object_with_name(name) == false)
      {
        const double  delta_t = this->get_delta_t(time);
        const double &det_F   = get_det_F();

        const SymmetricTensor<4, dim> dQ_t_dC =
          (1.0 / (1.0 + delta_t / this->get_tau_v())) *
          (delta_t / this->get_tau_v()) *
          ((2.0 / dim) * std::pow(det_F, 2.0 / dim - 1.0) *
             outer_product(get_C_inv(), get_ddet_F_dC()) +
           std::pow(det_F, 2.0 / dim) * get_dC_inv_dC());

        cache.add_unique_copy(name, dQ_t_dC);
      }

    return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<4, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_dC_inv_dC() const
  {
    const std::string name("dC_inv_dC");
    if (cache.stores_object_with_name(name) == false)
      {
        const SymmetricTensor<2, dim> &C_inv = get_C_inv();
        SymmetricTensor<4, dim>        dC_inv_dC;

        for (unsigned int A = 0; A < dim; ++A)
          for (unsigned int B = A; B < dim; ++B)
            for (unsigned int C = 0; C < dim; ++C)
              for (unsigned int D = C; D < dim; ++D)
                dC_inv_dC[A][B][C][D] -=               //
                  0.5 * (C_inv[A][C] * C_inv[B][D]     //
                         + C_inv[A][D] * C_inv[B][C]); //

        cache.add_unique_copy(name, dC_inv_dC);
      }

    return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<2, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_d_tr_C_dC() const
  {
    const std::string name("d_tr_C_dC");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, Physics::Elasticity::StandardTensors<dim>::I);

    return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<2, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_ddet_F_dC() const
  {
    const std::string name("ddet_F_dC");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, 0.5 * get_det_F() * get_C_inv());

    return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<2, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_dlog_det_F_dC() const
  {
    const std::string name("dlog_det_F_dC");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, 0.5 * get_C_inv());

    return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
  }

  template <int dim>
  const Tensor<1, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dH() const
  {
    const std::string name("dH_dot_C_inv_dot_H_dH");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, 2.0 * get_C_inv_dot_H());

    return cache.template get_object_with_name<Tensor<1, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<2, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dC() const
  {
    const std::string name("dH_dot_C_inv_dot_H_dC");
    if (cache.stores_object_with_name(name) == false)
      {
        const Tensor<1, dim> C_inv_dot_H = get_C_inv_dot_H();
        cache.add_unique_copy(
          name, -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H)));
      }

    return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<4, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_d2log_det_F_dC_dC() const
  {
    const std::string name("d2log_det_F_dC_dC");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, 0.5 * get_dC_inv_dC());

    return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<4, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_d2det_F_dC_dC() const
  {
    const std::string name("d2det_F_dC_dC");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name,
                            0.5 * (outer_product(get_C_inv(), get_ddet_F_dC()) +
                                   get_det_F() * get_dC_inv_dC()));

    return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<2, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dH_dH()
    const
  {
    const std::string name("d2H_dot_C_inv_dot_H_dH_dH");
    if (cache.stores_object_with_name(name) == false)
      cache.add_unique_copy(name, 2.0 * get_C_inv());

    return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
  }

  template <int dim>
  const Tensor<3, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dH()
    const
  {
    const std::string name("d2H_dot_C_inv_dot_H_dC_dH");
    if (cache.stores_object_with_name(name) == false)
      {
        const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H();
        const SymmetricTensor<2, dim> &C_inv       = get_C_inv();

        Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH;
        for (unsigned int A = 0; A < dim; ++A)
          for (unsigned int B = 0; B < dim; ++B)
            for (unsigned int C = 0; C < dim; ++C)
              d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -=
                C_inv[A][C] * C_inv_dot_H[B] + //
                C_inv_dot_H[A] * C_inv[B][C];  //

        cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dH);
      }

    return cache.template get_object_with_name<Tensor<3, dim>>(name);
  }

  template <int dim>
  const SymmetricTensor<4, dim> &
  Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dC()
    const
  {
    const std::string name("d2H_dot_C_inv_dot_H_dC_dC");
    if (cache.stores_object_with_name(name) == false)
      {
        const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H();
        const SymmetricTensor<2, dim> &C_inv       = get_C_inv();

        SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC;
        for (unsigned int A = 0; A < dim; ++A)
          for (unsigned int B = A; B < dim; ++B)
            for (unsigned int C = 0; C < dim; ++C)
              for (unsigned int D = C; D < dim; ++D)
                d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] +=
                  0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] +
                         C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] +
                         C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] +
                         C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]);

        cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dC);
      }

    return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
  }
  // @sect4{Rheological experiment parameters}

  //  @p RheologicalExperimentParameters 类是用来驱动数值实验的，这些实验将在我们已经实现了构成法则的耦合材料上进行。

  class RheologicalExperimentParameters : public ParameterAcceptor
  {
  public:
    RheologicalExperimentParameters();

    // 这些是要模拟的流变学试样的尺寸。它们有效地定义了我们虚拟实验的测量点。

    double sample_radius = 0.01;
    double sample_height = 0.001;

    // 三个稳态负载参数分别是

    // - 轴向拉伸。

    // -- 剪切应变振幅，和

    // - 轴向磁场强度。

    double lambda_2 = 0.95;
    double gamma_12 = 0.05;
    double H_2      = 60.0e3;

    // 此外，随时间变化的流变学负载条件的参数为

    // --加载周期的频率。

    // - 负载周期的数量，以及

    // - 每个周期的离散时间步数。

    double       frequency         = 1.0 / (2.0 * numbers::PI);
    unsigned int n_cycles          = 5;
    unsigned int n_steps_per_cycle = 2500;

    // 我们还声明了一些不言自明的参数，这些参数与用速率依赖型和速率非依赖型材料进行的实验所产生的输出数据有关。

    bool        output_data_to_file = true;
    std::string output_filename_rd  = "experimental_results-rate_dependent.csv";
    std::string output_filename_ri =
      "experimental_results-rate_independent.csv";

    // 接下来的几个函数将计算与时间有关的实验参数...

    double
    start_time() const;

    double
    end_time() const;

    double
    delta_t() const;

    // ...... 而下面两个则规定了任何时候的机械和磁力负载......

    Tensor<1, 3>
    get_H(const double time) const;

    Tensor<2, 3>
    get_F(const double time) const;

    // ...... 而这最后一个是将实验的状态输出到控制台。

    bool
    print_status(const int step_number) const;

    bool initialized = false;
  };

  RheologicalExperimentParameters::RheologicalExperimentParameters()
    : ParameterAcceptor("/Coupled Constitutive Laws/Rheological Experiment/")
  {
    add_parameter("Experimental sample radius", sample_radius);
    add_parameter("Experimental sample radius", sample_height);

    add_parameter("Axial stretch", lambda_2);
    add_parameter("Shear strain amplitude", gamma_12);
    add_parameter("Axial magnetic field strength", H_2);

    add_parameter("Frequency", frequency);
    add_parameter("Number of loading cycles", n_cycles);
    add_parameter("Discretisation for each cycle", n_steps_per_cycle);

    add_parameter("Output experimental results to file", output_data_to_file);
    add_parameter("Output file name (rate dependent constitutive law)",
                  output_filename_rd);
    add_parameter("Output file name (rate independent constitutive law)",
                  output_filename_ri);

    parse_parameters_call_back.connect([&]() -> void { initialized = true; });
  }

  double
  RheologicalExperimentParameters::start_time() const
  {
    return 0.0;
  }

  double
  RheologicalExperimentParameters::end_time() const
  {
    return n_cycles / frequency;
  }

  double
  RheologicalExperimentParameters::delta_t() const
  {
    return (end_time() - start_time()) / (n_steps_per_cycle * n_cycles);
  }

  bool
  RheologicalExperimentParameters::print_status(const int step_number) const
  {
    return (step_number % (n_cycles * n_steps_per_cycle / 100)) == 0;
  }

  // 施加的磁场总是与流变仪转子的旋转轴对齐。

  Tensor<1, 3>
  RheologicalExperimentParameters::get_H(const double) const
  {
    return Tensor<1, 3>({0.0, 0.0, H_2});
  }

  // 根据流变仪和样品的几何形状、采样点和实验参数，计算出应用的变形（梯度）。根据介绍中记录的位移曲线，变形梯度可以用直角坐标表示为
  // @f[
  //  \mathbf{F} = \begin{bmatrix}
  //     \frac{\cos\left(\alpha\right)}{\sqrt{\lambda_{3}}}
  //  & -\frac{\sin\left(\alpha\right)}{\sqrt{\lambda_{3}}}
  //  & -\tau R \sqrt{\lambda_{3}} \sin\left(\Theta + \alpha\right)
  //  \\  \frac{\sin\left(\alpha\right)}{\sqrt{\lambda_{3}}}
  //  & \frac{\cos\left(\alpha\right)}{\sqrt{\lambda_{3}}}
  //  & -\tau R \sqrt{\lambda_{3}} \cos\left(\Theta + \alpha\right)
  //  \\  0 & 0 & \lambda_{3}
  //  \end{bmatrix}
  //  @f] 。

  Tensor<2, 3>
  RheologicalExperimentParameters::get_F(const double time) const
  {
    AssertThrow((sample_radius > 0.0 && sample_height > 0.0),
                ExcMessage("Non-physical sample dimensions"));
    AssertThrow(lambda_2 > 0.0,
                ExcMessage("Non-physical applied axial stretch"));

    const double sqrt_lambda_2     = std::sqrt(lambda_2);
    const double inv_sqrt_lambda_2 = 1.0 / sqrt_lambda_2;

    const double alpha_max =
      std::atan(std::tan(gamma_12) * sample_height /
                sample_radius); // Small strain approximation
    const double A       = sample_radius * alpha_max;
    const double w       = 2.0 * numbers::PI * frequency; // in rad /s
    const double gamma_t = A * std::sin(w * time);
    const double tau_t =
      gamma_t /
      (sample_radius * sample_height); // Torsion angle per unit length
    const double alpha_t = tau_t * lambda_2 * sample_height;

    Tensor<2, 3> F;
    F[0][0] = inv_sqrt_lambda_2 * std::cos(alpha_t);
    F[0][1] = -inv_sqrt_lambda_2 * std::sin(alpha_t);
    F[0][2] = -tau_t * sample_radius * sqrt_lambda_2 * std::sin(alpha_t);
    F[1][0] = inv_sqrt_lambda_2 * std::sin(alpha_t);
    F[1][1] = inv_sqrt_lambda_2 * std::cos(alpha_t);
    F[1][2] = tau_t * sample_radius * sqrt_lambda_2 * std::cos(alpha_t);
    F[2][0] = 0.0;
    F[2][1] = 0.0;
    F[2][2] = lambda_2;

    AssertThrow((F[0][0] > 0) && (F[1][1] > 0) && (F[2][2] > 0),
                ExcMessage("Non-physical deformation gradient component."));
    AssertThrow(std::abs(determinant(F) - 1.0) < 1e-6,
                ExcMessage("Volumetric Jacobian is not equal to unity."));

    return F;
  }
  // @sect4{Rheological experiment: Parallel plate rotational rheometer}

  // 这是将驱动数值实验的函数。

  template <int dim>
  void
  run_rheological_experiment(
    const RheologicalExperimentParameters &experimental_parameters,
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
      &material_hand_calculated,
    Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
      &               material_assisted_computation,
    TimerOutput &     timer,
    const std::string filename)
  {
    // 我们可以利用手工实现的构成法，将我们用它达到的结果与用AD或SD得到的结果进行比较。通过这种方式，我们可以验证它们产生了相同的结果（这表明要么两种实现方式都有很大的可能性是正确的，要么就是它们都有相同的缺陷而不正确）。无论哪种方式，对于完全自我实现的变体来说，这都是一个很好的理智检查，当发现结果之间的差异时，当然可以作为一种调试策略）。)

    const auto check_material_class_results =
      [](const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &to_verify,
         const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &blessed,
         const double tol = 1e-6) {
        (void)to_verify;
        (void)blessed;
        (void)tol;

        Assert(std::abs(blessed.get_psi() - to_verify.get_psi()) < tol,
               ExcMessage("No match for psi. Error: " +
                          Utilities::to_string(std::abs(blessed.get_psi() -
                                                        to_verify.get_psi()))));

        Assert((blessed.get_B() - to_verify.get_B()).norm() < tol,
               ExcMessage("No match for B. Error: " +
                          Utilities::to_string(
                            (blessed.get_B() - to_verify.get_B()).norm())));
        Assert((blessed.get_S() - to_verify.get_S()).norm() < tol,
               ExcMessage("No match for S. Error: " +
                          Utilities::to_string(
                            (blessed.get_S() - to_verify.get_S()).norm())));

        Assert((blessed.get_DD() - to_verify.get_DD()).norm() < tol,
               ExcMessage("No match for BB. Error: " +
                          Utilities::to_string(
                            (blessed.get_DD() - to_verify.get_DD()).norm())));
        Assert((blessed.get_PP() - to_verify.get_PP()).norm() < tol,
               ExcMessage("No match for PP. Error: " +
                          Utilities::to_string(
                            (blessed.get_PP() - to_verify.get_PP()).norm())));
        Assert((blessed.get_HH() - to_verify.get_HH()).norm() < tol,
               ExcMessage("No match for HH. Error: " +
                          Utilities::to_string(
                            (blessed.get_HH() - to_verify.get_HH()).norm())));
      };

    // 我们将把材料的构成性响应输出到文件中进行后处理，所以在这里我们声明一个`stream`，它将作为这个输出的缓冲区。我们将使用一个简单的CSV格式来输出结果。

    std::ostringstream stream;
    stream
      << "Time;Axial magnetic field strength [A/m];Axial magnetic induction [T];Shear strain [%];Shear stress [Pa]\n";

    // 使用DiscreteTime类，我们使用一个固定的时间步长来迭代每个时间段。

    for (DiscreteTime time(experimental_parameters.start_time(),
                           experimental_parameters.end_time() +
                             experimental_parameters.delta_t(),
                           experimental_parameters.delta_t());
         time.is_at_end() == false;
         time.advance_time())
      {
        if (experimental_parameters.print_status(time.get_step_number()))
          std::cout << "Timestep = " << time.get_step_number()
                    << " @ time = " << time.get_current_time() << "s."
                    << std::endl;

        // 我们获取并计算在这个时间步长中应用于材料的负载...

        const Tensor<1, dim> H =
          experimental_parameters.get_H(time.get_current_time());
        const Tensor<2, dim> F =
          experimental_parameters.get_F(time.get_current_time());
        const SymmetricTensor<2, dim> C = Physics::Elasticity::Kinematics::C(F);

        // ...然后我们更新材料的状态...

        {
          TimerOutput::Scope timer_section(timer, "Hand calculated");
          material_hand_calculated.update_internal_data(C, H, time);
          material_hand_calculated.update_end_of_timestep();
        }

        {
          TimerOutput::Scope timer_section(timer, "Assisted computation");
          material_assisted_computation.update_internal_data(C, H, time);
          material_assisted_computation.update_end_of_timestep();
        }

        // ...并测试两者之间的差异。

        check_material_class_results(material_hand_calculated,
                                     material_assisted_computation);

        if (experimental_parameters.output_data_to_file)
          {
            // 接下来我们要做的是收集一些结果进行后处理。所有的数量都在
            // "当前配置 "中（而不是
            // "参考配置"，所有由构成法则计算的数量都在这个框架中）。

            const Tensor<1, dim> h =
              Physics::Transformations::Covariant::push_forward(H, F);
            const Tensor<1, dim> b =
              Physics::Transformations::Piola::push_forward(
                material_hand_calculated.get_B(), F);
            const SymmetricTensor<2, dim> sigma =
              Physics::Transformations::Piola::push_forward(
                material_hand_calculated.get_S(), F);
            stream << time.get_current_time() << ";" << h[2] << ";" << b[2]
                   << ";" << F[1][2] * 100.0 << ";" << sigma[1][2] << "\n";
          }
      }

    // 最后，我们将应变应力和磁载荷历史输出到文件中。

    if (experimental_parameters.output_data_to_file)
      {
        std::ofstream output(filename);
        output << stream.str();
      }
  }
  // @sect4{The CoupledConstitutiveLaws::run() function}

  // 这个驱动函数的目的是读取文件中的所有参数，并在此基础上创建每个构成法则的代表性实例，并调用函数对其进行流变学实验。

  void
  run(int argc, char *argv[])
  {
    using namespace dealii;

    constexpr unsigned int dim = 3;

    const ConstitutiveParameters          constitutive_parameters;
    const RheologicalExperimentParameters experimental_parameters;

    std::string parameter_file;
    if (argc > 1)
      parameter_file = argv[1];
    else
      parameter_file = "parameters.prm";
    ParameterAcceptor::initialize(parameter_file, "used_parameters.prm");

    // 我们开始实际工作，使用我们与速率无关的构成法配置和运行实验。这里的自动可微调数类型是硬编码的，但是通过一些巧妙的模板设计，可以在运行时选择使用哪种框架（例如，通过参数文件选择）。我们将同时用完全手工实现的反面材料法进行实验，并检查它与我们的辅助实现的计算结果。

    {
      TimerOutput timer(std::cout,
                        TimerOutput::summary,
                        TimerOutput::wall_times);
      std::cout
        << "Coupled magnetoelastic constitutive law using automatic differentiation."
        << std::endl;

      constexpr Differentiation::AD::NumberTypes ADTypeCode =
        Differentiation::AD::NumberTypes::sacado_dfad_dfad;

      Magnetoelastic_Constitutive_Law<dim> material(constitutive_parameters);
      Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode> material_ad(
        constitutive_parameters);

      run_rheological_experiment(experimental_parameters,
                                 material,
                                 material_ad,
                                 timer,
                                 experimental_parameters.output_filename_ri);

      std::cout << "... all calculations are correct!" << std::endl;
    }

    // 接下来我们对与速率相关的构成法则做同样的处理。如果SymEngine被设置为使用LLVM即时编译器，则默认选择最高性能的选项，该编译器（结合一些积极的编译标志）产生所有可用选项中最快的代码评估路径。作为后备措施，所谓的
    // "lambda
    // "优化器（它只需要一个兼容C++11的编译器）将被选中。同时，我们将要求CAS进行普通子表达式的消除，以尽量减少评估过程中使用的中间计算的数量。我们将记录在SD实现的构造器内执行
    // "初始化 "步骤所需的时间，因为这正是上述转换发生的地方。

    {
      TimerOutput timer(std::cout,
                        TimerOutput::summary,
                        TimerOutput::wall_times);
      std::cout
        << "Coupled magneto-viscoelastic constitutive law using symbolic differentiation."
        << std::endl;

#ifdef DEAL_II_SYMENGINE_WITH_LLVM
      std::cout << "Using LLVM optimizer." << std::endl;
      constexpr Differentiation::SD::OptimizerType optimizer_type =
        Differentiation::SD::OptimizerType::llvm;
      constexpr Differentiation::SD::OptimizationFlags optimization_flags =
        Differentiation::SD::OptimizationFlags::optimize_all;
#else
      std::cout << "Using lambda optimizer." << std::endl;
      constexpr Differentiation::SD::OptimizerType optimizer_type =
        Differentiation::SD::OptimizerType::lambda;
      constexpr Differentiation::SD::OptimizationFlags optimization_flags =
        Differentiation::SD::OptimizationFlags::optimize_cse;
#endif

      Magnetoviscoelastic_Constitutive_Law<dim> material(
        constitutive_parameters);

      timer.enter_subsection("Initialize symbolic CL");
      Magnetoviscoelastic_Constitutive_Law_SD<dim> material_sd(
        constitutive_parameters, optimizer_type, optimization_flags);
      timer.leave_subsection();

      run_rheological_experiment(experimental_parameters,
                                 material,
                                 material_sd,
                                 timer,
                                 experimental_parameters.output_filename_rd);

      std::cout << "... all calculations are correct!" << std::endl;
    }
  }

} // namespace CoupledConstitutiveLaws

} // namespace Step71
// @sect3{The main() function}

// 主函数只调用两组要执行的例子的驱动函数。

int
main(int argc, char *argv[])
{
  Step71::SimpleExample::run();
  Step71::CoupledConstitutiveLaws::run(argc, argv);

  return 0;
}
