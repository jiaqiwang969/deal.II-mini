//include/deal.II-translator/lac/solver_0.txt
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

#ifndef dealii_solver_h
#define dealii_solver_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_memory.h>

// Ignore deprecation warnings for auto_ptr.
#include <boost/signals2.hpp>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename number>
class Vector;
#endif

/**
 * 一个用于迭代线性求解器的基类。该类提供了与内存池和确定求解器是否收敛的对象的接口。
 *
 *  <h3>Requirements common to derived solver classes</h3>
 * 一般来说，迭代求解器不依赖于矩阵的任何特殊结构或存储格式。相反，它们只要求矩阵和向量定义某些操作，如矩阵-向量积，或向量之间的标量积。因此，这个类以及实现具体线性求解器的派生类及其成员函数是以矩阵和向量的类型为模板的。然而，矩阵或向量类型必须满足一些共同的要求，才有资格成为这个层次结构中的求解器的可接受类型。这些要求列举如下。
 * 我们下面展示的类不是任何具体的类。相反，它们旨在形成一个具体类必须符合的
 * "签名"。请注意，这个库中的矩阵和向量类当然符合这个接口；因此，SparseMatrix和Vector是这些类的好例子，因为它们提供了成员函数的必要签名（尽管它们也提供了许多解算器事实上不需要的接口
 *
 * - 例如，元素访问）。) 此外，你可能想看看 step-20 、 step-22 或LinearSolvers命名空间中的一些类，以了解如何定义可以作为线性求解器的线性运算符的类似矩阵的类。
 * 具体来说，可以传递给线性求解器的矩阵和向量类需要提供以下接口。
 *
 * @code
 * class Matrix
 * {
 * public:
 *   // Application of matrix to vector src. Write result into dst.
 *   void vmult (VectorType       &dst,
 *               const VectorType &src) const;
 *
 *   // Application of transpose to a vector. This function is,
 *   // however, only used by some iterative methods.
 *   void Tvmult (VectorType       &dst,
 *                const VectorType &src) const;
 * };
 *
 *
 *
 *
 * class Vector
 * {
 * public:
 *   // Define value type of the entries
 *   using value_type = double;
 *
 *   // Resize the current object to have the same size and layout as
 *   // the model_vector argument provided. The second argument
 *   // indicates whether to clear the current object after resizing.
 *   // The second argument must have a default value equal to false.
 *   void reinit (const Vector &model_vector,
 *                const bool  leave_elements_uninitialized = false);
 *
 *   // Inner product between the current object and the argument.
 *   double operator (const Vector &v) const;
 *
 *   // Set all the vector entries to a constant scalar.
 *   Vector & operator = (const double a);
 *
 *   // Deep copy of the vector.
 *   // Important if Vector contains pointers to data to duplicate data.
 *   Vector & operator = (const Vector &x);
 *
 *   // Addition of vectors
 *   void add (const Vector &x);
 *
 *   // Scaled addition of vectors
 *   void add (const double  a,
 *             const Vector &x);
 *
 *   // Scaled addition of vectors
 *   void sadd (const double  a,
 *              const double  b,
 *              const Vector &x);
 *
 *   // Scaled assignment of a vector
 *   void equ (const double  a,
 *             const Vector &x);
 *
 *   // Combined scaled addition of vector x into the current object and
 *   // subsequent inner product of the current object with v.
 *   double add_and_dot (const double  a,
 *                       const Vector &x,
 *                       const Vector &v);
 *
 *   // Multiply the elements of the current object by a fixed value.
 *   Vector & operator= (const double a);
 *
 *   // Return the l2 norm of the vector.
 *   double l2_norm () const;
 * };
 * @endcode
 *
 * 此外，对于某些求解器来说，必须有一个全局函数<tt>swap(VectorType
 * &a, VectorType &b)</tt>来交换这两个向量的值。
 * 最后，求解器还期望有一个GrowingVectorMemory @<VectorType@>.
 * 的实例，这些实例由deal.II库为内置的向量类型提供，但必须为用户提供的向量类明确添加。否则，链接器会抱怨说它找不到发生在
 * @p SolverBase  类中的 GrowingVectorMemory
 * 的构造函数和析构函数。
 *
 *
 * @code
 * // Definition and implementation of vector class
 * class UserVector { ... };
 *
 * // Create explicit instantiation for the vector class. If your project
 * // consists of multiple files, including header files, this instantiation
 * // must be put in a <code>.cc</code> file in order to instantiate only
 * // once.
 * #include <deal.II/lac/vector_memory.templates.h>
 *
 * template class VectorMemory<UserVector>;
 * template class GrowingVectorMemory<UserVector>;
 * @endcode
 *
 * 所用的预处理程序必须具有与矩阵相同的接口，也就是说，特别是它们必须提供一个成员函数
 * @p vmult ，表示预处理程序的应用。
 *
 *  <h3>AdditionalData</h3> 一些求解器需要额外的数据，比如 @p
 * SolverRichardson 类的阻尼参数 @p omega 或 @p SolverGMRES.
 * 的最大临时向量数。
 * 为了有一个标准化的求解器构造方式，每个求解器类都有一个<tt>struct
 * AdditionalData</tt>作为成员，并且所有求解器类的构造函数都接受这样一个参数。有些求解器不需要额外的数据，或者在当前不需要。对于这些求解器，结构
 * @p AdditionalData
 * 是空的，调用构造函数时可以不给额外的结构作为参数，因为默认
 * @p AdditionalData被设置。
 * 有了这个，创建一个求解器看起来就像下面的一个块。
 *
 * @code
 * // GMRES with restart every 50 iterations
 * SolverGMRES solver_gmres (solver_control, vector_memory,
 *                           SolverGMRES::AdditionalData(50));
 *
 * // Richardson with omega=0.8
 * SolverRichardson solver_richardson (solver_control, vector_memory,
 *                                     SolverGMRES::AdditionalData(0.8));
 *
 * // CG with default AdditionalData
 * SolverCG solver_cg (solver_control, vector_memory);
 * @endcode
 *
 * 对所有求解器使用统一的构造参数列表，支持 @p
 * SolverSelector类；统一的接口使我们能够不变地使用这个类，即使某个求解器的参数类型数量发生了变化，仍然可以以简单的方式为每个可能使用的求解器向
 * @p SolverSelector  对象提供这些附加数据。
 *
 *  <h3>Observing the progress of linear solver iterations</h3>
 * SolverBase类是所有迭代求解器（如SolverCG、SolverGMRES等）的基类，它提供了一些设施，实际的求解器实现通过这些设施来确定迭代是否已经收敛、尚未收敛或已经失败。通常，这是用一个SolverControl类型的对象来完成的，该对象被传递给求解器类的构造器，并从它们传递给该基类的构造器。每一个解决线性问题的教程程序（从
 * step-3
 * 开始）都使用这个方法，并且在那里有详细的描述。然而，底层机制更为普遍，允许许多其他用途，以观察线性求解器的迭代进展。
 * 基本的方法是，迭代求解器在每次迭代结束时调用一个<i>signal</i>来确定解是否收敛。信号是一个类，从概念上讲，它有一个指向函数的列表，每次信号被调用时，这些函数都会被调用。在信号的语言中，被调用的函数被称为<i>slots</i>，人们可以在一个信号上附加任何数量的槽。我们在这里使用的信号和槽的实现是BOOST.signals2库中的一个）。一些细节可以说明下面发生了什么。
 *
 *
 *
 * - 在现实中，信号对象并不存储指向函数的指针，而是作为槽的函数对象。每个槽必须符合一个特定的签名：在这里，它是一个可以用三个参数（当前线性迭代的编号、当前残差和当前迭代；更多的细节在connect()函数的文档中讨论）来调用的对象。一个指向具有该参数列表的函数的指针就满足了要求，但你也可以传递一个成员函数，其 <code>this</code> 参数已经用lambda函数绑定（见下面的例子）。
 *
 *
 *
 * - 每个槽将返回一个值，表示迭代是否应该继续，是否应该因为成功而停止，或者因为失败而停止。因此，槽的返回类型是 SolverControl::State. ，所有槽的返回值在返回给调用信号的迭代求解器之前必须进行合并。其工作方式是，如果至少有一个槽返回 SolverControl::failure, ，那么合并后的值就是 SolverControl::failure; ，否则，如果至少有一个槽返回 SolverControl::iterate, ，那么这将是信号的返回值；最后，只有当所有槽返回 SolverControl::success 时，信号的返回值才是 SolverControl::success.  。
 *
 *
 *
 * - 当然，也有可能某个槽连接到信号只是为了观察解决方案或其特定部分的收敛情况，但对是否应该继续迭代没有特别的意见。在这种情况下，槽应该只是返回 SolverControl::success, ，根据上面的规则，这是所有返回值中最弱的一个。
 * 鉴于这一切，现在应该很明显，SolverControl对象如何适合这个方案：当一个SolverControl对象被传递给当前类的构造函数时，我们只需将该对象的
 * SolverControl::check()
 * 函数作为槽连接到我们在这里维护的信号。换句话说，由于SolverBase对象总是使用SolverControl对象构造的，所以总是至少有一个槽与信号相关，即决定收敛的槽。
 * 另一方面，使用connect()成员函数，可以将任何数量的其他槽连接到信号上，以观察你想要观察的东西。connect()函数还返回一个描述从信号到槽的连接的对象，然后相应的BOOST函数允许你在需要时断开槽的连接。
 * 一个例子可以说明这些问题。在 step-3
 * 的教程程序中，让我们在主类中添加一个成员函数如下。
 *
 * @code
 * SolverControl::State
 * Step3::write_intermediate_solution (
 * const unsigned int    iteration,
 * const double          , //check_value
 * const Vector<double> &current_iterate) const
 * {
 * DataOut<2> data_out;
 * data_out.attach_dof_handler (dof_handler);
 * data_out.add_data_vector (current_iterate, "solution");
 * data_out.build_patches ();
 *
 * std::ofstream output ("solution-"
 *                       + Utilities::int_to_string(iteration,4)
 *                       + ".vtu");
 * data_out.write_vtu (output);
 *
 * return SolverControl::success;
 * }
 * @endcode
 * 该函数满足了成为上面讨论的信号槽所需的签名，但例外的是它是一个成员函数，因此需要一个
 * <code>this</code>
 * 指针。该函数所做的是将给定的向量作为最后一个参数，并将其写入一个VTU格式的文件中，文件名由迭代的次数得出。
 * 这个函数可以通过修改 <code>Step3::solve()</code>
 * 函数来连接到CG求解器中，如下所示。
 *
 * @code
 * void Step3::solve ()
 * {
 * SolverControl           solver_control (1000, 1e-12);
 * SolverCG<>              solver (solver_control);
 *
 * solver.connect ([this](const unsigned int iteration,
 *                        const double check_value,
 *                        const Vector<double>current_iterate){
 *                   this->write_intermediate_solution(
 *                     iteration, check_value, current_iterate);
 *                 });
 * solver.solve (system_matrix, solution, system_rhs,
 *               PreconditionIdentity());
 * }
 * @endcode
 * 这里使用lambda函数确保我们将有三个参数加上
 * <code>this</code>
 * 指针的成员函数转换为只取三个参数的函数，方法是将该函数的隐含
 * <code>this</code> 参数固定为当前函数中的 <code>this</code>
 * 的指针。
 * 众所周知，CG方法是一种平滑迭代（与更常用的Jacobi或SSOR迭代是平滑器的方式相同）。因此，上面的代码可以观察到解决方案在每次迭代中是如何变得越来越平滑的。这最好是通过用随机分布的数字初始化解向量来观察
 * $[-1,1]$ ，使用如下代码
 *
 * @code
 * for (unsigned int i=0; i<solution.size(); ++i)
 *   solution(i) = 2.*rand()/RAND_MAX-1;
 * @endcode
 * 利用这一点，槽将产生文件，当可视化时，在迭代0到5的过程中看起来像这样。
 * <table> <tr> <td>
 @image html "cg-monitor-smoothing-0.png"
 * </td> <td>
 @image html "cg-monitor-smoothing-1.png"
 * </td> <td>
 @image html "cg-monitor-smoothing-2.png"
 * </td> </tr> <tr> <td>
 @image html "cg-monitor-smoothing-3.png"
 * </td> <td>
 @image html "cg-monitor-smoothing-4.png"
 * </td> <td>
 @image html "cg-monitor-smoothing-5.png"
 * </td> </tr> </table>
 *
 *
 * @ingroup Solvers
 *
 */
template <class VectorType = Vector<double>>
class SolverBase : public Subscriptor
{
public:
  /**
   * 底层向量类型的一个别名
   *
   */
  using vector_type = VectorType;

  /**
   * 构造函数。接受一个评估收敛条件的控制对象，以及一个允许求解器为临时对象分配内存的对象。
   * 在这两个对象中，都存储了一个引用，所以用户有责任保证这两个参数的寿命至少与求解器对象的寿命一样长。
   *
   */
  SolverBase(SolverControl &           solver_control,
             VectorMemory<VectorType> &vector_memory);

  /**
   * 构造函数。接受一个控制对象，该对象评估收敛的条件。与其他构造函数不同，该构造函数指定一个GrowingVectorMemory类型的内部对象来分配内存。
   * 对控制对象的引用被存储，所以用户有责任保证该参数的寿命至少与求解器对象的寿命一样长。
   *
   */
  SolverBase(SolverControl &solver_control);

  /**
   * 连接一个将在迭代求解器中被定期调用的函数对象。这个函数用于将监视器附加到迭代求解器上，以确定何时发生收敛，或者只是观察迭代的进度。更多信息请参见该类的文档。
   * @param  槽
   * 这里指定的函数对象将在每次调用时接收当前迭代的编号、用于检查收敛的值（通常是当前迭代相对于要解决的线性系统的残差）以及当前迭代的当前最佳可用猜测。请注意，有些求解器并不是在每次迭代中都更新近似解，而只是在确定收敛或失败后才更新（GMRES就是一个例子）；在这种情况下，作为信号的最后一个参数传递的向量只是在信号被调用时的最佳近似值，而不是在信号的返回值表明迭代应该终止时将返回的向量。函数对象必须返回一个
   * SolverControl::State
   * 值，表明迭代是否应该继续，已经失败，或者已经成功。然后，所有连接的函数的结果将被结合起来，以确定迭代应该发生什么。
   * @return
   * 一个连接对象，表示从信号到函数对象的连接。它可以用来再次断开函数对象与信号的连接。关于连接管理的更多信息，请参见BOOST
   * Signals2库的文档。
   *
   */
  boost::signals2::connection
  connect(
    const std::function<SolverControl::State(const unsigned int iteration,
                                             const double       check_value,
                                             const VectorType &current_iterate)>
      &slot);



protected:
  /**
   * 一个静态的向量内存对象，只要没有给构造函数提供这样的对象，就可以使用。
   *
   */
  mutable GrowingVectorMemory<VectorType> static_vector_memory;

  /**
   * 对一个为辅助向量提供内存的对象的引用。
   *
   */
  VectorMemory<VectorType> &memory;

private:
  /**
   * 一个类，其operator()结合了两个状态，表明我们应该继续迭代还是停止，并返回一个占优势的状态。其规则是
   *
   *
   *
   *
   *
   * - 如果两个状态中的一个表示失败，则返回失败。
   *
   *
   *
   *
   *
   * - 否则，如果两个状态中的一个表示要继续迭代，那么就继续迭代。
   *
   *
   *
   *
   *
   *
   * - 否则，返回成功。
   *
   */
  struct StateCombiner
  {
    using result_type = SolverControl::State;

    SolverControl::State
    operator()(const SolverControl::State state1,
               const SolverControl::State state2) const;

    template <typename Iterator>
    SolverControl::State
    operator()(const Iterator begin, const Iterator end) const;
  };

protected:
  /**
   * 一个信号，迭代求解器可以在每次迭代结束时（或以其他周期性方式）执行，以了解我们是否应该继续迭代。该信号可以调用一个或多个槽，每个槽都会自己做出这个判断，而所有槽的结果（函数调用）将由StateCombiner对象决定。
   * 传递给信号的参数是：（i）当前迭代的编号；（ii）用于确定收敛的值（通常是残差，但在其他情况下可以使用其他量，只要它们在迭代接近线性系统的解时收敛为零）；以及（iii）对应于信号被调用时的解的当前最佳猜测的矢量。请注意，有些求解器并不是在每次迭代中都更新近似解，而只是在确定收敛或失败后才更新（GMRES就是一个例子）；在这种情况下，作为信号的最后一个参数传递的向量只是在信号被调用时的最佳近似值，而不是在信号的返回值表明迭代应该终止时将返回的向量。
   *
   */
  boost::signals2::signal<
    SolverControl::State(const unsigned int iteration,
                         const double       check_value,
                         const VectorType & current_iterate),
    StateCombiner>
    iteration_status;
};



 /*-------------------------------- Inline functions ------------------------*/ 


template <class VectorType>
inline SolverControl::State
SolverBase<VectorType>::StateCombiner::
operator()(const SolverControl::State state1,
           const SolverControl::State state2) const
{
  if ((state1 == SolverControl::failure) || (state2 == SolverControl::failure))
    return SolverControl::failure;
  else if ((state1 == SolverControl::iterate) ||
           (state2 == SolverControl::iterate))
    return SolverControl::iterate;
  else
    return SolverControl::success;
}


template <class VectorType>
template <typename Iterator>
inline SolverControl::State
SolverBase<VectorType>::StateCombiner::operator()(const Iterator begin,
                                                  const Iterator end) const
{
  Assert(begin != end,
         ExcMessage("You can't combine iterator states if no state is given."));

  // combine the first with all of the following states
  SolverControl::State state = *begin;
  Iterator             p     = begin;
  ++p;
  for (; p != end; ++p)
    state = this->operator()(state, *p);

  return state;
}


template <class VectorType>
inline SolverBase<VectorType>::SolverBase(
  SolverControl &           solver_control,
  VectorMemory<VectorType> &vector_memory)
  : memory(vector_memory)
{
  // connect the solver control object to the signal. SolverControl::check
  // only takes two arguments, the iteration and the check_value, and so
  // we simply ignore the third argument that is passed in whenever the
  // signal is executed
  connect([&solver_control](const unsigned int iteration,
                            const double       check_value,
                            const VectorType &) {
    return solver_control.check(iteration, check_value);
  });
}



template <class VectorType>
inline SolverBase<VectorType>::SolverBase(SolverControl &solver_control)
  : // use the static memory object this class owns
  memory(static_vector_memory)
{
  // connect the solver control object to the signal. SolverControl::check
  // only takes two arguments, the iteration and the check_value, and so
  // we simply ignore the third argument that is passed in whenever the
  // signal is executed
  connect([&solver_control](const unsigned int iteration,
                            const double       check_value,
                            const VectorType &) {
    return solver_control.check(iteration, check_value);
  });
}



template <class VectorType>
inline boost::signals2::connection
SolverBase<VectorType>::connect(
  const std::function<SolverControl::State(const unsigned int iteration,
                                           const double       check_value,
                                           const VectorType & current_iterate)>
    &slot)
{
  return iteration_status.connect(slot);
}



DEAL_II_NAMESPACE_CLOSE

#endif


