//include/deal.II-translator/base/exceptions_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_exceptions_h
#define dealii_exceptions_h

#include <deal.II/base/config.h>

#include <exception>
#include <ostream>
#include <string>
#include <type_traits>

#ifdef DEAL_II_WITH_CUDA
#  include <cusolverSp.h>
#  include <cusparse.h>
#endif


DEAL_II_NAMESPACE_OPEN


/**
 * 这个类是所有异常类的基类。不要直接使用它的方法和变量，因为它的接口和机制可能会被改变。而是使用<tt>DeclException</tt>宏系列创建新的异常类。
 * 参见
 * @ref Exceptions
 * 模块以了解关于这个类的更多细节，以及从它派生的类可以做什么。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
class ExceptionBase : public std::exception
{
public:
  /**
   * 默认构造函数。
   *
   */
  ExceptionBase();

  /**
   * 复制构造函数。
   *
   */
  ExceptionBase(const ExceptionBase &exc);

  /**
   * 解构器。
   *
   */
  virtual ~ExceptionBase() noexcept override;

  /**
   * 拷贝操作符。这个操作符被删除，因为异常对象是不可复制的。
   *
   */
  ExceptionBase
  operator=(const ExceptionBase &) = delete;

  /**
   * 设置异常出现的文件名和行，以及被违反的条件和作为char指针的异常名称。这个函数还可以填充堆栈跟踪。
   *
   */
  void
  set_fields(const char *file,
             const int   line,
             const char *function,
             const char *cond,
             const char *exc_name);


  /**
   * 覆盖标准函数，返回错误的描述。
   *
   */
  virtual const char *
  what() const noexcept override;

  /**
   * 获取异常名称。
   *
   */
  const char *
  get_exc_name() const;

  /**
   * 打印出错误信息的一般部分。
   *
   */
  void
  print_exc_data(std::ostream &out) const;

  /**
   * 打印关于发生的异常的更具体的信息。
   * 在你自己的异常类中重载这个函数。
   *
   */
  virtual void
  print_info(std::ostream &out) const;

  /**
   * 打印一个堆栈跟踪，如果之前有记录的话，到给定的流中。
   *
   */
  void
  print_stack_trace(std::ostream &out) const;

protected:
  /**
   * 这个异常所发生的文件的名称。
   *
   */
  const char *file;

  /**
   * 该文件中的行号。
   *
   */
  unsigned int line;

  /**
   * 函数的名称，漂亮的打印。
   *
   */
  const char *function;

  /**
   * 被违反的条件，作为一个字符串。
   *
   */
  const char *cond;

  /**
   * 异常的名称和调用序列。
   *
   */
  const char *exc;

  /**
   * 回溯到问题发生的位置，如果系统支持的话。
   *
   */
  mutable char **stacktrace;

  /**
   * 存储在前一个变量中的堆栈跟踪帧的数量。
   * 如果系统不支持堆栈跟踪，则为零。
   *
   */
  int n_stacktrace_frames;

#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  /**
   * 包含原始堆栈跟踪的指针数组
   *
   */
  void *raw_stacktrace[25];
#endif

private:
  /**
   * 产生c_string的内部函数。由what()调用。
   *
   */
  void
  generate_message() const;

  /**
   * 一个指向将被what()打印的c_string的指针。它是由generate_message()填充的。
   *
   */
  mutable std::string what_str;
};

#ifndef DOXYGEN

/**
 * 声明一个从ExceptionBase派生出来的没有参数的异常类。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，它们没有用一个字符串作为前缀，这很可能使它们成为deal.II独有的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException0(Exception0)                \
    class Exception0 : public dealii::ExceptionBase \
    {}


/**
 * 声明一个从 ExceptionBase
 * 派生的异常类，它可以接受一个运行时参数，但如果在你想抛出异常的地方没有给出参数，它就会简单地恢复到通过这个宏声明异常类时提供的默认文本。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它们没有用一个字符串作为前缀，这可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclExceptionMsg(Exception, defaulttext)    \
    class Exception : public dealii::ExceptionBase    \
    {                                                 \
    public:                                           \
      Exception(const std::string &msg = defaulttext) \
        : arg(msg)                                    \
      {}                                              \
      virtual ~Exception() noexcept                   \
      {}                                              \
      virtual void                                    \
      print_info(std::ostream &out) const override    \
      {                                               \
        out << "    " << arg << std::endl;            \
      }                                               \
                                                      \
    private:                                          \
      const std::string arg;                          \
    }

/**
 * 声明一个从ExceptionBase派生出来的异常类，并增加一个参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException1(Exception1, type1, outsequence) \
    class Exception1 : public dealii::ExceptionBase      \
    {                                                    \
    public:                                              \
      Exception1(type1 const &a1)                        \
        : arg1(a1)                                       \
      {}                                                 \
      virtual ~Exception1() noexcept                     \
      {}                                                 \
      virtual void                                       \
      print_info(std::ostream &out) const override       \
      {                                                  \
        out << "    " outsequence << std::endl;          \
      }                                                  \
                                                         \
    private:                                             \
      type1 const arg1;                                  \
    }


/**
 * 声明一个从ExceptionBase派生出来的异常类，有两个额外的参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException2(Exception2, type1, type2, outsequence) \
    class Exception2 : public dealii::ExceptionBase             \
    {                                                           \
    public:                                                     \
      Exception2(type1 const &a1, type2 const &a2)              \
        : arg1(a1)                                              \
        , arg2(a2)                                              \
      {}                                                        \
      virtual ~Exception2() noexcept                            \
      {}                                                        \
      virtual void                                              \
      print_info(std::ostream &out) const override              \
      {                                                         \
        out << "    " outsequence << std::endl;                 \
      }                                                         \
                                                                \
    private:                                                    \
      type1 const arg1;                                         \
      type2 const arg2;                                         \
    }


/**
 * 声明一个从ExceptionBase派生出来的异常类，有三个附加参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中独一无二。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException3(Exception3, type1, type2, type3, outsequence) \
    class Exception3 : public dealii::ExceptionBase                    \
    {                                                                  \
    public:                                                            \
      Exception3(type1 const &a1, type2 const &a2, type3 const &a3)    \
        : arg1(a1)                                                     \
        , arg2(a2)                                                     \
        , arg3(a3)                                                     \
      {}                                                               \
      virtual ~Exception3() noexcept                                   \
      {}                                                               \
      virtual void                                                     \
      print_info(std::ostream &out) const override                     \
      {                                                                \
        out << "    " outsequence << std::endl;                        \
      }                                                                \
                                                                       \
    private:                                                           \
      type1 const arg1;                                                \
      type2 const arg2;                                                \
      type3 const arg3;                                                \
    }


/**
 * 声明一个从ExceptionBase派生出来的异常类，有四个附加参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
    class Exception4 : public dealii::ExceptionBase                           \
    {                                                                         \
    public:                                                                   \
      Exception4(type1 const &a1,                                             \
                 type2 const &a2,                                             \
                 type3 const &a3,                                             \
                 type4 const &a4)                                             \
        : arg1(a1)                                                            \
        , arg2(a2)                                                            \
        , arg3(a3)                                                            \
        , arg4(a4)                                                            \
      {}                                                                      \
      virtual ~Exception4() noexcept                                          \
      {}                                                                      \
      virtual void                                                            \
      print_info(std::ostream &out) const override                            \
      {                                                                       \
        out << "    " outsequence << std::endl;                               \
      }                                                                       \
                                                                              \
    private:                                                                  \
      type1 const arg1;                                                       \
      type2 const arg2;                                                       \
      type3 const arg3;                                                       \
      type4 const arg4;                                                       \
    }


/**
 * 声明一个从ExceptionBase派生出来的异常类，有五个额外的参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中独一无二。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有用
 * <code>DEAL</code> or <code>deal</code>
 * 作为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException5(                                       \
    Exception5, type1, type2, type3, type4, type5, outsequence) \
    class Exception5 : public dealii::ExceptionBase             \
    {                                                           \
    public:                                                     \
      Exception5(type1 const &a1,                               \
                 type2 const &a2,                               \
                 type3 const &a3,                               \
                 type4 const &a4,                               \
                 type5 const &a5)                               \
        : arg1(a1)                                              \
        , arg2(a2)                                              \
        , arg3(a3)                                              \
        , arg4(a4)                                              \
        , arg5(a5)                                              \
      {}                                                        \
      virtual ~Exception5() noexcept                            \
      {}                                                        \
      virtual void                                              \
      print_info(std::ostream &out) const override              \
      {                                                         \
        out << "    " outsequence << std::endl;                 \
      }                                                         \
                                                                \
    private:                                                    \
      type1 const arg1;                                         \
      type2 const arg2;                                         \
      type3 const arg3;                                         \
      type4 const arg4;                                         \
      type5 const arg5;                                         \
    }

#else  /*ifndef DOXYGEN*/ 

// Dummy definitions for doxygen:

/**
 * 声明一个从ExceptionBase派生出来的没有参数的异常类。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有用
 * <code>DEAL</code> or <code>deal</code>
 * 作为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException0(Exception0) \
     /** @ingroup Exceptions */        \
    static dealii::ExceptionBase &Exception0()

/**
 * 声明一个从 ExceptionBase
 * 派生的异常类，它可以接受一个运行时参数，但如果在你想抛出异常的地方没有给出参数，它就会简单地恢复到通过这个宏声明异常类时提供的默认文本。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它们没有用一个字符串作为前缀，这可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclExceptionMsg(Exception, defaulttext) \
     /** @ingroup Exceptions */                      \
     /** @dealiiExceptionMessage{defaulttext} */     \
    static dealii::ExceptionBase &Exception()

/**
 * 声明一个从ExceptionBase派生出来的异常类，并增加一个参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中独一无二。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有用
 * <code>DEAL</code> or <code>deal</code>
 * 作为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException1(Exception1, type1, outsequence) \
     /** @ingroup Exceptions */                            \
     /** @dealiiExceptionMessage{outsequence} */           \
    static dealii::ExceptionBase &Exception1(type1 arg1)


/**
 * 声明一个从ExceptionBase派生出来的异常类，有两个额外的参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果它们的前缀不是
 * <code>DEAL</code> or <code>deal</code>
 * ，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException2(Exception2, type1, type2, outsequence) \
     /** @ingroup Exceptions */                                   \
     /** @dealiiExceptionMessage{outsequence} */                  \
    static dealii::ExceptionBase &Exception2(type1 arg1, type2 arg2)


/**
 * 声明一个从ExceptionBase派生出来的异常类，有三个额外的参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中独一无二。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException3(Exception3, type1, type2, type3, outsequence) \
     /** @ingroup Exceptions */                                          \
     /** @dealiiExceptionMessage{outsequence} */                         \
    static dealii::ExceptionBase &Exception3(type1 arg1, type2 arg2, type3 arg3)


/**
 * 声明一个从ExceptionBase派生出来的异常类，有四个附加参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果它们的前缀不是
 * <code>DEAL</code> or <code>deal</code>
 * ，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException4(Exception4, type1, type2, type3, type4, outsequence) \
     /** @ingroup Exceptions */                                                 \
     /** @dealiiExceptionMessage{outsequence} */                                \
    static dealii::ExceptionBase &Exception4(type1 arg1,                      \
                                             type2 arg2,                      \
                                             type3 arg3,                      \
                                             type4 arg4)


/**
 * 声明一个从ExceptionBase派生出来的异常类，有五个额外的参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，这很可能使它们在deal.II中独一无二。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define DeclException5(                                       \
    Exception5, type1, type2, type3, type4, type5, outsequence) \
     /** @ingroup Exceptions */                                   \
     /** @dealiiExceptionMessage{outsequence} */                  \
    static dealii::ExceptionBase &Exception5(                   \
      type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5)

#endif  /*ifndef DOXYGEN*/ 


/**
 * 声明一些重复出现的异常。这样，你就可以简单地使用这些异常，而不必在你的类中本地声明它们。声明这些异常的名字空间后来被包含在全局名字空间中，由
 *
 * @code
 * using namespace StandardExceptions;
 * @endcode
 *
 *
 *
 * @ingroup Exceptions
 *
 *
 */
namespace StandardExceptions
{
  /**
   * @addtogroup  Exceptions
   *
   */
  //@{

  /**
   * 表示除以0的异常情况。
   *
   */
  DeclExceptionMsg(ExcDivideByZero,
                   "A piece of code is attempting a division by zero. This is "
                   "likely going to lead to results that make no sense.");

  /**
   * 如果一个数字不是有限的，就会产生异常。
   * 这个异常应该被用来捕获不是除以零的算术运算所产生的无限数或非数的结果（使用ExcDivideByZero来处理这些结果）。
   * 这个异常使用 std::complex
   * 作为它的参数，以确保我们可以对所有标量参数（实数或复数）使用它。
   *
   */
  DeclException1(
    ExcNumberNotFinite,
    std::complex<double>,
    << "In a significant number of places, deal.II checks that some intermediate "
    << "value is a finite number (as opposed to plus or minus infinity, or "
    << "NaN/Not a Number). In the current function, we encountered a number "
    << "that is not finite (its value is " << arg1 << " and therefore "
    << "violates the current assertion).\n\n"
    << "This may be due to the fact that some operation in this function "
    << "created such a value, or because one of the arguments you passed "
    << "to the function already had this value from some previous "
    << "operation. In the latter case, this function only triggered the "
    << "error but may not actually be responsible for the computation of "
    << "the number that is not finite.\n\n"
    << "There are two common cases where this situation happens. First, your "
    << "code (or something in deal.II) divides by zero in a place where this "
    << "should not happen. Or, you are trying to solve a linear system "
    << "with an unsuitable solver (such as an indefinite or non-symmetric "
    << "linear system using a Conjugate Gradient solver); such attempts "
    << "oftentimes yield an operation somewhere that tries to divide "
    << "by zero or take the square root of a negative value.\n\n"
    << "In any case, when trying to find the source of the error, "
    << "recall that the location where you are getting this error is "
    << "simply the first place in the program where there is a check "
    << "that a number (e.g., an element of a solution vector) is in fact "
    << "finite, but that the actual error that computed the number "
    << "may have happened far earlier. To find this location, you "
    << "may want to add checks for finiteness in places of your "
    << "program visited before the place where this error is produced. "
    << "One way to check for finiteness is to use the 'AssertIsFinite' "
    << "macro.");

  /**
   * 试图分配一个新的对象，由于缺乏可用的内存而失败。
   *
   */
  DeclException1(ExcOutOfMemory,
                 std::size_t,
                 "Your program tried to allocate some memory but this "
                 "allocation failed. Typically, this either means that "
                 "you simply do not have enough memory in your system, "
                 "or that you are (erroneously) trying to allocate "
                 "a chunk of memory that is simply beyond all reasonable "
                 "size, for example because the size of the object has "
                 "been computed incorrectly."
                 "\n\n"
                 "In the current case, the request was for "
                   << arg1 << " bytes.");

  /**
   * 一个内存处理程序到达一个点，所有分配的对象都应该被释放。由于这个异常被抛出，一些对象仍然被分配。
   *
   */
  DeclException1(ExcMemoryLeak,
                 int,
                 << "Destroying memory handler while " << arg1
                 << " objects are still allocated.");

  /**
   * 读取或写入一个文件时发生错误。
   *
   */
  DeclExceptionMsg(ExcIO,
                   "An input/output error has occurred. There are a number of "
                   "reasons why this may be happening, both for reading and "
                   "writing operations."
                   "\n\n"
                   "If this happens during an operation that tries to read "
                   "data: First, you may be "
                   "trying to read from a file that doesn't exist or that is "
                   "not readable given its file permissions. Second, deal.II "
                   "uses this error at times if it tries to "
                   "read information from a file but where the information "
                   "in the file does not correspond to the expected format. "
                   "An example would be a truncated file, or a mesh file "
                   "that contains not only sections that describe the "
                   "vertices and cells, but also sections for additional "
                   "data that deal.II does not understand."
                   "\n\n"
                   "If this happens during an operation that tries to write "
                   "data: you may be trying to write to a file to which file "
                   "or directory permissions do not allow you to write. A "
                   "typical example is where you specify an output file in "
                   "a directory that does not exist.");

  /**
   * 发生了一个打开命名文件的错误。
   * 构造函数需要一个类型为 <tt>std::string</tt>
   * 的单一参数来命名文件。
   *
   */
  DeclException1(ExcFileNotOpen,
                 std::string,
                 << "Could not open file " << arg1
                 << "."
                    "\n\n"
                    "If this happens during an operation that tries to read "
                    "data: you may be "
                    "trying to read from a file that doesn't exist or that is "
                    "not readable given its file permissions."
                    "\n\n"
                    "If this happens during an operation that tries to write "
                    "data: you may be trying to write to a file to which file "
                    "or directory permissions do not allow you to write. A "
                    "typical example is where you specify an output file in "
                    "a directory that does not exist.");

  /**
   * 表示库或应用程序的一部分尚未实现的异常。在许多情况下，这只表明还没有太多需要的东西，而不是说这很难实现。因此，相当值得去看一看相应的地方，看看是否可以不费吹灰之力就能实现它。
   *
   */
  DeclExceptionMsg(ExcNotImplemented,
                   "You are trying to use functionality in deal.II that is "
                   "currently not implemented. In many cases, this indicates "
                   "that there simply didn't appear much of a need for it, or "
                   "that the author of the original code did not have the "
                   "time to implement a particular case. If you hit this "
                   "exception, it is therefore worth the time to look into "
                   "the code to find out whether you may be able to "
                   "implement the missing functionality. If you do, please "
                   "consider providing a patch to the deal.II development "
                   "sources (see the deal.II website on how to contribute).");

  /**
   * 这种异常通常表明，程序员认为在算法的某一点上必须满足的某些条件没有得到满足。这可能是由于上面的一些编程错误，由于对算法的修改没有保留这个断言，或者由于程序员所做的假设根本就不成立（也就是说，虽然这里没有错误，但却抛出了这个异常）。在库中，当我们写了某种复杂的算法，并且还不确定我们是否得到了正确的算法时，这种异常是最常用的；然后我们在算法的每个部分之后加入断言，检查一些应该成立的条件，如果不成立就抛出一个异常。
   * 即使我们确信实现是正确的，我们通常也会留下这些断言，因为如果后来有人改变或扩展了算法，这些异常会向他们表明他们是否违反了算法中后面使用的假设。此外，有时会发生这样的情况，即算法在非常罕见的角落情况下不工作。那么这些情况迟早会被异常所困住，这样一来，算法也可以针对这些情况进行修正。
   *
   */
  DeclExceptionMsg(ExcInternalError,
                   "This exception -- which is used in many places in the "
                   "library -- usually indicates that some condition which "
                   "the author of the code thought must be satisfied at a "
                   "certain point in an algorithm, is not fulfilled. An "
                   "example would be that the first part of an algorithm "
                   "sorts elements of an array in ascending order, and "
                   "a second part of the algorithm later encounters an "
                   "element that is not larger than the previous one."
                   "\n\n"
                   "There is usually not very much you can do if you "
                   "encounter such an exception since it indicates an error "
                   "in deal.II, not in your own program. Try to come up with "
                   "the smallest possible program that still demonstrates "
                   "the error and contact the deal.II mailing lists with it "
                   "to obtain help.");

  /**
   * 这个异常用于可能不被调用的函数（即在纯函数中），但不能被声明为纯函数，因为无论如何都要使用该类，尽管只有在使用派生类时才可能调用相应的函数。
   *
   */
  DeclExceptionMsg(
    ExcPureFunctionCalled,
    "You (or a place in the library) are trying to call a "
    "function that is declared as a virtual function in a "
    "base class but that has not been overridden in your "
    "derived class."
    "\n\n"
    "This exception happens in cases where the base class "
    "cannot provide a useful default implementation for "
    "the virtual function, but where we also do not want "
    "to mark the function as abstract (i.e., with '=0' at the end) "
    "because the function is not essential to the class in many "
    "contexts. In cases like this, the base class provides "
    "a dummy implementation that makes the compiler happy, but "
    "that then throws the current exception."
    "\n\n"
    "A concrete example would be the 'Function' class. It declares "
    "the existence of 'value()' and 'gradient()' member functions, "
    "and both are marked as 'virtual'. Derived classes have to "
    "override these functions for the values and gradients of a "
    "particular function. On the other hand, not every function "
    "has a gradient, and even for those that do, not every program "
    "actually needs to evaluate it. Consequently, there is no "
    "*requirement* that a derived class actually override the "
    "'gradient()' function (as there would be had it been marked "
    "as abstract). But, since the base class cannot know how to "
    "compute the gradient, if a derived class does not override "
    "the 'gradient()' function and it is called anyway, then the "
    "default implementation in the base class will simply throw "
    "an exception."
    "\n\n"
    "The exception you see is what happens in cases such as the "
    "one just illustrated. To fix the problem, you need to "
    "investigate whether the function being called should indeed have "
    "been called; if the answer is 'yes', then you need to "
    "implement the missing override in your class.");

  /**
   * 如果发现某些对象未被初始化，就会使用这个异常。
   *
   */
  DeclException0(ExcNotInitialized);

  /**
   * 该对象处于一个不适合该操作的状态。
   *
   */
  DeclException0(ExcInvalidState);

  /**
   * 如果一个功能在给定的维度上不可能实现，就会引发这个异常。主要用于抛出1D中的函数调用。
   * 构造函数接收一个<tt>int</tt>，表示维度。
   *
   */
  DeclException1(ExcImpossibleInDim,
                 int,
                 << "You are trying to execute functionality that is "
                 << "impossible in " << arg1
                 << "d or simply does not make any sense.");

  /**
   * 如果一个功能在给定的维度和空间维度的组合中是不可能的，就会引发这个异常。
   * 构造函数接收两个<tt>int</tt>，表示维度和空间维度。
   *
   */
  DeclException2(ExcImpossibleInDimSpacedim,
                 int,
                 int,
                 << "You are trying to execute functionality that is "
                 << "impossible in dimensions <" << arg1 << "," << arg2
                 << "> or simply does not make any sense.");


  /**
   * 一个数字是零，但它不应该在这里出现。
   *
   */
  DeclExceptionMsg(ExcZero,
                   "In a check in the code, deal.II encountered a zero in "
                   "a place where this does not make sense. See the condition "
                   "that was being checked and that is printed further up "
                   "in the error message to get more information on what "
                   "the erroneous zero corresponds to.");

  /**
   * 在这个成员函数被调用之前，对象应该已经被填满了东西。
   *
   */
  DeclExceptionMsg(ExcEmptyObject,
                   "The object you are trying to access is empty but it makes "
                   "no sense to attempt the operation you are trying on an "
                   "empty object.");

  /**
   * 当两个对象的大小被认为是相等的，但却不相等时，就会引发这个异常。
   * 构造函数的参数是第一个和第二个大小，都是<tt>int</tt>类型。
   *
   */
  DeclException2(ExcDimensionMismatch,
                 std::size_t,
                 std::size_t,
                 << "Dimension " << arg1 << " not equal to " << arg2 << ".");

  /**
   * 第一尺寸应该等于第二尺寸或第三尺寸，但它都不是。
   *
   */
  DeclException3(ExcDimensionMismatch2,
                 int,
                 int,
                 int,
                 << "Dimension " << arg1 << " neither equal to " << arg2
                 << " nor to " << arg3 << ".");

  /**
   * 这个异常表示一个索引不在预期范围内。  例如，可能是你试图访问一个不存在的向量的一个元素。    构造函数需要三个<tt>int</tt>参数，即  <ol>   <li>  违规索引  <li>  下限  <li>  上限加一  </ol>  。
   *
   */
  DeclException3(
    ExcIndexRange,
    int,
    int,
    int,
    << "Index " << arg1 << " is not in the half-open range [" << arg2 << ","
    << arg3 << ")."
    << (arg2 == arg3 ?
          " In the current case, this half-open range is in fact empty, "
          "suggesting that you are accessing an element of an empty "
          "collection such as a vector that has not been set to the "
          "correct size." :
          ""));

  /**
   * 这个异常表示某个指数不在预期范围内。  例如，可能是你试图访问一个不存在的向量的一个元素。    构造函数需要三个参数，即  <ol>   <li>  违规索引  <li>  下限  <li>  上限加一  </ol>  这个通用异常与ExcIndexRange不同，允许指定索引的类型。
   *
   */
  template <typename T>
  DeclException3(
    ExcIndexRangeType,
    T,
    T,
    T,
    << "Index " << arg1 << " is not in the half-open range [" << arg2 << ","
    << arg3 << ")."
    << (arg2 == arg3 ?
          " In the current case, this half-open range is in fact empty, "
          "suggesting that you are accessing an element of an empty "
          "collection such as a vector that has not been set to the "
          "correct size." :
          ""));

  /**
   * 一个数字太小。
   *
   */
  DeclException2(ExcLowerRange,
                 int,
                 int,
                 << "Number " << arg1 << " must be larger than or equal "
                 << arg2 << ".");

  /**
   * 上述ExcLowerRange的通用异常定义。
   *
   */
  template <typename T>
  DeclException2(ExcLowerRangeType,
                 T,
                 T,
                 << "Number " << arg1 << " must be larger than or equal "
                 << arg2 << ".");

  /**
   * 这个异常表示第一个参数应该是第二个参数的整数倍，但却不是。
   *
   */
  DeclException2(ExcNotMultiple,
                 int,
                 int,
                 << "Division " << arg1 << " by " << arg2
                 << " has remainder different from zero.");

  /**
   * 如果你访问的迭代器有损坏的数据，就会抛出这个异常。
   * 例如，可能是它所指的容器在迭代器所指的位置上没有一个入口。
   * 通常情况下，这将是deal.II的一个内部错误，因为增量和减量运算符不应该产生一个无效的迭代器。
   *
   */
  DeclExceptionMsg(ExcInvalidIterator,
                   "You are trying to use an iterator, but the iterator is "
                   "in an invalid state. This may indicate that the iterator "
                   "object has not been initialized, or that it has been "
                   "moved beyond the end of the range of valid elements.");

  /**
   * 如果你递增或递减的迭代器已经处于最终状态，就会抛出这个异常。
   *
   */
  DeclExceptionMsg(ExcIteratorPastEnd,
                   "You are trying to use an iterator, but the iterator is "
                   "pointing past the end of the range of valid elements. "
                   "It is not valid to dereference the iterator in this "
                   "case.");

  /**
   * 这个异常是围绕<tt>DeclException0</tt>宏中的一个设计缺陷进行的：通过DeclException0声明的异常不允许指定异常发生时显示的信息，而其他异常则允许与给定参数一起显示一个文本。
   * 当抛出这个异常时，你可以给一个消息作为
   * <tt>std::string</tt> 作为异常的参数，然后显示出来。
   * 当然，这个参数可以在运行时构建，例如包括一个无法打开的文件的名称，或者任何其他你可能想从不同的部分组合起来的文本。
   *
   */
  DeclException1(ExcMessage, std::string, << arg1);

  /**
   * 带有鬼魂元素的平行向量是只读向量。
   *
   */
  DeclExceptionMsg(ExcGhostsPresent,
                   "You are trying an operation on a vector that is only "
                   "allowed if the vector has no ghost elements, but the "
                   "vector you are operating on does have ghost elements. "
                   "Specifically, vectors with ghost elements are read-only "
                   "and cannot appear in operations that write into these "
                   "vectors."
                   "\n\n"
                   "See the glossary entry on 'Ghosted vectors' for more "
                   "information.");

  /**
   * 我们的一些数值类允许使用赋值运算符<tt>=</tt>将所有条目设置为零。
   * 在许多情况下，这个赋值运算符对参数0是有意义的<b>only</b>。在其他情况下，会抛出这种异常。
   *
   */
  DeclExceptionMsg(ExcScalarAssignmentOnlyForZeroValue,
                   "You are trying an operation of the form 'vector = C', "
                   "'matrix = C', or 'tensor = C' with a nonzero scalar value "
                   "'C'. However, such assignments are only allowed if the "
                   "C is zero, since the semantics for assigning any other "
                   "value are not clear. For example: one could interpret "
                   "assigning a matrix a value of 1 to mean the matrix has a "
                   "norm of 1, the matrix is the identity matrix, or the "
                   "matrix contains only 1s. Similar problems exist with "
                   "vectors and tensors. Hence, to avoid this ambiguity, such "
                   "assignments are not permitted.");

  /**
   * 这个函数需要对LAPACK库的支持。
   *
   */
  DeclExceptionMsg(
    ExcNeedsLAPACK,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use LAPACK, but cmake did not "
    "find a valid LAPACK library.");

  /**
   * 这个函数需要对MPI库的支持。
   *
   */
  DeclExceptionMsg(
    ExcNeedsMPI,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use MPI.");

  /**
   * 该函数需要对FunctionParser库的支持。
   *
   */
  DeclExceptionMsg(
    ExcNeedsFunctionparser,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use the function parser which "
    "relies on the muparser library, but cmake did not "
    "find a valid muparser library on your system and also did "
    "not choose the one that comes bundled with deal.II.");

  /**
   * 该函数需要对Assimp库的支持。
   *
   */
  DeclExceptionMsg(
    ExcNeedsAssimp,
    "You are attempting to use functionality that is only available "
    "if deal.II was configured to use Assimp, but cmake did not "
    "find a valid Assimp library.");

#ifdef DEAL_II_WITH_CUDA
  /**
   * 如果在CUDA内核中发生错误，会引发这个异常。
   * 构造函数接收一个<tt>char*</tt>，即cudaGetErrorString的输出。
   *
   */
  DeclException1(ExcCudaError, const char *, << arg1);
  /**
   * 如果在cuSPARSE函数中发生错误，就会引发这个异常。
   *
   */
  DeclException1(ExcCusparseError,
                 std::string,
                 << "There was an error in a cuSPARSE function: " << arg1);
#endif
  //@}

  /**
   * 这个函数需要对Exodus II库的支持。
   *
   */
  DeclExceptionMsg(
    ExcNeedsExodusII,
    "You are attempting to use functionality that is only available if deal.II "
    "was configured to use Trilinos' SEACAS library (which provides ExodusII), "
    "but cmake did not find a valid SEACAS library.");

#ifdef DEAL_II_WITH_MPI
  /**
   * MPI错误的异常。这个异常只有在 <code>deal.II</code>
   * 编译时支持MPI时才会被定义。这个异常应该和
   * <code>AssertThrow</code>
   * 一起使用，以检查MPI函数的错误代码。比如说。
   * @code
   * const int ierr = MPI_Isend(...);
   * AssertThrow(ierr == MPI_SUCCESS, ExcMPI(ierr));
   * @endcode
   * 或者，使用方便的宏  <code>AssertThrowMPI</code>  。
   * @code
   * const int ierr = MPI_Irecv(...);
   * AssertThrowMPI(ierr);
   * @endcode
   * 如果断言失败，那么错误代码将被用于利用
   * <code>MPI_Error_string</code>
   * 函数向屏幕打印一个有用的信息。
   * @ingroup Exceptions
   *
   */
  class ExcMPI : public dealii::ExceptionBase
  {
  public:
    ExcMPI(const int error_code);

    virtual void
    print_info(std::ostream &out) const override;

    const int error_code;
  };
#endif // DEAL_II_WITH_MPI



#ifdef DEAL_II_TRILINOS_WITH_SEACAS
  /**
   * 用于ExodusII错误的例外。只有在 <code>deal.II</code>
   * 编译时支持SEACAS的情况下才会定义这个异常，SEACAS可以通过Trilinos获得。这个函数应该和便利宏AssertThrowExodusII一起使用。
   * @ingroup Exceptions
   *
   */
  class ExcExodusII : public ExceptionBase
  {
  public:
    /**
     * 构造器。          @param  error_code
     * 由ExodusII函数返回的错误代码。
     *
     */
    ExcExodusII(const int error_code);

    /**
     * 打印错误的描述到给定的流中。
     *
     */
    virtual void
    print_info(std::ostream &out) const override;

    /**
     * 存储错误代码。
     *
     */
    const int error_code;
  };
#endif // DEAL_II_TRILINOS_WITH_SEACAS
}  /*namespace StandardExceptions*/ 



/**
 * 在这个命名空间中，与Assert和AssertThrow机制有关的函数被声明。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
namespace deal_II_exceptions
{
  namespace internals
  {
    /**
     * 将这个变量设置为false将禁用deal.II的异常机制来中止问题。Assert()宏将抛出异常，而AssertNothrow()宏将只打印错误信息。这个变量不应该被直接改变。使用
     * disable_abort_on_exception() 来代替。
     *
     */
    extern bool allow_abort_on_exception;
  } // namespace internals

  /**
   * 设置一个字符串，在表示触发<tt>Assert</tt>语句的信息输出时打印。这个字符串在通常的输出之外被打印出来，可能会指示出一些除非我们使用调试器否则不容易得到的信息。例如，对于集群计算机上的分布式程序，所有进程的输出都被重定向到同一个控制台窗口。在这种情况下，将程序运行的主机名称作为额外的名称是很方便的，这样就可以看到异常发生在程序的哪个实例中。
   * 参数所指向的字符串是复制的，所以在调用这个函数后不需要再存储。
   * 以前设置的附加输出被给这个函数的参数所取代。
   * @see  异常情况
   *
   */
  void
  set_additional_assert_output(const char *const p);

  /**
   * 当异常发生时，调用此函数可以禁止与其他输出一起打印堆栈跟踪记录。大多数时候，你会想看到这样的堆栈跟踪；然而，如果想在不同的机器和系统中比较程序的输出，抑制它是很有用的，因为堆栈跟踪显示的内存地址和库名/路径取决于机器的确切设置。
   * @see  异常情况
   *
   */
  void
  suppress_stacktrace_in_exceptions();

  /**
   * 当使用Assert()宏创建一个异常时，调用这个函数可以关闭对
   * <tt>std::abort()</tt>
   * 的使用。取而代之的是，异常将使用'throw'抛出，因此如果需要的话，它可以被捕获。一般来说，当Assert()被调用时，你想中止程序的执行，但如果你想记录所有创建的异常，或者你想测试一个断言是否正常工作，则需要关闭它。例如，在回归测试中就是这样做的。请注意，一些致命的错误仍然会调用abort()，例如，在异常处理过程中捕获的异常。
   * @see  异常现象
   *
   */
  void
  disable_abort_on_exception();

  /**
   * 这个命名空间中的函数与Assert和AssertThrow机制有关，但仅用于内部目的，不能在异常处理和抛出机制之外使用。
   * @ingroup Exceptions
   *
   */
  namespace internals
  {
    /**
     * 通过打印 @p exc 提供的错误信息并调用
     * <tt>std::abort()</tt>. 中止程序。
     *
     */
    [[noreturn]] void
    abort(const ExceptionBase &exc) noexcept;

    /**
     * 一个描述如何处理issue_error_noreturn中的异常的枚举。
     *
     */
    enum ExceptionHandling
    {
      /**
       * 通过调用 <code>std::abort</code> 中止程序，除非
       * deal_II_exceptions::disable_abort_on_exception
       * 已经被调用：在这种情况下，程序将抛出一个异常。
       *
       */
      abort_or_throw_on_exception,
      /**
       * 正常抛出异常。
       *
       */
      throw_on_exception
    };

    /**
     * 这个例程为<tt>Assert</tt>和<tt>AssertThrow</tt>宏中使用的异常生成机制做主要工作：正如其名称所暗示的，这个函数要么以抛出一个异常结束（如果
     * @p handling  ]是throw_on_exception，或者 @p handling
     * 是try_abort_exception且
     * deal_II_exceptions::disable_abort_on_exception
     * 是false），或者调用<tt>abort</tt>（如果 @p handling
     * 是try_abort_exception且
     * deal_II_exceptions::disable_abort_on_exception 是true）。
     * 实际的异常对象（最后一个参数）通常是一个在原地创建的未命名的对象；因为我们修改了它，所以我们不能通过常量引用来获取它，而且临时变量不会与非常量引用绑定。
     * 所以要用一个模板化的类型来取它的值（=copy
     * it），以避免切分
     *
     *
     *
     *
     *
     *
     * - 反正对性能的影响是相当小的。          @ref ExceptionBase
     *
     */
    template <class ExceptionType>
    [[noreturn]] void
    issue_error_noreturn(ExceptionHandling handling,
                         const char *      file,
                         int               line,
                         const char *      function,
                         const char *      cond,
                         const char *      exc_name,
                         ExceptionType     e) {
      // Fill the fields of the exception object
      e.set_fields(file, line, function, cond, exc_name);

      switch (handling)
        {
          case abort_or_throw_on_exception:
            {
              if (dealii::deal_II_exceptions::internals::
                    allow_abort_on_exception)
                internals::abort(e);
              else
                {
                  // We are not allowed to abort, so just throw the error:
                  throw e;
                }
            }
          case throw_on_exception:
            throw e;
          // this function should never return (and AssertNothrow can);
          // something must have gone wrong in the error handling code for us
          // to get this far, so throw an exception.
          default:
            throw ::dealii::StandardExceptions::ExcInternalError();
        }
    }

    /**
     * 内部函数，完成issue_error_nothrow的工作。
     *
     */
    void do_issue_error_nothrow(const ExceptionBase &e) noexcept;

    /**
     * 异常生成机制，以防我们必须不抛出。
     * @ref ExceptionBase
     * @note
     * 这个函数是用一个模板定义的，原因与issue_error_noreturn()相同。
     *
     */
    template <class ExceptionType>
    void
    issue_error_nothrow(const char *  file,
                        int           line,
                        const char *  function,
                        const char *  cond,
                        const char *  exc_name,
                        ExceptionType e) noexcept
    {
      static_assert(std::is_base_of<ExceptionBase, ExceptionType>::value,
                    "The provided exception must inherit from ExceptionBase.");
      // Fill the fields of the exception object
      e.set_fields(file, line, function, cond, exc_name);
      // avoid moving a bunch of code into the header by dispatching to
      // another function:
      do_issue_error_nothrow(e);
    }
#ifdef DEAL_II_WITH_CUDA
    /**
     * 返回一个给定错误代码的字符串。这类似于cudaGetErrorString函数，但在cuSPARSE中没有相应的函数。
     *
     */
    std::string
    get_cusparse_error_string(const cusparseStatus_t error_code);

    /**
     * 返回一个给定错误代码的字符串。这类似于cudaGetErrorString函数，但在cuSOLVER中没有相应的函数。
     *
     */
    std::string
    get_cusolver_error_string(const cusolverStatus_t error_code);
#endif
  }  /*namespace internals*/ 

}  /*namespace deal_II_exceptions*/ 



/**
 * 一个宏，作为调试模式错误检查的异常机制的主要程序。它断言某个条件得到满足，否则会发出错误并中止程序。
 * 更详细的描述可以在 @ref
 * Exceptions 模块中找到。它在  step-5  和  step-6
 * 中首次使用。更多信息也请参见<tt>ExceptionBase</tt>类。
 *
 *
 * @note  仅在DEBUG模式下活动
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，它们没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中包括
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#ifdef DEBUG
#  ifdef DEAL_II_HAVE_BUILTIN_EXPECT
#    define Assert(cond, exc)                                            \
      {                                                                  \
        if (__builtin_expect(!(cond), false))                            \
          ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
            ::dealii::deal_II_exceptions::internals::                    \
              abort_or_throw_on_exception,                               \
            __FILE__,                                                    \
            __LINE__,                                                    \
            __PRETTY_FUNCTION__,                                         \
            #cond,                                                       \
            #exc,                                                        \
            exc);                                                        \
      }
#  else  /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/ 
#    define Assert(cond, exc)                                            \
      {                                                                  \
        if (!(cond))                                                     \
          ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
            ::dealii::deal_II_exceptions::internals::                    \
              abort_or_throw_on_exception,                               \
            __FILE__,                                                    \
            __LINE__,                                                    \
            __PRETTY_FUNCTION__,                                         \
            #cond,                                                       \
            #exc,                                                        \
            exc);                                                        \
      }
#  endif  /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/ 
#else
#  define Assert(cond, exc) \
    {}
#endif



/**
 * 上面的<tt>Assert</tt>宏的一个变种，只要不调用disable_abort_on_exception，就会表现出相同的运行时行为。
 * 然而，如果disable_abort_on_exception被调用，这个宏只是打印出将被抛出的异常到deallog，并继续正常运行而不抛出异常。
 * 更详细的描述可以在 @ref
 * Exceptions
 * 模块中找到，在页面底部关于角落案例的讨论中。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，这些定义的前缀不是一个可能使它们对deal.II独一无二的字符串。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @note  仅在DEBUG模式下有效。
 *
 * @ingroup Exceptions
 *
 *
 */
#ifdef DEBUG
#  ifdef DEAL_II_HAVE_BUILTIN_EXPECT
#    define AssertNothrow(cond, exc)                                    \
      {                                                                 \
        if (__builtin_expect(!(cond), false))                           \
          ::dealii::deal_II_exceptions::internals::issue_error_nothrow( \
            __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc); \
      }
#  else  /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/ 
#    define AssertNothrow(cond, exc)                                    \
      {                                                                 \
        if (!(cond))                                                    \
          ::dealii::deal_II_exceptions::internals::issue_error_nothrow( \
            __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc); \
      }
#  endif  /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/ 
#else
#  define AssertNothrow(cond, exc) \
    {}
#endif

/**
 * 一个宏，作为动态错误检查的异常机制的主要程序。它断言某个条件得到满足，否则通过C++
 * @p throw  机制抛出一个异常。这个异常可以通过 @p catch
 * 子句来捕获，如 step-6 和以下所有教程程序所示。
 * 更详细的描述可以在 @ref
 * Exceptions 模块中找到。它在  step-9  和  step-13
 * 中首次使用。更多信息也请参见<tt>ExceptionBase</tt>类。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，这些定义没有用一个可能使它们对deal.II唯一的字符串作为前缀。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有用
 * <code>DEAL</code> or <code>deal</code>
 * 作为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @note  在DEBUG和RELEASE模式下都有效。
 *
 * @ingroup Exceptions
 *
 *
 */
#ifdef DEAL_II_HAVE_BUILTIN_EXPECT
#  define AssertThrow(cond, exc)                                       \
    {                                                                  \
      if (__builtin_expect(!(cond), false))                            \
        ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
          ::dealii::deal_II_exceptions::internals::throw_on_exception, \
          __FILE__,                                                    \
          __LINE__,                                                    \
          __PRETTY_FUNCTION__,                                         \
          #cond,                                                       \
          #exc,                                                        \
          exc);                                                        \
    }
#else  /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/ 
#  define AssertThrow(cond, exc)                                       \
    {                                                                  \
      if (!(cond))                                                     \
        ::dealii::deal_II_exceptions::internals::issue_error_noreturn( \
          ::dealii::deal_II_exceptions::internals::throw_on_exception, \
          __FILE__,                                                    \
          __LINE__,                                                    \
          __PRETTY_FUNCTION__,                                         \
          #cond,                                                       \
          #exc,                                                        \
          exc);                                                        \
    }
#endif  /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/ 

/**
 * 尺寸不匹配的特殊断言。
 * 因为这被经常使用，并且总是重复参数，所以我们为ExcDimensionMismatch引入了这个特殊的断言，以保持用户代码的简短。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，这些定义没有用一个可能使它们对deal.II唯一的字符串作前缀。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，只要在所有其他deal.II的头文件中包含
 * <code>deal.II/base/undefine_macros.h</code> ，就可以不以
 * <code>DEAL</code> or <code>deal</code> 为前缀了。
 *
 *
 * @ingroup Exceptions
 *
 */
#define AssertDimension(dim1, dim2)                                            \
  Assert(static_cast<typename ::dealii::internal::argument_type<void(          \
             typename std::common_type<decltype(dim1),                         \
                                       decltype(dim2)>::type)>::type>(dim1) == \
           static_cast<typename ::dealii::internal::argument_type<void(        \
             typename std::common_type<decltype(dim1),                         \
                                       decltype(dim2)>::type)>::type>(dim2),   \
         dealii::ExcDimensionMismatch((dim1), (dim2)))


/**
 * 一个测试<tt>vec</tt>是否有大小<tt>dim1</tt>的断言，并且向量的每个条目本身是一个大小为<tt>dim2</tt>的数组。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，这些定义的前缀不是一个字符串，而这个字符串很可能使它们成为deal.II的唯一。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有前缀
 * <code>DEAL</code> or <code>deal</code>
 * ，可以在所有其他deal.II的头文件中加入头文件
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#define AssertVectorVectorDimension(VEC, DIM1, DIM2) \
  AssertDimension(VEC.size(), DIM1);                 \
  for (const auto &subvector : VEC)                  \
    {                                                \
      (void)subvector;                               \
      AssertDimension(subvector.size(), DIM2);       \
    }

namespace internal
{
  // Workaround to allow for commas in template parameter lists
  // in preprocessor macros as found in
  // https://stackoverflow.com/questions/13842468/comma-in-c-c-macro
  template <typename T>
  struct argument_type;
  template <typename T, typename U>
  struct argument_type<T(U)>
  {
    using type = U;
  };
} // namespace internal

/**
 * 一个测试给定索引是否在半开范围内的断言
 * <code>[0,range)</code>
 * 。如果断言失败，它会抛出一个异常对象
 * <code>ExcIndexRange(index,0,range)</code>  。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，这些定义的前缀不是一个可能使它们在deal.II中独一无二的字符串。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 */
#define AssertIndexRange(index, range)                                         \
  Assert(                                                                      \
    static_cast<typename ::dealii::internal::argument_type<void(               \
        typename std::common_type<decltype(index), decltype(range)>::type)>::  \
                  type>(index) <                                               \
      static_cast<typename ::dealii::internal::argument_type<void(             \
        typename std::common_type<decltype(index), decltype(range)>::type)>::  \
                    type>(range),                                              \
    dealii::ExcIndexRangeType<typename ::dealii::internal::argument_type<void( \
      typename std::common_type<decltype(index), decltype(range)>::type)>::    \
                                type>((index), 0, (range)))

/**
 * 一个检查一个数字是否是有限的断言。我们明确地将这个数投到
 * std::complex 中，以匹配异常的签名（关于为什么使用
 * std::complex 的解释，见那里），并满足 std::complex
 * 没有隐含转换的事实。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，这些定义的前缀不是一个可能使它们对deal.II唯一的字符串。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#define AssertIsFinite(number)               \
  Assert(dealii::numbers::is_finite(number), \
         dealii::ExcNumberNotFinite(std::complex<double>(number)))

#ifdef DEAL_II_WITH_MPI
/**
 * 一个断言，检查MPI函数返回的错误代码是否等于
 * <code>MPI_SUCCESS</code>
 * 。如果检查失败，那么就会抛出一个ExcMPI类型的异常，给定的错误代码作为参数。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，这些定义的前缀不是一个可能使它们对deal.II独一无二的字符串。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有
 * <code>DEAL</code> or <code>deal</code>
 * 的前缀，可以在所有其他deal.II的头文件中包括
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @note  只有在deal.II被编译为MPI时才有效。
 *
 * @ingroup Exceptions
 *
 */
#  define AssertThrowMPI(error_code) \
    AssertThrow(error_code == MPI_SUCCESS, dealii::ExcMPI(error_code))
#else
#  define AssertThrowMPI(error_code) \
    {}
#endif // DEAL_II_WITH_MPI

#ifdef DEAL_II_WITH_CUDA
/**
 * 一个断言，检查调用CUDA例程产生的错误代码是否等于cudaSuccess。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，这些定义没有用一个字符串作为前缀，而这个字符串很可能使它们成为deal.II的唯一。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  ifdef DEBUG
#    define AssertCuda(error_code)      \
      Assert(error_code == cudaSuccess, \
             dealii::ExcCudaError(cudaGetErrorString(error_code)))
#  else
#    define AssertCuda(error_code) \
      {                            \
        (void)(error_code);        \
      }
#  endif

/**
 * AssertCuda的非抛出式等价物。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，这些定义没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果不以
 * <code>DEAL</code> or <code>deal</code>
 * 为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 */
#  ifdef DEBUG
#    define AssertNothrowCuda(error_code)      \
      AssertNothrow(error_code == cudaSuccess, \
                    dealii::ExcCudaError(cudaGetErrorString(error_code)))
#  else
#    define AssertNothrowCuda(error_code) \
      {                                   \
        (void)(error_code);               \
      }
#  endif

/**
 * 一个断言，检查内核是否被成功启动和执行。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，这些定义没有用一个字符串作为前缀，这很可能使它们对deal.II来说是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有
 * <code>DEAL</code> or <code>deal</code>
 * 的前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  ifdef DEBUG
#    define AssertCudaKernel()                                \
      {                                                       \
        cudaError_t local_error_code = cudaPeekAtLastError(); \
        AssertCuda(local_error_code);                         \
        local_error_code = cudaDeviceSynchronize();           \
        AssertCuda(local_error_code)                          \
      }
#  else
#    define AssertCudaKernel() \
      {}
#  endif

/**
 * 一个断言，检查调用cuSPARSE例程产生的错误代码是否等于CUSPARSE_STATUS_SUCCESS。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，这些定义的前缀不是一个可能使它们对deal.II唯一的字符串。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果它们的前缀不是
 * <code>DEAL</code> or <code>deal</code>
 * ，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  ifdef DEBUG
#    define AssertCusparse(error_code)                                      \
      Assert(                                                               \
        error_code == CUSPARSE_STATUS_SUCCESS,                              \
        dealii::ExcCusparseError(                                           \
          dealii::deal_II_exceptions::internals::get_cusparse_error_string( \
            error_code)))
#  else
#    define AssertCusparse(error_code) \
      {                                \
        (void)(error_code);            \
      }
#  endif

/**
 * AssertCusparse的非抛出式等价物。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理程序定义的例子，这些定义没有用一个字符串作为前缀，这很可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有用
 * <code>DEAL</code> or <code>deal</code>
 * 作为前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  ifdef DEBUG
#    define AssertNothrowCusparse(error_code)                               \
      AssertNothrow(                                                        \
        error_code == CUSPARSE_STATUS_SUCCESS,                              \
        dealii::ExcCusparseError(                                           \
          dealii::deal_II_exceptions::internals::get_cusparse_error_string( \
            error_code)))
#  else
#    define AssertNothrowCusparse(error_code) \
      {                                       \
        (void)(error_code);                   \
      }
#  endif

/**
 * 一个断言，检查调用cuSOLVER例程产生的错误代码是否等于CUSOLVER_STATUS_SUCCESS。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，它没有用一个字符串作为前缀，而这个字符串可能使它们在deal.II中是独一无二的。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有
 * <code>DEAL</code> or <code>deal</code>
 * 的前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  ifdef DEBUG
#    define AssertCusolver(error_code)                                      \
      Assert(                                                               \
        error_code == CUSOLVER_STATUS_SUCCESS,                              \
        dealii::ExcCusparseError(                                           \
          dealii::deal_II_exceptions::internals::get_cusolver_error_string( \
            error_code)))
#  else
#    define AssertCusolver(error_code) \
      {                                \
        (void)(error_code);            \
      }
#  endif

#endif

#ifdef DEAL_II_TRILINOS_WITH_SEACAS
/**
 * 断言，检查调用ExodusII例程产生的错误代码是否等于EX_NOERR（为零）。
 *
 *
 * @note
 * 这个和类似的宏名称是deal.II库中预处理器定义的例子，这些定义没有用一个可能使它们对deal.II唯一的字符串作前缀。因此，你的代码与其他库的接口有可能定义相同的名称，其结果将是名称冲突（见https://en.wikipedia.org/wiki/Name_collision）。我们可以
 * <code>\#undef</code>
 * 这个宏，以及所有其他由deal.II定义的宏，如果没有
 * <code>DEAL</code> or <code>deal</code>
 * 的前缀，可以在所有其他deal.II的头文件中加入
 * <code>deal.II/base/undefine_macros.h</code> 。
 *
 *
 * @ingroup Exceptions
 *
 *
 */
#  define AssertThrowExodusII(error_code) \
    AssertThrow(error_code == 0, ExcExodusII(error_code));
#endif // DEAL_II_TRILINOS_WITH_SEACAS

using namespace StandardExceptions;

DEAL_II_NAMESPACE_CLOSE

#endif


