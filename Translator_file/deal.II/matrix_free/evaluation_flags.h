//include/deal.II-translator/matrix_free/evaluation_flags_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_flags_h
#define dealii_matrix_free_evaluation_flags_h

#include <deal.II/base/config.h>


DEAL_II_NAMESPACE_OPEN



/**
 *
 * @brief The namespace for the EvaluationFlags enum
 * 这个命名空间包含了FEEvaluation中使用的EvaluationFlags枚举，用来控制数值、梯度等的评估和整合。
 *
 *
 */
namespace EvaluationFlags
{
  /**
   * @brief The EvaluationFlags enum  这个枚举包含一组由
   * FEEvaluation::integrate(),   FEEvaluation::evaluate()
   * 和其他方面使用的标志，以确定是否正在使用数值、梯度、豫度或它们的组合。
   *
   */
  enum EvaluationFlags
  {
    /**
     * 不要使用或计算任何东西。
     *
     */
    nothing = 0,
    /**
     * 使用或评估数值。
     *
     */
    values = 0x1,
    /**
     * 使用或评估梯度。
     *
     */
    gradients = 0x2,
    /**
     * 使用或评估赫西恩。
     *
     */
    hessians = 0x4
  };


  /**
   * 全局运算符，它返回一个对象，其中所有的位都被设置为第一个或第二个参数中的位。这个操作符的存在是因为如果它不存在，那么bit-or <tt>操作符|</tt>的结果将是一个整数，当我们试图将它赋值给UpdateFlags类型的对象时，又会引发编译器警告。
   * @ref EvaluationFlags
   *
   */
  inline EvaluationFlags
  operator|(const EvaluationFlags f1, const EvaluationFlags f2)
  {
    return static_cast<EvaluationFlags>(static_cast<unsigned int>(f1) |
                                        static_cast<unsigned int>(f2));
  }



  /**
   * 全局操作符，它将第二个参数的位也设置在第一个参数中。
   * @ref EvaluationFlags
   *
   */
  inline EvaluationFlags &
  operator|=(EvaluationFlags &f1, const EvaluationFlags f2)
  {
    f1 = f1 | f2;
    return f1;
  }


  /**
   * 全局操作符，它返回一个对象，其中所有位都被设置在第一个和第二个参数中。这个操作符的存在是因为如果它不存在，那么位和<tt>操作符&</tt>的结果将是一个整数，当我们试图将其分配给UpdateFlags类型的对象时，会引发编译器警告。
   * @ref EvaluationFlags
   *
   */
  inline EvaluationFlags operator&(const EvaluationFlags f1,
                                   const EvaluationFlags f2)
  {
    return static_cast<EvaluationFlags>(static_cast<unsigned int>(f1) &
                                        static_cast<unsigned int>(f2));
  }


  /**
   * 全局操作符，如果第一个参数中的所有位没有在第二个参数中设置，则将其清除。
   * @ref EvaluationFlags
   *
   */
  inline EvaluationFlags &
  operator&=(EvaluationFlags &f1, const EvaluationFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }

} // namespace EvaluationFlags


DEAL_II_NAMESPACE_CLOSE

#endif


