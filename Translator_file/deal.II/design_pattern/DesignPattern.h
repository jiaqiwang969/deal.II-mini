
/**
 * 参考：李建忠老师C++设计模式教程
 * 
 * https://cdn.mathpix.com/snip/images/UtopqRmgFLWYUjUz_Bex-aRnuLh28yyHUA1WEqcnx1c.original.fullsize.png
 *
 * # 基本原则 -> 面对变化，提高复用（二进制级，非源代码级） -> 寻找并划分变化点，在合适的时机使用设计模式 -> 重构到模式
 * - 原则1: 依赖倒置 (DIP)
 * -# 高层模块(稳定)不应该依赖于低层模块(变化)，二者都应该依赖
 * 于抽象(稳定)。 
 * 
 * -# 抽象(稳定)不应该依赖于实现细节(变化)，实现细节应该依赖于
 * 抽象(稳定)。
 *
 * - 原则2: 开放封闭原则 (OCP)
 * -# 对扩展开放，对更改封闭。 
 * 
 * -# 类模块应该是可扩展的，但是不可修改。
 *
 * - 原则3: 单一职责原则 (SRP)
 * -# 一个类应该仅有一个引起它变化的原因。 
 * 
 * -# 变化的方向隐含着类的责任。
 *
 * - 原则4: Liskov 替换原则 (LSP)
 * -# 子类必须能够替换它们的基类(IS-A)。 
 * 
 * -# 继承表达类型抽象。
 *
 * - 原则5: 接口隔离原则 (ISP)
 * -# 不应该强迫客户程序依赖它们不用的方法。 
 * 
 * -# 接口应该小而完备。
 *
 * - 原则6: 优先使用对象组合，而不是类继承
 * -# 类继承通常为 "白箱复用" ，对象组合通常为 "黑箱复用" 。 
 * 
 * -# 继承在某种程度上破坏了封装性，子类父类耦合度高。 
 * 
 * -# 而对象组合则只要求被组合的对象具有良好定义的接口，耦合度低。
 */

/**
 * - 原则7: 封装变化点
 * -# 使用封装来创建对象之间的分界层，让设计者可以在分界层的一侧进行修改，而不会对另一侧产生不良的影响，从而实现层次间的松耦合。
 * 
 * -# 使用封装来创建对象之间的分界层，让设计者可以在分界层的一侧进行修改，而不会对另一侧产生不良的影响，从而实现层次间的松耦合。
 */

/**
 * - 原则8: 针对接口编程，而不是针对实现编程
 * -# 不将变量类型声明为某个特定的具体类，而是声明为某个接口。
 * 
 * -# 客户程序无需获知对象的具体类型，只需要知道对象所具有的接口。
 * 
 * -# 减少系统中各部分的依赖关系，从而实现 "高内聚、松耦合" 的类型设计方案。
 */

/**
 * # 分类
 *
 * 从目的来看: 
 *   - 创建型 ( Creational ) 模式 ：将对象的部分创建工作延迟到子类或者其他对象，从市应对需求变化为对象仓建时具体类型实现引来的冲击。
 *
 *   - 结构型 ( Structural ) 模式 ：通过类继承或者对象组合获得更灵活的结构，从而应对需求变化为对象的结构带来的冲击。
 *
 *   - 行为型 ( Behavioral ) 模式 ：通过类继承或者对象组合来划分类与对象间的职责，从而应对需求变化为多个交互的对象带来的冲击。
 *
 *
 * 从范围来看 : 
 *   - 类模式：处理类与子类的静态关系。
 *
 *   - 对象模式：处理对象的动态关系。
 */

/**
 * # 从封装变化角度对模式分类
 * - 组件协作: 晚期绑定
 * Template Method Observer / Event Strategy
 *
 * - 对象性能:
 * Singleton Flyweight
 *
 * - 数据结构:
 * Composite Iterator Chain of Resposibility
 *
 * - 单一职责:
 * Decorator Bridge
 *
 * - 接口隔离:
 * Facade Proxy Mediator Adapter
 *
 * - 行为变化:
 * Command Visitor
 *
 * - 对象创建:
 * Factory Method Abstract Factory Prototype Builder
 *
 * - 状态变化:
 * Memento State
 *
 * - 领域问题:
 * Interpreter
 */

/**
 * # 重构关键技法：
 * - 静态 -> 动态
 * - 早绑定 -> 晚绑定
 * - 继承 -> 组合
 * - 编译时依赖 -> 运行时依赖
 * - 紧耦合 -> 松耦合
 */

/**
 * 
 *# 代码的坏味道 （马丁福勒）
 *- 神秘命名(Mysterious Name)
 *    -# 整洁代码最重要的一环就是用好的名字
 *        改变函数声明
 *        字段改名
 *        变量改名
 *    -# 如果你想不出一个好名字，说明背后很可能潜藏这更深的设计问题
 *- 重复代码(Duplicated Code)
 *    -# 如果你在一个以上的地点看到相同的代码结构，那么可以肯定：设法将它们合二为一，程序会变得更好
 *        提炼函数
 *        移动语句
 *        函数上移
 *- 过长函数(Long Function)
 *    -# 关键不在于函数的长度，而在于函数“做什么”和“如何做”之间的语义距离
 *    -# 如果函数有大量的参数和临时变量
 *        以查询取代临时变量
 *        引入参数对象
 *        保持对象完整
 *    -# 寻找注释，它们通常能指出代码用途和实现手法之间的语义距离
 *        提炼函数
 *    -# 条件表达式和循环常常也是提炼的信号
 *        分解条件表达式
 *        以多态取代条件表达式
 *- 全局数据(Global Data)
 *    -# 全局数据的问题在于，从代码库的任何一个角落都可以修改它，而且没有任何机制可以探测出到底哪段代码做出了修改
 *        封装变量
 *- 可变数据(Mutable Data)
 *    -# 在一处更新数据，却没有意识到软件的另一处获取的数据也发生了改变，功能失效了
 *    -# 数据永不改变——函数式编程
 *    -# 确保所有数据更新操作通过某几个函数进行
 *        封装变量
 *    -# 将某个多用途的变量拆分，避免危险的更新操作
 *        拆分变量
 *    -# 将没有副作用的代码与执行数据更新操作的代码分开
 *        移动语句
 *        提炼函数
 *- 设计良好的API
 *    -# 将查询函数和修改函数分离
 *    -# 移除设值函数
 *- 过长参数列表(Long Parameter List)
 *    -# 函数的参数应该总结该函数的可变性
 *        保持对象完整
 *        引入参数对象
 *        移除标记参数
 *        函数组合成类
 *- 冗赘的元素(Lazy Element)
 *    -# 这种程序元素（如类和函数）的结构与其实现代码相比不再有额外的意义
 *        内联函数
 *        内联类
 *        折叠继承体系
 *- 夸夸其谈通用性(Speculative Generality)
 *    -# 当前没有被用到，期待未来可能用到的复杂设计，应该移除它
 *        内联函数
 *        内联类
 *        改变函数声明
 *        折叠继承体系
 *- 临时字段(Temporary Field)
 *    -# 某个字段就像孤儿一样只为某种临时目的存在
 *        提炼类
 *        搬移函数
 *        引入特例
 *- 过长的消息链(Message Chains)
 *    -# 中间关系太多
 *        隐藏委托关系
 *        提炼函数&搬移函数
 *- 中间人(Middle Man)
 *    -# 过度引用委托
 *        移除中间人
 *        内联函数
 *        以委托取代超类/以委托取代子类
 *- 发散式变化(Divergent Change)
 *    -# 每次应该只关心一个上下文
 *    -# 两个不同的模块其上下文发生了重叠
 *        拆分阶段
 *        搬移函数
 *        提炼函数
 *        提炼类
 *- 散弹式修改(Shotgun Surgery)
 *    -# 某个模块的上下文四散各处
 *        搬移函数
 *        搬移字段
 *        函数组合成类
 *        函数组合成变换
 *        拆分阶段
 *        内联函数
 *        内联类
 *- 依恋情节(Feature Envy)
 *    -# 所谓模块化，就是力求将代码分出区域，最大化区域内部的交互，最小化跨区域的交互
 *    -# 一个函数跟另一个模块中的东西打交道特别多
 *        搬移函数
 *        提炼函数
 *- 数据泥团(Data Clumps)
 *    -# 数据成群结队，出现在了不该在的位置
 *        提炼类
 *        引入参数对象
 *        保持对象完整
 *    -# 删掉某个数据，如果这么做，其他数据也不再有意义，那么应该为它们产生一个新对象。
 *- 基本类型偏执(Primitive Obsession)
 *    -# 只用字符串等并不能表示某种特定数据结构，它应该有一个专门的类型
 *        以对象取代基本类型
 *        以子类取代类型码
 *        以多态取代条件表达式
 *        提炼类
 *        引入参数对象
 *- 重复的switch(Repeated Switches)
 *    -# 在不同地方反复使用同样的switch逻辑
 *        以多态取代条件表达式
 *- 循环语句(Loops)
 *    -# 管道操作（如filter和map）可以帮助我们更快地看清被处理的元素以及处理它们的动作
 *        以管道取代循环
 *- 内幕交易(Insider Trading)
 *    -# 两个模块总是在窃窃私语
 *        搬移函数&搬移字段
 *        隐藏委托关系
 *        以委托取代子类/以委托取代超类
 *- 过大的类(Large Class)
 *    -# 单个类做了太多事情
 *        提炼类
 *    -# 类中的数个变量有着相同的前缀或后缀
 *        提炼超类/以子类取代类型码
 *    -# 异曲同工的类(Alternative Classes with Different Interfaces)
 *        改变函数声明
 *        搬移函数&提炼超类
 *- 纯数据类(Data Class)
 *    -# 数据一定和行为相关联，纯数据类常常意味着行为被放在了错误的地方
 *        封装记录
 *        移除设值函数
 *        搬移函数
 *        提炼函数
 *        拆分阶段
 *- 被拒绝的遗赠(Refused Bequest)
 *    -# 子类不需要所有的继承
 *        以委托取代子类
 *        以委托取代超类
 *- 注释(Comments)
 *    -# 当你感觉需要撰写注释时，请先尝试重构
 *        提炼函数
 *        改变函数声明
 *        引入断言
 *    -# 并不是说不应该写注释！
 *
 *
 *
 *
 *# 重构手法
 *- 第一组重构
 *    -# 提炼函数(Extract Function)
 *    -# 内联函数(Inline Function)
 *    -# 提炼变量(Extract Variable)
 *    -# 内联变量(Inline Variable)
 *    -# 改变函数声明(Change Function Declaration)
 *    -# 封装变量(Encapsulate Variable)
 *    -# 引入参数对象(Introduce Parameter Object)
 *    -# 函数组合成类(Combine Functions into Class)
 *    -# 函数组合成变换(Combine Functions into Transform)
 *    -# 拆分阶段(Split Phase)
 *- 封装
 *    -# 封装记录(Encapsulate Record)
 *    -# 封装集合(Encapsulate Collection)
 *    -# 以对象取代基本类型(Replace Primitive with Object)
 *    -# 以查询取代临时变量(Replace Temp with Query)
 *    -# 提炼类(Extract Class)
 *    -# 内联类(Inline Class)
 *    -# 隐藏委托关系(Hide Delegate)
 *    -# 移除中间人(Remove Middle Man)
 *    -# 替换算法(Substitute Algorithm)
 *- 处理继承关系
 *    -# 函数上移(Pull Up Method)
 *    -# 字段上移(Pull Up Field)
 *    -# 构造函数本体上移(Pull Up Constructor Body)
 *    -# 函数下移(Push Down Method)
 *    -# 字段下移(Push Down Field)
 *    -# 以子类取代类型码(Replace Type Code with Subclasses)
 *    -# 移除子类(Remove Subclass)
 *    -# 提炼超类(Extract Superclass)
 *    -# 折叠继承体系(Collapse Hierarchy)
 *    -# 以委托取代子类(Replace Subclass with Delegate)
 *    -# 以委托取代超类(Replace Superclass with Delegate)
 *- 搬移特性
 *    -# 搬移函数(Move Function)
 *    -# 搬移字段(Move Field)
 *    -# 搬移语句到函数(Move Statements into Function)
 *    -# 搬移语句到调用者(Move Satements to Callers)
 *    -# 以函数调用取代内联代码(Replace Inline Code with Function Call)
 *    -# 移动语句(Slide Statements)
 *    -# 拆分循环(Split Loop)
 *    -# 以管道取代循环(Replace Loop with Pipeline)
 *    -# 移除死代码(Remove Dead Code)
 *- 重新组织数据
 *    -# 拆分变量(Split Variable)
 *    -# 字段改名(Rename Field)
 *    -# 以查询取代派生变量(Replace Derived Variable with Query)
 *    -# 将引用对象改为值对象(Change Reference to Value)
 *    -# 将值对象改为引用对象(Change Value to Reference)
 *- 简化条件逻辑
 *    -# 分解条件表达式(Decompose Conditional)
 *    -# 合并条件表达式(Consolidate Conditional Expression)
 *    -# 以卫语句取代嵌套条件表达式(Replace Nested Conditional with Guard Clauses)
 *    -# 以多态取代条件表达式(Replace Conditional with Polymorphsim)
 *    -# 引入特例(Introduce Special Case)
 *    -# 引入断言(Introduce Assertion)
 *- 重构API
 *    -# 将查询函数和修改函数分离(Separate Query from Modifier)
 *    -# 函数参数化(Parameterize Function)
 *    -# 移除标记参数(Remove Flag Argument)
 *    -# 保持对象完整(Preserve Whole Object)
 *    -# 以查询取代参数(Replace Parameter with Query)
 *    -# 以参数取代查询(Replace Query with Parameter)
 *    -# 移除设值函数(Remove Setting Method)
 *    -# 以工厂函数取代构造函数(Replace Constructor with Factory Function)
 *    -# 以命令取代函数(Replace Function with Command)
 *    -# 以函数取代命令(Replace Command with Function)
 *
 */

namespace DesignPattern
{

}
