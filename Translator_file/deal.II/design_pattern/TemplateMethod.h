#include <iostream>

namespace DesignPattern
{
    /**
     * TemplateMethod < 组件协作
     *
     * - 模式动机: 稳定的操作结构，子步骤变化
     *
     * 在超类中定义了算法的骨架，但允许子类在不改变算法结构的情况下覆盖算法的特定步骤。
     *
     * 1. 在软件构建过程中，对于某一项任务，它常常有稳定的整体操作 结构，但各个子步骤却有很多改变的需求，或者由于固有的原因 ( 比如框架与应用之间的关系）而无法和任务的整体结构同时实现。
     *
     * 2. 如何在确定稳定操作结构的前提下，来灵活应对各个子步骤的变 化或者晚期实现需求 ?
     */

    /**
     * - 模式定义:
     *
     * 定义一个操作中的算法的骨架
     * (稳定)，而将一些步骤延迟 (变化)到子类中。 Template
     * Method使得子类可以不改变
     * (复用)一个算法的结构即可重定义(override
     * 重写)该算法的某些特定步骤。 --《设计模式》 GOF
     *
     * 具体的子类应该实现这些操作，但要保持模板方法本身的完整性。
     *
     * 不要调用我，让我来调用你！需要虚函数晚绑定的机制来调用你。虚函数背后就是虚函数表上挂了一个函数指针。
     */

    /**
     * 设计FEM库：这种方法在各种算法在某些步骤上有差异但在全局上没有差异的情况下很有用。可以用这种模式设计策略，提供一类可由方案改变的算法，或者将这种模式应用于线性求解器的设计，可以使它们独立于矩阵和向量及其操作。 by Pooyan Dadvand
     * 
     */

    /**
     * 个人感悟：
     * 
     * 这个模式可以用来对大型的工程问题做变量分析，专门通过虚函数将变化的部分，隔离出来，每次运行，稳定部分的编译不会发生变化。
     * 
     * 非常适合做 不同方法 对整体效果的评估，比如论文里面，某些方法的效果对比。
     * 
     * 对于大厂开发来说，应用人员被隔离了，他们看不到也偷不走核心技术。
     */

    namespace TemplateMethod
    {

        class AbstractClass
        {
            /**
             * 模板方法定义了一个算法的骨架。
             */
        public:
            void TemplateMethod() const
            {
                this->BaseOperation1();      // defined by Lib
                this->RequiredOperations1(); // defined by App
                this->BaseOperation2();      // defined by Lib
                this->Hook1();               // defined by App
                this->RequiredOperation2();  // defined by App
                this->BaseOperation3();      // defined by Lib
                this->Hook2();               // defined by App
            }
            /**
             * 这些操作已经有了实现。
             */
        protected:
            void BaseOperation1() const
            {
                std::cout << "AbstractClass says: I am doing the bulk of the work\n";
            }
            void BaseOperation2() const
            {
                std::cout << "AbstractClass says: But I let subclasses override some operations\n";
            }
            void BaseOperation3() const
            {
                std::cout << "AbstractClass says: But I am doing the bulk of the work anyway\n";
            }
            /**
             * 这些操作必须在子类中实现。
             */
            virtual void RequiredOperations1() const = 0; // defined by App
            virtual void RequiredOperation2() const = 0;  // defined by App
            /**
             * 这些是
             * "钩子"。子类可以覆盖它们，但这不是强制性的，因为这些钩子已经有默认的（但空的）实现。钩子在算法的一些关键地方提供了额外的扩展点。
             */
            virtual void Hook1() const {} // defined by App
            virtual void Hook2() const {} // defined by App
        };
        /**
         * 具体类必须实现基类的所有抽象操作。他们也可以用一个默认的实现来覆盖一些操作。
         */
        class ConcreteClass1 : public
        {
        protected:
            void RequiredOperations1() const override
            {
                std::cout << "ConcreteClass1 says: Implemented Operation1\n";
            }
            void RequiredOperation2() const override
            {
                std::cout << "ConcreteClass1 says: Implemented Operation2\n";
            }
        };
        /**
         * 通常情况下，具体类只覆盖基类的一部分操作。
         */
        class ConcreteClass2 : public AbstractClass
        {
        protected:
            void RequiredOperations1() const override
            {
                std::cout << "ConcreteClass2 says: Implemented Operation1\n";
            }
            void RequiredOperation2() const override
            {
                std::cout << "ConcreteClass2 says: Implemented Operation2\n";
            }
            void Hook1() const override
            {
                std::cout << "ConcreteClass2 says: Overridden Hook1\n";
            }
        };
        /**
         * 客户端代码调用模板方法来执行该算法。客户端代码不必知道它所处理的对象的具体类别，只要它通过其基类的接口与对象一起工作即可。
         */
        void ClientCode(AbstractClass *class_) // defined by Lib
        {
            // ...
            class_->TemplateMethod();
            // ...
        }

    }

}

// defined by App

/**
 * 
 * int main() {
 *   std::cout << "Same client code can work with different subclasses:\n";
 *   ConcreteClass1 *concreteClass1 = new ConcreteClass1;
 *   ClientCode(concreteClass1);
 *   std::cout << "\n";
 *   std::cout << "Same client code can work with different subclasses:\n";
 *   ConcreteClass2 *concreteClass2 = new ConcreteClass2;
 *   ClientCode(concreteClass2);
 *   delete concreteClass1;
 *   delete concreteClass2;
 *   return 0;
 * }
 * 
 */
