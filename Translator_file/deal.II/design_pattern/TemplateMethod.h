#include <iostream>

/**
 * Template Method Design Pattern
 *
 * 意图: 在超类中定义了算法的骨架，但允许子类在不改变算法结构的情况下覆盖算法的特定步骤。
 */

namespace DesignPattern
{

    namespace TemplateMethod
    {
        /**
 * 抽象类定义了一个模板方法，它包含了某种算法的骨架，由对（通常）抽象的原始操作的调用组成。
 * 具体的子类应该实现这些操作，但要保持模板方法本身的完整性。
 */
        class AbstractClass
        {
            /**
   * 模板方法定义了一个算法的骨架。
   */
        public:
            void TemplateMethod() const
            {
                this->BaseOperation1();
                this->RequiredOperations1();
                this->BaseOperation2();
                this->Hook1();
                this->RequiredOperation2();
                this->BaseOperation3();
                this->Hook2();
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
            virtual void RequiredOperations1() const = 0;
            virtual void RequiredOperation2() const = 0;
            /**
   *  这些是 "钩子"。子类可以覆盖它们，但这不是强制性的，因为这些钩子已经有默认的（但空的）实现。钩子在算法的一些关键地方提供了额外的扩展点。
   */
            virtual void Hook1() const {}
            virtual void Hook2() const {}
        };
        /**
 * 具体类必须实现基类的所有抽象操作。他们也可以用一个默认的实现来覆盖一些操作。
 */
        class ConcreteClass1 : public AbstractClass
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
        void ClientCode(AbstractClass *class_)
        {
            // ...
            class_->TemplateMethod();
            // ...
        }

    }
}

// int main() {
//   std::cout << "Same client code can work with different subclasses:\n";
//   ConcreteClass1 *concreteClass1 = new ConcreteClass1;
//   ClientCode(concreteClass1);
//   std::cout << "\n";
//   std::cout << "Same client code can work with different subclasses:\n";
//   ConcreteClass2 *concreteClass2 = new ConcreteClass2;
//   ClientCode(concreteClass2);
//   delete concreteClass1;
//   delete concreteClass2;
//   return 0;
// }
