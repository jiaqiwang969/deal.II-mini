#include <array>
#include <iostream>
#include <string>

/**
  * 
  * Visitor
  *
  * - 模式动机: 让你把算法与它们所操作的对象分开。
  */

/**
  * - 模式定义:
  * 
  *   访客接口声明了一组与组件类对应的访问方法。一个访问方法的签名允许访问者识别它所处理的组件的确切类别。
  */
namespace DesignPattern
{
    namespace BehavioralDesign
    {
        namespace Visitor
        {

            class ConcreteComponentA;
            class ConcreteComponentB;

            class Visitor
            {
            public:
                virtual void VisitConcreteComponentA(const ConcreteComponentA *element) const = 0;
                virtual void VisitConcreteComponentB(const ConcreteComponentB *element) const = 0;
            };

            /**
         * 组件接口声明了一个`accept'方法，该方法应将基础访问者接口作为参数。
         */

            class Component
            {
            public:
                virtual ~Component() {}
                virtual void Accept(Visitor *visitor) const = 0;
            };

            /**
         * 每个具体的组件都必须以这样的方式实现`接受'方法，即调用与组件类相对应的访问者的方法。
         */
            class ConcreteComponentA : public Component
            {
                /**
             * 注意，我们正在调用`visitConcreteComponentA`，它与当前的类名相匹配。这样，我们让访问者知道它所工作的组件的类。
             */
            public:
                void Accept(Visitor *visitor) const override
                {
                    visitor->VisitConcreteComponentA(this);
                }
                /**
             * 具体组件可能有一些特殊的方法，这些方法在其基类或接口中并不存在。访问者仍然能够使用这些方法，因为它知道该组件的具体类。
             */
                std::string ExclusiveMethodOfConcreteComponentA() const
                {
                    return "A";
                }
            };

            class ConcreteComponentB : public Component
            {
                /**
             * Same here: visitConcreteComponentB => ConcreteComponentB
             */
            public:
                void Accept(Visitor *visitor) const override
                {
                    visitor->VisitConcreteComponentB(this);
                }
                std::string SpecialMethodOfConcreteComponentB() const
                {
                    return "B";
                }
            };

            /**
         * 具体的Visitor实现了同一算法的多个版本，它可以与所有的具体组件类一起工作。
         * 在与复杂的对象结构（如复合树）一起使用时，你可以体验到访问者模式的最大好处。在这种情况下，当在结构的各个对象上执行访问者的方法时，存储一些算法的中间状态可能会有帮助。
         */
            class ConcreteVisitor1 : public Visitor
            {
            public:
                void VisitConcreteComponentA(const ConcreteComponentA *element) const override
                {
                    std::cout << element->ExclusiveMethodOfConcreteComponentA() << " + ConcreteVisitor1\n";
                }

                void VisitConcreteComponentB(const ConcreteComponentB *element) const override
                {
                    std::cout << element->SpecialMethodOfConcreteComponentB() << " + ConcreteVisitor1\n";
                }
            };

            class ConcreteVisitor2 : public Visitor
            {
            public:
                void VisitConcreteComponentA(const ConcreteComponentA *element) const override
                {
                    std::cout << element->ExclusiveMethodOfConcreteComponentA() << " + ConcreteVisitor2\n";
                }
                void VisitConcreteComponentB(const ConcreteComponentB *element) const override
                {
                    std::cout << element->SpecialMethodOfConcreteComponentB() << " + ConcreteVisitor2\n";
                }
            };
            /**
         * 客户端代码可以在任何元素的集合上运行访问者操作，而不需要弄清楚它们的具体类别。接受操作会引导对访问者对象中适当的操作的调用。
         */
            void ClientCode(std::array<const Component *, 2> components, Visitor *visitor)
            {
                // ...
                for (const Component *comp : components)
                {
                    comp->Accept(visitor);
                }
                // ...
            }

        }
    }
}
// int main() {
//   std::array<const Component *, 2> components = {new ConcreteComponentA, new ConcreteComponentB};
//   std::cout << "The client code works with all visitors via the base Visitor interface:\n";
//   ConcreteVisitor1 *visitor1 = new ConcreteVisitor1;
//   ClientCode(components, visitor1);
//   std::cout << "\n";
//   std::cout << "It allows the same client code to work with different types of visitors:\n";
//   ConcreteVisitor2 *visitor2 = new ConcreteVisitor2;
//   ClientCode(components, visitor2);

//   for (const Component *comp : components) {
//     delete comp;
//   }
//   delete visitor1;
//   delete visitor2;

//   return 0;
// }
