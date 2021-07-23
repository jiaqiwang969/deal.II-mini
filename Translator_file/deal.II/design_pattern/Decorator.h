#include <iostream>
#include <string>

/**
 * Decorator
 *
 * - 模式动机:
 * 
 *   让你把新的行为附加到对象上，把这些对象放在包含行为的特殊包装对象中。
 *
 *   一般有两种方式可以实现给一个类或对象增加行为： 
 * 
 *   -- 继承机制，使用继承机制是给现有类添加功能的一种有效途径，通过继承一个现有类可以使得子类在拥有自身方法的同时还拥有父类的方法。但是这种方法是静态的，用户不能控制增加行为的方式和时机。
 *
 *   -- 关联机制，即将一个类的对象嵌入另一个对象中，由另一个对象来决定是否调用嵌入对象的行为以便扩展自己的行为，我们称这个嵌入的对象为装饰器(Decorator).
 * 
 * 装饰模式以对客户透明的方式动态地给一个对象附加上更多的责任，换言之，客户端并不会觉得对象在装饰前和装饰后有什么不同。装饰模式可以在不需要创造更多子类的情况下，将对象的功能加以扩展。这就是装饰模式的模式动机。
 */
/**
 * - 模式定义: 
 * 
 *   装饰模式(Decorator Pattern)：动态地给一个对象增加一些额外的职责(Responsibility)，就增加对象功能来说，装饰模式比生成子类实现更为灵活。其别名也可以称为包装器(Wrapper)，与适配器模式的别名相同，但它们适用于不同的场合。根据翻译的不同，装饰模式也有人称之为“油漆工模式”，它是一种对象结构型模式。
 * 
 *   基础组件接口定义了可由装饰器改变的操作。
 */
namespace DesignPattern
{

    namespace Decorator
    {
        class Component
        {
        public:
            virtual ~Component() {}
            virtual std::string Operation() const = 0;
        };
        /**
         * 具体组件提供操作的默认实现。这些类可能有几种变化。
         */
        class ConcreteComponent : public Component
        {
        public:
            std::string Operation() const override
            {
                return "ConcreteComponent";
            }
        };
        /**
         * 基础装饰器类遵循与其他组件相同的接口。这个类的主要目的是为所有具体的装饰器定义包装接口。包裹代码的默认实现可能包括一个用于存储被包裹组件的字段和初始化它的方法。
         */
        class Decorator : public Component
        {
            /**
             * @var Component
             */
        protected:
            Component *component_;

        public:
            Decorator(Component *component) : component_(component)
            {
            }
            /**
             * The Decorator delegates all work to the wrapped component.
             */
            std::string Operation() const override
            {
                return this->component_->Operation();
            }
        };
        /**
         * Concrete Decorators call the wrapped object and alter its result in
         * some way.
         */
        class ConcreteDecoratorA : public Decorator
        {
            /**
             * 装饰器可以调用操作的父级实现，而不是直接调用被包装的对象。这种方法简化了装饰器类的扩展。
             */
        public:
            ConcreteDecoratorA(Component *component) : Decorator(component)
            {
            }
            std::string Operation() const override
            {
                return "ConcreteDecoratorA(" + Decorator::Operation() + ")";
            }
        };
        /**
         * 装饰器可以在调用包装对象之前或之后执行其行为。
         */
        class ConcreteDecoratorB : public Decorator
        {
        public:
            ConcreteDecoratorB(Component *component) : Decorator(component)
            {
            }

            std::string Operation() const override
            {
                return "ConcreteDecoratorB(" + Decorator::Operation() + ")";
            }
        };
        /**
         * 客户端代码与所有使用组件接口的对象一起工作。这样，它就可以保持独立于它所使用的组件的具体类。
         */
        void ClientCode(Component *component)
        {
            // ...
            std::cout << "RESULT: " << component->Operation();
            // ...
        }

    }
}
