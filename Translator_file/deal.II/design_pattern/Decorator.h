#include <iostream>
#include <string>

/**
 * Decorator Design Pattern
 *
 * 意图: 让你把新的行为附加到对象上，把这些对象放在包含行为的特殊包装对象中。
 */
/**
 * 基础组件接口定义了可由装饰器改变的操作。
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
 * Concrete Decorators call the wrapped object and alter its result in some way.
 */
        class ConcreteDecoratorA : public Decorator
        {
            /**
   *  装饰器可以调用操作的父级实现，而不是直接调用被包装的对象。这种方法简化了装饰器类的扩展。
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