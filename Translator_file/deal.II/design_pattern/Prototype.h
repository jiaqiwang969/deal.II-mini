#include <iostream>
#include <string>
#include <unordered_map>

using std::string;

/**
 * Prototype Design Pattern
 *  意图: 让你复制现有的对象而不使你的代码依赖于它们的类。
 */
namespace DesignPattern
{

    namespace Prototype
    {

        enum Type
        {
            PROTOTYPE_1 = 0,
            PROTOTYPE_2
        };

        /**
 * 具有克隆能力的实例类。我们将看到不同类型的字段的值是如何被克隆的。
 */

        class Prototype
        {
        protected:
            string prototype_name_;
            float prototype_field_;

        public:
            Prototype() {}
            Prototype(string prototype_name)
                : prototype_name_(prototype_name)
            {
            }
            virtual ~Prototype() {}
            virtual Prototype *Clone() const = 0;
            virtual void Method(float prototype_field)
            {
                this->prototype_field_ = prototype_field;
                std::cout << "Call Method from " << prototype_name_ << " with field : " << prototype_field << std::endl;
            }
        };

        /**
 * ConcretePrototype1是Prototype的一个子类，实现了Clone方法 在这个例子中，Prototype类的所有数据成员都在堆栈中。如果你的属性中有指针，例如：String* name_ ，你将需要实现Copy-Constructor以确保你有一个来自clone方法的深度拷贝。
 */

        class ConcretePrototype1 : public Prototype
        {
        private:
            float concrete_prototype_field1_;

        public:
            ConcretePrototype1(string prototype_name, float concrete_prototype_field)
                : Prototype(prototype_name), concrete_prototype_field1_(concrete_prototype_field)
            {
            }

            /**
   * 注意，Clone方法返回一个指向新的ConcretePrototype1副本的指针。因此，客户端（调用Clone方法的人）有责任释放该内存。如果你有智能指针的知识，你可能更喜欢在这里使用unique_pointer。
   */
            Prototype *Clone() const override
            {
                return new ConcretePrototype1(*this);
            }
        };

        class ConcretePrototype2 : public Prototype
        {
        private:
            float concrete_prototype_field2_;

        public:
            ConcretePrototype2(string prototype_name, float concrete_prototype_field)
                : Prototype(prototype_name), concrete_prototype_field2_(concrete_prototype_field)
            {
            }
            Prototype *Clone() const override
            {
                return new ConcretePrototype2(*this);
            }
        };

        /**
 * 在PrototypeFactory中，你有两个具体的原型，每个具体的原型类都有一个，所以每次你想创建一个子弹，你可以使用现有的原型并克隆这些原型。
 */

        class PrototypeFactory
        {
        private:
            std::unordered_map<Type, Prototype *, std::hash<int>> prototypes_;

        public:
            PrototypeFactory()
            {
                prototypes_[Type::PROTOTYPE_1] = new ConcretePrototype1("PROTOTYPE_1 ", 50.f);
                prototypes_[Type::PROTOTYPE_2] = new ConcretePrototype2("PROTOTYPE_2 ", 60.f);
            }

            /**
   * 要注意释放所有分配的内存。同样，如果你有智能指针，那么在这里使用它将会更好。
   */

            ~PrototypeFactory()
            {
                delete prototypes_[Type::PROTOTYPE_1];
                delete prototypes_[Type::PROTOTYPE_2];
            }

            /**
   * Notice here that you just need to specify the type of the prototype you
   * want and the method will create from the object with this type.
   */
            Prototype *CreatePrototype(Type type)
            {
                return prototypes_[type]->Clone();
            }
        };

        void Client(PrototypeFactory &prototype_factory)
        {
            std::cout << "Let's create a Prototype 1\n";

            Prototype *prototype = prototype_factory.CreatePrototype(Type::PROTOTYPE_1);
            prototype->Method(90);
            delete prototype;

            std::cout << "\n";

            std::cout << "Let's create a Prototype 2 \n";

            prototype = prototype_factory.CreatePrototype(Type::PROTOTYPE_2);
            prototype->Method(10);

            delete prototype;
        }

    }
}