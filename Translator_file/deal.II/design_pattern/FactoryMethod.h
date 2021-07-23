#include <iostream>
#include <string>

namespace DesignPattern
{
    namespace CreationalDesign
    {
        /**
         * FactoryMethod
         *
         * - 模式动机：
         *
         * 为在父类中创建对象提供一个接口，但允许子类改变将被创建的对象的类型。
         *
         * 现在对该系统进行修改，不再设计一个按钮工厂类来统一负责所有产品的创建，而是将具体按钮的创建过程交给专门的工厂子类去完成，我们先定义一个抽象的按钮工厂类，再定义具体的工厂类来生成圆形按钮、矩形按钮、菱形按钮等，它们实现在抽象按钮工厂类中定义的方法。这种抽象化的结果使这种结构可以在不修改具体工厂类的情况下引进新的产品，如果出现新的按钮类型，只需要为这种新类型的按钮创建一个具体的工厂类就可以获得该新按钮的实例，这一特点无疑使得工厂方法模式具有超越简单工厂模式的优越性，更加符合“开闭原则”。
         *
         * - 模式定义：
         *
         * 工厂方法模式(Factory Method
         * Pattern)又称为工厂模式，也叫虚拟构造器(Virtual
         * Constructor)模式或者多态工厂(Polymorphic
         * Factory)模式，它属于类创建型模式。在工厂方法模式中，工厂父类负责定义创建产品对象的公共接口，而工厂子类则负责生成具体的产品对象，这样做的目的是将产品类的实例化操作延迟到工厂子类中完成，即通过工厂子类来确定究竟应该实例化哪一个具体产品类。
         */

        namespace FactoryMethod
        {

            class Product
            {
            public:
                virtual ~Product() {}
                virtual std::string Operation() const = 0;
            };

            /**
             * 具体的产品提供了产品接口的各种实现。
             */
            class ConcreteProduct1 : public Product
            {
            public:
                std::string Operation() const override
                {
                    return "{Result of the ConcreteProduct1}";
                }
            };
            class ConcreteProduct2 : public Product
            {
            public:
                std::string Operation() const override
                {
                    return "{Result of the ConcreteProduct2}";
                }
            };

            /**
             * Creator类声明了工厂方法，它应该返回一个产品类的对象。创造者的子类通常提供这个方法的实现。
             */

            class Creator
            {
                /**
                 * 请注意，创造者也可能提供一些工厂方法的默认实现。
                 */
            public:
                virtual ~Creator(){};
                virtual Product *FactoryMethod() const = 0;
                /**
                 * 还要注意的是，尽管它的名字叫
                 * "创造者"，但它的主要责任不是创造产品。通常，它包含一些依赖于产品对象的核心业务逻辑，这些对象由工厂方法返回。子类可以通过重写工厂方法并从中返回不同类型的产品来间接地改变该业务逻辑。
                 */

                std::string SomeOperation() const
                {
                    // Call the factory method to create a Product object.
                    Product *product = this->FactoryMethod();
                    // Now, use the product.
                    std::string result = "Creator: The same creator's code has just worked with " + product->Operation();
                    delete product;
                    return result;
                }
            };

            /**
             * 具体的创建者覆盖工厂方法，以改变生成的产品的类型。
             */
            class ConcreteCreator1 : public Creator
            {
                /**
                 * 请注意，方法的签名仍然使用抽象产品类型，尽管具体产品实际上是从该方法中返回的。这样，创造者就可以保持独立于具体的产品类型。
                 */
            public:
                Product *FactoryMethod() const override
                {
                    return new ConcreteProduct1();
                }
            };

            class ConcreteCreator2 : public Creator
            {
            public:
                Product *FactoryMethod() const override
                {
                    return new ConcreteProduct2();
                }
            };

            /**
             * 客户端代码与一个具体的创造者的实例一起工作，尽管是通过它的基础接口。只要客户端继续通过基接口与创造者一起工作，你就可以把任何创造者的子类传递给它。
             */
            void ClientCode(const Creator &creator)
            {
                // ...
                std::cout << "Client: I'm not aware of the creator's class, but it still works.\n"
                          << creator.SomeOperation() << std::endl;
                // ...
            }

        }
    }
}
/**
 * 应用程序根据配置或环境来挑选创建者的类型。
 */

/**
 * 应用程序根据配置或环境来挑选创建者的类型。
 */

// int main() {
//   std::cout << "App: Launched with the ConcreteCreator1.\n";
//   Creator* creator = new ConcreteCreator1();
//   ClientCode(*creator);
//   std::cout << std::endl;
//   std::cout << "App: Launched with the ConcreteCreator2.\n";
//   Creator* creator2 = new ConcreteCreator2();
//   ClientCode(*creator2);

//   delete creator;
//   delete creator2;
//   return 0;
// }
