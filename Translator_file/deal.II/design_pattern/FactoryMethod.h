#include <iostream>
#include <string>

/**
 * 工厂设计方法
 *
 * 意图：为在父类中创建对象提供一个接口，但允许子类改变将被创建的对象的类型。
 */

/**
 * 产品接口声明了所有具体产品必须 实现的操作。
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
     * 还要注意的是，尽管它的名字叫 "创造者"，但它的主要责任不是创造产品。通常，它包含一些依赖于产品对象的核心业务逻辑，这些对象由工厂方法返回。子类可以通过重写工厂方法并从中返回不同类型的产品来间接地改变该业务逻辑。
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
   *  请注意，方法的签名仍然使用抽象产品类型，尽管具体产品实际上是从该方法中返回的。这样，创造者就可以保持独立于具体的产品类型。
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