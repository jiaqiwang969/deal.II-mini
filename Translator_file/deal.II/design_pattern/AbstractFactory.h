#include <iostream>
#include <string>

/**
 * Abstract Factory Design Pattern
 *
 * 意图: 让你产生相关对象的系列，而不需要指定它们的具体类别。
 */

/**
 * 一个产品系列的每个不同产品都应该有一个基本接口。该产品的所有变体都必须实现这个接口。
 */
class AbstractProductA
{
public:
    virtual ~AbstractProductA(){};
    virtual std::string UsefulFunctionA() const = 0;
};

/**
 * Concrete产品是由相应的Concrete工厂创造的。
 */
class ConcreteProductA1 : public AbstractProductA
{
public:
    std::string UsefulFunctionA() const override
    {
        return "The result of the product A1.";
    }
};

class ConcreteProductA2 : public AbstractProductA
{
    std::string UsefulFunctionA() const override
    {
        return "The result of the product A2.";
    }
};

/**
 * 这里是另一个产品的基本接口。所有的产品都可以相互作用，但只有在相同的具体变体的产品之间才能进行适当的互动。
 */
class AbstractProductB
{
    /**
   * 产品B能够做它自己的事情...
   */
public:
    virtual ~AbstractProductB(){};
    virtual std::string UsefulFunctionB() const = 0;
    /**
   *   ...但它也可以与ProductA协作。
   * 抽象工厂确保它所创建的所有产品都具有相同的变体，因此是兼容的。
   */
    virtual std::string AnotherUsefulFunctionB(const AbstractProductA &collaborator) const = 0;
};

/**
 * Concrete产品是由相应的Concrete工厂创造的。
 */
class ConcreteProductB1 : public AbstractProductB
{
public:
    std::string UsefulFunctionB() const override
    {
        return "The result of the product B1.";
    }
    /**
   * 变体 "产品B1 "只能够与变体 "产品A1 "正常工作。然而，它接受AbstractProductA的任何实例作为参数。
   */
    std::string AnotherUsefulFunctionB(const AbstractProductA &collaborator) const override
    {
        const std::string result = collaborator.UsefulFunctionA();
        return "The result of the B1 collaborating with ( " + result + " )";
    }
};

class ConcreteProductB2 : public AbstractProductB
{
public:
    std::string UsefulFunctionB() const override
    {
        return "The result of the product B2.";
    }
    /**
   * 变体 "Product B2 "只能与变体 "Product A2 "正常工作。然而，它接受AbstractProductA的任何实例作为参数。
   */
    std::string AnotherUsefulFunctionB(const AbstractProductA &collaborator) const override
    {
        const std::string result = collaborator.UsefulFunctionA();
        return "The result of the B2 collaborating with ( " + result + " )";
    }
};

/**
 * 抽象工厂接口声明了一组返回不同抽象产品的方法。这些产品被称为一个家族，并通过一个高级的主题或概念进行关联。一个家族的产品通常能够在它们之间进行协作。一个家族的产品可能有几个变体，但一个变体的产品与另一个变体的产品是不兼容的。
 */
class AbstractFactory
{
public:
    virtual AbstractProductA *CreateProductA() const = 0;
    virtual AbstractProductB *CreateProductB() const = 0;
};

/**
 * Concrete工厂生产属于单一变体的产品系列。工厂保证产生的产品是兼容的。请注意，具体工厂的方法的签名会返回一个抽象产品，而在方法内部会实例化一个具体产品。
 */
class ConcreteFactory1 : public AbstractFactory
{
public:
    AbstractProductA *CreateProductA() const override
    {
        return new ConcreteProductA1();
    }
    AbstractProductB *CreateProductB() const override
    {
        return new ConcreteProductB1();
    }
};

/**
 * 每个Concrete工厂都有一个相应的产品变体。
 */
class ConcreteFactory2 : public AbstractFactory
{
public:
    AbstractProductA *CreateProductA() const override
    {
        return new ConcreteProductA2();
    }
    AbstractProductB *CreateProductB() const override
    {
        return new ConcreteProductB2();
    }
};

/**
 * 客户端代码只通过抽象类型与工厂和产品一起工作。AbstractFactory和AbstractProduct。这让你可以将任何工厂或产品的子类传递给客户端代码而不破坏它。
 */

void ClientCode(const AbstractFactory &factory)
{
    const AbstractProductA *product_a = factory.CreateProductA();
    const AbstractProductB *product_b = factory.CreateProductB();
    std::cout << product_b->UsefulFunctionB() << "\n";
    std::cout << product_b->AnotherUsefulFunctionB(*product_a) << "\n";
    delete product_a;
    delete product_b;
}