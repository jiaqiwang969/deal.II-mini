#include <iostream>
#include <string>
#include <vector>

/**
 * Builder Design Pattern
 *
 * 意图: 让你一步一步地构建复杂的对象。该模式允许你使用相同的构造代码产生不同类型和表现形式的对象。
 */

/**
 * 只有当你的产品相当复杂并需要大量配置时，使用构建器模式才有意义。
 * 与其他创建模式不同，不同的具体构建器可以产生不相关的产品。换句话说，各种构建器的结果可能并不总是遵循相同的接口。
 */

class Product1
{
public:
    std::vector<std::string> parts_;
    void ListParts() const
    {
        std::cout << "Product parts: ";
        for (size_t i = 0; i < parts_.size(); i++)
        {
            if (parts_[i] == parts_.back())
            {
                std::cout << parts_[i];
            }
            else
            {
                std::cout << parts_[i] << ", ";
            }
        }
        std::cout << "\n\n";
    }
};

/**
 * 生成器接口指定了用于创建产品对象不同部分的方法。
 */
class Builder
{
public:
    virtual ~Builder() {}
    virtual void ProducePartA() const = 0;
    virtual void ProducePartB() const = 0;
    virtual void ProducePartC() const = 0;
};
/**
 * 具体的Builder类遵循Builder接口并提供具体的构建步骤的实现。你的程序可能有几个不同的Builder，实现方式也不同。
 */
class ConcreteBuilder1 : public Builder
{
private:
    Product1 *product;

    /**
     * 一个新的构建器实例应该包含一个空白的产品对象，它被用于进一步的装配。
     */
public:
    ConcreteBuilder1()
    {
        this->Reset();
    }

    ~ConcreteBuilder1()
    {
        delete product;
    }

    void Reset()
    {
        this->product = new Product1();
    }
    /**
     * All production steps work with the same product instance.
     */

    void ProducePartA() const override
    {
        this->product->parts_.push_back("PartA1");
    }

    void ProducePartB() const override
    {
        this->product->parts_.push_back("PartB1");
    }

    void ProducePartC() const override
    {
        this->product->parts_.push_back("PartC1");
    }

    /**
     * Concrete建造商应该提供他们自己的方法来检索结果。这是因为各种类型的构建器可能会创建完全不同的产品，它们不遵循相同的接口。因此，这种方法不能在基本的Builder接口中声明（至少在静态类型的编程语言中）。请注意，PHP是一种动态类型的语言，这种方法可以在基础接口中出现。然而，为了清楚起见，我们不会在那里声明它。
     * 通常情况下，在将最终结果返回给客户之后，一个构建器实例被期望准备好开始生产另一个产品。这就是为什么通常的做法是在`getProduct`方法体的末尾调用重置方法。然而，这种行为并不是强制性的，你可以让你的构建器在处理之前的结果之前等待来自客户端代码的明确的重置调用。
     */

    /**
     * 请注意这里的内存所有权。一旦你调用GetProduct，这个函数的用户就有责任释放这段内存。这里可能是一个更好的选择，即使用智能指针来避免内存泄漏。
     */

    Product1 *GetProduct()
    {
        Product1 *result = this->product;
        this->Reset();
        return result;
    }
};

/**
 * 主任只负责按特定顺序执行建造步骤。当按照特定的顺序或配置生产产品时，它很有帮助。严格来说，Director类是可选的，因为客户可以直接控制建造者。
 */
class Director
{
    /**
     * @var Builder
     */
private:
    Builder *builder;
    /**
     * Director与客户代码传递给它的任何构建器实例一起工作。这样，客户代码可以改变新组装产品的最终类型。
     */

public:
    void set_builder(Builder *builder)
    {
        this->builder = builder;
    }

    /**
     * 主任可以使用相同的建造步骤建造多个产品变体。
     */

    void BuildMinimalViableProduct()
    {
        this->builder->ProducePartA();
    }

    void BuildFullFeaturedProduct()
    {
        this->builder->ProducePartA();
        this->builder->ProducePartB();
        this->builder->ProducePartC();
    }
};
/**
 * 客户端代码创建一个构建器对象，将其传递给主任，然后启动构建过程。最终结果从建造者对象中获取。
 */
/**
 * 为了简单起见，我使用了原始指针，但是你可能更喜欢在这里使用智能指针。
 */
void ClientCode(Director &director)
{
    ConcreteBuilder1 *builder = new ConcreteBuilder1();
    director.set_builder(builder);
    std::cout << "Standard basic product:\n";
    director.BuildMinimalViableProduct();

    Product1 *p = builder->GetProduct();
    p->ListParts();
    delete p;

    std::cout << "Standard full featured product:\n";
    director.BuildFullFeaturedProduct();

    p = builder->GetProduct();
    p->ListParts();
    delete p;

    // Remember, the Builder pattern can be used without a Director class.
    std::cout << "Custom product:\n";
    builder->ProducePartA();
    builder->ProducePartC();
    p = builder->GetProduct();
    p->ListParts();
    delete p;

    delete builder;
}