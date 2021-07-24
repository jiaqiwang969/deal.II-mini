#include <iostream>
#include <string>
namespace Pattern
{

        /**
         * AbstractFactory
         *
         * - 模式动机: 让你产生相关对象的系列，而不需要指定它们的具体类别。
         *
         * 在工厂方法模式中具体工厂负责生产具体的产品，每一个具体工厂对应一种具体产品，工厂方法也具有唯一性，一般情况下，一个具体工厂中只有一个工厂方法或者一组重载的工厂方法。但是有时候我们需要一个工厂可以提供多个产品对象，而不是单一的产品对象。
         *
         * 为了更清晰地理解工厂方法模式，需要先引入两个概念：
         *
         * - 产品等级结构:
         *
         * 产品等级结构即产品的继承结构，如一个抽象类是电视机，其子类有海尔电视机、海信电视机、TCL电视机，则抽象电视机与具体品牌的电视机之间构成了一个产品等级结构，抽象电视机是父类，而具体品牌的电视机是其子类。
         *
         * - 产品族:
         * 在抽象工厂模式中，产品族是指由同一个工厂生产的，位于不同产品等级结构中的一组产品，如海尔电器工厂生产的海尔电视机、海尔电冰箱，海尔电视机位于电视机产品等级结构中，海尔电冰箱位于电冰箱产品等级结构中。
         *
         * 当系统所提供的工厂所需生产的具体产品并不是一个简单的对象，而是多个位于不同产品等级结构中属于不同类型的具体产品时需要使用抽象工厂模式。
         *
         * 抽象工厂模式是所有形式的工厂模式中最为抽象和最具一般性的一种形态。
         *
         * 抽象工厂模式与工厂方法模式最大的区别在于，工厂方法模式针对的是一个产品等级结构，而抽象工厂模式则需要面对多个产品等级结构，一个工厂等级结构可以负责多个不同产品等级结构中的产品对象的创建。
         *
         * 当一个工厂等级结构可以创建出分属于不同产品等级结构的一个产品族中的所有对象时，抽象工厂模式比工厂方法模式更为简单、有效率。
         */

        /**
         * - 模式定义: 抽象工厂模式(Abstract Factory Pattern)：提供一个创建一系列相关或相互依赖对象的接口，而无须指定它们具体的类。抽象工厂模式又称为Kit模式，属于对象创建型模式。
         *
         * 一个产品系列的每个不同产品都应该有一个基本接口。该产品的所有变体都必须实现这个接口。
         */

        namespace AbstractFactory
        {

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
                 * ...但它也可以与ProductA协作。
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
                 * 变体 "产品B1 "只能够与变体 "产品A1
                 * "正常工作。然而，它接受AbstractProductA的任何实例作为参数。
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
                 * 变体 "Product B2 "只能与变体 "Product A2
                 * "正常工作。然而，它接受AbstractProductA的任何实例作为参数。
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

        }
}
