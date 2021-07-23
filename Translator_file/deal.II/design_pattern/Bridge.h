#include <iostream>
#include <string>

/**
 * Bridge Design Pattern
 *
 * 意图: 让你把一个大类或一组密切相关的类分成两个独立的层次--抽象和实现--它们可以互相独立开发。
 *
 *               A
 *            /     \                        A         N
 *          Aa      Ab        ===>        /     \     / \
 *         / \     /  \                 Aa(N) Ab(N)  1   2
 *       Aa1 Aa2  Ab1 Ab2
 */

/**
 * 实现定义了所有实现类的接口。它不一定要与抽象的接口相匹配。事实上，这两个接口可以完全不同。通常情况下，实现接口只提供原始操作，而抽象接口在这些原始操作的基础上定义了更高层次的操作。
 */

namespace Bridge
{
    class Implementation
    {
    public:
        virtual ~Implementation() {}
        virtual std::string OperationImplementation() const = 0;
    };

    /**
 * 每个具体的实现对应于一个特定的平台，并使用该平台的API实现实现接口。
 */
    class ConcreteImplementationA : public Implementation
    {
    public:
        std::string OperationImplementation() const override
        {
            return "ConcreteImplementationA: Here's the result on the platform A.\n";
        }
    };
    class ConcreteImplementationB : public Implementation
    {
    public:
        std::string OperationImplementation() const override
        {
            return "ConcreteImplementationB: Here's the result on the platform B.\n";
        }
    };

    /**
 * 抽象（Abstraction）定义了两个类层次的 "控制 "部分的接口。它维护一个对实现层次结构的对象的引用，并将所有的实际工作委托给这个对象。
 */

    class Abstraction
    {
        /**
   * @var Implementation
   */
    protected:
        Implementation *implementation_;

    public:
        Abstraction(Implementation *implementation) : implementation_(implementation)
        {
        }

        virtual ~Abstraction()
        {
        }

        virtual std::string Operation() const
        {
            return "Abstraction: Base operation with:\n" +
                   this->implementation_->OperationImplementation();
        }
    };
    /**
 * 你可以扩展抽象而不改变实现类。
 */
    class ExtendedAbstraction : public Abstraction
    {
    public:
        ExtendedAbstraction(Implementation *implementation) : Abstraction(implementation)
        {
        }
        std::string Operation() const override
        {
            return "ExtendedAbstraction: Extended operation with:\n" +
                   this->implementation_->OperationImplementation();
        }
    };

    /**
 * 除了在初始化阶段，一个抽象对象与一个特定的实现对象相联系，客户机代码应该只依赖于抽象类。这样一来，客户端代码就可以支持任何抽象-实现的组合。
 */
    void ClientCode(const Abstraction &abstraction)
    {
        // ...
        std::cout << abstraction.Operation();
        // ...
    }
    /**
 * 客户端代码应该能够与任何预先配置的抽象-实现组合一起工作。
 */

}