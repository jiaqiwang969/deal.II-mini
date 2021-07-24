#include <algorithm>
#include <iostream>
#include <list>
#include <string>

namespace Pattern
{

    /**
     * Composite
     *
     * - 模式动机:
     *
     * 让你把对象组成树状结构，然后像处理单个对象一样处理这些结构。
     */
    /**
     * - 模式定义：
     *
     * 基层的Component类声明了一个组合的简单和复杂对象的共同操作。
     */

    namespace Composite
    {
        class Component
        {
            /**
             * @var Component
             */
        protected:
            Component *parent_;
            /**
             * 作为选择，基础组件可以声明一个接口，用于设置和访问树状结构中组件的父级。它也可以为这些方法提供一些默认实现。
             */
        public:
            virtual ~Component() {}
            void SetParent(Component *parent)
            {
                this->parent_ = parent;
            }
            Component *GetParent() const
            {
                return this->parent_;
            }
            /**
             * 在某些情况下，在基础组件类中直接定义儿童管理操作将是有益的。这样，你就不需要向客户端代码暴露任何具体的组件类，即使是在对象树装配期间。缺点是，这些方法对于叶级组件来说将是空的。
             */
            virtual void Add(Component *component) {}
            virtual void Remove(Component *component) {}
            /**
             * 你可以提供一个方法，让客户端代码弄清一个组件是否可以产生子类。
             */
            virtual bool IsComposite() const
            {
                return false;
            }
            /**
             * 基组件可以实现一些默认的行为，或者将其留给具体的类（通过将包含该行为的方法声明为
             * "抽象的"）。
             */
            virtual std::string Operation() const = 0;
        };
        /**
         * 叶子类代表一个组合的终端对象。一个叶子不能有任何孩子。
         * 通常情况下，是叶子对象做实际工作，而复合对象只委托给它们的子构件。
         */
        class Leaf : public Component
        {
        public:
            std::string Operation() const override
            {
                return "Leaf";
            }
        };
        /**
         * 复合类代表了可能有子类的复杂组件。通常情况下，复合对象将实际工作委托给它们的子代，然后将结果
         * "汇总"。
         */
        class Composite : public Component
        {
            /**
             * @var \SplObjectStorage
             */
        protected:
            std::list<Component *> children_;

        public:
            /**
             * 一个复合对象可以将其他组件（包括简单的或复杂的）加入或移出其子列表。
             */
            void Add(Component *component) override
            {
                this->children_.push_back(component);
                component->SetParent(this);
            }
            /**
             * 请记住，这个方法删除了列表的指针，但没有释放内存，你应该手动操作，或者最好使用智能指针。
             */
            void Remove(Component *component) override
            {
                children_.remove(component);
                component->SetParent(nullptr);
            }
            bool IsComposite() const override
            {
                return true;
            }
            /**
             * 复合体以一种特殊的方式执行其主要逻辑。它递归地遍历其所有的子代，收集和计算其结果。由于复合体的子代将这些调用传递给他们的子代，以此类推，整个对象树被遍历的结果。
             */
            std::string Operation() const override
            {
                std::string result;
                for (const Component *c : children_)
                {
                    if (c == children_.back())
                    {
                        result += c->Operation();
                    }
                    else
                    {
                        result += c->Operation() + "+";
                    }
                }
                return "Branch(" + result + ")";
            }
        };
        /**
         * The client code works with all of the components via the base
         * interface.
         */
        void ClientCode(Component *component)
        {
            // ...
            std::cout << "RESULT: " << component->Operation();
            // ...
        }

        /**
         * 由于儿童管理操作是在基础组件类中声明的，客户代码可以与任何组件一起工作，不管是简单的还是复杂的，都不需要依赖它们的具体类。
         */
        void ClientCode2(Component *component1, Component *component2)
        {
            // ...
            if (component1->IsComposite())
            {
                component1->Add(component2);
            }
            std::cout << "RESULT: " << component1->Operation();
            // ...
        }

    }
}


/**
 * 这样，客户端代码就可以支持简单的叶子组件...
 */
