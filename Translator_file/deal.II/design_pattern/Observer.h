

#include <iostream>
#include <list>
#include <string>

namespace DesignPattern
{
    namespace BehavioralDesign
    {

        namespace Observer
        {
            /**
             * Observer
             *
             * - 模式动机:
             *
             * 让你定义一个订阅机制来通知多个对象关于他们正在观察的对象发生的任何事件。
             *
             * 建立一种对象与对象之间的依赖关系，一个对象发生改变时将自动通知其他对象，其他对象将相应做出反应。在此，发生改变的对象称为观察目标，而被通知的对象称为观察者，一个观察目标可以对应多个观察者，而且这些观察者之间没有相互联系，可以根据需要增加和删除观察者，使得系统更易于扩展，这就是观察者模式的模式动机。
             *
             * 请注意，有很多不同的术语与这种模式有类似的含义。只要记住，主体也被称为发布者，而观察者通常被称为订阅者，反之亦然。另外，动词
             * "观察"、"倾听 "或 "跟踪 "通常也是同一个意思。
             */

            /**
             * - 模式定义:
             *
             * 观察者模式(Observer
             * Pattern)：定义对象间的一种一对多依赖关系，使得每当一个对象状态发生改变时，其相关依赖对象皆得到通知并被自动更新。观察者模式又叫做发布-订阅（Publish/Subscribe）模式、模型-视图（Model/View）模式、源-监听器（Source/Listener）模式或从属者（Dependents）模式。
             *
             * 观察者模式是一种对象行为型模式。
             */

            class IObserver
            {
            public:
                virtual ~IObserver(){};
                virtual void Update(const std::string &message_from_subject) = 0;
            };

            class ISubject
            {
            public:
                virtual ~ISubject(){};
                virtual void Attach(IObserver *observer) = 0;
                virtual void Detach(IObserver *observer) = 0;
                virtual void Notify() = 0;
            };

            /**
             * 主体拥有一些重要的状态，并在状态改变时通知观察者。
             */

            class Subject : public ISubject
            {
            public:
                virtual ~Subject()
                {
                    std::cout << "Goodbye, I was the Subject.\n";
                }

                /**
             * 订阅管理方法。
             */
                void Attach(IObserver *observer) override
                {
                    list_observer_.push_back(observer);
                }
                void Detach(IObserver *observer) override
                {
                    list_observer_.remove(observer);
                }
                void Notify() override
                {
                    std::list<IObserver *>::iterator iterator = list_observer_.begin();
                    HowManyObserver();
                    while (iterator != list_observer_.end())
                    {
                        (*iterator)->Update(message_);
                        ++iterator;
                    }
                }

                void CreateMessage(std::string message = "Empty")
                {
                    this->message_ = message;
                    Notify();
                }
                void HowManyObserver()
                {
                    std::cout << "There are " << list_observer_.size() << " observers in the list.\n";
                }

                /**
             * 通常情况下，订阅逻辑只是Subject真正能做的一小部分。主题通常持有一些重要的业务逻辑，每当有重要的事情即将发生（或发生后），就会触发一个通知方法。
             */
                void SomeBusinessLogic()
                {
                    this->message_ = "change message message";
                    Notify();
                    std::cout << "I'm about to do some thing important\n";
                }

            private:
                std::list<IObserver *> list_observer_;
                std::string message_;
            };

            class Observer : public IObserver
            {
            public:
                Observer(Subject &subject) : subject_(subject)
                {
                    this->subject_.Attach(this);
                    std::cout << "Hi, I'm the Observer \"" << ++Observer::static_number_ << "\".\n";
                    this->number_ = Observer::static_number_;
                }
                virtual ~Observer()
                {
                    std::cout << "Goodbye, I was the Observer \"" << this->number_ << "\".\n";
                }

                void Update(const std::string &message_from_subject) override
                {
                    message_from_subject_ = message_from_subject;
                    PrintInfo();
                }
                void RemoveMeFromTheList()
                {
                    subject_.Detach(this);
                    std::cout << "Observer \"" << number_ << "\" removed from the list.\n";
                }
                void PrintInfo()
                {
                    std::cout << "Observer \"" << this->number_ << "\": a new message is available --> " << this->message_from_subject_ << "\n";
                }

            private:
                std::string message_from_subject_;
                Subject &subject_;
                static int static_number_;
                int number_;
            };

            int Observer::static_number_ = 0;

            void ClientCode()
            {
                Subject *subject = new Subject;
                Observer *observer1 = new Observer(*subject);
                Observer *observer2 = new Observer(*subject);
                Observer *observer3 = new Observer(*subject);
                Observer *observer4;
                Observer *observer5;

                subject->CreateMessage("Hello World! :D");
                observer3->RemoveMeFromTheList();

                subject->CreateMessage("The weather is hot today! :p");
                observer4 = new Observer(*subject);

                observer2->RemoveMeFromTheList();
                observer5 = new Observer(*subject);

                subject->CreateMessage("My new car is great! ;)");
                observer5->RemoveMeFromTheList();

                observer4->RemoveMeFromTheList();
                observer1->RemoveMeFromTheList();

                delete observer5;
                delete observer4;
                delete observer3;
                delete observer2;
                delete observer1;
                delete subject;
            }

        }
    }
}

// int main() {
//   ClientCode();
//   return 0;
// }
