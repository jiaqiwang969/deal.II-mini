template <class Iter, class Pred>
int count(Iter begin, Iter end, Pred p)
{
    int cnt = 0;
    while (begin != end)
    {
        if (p(*begin))
            ++cnt;
        ++begin;
    }
    return cnt;
}