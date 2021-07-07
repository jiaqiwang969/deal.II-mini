#include <cstdio>
#include <complex>

void mandelbrot(int image[], int xdim, int ydim, int max_iter)
{
    for (int y = 0; y < ydim; ++y)
    {
        for (int x = 0; x < xdim; ++x)
        { // <<<<< Breakpoint here
            std::complex<float> xy(-2.05 + x * 3.0 / xdim, -1.5 + y * 3.0 / ydim);
            std::complex<float> z(0, 0);
            int count = max_iter;
            for (int i = 0; i < max_iter; ++i)
            {
                z = z * z + xy;
                if (std::abs(z) >= 2)
                {
                    count = i;
                    break;
                }
            }
            image[y * xdim + x] = count;
        }
    }
}

int main()
{
    const int xdim = 500;
    const int ydim = 500;
    const int max_iter = 100;
    int image[xdim * ydim] = {0};
    mandelbrot(image, xdim, ydim, max_iter);
    for (int y = 0; y < ydim; y += 10)
    {
        for (int x = 0; x < xdim; x += 5)
        {
            putchar(image[y * xdim + x] < max_iter ? '.' : '#');
        }
        putchar('\n');
    }
    return 0;
}