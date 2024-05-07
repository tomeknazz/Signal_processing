#include "pybind11/pybind11.h"
#include <matplot/matplot.h>
#include <cmath>
#include <pybind11/numpy.h>


float add(float i, float j) {
    return i + j;
}

void sin_plot(double freq, double time)
{
    using namespace matplot;
    fplot([freq](double x) { return sin(2 * pi * freq * x); },
        std::array<double, 2>{0, time});
    grid(on);
    xlabel("t");
    ylabel("y");

    auto ax = gca();
    ax->x_axis().tick_values(iota(0, 0.5, time));


    show();
}

void test_plot()
{
    using namespace matplot;
    auto [X, Y] = meshgrid(iota(-3, .125, 3));
    auto Z = peaks(X, Y);
    auto C = transform(X, Y, [](double x, double y) { return x * y; });
    surfc(X, Y, Z, C);
    colorbar();

    show();
}


namespace py = pybind11;

PYBIND11_MODULE(Signal, m)
{
    m.doc() = "Signal processing";
    m.def("add", &add, "A function which adds two numbers");
    m.def("plot", &sin_plot, "A function which plots a graph");
    m.def("test_plot", &test_plot, "test");
}