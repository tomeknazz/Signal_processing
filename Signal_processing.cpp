#include "pybind11/pybind11.h"
#include <matplot/matplot.h>
#include <cmath>
#include <pybind11/numpy.h>
#include <complex>
#include <vector>

#define M_PI 3.14159265358979323846

float add(float i, float j) {
    return i + j;
}

std::vector<std::complex<double>> dft(const std::vector<double>& x) {
    int N = x.size();
    std::vector<std::complex<double>> X(N, { 0.0, 0.0 }); // Wynikowa transformata

    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * M_PI * k * n / N;
            std::complex<double> exp_term(std::cos(angle), std::sin(angle));
            X[k] += x[n] * exp_term;
        }
    }

    return X;
}

void plot_complex(const std::vector<std::complex<double>>& data) {
    using namespace matplot;

    std::vector<double> real_part, imag_part;
    for (const auto& val : data) {
        real_part.push_back(val.real());
        imag_part.push_back(val.imag());
    }

    auto t = linspace(0, real_part.size() - 1, real_part.size());

    subplot(2, 1, 1);
    plot(t, real_part);
    title("Real Part");

    subplot(2, 1, 2);
    plot(t, imag_part);
    title("Imaginary Part");

    show();
}

void sin_transform_plot(double freq, double time)
{
    using namespace matplot;
    auto t = linspace(0, time, 1000); // Generowanie wektora czasu
    std::vector<double> y_values; // Wektor przechowujący wartości sygnału
    for (auto x : t) {
        y_values.push_back(sin(2 * M_PI * freq * x));
    }

    // Obliczanie wyniku DFT
    auto dft_result = dft(y_values);

    // Rysowanie wykresu po przekształceniu DFT
    plot_complex(dft_result);
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
    m.def("test_plot", &test_plot, "test");
    m.def("sin_transform_plot", &sin_transform_plot, "zmiana sygnalu sinusxd");
}
