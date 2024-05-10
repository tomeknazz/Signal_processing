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

using namespace std;
//spectrum - liczba zespolona ofc ;)) jutro do dokoncze raczej


vector<double> dft(const vector<double>& x) 
{
    int N = x.size();
    vector<double> spectrum(N, 0.0);
    for (int k = 0; k < N; ++k) {
        complex<double> sum{ 0.0, 0.0 };
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * M_PI * k * n / N;
            std::complex<double> exp_term(cos(angle),sin(angle));
            sum += x[n] * exp_term;
        }
        spectrum[k] = abs(sum); 
    }

    return spectrum;
}

void plot_spectrum(const std::vector<double>& spectrum) {
    using namespace matplot;
    auto N = spectrum.size();
    auto f = linspace(0, 0.5, N); // Wektor częstotliwości
    plot(f, spectrum);
    xlabel("Frequency");
    ylabel("Amplitude");
    title("Amplitude Spectrum");
    show();
}

void sin_transform_plot(double amplitude, double freq, double time)
{
    using namespace matplot;
    auto t = linspace(0, time, 1000); 
    vector<double> y_values; 
    for (auto x : t) {
        y_values.push_back(amplitude * sin(2 * M_PI * freq * x));
    }
    auto spectrum = dft(y_values);
    plot_spectrum(spectrum);
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
    m.def("sin_transform_plot", &sin_transform_plot, "zmiana sygnalu dft");
}
