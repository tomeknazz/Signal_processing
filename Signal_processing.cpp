#include "pybind11/pybind11.h"
#include <matplot/matplot.h>
#include <cmath>
#include <pybind11/numpy.h>
#include <complex>
#include <vector>

#define M_PI 3.14159265358979323846

void sin_transform_plot(double amplitude, double frequency)
{
    using namespace matplot;
    using namespace std;
    int N = 1000; 
    double duration = 2 * M_PI;
    double sampling_rate = N / duration;

    vector<double> time;
    for (int i = 0; i < N; ++i) {
        time.push_back(i / sampling_rate);
    }

    vector<double> signal; // wartosci sinusa
    for (int i = 0; i < time.size(); ++i) {
        double t = time[i];
        signal.push_back(amplitude * sin(2 * M_PI * frequency * t));
    }

    
    vector<complex<double>> dft_result; // tu oblczam dft
    for (int k = 0; k < N; ++k) 
    {
        complex<double> sum(0.0, 0.0);

        for (int n = 0; n < N; ++n) {
            sum += signal[n] * exp(-2.0 * M_PI * complex<double>(0, 1) * static_cast<double>(n * k) / static_cast<double>(N)); // rownanie wrzucone na wiki 
        }
        dft_result.push_back(sum / static_cast<double>(N));
    }

    vector<double> frequencies;
    for (int k = 0; k < N; ++k) 
    {
        double freq;
        if (k <= N / 2) {
            freq = static_cast<double>(k) * sampling_rate / N;
        }
        else {
            freq = -(sampling_rate - static_cast<double>(k)) * sampling_rate / N;
        }
        frequencies.push_back(freq);
    }

    vector<double> real_parts;
    vector<double> imaginary_parts;
    for (int k = 0; k < N; ++k) {
        real_parts.push_back(abs(dft_result[k].real()) * 2);
        imaginary_parts.push_back(abs(dft_result[k].imag()) * 2);
    }

    plot(frequencies, real_parts);
    plot(frequencies, imaginary_parts);
    title("DFT Results");
    xlabel("Frequency (Hz)");
    ylabel("Amplitude");

    show();
}

namespace py = pybind11;

PYBIND11_MODULE(Signal, m)
{
    m.doc() = "Signal processing";
    m.def("sin_transform_plot", &sin_transform_plot, "Transforms and plots sinusoidal signal using DFT");
}
