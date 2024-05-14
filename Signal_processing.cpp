#include "pybind11/pybind11.h"
#include <matplot/matplot.h>
#include <cmath>
#include <pybind11/numpy.h>
#include <complex>
#include <vector>
#include <string>
#include <iostream>
using namespace std;
using namespace matplot;


void inverse_dft(const vector<complex<double>>& spectrum, vector<double> time)
{
	int N = spectrum.size();
	vector<double> signal;
	signal.resize(N);

	for (int n = 0; n < N; ++n)
	{
		complex<double> sum(0.0, 0.0);
		for (int k = 0; k < N; ++k)
		{
			sum += spectrum[k] * exp(2.0 * pi * complex<double>(0, 1) * static_cast<double>(n * k) / static_cast<double>(N));
		}
		signal[n] = sum.real();
	}
	figure();
	plot(time, signal);
	title("Inverse DFT Signal");
	xlabel("Time");
	ylabel("Amplitude");
	show();
}

void dft(double amplitude, double frequency, const vector<double>& signal, const vector<double>& time,double sampling_rate)
{
	using namespace matplot;
	using namespace std;
	int N = signal.size();

	vector<complex<double>> dft_result;
	for (int k = 0; k < N; ++k)
	{
		complex<double> sum(0.0, 0.0);
		for (int n = 0; n < N; ++n)
		{
			sum += signal[n] * exp(-2.0 * pi * complex<double>(0, 1) * static_cast<double>(n * k) / static_cast<double>(N));
		}
		dft_result.push_back(sum / static_cast<double>(N));
	}

	vector<double> frequencies;
	for (int k = 0; k < N; ++k)
	{
		double freq = static_cast<double>(k) * sampling_rate / N;
		frequencies.push_back(freq);
	}

	vector<double> magnitude;
	for (int k = 0; k < N; ++k)
	{
		magnitude.push_back(abs(dft_result[k]) * 2.02);
	}

	plot(frequencies, magnitude);
	title("DFT Results");
	xlabel("Frequency (Hz)");
	ylabel("Magnitude");
	show();
	inverse_dft(dft_result, time);
}

void create_signal_for_dft(double amplitude, double frequency)
{
	int N = 1000;
	double duration = 2 * pi;
	double sampling_rate = N / duration;

	vector<double> time;
	for (int i = 0; i < N; ++i)
	{
		time.push_back(i / sampling_rate);
	}

	vector<double> signal;
	for (int i = 0; i < time.size(); ++i)
	{
		double t = time[i];
		signal.push_back(amplitude * cos(2 * pi * frequency * t));
	}

	dft(amplitude, frequency, signal, time, sampling_rate);
}

void square_wave(double amplitude, double frequency) {
	using namespace matplot;

	auto X = linspace(0, 3, 5000);
	std::vector<std::vector<double>> Y(2);
	Y[0] = transform(X, [frequency](double x) { return sin(2 * pi * frequency * x); });


	// Modify Y to create a square wave
	for (size_t i = 0; i < Y[0].size(); ++i) {
		Y[0][i] = (Y[0][i] > 0) ? amplitude : -1.0 * amplitude;  // Assuming amplitude of 1
	}

	figure();
	stairs(X, Y);  // Plot the square wave

	// Adjust y-axis limits to provide space between the wave and the edges of the plot
	auto ax = gca();
	double ylim_min = *std::min_element(Y[0].begin(), Y[0].end()) - 0.5;
	double ylim_max = *std::max_element(Y[0].begin(), Y[0].end()) + 0.5;
	ax->ylim({ ylim_min, ylim_max });

	ax->x_axis().ticklabels({ "0", "pi/4", "pi/2", "3pi/4", "pi" });


	show();
}

void sawtooth_wave(double amplitude, double frequency)
{
	using namespace matplot;

	auto X = linspace(0, 3, 5000);
	auto Y = transform(X, [frequency, amplitude](double x) { return amplitude * (2 * (x * frequency - floor(1 / 2 + x * frequency)) - 1); });
	plot(X, Y);

	auto ax = gca();
	double ylim_min = *std::min_element(Y.begin(), Y.end()) - 0.5;
	double ylim_max = *std::max_element(Y.begin(), Y.end()) + 0.5;
	ax->ylim({ ylim_min, ylim_max });

	show();

}

PYBIND11_MODULE(Signal, m)
{
	m.doc() = "Signal processing";
	m.def("DFT", &create_signal_for_dft, "Transforms and plots sinusoidal signal using DFT");
	m.def("square_wave", &square_wave, "Plots square wave signal");
	m.def("sawtooth_wave", &sawtooth_wave, "Plots sawtooth wave signal");
}
