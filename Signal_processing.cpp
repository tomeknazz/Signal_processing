#include "pybind11/pybind11.h"
#include <matplot/matplot.h>
#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <random>

using namespace std;
using namespace matplot;

constexpr double plot_x_size = 3.0;
constexpr int samples = 5000;

void peak_in_signal()
{
	srand(time(NULL));
	double amplitude = rand() % 1000 + 1;
	double duration = 2 * pi;
	double sampling_rate = samples / duration;

	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> distribution(0, 1); // Szum
	vector<double> time(samples);
	vector<double> signal(samples);

	for (int i = 0; i < samples; ++i)
	{
		time[i] = i / sampling_rate;
		signal[i] = distribution(gen) * amplitude;
	}
	plot(time, signal);
	show();

	double peak_value = signal[0]; 
	size_t peak_index = 0;
	for (size_t i = 1; i < signal.size(); ++i)
	{
		if (signal[i] > peak_value)
		{
			peak_value = signal[i];
			peak_index = i;
		}
	}
	cout << "peak value is :" << peak_value << endl;
}

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

void dft(const vector<double>& signal, const vector<double>& time, double sampling_rate)
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
	double duration = 2 * pi;
	double sampling_rate = samples / duration;

	vector<double> time;
	for (int i = 0; i < samples; ++i)
	{
		time.push_back(i / sampling_rate);
	}

	vector<double> signal;
	for (int i = 0; i < time.size(); ++i)
	{
		double t = time[i];
		signal.push_back(amplitude * cos(2 * pi * frequency * t));
	}

	dft(signal, time, sampling_rate);
}

void square_wave(double amplitude, double frequency) {


	auto X = linspace(0, plot_x_size, samples);
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

	show();
}

void sawtooth_wave(double amplitude, double frequency)
{
	auto X = linspace(0, plot_x_size, samples);
	auto Y = transform(X, [frequency, amplitude](double x) { return amplitude * (2 * (x * frequency - floor(1 / 2 + x * frequency)) - 1); });
	plot(X, Y);

	auto ax = gca();
	double ylim_min = *std::min_element(Y.begin(), Y.end()) - 0.5;
	double ylim_max = *std::max_element(Y.begin(), Y.end()) + 0.5;
	ax->ylim({ ylim_min, ylim_max });

	show();

}

void sin_wave(double amplitude, double frequency)
{
	auto X = linspace(0, plot_x_size, samples);
	auto Y = transform(X, [frequency, amplitude](double x) { return amplitude * sin(2 * pi * frequency * x); });
	plot(X, Y);

	auto ax = gca();
	double ylim_min = *std::min_element(Y.begin(), Y.end()) - 0.5;
	double ylim_max = *std::max_element(Y.begin(), Y.end()) + 0.5;
	ax->ylim({ ylim_min, ylim_max });

	show();
}

void cos_wave(double amplitude, double frequency)
{
	auto X = linspace(0, plot_x_size, samples);
	auto Y = transform(X, [frequency, amplitude](double x) { return amplitude * cos(2 * pi * frequency * x); });
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
	m.def("sin_wave", &sin_wave, "Plots sin wave signal");
	m.def("cos_wave", &cos_wave, "Plots cos wave signal");
	m.def("peak", &peak_in_signal, "Plots cos wave signal");
}
