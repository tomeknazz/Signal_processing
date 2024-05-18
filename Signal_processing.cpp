#include "pybind11/pybind11.h"
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

using namespace std;
using namespace matplot;

constexpr double plot_x_size = 3.0;
constexpr int samples = 6000;

void peak_in_signal(const vector<double>& signal)
{
	double peak_value = signal[0];
	double trough_value = signal[0];
	size_t peak_index = 0;
	size_t trough_index = 0;
	for (size_t i = 1; i < signal.size(); ++i)
	{
		if (signal[i] > peak_value)
		{
			peak_value = signal[i];
			peak_index = i;
		}
		else if (signal[i] < trough_value)
		{
			trough_value = signal[i];
			trough_index = i;
		}
	}
	cout << "Peak value is: " << peak_value << " at time: " << peak_index / (samples / (2 * pi)) << endl;
	cout << "Trough value is: " << trough_value << " at time: " << trough_index / (samples / (2 * pi)) << endl;
	cin.get();
}

void random_signal()
{
	srand(time(NULL));
	double amplitude = rand() % 1500 + 1;
	vector<double> time(samples);
	vector<double> signal(samples);

	for (int i = 0; i < samples; ++i)
	{
		time[i] = i * (2 * pi) / samples;
		signal[i] = (static_cast<double>(rand()) / RAND_MAX - 0.5) * amplitude * (rand() % 849 + 2);
	}
	plot(time, signal);
	peak_in_signal(signal);
}

void inverse_dft(const vector<complex<double>>& spectrum, const vector<double>& time)
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

void dft(const vector<double>& signal, const vector<double>& time, const double sampling_rate)
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
	inverse_dft(dft_result, time);
}

void plot_function(const vector<double>& x, const vector<double>& y)
{
	figure();
	plot(x, y);
	auto ax = gca();
	double ylim_min = *std::min_element(y.begin(), y.end()) - 0.5;
	double ylim_max = *std::max_element(y.begin(), y.end()) + 0.5;
	ax->ylim({ ylim_min, ylim_max });
	title("Signal");
	xlabel("Time");
	ylabel("Amplitude");
	show();
}

void display_data_from_python(const vector<double>& v, const int samplerate) {
	auto x = linspace(0, plot_x_size, samplerate);
	plot_function(x, v);
}

void create_signal_for_dft(const double amplitude, const double frequency)
{
	constexpr double duration = 2 * pi;
	constexpr double sampling_rate = samples / duration;

	vector<double> time;
	for (int i = 0; i < samples; ++i)
	{
		time.push_back(i / sampling_rate);
	}

	vector<double> signal;
	for (int i = 0; i < time.size(); ++i)
	{
		const double t = time[i];
		signal.push_back(amplitude * cos(2 * pi * frequency * t)+5*amplitude*sin(2*pi*5*frequency*t));
	}

	dft(signal, time, sampling_rate);
}

void square_wave(const double amplitude, double frequency) {
	const auto x = linspace(0, plot_x_size, samples);
	auto y = transform(x, [frequency](const double x) { return sin(2 * pi * frequency * x); });

	for (size_t i = 0; i < y.size(); ++i) {
		y[i] = (y[i] > 0) ? amplitude : -1.0 * amplitude;  // Assuming amplitude of 1
	}

	plot_function(x, y);
}

void sawtooth_wave(double amplitude, double frequency)
{
	const auto x = linspace(0, plot_x_size, samples);
	const auto y = transform(x, [frequency, amplitude](const double x) { return amplitude * (2 * (x * frequency - floor(1 / 2 + x * frequency)) - 1); });
	plot_function(x, y);
}

void sin_wave(double amplitude, double frequency)
{
	const auto x = linspace(0, plot_x_size, samples);
	const auto y = transform(x, [frequency, amplitude](const double x) { return amplitude * sin(2 * pi * frequency * x); });
	plot_function(x, y);
}

void cos_wave(double amplitude, double frequency)
{
	const auto x = linspace(0, plot_x_size, samples);
	const auto y = transform(x, [frequency, amplitude](const double x) { return amplitude * cos(2 * pi * frequency * x); });
	plot_function(x, y);
}

PYBIND11_MODULE(Signal, m)
{
	m.doc() = "Signal processing";
	m.def("DFT", &create_signal_for_dft, "Transforms and plots sinusoidal signal using DFT");
	m.def("square_wave", &square_wave, "Plots square wave signal");
	m.def("sawtooth_wave", &sawtooth_wave, "Plots sawtooth wave signal");
	m.def("sin_wave", &sin_wave, "Plots sin wave signal");
	m.def("cos_wave", &cos_wave, "Plots cos wave signal");
	m.def("peak", &random_signal, "Finds peak signal in random signal");
	m.def("load_vector", &display_data_from_python, "Loading data from python numpy array to c++ vector");
}
