import sys
from scipy.io import wavfile
sys.path.append('out/build/x64-debug')
import Signal as sgn


def load_wav(samples, filename="Classical.wav"):
    samplerate, data = wavfile.read(filename)
    left_channel = []
    for i in range(0, samples):
        left_channel.append(data[i][0])
    sgn.load_vector(left_channel, samplerate)


def menu():
    while True:
        print("1. Visualise .wav file")
        print("2. Visualise sine wave")
        print("3. Visualise cosine wave")
        print("4. Visualise square wave")
        print("5. Visualise sawtooth wave")
        print("6. Generate random signal and find peak value")
        print("7. Discrete Fourier Transform")
        print("8. Exit")

        option = int(input("Option: "))
        if 1 <= option <= 8:
            return option
        else:
            print("Invalid option")


while True:
    match menu():
        case 1:
            load_wav(int(input("Number of samples: ")))
        case 2:
            sgn.sin_wave(float(input("Amplitude: ")), float(input("Frequency: ")))
        case 3:
            sgn.cos_wave(float(input("Amplitude: ")), float(input("Frequency: ")))
        case 4:
            sgn.square_wave(float(input("Amplitude: ")), float(input("Frequency: ")))
        case 5:
            sgn.sawtooth_wave(float(input("Amplitude: ")), float(input("Frequency: ")))
        case 6:
            sgn.peak()
        case 7:
            sgn.DFT(float(input("Amplitude: ")), float(input("Frequency: ")))
        case 8:
            exit(0)
