import sys
from scipy.io import wavfile
import numpy as np
sys.path.append('out/build/x64-debug')
import Signal as sgn
samplerate, data = wavfile.read('Classical.wav')
#sgn.square_wave(1.0, 4.0)
#sgn.DFT(100.0, 1000.0)
#s.sin_wave(1.0, 4.0)
#lista = [1.0,2.1,1.0,2.1,1.0,2.1,1.0,2.1,0.0]
left_channel=[]
#for i in range(0,10000):
#    left_channel.append(data[i][0])

#sgn.load_vector(left_channel,samplerate)