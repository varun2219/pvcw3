from flask.helpers import url_for
from werkzeug.utils import redirect
from application import app
from flask import render_template, request
import json
import plotly
import scipy
from scipy import signal
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np


@app.route('/')

def home():

    return render_template( 'index.html' )

@app.route( '/index', methods=['POST'] )

def index():
    x = [ str(x) for x in request.form.values( ) ]
    a = x[0]
    a = float(a)

    text_choice  = x[1]
    if(text_choice=='triangular') : 
        
        # triangular wave
        frequency = 2
        time = np.linspace(0,1,500, endpoint=False)
        wave = np.abs(signal.sawtooth(2 * np.pi * frequency * time))
    
        # plot wave
        fig1 = go.Figure(data=go.Scatter(x = time, y = wave))
    
        # generate fourier
        trans_signal = np.zeros_like(wave, dtype=np.complex128)
        wave = wave.copy().astype(np.complex128)
        n = len(wave)
        sf = np.fmod(np.arange(n) + np.fix(n / 2), n).astype(int)
        n_root = np.sqrt(n)
        a = np.remainder(a, 4.0)

        # Special cases
        if a == 0.0:
            return wave
        if a == 2.0:
            return np.flipud(wave)
        if a == 1.0:
            trans_signal[sf] = np.fft.fft(wave[sf]) / n_root
            return trans_signal
        if a == 3.0:
            trans_signal[sf] = np.fft.ifft(wave[sf]) * n_root
            return trans_signal

        # reduce to interval 0.5 < a < 1.5
        if a > 2.0:
            a = a - 2.0
            wave = np.flipud(wave)
        if a > 1.5:
            a = a - 1
            wave[sf] = np.fft.fft(wave[sf]) / n_root
        if a < 0.5:
            a = a + 1
            wave[sf] = np.fft.ifft(wave[sf]) * n_root

        # the general case for 0.5 < a < 1.5
        alpha = a * np.pi / 2
        tana2 = np.tan(alpha / 2)
        sina = np.sin(alpha)
    
        # store output of sin_intercp into dummy
        n = len(wave)
        y = np.zeros(2 * n - 1, dtype=wave.dtype)
        y[:2 * n:2] = wave
        xint = scipy.signal.fftconvolve(
          y[:2 * n],
           np.sinc(np.arange(-(2 * n - 3), (2 * n - 2)).T / 2),
        )
        dummy = xint[2 * n - 3: -2 * n + 3]
        wave = np.hstack((np.zeros(n - 1), dummy, np.zeros(n - 1))).T

        # chirp pre-multiplication
        chrp = np.exp(-1j * np.pi / n * tana2 / 4 *
                  np.arange(-2 * n + 2, 2 * n - 1).T ** 2)
        wave = chrp * wave

        # chirp convolution
        c = np.pi / n / sina / 4
        trans_signal = scipy.signal.fftconvolve(
           np.exp(1j * c * np.arange(-(4 * n - 4), 4 * n - 3).T ** 2),
           wave
        )
        trans_signal = trans_signal[4 * n - 4:8 * n - 7] * np.sqrt(c / np.pi)

        # chirp post multiplication
        trans_signal = chrp * trans_signal

        # normalizing constant
        trans_signal = np.exp(-1j * (1 - a) * np.pi / 4) * trans_signal[n - 1:-n + 1:2]
    
        # extracting real and imaginary part
        real_part = trans_signal.real
        imaginary_part = trans_signal.imag
        df = pd.DataFrame({'time':time, 'real': real_part, 'imag': imaginary_part })
        fig2 = px.line(df, x='time', y=['real', 'imag'])
     
        graph1JSON = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig2, cls=plotly.utils.PlotlyJSONEncoder)
    
        return render_template('index.html', title = "Home", graph1JSON=graph1JSON, graph2JSON=graph2JSON)

    if(text_choice=='sin') : 
        frequency = 2
        time = np.linspace(0,1,500, endpoint=False)
        wave = np.sin(time)
    
        # plot sine wave
        fig1 = go.Figure(data=go.Scatter(x = time, y = wave))
 
        # generate fourier
        trans_signal = np.zeros_like(wave, dtype=np.complex128)
        wave = wave.copy().astype(np.complex128)
        n = len(wave)
        sf = np.fmod(np.arange(n) + np.fix(n / 2), n).astype(int)
        n_root = np.sqrt(n)
        a = np.remainder(a, 4.0)

        # Special cases
        if a == 0.0:
            return wave
        if a == 2.0:
            return np.flipud(wave)
        if a == 1.0:
            trans_signal[sf] = np.fft.fft(wave[sf]) / n_root
            return trans_signal
        if a == 3.0:
            trans_signal[sf] = np.fft.ifft(wave[sf]) * n_root
            return trans_signal

        # reduce to interval 0.5 < a < 1.5
        if a > 2.0:
            a = a - 2.0
            wave = np.flipud(wave)
        if a > 1.5:
            a = a - 1
            wave[sf] = np.fft.fft(wave[sf]) / n_root
        if a < 0.5:
            a = a + 1
            wave[sf] = np.fft.ifft(wave[sf]) * n_root

        # the general case for 0.5 < a < 1.5
        alpha = a * np.pi / 2
        tana2 = np.tan(alpha / 2)
        sina = np.sin(alpha)
    
        # store output of sin_intercp into dummy
        n = len(wave)
        y = np.zeros(2 * n - 1, dtype=wave.dtype)
        y[:2 * n:2] = wave
        xint = scipy.signal.fftconvolve(
          y[:2 * n],
           np.sinc(np.arange(-(2 * n - 3), (2 * n - 2)).T / 2),
        )
        dummy = xint[2 * n - 3: -2 * n + 3]
        wave = np.hstack((np.zeros(n - 1), dummy, np.zeros(n - 1))).T

        # chirp pre-multiplication
        chrp = np.exp(-1j * np.pi / n * tana2 / 4 *
                  np.arange(-2 * n + 2, 2 * n - 1).T ** 2)
        wave = chrp * wave

        # chirp convolution
        c = np.pi / n / sina / 4
        trans_signal = scipy.signal.fftconvolve(
           np.exp(1j * c * np.arange(-(4 * n - 4), 4 * n - 3).T ** 2),
           wave
        )
        trans_signal = trans_signal[4 * n - 4:8 * n - 7] * np.sqrt(c / np.pi)

        # chirp post multiplication
        trans_signal = chrp * trans_signal

        # normalizing constant
        trans_signal = np.exp(-1j * (1 - a) * np.pi / 4) * trans_signal[n - 1:-n + 1:2]
    
        # extracting real and imaginary part
        real_part = trans_signal.real
        imaginary_part = trans_signal.imag
        df = pd.DataFrame({'time':time, 'real': real_part, 'imag': imaginary_part })
        fig2 = px.line(df, x='time', y=['real', 'imag'])
  
        graph1JSON = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig2, cls=plotly.utils.PlotlyJSONEncoder)
    
        return render_template('index.html', title = "Home", graph1JSON=graph1JSON, graph2JSON=graph2JSON)

    if(text_choice == 'rectangle') :
        
        frequency = 2
        time = np.linspace(0,1,500, endpoint=False)
        wave = signal.square(2 * np.pi * frequency * time)
    
        # plot rectangle wave
        fig1 = go.Figure(data=go.Scatter(x = time, y = wave))
 
        # generate fourier
        trans_signal = np.zeros_like(wave, dtype=np.complex128)
        wave = wave.copy().astype(np.complex128)
        n = len(wave)
        sf = np.fmod(np.arange(n) + np.fix(n / 2), n).astype(int)
        n_root = np.sqrt(n)
        a = np.remainder(a, 4.0)

        # Special cases
        if a == 0.0:
            return wave
        if a == 2.0:
            return np.flipud(wave)
        if a == 1.0:
            trans_signal[sf] = np.fft.fft(wave[sf]) / n_root
            return trans_signal
        if a == 3.0:
            trans_signal[sf] = np.fft.ifft(wave[sf]) * n_root
            return trans_signal

        # reduce to interval 0.5 < a < 1.5
        if a > 2.0:
            a = a - 2.0
            wave = np.flipud(wave)
        if a > 1.5:
            a = a - 1
            wave[sf] = np.fft.fft(wave[sf]) / n_root
        if a < 0.5:
            a = a + 1
            wave[sf] = np.fft.ifft(wave[sf]) * n_root

        # the general case for 0.5 < a < 1.5
        alpha = a * np.pi / 2
        tana2 = np.tan(alpha / 2)
        sina = np.sin(alpha)
    
        # store output of sin_intercp into dummy
        n = len(wave)
        y = np.zeros(2 * n - 1, dtype=wave.dtype)
        y[:2 * n:2] = wave
        xint = scipy.signal.fftconvolve(
          y[:2 * n],
           np.sinc(np.arange(-(2 * n - 3), (2 * n - 2)).T / 2),
        )
        dummy = xint[2 * n - 3: -2 * n + 3]
        wave = np.hstack((np.zeros(n - 1), dummy, np.zeros(n - 1))).T

        # chirp pre-multiplication
        chrp = np.exp(-1j * np.pi / n * tana2 / 4 *
                  np.arange(-2 * n + 2, 2 * n - 1).T ** 2)
        wave = chrp * wave

        # chirp convolution
        c = np.pi / n / sina / 4
        trans_signal = scipy.signal.fftconvolve(
           np.exp(1j * c * np.arange(-(4 * n - 4), 4 * n - 3).T ** 2),
           wave
        )
        trans_signal = trans_signal[4 * n - 4:8 * n - 7] * np.sqrt(c / np.pi)

        # chirp post multiplication
        trans_signal = chrp * trans_signal

        # normalizing constant
        trans_signal = np.exp(-1j * (1 - a) * np.pi / 4) * trans_signal[n - 1:-n + 1:2]
    
        # extracting real and imaginary part
        real_part = trans_signal.real
        imaginary_part = trans_signal.imag
        df = pd.DataFrame({'time':time, 'real': real_part, 'imag': imaginary_part })
        fig2 = px.line(df, x='time', y=['real', 'imag'])
   
        graph1JSON = json.dumps(fig1, cls=plotly.utils.PlotlyJSONEncoder)
        graph2JSON = json.dumps(fig2, cls=plotly.utils.PlotlyJSONEncoder)
    
        return render_template('index.html', title = "Home", graph1JSON=graph1JSON, graph2JSON=graph2JSON)