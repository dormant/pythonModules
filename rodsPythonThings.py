# Various plots

import matplotlib.pyplot as plt
import numpy as np
import obspy
import rodsPlotTfr as rpt
import scipy.ndimage as ndimage
from obspy.core import UTCDateTime, Stream, Trace
from obspy.signal import util
from scipy.signal import hilbert
from obspy.imaging.cm import obspy_sequential, obspy_divergent
from matplotlib.ticker import NullFormatter
from scipy import signal
from scipy.fft import fftshift
import math
import statistics as stat

plt.rcParams.update({'font.size': 16})



def plotStringThing( st, pre, dur ):

    tr1 = st[0]
    tr2 = st[1]

    # Plot
    fig, ax1 = plt.subplots()

    data = tr1.data
    data = np.abs( data )
    data = np.sqrt( data )
    npts = tr1.stats.npts
    dt = tr1.stats.delta
    samprate = tr1.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre
    ax1.plot(t, data, color=(0.8,0.8,0.8) )
    ax1.xaxis.set_major_formatter(NullFormatter())
    ax1.yaxis.set_major_formatter(NullFormatter())
    ax1.set_xlim(-1*pre, dur-pre)
    ax1.set_ylim(ymin=0)

    data = tr2.data
    npts = tr2.stats.npts
    dt = tr2.stats.delta
    samprate = tr2.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre
    ax2 = ax1.twinx()
    ax2.plot(t, data, color='r')
    ax2.xaxis.set_major_formatter(NullFormatter())
    ax2.yaxis.set_major_formatter(NullFormatter())
    ax2.set_xlim(-1*pre, dur-pre)

    plt.subplots_adjust(left=0.01, right=0.99, top=0.8, bottom=0.05)

    return fig


def plotZforai( st, pre, plotFmin=2.0, plotFmax=100 ):

    tr = st[0]

    # Data
    data = tr.data
    npts = tr.stats.npts
    dt = tr.stats.delta
    samprate = tr.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre

    # Spectrogram
    spec = calcSgramSt( st )

    # Plot
    cmap='jet'
    fig = plt.figure()
    ax = plt.gca()
    ax.pcolormesh( ttt, fff, sss, shading='gouraud', cmap=cmap )
    ax.set_yscale( 'log' )
    ax.set_ylim(plotFmin, plotFmax)
    #ax.set_xlim( tLimits[0]+pre+0.1, tLimits[1]+pre )
    ax.xaxis.set_major_formatter(NullFormatter())

    return fig



def plotZSpectrum( st, pre, plotSpec, plotFscale, plotZscale, plotFmin=0.5, plotFmax=100 ):

    tr = st[0]

    # Data
    data = tr.data
    npts = tr.stats.npts
    dt = tr.stats.delta
    samprate = tr.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre

    # Calculate fft
    p = 20*np.log10(np.abs(np.fft.rfft(data)))
    f = np.linspace(0, samprate/2, len(p))

    # Calculate power spectrum
    #ps = np.abs(np.fft.fft(data))**2
    #freqs = np.fft.fftfreq(data.size, dt )
    #idx = np.argsort(freqs)


    # Positions for subplots
    x1 = 0.1
    x2 = 0.5
    x3 = 0.9
    y1 = 0.1
    y2 = 0.9

    fig = plt.figure()

    # Plot signal
    ax_sig = fig.add_axes([ x2, y1, x3-x2-0.02, y2-y1 ] )
    ax_sig.plot( t, data, 'k', linewidth=0.5 )
    datamax = max(max(data),-min(data))
    ax_sig.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig.set_xlim(t[0], t[-1])
    ax_sig.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    tLabel = '.'.join( ['MV',tr.stats.station,tr.stats.location,tr.stats.channel] )
    ax_sig.text( xLabel, yLabel, tLabel )
    #ax_sig.xaxis.tick_top()
    yLimits = ax_sig.get_ylim()
    ax_sig.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig.vlines( 0.0, yLimits[0], yLimits[1], colors='g', linewidth=0.2 )

    # Plot spectrum
    ax_sig = fig.add_axes([ x1, y1, x2-x1-0.02, y2-y1 ] )
    ax_sig.loglog( f, p, 'k', linewidth=0.5 )
    ax_sig.set_xlim(1, 100)
    ax_sig.set_ylim(10, 1000)
    ax_sig.yaxis.set_major_formatter(NullFormatter())
    #ax_sig.xaxis.tick_top()
    yLimits = ax_sig.get_ylim()
    ax_sig.grid( visible=True, which='major', axis='both', linewidth=0.1 )
    ax_sig.vlines( 0.0, yLimits[0], yLimits[1], colors='g', linewidth=0.2 )


    return fig



def plotZManyWays( st, pre, plotSpec, plotFscale, plotZscale, plotFmin=0.5, plotFmax=100 ):

    tr = st[0]

    # Data
    data = tr.data
    npts = tr.stats.npts
    dt = tr.stats.delta
    samprate = tr.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre

    # Envelope
    tr_filt = tr.copy()
    tr_filt.filter( 'highpass', freq=1.0, corners=2, zerophase=True )
    data_envelope = obspy.signal.filter.envelope(tr_filt.data)
    kernel_size = 50
    kernel = np.ones(kernel_size) / kernel_size
    data_envelope = np.convolve(data_envelope, kernel, mode='same')

    # Spectrogram
    fft_zero_pad_fac = 4
    #nf =  128
    nf =  64
    w0 = 16.0
    if plotFmax > fnyq:
        plotFmax = fnyq
    cmap='jet'
    plot_args=['k', 'k']
    if fft_zero_pad_fac == 0:
        nfft = npts
    else:
        nfft = util.next_pow_2(npts) * fft_zero_pad_fac
    f_lin = np.linspace(0, 0.5 / dt, nfft // 2 + 1)
    _w = np.zeros((1, nf, npts), dtype=complex)
    _w[0] = rpt.cwt(data, dt, w0, plotFmin, plotFmax, nf)
    ntr = 1
    spec = np.zeros((1, nfft // 2 + 1), dtype=complex)
    spec[0] = np.fft.rfft(data, n=nfft) * dt

    # Spectrogram 2
    #fff, ttt, Sxx = signal.spectrogram( data, samprate )
    #if mode == 'sqrt':
    #    Sxx = np.sqrt( np.abs(Sxx) )


    if plotZscale == 'amp':
        _tfr = np.abs(_w)
        spec = np.abs(spec)
    elif plotZscale == 'power':
        _tfr = np.abs(_w) ** 2
        spec = np.abs(spec) ** 2
    elif plotZscale == 'log':
        _tfr = np.log( np.abs(_w) )
        spec = np.log( np.abs(spec) )
    elif plotZscale == 'sqrt':
        _tfr = np.sqrt( np.abs(_w) )
        spec = np.sqrt( np.abs(spec) )


    # Positions for subplots
    x1 = 0.1
    if plotSpec:
        x2 = 0.7
        x3 = 0.9
    else:
        x2 = 0.9
    y1 = 0.1
    y2 = 0.3
    y3 = 0.7
    y4 = 0.9

    fig = plt.figure()

    ax_sig = fig.add_axes([ x1, y3, x2-x1, y4-y3 ] )
    ax_sig.plot( t, data, 'k', linewidth=0.5 )
    datamax = max(max(data),-min(data))
    ax_sig.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig.set_xlim(t[0], t[-1])
    ax_sig.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    tLabel = '.'.join( ['MV',tr.stats.station,tr.stats.location,tr.stats.channel] )
    ax_sig.text( xLabel, yLabel, tLabel )
    ax_sig.xaxis.tick_top()
    yLimits = ax_sig.get_ylim()
    ax_sig.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig.vlines( 0.0, yLimits[0], yLimits[1], colors='g', linewidth=0.2 )

    ax_env = fig.add_axes([ x1, y1, x2-x1, y2-y1 ] )
    ax_env.plot( t, data_envelope, 'k', linewidth=0.5 )
    datamax = max(data_envelope)
    ax_env.set_ylim( 0.0, 1.05 * datamax)
    ax_env.set_xlim(t[0], t[-1])
    ax_env.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.1*datamax
    tLabel = 'Envelope'
    ax_env.set( xlabel='Time (seconds)' )
    ax_env.text( xLabel, yLabel, tLabel )
    yLimits = ax_env.get_ylim()
    ax_env.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_env.vlines( 0.0, yLimits[0], yLimits[1], colors='g', linewidth=0.2 )

    ax_tfr = fig.add_axes([ x1, y2, x2-x1, y3-y2 ])
    x, y = np.meshgrid(
            t, np.logspace(np.log10(plotFmin), np.log10(plotFmax),
            _tfr[0].shape[0]))
    img_tfr = rpt._pcolormesh_same_dim(ax_tfr, x, y, _tfr[0], cmap=cmap)
    img_tfr.set_rasterized(True)
    ax_tfr.set_yscale( plotFscale )
    if plotFscale == 'linear':
        ax_tfr.set_ylim(0.0, plotFmax-1)
    elif plotFscale == 'log':
        ax_tfr.set_ylim(plotFmin, plotFmax)
    ax_tfr.set_ylim(plotFmin, 100.0)
    ax_tfr.set_xlim(t[0], t[-1])
    ax_tfr.xaxis.set_major_formatter(NullFormatter())
    ax_tfr.set( ylabel='Frequency (Hz)' )
    clim = _tfr.max()
    img_tfr.set_clim(0., clim)
    #ax_tfr.grid( visible=True, which='major', axis='both', linewidth=0.2, color='w' )

    # For spectrogram 2
    #ax_sgram = fig.add_axes([ x1, y2, x2-x1, y3-y2 ])
    #ax_sgram.pcolormesh( ttt, fff, Sxx, shading='gouraud', cmap=cmap )
    #ax_sgram.set_yscale("log")
    #ax_sgram.set_ylim(fmin, 100.0)
    #ax_sgram.set_xlim(ttt[0], ttt[-1])
    #ax_sgram.xaxis.set_major_formatter(NullFormatter())
    #ax_sgram.set( ylabel='Frequency (Hz)' )


    if plotSpec:
        ax_spec = fig.add_axes([x2, y2, x3-x2, y3-y2])
        ax_spec.semilogy(spec[0], f_lin, plot_args[1], linewidth=0.5)
        ax_spec.set_ylim(plotFmin, 100.0)
        ax_spec.set_xlim(0.0, max(spec[0]))
        ax_spec.xaxis.set_major_formatter(NullFormatter())
        ax_spec.yaxis.tick_right()
        ax_spec.grid( visible=True, which='major', axis='both', linewidth=0.1 )

    return fig



def plot3CPartMot( st, pre, twin ):

    for itr in range(3):
        tr = st[itr]
        if tr.stats.channel[-1] == "Z":
            tr1 = tr
        elif tr.stats.channel[-1] == "N" or tr.stats.channel[-1] == "1":
            tr2 = tr
        elif tr.stats.channel[-1] == "E" or tr.stats.channel[-1] == "2":
            tr3 = tr

    data1 = tr1.data
    data2 = tr2.data
    data3 = tr3.data
    npts = tr1.stats.npts
    dt = tr1.stats.delta
    samprate = tr1.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre
    tlim1plus = 0.1

    mask = ~((t >= 0.0) & (t <= twin))
    unmask = np.logical_not(mask)

    # Positions for subplots
    x1 = 0.1
    x2 = 0.7
    x3 = 0.9
    y1 = 0.09
    y2 = 0.36
    y3 = 0.63
    y4 = 0.9

    fig = plt.figure()

    data1max = max(max(data1),-min(data1))
    data2max = max(max(data2),-min(data2))
    data3max = max(max(data3),-min(data3))
    datamax = max(data1max,data2max,data3max)

    ax_sig1 = fig.add_axes([ x1, y3, x2-x1, y4-y3 ] )
    ax_sig1.plot( np.ma.array(t, mask=unmask), np.ma.array(data1, mask=unmask), 'k', linewidth=0.5 )
    ax_sig1.plot( np.ma.array(t, mask=mask), np.ma.array(data1, mask=mask), 'r', linewidth=0.5 )
    ax_sig1.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig1.set_xlim(t[0]+tlim1plus, t[-1])
    ax_sig1.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    t1Label = '.'.join( ['MV',tr1.stats.station,tr1.stats.location,tr1.stats.channel] )
    ax_sig1.text( xLabel, yLabel, t1Label )
    ax_sig1.xaxis.tick_top()
    yLimits = ax_sig1.get_ylim()
    ax_sig1.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig1.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )

    ax_sig2 = fig.add_axes([ x1, y2, x2-x1, y3-y2 ] )
    ax_sig2.plot( np.ma.array(t, mask=unmask), np.ma.array(data2, mask=unmask), 'k', linewidth=0.5 )
    ax_sig2.plot( np.ma.array(t, mask=mask), np.ma.array(data2, mask=mask), 'r', linewidth=0.5 )
    ax_sig2.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig2.set_xlim(t[0]+tlim1plus, t[-1])
    ax_sig2.xaxis.set_major_formatter(NullFormatter())
    ax_sig2.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    t1Label = '.'.join( ['MV',tr1.stats.station,tr2.stats.location,tr2.stats.channel] )
    ax_sig2.text( xLabel, yLabel, t1Label )
    yLimits = ax_sig2.get_ylim()
    ax_sig2.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig2.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )

    ax_sig3 = fig.add_axes([ x1, y1, x2-x1, y2-y1 ] )
    ax_sig3.plot( np.ma.array(t, mask=unmask), np.ma.array(data3, mask=unmask), 'k', linewidth=0.5 )
    ax_sig3.plot( np.ma.array(t, mask=mask), np.ma.array(data3, mask=mask), 'r', linewidth=0.5 )
    ax_sig3.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig3.set_xlim(t[0]+tlim1plus, t[-1])
    ax_sig3.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    t1Labe3 = '.'.join( ['MV',tr1.stats.station,tr3.stats.location,tr3.stats.channel] )
    ax_sig3.text( xLabel, yLabel, t1Label )
    yLimits = ax_sig3.get_ylim()
    ax_sig3.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig3.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )

    ax_pm1 = fig.add_axes([ x2, y3, x3-x2, y4-y3 ] )
    ax_pm1.set_aspect('equal', adjustable='box')
    ax_pm1.plot( np.ma.array(data2,mask=mask), np.ma.array(data1,mask=mask), 'r', linewidth=0.5 )
    ax_pm1.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_pm1.set_xlim( -1.05 * datamax, 1.05 * datamax)
    ax_pm1.yaxis.set_major_formatter(NullFormatter())
    ax_pm1.xaxis.set_major_formatter(NullFormatter())
    xLabel = -1 * datamax
    yLabel = datamax - 0.2*datamax
    tLabel = tr1.stats.channel
    ax_pm1.text( xLabel, yLabel, tLabel, va='bottom', ha='left' )
    xLabel = datamax
    yLabel = -1 * datamax + 0.2*datamax
    tLabel = tr2.stats.channel
    ax_pm1.text( xLabel, yLabel, tLabel, va='top', ha='right' )
    ax_pm1.grid( visible=True, which='major', axis='both', linewidth=0.1 )
    ax_pm1.xaxis.tick_top()

    ax_pm2 = fig.add_axes([ x2, y2, x3-x2, y3-y2 ] )
    ax_pm2.set_aspect('equal', adjustable='box')
    ax_pm2.plot( np.ma.array(data3,mask=mask), np.ma.array(data1,mask=mask), 'r', linewidth=0.5 )
    ax_pm2.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_pm2.set_xlim( -1.05 * datamax, 1.05 * datamax)
    ax_pm2.yaxis.set_major_formatter(NullFormatter())
    ax_pm2.xaxis.set_major_formatter(NullFormatter())
    xLabel = -1 * datamax
    yLabel = datamax - 0.2*datamax
    tLabel = tr1.stats.channel
    ax_pm2.text( xLabel, yLabel, tLabel, va='bottom', ha='left' )
    xLabel = datamax
    yLabel = -1 * datamax + 0.2*datamax
    tLabel = tr3.stats.channel
    ax_pm2.text( xLabel, yLabel, tLabel, va='top', ha='right' )
    ax_pm2.grid( visible=True, which='major', axis='both', linewidth=0.1 )

    ax_pm3 = fig.add_axes([ x2, y1, x3-x2, y2-y1 ] )
    ax_pm3.set_aspect('equal', adjustable='box')
    ax_pm3.plot( np.ma.array(data2,mask=mask), np.ma.array(data3,mask=mask), 'r', linewidth=0.5 )
    ax_pm3.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_pm3.set_xlim( -1.05 * datamax, 1.05 * datamax)
    ax_pm3.yaxis.set_major_formatter(NullFormatter())
    ax_pm3.xaxis.set_major_formatter(NullFormatter())
    xLabel = -1 * datamax
    yLabel = datamax - 0.2*datamax
    tLabel = tr3.stats.channel
    ax_pm3.text( xLabel, yLabel, tLabel, va='bottom', ha='left' )
    xLabel = datamax
    yLabel = -1 * datamax + 0.2*datamax
    tLabel = tr2.stats.channel
    ax_pm3.text( xLabel, yLabel, tLabel, va='top', ha='right' )
    ax_pm3.grid( visible=True, which='major', axis='both', linewidth=0.1 )

    return fig



def plot3CManyWays( st, pre, plotFscale, plotZscale, plotFmin=0.5, plotFmax=100.0 ):

    for itr in range(3):
        tr = st[itr]
        if tr.stats.channel[-1] == "Z":
            tr1 = tr
        elif tr.stats.channel[-1] == "N" or tr.stats.channel[-1] == "1":
            tr2 = tr
        elif tr.stats.channel[-1] == "E" or tr.stats.channel[-1] == "2":
            tr3 = tr

    data1 = tr1.data
    data2 = tr2.data
    data3 = tr3.data
    npts = tr1.stats.npts
    dt = tr1.stats.delta
    samprate = tr1.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre
    cmap='jet'
    tlim1plus = 0.1

    # Envelope
    kernel_size = 50
    kernel = np.ones(kernel_size) / kernel_size
    tr_filt = tr1.copy()
    tr_filt.filter( 'highpass', freq=1.0, corners=2, zerophase=True )
    data_envelope = obspy.signal.filter.envelope(tr_filt.data)
    data1_envelope = np.convolve(data_envelope, kernel, mode='same')
    tr_filt = tr2.copy()
    tr_filt.filter( 'highpass', freq=1.0, corners=2, zerophase=True )
    data_envelope = obspy.signal.filter.envelope(tr_filt.data)
    data2_envelope = np.convolve(data_envelope, kernel, mode='same')
    tr_filt = tr3.copy()
    tr_filt.filter( 'highpass', freq=1.0, corners=2, zerophase=True )
    data_envelope = obspy.signal.filter.envelope(tr_filt.data)
    data3_envelope = np.convolve(data_envelope, kernel, mode='same')

    # Positions for subplots
    x1 = 0.1
    x2 = 0.5
    x3 = 0.9
    y1 = 0.1
    y2 = 0.3
    y3 = 0.5
    y4 = 0.7
    y5 = 0.9

    fig = plt.figure()

    data1max = max(max(data1),-min(data1))
    data2max = max(max(data2),-min(data2))
    data3max = max(max(data3),-min(data3))
    datamax = max(data1max,data2max,data3max)

    ax_sig1 = fig.add_axes([ x1, y4, x2-x1, y5-y4 ] )
    ax_sig1.plot( t, data1, 'r', linewidth=0.5 )
    ax_sig1.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig1.set_xlim(t[0]+tlim1plus, t[-1])
    ax_sig1.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    t1Label = '.'.join( ['MV',tr1.stats.station,tr1.stats.location,tr1.stats.channel] )
    ax_sig1.text( xLabel, yLabel, t1Label )
    ax_sig1.xaxis.tick_top()
    yLimits = ax_sig1.get_ylim()
    ax_sig1.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig1.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )

    ax_sig2 = fig.add_axes([ x1, y3, x2-x1, y4-y3 ] )
    ax_sig2.plot( t, data2, 'g', linewidth=0.5 )
    ax_sig2.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig2.set_xlim(t[0]+tlim1plus, t[-1])
    ax_sig2.xaxis.set_major_formatter(NullFormatter())
    ax_sig2.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    t2Label = '.'.join( ['MV',tr2.stats.station,tr2.stats.location,tr2.stats.channel] )
    ax_sig2.text( xLabel, yLabel, t2Label )
    yLimits = ax_sig2.get_ylim()
    ax_sig2.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig2.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )

    ax_sig3 = fig.add_axes([ x1, y2, x2-x1, y3-y2 ] )
    ax_sig3.plot( t, data3, 'b', linewidth=0.5 )
    ax_sig3.set_ylim( -1.05 * datamax, 1.05 * datamax)
    ax_sig3.set_xlim(t[0]+tlim1plus, t[-1])
    ax_sig3.xaxis.set_major_formatter(NullFormatter())
    ax_sig3.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.2*datamax
    t3Label = '.'.join( ['MV',tr3.stats.station,tr3.stats.location,tr3.stats.channel] )
    ax_sig3.text( xLabel, yLabel, t3Label )
    yLimits = ax_sig3.get_ylim()
    ax_sig3.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_sig3.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )

    ax_env = fig.add_axes([ x1, y1, x2-x1, y2-y1 ] )
    datamax = max( max(data1_envelope), max(data2_envelope), max(data3_envelope) )
    ax_env.plot( t, data1_envelope, 'r', linewidth=0.5 )
    ax_env.plot( t, data2_envelope, 'g', linewidth=0.5 )
    ax_env.plot( t, data3_envelope, 'b', linewidth=0.5 )
    ax_env.set_ylim( 0.0, 1.05 * datamax)
    ax_env.set_xlim(t[0]+tlim1plus, t[-1])
    ax_env.yaxis.set_major_formatter(NullFormatter())
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = datamax - 0.1*datamax
    eLabel = 'Envelope'
    ax_env.text( xLabel, yLabel, eLabel )
    yLimits = ax_env.get_ylim()
    ax_env.grid( visible=True, which='major', axis='x', linewidth=0.1 )
    ax_env.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )
    ax_env.set( xlabel='Time (seconds)' )

    ax_envr = fig.add_axes([ x2, y1, x3-x2, y2-y1 ] )
    data_env_hor = np.add( np.square(data2_envelope), np.square(data3_envelope) )
    data_env_hor = np.sqrt(data_env_hor)
    data_env_ratio = np.divide( data1_envelope, data_env_hor )
    ax_envr.semilogy( t, data_env_ratio, 'k', linewidth=0.5 )
    ax_envr.set_xlim(t[0]+tlim1plus, t[-1])
    ax_envr.set_ylim(0.08,12.0)
    xLabel = t[0] + 0.01*(t[-1]-t[0])
    yLabel = 5.0
    eLabel = 'Envelope ratio (V/H)'
    ax_envr.text( xLabel, yLabel, eLabel )
    yLimits = ax_envr.get_ylim()
    ax_envr.grid( visible=True, which='major', axis='both', color='g', linewidth=0.1 )
    ax_envr.vlines( 0.0, yLimits[0], yLimits[1], colors='g', linewidth=0.2 )
    ax_envr.yaxis.tick_right()
    ax_envr.set( xlabel='Time (seconds)' )

    tLimits = ax_envr.get_xlim()

    fff, ttt, sss = calcSgram( data1, samprate, plotFscale )
    if plotZscale =='sqrt':
        sss = np.sqrt( np.abs( sss ) )
    elif plotZscale =='log':
        sss = np.log( np.abs( sss ) )
    elif plotZscale =='power':
        sss = np.abs( sss ) ** 2
    elif plotZscale =='amp':
        sss = np.abs( sss ) 
    ax_sgram1 = fig.add_axes([ x2, y4, x3-x2, y5-y4 ] )
    ax_sgram1.pcolormesh( ttt, fff, sss, shading='gouraud', cmap=cmap )
    ax_sgram1.set_yscale( plotFscale )
    if plotFscale == 'linear':
        ax_sgram1.set_ylim(0.0, plotFmax-1)
    elif plotFscale == 'log':
        ax_sgram1.set_ylim(plotFmin, plotFmax)
    ax_sgram1.set_xlim( tLimits[0]+pre+0.1, tLimits[1]+pre )
    ax_sgram1.xaxis.set_major_formatter(NullFormatter())
    ax_sgram1.yaxis.tick_right()
    #ax_sgram1.grid( visible=True, color='w' )

    fff, ttt, sss = calcSgram( data2, samprate, plotFscale )
    if plotZscale =='sqrt':
        sss = np.sqrt( np.abs( sss ) )
    elif plotZscale =='log':
        sss = np.log( np.abs( sss ) )
    elif plotZscale =='amp':
        sss = np.abs( sss )
    elif plotZscale =='sqrt':
        sss = np.sqrt( np.abs( sss ) )
    ax_sgram2 = fig.add_axes([ x2, y3, x3-x2, y4-y3 ] )
    ax_sgram2.pcolormesh( ttt, fff, sss, shading='gouraud', cmap=cmap )
    ax_sgram2.set_yscale( plotFscale )
    if plotFscale == 'linear':
        ax_sgram2.set_ylim(0.0, plotFmax-1)
    elif plotFscale == 'log':
        ax_sgram2.set_ylim(plotFmin, plotFmax)
    ax_sgram2.set_xlim( tLimits[0]+pre+0.1, tLimits[1]+pre )
    ax_sgram2.xaxis.set_major_formatter(NullFormatter())
    ax_sgram2.yaxis.tick_right()
    ax_sgram2.yaxis.set_label_position("right")
    ax_sgram2.set( ylabel='Frequency (Hz)' )
    ax_sgram2.spines['top'].set_color('w')
    #ax_sgram2.grid( visible=True, color='w' )

    fff, ttt, sss = calcSgram( data3, samprate, plotFscale )
    if plotZscale =='sqrt':
        sss = np.sqrt( np.abs( sss ) )
    elif plotZscale =='log':
        sss = np.log( np.abs( sss ) )
    elif plotZscale =='amp':
        sss = np.abs( sss )
    elif plotZscale =='sqrt':
        sss = np.sqrt( np.abs( sss ) )
    ax_sgram3 = fig.add_axes([ x2, y2, x3-x2, y3-y2 ] )
    ax_sgram3.pcolormesh( ttt, fff, sss, shading='gouraud', cmap=cmap )
    ax_sgram3.set_yscale( plotFscale )
    if plotFscale == 'linear':
        ax_sgram3.set_ylim(0.0, plotFmax-1)
    elif plotFscale == 'log':
        ax_sgram3.set_ylim(plotFmin, plotFmax)
    ax_sgram3.set_xlim( tLimits[0]+pre+0.1, tLimits[1]+pre )
    ax_sgram3.xaxis.set_major_formatter(NullFormatter())
    ax_sgram3.yaxis.tick_right()
    ax_sgram3.spines['top'].set_color('w')
    #ax_sgram3.grid( visible=True, color='w' )

    return fig


def plotLahar( st, pre, plotRms, plotFscale, plotZscale, plotFmin=1.0, plotFmax=50.0 ):

    plt.rcParams.update({'font.size': 14})

    nTrace = len( st )

    stRMS = streamRMS( st, 0.1 )

    # Positions for subplots
    x1 = 0.1
    x2 = 0.9
    nYPanes = 2 * nTrace
    y = []
    for iy in range(nYPanes):
        y.append( 0.8*iy/nYPanes + 0.1 )
    y.append( 0.9)

    fig = plt.figure()

    iPane = nYPanes - 1
    for iTrace in range( nTrace ):

        ax = fig.add_axes( [ x1, y[iPane], x2-x1, y[iPane+1]-y[iPane] ] )
        iPane -=  1
        tr = st[iTrace]
        data = tr.data
        npts = tr.stats.npts
        dt = tr.stats.delta
        samprate = tr.stats.sampling_rate
        fnyq = samprate/2.0
        t = np.arange(0, npts / samprate, 1 / samprate)
        t = t - pre
        cmap='jet'
        if plotFmax > fnyq:
            plotFmax = fnyq
        tlim1plus = 0.1
        datamax = max(max(data),-min(data))

        if plotRms:
            ax.plot( t, data, '0.5', linewidth=0.5 )
        else:
            ax.plot( t, data, 'k', linewidth=0.5 )
        ax.set_ylim( -1.05 * datamax, 1.05 * datamax)
        ax.set_xlim(t[0]+tlim1plus, t[-1])
        ax.yaxis.set_major_formatter(NullFormatter())
        xLabel = t[0] + 0.01*(t[-1]-t[0])
        yLabel = datamax - 0.3*datamax
        tLabel = '.'.join( ['MV',tr.stats.station,tr.stats.location,tr.stats.channel] )
        ax.text( xLabel, yLabel, tLabel )
        if plotRms:
            xLabel = t[-1] - 0.01*(t[-1]-t[0])
            tLabel = 'RMS'
            ax.text( xLabel, yLabel, tLabel, color='r', ha='right' )
        ax.xaxis.tick_top()
        yLimits = ax.get_ylim()
        ax.grid( visible=True, which='major', axis='x', linewidth=0.1 )
        ax.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )
        ax.tick_params(axis='x', direction='in')
        if iTrace == nTrace -1:
            ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=True, left=False, right=False)
        else:
            ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=True, left=False, right=False)
        tLimits = ax.get_xlim()

        # Add RMS to same axes
        if plotRms:
            ax2 = ax.twinx()
            tr2 = stRMS[iTrace]
            data2 = tr2.data
            npts2 = tr2.stats.npts
            samprate2 = tr2.stats.sampling_rate
            t2 = np.arange(0, npts2 / samprate2, 1 / samprate2)
            t2 = t2 - pre
            datamax2 = max(max(data2),-min(data2))
            ax2.plot( t2, data2, 'r', linewidth=1.0 )
            ax2.set_ylim( 0.0, 1.05 * datamax2)
            #ax2.set_ylim( -1.05 * datamax2, 1.05 * datamax2)
            ax2.set_xlim(t[0]+tlim1plus, t[-1])
            ax2.yaxis.set_major_formatter(NullFormatter())

        ax = fig.add_axes( [ x1, y[iPane], x2-x1, y[iPane+1]-y[iPane] ] )
        fff, ttt, Sxx = signal.spectrogram( data, samprate )
        if plotZscale == 'sqrt':
            Sxx = np.sqrt( np.abs(Sxx) )
        elif plotZscale =='log':
            Sxx = np.log( np.abs( Sxx ) )
        elif plotZscale =='power':
            Sxx = np.abs( Sxx ) ** 2
        elif plotZscale =='amp':
            Sxx = np.abs( Sxx ) 
        ax.pcolormesh( ttt, fff, Sxx, shading='gouraud', cmap=cmap )
        ax.set_yscale( plotFscale )
        if plotFscale == 'linear':
            ax.set_ylim(0.0, plotFmax )
        elif plotFscale == 'log':
            ax.set_ylim(plotFmin, plotFmax )
        ax.set_xlim( tLimits[0]+pre+0.1, tLimits[1]+pre )
        if iTrace == nTrace -1:
            ax.tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=True, left=False, right=False)
            ax.set( xlabel='Time (seconds)' )
        else:
            ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=True, left=False, right=False)

        iPane -=  1

    return fig

def plotRockfall( st, pre, datimEventString ):

    plt.rcParams.update({'font.size': 14})

    nTrace = len( st )

    st.filter( "bandpass", freqmin=2.0, freqmax=20.0, corners=8, zerophase=True )
    #st.filter( "bandpass", freqmin=2.0, freqmax=10.0, corners=8, zerophase=True )

    stRMS = streamRMS( st, 0.04 )
    #stRMS = streamRMS( st, 0.03 )

    # Positions for subplots
    tlim1plus = 0.1
    x1 = 0.1
    x2 = 0.9
    y = []
    for iy in range(nTrace):
        y.append( 0.8*iy/nTrace + 0.1 )
    y.append( 0.9)

    fig = plt.figure()
    fileRockfallAmps = ''.join([ datimEventString, '--Rockfall-amplitudes.txt' ])
    file1 = open(fileRockfallAmps, "w") 

    iPane = nTrace - 1
    for iTrace in range( nTrace ):

        ax = fig.add_axes( [ x1, y[iPane], x2-x1, y[iPane+1]-y[iPane] ] )
        iPane -=  1
        tr = st[iTrace]
        data = tr.data
        npts = tr.stats.npts
        dt = tr.stats.delta
        samprate = tr.stats.sampling_rate
        fnyq = samprate/2.0
        t = np.arange(0, npts / samprate, 1 / samprate)
        t = t - pre
        datamax = max(max(data),-min(data))

        #clipped = isClipped( tr )

        ax.plot( t, data, '0.5', linewidth=0.5 )
        ax.set_ylim( -1.05 * datamax, 1.05 * datamax)
        ax.set_xlim(t[0]+tlim1plus, t[-1])
        ax.yaxis.set_major_formatter(NullFormatter())
        xLabel = t[0] + 0.01*(t[-1]-t[0])
        yLabel = datamax - 0.3*datamax
        tLabel = '.'.join( ['MV',tr.stats.station,tr.stats.location,tr.stats.channel] )
        ax.text( xLabel, yLabel, tLabel )
        xLabel = t[-1] - 0.01*(t[-1]-t[0])
        tLabel = 'RMS'
        ax.text( xLabel, yLabel, tLabel, color='r', ha='right' )
        ax.xaxis.tick_top()
        yLimits = ax.get_ylim()
        ax.grid( visible=True, which='major', axis='x', linewidth=0.1 )
        ax.vlines( 0.0, yLimits[0], yLimits[1], colors='k', linewidth=0.2 )
        ax.tick_params(axis='x', direction='in')
        if iTrace == nTrace -1:
            ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=True, left=False, right=False)
        else:
            ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=True, left=False, right=False)
        tLimits = ax.get_xlim()

        # Add RMS to same axes
        if npts > 0:
            ax2 = ax.twinx()
            tr2 = stRMS[iTrace]
            data2 = tr2.data
            npts2 = tr2.stats.npts
            samprate2 = tr2.stats.sampling_rate
            t2 = np.arange(0, npts2 / samprate2, 1 / samprate2)
            t2 = t2 - pre
            t2 = t2 + (0.5/samprate2)
            data2max = max(data2)
            whereDataMax = np.where( data2 == data2max )[0][0]
            data2Pre = data2[t2 < 0]
            data2Pos = data2[t2 >= 0]
            data2meanPre = np.mean( data2Pre )
            data2meanPos = np.mean( data2Pos )
            datamed2 = stat.median(data2)

            if data2meanPos > data2meanPre:
                ax2.plot( t2, data2, 'b', linewidth=1.0 )
            else:
                whereDataMax = 0
                ax2.plot( t2, data2, 'r', linewidth=1.0 )

            ax2.set_ylim( 0.0, 1.05 * data2max)
            ax2.set_xlim(t[0]+tlim1plus, t[-1])
            ax2.yaxis.set_major_formatter(NullFormatter())

            #if clipped:
            #    clipStatus = 'clipped'
            #else:
            #    clipStatus = '       '

            file1.write( ' '.join([ tr.stats.station, str(data2max), str(whereDataMax), "\n"]) ) 
            #file1.write( '  '.join( [tr.stats.station, '  '] ) ) 
            #file1.write( "{:7.2f}   {:2d}".format(data2max, whereDataMax) )
            #file1.write( '  '.join([ '  ', clipStatus, "\n"]) ) 

    file1.close()

    return fig

def calcSgram( data, fs, yaxis='linear' ):

    #data = tr.data
    #npts = tr.stats.npts
    #dt = tr.stats.delta
    #fs= tr.stats.sampling_rate
    npts = len( data )
    fnyq = fs/2.0
    nfft = util.next_pow_2(npts) * 4

    if yaxis == 'linear':
        fff, ttt, sss = signal.spectrogram( data, fs=fs, window='hamming', nperseg=256, noverlap=128, nfft=nfft, detrend='constant' )
    elif yaxis == 'log':
        fff, ttt, sss = signal.spectrogram( data, fs=fs, window='hamming', nperseg=256, noverlap=128, nfft=nfft, detrend='constant' )

    return fff, ttt, sss


def streamFiddle( st, what ):
    # Replaces stream with modified values
    st2 = Stream()
    nTrace = len( st )
    for itr in range(nTrace):
        tr = st[itr]
        tr2 = tr
        data = tr.data
        if what == 'abs':
            data = np.abs( data )
        elif what == 'sqrt':
            data = np.sqrt( data )
        elif what == 'env':
            analytic_signal = hilbert(data)
            data = np.abs(analytic_signal)
            data = util.smooth( data, 50 )
        elif what == '':
            data = np.log( data )
        tr2.data = data
        st2.append( tr2 )
    return st2



def streamFiddle3C( st, what ):
    # Replaces 3C streams with modified values
    # DOES NOT WORK IF ANY NON-3C DATA IS PASSED
    st2 = Stream()
    nTrace = len( st )
    for itr in xrange(0,nTrace,3):
        tr1 = st[itr]
        trnew = tr
        tr2 = st(itr+1]
        tr3 = st(itr+2]
        data1 = tr1.data
        data2 = tr2.data
        data3 = tr3.data
        if what == 'vec':
            data = np.sqrt( np.square(data1) + np.square(data2) + np.square(data3) )
        trnew.data = data
        st2.append( trnew )
    return st2



def streamRMS( st, windowMinutes ):

    # Calculates RMS for each trace in a stream

    st2 = Stream()

    nTrace = len( st )

    for itr in range(nTrace):

        tr = st[itr]

        t = tr.stats.starttime
        dataRMS = []
        timeRMS = []

        while t < tr.stats.endtime:
            trace_tmp = tr.slice( t, t + windowMinutes*60 )
            dataValue = np.sqrt( np.mean( np.square( trace_tmp.detrend("demean") ) ) )
            dataRMS.append( dataValue )
            t += windowMinutes * 60
            timeRMS.append( t )

        tr2 = Trace( data=np.array(dataRMS) )
        tr2.stats.delta = windowMinutes*60
        tr2.stats.npts = len(dataRMS)
        tr2.stats.sampling_rate = 1.0/(windowMinutes*60)
        tr2.stats.starttime  = tr.stats.starttime + windowMinutes*60
        tr2.id = tr.id

        st2.append( tr2 )

    return st2


def isClipped( tr ):

    # Tests if a waveform is clipped
    clipped = False

    sta = tr.stats.station
    cha = tr.stats.channel
    data = tr.data

    if sta == 'MSS1' and cha == 'SHZ':

        numClipTop = (data >= 31635).sum()
        numClipBot = (data <= -32500).sum()

        if numClipTop > 0 or numClipBot > 0:
            clipped = True

    return clipped




def calcSgramSt( st ):

    tr = st[0]

    # Data
    data = tr.data
    npts = tr.stats.npts
    dt = tr.stats.delta
    samprate = tr.stats.sampling_rate
    fnyq = samprate/2.0
    t = np.arange(0, npts / samprate, 1 / samprate)
    t = t - pre

    # Spectrogram
    fft_zero_pad_fac = 4
    nf =  64
    w0 = 16.0
    cmap='jet'
    plot_args=['k', 'k']
    if fft_zero_pad_fac == 0:
        nfft = npts
    else:
        nfft = util.next_pow_2(npts) * fft_zero_pad_fac
    f_lin = np.linspace(0, 0.5 / dt, nfft // 2 + 1)
    _w = np.zeros((1, nf, npts), dtype=complex)
    _w[0] = rpt.cwt(data, dt, w0, plotFmin, plotFmax, nf)
    ntr = 1
    spec = np.zeros((1, nfft // 2 + 1), dtype=complex)
    spec[0] = np.fft.rfft(data, n=nfft) * dt

    _tfr = np.sqrt( np.abs(_w) )
    spec = np.sqrt( np.abs(spec) )


    return spec 
