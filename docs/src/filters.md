# filters.jl

This unit implements banks of band-pass filters interfacing the [DSP.jl](https://github.com/JuliaDSP/DSP.jl) package. Those banks
are used to create time-frequency analytic signal representations
via the Hilbert transform.

At the bottom of this page you can find [notes on DSP package useful functions](@ref).

### FilterDesign
```julia
FilterDesign = Union{ZeroPoleGain, FIRWindow}
```

Those are possible filter designs implemented in the DSP package.
For design methods see the  [DSP documentation](https://juliadsp.github.io/DSP.jl/stable/filters/).

```@docs
filterbank
```

### notes on DSP package useful functions
```julia
using DSP, Plots

### create a time-series with dominant frequency at 10Hz
f, sr, t = 10, 128, 512
v=sinusoidal(1., f, sr, t, 0)+randn(512)
plot(v)

### IIR filters
filterBP = digitalfilter(Bandpass(8/(sr/2), 12/(sr/2)), Chebyshev2(8, 10))
filterLP = digitalfilter(Lowpass(13/(sr/2)), Butterworth(4))
filterHP = digitalfilter(Highpass(7/(sr/2)), Butterworth(4))
# forward filter, unlinear phase response
z=filt(filterBP, v)
plot([v, z])
# forward-backward filter, linear phase response
z=filtfilt(filterBP, v)
plot([v, z])

### FIR filters
responsetype = Bandpass(8, 12; fs=sr)
filtkind = FIRWindow(hanning(64))
filterBP = digitalfilter(responsetype, filtkind)
# forward filter, unlinear phase response
z=filt(digitalfilter(responsetype, filtkind), v)
plot([v, z])
# forward-backward filter, linear phase response
z=filtfilt(digitalfilter(responsetype, filtkind), v)
plot([v, z])

### resampling
z1=resample(x, 0.5)   # downsample by 2
z2=resample(x, 2)     # upsample by 2
plot([x, z1, z2])

### estimate peak frequency using ESPRIT algorithm
peakv=esprit(v, 200, 2, sr)[1]
# this may fails.
# If the approximate location of the peak is sought,
first band-pass filter the data
filter = digitalfilter(Bandpass(7/(sr/2), 13/(sr/2)), Chebyshev2(8, 10))
z = filtfilt(filter, v)
peakz=esprit(z, 200, 2, sr)[1]

### mean power frequency
meanfreq(v, t)

### root mean square
rms(v)

### cross-correlation
y=sinusoidal(1., f, sr, t, 0)+randn(512)
plot(xcorr(v, y))

### convolution (using the FFT algorithm)
plot(conv(v, y))

### closest product of 2, 3, 5, and 7 greater than or equal to n.
# This is used to find a window length for which FFTW will be able
# to find an efficient FFT plan
nextfastfft(1031)

### estimate the delay by locating the peak of the cross-correlation.
# The output delay will be positive when v is delayed with respect y,
# negative if advanced, 0 otherwise.
finddelay(v, y)
finddelay([0, 0, 1, 2, 3], [1, 2, 3]) # return 2
finddelay([1, 2, 3], [0, 0, 1, 2, 3]) # return -2

### shift elements of signal v in time by a given amount s of samples
### and fill the spaces with zeros.
# For circular shifting, use circshift.
shiftsignal([1, 2, 3], 2)  # return [0, 0, 1]
shiftsignal([1, 2, 3], -2) # return [3, 0, 0]

### use finddelay() and shiftsignal() to time align v to y.
# Also return the delay of v with respect to y.
alignsignals([0, 0, 1, 2, 3], [1, 2, 3]) # return ([1, 2, 3, 0, 0], 2)
alignsignals([1, 2, 3], [0, 0, 1, 2, 3]) # return ([0, 0, 1], -2)
```
