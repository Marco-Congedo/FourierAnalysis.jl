#   Unit "filters" of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2023,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements filtering procedures interfacing the DSP.jl package


FilterDesign = Union{ZeroPoleGain, FIRWindow}

#######################################################################
"""
```julia
function filterbank(x         :: Vector{T},
                    sr        :: Int,
                    bandwidth :: IntOrReal    = 2;
                filtkind      :: FilterDesign = Butterworth(2),
                fmin          :: IntOrReal    = bandwidth,
                fmax          :: IntOrReal    = sr÷2,
                ⏩           :: Bool         = true) where T<:Real
```

Pass signal vector `x` throught a bank of band-pass filters,
given sampling rate `sr` and `bandwidth` parameter.

The kind of filters is specified by optional keyword argument `filtkind`,
of type [FilterDesign](@ref), using the DSP package.
By default `filtkind` is a forward-backward (linear phase response)
Butterworth IIR filter of order 2 in each direction
(hence, 4th order total). See [notes on DSP package useful functions](@ref)
for tips on how to design other IIR and FIR filters.

Return 2-tuple `(f, Y)`, where `f` is a vector holding the center
frequencies of the filter bank band-pass regions
and `Y` a matrix holding in the ``i^{th}`` column the signal `x`
band-pass iltered by the ``i^{th}`` band-pass filter.
Hence, ```size(Y, 1)=length(x)``` and ```size(Y, 2)=length(f)```.

The filter bank is designed by means of argument `bandwidth`
and optional keyword arguments `fmin` and `fmax`.
These three arguments can be passed either as integers or real numbers.
All band-pass regions have bandwidth equal to the `bandwidth` argument
and overlap with adjacent band-pass regions.
By default the lower limit of the first band-pass region is set at
`bandwidth` Hz, and successive band-pass regions are defined up to,
but excluding, the Nyquist frequency (``sr/2``).
If `fmin` is specified (in Hz), the center frequency of the first band-pass
region is set as close as possible, but not below, `fmin`.
If `fmax` is specified (in Hz), the center frequency of the last band-pass
region is set as close as possible, but not above, `fmax`.

Here are some examples of filter bank definitions given different arguments
(`sr`=16 in all examples).

| bandwidth | fmin | fmax | center frequencies |       band-pass regions       |
|:---------:|:----:|:----:|:------------------:|:-----------------------------:|
|     4     |   -  |  -   |         4          |              2-6              |
|     2     |   -  |  -   | 2, 3, 4, 5, 6      |    1-3, 2-4, 3-5, 4-6, 5-7    |
|     2     |   3  |  7   | 3, 4, 5, 6         |    2-4, 3-5, 4-6, 5-7         |
|     1     |   3  |  5   | 3, 3.5, 4, 4.5, 5  | 2.5-3.5, 3-4, 3.5-4.5, 4-5, 4.5-5.5 |
|    1.1    |   3  |  5   | 2.75, 3.3, 3.85, 4.4, 4.95 | 2.2-3.3, 2.75-8.85, 3.3-4.4, 3.85-4.95, 4.4-5.5 |
|    1.9    |   3  |  5   | 2.85, 3.8, 4.75    | 1.9-3.8, 2.85-4.75, 3.8-5.7 |

!!! note "Nota Bene"
    At least `sr` samples should be trimmed at the beginning and end of the
    output signal `Y`, as those samples are severely distorted by the filtering
    process.

    If keyword optional argument ⏩ is true (default), the filtering is
    multi-threaded across band-pass filters. See [Threads](@ref).

This function is called by the following functions operating
on time-frequency reprsentations: [`TFanalyticsignal`](@ref),
[`TFamplitude`](@ref), [`TFphase`](@ref), [`meanAmplitude`](@ref),
[`concentration`](@ref), [`meanDirection`](@ref), [`comodulation`](@ref),
[`coherence`](@ref).


**Examples**:
```julia
using FourierAnalysis, DSP, Plots
# generate a sinusoidal + noise
f, sr, t = 8, 128, 512
v=sinusoidal(1., f, sr, t, 0)
x=v+randn(t)
flabels, Y=filterbank(x, 128)
flabels, Y=filterbank(x, 128; fmin=4, fmax=32)
flabels, Y=filterbank(x, 128, 4; fmin=4, fmax=32)
flabels, Y=filterbank(x, 128, 4;
                      filtkind=Chebyshev2(8, 10),
                      fmin=4,
                      fmax=16)
# trick for plotting the signal filtered in the band-pass regions
for i=1:size(Y, 2) Y[:, i].+=convert(eltype(Y), i)*1.5 end
mylabels=Array{String}(undef, 1, length(flabels))
for i=1:length(flabels) mylabels[1, i]=string(flabels[i])*" Hz" end
plot(Y; c=:grey, labels=mylabels)
```
"""
function filterbank(x         :: Vector{T},
                    sr        :: Int,
                    bandwidth :: IntOrReal    = 2;
                filtkind      :: FilterDesign = Butterworth(2),
                fmin          :: IntOrReal    = bandwidth,
                fmax          :: IntOrReal    = sr÷2,
                ⏩           :: Bool         = true) where T<:Real

    df, dm, ff, Bp=digitalfilter, filtkind, filtfilt, Bandpass # aliases for DSP.jl stuff
    fb=_ffilterbank(sr, bandwidth; fmin=fmin, fmax=fmax) # get limits for band-passes
    nf=length(fb)-2 # eliminate last frequency, which may be the Nyquist frequency
    flt=[df(Bp(fb[f], fb[f+2]; fs=sr), dm) for f=1:nf] # f band-pass filters
    x_=[zeros(T, sr); x; zeros(T, sr)]     # Before filtering, `x` is padded with `sr` zeros at both sides in order to symmetrize the distortion of DSP.jl butterworth filters.
    Y=Matrix{T}(undef, length(x_), nf)
    _thread(⏩, nf) ? (@threads for f=1:nf Y[:, f]=ff(flt[f], x_) end) :
                               (for f=1:nf Y[:, f]=ff(flt[f], x_) end)
    return fb[2:end-1], Y[sr+1:end-sr, :] # trim the padded segments in Y
end


##################################################################
# internal function: Create a vector holding necessary information
# to build a filter bank, given sampling rate `rt`, the `bandwidth`
# and optional keyword parameters `fmin` and `fmax`.
# `bandwidth`, `fmax` and `fmin` can be either Integers or Floats.
# Let y be the output vector, a filter bank with length(y)-2 filters
# can be constructed with band-pass regions (in Hz) going from y[i-1] to y[i+1],
# where y[i+1]-y[i-1] = `bandwidth` for all i in [2:end-1].
# All the length(y)-2 filters are centered at y[i].
# By default the filter bank is defined up to, but exclude, the Nyquist frequency.
# If `fmin` is specified (in Hz), the lowest center frequency is >= `fmin`.
# If `fmax` is specified (in Hz), the highest center frequency is <= `fmax`.
#
# Examples:
#
# sr=16
# bw=4
# _ffilterbank(sr, bw) # [2.0, 4.0, 6.0]
#
# bw=2
# _ffilterbank(sr, bw) # [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
# _ffilterbank(sr, bw; fmin=3, fmax=7) # [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
#
# bw=1
# _ffilterbank(sr, bw; fmin=3, fmax=5) # [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
#
# bw=1.1
# _ffilterbank(sr, bw; fmin=3, fmax=5) # [2.2, 2.75, 3.3, 3.85, 4.4, 4.95, 5.5]
#
# bw=1.9
# _ffilterbank(sr, bw; fmin=3, fmax=5) # [1.9, 2.85, 3.8, 4.75, 5.7]
function _ffilterbank(sr        :: Int,
                      bandwidth :: IntOrReal;
                  fmin          :: IntOrReal = bandwidth,
                  fmax          :: IntOrReal = sr÷2)
    f=collect(range(0, step=bandwidth/2, stop=sr÷2))
    popfirst!(f)
    a=max(1, (findmin([abs(g-fmin) for g∈f])[2])-1)
    b=min(length(f), (findmin([abs(g-fmax)  for g∈f])[2])+1)
    f[b]>=sr÷2 ? (return f[a:b-1]) : (return f[a:b])
end
