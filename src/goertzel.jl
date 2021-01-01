#   Unit "gortzel" of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2020, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements Goertzel's algorithms for estimating
#   DFT complex coefficients at a specific frequency
#   in the time domain for FFT analysis.
#   NB: The result is practically identical to the output of an FFT,
#   but requieres only `t` operations, where `t` is the window length
#   and is not restricted to a discrete Fourier Frequency.
#   Since the complexity of FFT algirthms is t*log(t), Goertzels' algorithm
#   are advantageous when only one or a few coefficients are needed.
#   Smart block on-line implementations are also possible.
#
#   See: Goertzel G (1958) An algorithm for the evaluation of finite trigonometric series
#   The American Mathematical Monthly, 65(1), 34-35.
#   https://pdfs.semanticscholar.org/a5e4/d0faf65627374b1ac82c3c79006d010173c9.pdf
#   See also
#   https://www.st.com/content/ccc/resource/technical/document/design_tip/group0/20/06/95/0b/c3/8d/4a/7b/DM00446805/files/DM00446805.pdf/jcr:content/translations/en.DM00446805.pdf

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #

"""
```julia
function goertzel( x   :: Vector{T},
                   f   :: IntOrReal,
                   sr  :: Int,
                   wl  :: Int,
               startAT :: Int = 1) where T<:Union{Real, Complex}
```

Given a time series as input vector `x` sampled at `sr` sampling rate,
return the DFT complex coefficient at the discrete Fourier frequency which is
the closest to `f`. The coefficient is computed for the time window starting
at `startAt` (one-based) and of lenght `wl`.

If `startAT`=1 (default) the first `wl` samples of the vector are considered.

`wl` does not need to be power of 2.

**See also**: [`goertzel_fast`](@ref), [`goertzel2`](@ref).

**Examples**:
```julia
using FourierAnalysis
sr, t, f, a = 128, 128, 5, 10
v=sinusoidal(a, f, sr, t, 0)
c=goertzel(v, f, sr, t) # should output 0+aim
```
"""
function goertzel( x   :: Vector{T},
                   f   :: IntOrReal,
                   sr  :: Int,
                   wl  :: Int,
               startAT :: Int = 1) where T<:Union{Real, Complex}
    p, a = 2π*round(UInt16, (f*wl)/sr)/wl, 2/wl
    s, c = sincos(p)
    d = 2*c
    return goertzel_fast(x, wl, a, c, s, d, startAT)
end


"""
```julia
function goertzel_fast( x  :: Vector{T},
                        wl :: Int,
                        a  :: Real,
                        c  :: Real,
                        s  :: Real,
                        d  :: Real,
                    startAT:: Int = 1) where T<:Union{Real, Complex}
```

Fast version of the [`goertzel`](@ref) function to be preferred if the function
is invoked repeatedly and speed is of concern. The user provides as arguments:

- a time series as input vector `x`,
- `wl`, the number of samples in `x` (or the window length),
- `a` = ``2/wl``,
- `c` = ``cos(p)`` where ``p``=`2*π*round(UInt16, (f*wl)/sr)/wl`, ``f`` is the desired frequency and ``sr`` the sampling rate,
- `s` = ``sin(p)``
- `d` = 2*`c`
- `startAt` as in the [`goertzel`](@ref) function.

**See also**: [`goertzel`](@ref), [`goertzel2`](@ref).
"""
function goertzel_fast( x  :: Vector{T},
                        wl :: Int,
                        a  :: Real,
                        c  :: Real,
                        s  :: Real,
                        d  :: Real,
                    startAT:: Int = 1) where T<:Union{Real, Complex}
    z1, z2 = T(0), T(0)
    @inbounds for i = startAT : startAT + wl-1
     z0 = muladd(d, z1, x[i]-z2)
     z2 = z1
     z1 = z0
    end
    return muladd(z1, c, -z2)*a + *(z1, s, a)im
end


"""
```julia
function goertzel2( x  :: Vector{T},
                    f  :: IntOrReal,
                    sr :: Int,
                    wl :: Int,
                startAT:: Int = 1) where T<:Union{Real, Complex}
```

Like the [`goertzel`](@ref) function, but allows estimating the DFT coefficient
in the whole positive real line up to the Nyquist frequency.
This is useful when the DFT coefficient is sought
for a frequency which is not a discrete Fourier Frequency.

Using an appropriate taper this function allows almost exact recovery of
both amplitude and phase at all frequencies in case of one sinusoid,
whereas the `goertzel` function can do so only at exact Fourier discrete
frequencies.

**See also**: [`goertzel_fast`](@ref), [`goertzel`](@ref).

**Examples**:
```julia
using FourierAnalysis
sr, t, f, a = 128, 128, 5, 10
ftrue=f+0.15
v=sinusoidal(a, ftrue, sr, t, 0 )
c=goertzel(v, f, sr, t) # should be 0+aim, but clearly fails
d=goertzel2(v, ftrue, sr, t) # this get closer
```
"""
function goertzel2( x  :: Vector{T},
                    f  :: IntOrReal,
                    sr :: Int,
                    wl :: Int,
                startAT:: Int = 1) where T<:Union{Real, Complex}
    q = (2π*f*wl)/sr
    p = q/wl
    a = 2/wl
    s, c = sincos(p)
    d = 2*c
    z = goertzel_fast(x, wl, a, c, s, d, startAT)
    s2, c2 = sincos(q)
    return muladd(real(z), c2, imag(z)*s2) + muladd(-real(z), s2, imag(z)*c2)im
end
