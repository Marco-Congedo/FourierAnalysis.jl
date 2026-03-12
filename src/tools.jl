#   Unit "tools" of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2025,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements useful functions for Fourier analysis.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #

"""
```julia
function sinusoidal(a :: IntOrReal,
                    f :: IntOrReal,
                   sr :: Int,
                    t :: Int,
                    θ :: IntOrReal = 0.;
                DClevel = 0.)
```

Generate a sinusoidal wave with peak amplitude `a`, frequency `f`,
sampling rate `sr`, duration (in samples) `t`, angle `θ`
(θ=0 makes a sine, θ=π/2 makes a cosine)
and optional keyword argument `DC` (float), the DC level defaulting to zero.
It is adopted the convention that a sine wave starts at zero.

**Examples**:
```julia
using FourierAnalysis, Plots

# create and plot a sinusoidal wave of 128 samples with
# peak amplitude 1, frequency 12Hz, sr=64, phase=π/2
v=sinusoidal(1., 12, 64, 128, π/2)
plot(v)

# estimate amplitude of a sinusoidal wave using Goertzel algorithm
f, sr, t = 32, 128, 128
v=sinusoidal(3., f, sr, t, 0)
c=goertzel(v, f, sr, t) # c should be equal to 0+3.0im
```
"""
function sinusoidal(a  :: IntOrReal,
                    f  :: IntOrReal,
                    sr :: Int,
                    t  :: Int,
                    θ  :: IntOrReal = 0.;
                DC = 0.)
    Δi = inv(sr)
    f2π = f * 2π
    ν=Vector([DC+(a*sin(f2π * i * Δi + θ)) for i=0:t-1])
end


"""
```julia
function fres(sr :: Int,
              wl :: Int)
```

FFT **f**requency **res**olution given sampling rate `sr` and window length `wl`.

**See also**: [`f2b`](@ref), [`b2f`](@ref), [`fdf`](@ref), [`brange`](@ref).

**Examples**:
```julia
using FourierAnalysis
fres(1024, 2048) # return 0.5
```
"""
fres(sr :: Int,
     wl :: Int) = sr/wl


"""
```julia
function f2b(f :: IntOrReal,
            sr :: Int,
            wl :: Int;
        DC :: Bool = false)
```

**f**requency **to b**in. Return the bin (position) in a real-FFT vector best matching
a frequency `f` (in Hz), given sampling rate `sr` and window length `wl`.
The frequency can be given either as an integer or as a real number.

If the requested `f` is exactly in between two Fourier discrete frequencies,
then the smallest of the two equidistant frequencies is returned.

The FFT vector is assumed to be 1-based (as always in Julia).
If `DC` is false, the first discrete frequency is assumed to be at bin 1,
otherwise the DC is assumed to be at bin 1 and the first discrete frequency
at bin 2.

If `DC` is false return 0 for frequencies inferior to half the frequency
resolution.

**See also**: [`fres`](@ref), [`b2f`](@ref), [`fdf`](@ref), [`brange`](@ref).

**Examples**:
```julia
using FourierAnalysis
f2b(10, 512, 1024) # return 20
```
"""
f2b(f  :: IntOrReal,
    sr :: Int,
    wl :: Int;
    DC :: Bool = false) =
    !DC*fres(sr, wl)<=f<=sr/2 ? (return round(Int, (f*(wl/sr))-eps()) + DC) :
    @error 📌*", call to f2b function; invalid frequency argument. The frequency must be comprised between $(!DC*fres(sr, wl)) and $(sr/2) (half the sampling rate)" f
#¤ if DC is true the first frequency is the DC level, hence !DC*fres(sr, wl)
#¤ is equal to zero since !DC=0, otherwise the first frequency is fres(sr, wl)


"""
```julia
function b2f(bin :: Int,
              sr :: Int,
              wl :: Int;
        DC :: Bool = false)
```

**b**in **to f**requency. Return the closest discrete Fourier frequency (in Hz)
that corresponds to a bin (position) in a real-FFT vector,
given sampling rate `sr` and window length `wl`.
The FFT vector is assumed to be 1-based, as always in Julia.

If `DC` is false, the first discrete frequency is assumed to be at bin 1,
otherwise the DC level is assumed to be at bin 1 and the first discrete
frequency at bin 2.

**See also**: [`f2b`](@ref), [`fres`](@ref), [`fdf`](@ref), [`brange`](@ref).

**Examples**:
```julia
using FourierAnalysis
f2b(20, 512, 1024) # return 40
f2b(10, 128, 128) # return 10
```
"""
b2f(bin :: Int,
    sr  :: Int,
    wl  :: Int;
    DC :: Bool = false) =
    0<bin<=wl÷2+DC ? (return (bin-DC)*(sr/wl) ) :
    @error 📌*", call to b2f function; invalid bin argument. The bin must be comprised between 1 and half the window length (+ 1 if DC=true)" bin
#¤ if DC is true there are wl÷2+1 bins, otherwise there are wl÷2 bins


"""
```julia
function fdf(sr :: Int,
             wl :: Int;
          DC :: Bool = false)
```
Return a vector with all **F**ourier **d**iscrete **f**requencies for a real-FFT,
given sampling rate `sr` and window length `wl`.
If `DC` is false, the first discrete frequency starts at bin (position) 1 and
the length of the vector is ``wl÷2`` (integer division), otherwise the DC level is at position 1.
and the length of the vector is ``(wl÷2)+1``.

**See also**: [`f2b`](@ref), [`fres`](@ref), [`b2f`](@ref), [`brange`](@ref).

**Examples**:
```julia
using FourierAnalysis
fdf(8, 16)
# return the 8-element Array{Float64,1}:
# [0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4.0]
```
"""
fdf(sr :: Int,
    wl :: Int;
    DC :: Bool = false) = [i*(sr/wl) for i=!DC:wl÷2]
#¤ if DC is true !DC=0, hence !DC:wl÷2 starts at zero


"""
```julia
function brange(wl :: Int;
             DC :: Bool = false)
```
Return a range of bins for a real-FFT vector covering all Fourier discrete
frequencies given window length `wl`.

If `DC` is false, the range is ``1:(wl÷2)`` (integer division),
otherwise it is ``1:(wl÷2)+1``.

**See also**: [`f2b`](@ref), [`fres`](@ref), [`b2f`](@ref), [`fdf`](@ref).

**Examples**:
```julia
using FourierAnalysis
brange(0.5, 8) # return 1:4
```
"""
brange(wl::Int;
        DC::Bool=false) = (1:wl÷2+DC)


"""
```julia
function bbands(sr :: Int,
                wl :: Int,
         bandwidth :: IntOrReal;
    DC :: Bool = false)
```

Return a vector of integers holding the limits of all `bandwidth`-spaced
band-pass regions of a real-FFT, in bins of discrete Fourier frequencies,
from one to ``wl÷2`` (integer division).

This is used by function [`bands`](@ref).

To know the frequencies in Hz to which these bins correspond, call
[`fbands`](@ref).

**See**: [`bands`](@ref).

**See also**: [`fbands`](@ref).

**Examples**:
```julia
using FourierAnalysis
bbands(128, 256, 16) # return [1, 32, 64, 96, 128]
fbands(128, 256, 16) # return [0.5, 16.0, 32.0, 48.0, 64.0]

bbands(128, 256, 16; DC=true) # return [2, 33, 65, 97, 129]
fbands(128, 256, 16; DC=true) # return [0.5, 16.0, 32.0, 48.0, 64.0]

bbands(128, 128, 16) # return [1, 16, 32, 48, 64]
fbands(128, 128, 16) # return [1.0, 16.0, 32.0, 48.0, 64.0]
```
"""
function bbands(sr        :: Int,
                wl        :: Int,
                bandwidth :: IntOrReal;
            DC        :: Bool=false)
    fr=fres(sr, wl)
    if bandwidth<fr
        @error 📌*", call to function `bbands` or `bandsaverage`: bandwidth cannot be smaller than the FFT frequency resolution" bandwidth fr
        return
    end
    b=collect(StepRange(0, f2b(bandwidth, sr, wl; DC=DC), wl÷2))
    b[1]=1
    if b[2]==1 popfirst!(b) end
    return DC ? b.+=1 : b
end

"""
```julia
function fbands(sr :: Int,
                wl :: Int,
         bandwidth :: IntOrReal;
      DC :: Bool = false)
```
Return a vector of Frequencies (in Hz) to which the bins created by a call
to function [`bbands`](@ref) with the same arguments correspond.

**See**: [`bbands`](@ref).

**See also**: [`bands`](@ref).
"""
fbands(sr        :: Int,
       wl        :: Int,
       bandwidth :: IntOrReal;
    DC :: Bool = false) = b2f.(bbands(sr, wl, bandwidth; DC=DC), sr, wl; DC=DC)


"""
```julia
(1)
function decibel(S :: Union{Real, AbstractArray{T}}) where T<:Real

(2)
function decibel(S1 :: Union{Real, AbstractArray{T}},
            S2 :: Union{Real, AbstractArray{T}}) where T<:Real
```

Convert (1) a measure `S`, or (2) a ratio between two measures `S1`./`S2`
into deciBels.

Input measures can be real numbers or real arrays of any dimensions.

For array input, the ratio and the conversion is computed element-wise.

**Examples**:
```julia
using FourierAnalysis
v=sinusoidal(3., 1, 128, 256, 0)
s=spectra(v, 128, 256; func=decibel) # compute the spectra in dB
s.y # show the spectra

decibel(s.y)

decibel(10.0)

N=abs.(randn(3, 3))
decibel(N)
```
"""
decibel(S::Union{Real, AbstractArray{T}}) where T<:Real = 10*log10.(S)
decibel(S1::Union{Real, AbstractArray{T}},
   S2::Union{Real, AbstractArray{T}}) where T<:Real = 10*log10.(S1./S2)


"""
```julia
(1)
function amplitude(c::Complex;
                        func::Function=identity) = func(abs(c))

(2)
function amplitude(A::AbstractArray{T};
                        func::Function=identity) where T<:Complex

(3)
function amplitude(A::TFAnalyticSignal;
                        func::Function=identity)

(4)
function amplitude(𝐀::TFAnalyticSignalVector;
                        func::Function=identity)
```
(1)

Return the amplitude (modulus) of a complex number.
This corresponds to Julia's
[abs](https://docs.julialang.org/en/v1/base/math/#Base.abs) function.
It is here provided for syntactic consistency with the following methods.

(2)

Return the amplitude of a complex array `Z`.
Typically, `Z` holds analytic signal, in which case the output is
the analytic (instantaneous) amplitude (also known as envelope).
The output is a real array of the same size as `Z`.

(3)

Return a real matrix with the analytic (instantaneous) amplitude of the
[TFAnalyticSignal](@ref) object `Z`. The output is of the same
size as the data field `Z.y`.

(4)

As (3), but return a vector of amplitude matrices for all
[TFAnalyticSignal](@ref) objects in 𝐀

~

In all methods if a function is provided by the optional keyword
argument `func`, it is applied element-wise to the output. For example,
- passing `func=x->x^2` will return the power,
- passing `func=x->log(x^2)` will return the log-power,
- passing `func=x->decibel(x^2)` will return the power in deciBels.

**See**: [TFAnalyticSignal](@ref).

**Examples**:
```julia
using FourierAnalysis, Plots
x=sinusoidal(10, 2, 128, t*4, 0).*sinusoidal(10, 1, 128, t*4, 0)

# amplitude and phase of a vector using analytic signal standard method
y=analyticsignal(x)
a=amplitude(y)
ϕ=phase(y, func=x->(x+π)/2π*50)
plot([x, a, ϕ]; labels=["signal", "amplitude", "phase"])

# see what happen if `x` contains energy in frequencies below sr/wl Hz
# (see documentation of `analyticSignal` function)
y=analyticsignal(x, 64)
a=amplitude(y)
ϕ=phase(y, func=x->(x+π)/2π*50)
plot([x, a, ϕ]; labels=["signal", "amplitude", "phase"])

# unwrapped phase
# the line below will do nothing as argument `unwrapdims` is 0 by default
ϕ2=unwrapPhase(phase(y))
# this will do the job
ϕ2=unwrapPhase(phase(y); unwrapdims=1)
plot([x, a, ϕ2./25]; labels=["signal", "amplitude", "unwr. phase"])

# amplitude from analytic signal of a data matrix holding multiple series
X=randn(t, 4)
Y=analyticsignal(X)
A=amplitude(Y)
plot(A[:, 1:2])

# phase
𝛷=phase(Y)
plot(𝛷[:, 1:1])

# unwrapped phase
𝛷2=unwrapPhase(𝛷; unwrapdims=1)
plot(𝛷2)

# phase represented in [-1, 1]
𝛷=phase(Y, func=x->(x+π)/2π)
plot(𝛷[:, 1:1])

# sine of the phase
𝛷=phase(Y, func=sin)
plot(𝛷[:, 1:1])

# get Amplitude and phase from analytic Signal
A, 𝛷=polar(Y)
A
𝛷
```
"""
amplitude(c::Complex;
            func::Function=identity) = func(abs(c))

amplitude(A::AbstractArray{T};
            func::Function=identity) where T<:Complex = @.func(abs(A))

amplitude(A::TFAnalyticSignal;
            func::Function=identity) = amplitude(A.y, func)

amplitude(𝐀::TFAnalyticSignalVector;
            func::Function=identity) = [amplitude(A.y, func) for A ∈ 𝐀]

# Standard atan function in other languages
atan2(z::Complex) = atan(imag(z), real(z))

"""
```julia
(1)
function phase(z::Complex; func::Function=identity)

(2)
function phase(Z::AbstractArray{T};
                        unwrapdims::Int=0,
                        func::Function=identity) where T<:Complex

(3)
function phase(Z::TFAnalyticSignal;
                        unwrapped::Bool=false,
                        func::Function=identity)

(4)
function phase(𝐙::TFAnalyticSignalVector;
                        unwrapped::Bool=false,
                        func::Function=identity)

```
(1)

Return the phase (argument) of a complex number.
This corresponds to a standard
[atan2](https://en.wikipedia.org/wiki/Atan2) function.
It is here provided for syntactic consistency with the following methods.

(2)

Return the phase of a complex array `Z`.
Typically, `Z` holds analytic signal, in which case the output is
the analytic (instantaneous) phase.
The output is a real array of the same size as `Z`.

If optional keyword argument `unwrapdims` is > 0, return the unwrapped phase
along the `unwrapdims` dimension of the array. For example, if `Z` is a matrix,
passing `unwrapdims=1` unwrap the phase indipendently along its columns.

(3)

Return a real matrix with the analytic (instantaneous) phase of the
[TFAnalyticSignal](@ref) object `Z`. The output is of the same
size as the data field `Z.y`.

If optional keyword argument `unwrapped` is true, return the unwrapped phase
along the time dimension of the analytic signal (dims=2).

(4)

As (3), but return a vector of phase matrices for all
[TFAnalyticSignal](@ref) objects in 𝚯.

~

In all methods by default the phase is returned in [−π, π].
If a function is provided by the optional keyword argument `func`,
it is applied to the phase. For example
- passing `func=x->x+π` will return the phase in [0, 2π],
- passing `func=x->x/π` will return the phase in [-1, 1],
- passing `func=sin` will return the sine of the phase.

!!! note "Nota Bene"
    If in method (2) `unwrapdims` is >0 or in method (3) and (4)
    `unwrapped` is true, the function `func` is applied to the unwrapped phase.

**See**: [`unwrapPhase`](@ref), [TFAnalyticSignal](@ref).

**Examples**: see examples of [`amplitude`](@ref).
"""
phase(z::Complex; func::Function=identity) = func(atan2(z))

phase(Z::AbstractArray{T};
        unwrapdims::Int=0,
        func::Function=identity) where T<:Complex =
    unwrapdims>0 ? func.(unwrapPhase(Z, unwrapdims=unwrapdims)) : @.func(atan2(Z))  # see Base.Broadcast.@__dot__

phase(Z::TFAnalyticSignal;
        unwrapped::Bool=false, func::Function=identity) =
    phase(Z.y; unwrapdims=2*unwrapped, func=func) # unwrap==true -> unwrapdims=2, else unwrapdims=0

phase(𝐙::TFAnalyticSignalVector; unwrapped::Bool=false, func::Function=identity) =
    [phase(Z; unwrapped=unwrapped, func=func) for Z ∈ 𝐙]

"""
```julia
(1)
function polar(c::Complex)

(2)
function polar(Z::AbstractArray{T}) where T<:Complex

(3)
function polar(Z::TFAnalyticSignal)

```

(1)

Return the amplitude (modulus) and phase (argument)
of a complex number as a 2-tuple.

(2)

Return the amplitude and phase of a complex array `Z`.
Typically, `Z` holds analytic signal, in which case return
the analytic (instantaneous) amplitude and phase.
The output is a tuple of two real arrays of the same size as data field `Z.y`.

(3)

Return the analytic (instantaneous) amplitude
and phase of the [TFAnalyticSignal](@ref) object `Z`.
The output is a tuple of two real arrays of the same size as data field `Z.y`.

~

In all methods the phase is returned in [−π, π].

**See**: [`amplitude`](@ref), [`phase`](@ref), [TFAnalyticSignal](@ref).

**Examples**: see examples of [`amplitude`](@ref).
"""
polar(c::Complex) = amp(c), atan2(c)

polar(Z::AbstractArray{T}) where T<:Complex = amplitude(Z), phase(Z)

polar(Z::TFAnalyticSignal) = polar(Z.y)


# instantaneous frequency from wikipedia
  #function frequency(Y::TimeFrequency)
    #  W=zeros(real(eltype(Y.Z)), size(Y.Z))
    #  @inbounds for j=2:size(Y.Z, 2), i=1:size(Y.Z, 1)
    #                W[i, j]=atan2(Y.Z[i, j]*conj(Y.Z[i, j-1]))/2π end
    #  return
 # end



"""
```julia
(1)
function unwrapPhase(Z::AbstractArray{T};
                                unwrapdims::Int=0) where T<:Complex

(2)
function unwrapPhase(ϴ::AbstractArray{T};
                                unwrapdims::Int=0) where T<:Real

(3)
unwrapPhase(ϴ::TFPhase) [constructor of a TFPhase object]

(4)
unwrapPhase(𝚯::TFPhaseVector) [constructor of a TFPhaseVector object]
```
(1)

If optional keyword argument `unwrapdims` is > 0, compute the phase
(argument) from a *complex* array and unwrap it along the `unwrapdims`
dimension, otherwise (default) return `Z`.
Typically, `Z` holds analytic signal.

(2)

If optional keyword argument `unwrapdims` is > 0, unwrap along the
`unwrapdims` dimension a *real* array holding phase data in [−π, π],
otherwise return `ϴ`.

(3)

Construct a [TFPhase](@ref) object by unwrapping its phase along the time
dimension and copying all other fields from the `ϴ` object. If `ϴ.func`
is different from the `identity` (do nothing) function,
return instead an error message.

(4)

As (3), but conctruct a [TFPhaseVector](@ref) holding
[TFPhase](@ref) objects in 𝚯 with the phase unwrapped.
`ϴ.func` must be the identity function for all ϴ ∈ 𝚯.

The unwrapped phase is defined as the cumulative sum
of the phase (along the relevant dimension)
once this is represented in [0, 2π].

**Examples**: see examples of [`amplitude`](@ref).
"""
unwrapPhase(Z::AbstractArray{T}; unwrapdims::Int=0) where T<:Complex =
    unwrapdims>0 ? accumulate(+, (atan2.(Z)).+π; dims=unwrapdims) : Z

unwrapPhase(ϴ::AbstractArray{T}; unwrapdims::Int=0) where T<:Real =
    unwrapdims>0 ? accumulate(+, ϴ.+π; dims=unwrapdims) : ϴ

unwrapPhase(ϴ::TFPhase) =
    ϴ.func==identity ? TFPhase(unwrapPhase(ϴ.y; unwrapdims=2), ϴ.bandwidth,
                               ϴ.flabels, ϴ.nonlinear,
                               ϴ.fsmoothing, ϴ.tsmoothing, true, ϴ.func) :
   @error 📌*", call to unwrapPhase constructor; I shall not unwrap a phase on which a function has been applied. The unwrapped phase object has not been created.)" ϴ.func


unwrapPhase(𝚯::TFPhaseVector) =
   sum(ϴ.func==identity for ϴ ∈ 𝚯)==length(𝚯) ? TFPhaseVector([
                    TFPhase(unwrapPhase(ϴ.y; unwrapdims=2), ϴ.bandwidth,
                    ϴ.flabels, ϴ.nonlinear, ϴ.fsmoothing, ϴ.tsmoothing,
                    true, ϴ.func) for ϴ ∈ 𝚯]) :
  @error 📌*", call to unwrapPhase constructor; I shall not unwrap a phase on which a function has been applied. The unwrapped phase object has not been created.)" ϴ.func



"""
```julia
function sameParams(𝐒        :: FDobjectsVector,
                    funcname :: String)
```


Return true if all objects in 𝐒 have the same `sr`, `wl`, `DC`, `taper`,
`func`(only for [SpectraVector](@ref) objects), `nonlinear` (only for
[CrossSpectraVector](@ref) and [CoherenceVector](@ref)) and `smoothing` fields,
otherwise print an error message pointing to the first field that is not
identical in all objects and return `Nothing`. This method applies to
all [FDobjectsVector](@ref) types, that is,
to [SpectraVector](@ref), [CrossSpectraVector](@ref) and
[CoherenceVector](@ref).

`funcname` is an optional string that the user can
provide. It is inserted into the error message to locate the part
of the code that generated the error. By defalut, "unknown" is used.
"""
sameParams(𝐒        :: FDobjectsVector,
           funcname :: String = "unknown") =
   if      !_allsame([S.sr for S ∈ 𝐒])
               @error 📌*", "*funcname*" function. All objects in argument of type $typeof(𝐒) must have the same sampling rate (.sr field)"
   elseif  !_allsame([S.wl for S ∈ 𝐒])
               @error 📌*", "*funcname*" function. All objects in argument of type $typeof(𝐒) must have the same window length (.wl field)"
   elseif  !_allsame([S.DC for S ∈ 𝐒])
               @error 📌*", "*funcname*" function. All objects in argument of type $typeof(𝐒) must all hold or not hold DC level (.DC field)"
   elseif  !_allsame([S.taper for S ∈ 𝐒])
               @error 📌*", "*funcname*". All objects in argument of type $typeof(𝐒) must have the same taper (.taper field)"
   elseif  𝐒 isa SpectraVector && !_allsame([S.func for S ∈ 𝐒])
               @error 📌*", "*funcname*". All objects in argument of type $typeof(𝐒) must have been subjected to the same function (.func field)"
   elseif  (𝐒 isa CrossSpectraVector || 𝐒 isa CoherenceVector) && !_allsame([S.nonlinear for S ∈ 𝐒])
               @error 📌*", "*funcname*". The objects in argument of type $typeof(𝐒) must be either all linear or all non-linear (.nonlinear field)"
   elseif  !_allsame([S.smoothing for S ∈ 𝐒])
               @error 📌*", "*funcname*". All the objects in argument of type $typeof(𝐒) must have been subjected to the same smoother (.smoothing field)"
   elseif      return true
   end


"""
```julia
function sameParams(𝒀        :: TFobjectsVector,
                    funcname :: String) =
```
Return true if all objects in 𝒀 have the same `bandwidth`, `nonlinear`,
`fsmoothing` and `tsmoothing` field, otherwise print an error message
pointing to the first field that is not identical in all objects and
return `Nothing`. This method applies to all [TFobjectsVector](@ref) types,
that is, [TFAnalyticSignalVector](@ref), [TFAmplitudeVector](@ref) and
[TFPhaseVector](@ref).

`funcname` has the same meaning as in the previous method.
"""
sameParams(𝒀        :: TFobjectsVector,
           funcname :: String = "unknown") =
   if      !_allsame([Y.bandwidth for Y ∈ 𝒀])
               @error 📌*", "*funcname*" function. All the objects in argument of type $typeof(𝒀) must all be definied with the same bandwidth"
   elseif  !isa(𝒀, TFAmplitudeVector) && !_allsame([Y.nonlinear for Y ∈ 𝒀])
               @error 📌*", "*funcname*" function. All the objects in argument of type $typeof(𝒀) must be definied as either all nonlinear or all linear"
   elseif  !_allsame([Y.fsmoothing for Y ∈ 𝒀])
               @error 📌*", "*funcname*" function. All the objects in argument of type $typeof(𝒀) must be definied with the same frequency smoother"
   elseif  !_allsame([Y.tsmoothing for Y ∈ 𝒀])
               @error 📌*", "*funcname*" function. All the objects in argument of type $typeof(𝒀) must be definied with the same time smoother"
   elseif      return true
   end


"""
```julia
function isLinear(𝒀::Union{FDobjectsVector, TFobjectsVector})
```
Return true if all objects in `𝒀` are linear.
By definition, [Spectra](@ref) and [TFAmplitude](@ref) objects are linear.
[CrossSpectra](@ref), [Coherence](@ref), [TFAnalyticSignal](@ref)
and [TFPhase](@ref) objects may be linear or non-linear.

**See**:[FDobjectsVector](@ref), [TFobjectsVector](@ref).
"""
isLinear(𝒀::FDobjectsVector) =
   𝒀 isa SpectraVector ? true : sum(Y.nonlinear for Y ∈ 𝒀)==0

isLinear(𝒀::TFobjectsVector) =
   𝒀 isa TFAmplitudeVector ? true : sum(Y.nonlinear for Y ∈ 𝒀)==0

"""
```julia
function isNonLinear(𝒀::Union{FDobjectsVector, TFobjectsVector})
```
Return true if all objects in `𝒀` are non-linear.
By definition, [Spectra](@ref) and [TFAmplitude](@ref) objects are linear.
[CrossSpectra](@ref), [Coherence](@ref), [TFAnalyticSignal](@ref)
and [TFPhase](@ref) objects may be linear or non-linear.

**See**:[FDobjectsVector](@ref), [TFobjectsVector](@ref).

"""
isNonLinear(𝒀::FDobjectsVector) =
   𝒀 isa SpectraVector ? false : sum(Y.nonlinear for Y ∈ 𝒀)==length(𝒀)

isNonLinear(𝒀::TFobjectsVector) =
   𝒀 isa TFAmplitudeVector ? false : sum(Y.nonlinear for Y ∈ 𝒀)==length(𝒀)

"""
```julia
(1)
function isUnwrapped(ϴ::TFPhase)

(2)
function isUnwrapped(𝚯::TFPhaseVector)
```
(1)
Return true if the TFPhase objects ϴ have the phase unwrapped.

(2)
Return true if all TFPhase objects in 𝚯 have the phase unwrapped.

**See**: [`unwrapPhase`](@ref), [TFPhase](@ref), [TFPhaseVector](@ref).
"""
isUnwrapped(ϴ::TFPhase) =  ϴ.unwrapped

isUnwrapped(𝚯::TFPhaseVector) =  sum(ϴ.unwrapped for ϴ ∈ 𝚯)==0


### Frequency-domain smoothing ###
###############################################################################

# Internal function: computes coefficients for smoothing
function _getSmoothCoeff(smoother::Smoother)
    if      smoother == hannSmoother
            return [0.25, 0.5, 0.25], 0.75, 0. # 0.75=c1+c2
    elseif  smoother == hammingSmoother
            return [0.23, 0.54, 0.23], 0.77, 0. # 0.77=c1+c2
    elseif  smoother == blackmanSmoother
            return [0.04, 0.25, 0.42, 0.25, 0.04], 0.71, 0.96 # 0.71=c1+c2+c3, 0.96=c1+c2+c3+c4
    end
end


# Internal function: apply frequency-domain smoothing to the `X` vector
# of real or complex numbers, or real or complex lower-triangular matrices
function __smooth(X, c, s, t)
    Y=similar(X)
    k=length(X)
    type=typeof(X)

    if      type isa Vector{LowerTriangular} # || type===Array{LowerTriangular{Complex{Float64},S} where S<:AbstractArray{Complex{Float64},2},1}
            T=LowerTriangular
    elseif  type isa Array{Hermitian} # || Array{Hermitian{Complex{Float64},S} where S<:AbstractArray{Complex{Float64},2},1}
            T=Hermitian
    else    T=eltype(X)
    end

    if length(c) == 3 # `smoother` = `hannSmoother` or `hammingSmoother`
        Y[1] = T(X[1]*(c[2]/s) + X[2]*(c[3]/s))
        @inbounds for i=2:k-1 Y[i] = T(X[i-1]*c[1] + X[i]*c[2] + X[i+1]*c[3]) end
        Y[k] = T(X[k-1]*(c[1]/s) + X[k]*(c[2]/s))
    else # `smoother` = `blackmanSmoother`, for which length(c) == 5
        Y[1] = T(X[1]*(c[3]/s) + X[2]*(c[4]/s) + X[3]*(c[5]/s))
        Y[2] = T(X[1]*(c[2]/t) + X[2]*(c[3]/t) + X[3]*(c[4]/t) + X[4]*(c[5]/t))
        @inbounds for i=3:k-2 Y[i] = T(X[i-2]*c[1] + X[i-1]*c[2] + X[i]*c[3] + X[i+1]*c[4] + X[i+2]*c[5]) end
        Y[k-1] = T(X[k-3]*(c[1]/t) + X[k-2]*(c[2]/t) + X[k-1]*(c[3]/t) + X[k]*(c[4]/t))
        Y[k] = T(X[k-2]*(c[1]/s) + X[k-1]*(c[2]/s) + X[k]*(c[3]/s))
    end
    return type(Y)
end


smooth(smoothing::Smoother, v::Vector{R}) where R<:RealOrComplex =
    return smoothing == noSmoother ? v : __smooth(v, _getSmoothCoeff(smoothing)...)

function smooth(smoothing::Smoother, X::Matrix{R};
                dims::Int=1) where R<:RealOrComplex
    if smoothing == noSmoother
        return X
    end

    if dims ∉ (1, 2)
        throw(ArgumentError(📌*", function smooth: the `dims` keyword must be 1 or 2"))
    end

    if size(X, dims) < 3 && smoothing∉(hannSmoother, hammingSmoother)
        throw(ArgumentError(📌*", function smooth: at least three points along dimension $dims must be available in order to apply a hannSmoother or a hammingSmoother"))
    end

    if size(X, dims) < 5 && smoothing==blackmanSmoother
        throw(ArgumentError(📌*", function smooth: at least five points along dimension $dims must be available in order to apply a blackmanSmoother"))
    end

    Y=similar(X)
    coeff = _getSmoothCoeff(smoothing)
    if dims==1
        for i=1:size(X, 2) Y[:, i]=__smooth(X[:, i], coeff...) end
    elseif dims==2
        for i=1:size(X, 1) Y[i, :]=__smooth(X[i, :], coeff...) end
    end
    return Y
end


# Internal function: smooth the output of `spectra` function; spectra of a Matrix (columns)
function _smooth(smoother::Smoother, S::Spectra, γ, c1, c2)
    Y=similar(S.y)
    @inbounds for i=1:size(S.y, 2) Y[:, i]=__smooth(S.y[:, i], γ, c1, c2) end
    return Spectra(Y, S.sr, S.wl, S.DC, S.taper, S.flabels, S.func, smoother)
end

# Internal function: smooth the output of `spectra` function; k-spectra of a Matrix Vector
_smooth(smoother::Smoother, 𝐒::SpectraVector, γ, c1, c2) =
    SpectraVector([_smooth(smoother, S, γ, c1, c2) for S ∈ 𝐒])

# Internal function: smooth the output of `cross-spectra` or `coherence` function
_smooth(smoother::Smoother, 𝙎::Union{CrossSpectra, Coherence}, γ, c1, c2) =
    typeof(𝙎)(__smooth(𝙎.y, γ, c1, c2), 𝙎.sr, 𝙎.wl, 𝙎.DC, 𝙎.taper, 𝙎.flabels, 𝙎.nonlinear, smoother, 𝙎.tril)

# Internal function: smooth the output of `cross-spectra` or `coherence` function; k-cross-spectra of a Matrix Vector
_smooth(smoother::Smoother, 𝓢::Union{CrossSpectraVector, CoherenceVector}, γ, c1, c2) =
    typeof(𝓢)([_smooth(smoother, 𝐒, γ, c1, c2) for 𝐒 ∈ 𝓢])


"""
```julia
(1)
function smooth(smoothing::Smoother,
                v::Vector{R}) where R<:RealOrComplex


(2)
function smooth(smoothing::Smoother,
                X::Matrix{R}; dims::Int=1) where R<:RealOrComplex

(3)
function smooth(smoother :: Smoother,
                S :: Union{FDobjects, FDobjectsVector})

(4)
function smooth(fsmoothing :: Smoother,
                tsmoothing :: Smoother,
                Y :: Union{TFobjects, TFobjectsVector})
```

Apply a smoothing function of type [Smoother](@ref) to
- (1) a vector of real or complex numbers,
- (2) a real of complex matrix along dimension `dims` (default=1),
- (3) a [FDobjects](@ref) or all objects in a [FDobjectsVector](@ref),
- (4) a [TFobjects](@ref) or all objects in a [TFobjectsVector](@ref).

Methods (1) and (2) are provided for low-level computations.
Methods (3) and (4) are constructors;
for all methods the output is always of the same type as the input.

Method (3) smooths across the frequency dimension:
- for [Spectra](@ref) objects this amounts to smoothing the column vectors in their `.y` field,
- for [CrossSpectra](@ref) and [Coherence](@ref) objects this amounts to smoothing adjacent matrices in their .y field.

Method (4) smooths across the frequency dimension, time dimension or both.
This amounts to smoothing across the column vectors (frequency) and/or row vectors
(time) in the `.y` field of the object.
A smoother must be specified for the frequency dimension
(`fsmoothing`) and for the time dimension (`tsmoothing`).
Either one may be `noSmoother`, but if the two are different from `noSmoother`,
then they must be the same. If smoothing is requested in both the frequency and
time dimension, then the data is smoothed first in the time then in the
frequency dimension.
For [TFPhase](@ref) objects, smoothing is allowed only if the phase is unwrapped.

This function allow smoothing frequency domain and time-frequency domain objects
after they have been created, however, smoothing can also be requested upon
creation. For example, see the documentation of [Spectra](@ref).

!!! note "Nota Bene"
    For methods (1), (2) and (3), if `Smoother` is `noSmoother`, then the input
    is returned unchanged. For method (4) this is the case if both `fsmoother`
    and `tsmoother` are `noSmoother`.

    The data input must hold in the concerned dimension at least three elements
    for applying an Hann or Hamming smoother and at least five elements for
    applying the Blackman smoother.

## Maths

Smoothing of a series ``x`` composed of ``k`` elements is carried out at element
``i`` such as

``x_{i}=ax_{i-2}+bx_{i-1}+cx_{i}+bx_{i+1}+ax_{i+2}``.

The coefficients are

| smoothing window  |  a   |  b   |  c   |
|:-----------------:|:----:|:----:|:----:|
| Hann              | 0    | 0.25 | 0.50 |
| Hamming           | 0    | 0.23 | 0.54 |
| Blackman          | 0.04 | 0.25 | 0.42 |

For 3-point smoothers, the first point is smoothed as

``x_{1}=\\frac{c}{b+c}x_{1} + \\frac{b}{b+c}x_{2}``

and the last (the ``k^{th}``) as

``x_{k}=\\frac{c}{b+c}x_{k} + \\frac{b}{b+c}x_{k-1}``.

For 5-point smoothers, the first point is smoothed as

``x_{1}=\\frac{c}{a+b+c}x_{1} + \\frac{b}{a+b+c}x_{2} + \\frac{a}{a+b+c}x_{3}``,

the second as

``x_{2}=\\frac{b}{a+2b+c}x_{1} + \\frac{c}{a+2b+c}x_{2} + \\frac{b}{a+2b+c}x_{3} + \\frac{a}{a+2b+c}x_{4}``,

the second to last as

``x_{k-1}=\\frac{a}{a+2b+c}x_{k-3} + \\frac{b}{a+2b+c}x_{k-2} + \\frac{c}{a+2b+c}x_{k-1} + \\frac{b}{a+2b+c}x_{k}``

and the last as

``x_{k}=\\frac{a}{a+b+c}x_{k-2} + \\frac{b}{a+b+c}x_{k-1} + \\frac{c}{a+b+c}x_{k}``.

**See**: [Smoother](@ref)

**Examples**:
```julia
using FourierAnalysis, Plots
sr, t, f, a = 128, 128, 10, 0.5
# create a sinusoidal superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# compute Amplitude Spectra
Σ=spectra(v, sr, t; func=√)
bar(Σ.y, labels="raw amplitude spectra")

#smooth spectra
Σ2=smooth(blackmanSmoother, Σ)
bar!(Σ2.y, labels="smoothed amplitude spectra")

# smooth cross-spectra (or coherence) matrices
X=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
S=crossSpectra(X, sr, t) # or coherence (X, sr, t)
# smooth the cross-spectra # or coherence
S2=smooth(blackmanSmoother, S)

# smooth time-frequency object
Y = TFanalyticsignal(v, sr, sr*4)
# smooth frequency
Z=smooth(blackmanSmoother, noSmoother, Y)
# plot amplitude of smoothed analytic signal
heatmap(Z, amplitude)

# smooth AS: smooth both frequency and time
E=smooth(blackmanSmoother, blackmanSmoother, Y)
# plot real part of smoothed analytic signal
heatmap(Z, real)
```
"""
smooth(smoothing::Smoother, S::Union{FDobjects, FDobjectsVector}) =
    return smoothing == noSmoother ? S : _smooth(smoothing, S, _getSmoothCoeff(smoothing)...)

# internal function: Smooth the TimeFrequency types; dim1 = time, dim2=frequency
function _smooth(smoother  :: Smoother,
                 Z         :: TFobjects,
                 frequency :: Bool,
                 time      :: Bool,
                 γ, c1, c2)

    cT=eltype(Z.y)
    if     !frequency & !time
        Z_=Z.y
    elseif !frequency & time
        Z_=Matrix{cT}(undef, size(Z.y))
        @inbounds for i=1:size(Z_, 1) Z_[i, :]=__smooth(Z.y[i, :], γ, c1, c2) end
    elseif frequency & !time
        Z_=Matrix{cT}(undef, size(Z.y))
        @inbounds for i=1:size(Z_, 2) Z_[:, i]=__smooth(Z.y[:, i], γ, c1, c2) end
    elseif frequency & time
        Z_=Matrix{cT}(undef, size(Z.y))
        @inbounds for i=1:size(Z_, 1) Z_[i, :]=__smooth(Z.y[i, :], γ, c1, c2) end
        @inbounds for i=1:size(Z_, 2) Z_[:, i]=__smooth(Z_[:, i], γ, c1, c2) end
    end

    type=typeof(Z)
    if     type===TFAnalyticSignal
        return type(Z_, Z.bandwidth, Z.flabels, Z.nonlinear, frequency ? smoother : noSmoother, time ? smoother : noSmoother)
    elseif type===TFAmplitude
        return type(Z_, Z.bandwidth, Z.flabels, frequency ? smoother : noSmoother, time ? smoother : noSmoother, Z.func)
    elseif type===TFPhase
        return type(Z_, Z.bandwidth, Z.flabels, Z.nonlinear, frequency ? smoother : noSmoother, time ? smoother : noSmoother, Z.unwrapped, Z.func)
    end
end

# Smooth the TFobjects types; dim1 = time, dim2=frequency
smooth(fsmoothing :: Smoother,
       tsmoothing :: Smoother,
       Z          :: TFobjects,
       funcname   :: String="") =
    if Z isa TFPhase && !Z.unwrapped
        @error 📌*", function "*funcname*": smoothing of TFPhase objects can be applied only if the phase is unwrapped"
        return
    elseif fsmoothing == tsmoothing == noSmoother
        return Z
    elseif fsmoothing==tsmoothing || (fsmoothing==noSmoother || tsmoothing==noSmoother)
        smoother=max(fsmoothing, tsmoothing) # return the smoother if they are equal or one of them is noSmoother
        return _smooth(smoother, Z, fsmoothing≠noSmoother, tsmoothing≠noSmoother, _getSmoothCoeff(smoother)...)
    else
        @error 📌*", function "*funcname*": if smoothing is requested, then the `Smoother` for frequency and time must be of the same kind"
    end

# Smooth the TFobjectsVector types; dim1 = time, dim2=frequency
smooth(fsmoothing :: Smoother,
       tsmoothing :: Smoother,
       𝐙          :: TFobjectsVector,
       funcname   :: String="") =
   if 𝐙 isa TFPhaseVector && isUnwrapped(𝐙)
       @error 📌*", function "*funcname*": smoothing of TFPhaseVector objects can be applied only if the phase is unwrapped for all TFPhase objects it holds"
       return
   elseif fsmoothing == tsmoothing == noSmoother
       return 𝐙
   elseif fsmoothing==tsmoothing || (fsmoothing==noSmoother || tsmoothing==noSmoother)
       smoother=max(fsmoothing, tsmoothing) # return the smoother if they are equal or one of them is noSmoother
       return typeof(𝐙)([_smooth(smoother, Z, fsmoothing≠noSmoother, tsmoothing≠noSmoother, _getSmoothCoeff(smoother)...) for Z ∈ 𝐙])
   else
       @error 📌*", function "*funcname*": if smoothing is requested, then the `Smoother` for frequency and time must be of the same kind"
   end
#    return smoother == noSmoother ? 𝐙 :
#           typeof(𝐙)([_smooth(smoother, Z, frequency, time, _getSmoothCoeff(smoother)...) for Z ∈ 𝐙])


### Frequency-domain data extraction and averaging ###
###############################################################################

# get an fInterval type and output a UnitRange whose start and stop are the
# corresonding Fourier discrete frequency bins
function _getfrange(S::FDobjects, frange::fInterval, funcname::String)
    if      isa(frange, IntOrReal)
            !S.DC*fres(S.sr, S.wl)<=frange<=S.sr/2 ? (return f2b(frange, S.sr, S.wl; DC=S.DC):f2b(frange, S.sr, S.wl; DC=S.DC)) :
            @error 📌*", "*funcname*" function passed invalid frange Int or Real argument. The frequency must be comprised between $(!S.DC*fres(S.sr, S.wl)) and $(S.sr/2) (half the sampling rate)" frange
            return 0
    elseif  isa(frange, Colon)
            return 1:size(S.y, 1)
    elseif  isa(frange, Tuple{IntOrReal, IntOrReal})
            # println(frange[1], ", ", frange[2])
            if      frange[1]>frange[2]
                        @error 📌*", "*funcname*" function passed invalid frange 2-tuple argument. The second element cannot be inferior to the first" frange[1] frange[2]
                        return 0
            elseif  frange[1]<!S.DC*fres(S.sr, S.wl) || frange[2]>S.sr/2
                        @error 📌*", "*funcname*" function passed invalid frange 2-tuple argument. The two frequencies must be comprised between $(!S.DC*fres(S.sr, S.wl)) and $(S.sr/2) (half the sampling rate)" frange[1] frange[2]
                        return 0
            else        return f2b(frange[1], S.sr, S.wl; DC=S.DC): f2b(frange[2], S.sr, S.wl; DC=S.DC)
            end
    end
end

# as _getfrange above, but for a FDobjectsVector object.
# By default, it is checked that the elements of the vector are homogeneous.
# If so return the UnitRange obtained on the first FDobject in the vector,
# otherwise return 0.
# if `check` is false, the check is skypped.
_getfrange(𝐒::FDobjectsVector, frange::fInterval, funcname::String; check::Bool=true) =
    check ? (sameParams(𝐒, funcname) ? _getfrange(𝐒[1], frange, funcname) : 0) :
    _getfrange(𝐒[1], frange, funcname)

_extract(S::Spectra, frange::UnitRange{Int64}) =
    _isnull(frange) && size(S.y, 2)==1 ? S.y[first(frange), 1] :
                                         copy(S.y[frange, :])

_extract(𝐒::SpectraVector, frange::UnitRange{Int64}) =
    _isnull(frange) && size(𝐒[1].y, 2)==1 ? [S[first(frange), 1] for S ∈ 𝐒] :
                                            [S[frange, :] for S ∈ 𝐒]

function _extract(S::Union{CrossSpectra, Coherence}, frange::UnitRange{Int64})
    mattype = S.tril ? LowerTriangular : Hermitian
    _isnull(frange) ? mattype(copy(S.y[first(frange)])) :
                      typeof(S.y)([mattype(S.y[i]) for i∈frange])
end

function _extract(𝐒::Union{CrossSpectraVector, CoherenceVector}, frange::UnitRange{Int64})
    mattype = 𝐒[1].tril ? LowerTriangular : Hermitian
    _isnull(frange) ? typeof(𝐒[1].y)([mattype(S.y[first(frange)]) for S ∈ 𝐒]) :
                      [typeof(𝐒[1].y)([mattype(S.y[i]) for i∈frange]) for S ∈ 𝐒]
end

_mean(S::Spectra, frange::UnitRange{Int64}) =
    size(S.y, 2)==1 ? mean(view(S.y, frange, 1)) :
                     [mean(view(S.y, frange, k)) for k=1:size(S.y, 2)]

_mean(S::Union{CrossSpectra, Coherence}, frange::UnitRange{Int64}) =
    eltype(S.y)(mean(view(S.y[i], :, :) for i∈frange))

_mean(𝐒::SpectraVector, frange::UnitRange{Int64}) =
    size(𝐒[1].y, 2)==1 ? [mean(view(S.y, frange, 1)) for S ∈ 𝐒] :
                         [[mean(view(S.y, frange, k)) for k=1:size(S.y, 2)] for S ∈ 𝐒]

_mean(𝐒::Union{CrossSpectraVector, CoherenceVector}, frange::UnitRange{Int64}) =
    typeof(𝐒[1].y)([_mean(S, frange) for S ∈ 𝐒])

"""
```julia
(1)
function extract(S :: FDobjects,
            frange :: fInterval)

(2)
function extract(𝐒 :: FDobjectsVector,
            frange :: fInterval;
        w :: Vector = [],
    check :: Bool   = true)

(3)
function extract(Y :: TFobjects,
            frange :: fInterval,
            trange :: tInterval)

(4)
function extract(𝒀 :: TFobjectsVector,
            frange :: fInterval,
            trange :: tInterval;
        w :: Vector = [],
    check :: Bool   = true)
```

**alias**: `extr`

Extract data in a frequency region from [FDobjects](@ref) and data in a
time-frequency region from [TFobjects](@ref). The frequency and time region
are indicated by `frange` and `trange`, which are of type [fInterval](@ref)
and [tInterval](@ref), respectively.

The input/output types of this function for a region with more then one
frequency and more than one sample is reported in the following table:

|method|        input object          |                   output                |
|:---:|:------------------------------|:---------------------------------------------|
|(1.1)| [Spectra](@ref)               | a real matrix with spectra in `frange` arranged in columns¹|
|(1.2)| [CrossSpectra](@ref)          | a vector of complex matrices holding the cross-spectra in `frange`²|
|(1.3)| [Coherence](@ref)             | a vector of real matrices holding the coherence in `frange`²|
|(2.1)| [SpectraVector](@ref)         | a vector of matrices of type (1.1)|
|(2.2)| [CrossSpectraVector](@ref)    | a vector of vectors of type (1.2)|
|(2.3)| [CoherenceVector](@ref)       | a vector of vectors of type (1.3)|
|(3.1)| [TFAnalyticSignal](@ref)      | a complex matrix holding the analytic signal in [`frange`, `trange`]|
|(3.2)| [TFAmplitude](@ref)           | a real matrix holding the amplitude in [`frange`, `trange`]|
|(3.3)| [TFPhase](@ref)               | a real matrices holding the phase in [`frange`, `trange`]|
|(4.1)| [TFAnalyticSignalVector](@ref)| a vector of matrices of type (3.1)|
|(4.2)| [TFAmplitudeVector](@ref)     | a vector of matrices of type (3.2)|
|(4.3)| [TFPhaseVector](@ref)         | a vector of matrices of type (3.3)|
Legend: ¹ *each column refers to a time-series on which the spectra have been computed.*
² *depending on how the objects has been created, the matrices may be either
Hermitian or LowerTriangular. See the documentation of [CrossSpectra](@ref) and [Coherence](@ref).

Note that, depending on the arguments, the type of the output may loose one or two dimensions.
For instance,
- if the [Spectra](@ref) object holds only one spectrum, (1.1) will output a column vector and (2.1) a vector of column vectors.
- if `frange` points to a single frequency, (1.1) will output a row vector and (2.1) a vector of row vectors.
- if both the above two conditions hold, (1.1) will output a real number and (2.1) a vector.
- if `frange` points to a single frequency, (1.2), (1.3) will output a matrix and (2.2), (2.3) a vector of matrices.
- If `frange` points to a single frequency band, (3.1), (3.2), (3.3) will output a row vector and (4.1), (4.2), (4.3) a vector of row vectors.
- If `trange` points to a single time sample, (3.1), (3.2), (3.3) will output a column vector and (4.1), (4.2), (4.3) a vector of column vectors.
- if both the above two conditions hold, (3.1), (3.2), (3.3) will output a number and (4.1), (4.2), (4.3) a vector.

Method (2) and (4) allows the following *optional keyword arguments*:

`w`, a ``k``-vector of non-negative integers or real numbers, where ``k`` is the
numbers of objects hold in the input [FDobjectsVector](@ref) or
[TFobjectsVector](@ref). `w` is a vector of weights for the regions extracted from
the input objects. By default, no weights are assigned.

`check`, a boolean. If it is true (default), it is checked that the non-data fields
of the input objects are all the same (for example, sampling rate, bandwidth, etc.).
Set it to false to improve speed.

**See also**: [`mean`](@ref).

**Examples**:
```julia
using FourierAnalysis

# example with univariate Spectra objects (one series -> one spectrum)
sr, t, f, a = 128, 256, 10, 1
# create a sinusoidal superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# compute univariate spectra
Σ=spectra(v, sr, t)
# spectra in between 8Hz and 12Hz
s=extract(Σ, (8, 12))
# spectra in between 8Hz and 12.5Hz
s=extract(Σ, (8, 12.5))
# spectra at 10Hz
s=extract(Σ, 10) # or s=extract(S, (10, 10))
# these two expressions are equivalent: s=extract(Σ, :), s=Σ.y

# example with multivariate spectra (several series -> several spectra)
Σ=spectra(hcat(v, v+randn(t*16)), sr, t)
# spectra in between 8Hz and 12Hz
S=extract(Σ, (8, 12))
# spectra at 10Hz
S=extract(Σ, 10)

# example with CrossSpectra objects (the same goes for Coherence objects)
X=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
Σ=crossSpectra(X, sr, t)
# cross-spectra in between 8Hz and 12Hz (Hermitian matrices)
S=extract(Σ, (8, 12))
Σ=crossSpectra(X, sr, t; tril=true)
# cross-spectra in between 8Hz and 12Hz (LowerTriangular matrices)
S=extract(Σ, (8, 12))

# example with multiple cross-spectra
X2=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
Σ=crossSpectra([X, X2], sr, t) # this is a CrossSpectraVector
S=extract(Σ, (8, 12); w=[0.4, 0.6])
# now S[1] holds the cross-spectra in range 8-12Hz for X
# and S[2] holds the cross-spectra in range 8-12Hz for X2

# example with time-frequency objects
# (work in the same way for TFAnalyticSignal, TFAmplitude and TFPhase)
Y = TFanalyticsignal(v, sr, t)
# analytic signal within frequencies 8Hz and 12Hz and time samples 1 to 64.
AS=extract(Y, (8, 12), (1, 64))

# all analytic signal within frequencies 8Hz and 12Hz.
AS=extract(Y, (8.0, 12), :) # accept integers and reals for frequencies

# all analytic signal within time samples 1 to 64.
AS=extract(Y, :, (1, 64))

# example with multiple time-frequency objects
# (notice how the type of the output changes)
Y = TFanalyticsignal([v, v+randn(t*16)], sr, t)
AS=extract(Y, (8, 12), (1, 64))
AS=extract(Y, (8), :)
AS=extract(Y, 8, 2)
```
"""
extract(S::FDobjects, frange::fInterval) =
    if (_frange=_getfrange(S, frange, "extract")) ≠ 0
        return _extract(S, _frange)
    end

extract(𝐒::FDobjectsVector, frange::fInterval;
        w::Vector=[], check::Bool=true) =
    if (_frange=_getfrange(𝐒, frange, "extract"; check=check)) ≠ 0
        return isempty(w) ? [_extract(S, _frange) for S ∈ 𝐒] :
                            [ω*_extract(S, _frange) for (ω, S) ∈ zip(w, 𝐒)]
    end


"""
```julia
(1)
function mean(S :: FDobjects,
         frange :: fInterval)

(2)
function mean(𝐒 :: FDobjectsVector,
         frange :: fInterval;
        w :: Vector = [],
    check :: Bool   = true)

(3)
function mean(Y :: TFobjects,
         frange :: fInterval,
         trange :: tInterval)

(4)
function mean(𝒀 :: TFobjectsVector,
         frange :: fInterval,
         trange :: tInterval;
          w :: Vector = [],
      check :: Bool   = true)
```

Return the mean of data in a frequency region from [FDobjects](@ref) and data in a
time-frequency region from [TFobjects](@ref). The frequency and time region
are indicated by `frange` and `trange`, which are of type [fInterval](@ref)
and [tInterval](@ref), respectively.

The complete input/output types for this function is reported in the following table:

|method|        input object           |                   output                |
|:----:|:------------------------------|:---------------------------------------------|
|(1.1)| [Spectra](@ref)                | a vector holding the mean spectra in `frange`¹|
|(1.2)| [CrossSpectra](@ref)           | a complex matrix holding the mean cross-spectra in `frange`²|
|(1.3)| [Coherence](@ref)              | a real matrix holding the mean coherence in `frange`²|
|(2.1)| [SpectraVector](@ref)          | a vector of vectors of type (1.1)|
|(2.2)| [CrossSpectraVector](@ref)     | a vector of matrices of type (1.2)|
|(2.3)| [CoherenceVector](@ref)        | a vector of matrices of type (1.3)|
|(3.1)| [TFAnalyticSignal](@ref)       | a complex number holding the mean analytic signal in [`frange`, `trange`]|
|(3.2)| [TFAmplitude](@ref)            | a real number holding the mean amplitude in [`frange`, `trange`]|
|(3.3)| [TFPhase](@ref)                | a real number holding the mean phase in [`frange`, `trange`]|
|(4.1)| [TFAnalyticSignalVector](@ref) | a vector of numbers of type (3.1)|
|(4.2)| [TFAmplitudeVector](@ref)      | a vector of numbers of type (3.2)|
|(4.3)| [TFPhaseVector](@ref)          | a vector of numbers of type (3.3)|
legend: ¹*each element of the vector refers to a time-series on which the spectra have been computed.*
² *depending on how the objects has been created, the matrices may be either
Hermitian or LowerTriangular.*

Method (2) and (4) allows the following *optional keyword arguments*:

`w`, a ``k``-vector of non-negative integers or real numbers, where ``k`` is the
numbers of objects hold in the input [FDobjectsVector](@ref) or
[TFobjectsVector](@ref). `w` is a vector of weights for the means extracted from
the input objects. By default, no weights are assigned.

`check`, a boolean. If it is true (default), it is checked that the non-data fields
of the input objects are all the same (for example, sampling rate, bandwidth, etc.).

**See also**: [`extract`](@ref).

**Examples**:
```julia
using FourierAnalysis, Plots

# example with univariate Spectra objects (one series -> one spectrum)
sr, t, f, a = 128, 256, 10, 1
# create a sinusoidal superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# compute the spectrum
Σ=spectra(v, sr, t)
# mean spectrum in between 8Hz and 12Hz
s=mean(Σ, (8, 12))
# mean spectrum in between 8Hz and 12.5Hz
s=mean(Σ, (8, 12.5))

# example with multivariate spectra (several series -> several spectra)
Σ=spectra(hcat(v, v+randn(t*16)), sr, t)
# mean spectra in between 8Hz and 12Hz
S=mean(Σ, (8, 12))
# mean spectra at 10Hz, i.e., the spectra at 10Hz
S=mean(Σ, 10)

# example with CrossSpectra objects (the same goes for Coherence objects)
X=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
Σ=crossSpectra(X, sr, t)
# mean cross-spectra in between 8Hz and 12Hz (an Hermitian matrix)
S=mean(Σ, (8, 12))
Σ=crossSpectra(X, sr, t; tril=true)
# mean cross-spectra in between 8Hz and 12Hz (a LowerTriangular matrix)
S=mean(Σ, (8.0, 12.0)) # accept integers and reals for frequencies

# example with multiple CrossSpectra objects
X2=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
Σ=crossSpectra([X, X2], sr, t) # this is a CrossSpectraVector
S=mean(Σ, (8, 12); w=[0.4, 0.6])
# now S[1] will hold the mean cross-spectrum in range 8-12Hz for X
# and S[2] will hold the mean cross-spectrum in range 8-12Hz for X2

# example with time-frequency objects
# (work in the same way for TFAnalyticSignal, TFAmplitude and TFPhase)
Y = TFanalyticsignal(v, sr, t)
# mean analytic signal within frequencies 8Hz and 12Hz and time samples 1 to 64.
as=mean(Σ, (8, 12), (1, 64))
# mean analytic signal within frequencies 8Hz and 12Hz.
as=mean(Σ, (8, 12), :)
# mean analytic signal within time samples 1 to 64.
as=mean(Σ, :, (1, 64))

# example with multiple time-frequency objects
Y = TFanalyticsignal([v, v+randn(t*16)], sr, t)
AS=mean(Y, (8, 12), (1, 64))
# get the mean across TFobjects of those means
m=mean(mean(Y, (8, 12), (1, 64)))
AS=mean(Y, (8), :)
AS=mean(Y, 8, 2)
```
"""
mean(S::FDobjects, frange::fInterval) =
    if (_frange=_getfrange(S, frange, "mean")) ≠ 0
        return _mean(S, _frange)
    end

mean(𝐒::FDobjectsVector, frange::fInterval;
     w::Vector=[], check::Bool=true) =
    if (_frange=_getfrange(𝐒, frange, "mean"; check=check)) ≠ 0
        return isempty(w) ? [_mean(S, _frange) for S ∈ 𝐒] :
                            [ω*_mean(S, _frange) for (ω, S) ∈ zip(w, 𝐒)]
    end

### Time-Frequency-domain data extraction and averaging ###
###############################################################################

function _getfrange(Y::TFobjects, frange::fInterval, funcname::String)
    # find frequencies in filterbanks
    hb=Y.bandwidth/2
    if      isa(frange, IntOrReal)
            (Y.flabels[1]-hb)<=frange<=(Y.flabels[end]+hb) ? (return max(1, (findmin([abs(f-frange) for f∈Y.flabels])[2])):min(length(Y.flabels), (findmin([abs(f-frange)  for f∈Y.flabels])[2]))) :
            @error 📌*", "*funcname*" function passed invalid frange Int or Real argument. The frequency must be comprised between $((Y.flabels[1]-hb)) and $((Y.flabels[end]+hb))" frange
            return 0
    elseif  isa(frange, Colon)
            return 1:size(Y.y, 1)
    elseif  isa(frange, Tuple{IntOrReal, IntOrReal})
        if      frange[1]>frange[2]
                    @error 📌*", "*funcname*" function passed invalid frange 2-tuple argument. The second element cannot be inferior to the first" frange[1] frange[2]
                    return 0
        elseif  frange[1]<(Y.flabels[1]-hb) || frange[2]>(Y.flabels[end]+hb)
                    @error 📌*", "*funcname*" function passed invalid frange 2-tuple argument. The two frequencies must be comprised between $(Y.flabels[1]-hb) and $(Y.flabels[end]+hb)" frange[1] frange[2]
                    return 0
        else        return max(1, (findmin([abs(f-frange[1]) for f∈Y.flabels])[2])):min(length(Y.flabels), (findmin([abs(f-frange[2])  for f∈Y.flabels])[2]))
        end
    end
end


function _gettrange(Y::TFobjects, trange::tInterval, funcname::String)
    if      isa(trange, Int)
            1<=trange<=size(Y.y, 2) ? (return trange:trange) :
            @error 📌*", "*funcname*" function passed invalid trange Int argument. The time sample must be comprised between 1 and $(size(Y.y, 2))" trange
            return 0
    elseif  isa(trange, Colon)
            return 1:size(Y.y, 2)
    elseif  isa(trange, Tuple{Int, Int})
        if      trange[1]>trange[2]
                    @error 📌*", "*funcname*" function passed invalid trange 2-tuple argument. The second element cannot be inferior to the first" trange[1] trange[2]
                    return 0
        elseif  trange[1]<1 || trange[2]>size(Y.y, 2)
                    @error 📌*", "*funcname*" function passed invalid trange 2-tuple argument. The time limits in samples must be comprised between 1 and $(size(Y.y, 2))" trange[1] trange[2]
                    return 0
        else        return trange[1]:trange[2]
        end
    end
end

_getftrange(Y::TFobjects, frange::fInterval, trange::tInterval, funcname::String) =
    _getfrange(Y, frange, funcname), _gettrange(Y, trange, funcname)


function extract(Y::TFobjects, frange::fInterval, trange::tInterval)
    _frange, _trange = _getftrange(Y, frange, trange, "extract")
    if _frange ≠ 0 && _trange ≠ 0
        _isnull(_frange, _trange) ? Y.y[first(_frange), first(_trange)] :
                                    copy(Y.y[_frange, _trange])
    end
end

function mean(Y::TFobjects, frange::fInterval, trange::tInterval)
    _frange, _trange = _getftrange(Y, frange, trange, "mean")
    if _frange ≠ 0 && _trange ≠ 0
        return mean(view(Y.y, _frange, _trange))
    end
end

_getftrange(𝒀::TFobjectsVector, frange::fInterval, trange::tInterval, funcname::String;
            check::Bool=true) =
    check ? (sameParams(𝒀, funcname) ? _getftrange(𝒀[1], frange, trange, funcname) : 0) :
    _getftrange(𝒀[1], frange, trange, funcname)

extract(𝒀::TFobjectsVector, frange::fInterval, trange::tInterval;
        w::Vector=[], check::Bool=true) =
    if      check && frange isa Colon && !_allsame([size(Y.y, 1) for Y ∈ 𝒀])
            @error 📌*", extract function. In order of using a column (:) as `frange` argument, the number of rows (frequencies) must be identical in all time-frequency objects"
    elseif  check && trange isa Colon && !_allsame([size(Y.y, 2) for Y ∈ 𝒀])
            @error 📌*", extract function. In order of using a column (:) as `trange` argument, the number of columns (time samples) must be identical in all time-frequency objects"
    elseif  ((_frange, _trange)=_getftrange(𝒀, frange, trange, "extract"; check=check)) ≠ (0, 0)
            if _isnull(_frange, _trange)
                isempty(w) ? [Y.y[first(_frange), first(_trange)] for Y ∈ 𝒀] :
                             [ω*Y.y[first(_frange), first(_trange)] for (ω, Y) ∈ zip(w, 𝒀)]
            else
                isempty(w) ? [Y.y[_frange, _trange] for Y ∈ 𝒀] :
                             [ω*Y.y[_frange, _trange] for (ω, Y) ∈ zip(w, 𝒀)]
            end
    end


mean(𝒀::TFobjectsVector, frange::fInterval, trange::tInterval;
     w::Vector=[], check::Bool=true) =
    if      check && frange isa Colon && !_allsame([size(Y.y, 1) for Y ∈ 𝒀])
            @error 📌*", extract function. In order of using a column (:) as `frange` argument, the number of rows (frequencies) must be identical in all time-frequency objects"
    elseif  check && trange isa Colon && !_allsame([size(Y.y, 2) for Y ∈ 𝒀])
            @error 📌*", extract function. In order of using a column (:) as `trange` argument, the number of columns (time samples) must be identical in all time-frequency objects"
    elseif  ((_frange, _trange)=_getftrange(𝒀, frange, trange, "mean"; check=check)) ≠ (0, 0)
            return isempty(w) ? [mean(view(Y.y, _frange, _trange)) for Y ∈ 𝒀] :
                                [ω*mean(view(Y.y, _frange, _trange)) for (ω, Y) ∈ zip(w, 𝒀)]
    end

extr=extract

### Band-pass Averages ###

"""
```julia
function bands(S :: Union{FDobjects, FDobjectsVector}
       bandwidth :: IntOrReal)
```

Return band-pass average of spectral, cross-spectral or coherence estimates
in equally spaced band-pass regions with the given `bandwidth`.
`bandwidth` can be given as an integer or as a real number. See [`bbands`](@ref)
for details on the definition of band-pass regions.

Band-pass average is not supported for time-frequency objects as for those
objects a similar averaging is natively avaiable using argument `bandwidth`
in their constructors.

The output of this function is as it follows:
- for univariate [Spectra](@ref) objects (i.e., those hodling one spectrum only), a real column vector,
- for multivariate [Spectra](@ref) objects, a real matrix,
- for [SpectraVector](@ref) objects, a vector of the above,
- for [CrossSpectra](@ref) and [Coherence](@ref) objects, a vector of Hermitian or LowerTriangular matrices, depending on how the object has been cosntructed,
- for [CrossSpectraVector](@ref) and [CoherenceVector](@ref) objects, a vector  of the above.

**See**: [`bbands`](@ref).

**Examples**:
```julia
using FourierAnalysis, Plots

# example with univariate Spectra objects (one series -> one spectrum)
sr, t, f, a = 128, 256, 10, 1
# create a sinusoidal superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# compute the spectrum
Σ=spectra(v, sr, t)
# mean spectra in 2Hz-band-pass regions
b=bands(Σ, 2)
plot(b)

# example with multivariate spectra (several series -> several spectra)
Σ=spectra(hcat(v, v+randn(t*16), v+randn(t*16) ), sr, t)
# mean spectra in 2Hz-band-pass regions for all time-series
b=bands(Σ, 2)
plot(b)
# plot mean spectra in 2Hz-band-pass regions for time-series 2 and 3 only
plot(bands(Σ, 2)[:, 2:3])

# example with CrossSpectra objects (the same goes for Coherence objects)
X=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
Σ=crossSpectra(X, sr, t)
# mean cross-spectra in 4Hz-band-pass regions
B=bands(Σ, 4)

# example with multiple CrossSpectra objects
X2=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
Σ=crossSpectra([X, X2], sr, t) # this is a CrossSpectraVector
# mean cross-spectra in 4Hz-band-pass regions for all cross-spectra objects
B=bands(Σ, 4)
```
"""
function bands(S::Spectra, bandwidth::IntOrReal)
    if (bands=bbands(S.sr, S.wl, bandwidth; DC=S.DC)) isa Nothing return S end
    r, n=length(bands)-1, size(S.y, 2)
    if n==1 return [mean(S.y[bands[b]:bands[b+1]]) for b=1:r]
    else
        Y=Matrix{eltype(S.y)}(undef, r, n)
        @inbounds for i=1:n Y[:, i]=[mean(S.y[bands[b]:bands[b+1], i]) for b=1:r] end
        return Y
    end
end # input Spectra, Output a vector or a Matrix


bands(𝐒::SpectraVector, bandwidth::IntOrReal) =
    [bands(S, bandwidth) for S ∈ 𝐒] # output a vector of vectors or matrices

function bands(S::Union{CrossSpectra, Coherence}, bandwidth::IntOrReal)
    if (bands=bbands(S.sr, S.wl, bandwidth; DC=S.DC)) isa Nothing return S end
    r=length(bands)-1
    mattype = S.tril ? LowerTriangular : Hermitian
    return typeof(S.y)([mattype(mean(S.y[i] for i=bands[b]:bands[b+1])) for b=1:r])
end # output: a vector of LowerTriangular (if tril=true) or Hermitian (if tril=false) matrices


bands(𝐒::Union{CrossSpectraVector, CoherenceVector}, bandwidth::IntOrReal) =
    [bands(S, bandwidth) for S ∈ 𝐒] # output: a vector of vectors of LowerTriangular (if tril=true) or Hermitian (if tril=false) matrices


# internal function: used by functions that may execute in multi-threading
# to decide to so so or not.
# The multi-threading is to be done for a for looping over `n` elements.
# Return true if `n` is large enough to justify multi-threading (n>thr, where
# `thr` is the number of threads Julia is instructed to use) and thr>1,
# otherwise return false.
function _thread(⏩, n)
    thr = nthreads() # of threads Julia is instructed to use
    ⏩ && n>=thr*2 && thr > 1 ? (return true) : (return false)
end

# check that all elements of a vector are equal. General function
# taken from https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal
_allsame(x) = all(y -> y == (first(x)), x)

_isnull(frange::UnitRange{Int64}) = first(frange)==last(frange)
_isnull(frange::UnitRange{Int64}, trange::UnitRange{Int64}) =
    (first(frange)==last(frange)) & (first(trange)==last(trange))
