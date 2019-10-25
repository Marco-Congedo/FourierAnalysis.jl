# tools.jl

This unit implements
- [functions](@ref) that are useful when employing the FFT,
- [specific methods](@ref): they apply to Julia's types or to single objects created by *FourierAnalysis*,
- [generic methods](@ref): they apply to entire categories of objects created by *FourierAnalysis*.

## functions

|         function         |           description             |
|:-------------------------|:----------------------------------|
| [`sinusoidal`](@ref) | generate a sinusoidal wave |
| [`fres`](@ref)       | FFT frequency resolution |
| [`f2b`](@ref)        | bin on a real-FFT vector best matching a frequency |
| [`b2f`](@ref)        | frequency (in Hz) that correspond to a bin on a real-FFT vector|
| [`fdf`](@ref)        | all Fourier discrete frequencies for a real-FFT |
| [`brange`](@ref)     | range of bins for a real-FFT vector covering all Fourier discrete frequencies |
| [`bbands`](@ref)     | limits of all bandwidth-spaced band-pass regions of a real-FFT vectors, in bins |
| [`fbands`](@ref)     | limits of all bandwidth-spaced band-pass regions of a real-FFT vectors, in frequencies (Hz) |
| [`dB`](@ref)         | convert a measure or a ratio between two measures into deciBels |

```@docs
sinusoidal
fres
f2b
b2f
fdf
brange
bbands
fbands
dB
```

## specific methods

|      function        |           description             |
|:---------------------|:----------------------------------|
| [`amplitude`](@ref)  | return the amplitude (modulus) of a complex number, complex array or [TFAnalyticSignal](@ref) object|
| [`phase`](@ref)      | return the phase (argument) of a complex number, complex array or  [TFAnalyticSignal](@ref) object|
| [`polar`](@ref)      | return the phase (argument) and amplitude (modulus) of a complex number, a complex array or [TFAnalyticSignal](@ref) object|
| [`unwrapPhase`](@ref)| compute and unwrap the phase of a complex array, unwrap a real array holding phase in [−π, π], or construct a [TFPhase](@ref) object with the phase unwrapped|
| [`isUnwrapped`](@ref)| return true if the phase of all objects in a [TFPhaseVector](@ref) is unwrapped|

```@docs
amplitude
phase
polar
unwrapPhase
isUnwrapped
```

## generic methods

**Generic methods** applying to [FDobjects](@ref), [FDobjectsVector](@ref), [TFobjects](@ref) and [TFobjectsVector](@ref):

|      function     |           description             |
|:------------------|:----------------------------------|
| [`smooth`](@ref)  | smooth the data across frequencies and/or across time    |
| [`extract`](@ref) | extract the data in a frequency or time-frequency region |
| [`mean`](@ref)    | compute the mean in a frequency or time-frequency region |
¤

**Generic methods** applying only to [FDobjects](@ref) and [FDobjectsVector](@ref):

|     function    |           description             |
|:----------------|:----------------------------------|
| [`bands`](@ref) |  Return band-pass average of spectral, cross-spectral or coherence estimates in equally spaced band-pass regions |
¤

**Generic methods** applying to [FDobjectsVector](@ref),
and [TFobjectsVector](@ref):

|       function       |           description             |
|:---------------------|:----------------------------------|
| [`sameParams`](@ref) | return true if the non-data fields of all objects in the vector have the same value |
| [`isLinear`](@ref)   | return true if all objects in the vector are linear |
| [`isNonLinear`](@ref)| return true if all objects in the vector are non-linear |
¤

```@docs
smooth
extract
mean
bands
sameParams
isLinear
isNonLinear
```
