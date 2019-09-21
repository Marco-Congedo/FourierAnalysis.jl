# fftw.jl

This unit implements a convenient interface to
[FFTW.jl](https://github.com/JuliaMath/FFTW.jl) and
[AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl).
It is not part of the main *FourierAnalysis* API,
since all its functions manage basic FFTW plans
for FFT and iFFT computations using default settings.
However, the use of this unit is important if massive FFT computations are to be performed.

In order to use effectively FFTW, see [window length in FFTW](@ref).

The following are constants used by FFTW as flags in creating FFT plans.
See [here](http://www.fftw.org/fftw3_doc/Planner-Flags.html) for details.

```
plan_estimate        = UInt32(64)      # (1 << 6)   very small search
plan_measure         = UInt32(0)       #            small search
plan_patient         = UInt32(32)      # (1 << 5)   large search
plan_exhaustive      = UInt32(8)       # (1 << 3)   very large search
plan_conserve_memory = UInt32(4)       # (1 << 2)
plan_wisdom_only     = UInt32(2097152) # (1 << 21)
plan_unaligned       = UInt32(2)       # (1 << 1)
plan_preserve_input  = UInt32(16)      # (1 << 4)
```

```@docs
Planner
```
