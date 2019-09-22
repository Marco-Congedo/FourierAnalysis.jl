#   Unit "hilbert" of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 5th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements Welch high-resolution Hilbert transform and
#   analytic signal estimations using FFTW

"""
```
   (1)
   function analyticsignal( X  :: Union{Vector{T}, Matrix{T}},
                            wl :: Int     = size(X, 1);
                     nonlinear :: Bool    =  false,
                     planner   :: Planner =  getplanner,
                     ⏩       :: Bool    =  true) where T<:Real

   (2)
   function analyticsignal( 𝐗      :: Vector{Matrix{T}},
                            wl     :: Int;
                        nonlinear  ::  Bool      =  false,
                        planner    ::  Planner   =  getplanner,
                        ⏩        ::  Bool      =  true) where T<:Real
```

(1)

Compute the analytic signal(AS) of vector `X` or of
all column vectors of matrix `X` via the FFT and iFFT procedure,
as explained in Marple(1999).
If `wl`=size(X, 1) (default), use the standard method passing to the FFT
and iFFT all samples in `X` altogether, whereas if `wl`<size(X, 1) a Welch-like
method is used (see below).

Return the analytic signal `𝑌`, a complex vector if `X` is a vector
or a complex matrix holding in its columns the analytic signal of the
columns of `X` if `X` is a matrix. `𝑌` has the same number of samples (rows) as `X`.
Contrarely to what is done in the [DSP](https://github.com/JuliaDSP/DSP.jl)
package, the DC level of the signal is removed, therefore,
if the input signal features a non-zero DC level,
the real part of the AS will be equal to the input signal with the
DC level removed. The imaginary part of `𝑌` is the Hilbert transform
of such no-DC `X`.

The Welch-like AS allows an efficient estimation of the AS for vectors and
matrices of any length, that is, even if they are very large; it proceeds
computing the AS on 50% sliding overlapping windows and forming the AS
by retaining the central half of each window. The number of points effectively
used to obtain the final estimation is ``wl÷2`` (integer division)
if `wl` is even, ``wl÷2+1`` otherwise. Using the central half of the
AS computed at each window eliminates edge effects in the range ``wl÷2:end-wl÷2``.
Before running the procedure, ``wl÷2`` zeros are padded at the beginning
and at the end of `X` and trimmed on exiting the function,
in order to reduce edge effects at the beginning
and end of `X`. This procedure does not eliminate these edge effects
completely, thus the first and last ``wl÷2`` samples of the AS estimation
should be discarded. In practice, one requires the AS of a larger data
segment and trims at least ``wl÷2`` samples at the beginning and end of
the estimation. This is done automatically by the [`TFanalyticsignal`](@ref)
function.

!!! note "Nota Bene"
    In order to avoid FFT computation of very long epochs,
    if `wl` > 2^14, then `wl` is set to 2^10.

    See [window length in FFTW](@ref) for setting efficiently argument `wl`.

    The input signal should be previously band-pass or high-pass filtered
    so as not to contain frequency components below the first discrete
    Fourier frequency obtained with windows of `wl` samples, that is,
    below sr/wl Hz, where sr is the sampling rate.

**Optional Keyword Arguments**:

`nonlinear`, if true, the analytic signal is normalized so that its amplitude
is ``1.0`` at all points. This allow non-linear univariate and bivariate estimations
(see [timefrequencyuni.jl](@ref) and [timefrequencybi.jl](@ref)).

`planner` is an instance of the [`Planner`](@ref) object, holding the forward
and backward FFTW plans used to compute the FFTs and the iFFTs.
By default the planner is computed, but it can be passed as an
argumet here if it is pre-computed. This is interesting if the
`analyticsignal` function is to be invoked repeatedly.

if `⏩` is true, the method is run in multi-threaded mode across the series in `X`
if the number of series is at least twice the number of threads Julia
is instructed to use. See [Threads](@ref).

(2)

Compute the analytic signal for all ``k``
multivariate data matrices given as a vector of matrices `𝐗`.
Return a vector of matrices hodling the corresponding analytic signals
as in method (1). The FFT and iFFT plans are computed only once.
The ``k`` matrices in `𝐗` may have different number of columns (i.e., different
number of series) and different number of rows (samples).
However, the number of rows must be larger than `wl` for all of them.
if `⏩` is true, this method run in multi-threaded mode across the
matrices in `𝐗` if the number of matrices is at least twice
the number of threads Julia is instructed to use, otherwise it
tries to run each analytic signal estimation in multi-threaded mode
as per method (1). See [Threads](@ref).

This function is called by the following functions operating
on time-frequency reprsentations: [`TFanalyticsignal`](@ref),
[`TFamplitude`](@ref), [`TFphase`](@ref), [`meanAmplitude`](@ref),
[`concentration`](@ref), [`meanDirection`](@ref), [`comodulation`](@ref),
[`coherence`](@ref).

**References**
Marple S.L. (1999)
Computing the Discrete-Time Analytic Signal via FFT.
IEEE Transactions on Signal Processing 47(9), 2600-3.

**Examples**:
```
using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, DSP
t=128; lab=["x", "real(y)", "imag(y)"]

# Analytic signal of one vector
x=sinusoidal(10, 2, 128, t, π/2; DC=10) # cosine
y=analyticsignal(x)
# note that real(y) is x without the DC level, i.e., x=real(y)+DC
plot([x, real(y), imag(y)]; labels=lab)

# make a check
s=sinusoidal(10, 2, 128, t, 0) # sine
norm(s-imag(y)) # should be zero

# Analytic Signal by DSP.jl
y2=hilbert(x)
norm(s-imag(y2)) # should be zero
# DSP.jl does not remove the DC level
# thus x=real(y2) in this case
plot([x, real(y2), imag(y2)]; labels=lab)

# Analytic signal of multiple vectors
x=hcat(x, sinusoidal(10, 3, 128, t, π/2; DC=10))
y=analyticsignal(x)

# welch-like analytic signal of one vector
# (note edge effects)
x=sinusoidal(10, 2, 128, t*4, π/2; DC=0)
y=analyticsignal(x, t)
plot([x, real(y), imag(y)]; labels=lab)

# Welch-like analytic signal of multiple vectors
x=hcat(x, sinusoidal(10, 3, 128, t*4, π/2; DC=0))
y=analyticsignal(x, t)
```
"""
function analyticsignal( X  :: Union{Vector{T}, Matrix{T}},
                         wl :: Int     = size(X, 1);
                  nonlinear :: Bool    =  false,
                  planner   :: Planner =  getplanner,
                  ⏩       :: Bool    =  true) where T<:Real

   # get all needed paramters (see below the _paramHT! function)
   X_, 𝚙, i𝚙, e, tₐₗₗ, n, wl½, two_wl⁻¹, cT, f, g, ζ = _paramHT!(X, wl, planner)

   𝑌=zeros(cT, tₐₗₗ, n) # Preallocate memory for the Analytic Signal of X

   # Add to `𝑌` the analytic signal of a vector `X_` or of the column vectors of
   # matrix `X_`, computed on the epoch(s) of duration `t` samples starting at
   # sample `at`. It does not return anything.
   # ARGUMENTS:
   # `n` (Int) is the number of columns in `X_`; it must be 1 if `X_` is a vector.
   # `wl½` = wl÷2
   # `two_wl⁻¹` = 2/wl
   # `e` (Int) the number of sliding windows
   # `𝚙` is the forward FFTW plan performing the FFT.
   # `i𝚙` is the backward FFTW plan performing the iFFT (Hilber transform).
   # `ζ`, =zeros(cT, wl-wl½_), a zero-vector append to the iFFT vectors, since the second half of the FFT is not computed
   # `f`, the lower limit (in samples) of the central half of each Hilbert transform to cumulate with respect to the FFT window
   # `g`, the upper limit (in samples) of the central half of each Hilbert transform to cumulate with respect to the FFT window
   function _analyticsignal!(at, 𝚗)
      y = 𝚙*X_[at:at+wl-1, 𝚗]        # FFT
      @inbounds begin                # transform
         y[1]=0.
         for i=2:wl½ y[i] *= two_wl⁻¹ end
      end
      e==1 ? 𝑌[:, 𝚗]+=i𝚙*vcat(y, ζ) : 𝑌[at+f-1:at+g-1, 𝚗]+=(i𝚙*vcat(y, ζ))[f:g]  # iFFT, cumulate
   end

   # add the Analytic Signal
   as!(𝚎, 𝚗)=_analyticsignal!((𝚎*wl½)+1, 𝚗)
   _thread(⏩, n) ? (for 𝚎=0:e-1 @threads for 𝚗=1:n as!(𝚎, 𝚗) end end) : (for 𝚗=1:n, 𝚎=0:e-1 as!(𝚎, 𝚗) end)

   # return a vector if `X` is a vector, a matrix if `X` is a matrix.
   # for the high-resolution version eliminate the first and last t½ samples that have been previously padded.
   𝑌_ = n==1 ? (e>1 ? 𝑌[:][wl½+1:end-wl½] : 𝑌[:]) : (e>1 ? 𝑌[wl½+1:end-wl½, 1:n] : 𝑌)
   if nonlinear @inbounds for i in eachindex(𝑌_) 𝑌_[i]/=abs(𝑌_[i]) end end
   return 𝑌_
end


function analyticsignal( 𝐗      :: Vector{Matrix{T}},
                         wl     :: Int;
                     nonlinear  ::  Bool      =  false,
                     planner    ::  Planner   =  getplanner,
                     ⏩        ::  Bool      =  true) where T<:Real

   cT, k=typeof(Complex(T(0))), length(𝐗)

   planner ≠ getplanner ? plan=planner : plan=Planner(plan_estimate, -1.0, wl, T, true)

   𝒀=Vector{Matrix{cT}}(undef, k)

   # function to add the analytic signal
   as(i, ⏩)=analyticsignal(𝐗[i], wl; nonlinear=nonlinear, planner=plan, ⏩=⏩)
   _thread(⏩, k) ? (@threads for i=1:k 𝒀[i]=as(i, false) end) : (for i=1:k 𝒀[i]=as(i, ⏩) end)

   return 𝒀
end


# Internal function: prepare all parameters that are needed for
# Analytic Signal and High-Resolution Analytic Signal estimations
#
# ARGUMENTS:
# `X` is the input data vector or matrix (n time-series)
# `inrange` is the range (UnitRange{Int}, e.g., 128:1280-1) of the vectors
#      or of rows of `X` on which the analysis is to be performed
# `wl` (Int) is the epoch length for FFT
# `plan` a FFTW.rFFTWPlan forward plan for computing the FFT
# `plan` a FFTW.cFFTWPlan backward plan for computing the FFT
# `planflag` flags for the computation of the plan by FFYW. See doc of `spectra`
# `plantime` (Int) is the maximum time allowed (in seconds) to optimize the FFTW plan
# `⏩` if true multi-threaded computations are requested
#
# RETURN
# `X_`: if a standard anaytic signal is requested (t==size(X, 1)) then X_=X
#       else X_ is the data in `X` with `t`-1 zeros prepended and appended
# `𝚙` return the forward FFTW plan. It is computed if argument `plan`=`computeplan`
# `i𝚙` return the backward FFTW plan. It is computed if argument `iplan`=`computeiplan`
# `e` (Int) the number of sliding windows for the validated `inrange`, always with 1-sample step
# `n` (Int) 1 if `X` is a Vector, the number of columns of `X` otherwise
# `from` the first sample in the validated `inrange`
# `wl½` `t`÷2 computed by the right aritmetic bit-wise operation `wl`>>1
# `wl½_`, the length of FFT vectors:  t½+1 if DC=true, t½ is DC=false
# `two_wl⁻¹` = 2/wl
# `cT` the complex type corresponding to type `T` in function definition (e.g., ComplexF64 if T=Float64)
# `thrn`: true if ⏩=true and # of threads is >1 and is >= 2*`n`, false otherwise
# `navg`, the number of analytic signal estimations for each sample. This is
#       different from `e` because only the central half of each Hilbert transform is used for each window
# `f`, the lower limit (in samples) of the central half of each Hilbert transform to cumulate with respect to the FFT window
# `g`, the upper limit (in samples) of the central half of each Hilbert transform to cumulate with respect to the FFT window
# `ζ`, =zeros(cT, wl-wl½_), a zero-vector append to the iFFT vectors, since the second half of the FFT is not computed
# NB if wl > 2^14 then t is set to 2^10. This affects the default behavior of all AS functions
function _paramHT!( X       :: AbstractArray{T},
                    wl      :: Int,
                    planner :: Planner) where T<:Real

    size(X, 1) < wl && @error 📌*", the number of samples in input matrix is smaller than the desired window length."
    wl > 2^14 ? wl = 2^10 : nothing
    n, wl½, two_wl⁻¹ = size(X, 2), (wl>>1), T(2/wl)
    if size(X, 1)>wl
        X_=[zeros(T, wl½, n); X; zeros(T, wl½, n)] else X_=X
        isodd(wl) && @error 📌*", for Welch-like Analytic Signal estimation `wl` must be even."
    end
    tₐₗₗ = size(X_, 1) # for high-resolution anaytic signal, this is not size(X, 1)
    cT = typeof(Complex(T(0)))
    if planner ≠ getplanner
       𝚙=planner.p
       i𝚙=planner.ip
    else
       𝚙=plan_rfft(zeros(T, wl), flags=plan_estimate, timelimit=-1.0)
       i𝚙=plan_bfft(zeros(cT, wl), flags=plan_estimate, timelimit=-1.0)
    end
    f=wl½÷2+1 # lower limit of central region to copy
    g=f+wl½-1 # upper limit of central region to copy
    e = (tₐₗₗ-wl)÷wl½+1 # number of 50% overlapping epochs
    ζ=zeros(cT, wl-wl½-1)

    return X_, 𝚙, i𝚙, e, tₐₗₗ, n, wl½, two_wl⁻¹, cT, f, g, ζ
end
