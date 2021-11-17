#   Unit "spectra" of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2021,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements Welch power spectral estimates using FFTW

#############################################################################
# Generic constructors of Spectra objects

# Enable construction giving only `y`, `sr`, `wl`, `DC` and `taper` argument.
# `flabels` is generated automatically, `smoothing` is set to `noSmoother`
# and `func` is set to `identity`
Spectra(y, sr, wl, DC, taper) =
    Spectra(y, sr, wl, DC, taper, fdf(sr, wl; DC=DC), identity, noSmoother)
#
# As above, but setting by default also `DC` to false
Spectra(y, sr, wl, taper) =
    Spectra(y, sr, wl, false, taper, fdf(sr, wl; DC=false), identity, noSmoother)
#
#
# Get the spectra from a CrossSpectra object and construct a Spectra object
# The `func` function can be applied element-wise to the spectra (see `spectra`).
Spectra(𝙎::CrossSpectra; func::Function=identity) =
    Spectra(func.(Real.([𝙎.y[i][j, j] for i=1:length(𝙎.y), j=1:size(𝙎.y[1], 2)])), 𝙎.sr, 𝙎.wl, 𝙎.DC, 𝙎.taper, 𝙎.flabels, func, 𝙎.smoothing)
#
# Get the spectra from a CrossSpectraVector and construct a SpectraVector
# The `func` function can be applied element-wise to the spectra (see `spectra`).
Spectra(𝙎::CrossSpectraVector; func::Function=identity) =
    SpectraVector([Spectra(S; func=func) for S in 𝙎])

#############################################################################
"""
```julia
(1)
function spectra( X     :: Union{Vector{T}, Matrix{T}},
                  sr    :: Int,
                  wl    :: Int;
            tapering  :: Union{Taper, TaperKind} = harris4,
            planner   :: Planner                 = getplanner,
            slide     :: Int                     = 0,
            DC        :: Bool                    = false,
            func      :: Function                = identity, # equal to x->x
            smoothing :: Smoother                = noSmoother,
            ⏩       :: Bool                    = true) where T<:Real

(2)
function spectra( 𝐗   :: Vector{Matrix{T}},
            < same argument sr, ..., ⏩ of method (1) > where T<:Real
```

(1)

Construct a [Spectra](@ref) object from real univariate or
multivariate data. Given sampling rate `sr` and epoch length `wl`,
compute the Welch power spectra of a vector (univariate) or of a data matrix
(multivariate) `X` of dimension ``t``x``n``, where ``t`` is the
number of samples (rows) and ``n`` is the number of time-series (columns).

The spectra are hold in the `.y` field of the created object.
If `X` is a vector, `.y` is a vector, whereas if `X` is a matrix,
`.y` is a matrix holding in its columns the spectra of the
signals in the columns of `X`. The size of the spectra depends on
the `DC` optional keyword argument (see below),
as reported in the documentation of the [Spectra](@ref) structure.

**Optional Keyword Arguments**:

`tapering`: this is a tapering object of type [Taper](@ref) or a tapering
window kind of type [TaperKind](@ref).
By default the *harris4* tapering window is used.
If no tapering is sought, pass as argument `tapering=rectangular`.
This same syntax is the most convenient way to specify all simple
tapering window, e.g., `tapering=hann`, `tapering=hamming`, etc.
For *discrete prolate spheroidal sequences (dpss)* multi-tapering,
use instead the [`slepians`](@ref) constructor, e.g.,
pass as argument something like ```tapering=slepians(sr, wl, 2)```.

`planner`: this is an instance of the [`Planner`](@ref) object, holding the
forward FFTW plan used to compute the FFTs.
By default the planner is computed by this method,
but it can be passed as an
argumet here if it is pre-computed. This is interesting if this
function is to be invoked repeatedly.

`slide`: this is the number of samples the windows slide on (Welch method).
By default the number of samples is chosen to allow 50% overlap.

`DC`: if true, the spectrum/a of the DC level is returned
in the first row of `y` (see the fields of the [Spectra](@ref)
object), otherwise (default) the rows in `y` start with the first positive
discrete frequency, that is, ``sr/wl`` Hz.

`func`: this is a function that operate element-wise to transfrom the power
spectra before returning them, including anonymous functions.
Common choices are:
- `func=sqrt` return the amplitude spectra,
- `func=log` return the log-power spectra,
- `func=decibel` return the power spectra in deciBels (see [`decibel`](@ref)),
By default no function is applied and the power spectra are returned.
If smoothing has been requested (see below), it is applied after the function.

`smoothing`: it applies a smoothing function of type [Smoother](@ref)
to the spectra across adjacent frequencies. By default no smoothing is applied.

`⏩`: if true (default), the method run in multi-threaded mode
across the series in `X` if the number of series is at least twice
the number of threads Julia is instructed to use. See [Threads](@ref).

(2)

Construct a [SpectraVector](@ref) object from a vector of real
univariate (vectors) or multivariate data (matrices).
Compute the spectra as per method (1)
for all ``k`` data vectors or data matrices in `𝐗`.

If `𝐗` holds data matrices, the ``k`` matrices in `𝐗` must have the same
number of columns (i.e., the same number of time series),
but may have any number of (at least `wl`) rows (samples).
All other arguments have
the same meaning as in method (1), with the following differences:

- `⏩`: if true (default), the method run in multi-threaded mode across the
    ``k`` data matrices if ``k`` is at least twice the number of threads Julia
    is instructed to use, otherwise this method attempts to run each spectra
    estimation in multi-threaded mode across series as per method (1).
    See [Threads](@ref).

- If a `planner` is not explicitly passed as an argument,
    the FFTW plan is computed once and applied for all spectral estimations.

**See**: [Spectra](@ref), [plot spectra](@ref).

**See also**: [`crossSpectra`](@ref), [`coherence`](@ref), [goertzel.jl](@ref).

**Examples**:
```julia
using FourierAnalysis

###################################################################

# (1)

# Check that correct amplitude spectra is returned at all discrete
# Fourier Frequency (using a rectangular taper).
# Sinusoids are created at all frequencies with amplitude 10 at the
# first frequency and then incrementing by 10 units along frequencies.
# NB when t is even, correct amplitude for the last frequency is obtained
# only if the sinusoidal has a particular phase.

sr, t, wl= 16, 32, 16
V=Matrix{Float64}(undef, t, wl)
for i=1:wl V[:, i]=sinusoidal(10*i, b2f(i, sr, t), sr, t, π/6) end

# using FFTW.jl only
using FFTW
P=plan_rfft(V, 1)*(2/t);
Σ=abs.(P*V)
using Plots
bar(Σ[brange(t, DC=true), :], labels="")

# using FourierAnalysis.jl
Σ2=spectra(V, sr, t; tapering=rectangular, func=√, DC=true)
using Plots
bar(Σ2.y[brange(t, DC=true), :], labels="")

#############################################################################

### Check amplitude spectra on long signals obtained by welch methods
# one sinusoidal is at an exact discrete Fourier Frequency and the other not
# Rectangular window
sr, t, f, a = 128, 128, 10, 0.5
v=sinusoidal(a, f, sr, t*16)+sinusoidal(a, f*3.5+0.5, sr, t*16)+randn(t*16);
Σ=spectra(v, sr, t; tapering=rectangular, func=√)
bar(Σ.y, labels="rectangular")

# harris4 window (default)
Σ2=spectra(v, sr, t; func=√)
bar!(Σ2.y, labels="harris4")

#smooth spectra
Σ3=smooth(blackmanSmoother, Σ2)
bar!(Σ3.y, labels="smoothed")

#############################################################################

function generateSomeData(sr::Int, t::Int; noise::Real=1.)
    # four sinusoids of length t samples and sr sampling rate
    # peak amplitude: 0.7, 0.6, 0.5, 0.4
    # frequency:        5,   7,  13,  27
    # phase:            0, π/4, π/2,   π
    v1=sinusoidal(0.7, 5,  sr, t, 0)
    v2=sinusoidal(0.6, 7,  sr, t, π/4)
    v3=sinusoidal(0.5, 13, sr, t, π/2)
    v4=sinusoidal(0.4, 27, sr, t, π)
    return hcat(v1, v2, v3, v4) + (randn(t, 4)*noise)
end

sr, wl, t = 128, 512, 8192
X=generateSomeData(sr, t)
# multivariate data matrix 8192x4

# compute spectra
S=spectra(X, sr, wl)

# check the spectrum of first series
S.y[:, 1]

# gather some plot attributes to get nice plots
using Plots.Measures
spectraArgs=(fmax = 32,
             left_margin = 2mm,
             bottom_margin = 2mm,
             xtickfont = font(11, "Times"),
             ytickfont = font(11, "Times"))
plot(S; spectraArgs...)
plot(S; xspace=2, spectraArgs...)

# use a custom simple taperig window
S=spectra(X, sr, wl; tapering=riesz)

# use Slepian's multi-tapering
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.5))

# construct with smoothing
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.5), smoothing=hannSmoother)

# compute Amplitude Spectra instead
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.5), func=√)

# plot Aplitude spectra
plot(S; ytitle="Amplitude", spectraArgs...)

# smooth the spectra a-posteriori
S=smooth(blackmanSmoother, S)

# extract spectra in range (8Hz to 12Hz)
e=extract(S, (8, 12))

# extract spectra in range (8Hz to 12Hz) for series 1 and 2
e=extract(S, (8, 12))[:, 1:2]

# extract the spectra at 10Hz only (1 bin per series)
e=extract(S, 10)

# average spectra in the 8Hz-12Hz range
bar(mean(S, (8.0, 12.0)))

# average across series of the average spectra in the 8Hz-12Hz range
mean(mean(S, (8.0, 12.0)))

# average spectrum across all frequencies for each series
bar(mean(S, :))

# average spectra in equally-spaced 2Hz-band-pass regions for all series
plot(bands(S, 2))

# average spectra in equally-spaced 2Hz-band-pass regions for series 1 and 2
plot(bands(S, 2)[:, 1:2])

# (2)

# generate 3 multivariate data matrices 8192x4
X=[generateSomeData(sr, t) for i=1:3]

# Now the call to the spectra function will generate 3 Spectra objects
S=spectra(X, sr, wl)
plot(S[1]; spectraArgs...)
plot(S[2]; spectraArgs...)
plot(S[3]; spectraArgs...)

# when you want to compute the spectra of many data matrices you may want
# to do it using a fast FFTW plan (wait up to 10s for computing the plan)
plan=Planner(plan_exhaustive, 10.0, wl)
S=spectra(X, sr, wl; planner=plan)

# how faster is this?
using BenchmarkTools
@benchmark(spectra(X, sr, wl))
@benchmark(spectra(X, sr, wl; planner=plan))


# average spectra in range (8Hz-12Hz) for all series of all objects
M=mean(S, (8, 12))

# plot the average spectrum across all series for the three objects
# using Julia standard mean function
plot(mean(S[1].y[:, i] for i=1:size(S[1].y, 2)))
plot!(mean(S[2].y[:, i] for i=1:size(S[2].y, 2)))
plot!(mean(S[3].y[:, i] for i=1:size(S[3].y, 2)))

# average spectra in range (4Hz-32.5Hz) across objects for the 4 series
plot(mean(mean(S, (4, 32.5))))

# extract spectra in range (8Hz-12Hz) for all series and all subjects
extract(S, (8, 12))

# if you enter en illegal range, nothing will be done and you will get
# an error in the REPL explaining what is wrong in your argument
extract(S, (0, 12))
mean(S, (1, 128))

# extract 4Hz-band-pass average spectra for all series and all objects
bands(S, 4)

# Apply smoothing in the spectra computations
S=spectra(X, sr, t; smoothing=blackmanSmoother)
plot(S[1]; spectraArgs...)
plot(S[2]; spectraArgs...)
plot(S[3]; spectraArgs...)

# plot spectra in in 1Hz band-pass regions for all series in S[1]
plot(bands(S[1], 1))

# use slepian multi-tapering
S=spectra(X, sr, wl; tapering=slepians(sr, wl, 1.))
plot(S[1]; spectraArgs...)

# average spectra across objects
plot(mean(s.y for s ∈ S))


```
"""
function spectra( X   :: Union{Vector{T}, Matrix{T}},
                  sr  :: Int,
                  wl  :: Int;
            tapering  :: Union{Taper, TaperKind} = harris4,
            planner   :: Planner                 = getplanner,
            slide     :: Int                     = 0,
            DC        :: Bool                    = false,
            func      :: Function                = identity, # equal to x->x
            smoothing :: Smoother                = noSmoother,
            ⏩       :: Bool                    = true) where T<:Real

    # function to set all parameters
    n, 𝚙, 𝜏, s, e, 𝜎, η, wl½, wl½_, cT = _paramSP!(X, tapering, wl, planner, slide, DC)

    # allocate memory and initialize the spectra :  𝚂 is a (wl½_ x n)-Matrix
    𝚂 = zeros(T, wl½_, n)

    # add to 𝚂 the power spectra of a vector or of the column vectors of matrix `X`,
    # computed on the epoch(s) of duration `t` samples starting at sample `at`.
    # ARGUMENTS:
    # `η` is a normalization Vector to be applied to the FFT coefficients, computed
    #    by function `_fftNormalization`. It is applied element-wise to the FFT coefficient vector(s).
    # `n` (Int) is the number of columns in `X`; it must be 1 if `X` is a vector.
    # `𝜏` is a tapering window (Vector), `rectangular` if no tapering is applied (rectangular taper)
    # `𝚙` is the FFTW plan performing the FFT.
    # `DC` (bool) if true the DC level(s) is(are) computed and stored in the first element (row) of 𝚂
    _spectra!(at, 𝚗) =
        𝜏.kind == rectangular ? ( 𝚂[:, 𝚗] += @. abs2( η*$((𝚙*X[at:at+wl-1, 𝚗])[1+!DC:end]) ) ) :
        ( 𝚂[:, 𝚗] += mean(abs2.( η.*(𝚙*(𝜏.y[:, i].*X[at:at+wl-1, 𝚗]))[1+!DC:end] ) for i=1:𝜏.n) )

    # function cumulating the avarege spectra
    sp!(𝚎, 𝚗) = _spectra!((𝚎*s)+1, 𝚗)

    # average spectra across sliding windows, threaded over `n` if `n` is large
    _thread(⏩, n) ? (for 𝚎=0:e-1 @threads for 𝚗=1:n sp!(𝚎, 𝚗) end end) : (for 𝚗=1:n, 𝚎=0:e-1 sp!(𝚎, 𝚗) end)

    # return spectra as a Vector if input data is a vector, as a matrix otherwise
    return smooth(smoothing, Spectra(func.(X isa Vector ? 𝚂[:] : 𝚂), sr, wl, DC, taperinfo(𝜏), fdf(sr, wl; DC=DC), func, smoothing))
end



function spectra( 𝐗   :: Vector{Matrix{T}},
                  sr  :: Int,
                  wl  :: Int;
            tapering  ::  Union{Taper, TaperKind} = harris4,
            planner   ::  Planner                =  getplanner,
            slide     ::  Int                    =  0,
            DC        ::  Bool                   =  false,
            func      ::  Function               =  identity, # equal to x->x
            smoothing ::  Smoother               =  noSmoother,
            ⏩       ::  Bool                   =  true) where T<:Real

    # set some parameters. k=# input data matrices, n=# time series
    k, n, wl½, thr = length(𝐗), size(𝐗[1], 2), (wl>>1), nthreads()

    # get the FFTW plan (for one time-series) if not provided
    planner ≠ getplanner ? plan=planner : plan = _Planner(plan_rfft(zeros(T, wl)))

    # get the tapering window if not provided
    tapering isa Taper ? 𝜏=tapering : 𝜏=taper(tapering, wl)

    # allocate memory for spectra : 𝗦 is a k-Vector of Spectra objects
    𝗦 = SpectraVector(undef, k)

    # function to compute spectra for the ith input data matrix `𝐗[i]`
    sp!(i, ⏩) = 𝗦[i]=spectra(𝐗[i],
                              sr,
                              wl;
                        tapering  = 𝜏,
                        planner   = plan,
                        slide     = slide,
                        DC        = DC,
                        func      = func,
                        smoothing = smoothing,
                        ⏩       = ⏩)

    # threads over k if k is large, otherwise let spectra function decide if to thread over n
    _thread(⏩, k) ? (@threads for i=1:k sp!(i, false) end) : (for i=1:k sp!(i, ⏩) end)

    return 𝗦
end




# Internal function: create a vector for normalizing element-wise the
# FFT coefficients for spectra estimations.
# With this normalization the `spectra` function
# returns amplitude spectra values equal to the peak-amplitude
# when taking as input a synusoid at any phase and any Fourier discrete
# frequency (with a rectangular tapering window).
# NB when wl is even, correct amplitude for the last frequency is obtained
# only if the sinusoidal has a phase=π/6.
#
# ARGUMENTS:
# `wl` = FFT epoch length (Int)
# `wl½` = t÷2 (Int)
# `wl⁻¹` = inv(t) (Float64)
# `two_wl⁻¹` = 2*inv(t)
# `DC` (Bool); is false the output is of size t÷2, otherwise it is of size t÷2+1
_fftNormalization(wl           :: Int,
                  wl½          :: Int,
                  wl⁻¹         :: Real,
                  two_wl⁻¹     :: Real,
                  DC           :: Bool) =
    DC ? [wl⁻¹; [two_wl⁻¹ for i=1:wl½]] : [two_wl⁻¹ for i=1:wl½]



# Internal function: prepare all parameters that are needed for Welch analysis
# in spectral, cross-spectral and coherence estimates
#
# ARGUMENTS:
# `X` is the input data vector or matrix (n time-series)
# `taper` is a a tapering window ( of enum type `TaperTD`)
# `inrange` is the range (UnitRange{Int}, e.g., 128:1280-1) of the vectors
#      or of rows of `X` on which the analysis is to be performed
# `wl` (Int) is the epoch length for FFT
# `plan` a FFTW.rFFTWPlan forward plan for computing the FFT
# `planflag` flags for the computation of the plan by FFYT. See doc of `spectra`
# `plantime` (Int) is the maximum time allowed (in seconds) to optimize the FFTW plan
# `slide` (Int) is the number of samples the windows slides on (Welch method)
# `⏩` if true multi-threaded computations are requested
#
# RETURN
# `n` (Int) 1 if `X` is a Vector, the number of columns of `X` otherwise
# `𝚙` return the forward FFTW plan. It is computed if argument `plan`=`computeplan`
# `𝜏` (Vector or matrix) the computed tapering window(s)
#       or 0 if no tapering is requested (rectangular taper)
# `n𝜏` the number of tapers. 1 for all tapers (𝜏 is a Vector), but for slepians (𝜏 is a matrix)
# `slide` (Int) checked to be valid
# `e` (Int) the number of sliding windows for the validated `inrange`
# `𝜎=√e` (Real) the squre toot of e, used for in-line computation of the avarage
#       across sliding windows as a normalization factor for the FFT coefficients
# `from` the first sample in the validated `inrange`
# `wl½` `wl`÷2 computed by the right aritmetic bit-wise operation `wl`>>1
# `wl½_`, the length of FFT vectors:  wl½+1 if DC=true, wl½ is DC=false
# `cT` the complex type corresponding to type `T` in function definition (e.g., ComplexF64 if T=Float64)
function _paramSP!( X        :: AbstractArray{T},
                    tapering :: Union{Taper, TaperKind},
                    wl       :: Int,
                    planner  :: Planner,
                    slide    :: Int,
                    DC       :: Bool) where T<:Real

    tₐₗₗ, wl½, wl⁻¹, n = size(X, 1), (wl>>1), inv(wl), size(X, 2)
    wl½_ = wl½ + DC
    tₐₗₗ < wl && @error 📌*", the number of samples in input matrix is smaller than the desired window length."
    planner ≠ getplanner ? 𝚙=planner.p :
        𝚙=plan_rfft(zeros(T, wl), flags=plan_estimate, timelimit=-1.0)
    tapering isa Taper ? 𝜏=tapering : 𝜏=taper(tapering, wl) # tapering window
    if slide <= 0 || slide>tₐₗₗ-wl slide = wl½ end # check the number of sliding samples
    e = (tₐₗₗ-wl)÷slide+1 # number of overlapping epochs
    𝜎 = T(√e)
    η = _fftNormalization(wl, wl½, wl⁻¹/𝜎, T(2/wl)/𝜎, DC)
    cT = typeof(Complex(T(0)))

    return n, 𝚙, 𝜏, slide, e, 𝜎, η, wl½, wl½_, cT
end


# ++++++++++++++++++++  Show override  +++++++++++++++++++ # (REPL output)
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, S::Spectra)
    r=size(S.y, 1)
    vec=S.y isa Vector
    vec ? c=1 : c=size(S.y, 2)
    l=length(S.flabels)
    println(io, titleFont, "▥ Spectra type; $r freq. x $c series")
    #println(io, "□  □    □      □        □           □", defaultFont)
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "sampling rate   (.sr): $(S.sr)")
    println(io, "epoch length    (.wl): $(S.wl)")
    println(io, "DC level        (.DC): $(S.DC)")
    println(io, "taper kind   (.taper): $(S.taper)")
    println(io, "freq. lab. (.flabels): $(l)-", typeof(S.flabels))
    vec ? println(io, "data             (.y): $(r)-", typeof(S.y)) :
          println(io, "data             (.y): $(r)x$(c)-", typeof(S.y))
    # println(io, S.y)
    println(io, "function      (.func): ", string(S.func))
    println(io, "smoother (.smoothing): ", string(S.smoothing))
    r≠l && @warn "number of frequency labels does not match the spectra data size" l r
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, 𝐒::SpectraVector)
    println(io, titleFont, "▥ ⋯ ▥ SpectraVector Type")
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "$(length(𝐒))-element Vector{Spectra}")
end
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
