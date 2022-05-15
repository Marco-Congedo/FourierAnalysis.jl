#   Unit "crossspectra" of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2022,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements Welch cross-spectral estimates using FFTW

#############################################################################
# Generic constructors of CrossSpectra objects:
#
# Enable construction giving only `y`, `sr`, `wl`, `DC`, `taper`
# and `non-linear` argument.
# `flabels` is generated automatically, `smoothing` is set to `noSmoother`
# and `tril` is set to false
CrossSpectra(y, sr, wl, DC, taper, nonlinear) =
        CrossSpectra(y, sr, wl, DC, taper, fdf(sr, wl; DC=DC), nonlinear, noSmoother, false)
#
# As above, but setting by default also `DC` and `nonlinear` to false
CrossSpectra(y, sr, wl, taper) =
    CrossSpectra(y, sr, wl, false, taper, fdf(sr, wl; DC=false), false, noSmoother, false)
#############################################################################

"""
```julia
(1)
function crossSpectra( X    :: Matrix{T},
                       sr   :: Int,
                       wl   :: Int;
                  tapering  :: Union{Taper, TaperKind} = harris4,
                  planner   :: Planner                 = getplanner,
                  slide     :: Int                     = 0,
                  DC        :: Bool                    = false,
                  nonlinear :: Bool                    = false,
                  smoothing :: Smoother                = noSmoother,
                  tril      :: Bool                    = false,
                  ⏩       :: Bool                    = true) where T<:Real

(2)
function crossSpectra( 𝐗    :: Vector{Matrix{T}},
              < same argument sr, ..., ⏩ of method (1) > where T<:Real

```

(1)

Construct a [CrossSpectra](@ref) objects from real multivariate data
using the Welch method.

Given sampling rate `sr` and epoch length `wl`, compute the cross-spectral
matrices of dimension ``n``x``n`` for a multivariate data matrix
`X` of dimension ``t``x``n``, where ``t`` is the number of samples (rows) and
`n>1` is the number of time-series (columns).

The cross-spectral matrices are hold in the `.y` vector field of the created
object. The length of `.y` depends upon
the `wl` argument and `DC` optional keyword argument (see below).

**Optional Keyword Arguments**:

`sr`, `wl`, `tapering`, `planner` and `slide` have the same meaning as for the
[`spectra`](@ref) function.

`DC`: if true the cross-spectral matrix of the DC level is returned
in the first position of `y` (see the fields of the [CrossSpectra](@ref)
object), otherwise (default) the matrices in `y` start with the first positive
discrete frequency, that is, ``sr/wl`` Hz.

`nonlinear`: if true, the amplitude information is eliminated from the
DFT coefficients, that is, they are normalized by their modulus before being
averaged.
This leads to non-linear estimates (Congedo, 2018;
Pascual-Marqui 2007) where the diagonal elements of the cross-spectral matrices
(the spectra) are 1.0 for all frequencies. By default, it is false.

`smoothing`: apply a smoothing function of type [Smoother](@ref)
to the cross-spectral matrices across frequencies.
By default no smoothing is applied.

`tril`: if false (default), the whole cross-spectra matrices will be computed,
otherwise only their lower triangular part (see below).

`⏩`: if true (default), the method is run in multi-threaded mode
across the series in `X` if the number of series is at least twice
the number of threads Julia is instructed to use. See [Threads](@ref).

**Return**

If `tril` is false (default), the output is of type `Array{Hermitian,1}`,
which is the `ℍVector` type used in package [PosDefManifold](https://github.com/Marco-Congedo/PosDefManifold.jl).
Since cross-spectral estimates are Hermitian
positive definite, they can be straightaway used as
argument to PosDefManifold's functions, e.g., for computing matrix moves on
[geodesics](https://marco-congedo.github.io/PosDefManifold.jl/dev/riemannianGeometry/#Geodesic-equations-1),
matrix [distances](https://marco-congedo.github.io/PosDefManifold.jl/dev/riemannianGeometry/#PosDefManifold.distance),
etc. and the the whole vector output to compute
matrix [means](https://marco-congedo.github.io/PosDefManifold.jl/dev/riemannianGeometry/#Means-1),
[spectral embedding](https://marco-congedo.github.io/PosDefManifold.jl/dev/riemannianGeometry/#PosDefManifold.spectralEmbedding)
and more.

If `tril` is true, the output is of type `Array{LowerTriangular,1}`,
which is the `𝕃Vector` type used in PosDefManifold, that is, only the
lower triangle of the cross-spectra is computed in order to save time
and memory.

(2)

Construct a [CrossSpectraVector](@ref) object from a vector of
real multivariate data matrices. Compute the cross-spectral matrices
using the Welch method as per method (1) for all ``k`` data
matrices in `𝐗`.

The ``k`` matrices in `𝐗` must have the same number of columns
(i.e., the same number of time-series), but may have any number of (at least
`wl`) rows (samples).
All other arguments have the same meaning as in method (1),
with the following difference:

`⏩`: if true (default), the method is run in multi-threaded mode across the
``k`` data matrices if ``k`` is at least twice the number of threads Julia
is instructed to use, otherwise this method attempts to run each cross-spectral
estimation in multi-threaded mode across series as per method (1).
See [Threads](@ref).

If a `planner` is not explicitly passed as an argument,
the FFTW plan is computed once and applied for all cross-spectral
estimations.

**Return**

If `tril` is false, the output is of type `Array{Array{Hermitian,1},1}`,
which is the `ℍVector₂` type used in
[PosDefManifold](https://github.com/Marco-Congedo/PosDefManifold.jl).

If `tril` is true, the output is of type `Array{Array{LowerTriangular,1},1}`,
which is the `𝕃Vector₂` type used in PosDefManifold.

**See**: [CrossSpectra](@ref).

**See also**: [`spectra`](@ref), [`coherence`](@ref).

**Examples**:
```julia
using FourierAnalysis, Plots, LinearAlgebra

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

# (1)

X=generateSomeData(sr, t) # multivariate data matrix 8192x4

# cross-spectra using default harris4 tapering window
S=crossSpectra(X, sr, wl)

# check the cross-spectral matrix at frequency 5Hz
S.y[f2b(5, sr, wl)]

# check only the diagonal part of this matrix as a vector
diag(Diagonal(S.y[f2b(5, sr, wl)]))

# cross-spectra using hann tapering window
S=crossSpectra(X, sr, wl; tapering=hann)

# using Slepian's multi-tapering
S=crossSpectra(X, sr, wl; tapering=slepians(sr, wl))

# compute non-linear cross-spectra
S=crossSpectra(X, sr, wl; tapering=slepians(sr, wl), nonlinear=true)

# compute only the lower triangle of cross-spectral matrices
S=crossSpectra(X, sr, wl; tapering=slepians(sr, wl), tril=true)

# smooth a-posteriori the cross-spectra
S2=smooth(blackmanSmoother, S)

# or compute cross-spectra already smoothed
S=crossSpectra(X, sr, wl;
               tapering=slepians(sr, wl), tril=true, smoothing=blackmanSmoother)

# mean cross-spectral matrix in 8Hz-12Hz range
M=mean(S, (8, 12)) # or also M=mean(S, (8.0, 12.0))

# extract all cross-spectral matrices in 8Hz-12Hz range
E=extract(S, (8, 12))

# cross-spectral matrices averaged in 2Hz band-pass regions
B=bands(S, 2)

# Get the spectra from a CrossSpectra object
PowerSpectra=Spectra(S)

# Get the amplitude spectra from a CrossSpectra object
AmpSpectra=Spectra(S, func=√)

# Get the log10-spectra from a CrossSpectra object
log10Spectra=Spectra(S, func=log10)

# plot the spectra (see recipes.jl)
plot(AmpSpectra; fmax=32, xspace=4, ytitle="Amplitude")

# (2)
# generate 3 multivariate data matrices 8192x4
X=[generateSomeData(sr, t) for i=1:3]

# The examples here below use exactly the same syntax as the previous method.
# However, since the input here is a vector of data matrices
# and not a single data matrix, the examples here below create a vector
# of the object created by the examples above.
# For example:

# cross-spectra using the default harris4 tapering window
# this creates a CrossSpectraVector object
S=crossSpectra(X, sr, wl)

# check the cross-spectral matrix at fr. 5Hz for the first CrossSpectra object
S[1].y[f2b(5, sr, wl)]

# check only the diagonal part of this matrix as a vector
diag(Diagonal(S[1].y[f2b(5, sr, wl)]))

# cross-spectra using Hamming's tapering window
S=crossSpectra(X, sr, wl; tapering=hamming)

# using Slepian's multi-tapering
S=crossSpectra(X, sr, wl; tapering=slepians(sr, wl))

# compute non-linear cross-spectra
S=crossSpectra(X, sr, wl; tapering=slepians(sr, wl), nonlinear=true)

# compute only the lower triangle of cross-spectral matrices
S=crossSpectra(X, sr, wl; tapering=slepians(sr, wl), tril=true)

# smooth a-posteriori all CrossSpectra objects in S
S2=smooth(blackmanSmoother, S)

# or compute them all already smoothed
S=crossSpectra(X, sr, wl; tapering=parzen, smoothing=blackmanSmoother)

# mean cross-spectral matrix in 8Hz-12Hz range for all CrossSpectra (CS) objects
M=mean(S, (8, 12)) # or also M=mean(S, (8.0, 12.0))

# grand-average mean of the above across all CS objects
meanM=mean(mean(S, (8, 12)))

# extract all cross-spectral matrices in 8Hz-12Hz range for all CS objects
E=extract(S, (8, 12))

# grand average of cross-spectral matrices in 8Hz-12Hz range for all CS objects
meanE=mean(extract(S, (8, 12)))

# cross-spectral matrices averaged in 2Hz band-pass regions for all CS objects
B=bands(S, 2)

# Get and plot the spectra from a CrossSpectra object
plot(Spectra(S[1]); fmax=32, xspace=4)

# Pre-compute a FFT planner and pass it as argument
# (this interesting if the function is to be called repeatedly).
plan=Planner(plan_exhaustive, 10.0, wl, eltype(X[1])) # wait 10s
S=crossSpectra(X, sr, wl; planner=plan)

# how faster is this?
using BenchmarkTools
@benchmark(crossSpectra(X, sr, wl))
@benchmark(crossSpectra(X, sr, wl; planner=plan))
...
...


```

"""
function crossSpectra( X    :: Matrix{T},
                       sr   :: Int,
                       wl   :: Int;
                  tapering  :: Union{Taper, TaperKind} = harris4,
                  planner   :: Planner                 = getplanner,
                  slide     :: Int                     = 0,
                  DC        :: Bool                    = false,
                  nonlinear :: Bool                    = false,
                  smoothing :: Smoother                = noSmoother,
                  tril      :: Bool                    = false,
                  ⏩       :: Bool                    = true) where T<:Real

    # function to set all parameters
    n, 𝚙, 𝜏, s, e, 𝜎, η, wl½, wl½_, cT = _paramSP!(X, tapering, wl, planner, slide, DC)

    # allocate memory for DFTs (𝙵)
    𝙵 = [zeros(cT, wl½_, n) for i=1:𝜏.n]

    # allocate memory for cross-spectra and initialize them with zeros (𝙎)
    𝙎 = Vector{LowerTriangular}([LowerTriangular(zeros(cT, n, n)) for i=1:wl½_])

    # compute the descrete Fourier transform (DFT) of a Vector or of the column vectors
    # of matrix `X` and put it into the columns of `𝙵`. The DFT is computed on the epoch
    # of duration `t` samples starting at sample `at`. It does not return anything.
    # `𝙵` is a `n𝜏`-Vector of DFT matrices, holding a DFT matrix (# discr. freq. x time-series),
    # where `n𝜏` is the number of tapering window (>1 only for Slepian tapers).
    # ARGUMENTS:
    # `η` is a calibration Vector to be applied to the FFT coefficients, computed
    #    by function `_fftNormalization`. It is applied element-wise to the FFT
    #    coefficient vector(s) to obtain spectral estimates in correct physical units.
    #    It is applied only if `nonlinear` is false.
    # `n` (Int) is the number of columns in `X`; it must be 1 if `X` is a vector.
    # `𝜏` is a tapering window, ``rectangular` if no tapering is applied
    # `𝚙` is the FFTW plan performing the FFT.
    # `DC` (bool): if true the DC level(s) is(are) computed and stored in the first matrix of 𝚂
    _dft!(at, 𝚗) =
        if nonlinear
            𝜏.kind == rectangular ? ( 𝙵[1][:, 𝚗] = ((𝚙*X[at:at+wl-1, 𝚗])[1+!DC:end]) ) : # no tapering
                                    ( for i=1:𝜏.n 𝙵[i][:, 𝚗] = (𝚙*(𝜏.y[:, i].*X[at:at+wl-1, 𝚗]))[1+!DC:end] end ) # with tapering
        else
            𝜏.kind == rectangular ? ( 𝙵[1][:, 𝚗] = η.*((𝚙*X[at:at+wl-1, 𝚗])[1+!DC:end]) ) : # no tapering
                                    ( for i=1:𝜏.n 𝙵[i][:, 𝚗] = η.*(𝚙*(𝜏.y[:, i].*X[at:at+wl-1, 𝚗]))[1+!DC:end] end ) # with tapering
        end

    # function computing DFTs (discrete Fourier transform)
    dft!(𝚎, 𝚗) = _dft!((𝚎*s)+1, 𝚗)

    # function cumulating the avarage cross-spectra. The mean is over tapering windows (>1 only for Slepians)
    cs!(𝚏) = for j=1:n, i=j:n 𝙎[𝚏][i, j] += mean(𝙵[l][𝚏, i] * conj(𝙵[l][𝚏, j]) for l=1:𝜏.n) end

    # average cross-spectra across sliding windows
    for 𝚎=0:e-1
        _thread(⏩, n) ? (@threads for 𝚗=1:n dft!(𝚎, 𝚗) end) : (for 𝚗=1:n dft!(𝚎, 𝚗) end)
        nonlinear ? (for i=1:length(𝙵) 𝙵[i]./=(abs.(𝙵[i])*𝜎) end) : nothing
        _thread(⏩, wl½_) ? (@threads for 𝚏=1:wl½_ cs!(𝚏) end) : (for 𝚏=1:wl½_ cs!(𝚏) end)
    end

    # if `tril`=false return the full matrices, otherwise the lower-triangular matrices
    # Typecasting output: 𝔸=Array; 𝕄=Matrix; ℍ=Hermitian; 𝕃=LowerTriangular
    H=Hermitian; V=Vector; M=Matrix; CS=CrossSpectra; sm=smoothing
    args=(sr, wl, DC, taperinfo(𝜏), fdf(sr, wl; DC=DC), nonlinear, sm, tril)
    return tril ? smooth(sm, CS(𝙎, args...)) :
                  smooth(sm, CS(V{H}([H(M(𝙎[𝚏]), :L) for 𝚏=1:wl½_]), args...))
end



function crossSpectra( 𝐗    :: Vector{Matrix{T}},
                       sr   :: Int,
                       wl   :: Int;
                  tapering  :: Union{Taper, TaperKind} = harris4,
                  planner   :: Planner                 = getplanner,
                  slide     :: Int                     = 0,
                  DC        :: Bool                    = false,
                  nonlinear :: Bool                    = false,
                  smoothing :: Smoother                = noSmoother,
                  tril      :: Bool                    = false,
                  ⏩       :: Bool                    = true) where T<:Real

    # set some parameters. k=# input data matrices, n=# time series
    k, n, wl½, thr, cT = length(𝐗), size(𝐗[1], 2), (wl>>1), nthreads(), typeof(Complex(T(0)))

    # get the FFTW plan (for one time-series)
    planner ≠ getplanner ? plan=planner : plan = _Planner(plan_rfft(zeros(T, wl)))

    # get the tapering window if not provided
    tapering isa Taper ? 𝜏=tapering : 𝜏=taper(tapering, wl)

    # allocate memory for cross-spectra
    𝓢 = CrossSpectraVector(undef, k)

    # define internal function to compute cross-spectra for the ith input data matrix `𝐗[i]`
    cs!(i, ⏩) = 𝓢[i]=crossSpectra(𝐗[i],
                                   sr,
                                   wl;
                               tapering  = 𝜏,
                               slide     = slide,
                               DC        = DC,
                               planner   = plan,
                               nonlinear = nonlinear,
                               smoothing = smoothing,
                               tril      = tril,
                               ⏩       = ⏩)

    # threads over k if k is large, otherwise let crossSpectra function decide if to thread over n
    _thread(⏩, k) ? (@threads for i=1:k cs!(i, false) end) : (for i=1:k cs!(i, ⏩) end)

    return 𝓢
end


# ++++++++++++++++++++  Show override  +++++++++++++++++++ # (REPL output)
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, S::CrossSpectra)
    r=size(S.y[1], 1)
    l=length(S.flabels)
    ll=length(S.y)
    println(io, titleFont, "▦ ⋯ ▦ CrossSpectra type; $l freq. x $(r)²-matrices")
    #println(io, "□  □    □      □        □           □", defaultFont)
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "sampling rate    (.sr): $(S.sr)")
    println(io, "epoch length     (.wl): $(S.wl)")
    println(io, "DC level         (.DC): $(S.DC)")
    println(io, "taper kind (.tapering): $(S.taper)")
    println(io, "freq. lab.  (.flabels): $(l)-", typeof(S.flabels))
    println(io, "data              (.y): $(l)-", typeof(S.y))
    # println(io, S.y)
    println(io, "non-linear(.nonlinear): ", string(S.nonlinear))
    println(io, "smoother  (.smoothing): ", string(S.smoothing))
    println(io, "triang. form   (.tril): ", string(S.tril))
    ll≠l && @warn "number of frequency labels does not match the CrossSpectra data size" l ll
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, 𝐒::CrossSpectraVector)
    println(io, titleFont, "▦ ⋯ ▦ CrossSpectraVector Type")
    println(io, titleFont, "  ⋯")
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "$(length(𝐒))-element Vector{CrossSpectra}")
end

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
