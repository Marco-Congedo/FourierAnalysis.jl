# crossspectra.jl

Cross-spectra objects created by *FourierAnalysis* are incapsulated in
the following structure:

## CrossSpectra

**Categories**: [data objects](@ref), [FDobjects](@ref).

```julia
struct CrossSpectra
    y           :: Union{Vector{LowerTriangular}, Vector{Hermitian}}
    sr          :: Int
    wl          :: Int
    DC          :: Bool
    taper       :: String
    flabels     :: Vector{T} where T<:Union{Real, Int}
    nonlinear   :: Bool
    smoothing   :: Smoother
    tril        :: Bool
end
```

**Fields**:

`y`: a real vector of `LowerTriangular` or `Hermitian` matrices
(see [Linear Algebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
Julia standard package), depending on the `tril` field, as explained below,
holding the *cross-spectral* matrices.

`sr`: the *sampling rate* of the data on which the cross-spectra have been
estimated.

`wl`: the FFT *window length* used for estimating the cross-spectra.

`DC`: if true, the first matrix holds the cross-spectra of the *DC level*,
otherwise it holds the cross-spectra of the first positive frequency.
Thus, if `DC` is false, the number of matrices in `y` is equal to ``wl÷2``
(integer division), otherwise it is equal to ``(wl÷2)+1`` (see [Overview](@ref)).
In all constructors it is false by default.

`taper`: the time-domain *tapering window* used for FFT computation,
as a string, with parameters in parentheses for Slepian's dpss.
See [tapers.jl](@ref).

`flabels`: a vector holding all Fourier discrete frequencies in Hz.
Those are the *frequency labels* for the matrices in `y`. If `DC` is true,
the first label is ``0``, otherwise it is the first positive frequency,
which is equal to the frequency resolution ``sr/wl``.

`nonlinear`: if true, the amplitude information has been eliminated from the
DFT coefficients, that is, the coefficients have been normalized by their
modulus before being averaged to obtain the cross-spectra.
This leads to non-linear estimates (Congedo, 2018;
Pascual-Marqui 2007) where the diagonal elements of the cross-spectral matrices
(the spectra) are 1.0 for all frequencies. In all constructors it is false
by default.

`smoothing`: a [Smoother](@ref) flag indicating
whether the cross-spectral matrices have been smoothed across adjacent
frequencies. If no smoothing has been applied, it is equal to `noSmoother`,
which is the default in all constructors.

`tril`: if false, the cross-spectral matrices in `y` are full `Hermitian` matrices,
otherwise they are `LowerTriangular` matrices holding only the lower triangles
of the cross-spectra. In all constructors it is false by default.

**Note**: In Julia the fields are accessed by the usual dot notation, e.g.,
you may verify that for *CrossSpectra* object `S`,
```length(S.flabels) == length(S.y)== (S.wl/2)+S.DC```.

A vector of *CrossSpectra* objects is of type [CrossSpectraVector](@ref).

**Methods for CrossSpectra and CrossSpectraVector objects**

|      method          | CrossSpectra | CrossSpectraVector |
|:---------------------|:------------:|:------------------:|
| [`bands`](@ref)      |     ✔        |         ✔         |
| [`extract`](@ref)    |     ✔        |         ✔         |
| [`mean`](@ref)       |     ✔        |         ✔         |
| [`smooth`](@ref)     |     ✔        |         ✔         |
| [`sameParams`](@ref) |              |         ✔         |
| [`isLinear`](@ref)   |              |         ✔         |
| [`isNonLinear`](@ref)|              |         ✔         |


**Generic Constructors**:

In order to construct a *CrossSpectra* object from multivariate
data using the Welch method, *FourierAnalysis*
provides two [`crossSpectra`](@ref) constuctors, which
is what you will use in practice most of the time.

Manual constructors are also possible, for which you have to provide
appropriate arguments. The default manual constructor of *CrossSpectra*
objects is

```julia
CrossSpectra(y, sr, wl, DC, taper, flabels, nonlinear, smoothing, tril)
```

Other constructors are also provided:

```julia
CrossSpectra(y, sr, wl, DC, taper, nonlinear)
```
enable construction giving only `y`, `sr`, `wl`, `DC`, `taper`
and `nonlinear` argument.
`flabels` is generated automatically, `smoothing` is set to `noSmoother`
and `tril` is set to false;

```julia
CrossSpectra(y, sr, wl, taper)
```

as above, but setting by default also `DC` and `nonlinear` to false.

**Constructors from data**:

```@docs
crossSpectra
```

**References**:

M. Congedo (2018),
[Non-Parametric Synchronization Measures used in EEG and MEG](https://hal.archives-ouvertes.fr/hal-01868538v2/document),
Technical Report. GIPSA-lab, CNRS, University Grenoble Alpes, Grenoble INP.

R. Pascual-Marqui (2007),
[Instantaneous and lagged measurements of linear and nonlinear dependence between groups of multivariate time series: frequency decomposition](https://arxiv.org/ftp/arxiv/papers/0711/0711.1455.pdf),
arXiv:0711.1455.
