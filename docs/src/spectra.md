# spectra.jl

Spectra objects created by *FourierAnalysis* are incapsulated in
the following structure:

## Spectra

**Categories**: [data objects](@ref), [FDobjects](@ref)

```
struct Spectra
    y           :: AbstractArray{T} where T<:Real
    sr          :: Int
    wl          :: Int
    DC          :: Bool
    taper       :: String
    flabels     :: Vector{T} where T<:Union{Real, Int}
    func        :: Function
    smoothing   :: Smoother
end
```

**Fields**:

`y`: a real vector holding one *spectrum*, or a real matrix holding several *spectra* in its columns.

`sr`: the *sampling rate* of the data on which the spectra have been estimated.

`wl`: the FFT *window length* used for estimating the spectra.

`DC`: if true, the first row of `y` holds the *DC level*, otherwise it holds the first positive frequency. Thus, if `DC` is false, the first dimension of `y` is equal to wl÷2 (integer division), otherwise it is equal to
(wl÷2)+1 (see [Overview](@ref)). In all constructors it is false
by default.

`taper`: the time-domain *tapering window* used for FFT computation, as a string, with parameters in parentheses for Slepian's dpss. See [tapers.jl](@ref).

`flabels`: a vector holding all Fourier discrete frequencies in Hz.
Those are the *frequency labels* for the rows of `y`. If `DC` is true,
the first label is ``0``, otherwise it is the first positive frequency,
which is equal to the frequency resolution sr/wl.

`func`: a *function* applied element-wise to the spectra.
In all constructors the default is the `identity` (do-nothing) function.

`smoothing`: a [Smoother::Enumerated Type](@ref). A flag indicating
whether the spectra have been smoothed across adjacent
frequencies. If no smoothing has been applied, it is equal to `noSmoother`,
which is the default in all constructors.

**Note**: In Julia the fields are accessed by the usual dot notation, e.g.,
you may verify that for *Spectra* object S, length(S.flabels) == size(S.y, 1) == (S.wl/2)+S.DC.

A vector of *Spectra* objects is of type [SpectraVector](@ref).

**Methods for Spectra and SpectraVector objects**

|      method          |   Spectra   | SpectraVector |
|:---------------------|:-----------:|:-------------:|
| [`bands`](@ref)      |     ✔      |      ✔      |
| [`extract`](@ref)    |     ✔      |      ✔      |
| [`mean`](@ref)       |     ✔      |      ✔      |
| [`smooth`](@ref)     |     ✔      |      ✔      |
| [`sameParams`](@ref) |            |      ✔      |


**Generic Constructors**:

In order to construct a *Spectra* object from univariate and multivariate
data, *FourierAnalysis* provides two [`spectra`](@ref) constuctors, which
is what you will use in practice most of the time.

Manual constructors are also possible, for which you have to provide
appropriate arguments. The default manual constructor of *Spectra* objects is

```
Spectra(y, sr, wl, DC, taper, flabels, func, smoothing).
```

Other constructors are also provided:

```
Spectra(y, sr, wl, DC, tapers)
```

generate the appropriate `flabels`, set `func` to `identity`
(do-nothing) and `smoothing` to `noSmoother`;

```
Spectra(y, sr, wl, taper)
```

is like the constructor here above, but also set `DC` to false;

```
Spectra(𝙎::CrossSpectra;
        func::Function=identity)
```

create a *Spectra* object extracting the spectra from a [CrossSpectra](@ref)
object. If a function is provided with the `func` argument,
this function is applied element-wise to the spectra.
For instance,
- `func=sqrt` will extract amplitude spectra,
- `func=log` will extract log-spectra,
- `func=dB` will extract spectra in deciBels (see [`dB`](@ref)).
By default the `identity` (do-nothing) function
is applied, thus (power) spectra are extracted;

!!! note "Nota Bene"
    If the CrossSpectra object is non-linear, the spectra are
    uniformly equal to 1.0. See [`crossSpectra`](@ref).

```
Spectra(𝙎::CrossSpectraVector;
        func::Function=identity)
```

create a [SpectraVector](@ref) object from a [CrossSpectraVector](@ref) object,
calling the constructor here above for all [`crossSpectra`](@ref) objects
hold by `𝙎`.

**Constructors from data**:

```@docs
spectra
```
