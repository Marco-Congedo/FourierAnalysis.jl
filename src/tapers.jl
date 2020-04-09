#   Unit "tapers" of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2020, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements several tapering windows in the time domain.
#
#   See F.J. Harris
#   "On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform"
#   Proc. IEEE, 66, 51-53, 1978
#   http://www.utdallas.edu/~cpb021000/EE%204361/Great%20DSP%20Papers/Harris%20on%20Windows.pdf
#
#   D. Slepian
#   "Prolate Spheroidal Wave Functions. Fourier Analysis, and Uncertainty—V: The Discrete Case"
#   The Bell System Technical Journal,VoL 57, No. 5. May-June 1978
#   https://ieeexplore-ieee-org.gaelnomade-1.grenet.fr/stamp/stamp.jsp?tp=&arnumber=6771595
#
#   D.J. Thomson
#   "Spectrum estimation and harmonic analysis."
#   Proc. IEEE 70: 1055-1096, 1982.
#   http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.471.1278&rep=rep1&type=pdf

#   Other resources :
#   https://pdfs.semanticscholar.org/752d/1a551b96559458064323eb3de7faaaef4c4e.pdf
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


# time-domain tapering windows.
@enum TaperKind begin
    rectangular = 1
    triangular  = 2
    hann        = 3
    hamming     = 4
    blackman    = 5
    harris4     = 6
    riesz       = 7
    parzen      = 8
    slepian     = 9
end

TaperData=Union{Vector{T}, Matrix{T}} where T<:Union{Real, Complex}

struct Taper
    y    :: TaperData # here is the taper
    kind :: TaperKind # enumerated type; see here above
    α    :: Real # 1-parameter tapers are supported
    n    :: Int # number of tapers; >1 only for slepians
end


### Time-domain tapering ###
###############################################################################
"""
```
function taper( kind  :: TaperKind,
                wl    :: Int;
            α       :: Real    = 2.,
            n       :: Int     = ceil(Int, 2*α)-1,
            padding :: Int     = 0,
            type    :: Type{T} = Float64) where T<:Union{Real, Complex}

```
Universal constructor of [Taper](@ref) objects, given a tapering window
`kind`, of type [TaperKind](@ref) and the window length `wl`.

Return a vector of length `wl` for all types of tapers, but for the dpss
(Slepian multi-tapers), for which return a matrix of size `wl` x `n`.

If optional keyword argument `padding` is >0, then the actual window length
will be `wl` + `padding`.

The `type` optional keyword argument can be use to specify the the type of
elements of the typering window. By default, this is the `Float64` type.

Optional keywords arguments `α` and `n` currently apply only for slepian
multi-tapering (*discrete prolate spheroidal sequences*):

`α` is the *half bandwidth* (hbw) parameter as per the DSP.jl package.
This unit is used in many dpss implementations and it is often reported
that "typical values" for `α` are 2, 2.5, 3, 3.5 or 4.
However, the optimal smoothing increases with the ratio between window length
and sampling rate, thus these values in absolute terms are not useful.
In fact, the larger the *hbw* and the higher `n`, the smoother the
spectra will be (variance reduction). In order to overcome these
difficulties, for Slepian's multitapering
*FourirAnalysis* implements the [`slepians`](@ref) constructor,
which allows the bandwidth parameter to be given in Hz and the
number of tapering windows to be chosen automatically.

`n` is the number of tapering windows. For slepian tapers this is
the number of the discrete prolate spheroidal sequences.
As in the DSP package, by default it is set to
`ceil(Int, 2*α)-1`, however, depending on how large this number is,
low eigenvalues may correspond to the last sequences, therefore
those should be discarded.

**See**: [plot tapering windows](@ref).

**See also**: [`slepians`](@ref), [`taperinfo`](@ref)

**Examples**:
```
using FourierAnalysis

## Use the constructor
sr, t, f, a = 128, 128, 10, 0.5
# create a sinusoidal superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# create a data matrix
X=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
# compute spectra using hamming tapering window
# we need to prepend 'FourierAnalysis.' since `hamming`
# is declared also in DSP.jl
H=taper(FourierAnalysis.hamming, t)
S=spectra(X, sr, t; tapering=H)
# you can obtain the same thing with
S=spectra(X, sr, t; tapering=H)
# which will create the hamming tapering window on the fly,
# thus calling explicitly the constructor is interesting
# only if you need to reuse the same tapering window many times.

## Plot tapering windows using the standard plot function
using Plots
tapers=[TaperKind(i) for i=1:8]
X=zeros(t, 8)
for i=1:8 X[:, i] = taper(tapers[i], t).y end
mylabels=Array{String}(undef, 1, 8)
for i=1:8 mylabels[1, i]=string(tapers[i]) end
plot(X; labels=mylabels)

## using the recipe declared in recipes.jl
plot(taper(parzen, 256))
plot(taper(slepian, 256, α=4, n=7))

```
"""
function taper( kind :: TaperKind,
                wl   :: Int;  # wl is the window length
            α        :: Real    = 2.,    # only for splepians;
            n        :: Int     = ceil(Int, 2*α)-1, #ceil : Int greater then or equal to
            padding  :: Int     = 0,
            type     :: Type{T} = Float64) where T<:Union{Real, Complex}

    if      kind ∉ (rectangular, slepian)
            wl₋₁     = wl-1
            wl₋₁½    = wl₋₁/2
            c(i::Int) = cos(2π*(i-1)/wl₋₁)
            c(i::Int, j::Int) = cos(2π*j*(i-1)/wl₋₁)
    end

    if      kind == rectangular v = ones(type, wl)
    elseif  kind == triangular  v = Vector{type}([1 - abs(2*(i-1)-wl₋₁)/wl₋₁ for i=1:wl]) # zero at the edge
    elseif  kind == hann        v = Vector{type}([0.5*(1-c(i)) for i=1:wl])
    elseif  kind == hamming     v = Vector{type}([0.54 - 0.46*c(i) for i=1:wl])
    elseif  kind == blackman    v = Vector{type}([0.42 + 0.50*-c(i) + 0.08*c(i, 2) for i=1:wl])
                                     v[1]=0.; v[end]=0.
    elseif  kind == harris4     v = Vector{type}([0.35875 - 0.48829*c(i) + 0.14128*c(i, 2) - 0.01168*c(i, 3) for i=1:wl])
    elseif  kind == riesz       v = Vector{type}([1 - ((i-1-wl₋₁½)/wl₋₁½)^2 for i=1:wl])
    elseif  kind == parzen
            v=Vector{type}(undef, wl)
            wl½, wl¼= wl÷2, wl÷4
            for i=1:wl¼        v[i] = 2*((1 - (abs(i-wl½)/wl½))^3) end
            for i=wl¼+1:wl½    v[i] = 1 - 6*(abs(i-wl½)/wl½ )^2 + 6*(abs(i-wl½)/wl½)^3 end
            for i=wl½+1:wl-wl¼ v[i] = 1 - 6*(abs(i-1-wl½)/wl½)^2 + 6*(abs(i-1-wl½)/wl½)^3 end
            for i=wl-wl¼+1:wl  v[i] = 2*((1 - (abs(i-1-wl½)/wl½))^3) end
    elseif  kind == slepian    v    = Matrix{type}(dpss(wl, α, n, zerophase=false))
    end

    # This has to be checked:

    # normalize to unit mean or unit mean of absolute values
    if      kind ∉ (rectangular, slepian) v ./= mean(v) end
    if      kind == slepian v./=mean(abs.(a) for a in v) end

    # if      kind == slepian v./=mean(v.^2; dims=1) end
    # if      kind ∉ (rectangular, slepian) v ./= mean(v) end

    # normalization in DSP.jl
    #    if      taper == slepian for i=1:size(v, 2) v[:, i]./= (norm(v[:, i])/wl) end end
    # as above but givin less and less weight to the sequences
    #    if      taper == slepian for i=1:size(v, 2) v[:, i]./= √i*(norm(v[:, i])/wl) end end

    p=padding
    p>0 && (kind == slepian ? v=[v; zeros(type, p, n)] : v=[v; zeros(type, p)])

    return kind == slepian ? Taper(v, kind, α, n) : Taper(v, kind, 0., 1)
end

"""
Construct a [Taper](@ref) objects holding Slepian's multi-tapering
*discrete prolate spheroidal sequences*,
given sampling rate `sr`, window length `wl` and
the `bandwidth` argument in Hz.
For EEG data, 1<=bandwidth<=2 is an adequate choice.

The 'half-bandwidth' parameter `α` used in the DSP package and in the
universal [Taper](@ref) constructor is set as

        `α=(bandwidth/2)*wl/sr`.

The optimal number of dpss is heuristically set to

        `n=max(1, trunc(Int, 2*α)-trunc(Int, log(2*α)))`.

The created object can be passed as argument
in constructors [`spectra`](@ref), [`crossSpectra`](@ref) and
[`coherence`](@ref).

**See**: [plot tapering windows](@ref).

**Examples**:
```
using FourierAnalysis
sr, t, f, a = 128, 128, 10, 0.5
# create a sinusoidal superimposed to white noise
v=sinusoidal(a, f, sr, t*16, 0) + randn(t*16)
# create a data matrix
X=broadcast(+, v, randn(t*16, 3))*randn(3, 3)
# compute spectra using slepian multi-tapering with bandwidth 1.5
H=slepians(sr, t, 2)
S=spectra(X, sr, t; tapering=H)

using Plots
plot(H)
plot(S)
```
"""
function slepians( sr    :: Int,
                   wl    :: Int,
               bandwidth :: Real = 1.5)
    α=(bandwidth/2)*wl/sr # α parameter
    # heuristic to eliminate eigenfunctions with small eigenvalues
    n=max(1, trunc(Int, 2*α)-trunc(Int, log(2*α)))
    return taper(slepian, wl, α=α, n=n)
end
###############################################################################

"""
    function taperinfo(taper::Taper)

Return the name of the
tapering window(s) encapsulated in the [Taper](@ref) object
as a string.

Only for Slepian's discrete prolate spheroidal sequences (dpss),
their parameters, namely, ``α`` (half-bandwidth)
and ``n`` (number of windows), are reported within parentheses as well.

**Examples**:
```
H=taper(hamming, 128*8)
taperinfo(H)

H=slepians(128, 128*8, 2)
taperinfo(H)
```
"""
taperinfo(taper::Taper) =
     taper.kind==slepian ? string(taper.kind)*"'s dpss (alpha=$(taper.α), n=$(taper.n))" :
                           string(taper.kind)


# ++++++++++++++++++++  Show override  +++++++++++++++++++ # (REPL output)
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, 𝜏::Taper)
println(io, titleFont, "⍓ Taper type; $(size(𝜏.y, 1))-samples")
#println(io, "□  □    □      □        □           □", defaultFont)
println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
println(io, "taper kind   (.kind): ", string(𝜏.kind))
println(io, "half-bandwidth  (.α): $(𝜏.α)")
println(io, "number of tapers(.n): $(𝜏.n)")
𝜏.kind==slepian ? println(io, "data            (.y): $(size(𝜏.y, 1))x$(size(𝜏.y, 2))-", typeof(𝜏.y)) :
                 println(io, "data            (.y): $(length(𝜏.y))-", typeof(𝜏.y))

end
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
