# tapers.jl

This unit implements eight tapering windows and the
discrete prolate spheroidal sequences, the latters via the
[DSP](https://github.com/JuliaDSP/DSP.jl) package.


## TaperKind
```
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
```

The Hann tapering window is also known as 'squared cosine' window and
the Riesz window is similar to window known as the 'cosine' window.

The design of tapering windows implies a trade-off between
the *equivalent noise bandwidth* (enb),
the energy of the *first sidelobe* (fsl) and the rate of
*sidelobe falloff* (slf). The characteristics of the implemented tapering windows are reported in the following table:

|    window   | notable points   |enb(bins)|fsl(dB)| slf(dB/octave) |
|:-----------:|:----------------:|:-------:|:-----:|:--------------:|
| rectangular | 1 everywhere     |1        |  -13  |    - 6         |
| triangular  | 0 at boundaries  |1.33     |  -26  |    -12         |
| hann        | 0 at boundaries  |1.50     |  -32  |    -18         |
| hamming     | >0 at boundaries |1.36     |  -43  |    - 6         |
| blackman    | 0 at boundaries  |1.73     |  -58  |    -18         |
| harris4     | >0 at boundaries |2        |  -92  |    - 6         |
| riesz       | 0 at boundaries  |1.2      |  -21  |    -12         |
| parzen      | >0 at boundaries |1.92     |  -53  |    -24         |

The harris4 tapering window features excellent first sidelobe (-92dB) and sidelob falloff (-6dB rate),
at the expenses of the highest equivalent noise bandwidth among all.
Since this latter parameter is not critical in many applications,
this window is employed as default by all *FourierAnalysis* constructors of objects in the frequency domain.

For reducing the variance of the spectral estimations, use the (Slepian) *discrete prolate spheroidal sequences (dpss)* multi-tapering (see [`slepians`](@ref)).

## Taper

Tapering windows in **FourirAnalysis** are encapsulated in the following structure:

```
struct Taper
    y    :: Union{Vector{T}, Matrix{T}} where T<:Real
    kind :: TaperKind
    α    :: Real
    n    :: Int
end
```

The fields of the structure are:

- `y`, a real vector holding the tapering window, but for Slepian multi-tapers, for which this is a matrix holding in its columns the dpss
- `kind`, the tapering window(s) as a [TaperKind](@ref)
- `α`, a parameter for the tapering window(s). This is needed only for dpss
- `n`, the number of tapering windows. It is >1 only for dpss.

If you need to construct *Taper* objects for single tapering windows, use the universal [`taper`](@ref) constructor.
For constructing *dpss* use the specialized constructor [`slepians`](@ref).

```@docs
taper
slepians
taperinfo
```


**References**

F.J. Harris (1978)
[On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform](http://www.utdallas.edu/~cpb021000/EE%204361/Great%20DSP%20Papers/Harris%20on%20Windows.pdf)
Proc. IEEE, 66, 51-53, 1978

D. Slepian (1978)
[Prolate Spheroidal Wave Functions. Fourier Analysis, and Uncertainty—V: The Discrete Case](https://ieeexplore-ieee-org.gaelnomade-1.grenet.fr/stamp/stamp.jsp?tp=&arnumber=6771595)
The Bell System Technical Journal,VoL 57, No. 5. May-June 1978

D.J. Thomson (1982)
[Spectrum estimation and harmonic analysis](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.471.1278&rep=rep1&type=pdf)
Proc. IEEE 70: 1055-1096, 1982.

[More resources](https://pdfs.semanticscholar.org/752d/1a551b96559458064323eb3de7faaaef4c4e.pdf)
