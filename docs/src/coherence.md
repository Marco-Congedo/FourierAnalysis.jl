# coherence.jl

*FourierAnalysis* estimates coherence both in
the *frequency domain*, in this unit,
and in the *time-frequency domain*,
in unit [timefrequencybi.jl](@ref).

### Notation

The following notation will be used:

|        symbol         |            meaning             | Julia function |
|:---------------------:|:------------------------------:|:--------------:|
| superscript ``^*``    |    complex conjugate           |     conj()     |
| superscript ``^H``    | complex conjugate transpose    |       '        |
| ``\mid \cdot \mid``   |           modulus              |     abs()      |
| ``\Bbb R(\cdot)``     |  real part of the argument     |     real()     |
| ``\Bbb C(\cdot)``     | imaginary part of the argument |     imag()     |
| ``\left<\cdot\right>``| an average across argument's elements | mean()|

## coherence estimates

**frequency domain (FD)**

In the FD, the input is a continuous multivariate data matrix
``X=[x_1 ... x_n]``,
where each column is a series (e.g., a time-series).
Coherence estimates take the form of a ``n``x``n`` coherence matrix ``C(f)``
for each discrete Fourier frequency ``f``. Those are
normalized versions of the ``n``x``n`` *cross-spectral matrices* ``S(f)``
at corresponding frequencies.

Cross-spectral matrices are estimated
via a [Welch procedure](https://en.wikipedia.org/wiki/Welch%27s_method),
that is, they are averaged across windows sliding
over ``X``. The windows may be overlapping or not. Let ``Y_l(f)`` be the discrete Fourier transform at discrete
Fourier frequency ``f`` of the ``l^{th}`` window of ``X``.
The cross-spectral matrix estimation at frequency ``f`` is then given by

``S(f)=\left<Y(f)Y(f)^H\right>``,

where the average is across a number of sliding windows.

Therefore, coherence matrices in the FD is a measure of the
*synchronization between all ``n`` series of ``X`` taken pair-wise*,
resulting in a symmetric matrix for each discrete Fourier frequency.

**time-frequency domain (TFD)**

There exist many time-frequency representation, e.g., those obtained
using [wavelets](https://en.wikipedia.org/wiki/Wavelet),
[short-time Fourier transform](https://en.wikipedia.org/wiki/Short-time_Fourier_transform), etc.
In *FourierAnalysis* time-frequency representations are obtained passing
a univariate signal through a
[filter bank](https://en.wikipedia.org/wiki/Filter_bank)
and computing the
[analytic signal](https://en.wikipedia.org/wiki/Analytic_signal)
of the filter output via the
[Hilbert transform](https://en.wikipedia.org/wiki/Hilbert_transform).
This results in a matrix with the band-pass regions of the filer bank in rows and the ``t`` time samples in columns (see [`filterbank`](@ref) and
[`analyticsignal`](@ref) for details).

For estimating coherence in the TFD, the input is composed of two sets of ``k``
univariate data vectors
``𝐱_a=\{x_{a1}, x_{a2},...,x_{ak}\}`` and ``𝐱_b=\{x_{b1}, x_{b2},...,x_{bk}\}``,
where each element of the set is to be understood as a single window holding
``t`` samples and the element in corresponding positions forms pairs.
Let ``Y_r(a)`` and ``Y_r(b)`` be the time-frequency analytic signal
of realizations ``x_{ar}`` and ``x_{ar}``, respectively, with ``r=1...k``.
Coherence estimates also have a time-frequency representation.
They are a normalized version of crossed analytic signals, which
is analogous to cross-spectra,and in a time-frequency plane is given by

``Z=\left<Y(a) \odot Y(b)^*\right>``,

where the average is across the ``k`` realizations.

Therefore, the coherence matrix in the TFD is a measure of the
*synchronization between corresponding elements of the set ``𝐱_a`` and ``𝐱_b``* for each point in the time-frequency plane.

### kinds of coherence

*FourierAnalysis* estimates several kinds of *linear* and *non-linear*
squared coherences, simply referred to as *coherence*.
For simplicity of notation, let ``i`` and ``j`` indicate two series
in the *frequency domain* and two sets of series in the
*time-frequency domain*, as explained in the previous section. Let then
``\left<s_{ij}\right>``
be an average cross-spectrum between
``i`` and ``j``, and ``\left<s_{i}\right>``, ``\left<s_{j}\right>``
be the auto-spectrum of ``i`` , ``j``. In the frequency domain those
are function of frequency, whereas in the time-frequency domain they are
functions of both time and frequency.

Finally, let ``w`` denote non-negative weights normalized so that their
average is 1. Those are weights for the windows on which the averages
are computed.
Setting all weights equal to 1, gives the unweighted version
of all measures, which is the only supported option in the FD, since
therein windows segmentation is arbitrary and weighting
windows is usually meaningless.

All kinds of coherences estimated in *FourierAnalysis* are summarized
in the following table:

|   kind         |                linear                       |              non-linear              |
|---------------:|:-------------------------------------------:|:------------------------------------:|
|    *total*     | ``\left<w\mid s_{ij} \mid^2\right>\big /\left<ws_is_j\right>`` | ``\left<w\mid s_{ij} \mid^2\right>``  |
|    *real*      | ``\left<w\mid \Bbb R(s_{ij}) \mid^2\right>\big/\left<ws_is_j\right>`` | ``\left<w\mid \Bbb R(s_{ij})\mid^2\right>`` |
|  *imaginary*   | ``\left<w\mid \Bbb C(s_{ij}) \mid^2\right>\big/\left<ws_is_j\right>`` | ``\left<w\mid \Bbb C(s_{ij}) \mid^2\right>`` |
| *instantaneous*| ``\left<w\mid \Bbb R(s_{ij}) \mid^2\right>\big/\left<w(s_is_j-\mid \Bbb C(s_{ij}) \mid^2)\right>`` | ``\left<w\mid \Bbb R(s_{ij}) \mid^2\right>\big/\left<w(1-\mid \Bbb C(s_{ij}) \mid^2)\right>`` |
|   *lagged*     | ``\left<w\mid \Bbb C(s_{ij}) \mid^2\right>\big/\left<w(s_is_j-\mid \Bbb R(s_{ij}) \mid^2)\right>`` | ``\left<w\mid \Bbb C(s_{ij}) \mid^2\right>\big/\left<w(1-\mid \Bbb R(s_{ij}) \mid^2)\right>`` |

The *linear* **total**  coherence is the classical *squared coherence* measure
and its *non-linear* counterpart is known as *phase-locking value* or
*phase coherence*.

The **real** coherence (Pascual-Marqui, 2007) and the
**imaginary** coherence (Nolte *et al.*, 2004), sum up to the total coherence.

The **lagged** coherence has been proposed by Pascual-Marqui (2007).
As a completion of the table, we here also define the **instantaneous**
coherence in an analogous way.

Corresponding *linear* and *non-linear* measures are the same, but are computed
on linear and non-linear cross-spectra in the FD and crossed analytic-signal in the TFD, respectively.
In fact, the non-linear expressions are obtained setting ``s_i=s_j=1``,
which is the case for non-linear cross-spectra and analytic signals. See
[CrossSpectra](@ref) and [TFAnalyticSignal](@ref).
For further explanations on the linearity of those objects, see Congedo (2018).

### interpreting coherence

The real part of the cross-spectrum, named the *co-spectrum*, describes the
*synchronous* synchronization, that is, in-phase synchronization or
out-of-phase synchronization. The imaginary part, named the
*quadrature spectrum*, describes the *asynchronous* synchronization,
that is, with a quarter of a cycle lead or lag (see Congedo, 2018).
Consequently,
- the *total* coherence is sensitive to all kinds of synchronization,
- the *real* coherence is sensitive only to synchronous synchronization,
- the *instantaneous* coherence is sensitive to all but asynchronous synchronization,
- the *imaginary* coherence is sensitive only to asynchronous synchronization,
- the *lagged* coherence is sensitive to all but synchronous synchronization.

## Coherence

**Categories**: [data objects](@ref), [FDobjects](@ref)

In the time-frequency domain coherence estimates are given as Julia Matrix
objects. In the frequency-domain they are encapsulated in the
following structure:

```
struct Coherence
   <same fields of the CrossSpectra structure>
end
```

This object has the same structure of the [CrossSpectra](@ref) object,
with the difference that its `.y` data field holds real matrices,
whereas for cross-spectra holds complex matrices.

By convention the diagonal elements of all FD coherence matrices are
filled with 1.0. Note that the FFT coefficients for the DC level
and the Nyquist frequency are real, therefore for those frequencies the
identity matrix is returned for the *imaginary* and *lagged* coherence.
Also, note that at these frequencies the *instantaneous* coherence is equal
to the *real* coherence. In the time-frequency domain,
we do not have discrete frequencies but a filter bank,
thus no special behavior appears at the edges.

A vector of *Coherence* objects is of type [CoherenceVector](@ref).

**Methods for Coherence and CoherenceVector objects**

|      method          | Coherence    | CoherenceVector    |
|:---------------------|:------------:|:------------------:|
| [`bands`](@ref)      |     ✔        |         ✔         |
| [`extract`](@ref)    |     ✔        |         ✔         |
| [`mean`](@ref)       |     ✔        |         ✔         |
| [`smooth`](@ref)     |     ✔        |         ✔         |
| [`sameParams`](@ref) |              |         ✔         |


**Generic Constructors**:

In order to construct a *Coherence* object in the frequency domain from
multivariate data, *FourierAnalysis* provides two [`coherence`](@ref)
constuctors from raw data two constuctors building coherence matrices
from [CrossSpectra](@ref) and [CrossSpectraVector](@ref) objects.
Those are what you will use in practice most of the time.

Manual constructors are also possible, for which you have to provide
appropriate arguments. Generic constructors follow exactly the same syntax
as the generic constructors of [CrossSpectra](@ref) objects.

In order to estimate coherence in the time-frequency domain from
sets of univariate data, *FourierAnalysis* provides two more
[`coherence`](@ref) functions.

**Constructors from data**:

```@docs
coherence
```

**References**:

M. Congedo (2018),
[Non-Parametric Synchronization Measures used in EEG and MEG](https://hal.archives-ouvertes.fr/hal-01868538),
Technical Report. GIPSA-lab, CNRS, University Grenoble Alpes, Grenoble INP.

G. Nolte *et al.* (2004),
[Identifying true brain interaction from EEG data using the imaginary part of coherency](https://www.researchgate.net/file.PostFileLoader.html?id=55e56d5260614bed268b45e5&assetKey=AS%3A273843090853888%401442300689008).
Clin Neurophysiol, 115(10), 2292-307.

R. Pascual-Marqui (2007),
[Instantaneous and lagged measurements of linear and nonlinear dependence between groups of multivariate time series: frequency decomposition](https://arxiv.org/ftp/arxiv/papers/0711/0711.1455.pdf),
arXiv:0711.1455.
