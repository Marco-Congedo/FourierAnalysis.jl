# timefrequencyuni.jl

This unit implements average time-frequency *univariate measures* based
on unit [timefrequency.jl](@ref) and [tools.jl](@ref).

These measures here implemented are weighted version of the measures
described in [Congedo (2018)](https://hal.archives-ouvertes.fr/hal-01868538/document).

They can be obtained from [TFAnalyticSignalVector](@ref) objects
or from raw data. Some of them can be obtained also from
[TFAmplitudeVector](@ref) objects.

The measures are always estimated as the average across several analytic signals
at all points in a time-frequency region or as the grand-average
across several analytic signals of all the points in a time-frequency region.
In the following an average process is denoted by angle brackets
``\left<\cdot\right>`` and is used to indicate
generically both averaging processes.
Analytic signal in a time-frequency region,
it does not matter if this is actually a single point,
a vector or a matrix, is denoted such as

``z=x+𝑖y=re^{𝑖𝜑}``,

where ``𝑖`` is the imaginary unit, ``r=\mid z \mid`` is the amplitude
(modulus) of ``z`` and ``𝜑=\textrm{ArcTan}\space (x/y)`` the phase (argument)
of ``z``.

Also, ``w`` denotes non-negative weights normalized so that their average is
``1.0``. ``w`` weights the different analytic signals on which the average
is computed. Setting all weights equal to ``1.0`` gives the unweighted version
of all measures.

Some of the measures come in a *linear* and *non-linear* form.
It is adopted throughout the convention of prepending 'phase'
to the name of a measure to signal it is non-linear.
The reason is that non-linear measures are not sensitive to amplitude,
but only to phase.
See [Congedo (2018)](https://hal.archives-ouvertes.fr/hal-01868538/document)
for a throughout discussion.

The implemented measures are:

#### (weighted) mean amplitude

``(w)MAmp=w\left<\mid z \mid\right>\big /\left<w\right>=w\left<r\right>\big /\left<w\right>``.

#### (weighted) concentration

``(w)Con=\mid\left<wz\right>\mid\big /\left<w\right>=\mid\left<wre^{𝑖𝜑}\right>\mid\big /\left<w\right>``.

#### (weighted) phase concentration

``(w)PCon=\mid\left<we^{𝑖𝜑}\right>\mid\big /\left<w\right>``.

This is the non-linear version of the *(weighted) concentration*.
In the litetrature it is also known as
*circular mean resultant length*, *inter-trial phase coherence*,
*inter-trial phase clustering* and
*phase coherence*, among other names (Congedo, 2018).

#### (weighted) mean direction

``(w)MDir=\textrm{ArcTan}\space \big((\left<wy\right>/\left<wx\right>)\big/ \left<w\right>\big)``

#### (weighted) phase mean direction

``(w)PMDir=\textrm{ArcTan}\space \big((\left<wy/r\right>/\left<wx/r\right>)\big/ \left<w\right>\big)``

This is the non-linear version of the *(weighted) phase mean direction*


```@docs
meanAmplitude
concentration
meanDirection
```
