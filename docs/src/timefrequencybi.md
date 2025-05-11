# timefrequencybi.jl

This unit implements average time-frequency *bivariate measures* based
on unit [timefrequency.jl](@ref) and [tools.jl](@ref).

These measures are bivariate extension of the
*mean amplitude*, *concentration* and *phase concentration*
univariate measures implemented in [timefrequencyuni.jl](@ref),
as explained in details in
[Congedo (2018)](https://hal.archives-ouvertes.fr/hal-01868538v2/document).
Please read the documentation of that unit first, as it contains
useful information on how time-frequency measures are estimated
and explain the notation used here below.

The implemented measures are:

#### (weighted) comodulation

``(w)Com=\left<wr_1r_2\right>\big /\sqrt{\left<wr_1^2\right>\left<wr_2^2\right>}``.

This is the bivariate extension of univariate [(weighted) mean amplitude](@ref)
measure.

#### (weighted) coherence

This is the bivariate extension of the
*concentration* and *phase concentration* measure,
yielding *coherence* and *phase coherence* estimates,
the latter also known as *phase-locking value*.

For a complete account on coherence estimations in the time-frequency
domain and the difference with their frequency domain counterpart
see the [coherence.jl](@ref) unit,
where the documentation of the [`coherence`](@ref) functions is also found.


```@docs
comodulation
```

For coherence estimation in the time-frequency domain see [`coherence`](@ref).
