#   Unit "plots" of the FourierAnalysis Package for julia language
#   v 0.2.0 - last update 20th of October 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements recipes for plotting data generated
#   by FourierAnalysis without having any dependency to plot
#   packages

########################################################

@recipe function f(::Type{Taper}, taper::Taper)
    title               --> taperinfo(taper)
    legend              --> false
    dpi                 --> 300
    background_color    --> :white
    gridcolor           --> :white
    taper.y
end



@recipe function f(::Type{Spectra}, S::Spectra;
            fmax   = S.sr/2,
            xspace = 4,
            ytitle = "Power (\\muV²)")

    b=f2b(fmax, S.sr, S.wl)
    xticstring=["$(xspace*i)" for i=1:fmax÷xspace]
    pushfirst!(xticstring, "$(fres(S.sr, S.wl))")
    xticrange=[1:S.wl/S.sr*xspace:b+1;]
    for i=2:length(xticrange) xticrange[i]=xticrange[i]-1 end

    legend         --> false
    dpi            --> 300
    xlabel         --> "Frequency (Hz)"
    xticks         --> (xticrange, xticstring)
    ylabel         --> ytitle
    delete!(plotattributes, :fmax)
    delete!(plotattributes, :xspace)
    delete!(plotattributes, :ytitle)
    S.y[1:(f2b(fmax, S.sr, S.wl; DC=S.DC)), :]
end


"""
Generate labels for the axes to be used in heatmap plots of [TFobjects](@ref).
On the y-axis the labels correspond to the center frequencies of the band-pass
filter bank used to obtain the time-frequency representation, that is,
to the .flabels field of the object.

See [plot time-frequency objects](@ref) for examples.
"""
function tfAxes(Y::TFobjects)
    function formatf(r::Real)
        s="$r"
        ss=split(s, ".")
        length(ss)>1 && parse(Float64, ss[2])==0.0 ? (return ss[1]) : (return s)
    end
    xs = [string(i) for i = 1:size(Y.y, 2)]
    ys = Vector{String}([formatf(f) for f in Y.flabels])
    return xs, ys
end
