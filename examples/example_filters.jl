#   Unit examples of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2021,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute filter banks.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


####################
# test Filters

using Plots, DSP, FourierAnalysis
x=sinusoidal(2, 8, 128, 512, 0)
f, Y=filterbank(x, 128)
plot(x)
plot!((Y[:, 7]))

# Show amplitude distorsions at the edges and the linear phase response
plot!((x-Y[:, 7]))

f, Y=filterbank(x, 128; filtkind=Butterworth(4))

x2=randn(512)+sinusoidal(2, 8, 128, 512, 0)
f, Y=filterbank(x2, 128)
plot(x)
plot!(x2)
plot!((Y[:, 7]))
plot!((x-Y[:, 7]))

# generate a sinusoidal + noise
f, sr, t = 8, 128, 512
v=sinusoidal(1., f, sr, t, 0)
x=v+randn(t)
flabels, Y=filterbank(x, 128)
flabels, Y=filterbank(x, 128; fmin=4, fmax=32)
flabels, Y=filterbank(x, 128, 4; fmin=4, fmax=32)
flabels, Y=filterbank(x, 128, 4;
                      filtkind=Chebyshev2(8, 10),
                      fmin=4,
                      fmax=16)
# trick for plotting the signal filtered in the band-pass regions
for i=1:size(Y, 2) Y[:, i].+=convert(eltype(Y), i)*1.5 end
mylabels=Array{String}(undef, 1, length(flabels))
for i=1:length(flabels) mylabels[1, i]=string(flabels[i])*" Hz" end
plot(Y; c=:grey, labels=mylabels)
