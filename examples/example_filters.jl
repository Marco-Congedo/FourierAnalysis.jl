#   Unit examples of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 24th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example show how to compute filter banks.

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
