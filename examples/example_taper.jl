#   Unit examples of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2021,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows hot to create and plot tapering window for FFT.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


############ plot tapering windows #############
using FourierAnalysis, Plots

# using the standard plot function
tapers=[TaperKind(i) for i=1:8]
X=zeros(t, 8)
for i=1:8 X[:, i] = taper(tapers[i], t).y end
mylabels=Array{String}(undef, 1, 8)
for i=1:8 mylabels[1, i]=string(tapers[i]) end
plot(X; labels=mylabels)

# discrete prolate spheroid sequences
H=taper(slepian, t, α=4, n=7)
# This below is a plot recipe declared in the `recipes.jl unit`
plot(H)


################################################
