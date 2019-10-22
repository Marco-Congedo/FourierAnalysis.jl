#   Unit examples of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 24th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute and plot time-frequency analytic signal,
#   amplitude and phase objects and how to extract inforation from them.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots, Plots.Measures

# add module for reading the two EEG text files to be used ater
push!(LOAD_PATH, @__DIR__)
using IOtxt

# Get EEG file names with complete path (they have extension .txt)
S=getFilesInDir(@__DIR__; ext=(".txt",))

# read the two EEG data files and put them in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])



####################
# test timefrequency

using Plots, FourierAnalysis
sr, t, bandwidht=128, 512, 2
h=taper(harris4, t)
x1=sinusoidal(10, 8, sr, t, 0)
x2=sinusoidal(10, 19, sr, t, 0)
x=Vector((x1+x2).*h.y+randn(t))
y1=sinusoidal(10, 6, sr, t, 0)
y2=sinusoidal(10, 16, sr, t, 0)
y=Vector((y1+y2).*h.y+randn(t))

# Vector of Time-Frequency Object for x and y
𝒀 = TFanalyticsignal([x, y], sr, sr*4)

# compute the mean in a TF region (8:12Hz and samples 1:128) for the first object
y=mean(𝒀[1], (8, 12), (1, 128))

# extract a TF region (8:12Hz and samples 1:128) for the first object
Z=extract(𝒀[1], (8, 12), (1, 128))

# compute the mean in a TF region (8:12Hz and samples 1:128) for the two objects
Y=mean(𝒀, (8, 12), (1, 128))
# do the same computation without checking homogeneity of the two objects in 𝒀
Y=mean(𝒀, (8, 12), (1, 128); check=false)

# extract the data in a TF region (8:12Hz and samples 1:128) for the two objects
E=extract(𝒀, (8, 12), (8, 12))
# do the same computation without checking homogeneity of the two objects in 𝒀
E=extract(𝒀, (8, 12), (8, 12); check=false)

# compute time-frequency object for vector x
Y=TFanalyticsignal(x, sr, t, bandwidht; fmax=32)

# gather useful attributes to obtain nice heatmpap plots
tfArgs=(right_margin = 2mm,
        top_margin = 2mm,
        xtickfont = font(10, "Times"),
        ytickfont = font(10, "Times"))

# plot the real part of the AS
heatmap(tfAxes(Y)..., real(Y.y); c=:pu_or, tfArgs...)

# ...the imaginary part of the AS
heatmap(tfAxes(Y)..., imag(Y.y); c=:bluesreds, tfArgs...)

# ...the amplitude of the AS
heatmap(tfAxes(Y)..., amplitude(Y.y); c=:amp, tfArgs...)


# or generate a TFAmplitude object
A=TFamplitude(Y)
# and plot it (with different colors)
heatmap(tfAxes(A)..., A.y; c=:fire, tfArgs...)

# ...the amplitude of the AS smoothed in the freq. dim.
heatmap(tfAxes(Y)...,
        amplitude(smooth(hannSmoother, noSmoother, Y).y);
        c=:amp, tfArgs...)

# ...the amplitude of the AS smoothed in freq. and time
heatmap(tfAxes(Y)...,
        amplitude(smooth(blackmanSmoother, blackmanSmoother, Y).y);
        c=:amp, tfArgs...)

# ...the phase of the AS
heatmap(tfAxes(Y)..., phase(Y.y);
        c=:bluesreds, tfArgs...)

#or, generate a TFPhase object
ϴ=TFphase(Y)
# and plot it
heatmap(tfAxes(ϴ)..., ϴ.y;
        c=:pu_or, tfArgs...)

# compute and plot phase in [0, 2π]
heatmap(tfAxes(Y)..., TFphase(Y; func=x->x+π).y;
        c=:amp, tfArgs...)

# compute and plot unwrapped phase
heatmap(tfAxes(Y)..., TFphase(Y; unwrapped=true).y;
        c=:amp, tfArgs...)

# smooth time-frequency analytic signal: smooth frequency
Z=smooth(blackmanSmoother, noSmoother, Y)
# plot amplitude of smoothed analytic signal
heatmap(tfAxes(Z)..., amplitude(Z.y); c=:amp, tfArgs...)

# not equivalently (!), you can create an amplitude object and smooth them:
# in this case the amplitude is smoothed, not the analytic signal
A=smooth(blackmanSmoother, noSmoother, TFamplitude(Y))
heatmap(tfAxes(A)..., A.y; c=:amp, tfArgs...)

# Smoothing raw phase estimates is unappropriate
# since the phase is a discontinous function, however it makes sense to smooth
# phase if the phase is unwrapped.
heatmap(tfAxes(Y)...,
        smooth(blackmanSmoother, noSmoother, TFphase(Y; unwrapped=true)).y;
        c=:amp, tfArgs...)

# smooth AS: smooth both frequency and time
E=smooth(blackmanSmoother, blackmanSmoother, Y)

# plot amplitude of smoothed analytic signal
heatmap(tfAxes(E)..., amplitude(E.y);
        c=:fire, tfArgs...)

# plot phase of smoothed analytic signal
heatmap(tfAxes(E)..., phase(E.y);
        c=:bluesreds, tfArgs...)

# not equivalently (!), create amplitude and phase objects and smooth them
A=smooth(blackmanSmoother, blackmanSmoother, TFamplitude(Y))
heatmap(tfAxes(A)..., A.y;
        c=:fire, tfArgs...)

ϴ=smooth(blackmanSmoother, blackmanSmoother, TFphase(Y, unwrapped=true))
heatmap(tfAxes(ϴ)..., ϴ.y;
        c=:pu_or, tfArgs...)

# smooth again
ϴ=smooth(blackmanSmoother, blackmanSmoother, ϴ)
heatmap(tfAxes(ϴ)..., ϴ.y;
        c=:pu_or, tfArgs...)
# and again ...

# you may also create all these objects already smoothed, for example
# create directly smoothed AS
Y=TFanalyticsignal(x, sr, t, bandwidht;
                   fmax=32,
                   fsmoothing=hannSmoother,
                   tsmoothing=hannSmoother)

# plot amplitude of smoothed analytic signal
heatmap(tfAxes(Y)..., amplitude(Y.y);
   c=:amp, tfArgs...)

# create directly smoothed Amplitude
A=TFamplitude(x, sr, t, bandwidht;
         fmax=32,
         fsmoothing=hannSmoother,
         tsmoothing=hannSmoother)

# plot smoothed amplitude
heatmap(tfAxes(A)..., A.y;
   c=:amp, tfArgs...)

# compute a TFAnalyticSignal object with non-linear AS
Y=TFanalyticsignal(x, sr, t, bandwidht; fmax=32, nonlinear=true)

# check that it is non-linear
Y.nonlinear

# check that the amplitude is now 1.0 everywhere
norm(amplitude(Y.y)-ones(eltype(Y.y), size(Y.y))) # must be zero

# plot non-linear phase
heatmap(tfAxes(Y)..., phase(Y.y);
   c=:bkr, tfArgs...)

# get the center frequencies of TFAmplitude object A
A.flabels

# extract the amplitude in a time-frequency region
extract(A, (2, 10), (1, 256))
# extract the amplitude in a time-frequency region at only one frequency
extract(A, 10, (1, 256))
# extract the amplitude at one temporal sample at one frequency
extract(A, 10, 12)
# or
extract(A, 10.0, 12)
# extract amplitude at one temporal sample in a frequency band-pass
extract(A, (10, 12), 12)
#or
extract(A, (10.0, 12.0), 12)
# extract amplitude at one temporal sample and all frequencies
extract(A, :, 12)
# NB: the extract function work in the same way
#     for objects TFAnayticSignal and TFPhase


# compute the mean in a time-frequency region:
mean(A, (2, 10), (1, 256))
# is equivalent to (but may be less efficient than)
mean(extract(A, (2, 10), (1, 256)))

# using column sign for extracting all time samples
extract(A, (2, 10), :)
# This:
extract(A, :, :)
# is equivalent to this:
A.y
# but if you don't need to extract all frequencies,
# use the extract function to control what frequencies will be extracted:
# This
extract(A, (4, 5), 10)
# is not equivalent to this
A.y[4:5, 10]
# since the `extract` function finds the correct rows corresponding
# to the sought frequencies (in Hz), while A.y[4:5, 10]
# just returns the element [4:5, 10] in the TF amplitude object


mean(A, 10, (1, 256))
mean(A, 10, 12)
# this is not equivalent to the following
A.y[10, 12]
# since the `mean` function also finds the correct rows corresponding
# to the sought frequencies (in Hz), while A.y[10, 12]
# just returns the element [10, 12] in the TF amplitude object

# Although the first center frequency in A is 2Hz, its
# band-pass region is 1-3Hz, therefore the frequency range 1:10 is accepted
mean(A, (1, 10), (1, 256))
# but this result in an error (see the REPL) since 0.5 is out of range
mean(A, (0.5, 10), (1, 256))

# using a colon sign for time range
a=mean(A, (1, 10), :)
# using an integer for time range (indicates one specific sample)
a=mean(A, (1, 10), 16)

# using a colon sign for frequency range
a=mean(A, :, (1, 16))
# using a real number for frequency range
a=mean(A, 8.5, (1, 16))

# This
a=mean(A, :, :)
# is equivalent to this
a=mean(A.y)


###############################################################

Pz=15
𝐱=[X1[:, Pz], X2[:, Pz]] # get the two times-series at electrode Pz
𝐘=TFanalyticsignal(𝐱, sr, t, bandwidht; fmax=32, nonlinear=false)
𝐀=TFamplitude(𝐘)
heatmap(tfAxes(𝐀[1])..., 𝐀[1].y; c=:amp, tfArgs...)
heatmap(tfAxes(𝐀[2])..., 𝐀[2].y; c=:pu_or, tfArgs...)

# plot the power over time from instantaneous amplitude
plot([sum(𝐀[2].y[:, t]) for t=1:size(𝐀[2].y, 2)])
# plot the smoothed power over time from the instantaneous amplitude
A=smooth(noSmoother, blackmanSmoother, 𝐀[2])
plot([sum(A.y[:, t]) for t=1:size(A.y, 2)])

# 𝐀 holds 3 TFAmplitude objects;
# The following will create a 3-Vector of matrices
# holding the amplitude in the TF region defined by frequency 8:12 and samples 1:128:
extract(𝐀, (8, 12), (1, 128))

# use of a column is allowed for the frequency range as all the amplitudes have the same number of frequencies:
extract(𝐀, :, (1, 128))

# however, the use of a column is not allowed for the time range as the number
# of samples is different in the three TFAmplitude objects
extract(𝐀, (8, 12), :)

# the same goes for the mean function
mean(𝐀, (8, 12), :)

# to obtain the grand-average across all objects:
a=mean(mean(𝐀, (8, 12), (1, 128)))
