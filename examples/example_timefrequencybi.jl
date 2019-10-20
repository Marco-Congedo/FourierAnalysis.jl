#   Unit examples of the FourierAnalysis Package for julia language
#   v 0.0.1 - last update 24th of September 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute univariate measures in
#   the time-frequency domain.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


using FourierAnalysis, FFTW, LinearAlgebra, Statistics, Plots

# add module for reading the two EEG text files to be used ater
push!(LOAD_PATH, @__DIR__)
using IOtxt

# Get EEG file names with complete path (they have extension .txt)
S=getFilesInDir(@__DIR__; ext=(".txt",))

# read the two EEG data files and put them in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])

#############################################################################
# Bivariate measures
# see Congedo, 2018: https://hal.archives-ouvertes.fr/hal-01868538/document

Pz=15
Fz=5
𝐱₁=[X1[:, Pz], X2[:, Pz]] # get the two times-series at electrode Pz
𝐱₂=[X1[:, Fz], X2[:, Fz]] # get the two times-series at electrode Fz
sr, wl, bandwidht=128, 512, 2
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=false)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=false)
𝐀₁=TFamplitude(𝐘₁)
𝐀₂=TFamplitude(𝐘₂)


# Comodulation Eq. 0.11 (Com)
# compute the Com averaging in a TF region from a TFAnalyticSignalVector object
# NB If you compute it from Analytic Signal or Amplitude objects, those
# must be linear (note the nonlinear=false above)
Com=comodulation(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=mean)
# compute the Com averaging in a TF region from a TFAmplitudeVector object
Com=comodulation(𝐀₁, 𝐀₂, (8, 12), (1, 512); mode=mean)
# compute the Com averaging in a TF region directly from data
# In this care you don't have to worry about linearity of the analytic signal
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean)
# You can compute comodulation from smoothed amplitude:
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht;
                 mode=mean,
                 fsmoothing=blackmanSmoother,
                 tsmoothing=blackmanSmoother)


# you can go faster pre-computing a FFTW plan.
# This is aslo useful when you have to call the comodulation function several times
plan=Planner(plan_patient, 5, wl, Float64, true)
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean, planner=plan)


# compute the Com in a TF region from a TFAnalyticSignalVector object
Com=comodulation(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=extract)
# compute the Com in a TF region from a TFAmplitudeVector object
Com=comodulation(𝐀₁, 𝐀₂, (8, 12), (1, 512); mode=extract)
# compute the Com in a TF region directly from data
Com=comodulation(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=extract)

# All these operations can be done also for coherence measures, for example
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=mean)
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=extract)

# Compute all 5 coherence types
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=extract, allkinds=true)


# phase coherence (phase-locking value)
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=true)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=true)
Coh=coherence(𝐘₁, 𝐘₂, (8, 12), (1, 512); mode=mean, nonlinear=true)

# or directly from data (no need to compute non-linear analytic signal in this case)
Coh=coherence(𝐱₁, 𝐱₂, sr, wl, (8, 12), (1, 512), bandwidht; mode=mean, nonlinear=true)

# and also for non-linear meausures
# compute non-linear analyticSignal
𝐘₁=TFanalyticsignal(𝐱₁, sr, wl, bandwidht; fmax=32, nonlinear=true)
𝐘₂=TFanalyticsignal(𝐱₂, sr, wl, bandwidht; fmax=32, nonlinear=true)
# although you are allowed to compute the amplitude of non-linear
# analytic sygnal this way, it does not make sense as the amplitude
# is 1.0 everywhere, as you can check
𝐀₁=TFamplitude(𝐘₁)
𝐀₂=TFamplitude(𝐘₂)

##################################################################
