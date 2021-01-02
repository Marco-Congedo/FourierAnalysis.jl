#   Unit examples of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2021,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to compute and plot spectra
#   and how to extract inforation from them.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


using FourierAnalysis, FFTW, LinearAlgebra, Statistics,
      Plots, Plots.PlotMeasures

# add module for reading the two EEG text files to be used ater
push!(LOAD_PATH, @__DIR__)
using IOtxt

# Get EEG file names with complete path (they have extension .txt)
S=getFilesInDir(@__DIR__; ext=(".txt",))

# read the two EEG data files and put them in a Matrix object
X1=readEEG(S[1])
X2=readEEG(S[2])

######## compute the FFT of one epoch of a time series ########
sr, t, f, a = 128, 128, 8, 2
v=sinusoidal(a, f, sr, t, 0.6)
plot(v) # plot the time series
tapering=hann

### The long way using FFTW.jl only:
# get a plan.  The normalization 2/t outputs correct peak amplitude
P=plan_rfft(v)*(2/t)

# get a tapering window
a=taper(tapering, t)

# plot the tapered signal
plot!(v.*a.y)

# do the FFT on tapered data
Ψ=P*(v.*a.y)

# get the amplitude spectrum
Σ=abs.(Ψ)
bar!(Σ[2:end])

# check the amplitude at frequency f and neighbooring bins
pos=f2b(f, sr, t, DC=true)

Σ[pos-1]
Σ[pos]
Σ[pos+1]

# The fast way using FourierAnalysis.jl:
Σ2=spectra(v, sr, t; tapering=tapering, func=√)
bar!(Σ2.y)

# equivalently:
Σ2=spectra(v, sr, t; tapering=taper(hann, t), func=√)
bar!(Σ2.y)

####################################################

# test spectra with odd samples and with DC
sr=65; t=65
f, a = fres(sr, t)*4, 2
v=sinusoidal(a, f, sr, t, 0.6; DC=3)
plot(v) # plot the time series
Σ2=spectra(v, sr, t; tapering=rectangular, DC=true, func=√)
bar!(Σ2.y)

####################################################
### Check amplitude spectra at all Fourier Frequency (rectangular taper) ###
# the amplitudes are in increment of 10 along frequencies
# NB when t is even, correct amplitude for the last frequency is obtained
# only if the sinusoidal has a phase=π/6
sr, t, = 16, 32
t½ = t÷2
V=Matrix{Float64}(undef, t, t½)
for i=1:t½ V[:, i]=sinusoidal(10*i, b2f(i, sr, t), sr, t, π/6) end

# plot the sinusoids in V
Z=deepcopy(V)
for i=2:size(Z, 2), j=1:size(Z, 1) Z[j, i]+=(10.0*i) end
plot(Z, label="")

# long way using FFTW only.jl
P=plan_rfft(V, 1)*(2/t)
Ψ=P*V
Ψ[end, :]*=0.5*2^0.25 # correction for Nyquist frequency
Σ=abs.(Ψ)
bar(Σ[brange(t, DC=true), :], labels="")

# fast way using FourierAnalysis.jl
Σ2=spectra(V, sr, t; tapering=rectangular, func=√, DC=true)
bar(Σ2.y[brange(t, DC=true), :], labels="")

#############################################################################
### Check amplitude spectra on long signals obtained by welch methods
# one sinusoidal is at an exact discrete Fourier Frequency and the other not
# Rectangular window
sr, t, f, a = 128, 128, 10, 0.5
v=sinusoidal(a, f, sr, t*16)+sinusoidal(a, f*3.5+0.5, sr, t*16)+randn(t*16)
Σ=spectra(v, sr, t; tapering=rectangular, func=√)
bar(Σ.y, labels="rectangular")

# harris4 window (default)
Σ2=spectra(v, sr, t; func=√)
bar!(Σ2.y, labels="harris4")

#smooth spectra
Σ3=smooth(blackmanSmoother, Σ2)
bar!(Σ3.y, labels="smoothed")


######## Check amplitude spectrum obtained with dpss (slepian tapers) ########
sr, t = 128, 128
t½ = t÷2
V=Matrix{Float64}(undef, t, t½)
for i=1:t½ V[:, i]=sinusoidal(1, b2f(i, sr, t), sr, t, 0 ) end
# plot half of the the first 7 sinusoids in V
plot(V[1:t½+1, 1:7], color=:gray, legend=false)

Σ=spectra(V, sr, t; tapering=taper(slepian, t; α=2, n=3), func=√)
plot(Σ.y, labels="")
#############################################################

#### Compute Welch amplitude spectra of EEG data ####
t, sr, slide, tapering = 1024, 128, 512, harris4

# gather some attributes to obtain nice spectra plots
spectraArgs=(fmax          = 48,
             left_margin   = 2mm,
             bottom_margin = 2mm,
             xtickfont     = font(10, "Times"),
             ytickfont     = font(10, "Times"))

S=spectra(X1, sr, t; tapering=tapering, slide=sr, func=sqrt)
plot(S; ytitle="Amplitude (\\muV)", spectraArgs...)

#smooth spectra
S2=smooth(blackmanSmoother, S)
plot(S2; ytitle="Amplitude (\\muV)", spectraArgs...)

# extract spectra in alpha range (8Hz to 12Hz) at all electrodes
e=extract(S, (8, 12))

# You can use any combination of integer and real numbers
e=extract(S, (8, 12.5))
e=extract(S, (8.5, 12))

# the following equivalently extract the spectra at 10Hz only (1 bin)
e=extract(S, 10)
e=extract(S, (10))
e=extract(S, (10, 10))

# extract spectra at all frequencies. Equivalent to S.y
e=extract(S, :)

# All these ways to indicate frequency ranges work with the mean function
# as well, e.g.,

# compute average spectra in alpha range for each time-series
bar(mean(S, (8, 12)))
# or
bar(mean(S, (8.0,12.0)))

# this actually extract the power at 10Hz for each time series
# and is equivalent to extract(S, 10)
mean(S, 10)
# thus this computes the average power at 10Hz across time-series
mean(extract(S, 10)) # or mean(mean(S, 10))
# and this computes the mean across time-series and across frequencies in range [8Hz, 12Hz]
mean(extract(S, (8, 12))) # or equivalently mean(mean(S, (8, 12)))

# average spectrum across all frequencies for each time-series
bar(mean(S, :))

# plot average spectra in 2Hz-band-pass regions for all time-series
plot(bands(S, 2))

# plot average spectra in in 2Hz-band-pass regions for series 15-19
plot(bands(S, 2)[:, 15:19])

# plot average spectra in in 2Hz-band-pass regions for series 19
plot(bands(S, 2)[:, 19])

# plot average spectra in alpha range for series 15-19
bar(mean(S, (8,12))[15:19])

# extract spectra in alpha range for series 15-19
e=extract(S, (8, 12))[:, 15:19]

# Compute several spectra altogether
𝐗=[X1, X2]
𝐒=spectra(𝐗, sr, t; tapering=tapering, slide=sr, func=√)
plot(𝐒[1]; ytitle="Amplitude (\\muV)", spectraArgs...)
plot(𝐒[2]; ytitle="Amplitude (\\muV)", spectraArgs...)

# do the same thing using a fast FFTW plan (wait up to 10s for computing the plan)
plan=Planner(plan_exhaustive, 10.0, t, eltype(𝐗[1]))
𝐒_fast=spectra(𝐗, sr, t; planner=plan, slide=sr, func=√)

# compute the average spectra in the alpha range for all time-series and all subjects
mean(𝐒, (8,12))

# compute the average spectra in the alpha range across subjects for all time-series
mean(mean(𝐒, (8,12)))
# compute and plot it
plot(mean(mean(𝐒, (8,12))))

# extract spectra in alpha range for all time-series and all subjects
extract(𝐒, (8, 12))

# if you enter en illegal range, nothing will be done and you will get
# an error in the REPL explaining what is wrong in your argument
extract(𝐒, (0, 12))
mean(𝐒, (0, 128))

# extract 4Hz-band-pass average spectra for all electrodes and all subjects
bands(𝐒, 4)

# Apply smoothing in the above spectra computations
𝐒=spectra(𝐗, sr, t;
    tapering=tapering, smoothing=blackmanSmoother, func=√)
plot(𝐒[1]; ytitle="Amplitude (\\muV)", spectraArgs...)
plot(𝐒[2]; ytitle="Amplitude (\\muV)", spectraArgs...)

# plot the average spectrum across all electrodes for the two files
# using Julia standard mean function
plot(mean(𝐒[1].y[:, i] for i=1:size(𝐒[1].y, 2)))
plot!(mean(𝐒[2].y[:, i] for i=1:size(𝐒[2].y, 2)))

# plot spectra in in 1Hz band-pass regions for all electrodes in 𝐒[1]
plot(bands(𝐒[1], 1))

# use slepian multi-tapering on EEG data
S_sl=spectra(X1, sr, t; tapering=slepians(sr, t, 1.25), func=√)
plot(S_sl; ytitle="Amplitude (\\muV)", spectraArgs...)

# use slepian multi-tapering and get smoothed spectra
S_sl=spectra(X1, sr, t;
    tapering=slepians(sr, t, 1.25), func=√, smoothing=blackmanSmoother)
plot(S_sl; xspace=8, ytitle="Amplitude (\\muV)", spectraArgs...)

# See how the variance of the spectra estimation is lower
# as compared to using a normal tapering window !
S=spectra(X1, sr, t;
    tapering=harris4, func=√, smoothing=blackmanSmoother)
plot(S; ytitle="Amplitude (\\muV)", spectraArgs...)
