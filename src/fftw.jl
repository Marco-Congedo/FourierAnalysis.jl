#   Unit "fftw.jl" of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019-2025,
#   Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This unit implements a convenient interface to FFTW.jl and AbstractFFTs.jl


# ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤

# flags to be passed as `planflag` parameter in several functions
# see www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags for details
const plan_estimate        = UInt32(64)      # very small search =(1 << 6)
const plan_measure         = UInt32(0)       # small search
const plan_patient         = UInt32(32)      # large search =(1 << 5)
const plan_exhaustive      = UInt32(8)       # very large search =(1 << 3)
const plan_conserve_memory = UInt32(4)       # = (1 << 2)
const plan_wisdom_only     = UInt32(2097152) # = (1 << 21)
const plan_unaligned       = UInt32(2)       # = (1 << 1)
const plan_preserve_input  = UInt32(16)      # do not copy the input data = (1 << 4)

"""
FFTW plans to be used by all functions in **Fourier Analysis**
are incapsulated in this structure. This is the only *operator object*
created by this package, all others being *data objects*.

Spectral and cross-spectral computations (frequency domain objects)
need only a forward plan (real FFT),
while objects based on the analytic signal (time-frequency domain objects)
need both a forward (`.p `) and a backward (`.ip `) plan (complex iFFT).

Argument `flags` must be one of the constants here above.
Basic usage of the flags involves a trade-off between the time needed to compute
the planner and the efficacy for computing the FFTs. The following flags
are sorted in ascending order of time needed to be computed and efficacy:
`plan_estimate`, `plan_measure`, `plan_patient`,
`plan_exhaustive`.
By default *FourierAnalysis* adopts `plan_estimate`,
which allows the quickest search, but the least optimal choice.
Other flags may be worth if several FFT computations are to be done
with the same setting.
See [here](http://www.fftw.org/fftw3_doc/Planner-Flags.html) for details.

Argument `timelimit` is the maximum time (in seconds) allowed to FFTW
for optimizing the plans. Setting this to `-1.0` amounts to imposing
no time limits.

FFTW plans are computed for a given window length `wl` and data type `type`.

`Planners` objects may be passed as an argument to constructors of
[FDobjects](@ref) and [TFobjects](@ref) by using one of the `Planner` constuctor
here below.

**Constructors**

```julia
Planner(flags :: UInt32,
    timelimit :: Union{Int, Float64},
           wl :: Int,
         type :: Type,
           bw :: Bool = false)
```

Use this to create a Planner object, passing as argument a FFTW `flags` constant
and `timelimit`, the window length `wl` and the type of the input data `type`,
which must be real, e.g., Float64.

If `bw` is false (default), a dummy backward plan is created,
otherwise a backward plan is created for the complex data type corresponding
to `type` (e.g., if `type` is Float64, it will created for ComplexF64 data.)

For example, suppose 𝐗 is a vector of many matrices of multivariate time-series
sampled at 128 samples per second and that we want to compute the spectra
for all of them using a 256-point FFT. We first create a plan by

```julia
p=Planner(plan_exhaustive, 10.0, 256, eltype(𝐗[1]))
```

Then we invoke the [Spectra](@ref) function passing the plan as argument:

```julia
𝐒=spectra(𝐗, sr, t; planner=p)
```

A shorter construction syntax is available when only the
forward plan is needed and the type of the data is Float64:

```julia
Planner(flags :: UInt32,
    timelimit :: Union{Int, Float64},
           wl :: Int)
```

For example, the line above could have been written more shortly as

```julia
p=Planner(plan_exhaustive, 10.0, 256)
```

In order to create also a backward plan you would use instead

```julia
p=Planner(plan_exhaustive, 10.0, 256, eltype(𝐗[1]), true)
```


"""
struct Planner
    flags     ::  UInt32            # see
    timelimit ::  Float64           # =for no time limit set to -1.0
    wl        ::  Int               # window length
    type      ::  Type              # element type of forward planner (e.g., Float64)
    p         ::  FFTW.rFFTWPlan    # forward plan (FFT)
    ip        ::  FFTW.cFFTWPlan    # backward plan (inverse FFT)
end


Planner(flags::UInt32, timelimit::Union{Int, Float64}, wl::Int, type::Type, bw::Bool=false) =
    bw ? Planner(flags, Float64(timelimit), wl, type, plan_rfft(zeros(type, wl);
                 flags=flags, timelimit=timelimit), plan_bfft(zeros(eltype(Complex(type(0))), wl); flags=flags, timelimit=timelimit)) :
         Planner(flags, Float64(timelimit), wl, type, plan_rfft(zeros(type, wl);
                 flags=flags, timelimit=timelimit), plan_bfft(zeros(ComplexF64, 2)))

# create a planner with only a forward plan defaulting to the Float64 type
Planner(flags::UInt32, timelimit::Union{Int, Float64}, wl::Int) =
    Planner(flags, Float64(timelimit), wl, Float64, plan_rfft(zeros(Float64, wl), flags=flags, timelimit=timelimit), plan_bfft(zeros(ComplexF64, 2)))

# example of full constructor:
# plan=Planner( plan_exhaustive, 10, 256, Float64,
#               plan_rfft(zeros(Float64, 256), plan_exhaustive, timelimit=10),
#               plan_bfft(zeros(ComplexF64, 256), plan_exhaustive, timelimit=10))

# fast constructor writing only the forward planner passed as argument.
# all other fields are dummy. This is used internally by spectra and
# cross-spectra functions
_Planner(p::FFTW.rFFTWPlan)=Planner(plan_estimate, -1.0, 0, Float64, p, plan_bfft(zeros(ComplexF64, 2)))

# dummy FFTW planner used for default settings. Only for internal use
getplanner=_Planner(plan_rfft(zeros(2)))


_getflags(flags::UInt32) =
 if     flags==UInt32(64) return "plan_estimate"
 elseif flags==UInt32(0)  return "plan_measure"
 elseif flags==UInt32(32) return "plan_patient"
 elseif flags==UInt32(8)  return "plan_exhaustive"
 elseif flags==UInt32(4)  return "plan_conserve_memory"
 elseif flags==UInt32(2097152) return "plan_wisdom_only"
 elseif flags==UInt32(2)  return "plan_unaligned"
 elseif flags==UInt32(16) return "plan_preserve_input"
 elseif @error 📌*", the FFTW flags does not exist"
 end

 # ++++++++++++++++++++  Show override  +++++++++++++++++++ # (REPL output)
function Base.show(io::IO, ::MIME{Symbol("text/plain")}, p::Planner)
    println(io, titleFont, "∿ Planner type")
    #println(io, "□  □    □      □        □           □", defaultFont)
    println(io, separatorFont, "⭒  ⭒    ⭒      ⭒        ⭒           ⭒", defaultFont)
    println(io, "FFTW flags        (.flags): "*_getflags(p.flags)*" ($(p.flags))")
    p.timelimit==-1.0 ? println(io, "FFTW timelimit(.timelimit): no limit") :
                        println(io, "FFTW timelimit(.timelimit): $(p.timelimit)")
    println(io, "window length        (.wl): $(p.wl)")
    println(io, "element type       (.type): $(p.type)")
    println(io, "forward plan          (.p): ")
    println(io, greyFont, "$(p.p)", defaultFont)
    println(io, "backward plan        (.ip): ")
    println(io, greyFont, "$(p.ip)", defaultFont)
end
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤ ¤
