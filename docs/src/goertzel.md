# goertzel.jl

This unit implements three **Goertzel's algorithms** for estimating DFT complex
coefficients at a specific frequency.

The result of Goertzel's algorithm for a given frequency is practically
identical to the output of an FFT for that frequency, but requieres only ``t`` operations, where ``t`` is the window length in samples.
Since the complexity of FFT algirthms is ``t*log(t)``, Goertzels' algorithm
are advantageous when only one or a few coefficients are needed.
Furthermore, the estimation of Goertzel's algorithms is not restricted
to a Fourier discrete frequency, but can be in the whole positive real line
up to the Nyquist frequency. Smart block on-line implementations are also possible.

**References**

Goertzel G (1958) [An algorithm for the evaluation of finite trigonometric series](https://pdfs.semanticscholar.org/a5e4/d0faf65627374b1ac82c3c79006d010173c9.pdf).
The American Mathematical Monthly, 65(1), 34-35.

See also [here](https://www.st.com/content/ccc/resource/technical/document/design_tip/group0/20/06/95/0b/c3/8d/4a/7b/DM00446805/files/DM00446805.pdf/jcr:content/translations/en.DM00446805.pdf).

```@docs
goertzel
goertzel_fast
goertzel2
```
