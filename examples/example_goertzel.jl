#   Unit examples of the FourierAnalysis Package for julia language
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# ? CONTENTS :
#   This example shows how to use Goertzel's algorithms.

#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #
#                                                                             #
#   ~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~~¤~  #


########## Test goertzel algorithms #############
using FourierAnalysis
sr, t, f, a = 128, 128, 5, 10
v=sinusoidal(a, f, sr, t, 0)
c=goertzel(v, f, sr, t) # should be 0+aim
d=goertzel2(v, f, sr, t) # should be 0+aim
abs2(c) # power, should be a²
abs2(d) # power, should be a²

ftrue=f+0.15
v=sinusoidal(a, ftrue, sr, t, 0 )
c=goertzel(v, f, sr, t) # should be 0+aim
d=goertzel2(v, ftrue, sr, t) # should be 0+aim
abs2(c) # power, should be a²
abs2(d) # power, should be a²
#################################################
