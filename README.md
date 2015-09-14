CBRF - A Matlab Toolbox (Cosine Bell Rate Function) for the computation of neural rate functions

This Matlab toolbox calculates a neural rate function based on cosine bell functions and the time-rescaling theorem. The cosine bell function give a continuous function representation to each spike event time, as Gaussians often do in other algorithms. In contrast to a Gaussian, the cosine bell function has compact support, which in this algorithm is adaptive and depends on local interspike intervals. More details of the CBRF algorithm are available in the Documentation folder.

The input data is assumed to be a set of spike times associated with N "repeat trials" (same stimulus). The rate function is constructed through two stages of spike time transformations: a "time A" transformation designed to account for instantaneous rates within each trial, and a "time B" transformation that attempts to homogenize the rate function by accounting for local rate variability across all trials, after all trials are merged. The ultimate goal is to obtain a set of "time B" event times
that, across trials, have the statistics of a unit rate homogeneous Poisson process, as would be the case (according to the time-rescaling theorem) if the two-stage transformation properly "smoothed out" the local rate variations. The real-time rate function is constructed from these two transformations. In short, the original time rate function can be viewed as a clock that, statistically speaking, ticks at irregular intervals. In time B, the clock ticks at intervals with a constant mean.

A Demo folder illustrating the algorithm on a set of 128 repeat trials associated with the response of a cat X retinal ganglion cell is provided.


