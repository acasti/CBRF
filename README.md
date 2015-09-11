CBRF - A Matlab Toolbox (Cosine Bell Rate Function) for the computation of neural rate functions

This Matlab toolbox calculates a neural rate function based on cosine bell functions and the time-rescaling theorem.
The cosine bell functions serve the purpose that Gaussians often do in giving a continuous function representation of
of event time. However, the CBRF algorithm is more complicated than averaging Gaussians across trials.  See the 
Documentation folder.

The input data is assumed to be a set of spike times associated with N "repeat trials" (same stimulus).  The rate function
is constructed through two stages of spike time transformations:  a "time A" transformation designed to account for 
instantaneous rates within each trial, and a "time B" transformation that attempts to homogenize the rate function by
accounting for local rate variability across all trials.  The ultimate goal is to obtain a set of "time B" event times 
that, across trials, have the statistics of a unit rate homogeneous Poisson process, as would be the case (according to
the time-rescaling theorem) if the two-stage transformation properly "smoothed out" the local rate variations.  The 
real-time rate function is constructed from these two transformations.

A Demo folder illustrating the algorithm on a set of 128 repeat trials associated with the response of a cat X retinal 
ganglion cell is provided.

