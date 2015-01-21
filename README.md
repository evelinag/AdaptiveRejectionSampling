# Derivative-free adaptive rejection sampling
Adaptive rejection sampling in F#

Adaptive rejection sampling (derivative-free). Based on matlab implementation from 
[PMTK3 library](https://github.com/probml/pmtk3).

Compared to the original implementation in matlab, this version computes 
probabilities in log space which avoids numerical problems. 
 
# Parameters

* `func` - Function which computes logarithm of a likelihood function
* `a`, `b`  - Starting points, a < b. 
* `domain` - Tuple of `func` boundaries. If the left boundary is -infinity, 
    derivative of `f` at `a` must be positive. If the right boundary is infinity,
    derivative of `f` at `b` must be negative. 
* `nSamples` - Total number of samples that the sampler should return.

# References
* Gilks W. 1992. Derivative-free adaptive rejection sampling for Gibbs sampling. In: Bernardo J.M., Berger J.O., Dawid A.P., and Smith A.F.M. (Eds.), Bayesian Statistics 4. Oxford University Press, xford, pp. 641–665.
* Gilks W.R. and Wild P. 1992. Adaptive rejection sampling for Gibbs sampling. Applied Statistics 41: 337–348.
