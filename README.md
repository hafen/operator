# Operator: Loess and STL Operator Matrix Computation for Time Series

This package contains methods to compute operator matrices for loess and STL decompositions.  It was written during graduate work at Purdue University and no development has occurred while I have been at PNNL.  I am making it public here because I never got around to putting it on CRAN.  

Some of its features:

- Model selection (ANOVA, Cp statistic)
- Prediction (with plot methods)
- Confidence intervals (with plot methods)
- "Blending" to lower-degree polynomial at endpoints
- Computation and plotting of frequency response functions for different parameter choices of each component in a decomposition, useful for diagnosing whether there is interference between seasonal and trend smoothing
- Experimental support for including an ARMA component in the mix

Note that these methods consume a lot of time and memory for very long series, and due to the locality of STL and Loess, one should consider applying them to only the portion of the series of interest.

For functions to apply and plot STL decompositions, see the [stl2](http://github.com/hafen/stl2) package.

## References

- [Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990). STL: A seasonal-trend decomposition procedure based on loess. *Journal of Official Statistics*, 6(1), 3-73.](http://cs.wellesley.edu/~cs315/Papers/stl%20statistical%20model.pdf)
- [Hafen, R. P. "Local regression models: Advancements, applications, and new methods." (2010).](http://search.proquest.com/docview/749923640)

## License

This software is released under the BSD license.  Please read the [license](https://github.com/hafen/operator/blob/master/LICENSE.md) document.

