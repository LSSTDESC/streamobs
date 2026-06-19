#!/usr/bin/env python
"""
Probabilistic samplers.
"""

import numpy as np
import scipy.stats

from streamobs.functions import (CubicSplineInterpolation,
                                 FileCubicSplineInterpolation,
                                 FileInterpolation,
                                 FileLinearDensityCubicSplineInterpolation,
                                 Interpolation,
                                 LinearDensityCubicSplineInterpolation,
                                 Sinusoid)


def sampler_factory(type_, **kwargs):
    """Create a sampler with given kwargs.

    Parameters
    ----------
    type_ : sampler type
    kwargs: passed to sampler init

    Returns:
    -------
    sampler : the sampler
    """
    if type_ == "uniform":
        sampler = UniformSampler(**kwargs)
    elif type_ == "gaussian":
        sampler = GaussianSampler(**kwargs)
    elif type_ == "sinusoid":
        sampler = SinusoidSampler(**kwargs)
    elif type_ in ("interpolation"):
        sampler = InterpolationSampler(**kwargs)
    elif type_ in ("file", "fileinterpolation"):
        sampler = FileInterpolationSampler(**kwargs)
    elif type_ in ("cubicsplineinterpolation"):
        sampler = CubicSplineInterpolationSampler(**kwargs)
    elif type_ in ("filecubicsplineinterpolation"):
        sampler = FileCubicSplineInterpolationSampler(**kwargs)
    elif type_ in ("lineardensitycubicsplineinterpolation"):
        sampler = LinearDensityCubicSplineInterpolationSampler(**kwargs)
    elif type_ in ("filelineardensitycubicsplineinterpolation"):
        sampler = FileLinearDensityCubicSplineInterpolationSampler(**kwargs)
    else:
        raise Exception(f"Unrecognized sampler: {type_}")

    return sampler


def inverse_transform_sample(vals, pdf, size, rng=None):
    """Perform inverse transform sampling

    Parameters
    ----------
    vals: value at which pdf is measured
    pdf : pdf value
    size : number of stars to sample
    rng : numpy.random.Generator, optional
        Random number generator used to draw the underlying uniform variates.
        A fresh, unseeded generator is created if omitted.

    Returns
    -------
    samples : samples of vals
    """
    rng = np.random.default_rng() if rng is None else rng

    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    fn = scipy.interpolate.interp1d(
        cdf, list(range(0, len(cdf))), bounds_error=False, fill_value=0.0
    )
    x_new = rng.uniform(size=np.rint(size).astype(int))
    index = np.rint(fn(x_new)).astype(int)
    return vals[index]


############################################################
# Probabilistic Samplers
############################################################


class Sampler(object):
    def __init__(self, *args, **kwargs):
        pass

    def pdf(self, values):
        return self._pdf(values)

    def sample(self, size, rng=None):
        pass


class ScipySampler(Sampler):
    def __init__(self, rv):
        self._rv = rv

    def pdf(self, x):
        return self._rv.pdf(x)

    def sample(self, size, rng=None):
        return self._rv.rvs(size=int(size), random_state=rng)


class UniformSampler(ScipySampler):
    """Sample from uniform distribution."""

    def __init__(self, xmin, xmax):
        """Uniform distribution sampled between xmin and xmax.

        Parameters
        ----------
        xmin : minimum value of sampled range.
        xmax : maximum value of sampled range.

        Returns
        -------
        self
        """
        rv = scipy.stats.uniform(loc=xmin, scale=xmax - xmin)
        super().__init__(rv)


class GaussianSampler(ScipySampler):
    """Sample from Gaussian."""

    def __init__(self, mu, sigma):
        """Gaussian distribution described by mean and sigma.

        Parameters
        ----------
        mu : mean of Gaussian
        sigma : sigma of Gaussian

        Returns
        -------
        self
        """
        rv = scipy.stats.norm(loc=mu, scale=sigma)
        super().__init__(rv)


class InterpolationSampler(Sampler):
    """Sample from interpolated function."""

    def __init__(self, xvals, yvals, **kwargs):
        """Interpolation sampler.

        Parameters
        ----------
        xvals : x-values of interpolation
        yvals : y-values of interpolation
        kwargs : passed to baseclass

        Returns
        -------
        sampler
        """
        self.interp = Interpolation(xvals, yvals, **kwargs)

    def sample(self, size, nsteps=1e5, rng=None):
        xvals = np.linspace(self.interp.xmin, self.interp.xmax, int(nsteps))
        pdf = self.interp(xvals)
        return inverse_transform_sample(xvals, pdf, size=size, rng=rng)


class FileInterpolationSampler(InterpolationSampler):

    def __init__(self, filename, columns=None):
        self.interp = FileInterpolation(filename, columns)


class SinusoidSampler(InterpolationSampler):

    def __init__(self, **kwargs):
        """Sample from Sinusoid function.

        Parameters
        ----------
        Passed to Sinusoid function.

        Returns
        -------
        sampler
        """
        self.interp = Sinusoid(**kwargs)


class CubicSplineInterpolationSampler(InterpolationSampler):
    def __init__(self, nodes, node_values):
        """Sample using Cubic Spline Interpolation method.

        Parameters
        ----------
        nodes : array-like
            The nodes for the interpolation.
        node_values : array-like
            The values at the nodes.

        Returns
        -------
        sampler
        """
        self.interp = CubicSplineInterpolation(nodes, node_values)


class FileCubicSplineInterpolationSampler(InterpolationSampler):
    def __init__(self, filename, stream_name=None, spline_type=None):
        """Sample using Cubic Spline Interpolation from a file.

        Parameters
        ----------
        filename : str
            The name of the file containing interpolation data.
        stream_name : str, optional
            The name of the stream within the file.
        type : str, optional
            Type of data stream.

        Returns
        -------
        sampler
        """
        self.interp = FileCubicSplineInterpolation(filename, stream_name, spline_type)


class LinearDensityCubicSplineInterpolationSampler(InterpolationSampler):
    def __init__(
        self, intensity_nodes, intensity_node_values, spread_nodes, spread_node_values
    ):
        """Sample using Linear Density Cubic Spline Interpolation method.

        Parameters
        ----------
        intensity_nodes : array-like
            The nodes for the intensity interpolation.
        intensity_node_values : array-like
            The values at the intensity nodes.
        spread_nodes : array-like
            The nodes for the spread interpolation.
        spread_node_values : array-like
            The values at the spread nodes.

        Returns
        -------
        sampler
        """
        self.interp = LinearDensityCubicSplineInterpolation(
            intensity_nodes, intensity_node_values, spread_nodes, spread_node_values
        )


class FileLinearDensityCubicSplineInterpolationSampler(InterpolationSampler):
    def __init__(self, filename, stream_name):
        """Sample using Linear Density Cubic Spline Interpolation from a file.

        Parameters
        ----------
        filename : str
            The name of the file containing interpolation data.
        stream_name : str
            The name of the stream within the file.

        Returns
        -------
        sampler
        """
        self.interp = FileLinearDensityCubicSplineInterpolation(filename, stream_name)
