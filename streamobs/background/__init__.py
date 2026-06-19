"""
Background generation for stellar stream analysis.

This subpackage provides two complementary pipelines:

- **Full injection** (:class:`BackgroundCatalogInjector`): runs a known
  catalog of true stars or galaxies through the complete
  :class:`~streamobs.observed.StreamInjector` pipeline.
- **Light generation** (:class:`LightBackgroundGenerator`): samples from
  precomputed binned color–magnitude distributions stored by
  :class:`BackgroundStorage`, giving fast per-pixel background realizations
  without re-running the full injection.

The top-level :class:`Background` class wraps both modes and lets users
generate stars-only, galaxies-only, or combined backgrounds.

Typical usage — light method with default bundled resources::

    from streamobs.surveys import Survey
    from streamobs.background import Background

    survey = Survey.load('lsst', release='yr5')
    bg = Background(survey, source_type='both', method='light')
    catalog = bg.generate(phi1_limits=(-20, 20), phi2_limits=(-2, 2), gc_frame=frame)
"""

from .background import Background
from .catalog_injector import BackgroundCatalogInjector
from .generator import LightBackgroundGenerator
from .resource_builder import BackgroundResourceBuilder
from .storage import BackgroundStorage

__all__ = [
    "Background",
    "BackgroundCatalogInjector",
    "BackgroundResourceBuilder",
    "BackgroundStorage",
    "LightBackgroundGenerator",
]
