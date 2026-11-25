"""SiteWindow: site-energy landscapes and deliverable-capacity limits for hydrogen in sp2 carbon."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("sitewindow")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

__all__ = ["__version__"]
