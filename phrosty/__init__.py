from importlib.metadata import version, PackageNotFoundError
__all__ = []



try:
    __version__ = version("phrosty")
except PackageNotFoundError:
    # package is not installed
    pass


# These were in Lauren's original __init__.py, but I don't think we want them--
#   I think we want the subpackages.
# from .photometry import *
# from .plotting import *
# from .utils import *
# from .imagesubtraction import *
# from .nn2 import *
# # This is broken:
# # from .sourcedetection import *
