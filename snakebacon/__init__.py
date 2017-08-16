from .agedepth import AgeDepthModel
from .records import CalibCurve, ChronRecord, ProxyRecord, DatedProxyRecord
from .records import read_14c, read_chron, read_proxy
from .utils import suggest_accumulation_rate

try:
    from .version import version as __version__
except ImportError:
    raise ImportError('snakebacon not properly installed. If you are running from the source directory, please instead '
                      'create a new virtual environment (using conda or virtualenv) and then install it in-place by '
                      'running: pip install -e .')
