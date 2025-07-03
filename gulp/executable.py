import os
import subprocess
import sys
try:
    from importlib import metadata
except ImportError: # for Python<3.8
    import importlib_metadata as metadata
ROOT_DIR = os.path.dirname(__file__)


def _program(name, args):
    return subprocess.call([os.path.join(ROOT_DIR, name)] + args, close_fds=False)


def gulp():
    suffix = '.exe' if os.name == 'nt' else ''
    raise SystemExit(_program('gulp' + suffix, sys.argv[1:]))