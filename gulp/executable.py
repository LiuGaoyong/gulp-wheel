import os
import subprocess
import sys

ROOT_DIR = os.path.dirname(__file__)


def _program(name, args):
    return subprocess.call(
        [os.path.join(ROOT_DIR, name)] + args,
        close_fds=False,
    )


def gulp():
    suffix = ".exe" if os.name == "nt" else ""
    raise SystemExit(_program("gulp" + suffix, sys.argv[1:]))
