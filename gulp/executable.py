import os
import subprocess
import sys
from tempfile import mkdtemp

ROOT_DIR = os.path.dirname(__file__)
__SUFFIX = ".exe" if os.name == "nt" else ""
EXE = os.path.join(ROOT_DIR, f"gulp{__SUFFIX}")
ENCODING_LST = encoding_lst = ["UTF-8", "GB2312", "GBK", "ISO-8859-1", "UTF-16"]


def gulp():
    raise SystemExit(
        subprocess.call(
            [EXE] + sys.argv[1:],
            close_fds=False,
        )
    )


def available() -> bool:
    test_file = os.path.join("fortran", "Examples", "example1.gin")
    test_file = os.path.join(ROOT_DIR, "..", test_file)
    if os.path.exists(EXE) and os.path.isfile(EXE):
        pass
    else:
        return False

    result = subprocess.run(
        [f"{EXE} < {test_file}"],
        capture_output=True,
        cwd=mkdtemp(),
        shell=True,
    )
    if result.returncode != 0:
        return False
    out, has_decode = result.stdout, False
    for encoding in encoding_lst:
        try:
            out = out.decode(encoding=encoding)
            has_decode = True
            break
        except Exception:
            has_decode = False
    if not has_decode:
        raise Exception(f"Cannot decode the output {out} of the subprocess")

    if "Julian Gale" in out and "Version = 6.2.0" in out:  # type: ignore
        return True
    else:
        return False


if __name__ == "__main__":
    assert available()
