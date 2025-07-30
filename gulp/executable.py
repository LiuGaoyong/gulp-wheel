import os
import subprocess
import sys
from tempfile import TemporaryDirectory

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


def available(verbose: bool = False) -> bool:
    with TemporaryDirectory() as tmpdir:
        test_file = os.path.join(tmpdir, "example1.gin")
        with open(test_file, "w") as f:
            f.write("""opti prop conp
title
alumina test file
end
cell
4.7602   4.7602  12.9933  90.000000  90.000000 120.0
frac
Al core 0.000000   0.000000   0.352160
Al shel 0.000000   0.000000   0.352160
O  core 0.306240   0.000000   0.250000
O  shel 0.306240   0.000000   0.250000
space
167
species
Al core  0.043
Al shel  2.957
O  core  0.513
O  shel -2.513
buckingham
Al shel O shel  2409.505 0.2649  0.00 0.0 10.0
O  shel O shel    25.410 0.6937 32.32 0.0 12.0
spring
Al 403.98
O   20.53
output xr example1
output marvin example1.mvn
""")
        if verbose:
            print(f"The GULP program: {EXE}")
            print(f"Test GULP program with file: {test_file}")
        if os.path.exists(EXE) and os.path.isfile(EXE):
            pass
        else:
            return False
        result = subprocess.run(
            [f"{EXE} < {test_file}"],
            capture_output=True,
            cwd=tmpdir,
            shell=True,
        )

    if verbose:
        print(f"Subprocess result: {result}")
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
    if verbose:
        print(out)

    if "Julian Gale" in out and "Version = 6.2.0" in out:  # type: ignore
        return True
    else:
        return False


if __name__ == "__main__":
    assert available(verbose=True)
