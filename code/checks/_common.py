"""Diagnostic I/O. Scientific input files are read-only; reports use outputs/checks."""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
CODE = ROOT / "code"
DATA = ROOT / "data" / "raw"
OUTPUT = ROOT / "outputs" / "checks"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    path = Path(path).resolve()
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.name


def jsonable(value):
    if hasattr(value, "tolist"):
        return jsonable(value.tolist())
    if isinstance(value, dict):
        return {jsonable(str(key)): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    if isinstance(value, Path):
        return relative(value)
    if isinstance(value, str):
        normalized = value.replace("\\", "/")
        prefix = ROOT.as_posix() + "/"
        if normalized.startswith(prefix):
            return normalized[len(prefix):]
        if re.match(r"^[A-Za-z]:/", normalized):
            return "external-origin/" + normalized.rsplit("/", 1)[-1]
    return value


def output_directory(path):
    path = Path(path).resolve()
    if path == DATA or DATA in path.parents:
        raise ValueError("Diagnostics cannot overwrite shipped data/raw inputs")
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_report(path, value):
    path = Path(path)
    output_directory(path.parent)
    path.write_text(json.dumps(jsonable(value), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return path


@contextmanager
def protect_inputs(paths):
    paths = [Path(path) for path in paths]
    before = {relative(path): digest(path) for path in paths}
    try:
        yield before
    finally:
        after = {relative(path): digest(path) for path in paths}
        if after != before:
            raise RuntimeError("A diagnostic changed a shipped scientific input")
