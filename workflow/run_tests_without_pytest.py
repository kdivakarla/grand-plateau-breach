#!/usr/bin/env python
"""Run the test suite without pytest installed.

STOPGAP. ``pytest`` is the supported runner and is declared in environment.yml,
but the current ``gpgn-318`` environment does not have it and installing into
it needs owner approval (CLAUDE.md Section 2.2). This script executes the real
test functions in tests/ by supplying the small part of the pytest API they use
(``approx``, ``mark.parametrize``, ``raises``, ``skip``). It does not
re-implement any assertion.

Once pytest is available, use it instead:

    pytest -q

Usage:
    python workflow/run_tests_without_pytest.py [substring-filter]
"""

from __future__ import annotations

import importlib.util
import sys
import traceback
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


class _Skipped(Exception):
    pass


class _Approx:
    # Make numpy hand `ndarray == _Approx` to us rather than comparing elementwise,
    # which is what pytest.approx does for array comparisons.
    __array_ufunc__ = None
    __array_priority__ = 1000

    def __init__(self, expected, rel=None, abs=None):
        self.expected, self.rel, self.abs = expected, rel, abs

    def _close(self, a, b) -> bool:
        tol = self.abs if self.abs is not None else 1e-6
        if self.rel is not None:
            tol = max(tol, self.rel * abs(b))
        return abs(a - b) <= tol

    def __eq__(self, other):
        import numpy as np

        e = np.asarray(self.expected, dtype=float)
        o = np.asarray(other, dtype=float)
        if e.shape != o.shape:
            return False
        if e.ndim == 0:
            return self._close(float(o), float(e))
        return all(self._close(float(x), float(y)) for x, y in zip(o.ravel(), e.ravel()))

    def __ne__(self, other):
        return not self.__eq__(other)


def _build_pytest_shim() -> types.ModuleType:
    mod = types.ModuleType("pytest")
    mod.approx = lambda expected, rel=None, abs=None: _Approx(expected, rel, abs)

    def skip(reason=""):
        raise _Skipped(reason)

    class _Raises:
        def __init__(self, exc):
            self.exc = exc

        def __enter__(self):
            return self

        def __exit__(self, t, v, tb):
            if t is None:
                raise AssertionError(f"expected {self.exc.__name__}")
            return issubclass(t, self.exc)

    def parametrize(argnames, argvalues):
        names = [a.strip() for a in argnames.split(",")]

        def deco(fn):
            fn._params = (names, list(argvalues))
            return fn

        return deco

    mod.skip = skip
    mod.raises = lambda exc: _Raises(exc)
    mod.mark = types.SimpleNamespace(parametrize=parametrize)
    mod.importorskip = lambda name, **k: __import__(name)
    return mod


sys.modules.setdefault("pytest", _build_pytest_shim())


def _load(path: Path) -> types.ModuleType:
    spec = importlib.util.spec_from_file_location(f"t_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv: list[str]) -> int:
    keep = argv[1] if len(argv) > 1 else ""
    files = sorted((ROOT / "tests").rglob("test_*.py"))
    passed = failed = skipped = 0
    failures: list[tuple[str, str]] = []

    for path in files:
        module = _load(path)
        rel = path.relative_to(ROOT)
        for name in sorted(vars(module)):
            if not name.startswith("test_"):
                continue
            fn = getattr(module, name)
            if not callable(fn):
                continue
            params = getattr(fn, "_params", None)
            calls = (
                [((), {})] if params is None
                else [((), dict(zip(params[0], v if isinstance(v, tuple) else (v,))))
                      for v in params[1]]
            )
            for args, kwargs in calls:
                label = f"{rel}::{name}" + (f"[{list(kwargs.values())}]" if kwargs else "")
                if keep and keep not in label:
                    continue
                try:
                    fn(*args, **kwargs)
                except _Skipped as exc:
                    skipped += 1
                    print(f"SKIP {label}: {exc}")
                except Exception:
                    failed += 1
                    failures.append((label, traceback.format_exc()))
                    print(f"FAIL {label}")
                else:
                    passed += 1
                    print(f"pass {label}")

    for label, tb in failures:
        print(f"\n{'=' * 70}\nFAILED {label}\n{tb}")
    print(f"\n{passed} passed, {failed} failed, {skipped} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
