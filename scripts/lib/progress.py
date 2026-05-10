"""Thin tqdm wrappers with consistent formatting."""
from tqdm import tqdm as _tqdm
from typing import Iterable, TypeVar

T = TypeVar("T")

_BAR_FMT = "{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"


def bar(iterable: Iterable[T], desc: str, total: int | None = None, **kw) -> Iterable[T]:
    return _tqdm(iterable, desc=desc, total=total, bar_format=_BAR_FMT,
                 dynamic_ncols=True, **kw)


def counter(desc: str, total: int | None = None, **kw) -> _tqdm:
    """Return a bare tqdm instance for manual .update() calls."""
    return _tqdm(total=total, desc=desc, bar_format=_BAR_FMT,
                 dynamic_ncols=True, **kw)
