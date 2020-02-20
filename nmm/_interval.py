from __future__ import annotations

from dataclasses import dataclass
from typing import Type


@dataclass
class Interval:
    start: int
    stop: int

    @classmethod
    def from_slice(cls: Type[Interval], sl: slice) -> Interval:
        if sl.step not in [None, 1]:
            raise ValueError("Slice has to have one-length steps.")
        if sl.start > sl.stop:
            raise ValueError("Start cannot be higher than stop.")
        return cls(start=sl.start, stop=sl.stop)
