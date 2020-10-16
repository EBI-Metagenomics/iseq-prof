from typing import Dict

__all__ = ["StrDB"]


class StrDB:
    __slots__ = ["_strings"]

    def __init__(self):
        self._strings: Dict[int, str] = {}

    def add(self, string: str) -> int:
        key = hash(string)
        self._strings[key] = string
        return key

    def get(self, key: int) -> str:
        return self._strings[key]
