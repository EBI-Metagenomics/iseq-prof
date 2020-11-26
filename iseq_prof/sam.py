from pathlib import Path
from typing import Dict, Optional, Union

from gff_io.interval import PyInterval
from numpy import asarray
from numpy import round as npy_round
from pandas import Series
from pysam import AlignedSegment, AlignmentFile

__all__ = ["SAMMap", "SAMMapItem"]


class SAMMapItem:
    def __init__(self, aseq: AlignedSegment):
        # 0-start coordinate
        reference = asarray(aseq.get_reference_positions(True), float)

        # 0-coordinate, half-open interval
        clip_start = aseq.query_alignment_start
        clip_end = aseq.query_alignment_end
        self._clip = PyInterval(clip_start, clip_end)

        reference = reference[self._clip.to_slice()]
        reference[:] = Series(reference).interpolate().values

        self._reference = asarray(npy_round(reference), int)

    def back_to_query(self, interval: PyInterval) -> Optional[PyInterval]:
        """
        Parameters
        ----------
        interval
            Python interval.

        Returns
        -------
        PyInterval
            Python interval mapped back to query.
        """
        interval = interval.offset(-self._clip.start)

        if interval.start >= len(self._reference):
            return None

        if interval.end <= 0:
            return None

        from_pos = max(interval.start, 0)
        to_pos = min(interval.end, len(self._reference)) - 1

        start = self._reference[from_pos]
        end = self._reference[to_pos] + 1

        return PyInterval(start, end)

    def full_query_interval(self) -> PyInterval:
        return PyInterval(self._reference[0], self._reference[-1] + 1)


class SAMMap:
    def __init__(self, filepath: Union[str, Path]):
        self._items: Dict[str, SAMMapItem] = {}

        filepath = Path(filepath)

        with AlignmentFile(filepath, "rb") as samfile:

            for aseq in samfile.fetch():
                if aseq.flag != 0:
                    continue

                assert aseq.query_name not in self._items
                self._items[aseq.query_name] = SAMMapItem(aseq)

    def back_to_query(
        self, query_name: str, interval: PyInterval
    ) -> Optional[PyInterval]:
        """
        Parameters
        ----------
        query_name
            Query name.
        interval
            Python interval.

        Returns
        -------
        PyInterval
            Python interval mapped back to query.
        """
        return self._items[query_name].back_to_query(interval)

    def full_query_interval(self, query_name: str) -> PyInterval:
        return self._items[query_name].full_query_interval()
