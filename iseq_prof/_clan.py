from pathlib import Path

import pandas as pd

__all__ = ["Clans"]


class Clans:
    def __init__(self, filepath: Path):
        df = pd.read_csv(filepath)
        self._prof_to_clan = {}
        for row in df.itertuples(False):
            self._prof_to_clan[row.prof_acc] = row.clan_id

    def belongs_to_clan(self, profile_acc: str) -> bool:
        profile_acc = profile_acc.partition(".")[0]
        return profile_acc in self._prof_to_clan

    def get_clan(self, profile_acc: str) -> str:
        profile_acc = profile_acc.partition(".")[0]
        return self._prof_to_clan.get(profile_acc, profile_acc)
