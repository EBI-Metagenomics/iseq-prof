# from __future__ import annotations

# from typing import Dict, Iterable, Set

# from .solut_space import Sample, SolutSpace


# class OSample:
#     __slots__ = ["organism_hash", "sample"]

#     def __init__(
#         self,
#         organism_hash: int,
#         sample: Sample,
#     ):
#         self.organism_hash = organism_hash
#         self.sample = sample

#     def __hash__(self) -> int:
#         return hash((self.organism_hash, self.sample))

#     def __eq__(self, you: OSample):  # type: ignore[override]
#         return self.organism_hash == you.organism_hash and self.sample == you.sample


# class DB:
#     __slots__ = ["_organisms", "_profiles", "_targets"]

#     def __init__(self):
#         self._organisms: Dict[int, str] = {}
#         self._profiles: Dict[int, str] = {}
#         self._targets: Dict[int, str] = {}

#     def add_organism(self, organism: str) -> int:
#         key = hash(organism)
#         self._organisms[key] = organism
#         return key

#     def add_profile(self, profile: str) -> int:
#         key = hash(profile)
#         self._profiles[key] = profile
#         return key

#     def add_target(self, target: str) -> int:
#         key = hash(target)
#         self._targets[key] = target
#         return key

#     def create_sample(self, organism: str, profile: str, target: str) -> OSample:
#         ohash = self.add_organism(organism)
#         phash = self.add_profile(profile)
#         thash = self.add_target(target)
#         return OSample(ohash, phash, thash)

#     def cross_create_samples(
#         self, organism: str, profiles: Iterable[str], targets: Iterable[str]
#     ):
#         ohash = self.add_organism(organism)
#         for profile in profiles:
#             phash = self.add_profile(profile)
#             for target in targets:
#                 thash = self.add_target(target)
#                 yield OSample(ohash, phash, thash)

#     def get_organism(self, sample: OSample) -> str:
#         return self._organisms[sample.organism_hash]

#     def get_profile(self, sample: OSample) -> str:
#         return self._profiles[sample.profile_hash]

#     def get_target(self, sample: OSample) -> str:
#         return self._targets[sample.target_hash]


# class OSolutSpace:
#     def __init__(self):

#         self._db = DB()
#         self._sample_space: Set[OSample] = set()
#         self._true_samples: Set[OSample] = set()
#         self._hits: Dict[OSample, float] = {}

#     def add_organism(self, name: str, solut_space: SolutSpace):
#         stype = SolutSpaceType(SampleType.PROF_TARGET, True)
#         sample_space, true_samples, ordered_hits = solut_space._get_samples(stype)

#         for sample in sample_space:
#             osample = self._convert(name, solut_space, sample)
#             self._sample_space.add(osample)

#         for sample in true_samples:
#             osample = self._convert(name, solut_space, sample)
#             self._true_samples.add(osample)

#         for sample in ordered_hits:
#             osample = self._convert(name, solut_space, sample)
#             self._hits[osample] = solut_space.hit_evalue(sample)

#     def hit_evalue(self, hit: OSample) -> float:
#         return self._hits[hit]

#     def _convert(self, name: str, solut_space: SolutSpace, sample: Sample):
#         prof = solut_space.profile(sample)
#         tgt = solut_space.target(sample)
#         return self._db.create_sample(name, prof, tgt)

#     @property
#     def _sorted_hits(self):
#         return [k for k, _ in sorted(self._hits.items(), key=lambda x: x[1])]

#     def per_profile(self, profile: str):
#         return (
#             [s for s in self._sample_space if self._db.get_profile(s) == profile],
#             [s for s in self._true_samples if self._db.get_profile(s) == profile],
#             [s for s in self._sorted_hits if self._db.get_profile(s) == profile],
#         )
