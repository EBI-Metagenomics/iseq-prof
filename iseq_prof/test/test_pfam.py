import hashlib
import os
from pathlib import Path

from assertpy import assert_that
from iseq_prof import pfam


def test_pfam_clan_default():
    clans = pfam.Clans()
    clan_ids = list(clans)

    assert_that(len(clan_ids)).is_equal_to(634)
    assert_that(clan_ids).contains("Globin")
    assert_that(clan_ids).contains("PHM_PNGase_F")

    assert_that(clans.get("PF07823")).is_equal_to("2H")
    assert_that(clans.get("PF07823.1")).is_equal_to("2H")
    assert_that(clans.get("PF07823.2")).is_equal_to("2H")
    assert_that(clans.get("ABCX")).is_equal_to(None)

    assert_that(clans.get("PF15446")).is_equal_to("zf-FYVE-PHD")


def test_pfam_latest_version():
    version = pfam.latest_version()
    v = version.split(".")
    a = int(v[0])
    b = int(v[1])
    assert_that(a).is_greater_than_or_equal_to(33)
    assert_that(b).is_greater_than_or_equal_to(1)


def test_pfam_clan_latest():
    clans = pfam.Clans()
    clan_ids = list(clans)

    default_clans = pfam.Clans()
    default_set = set(list(default_clans))
    frac = len(set(clan_ids) & default_set) / len(set(clan_ids) | default_set)

    assert_that(frac).is_greater_than(0.7)


def test_pfam_make_csv(tmp_path: Path):
    os.chdir(tmp_path)
    filepath = pfam.make_csv("33.1")
    expected = "35d6380008d91d124b81865330098b18a42c6f01d1a5ceba4e5f2fe61128aa84"
    with open(filepath, "rb") as f:
        got = hashlib.sha256(f.read()).hexdigest()
    assert_that(got).is_equal_to(expected)
