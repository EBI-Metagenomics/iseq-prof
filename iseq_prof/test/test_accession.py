import pytest
from assertpy import assert_that
from iseq_prof import Accession


@pytest.mark.net
def test_accession():

    acc = Accession("AE014075.1")
    assert_that(acc.name).is_equal_to("AE014075.1")
    assert_that(acc.taxonomy).is_equal_to(
        [
            "Bacteria",
            "Proteobacteria",
            "Gammaproteobacteria",
            "Enterobacterales",
            "Enterobacteriaceae",
            "Escherichia",
        ]
    )
    assert_that(acc.order).is_equal_to("Enterobacterales")
    assert_that(acc.domain).is_equal_to("Bacteria")
    assert_that(acc.organism).is_equal_to("Escherichia coli CFT073")
    assert_that(acc.phylum).is_equal_to("Proteobacteria")
    assert_that(acc.class_).is_equal_to("Gammaproteobacteria")

    assert_that(acc.molecule).is_equal_to("DNA")
    assert_that(str(acc)).is_equal_to("<AE014075.1>")
    assert_that(repr(acc)).is_equal_to('Accession("AE014075.1")')
