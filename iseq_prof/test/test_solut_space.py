from assertpy import assert_that
from iseq_prof.solut_space import SampleType, SolutSpaceType


def test_solut_space():
    sstype = SolutSpaceType.from_string("prof.drop_dupl")
    assert_that(sstype.sample_type).is_equal_to(SampleType.PROF)
    assert_that(sstype.drop_duplicates).is_equal_to(True)

    sstype = SolutSpaceType.from_string("prof.keep_dupl")
    assert_that(sstype.sample_type).is_equal_to(SampleType.PROF)
    assert_that(sstype.drop_duplicates).is_equal_to(False)

    sstype = SolutSpaceType.from_string("prof_target.drop_dupl")
    assert_that(sstype.sample_type).is_equal_to(SampleType.PROF_TARGET)
    assert_that(sstype.drop_duplicates).is_equal_to(True)

    sstype = SolutSpaceType.from_string("prof_target.keep_dupl")
    assert_that(sstype.sample_type).is_equal_to(SampleType.PROF_TARGET)
    assert_that(sstype.drop_duplicates).is_equal_to(False)

    sstype = SolutSpaceType.from_string("target.drop_dupl")
    assert_that(sstype.sample_type).is_equal_to(SampleType.TARGET)
    assert_that(sstype.drop_duplicates).is_equal_to(True)

    sstype = SolutSpaceType.from_string("target.keep_dupl")
    assert_that(sstype.sample_type).is_equal_to(SampleType.TARGET)
    assert_that(sstype.drop_duplicates).is_equal_to(False)
