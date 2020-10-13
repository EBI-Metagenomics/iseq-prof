from assertpy import assert_that
from iseq_prof import ConfusionMatrix
from numpy import argsort
from numpy.random import RandomState


def test_confusion_matrix():
    random = RandomState(8)

    ntrues = 1000
    nfalses = 207467

    true_samples = random.choice(ntrues + nfalses, ntrues, False)

    nsamples = 500
    samples = random.choice(ntrues + nfalses, nsamples, False)
    scores = random.randn(nsamples)
    idx = argsort(scores)

    cm = ConfusionMatrix(true_samples, nfalses, samples[idx])
    assert_that(cm.P).is_equal_to(ntrues)
    assert_that(cm.N).is_equal_to(nfalses)
    assert_that(cm.TP.shape[0]).is_equal_to(nsamples + 1)
    assert_that(cm.FP.shape[0]).is_equal_to(nsamples + 1)
    assert_that(cm.TN.shape[0]).is_equal_to(nsamples + 1)
    assert_that(cm.FN.shape[0]).is_equal_to(nsamples + 1)
    assert_that(cm.cutpoint(0.6)).is_equal_to(300)
