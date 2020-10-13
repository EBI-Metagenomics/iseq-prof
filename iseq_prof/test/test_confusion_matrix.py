from math import nan

from assertpy import assert_that
from iseq_prof import ConfusionMatrix
from numpy import argsort, errstate
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

    tol = 1e-7
    TPR = cm.TP / cm.P
    assert_that(TPR[0]).is_close_to(0.0, tol)
    assert_that(TPR[460]).is_close_to(0.005, tol)
    assert_that(TPR[-1]).is_close_to(0.006, tol)

    TNR = cm.TN / cm.N
    assert_that(TNR[0]).is_close_to(1.0, tol)
    assert_that(TNR[460]).is_close_to(0.997806880130334, tol)
    assert_that(TNR[-1]).is_close_to(0.9976188984272197, tol)

    with errstate(divide="ignore", invalid="ignore"):
        PPV = cm.TP / (cm.TP + cm.FP)
    assert_that(PPV[0]).is_close_to(nan, tol)
    assert_that(PPV[460]).is_close_to(0.010869565217391304, tol)
    assert_that(PPV[-1]).is_close_to(0.012, tol)

    NPV = cm.TN / (cm.TN + cm.FN)
    assert_that(NPV[0]).is_close_to(0.9952030777053442, tol)
    assert_that(NPV[460]).is_close_to(0.995216507136779, tol)
    assert_that(NPV[-1]).is_close_to(0.9952203955435237, tol)

    FNR = cm.FN / cm.P
    assert_that(FNR[0]).is_close_to(1.0, tol)
    assert_that(FNR[460]).is_close_to(0.995, tol)
    assert_that(FNR[-1]).is_close_to(0.994, tol)

    FPR = cm.FP / cm.N
    assert_that(FPR[0]).is_close_to(0.0, tol)
    assert_that(FPR[460]).is_close_to(0.002193119869666019, tol)
    assert_that(FPR[-1]).is_close_to(0.0023811015727802495, tol)

    with errstate(divide="ignore", invalid="ignore"):
        FDR = cm.FP / (cm.FP + cm.TP)
    assert_that(FDR[0]).is_close_to(nan, tol)
    assert_that(FDR[460]).is_close_to(0.9891304347826086, tol)
    assert_that(FDR[-1]).is_close_to(0.988, tol)

    FOR = cm.FN / (cm.FN + cm.TN)
    assert_that(FOR[0]).is_close_to(0.004796922294655749, tol)
    assert_that(FOR[460]).is_close_to(0.004783492863220949, tol)
    assert_that(FOR[-1]).is_close_to(0.004779604456476268, tol)

    ACC = (cm.TP + cm.TN) / (cm.P + cm.N)
    assert_that(ACC[0]).is_close_to(0.9952030777053442, tol)
    assert_that(ACC[460]).is_close_to(0.9930444626727492, tol)
    assert_that(ACC[-1]).is_close_to(0.9928621796255522, tol)

    for i in range(len(cm.sensitivity)):
        assert_that(TPR[i]).is_close_to(cm.sensitivity[i], tol)
        assert_that(TPR[i]).is_close_to(cm.recall[i], tol)
        assert_that(TPR[i]).is_close_to(cm.tpr[i], tol)

    for i in range(len(cm.specificity)):
        assert_that(TNR[i]).is_close_to(cm.specificity[i], tol)
        assert_that(TNR[i]).is_close_to(cm.selectivity[i], tol)
        assert_that(TNR[i]).is_close_to(cm.tnr[i], tol)

    for i in range(len(cm.precision)):
        assert_that(PPV[i]).is_close_to(cm.precision[i], tol)
        assert_that(PPV[i]).is_close_to(cm.ppv[i], tol)
        assert_that(TNR[i]).is_close_to(cm.tnr[i], tol)

    for i in range(len(cm.npv)):
        assert_that(NPV[i]).is_close_to(cm.npv[i], tol)

    for i in range(len(cm.fnr)):
        assert_that(FNR[i]).is_close_to(cm.miss_rate[i], tol)
        assert_that(FNR[i]).is_close_to(cm.fnr[i], tol)

    for i in range(len(cm.fpr)):
        assert_that(FPR[i]).is_close_to(cm.fallout[i], tol)
        assert_that(FPR[i]).is_close_to(cm.fpr[i], tol)

    for i in range(len(cm.fdr)):
        assert_that(FDR[i]).is_close_to(cm.fdr[i], tol)

    for i in range(len(cm.for_)):
        assert_that(FOR[i]).is_close_to(cm.for_[i], tol)

    for i in range(len(cm.accuracy)):
        assert_that(ACC[i]).is_close_to(cm.accuracy[i], tol)
