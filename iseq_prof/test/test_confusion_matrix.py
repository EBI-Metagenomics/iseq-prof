from math import nan
from pathlib import Path

from assertpy import assert_that
from iseq_prof import ConfusionMatrix
from numpy import argsort, errstate
from numpy.random import RandomState

TOL = 1e-7


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

    TOL = 1e-7
    TPR = cm.TP / cm.P
    assert_that(TPR[0]).is_close_to(0.0, TOL)
    assert_that(TPR[460]).is_close_to(0.005, TOL)
    assert_that(TPR[-1]).is_close_to(0.006, TOL)
    tmp = cm.TP / (cm.TP + cm.FN)
    for i in range(cm.tpr.shape[0]):
        assert_that(TPR[i]).is_close_to(tmp[i], TOL)

    TNR = cm.TN / cm.N
    assert_that(TNR[0]).is_close_to(1.0, TOL)
    assert_that(TNR[460]).is_close_to(0.997806880130334, TOL)
    assert_that(TNR[-1]).is_close_to(0.9976188984272197, TOL)
    tmp = cm.TN / (cm.TN + cm.FP)
    for i in range(cm.tnr.shape[0]):
        assert_that(TNR[i]).is_close_to(tmp[i], TOL)

    with errstate(divide="ignore", invalid="ignore"):
        PPV = cm.TP / (cm.TP + cm.FP)
    assert_that(PPV[0]).is_close_to(nan, TOL)
    assert_that(PPV[460]).is_close_to(0.010869565217391304, TOL)
    assert_that(PPV[-1]).is_close_to(0.012, TOL)

    NPV = cm.TN / (cm.TN + cm.FN)
    assert_that(NPV[0]).is_close_to(0.9952030777053442, TOL)
    assert_that(NPV[460]).is_close_to(0.995216507136779, TOL)
    assert_that(NPV[-1]).is_close_to(0.9952203955435237, TOL)

    FNR = cm.FN / cm.P
    assert_that(FNR[0]).is_close_to(1.0, TOL)
    assert_that(FNR[460]).is_close_to(0.995, TOL)
    assert_that(FNR[-1]).is_close_to(0.994, TOL)
    tmp = cm.FN / (cm.FN + cm.TP)
    for i in range(cm.fnr.shape[0]):
        assert_that(FNR[i]).is_close_to(tmp[i], TOL)

    FPR = cm.FP / cm.N
    assert_that(FPR[0]).is_close_to(0.0, TOL)
    assert_that(FPR[460]).is_close_to(0.002193119869666019, TOL)
    assert_that(FPR[-1]).is_close_to(0.0023811015727802495, TOL)
    tmp = cm.FP / (cm.FP + cm.TN)
    for i in range(cm.fpr.shape[0]):
        assert_that(FPR[i]).is_close_to(tmp[i], TOL)

    with errstate(divide="ignore", invalid="ignore"):
        FDR = cm.FP / (cm.FP + cm.TP)
    assert_that(FDR[0]).is_close_to(nan, TOL)
    assert_that(FDR[460]).is_close_to(0.9891304347826086, TOL)
    assert_that(FDR[-1]).is_close_to(0.988, TOL)

    FOR = cm.FN / (cm.FN + cm.TN)
    assert_that(FOR[0]).is_close_to(0.004796922294655749, TOL)
    assert_that(FOR[460]).is_close_to(0.004783492863220949, TOL)
    assert_that(FOR[-1]).is_close_to(0.004779604456476268, TOL)

    ACC = (cm.TP + cm.TN) / (cm.P + cm.N)
    assert_that(ACC[0]).is_close_to(0.9952030777053442, TOL)
    assert_that(ACC[460]).is_close_to(0.9930444626727492, TOL)
    assert_that(ACC[-1]).is_close_to(0.9928621796255522, TOL)

    for i in range(len(cm.sensitivity)):
        assert_that(TPR[i]).is_close_to(cm.sensitivity[i], TOL)
        assert_that(TPR[i]).is_close_to(cm.recall[i], TOL)
        assert_that(TPR[i]).is_close_to(cm.tpr[i], TOL)

    for i in range(len(cm.specificity)):
        assert_that(TNR[i]).is_close_to(cm.specificity[i], TOL)
        assert_that(TNR[i]).is_close_to(cm.selectivity[i], TOL)
        assert_that(TNR[i]).is_close_to(cm.tnr[i], TOL)

    for i in range(len(cm.precision)):
        assert_that(PPV[i]).is_close_to(cm.precision[i], TOL)
        assert_that(PPV[i]).is_close_to(cm.ppv[i], TOL)
        assert_that(TNR[i]).is_close_to(cm.tnr[i], TOL)

    for i in range(len(cm.npv)):
        assert_that(NPV[i]).is_close_to(cm.npv[i], TOL)

    for i in range(len(cm.fnr)):
        assert_that(FNR[i]).is_close_to(cm.miss_rate[i], TOL)
        assert_that(FNR[i]).is_close_to(cm.fnr[i], TOL)

    for i in range(len(cm.fpr)):
        assert_that(FPR[i]).is_close_to(cm.fallout[i], TOL)
        assert_that(FPR[i]).is_close_to(cm.fpr[i], TOL)

    for i in range(len(cm.fdr)):
        assert_that(FDR[i]).is_close_to(cm.fdr[i], TOL)

    for i in range(len(cm.for_)):
        assert_that(FOR[i]).is_close_to(cm.for_[i], TOL)

    for i in range(len(cm.accuracy)):
        assert_that(ACC[i]).is_close_to(cm.accuracy[i], TOL)


def test_confusion_matrix_pr_curve():
    random = RandomState(8)

    ntrues = 1000
    nfalses = 207467

    true_samples = random.choice(ntrues + nfalses, ntrues, False)

    nsamples = 500
    samples = random.choice(ntrues + nfalses, nsamples, False)
    scores = random.randn(nsamples)
    idx = argsort(scores)

    cm = ConfusionMatrix(true_samples, nfalses, samples[idx])
    pr = cm.pr_curve

    assert_that(pr.recall[420]).is_close_to(0.005, TOL)
    assert_that(pr.precision[420]).is_close_to(0.011876484560570071, TOL)
    assert_that(pr.auc).is_close_to(4.98054484430107e-05, TOL)


def test_confusion_matrix_roc_curve():
    random = RandomState(8)

    ntrues = 1000
    nfalses = 207467

    true_samples = random.choice(ntrues + nfalses, ntrues, False)

    nsamples = 500
    samples = random.choice(ntrues + nfalses, nsamples, False)
    scores = random.randn(nsamples)
    idx = argsort(scores)

    cm = ConfusionMatrix(true_samples, nfalses, samples[idx])
    roc = cm.roc_curve

    assert_that(roc.fpr[420]).is_close_to(0.00200031812288215, TOL)
    assert_that(roc.tpr[420]).is_close_to(0.005, TOL)
    assert_that(roc.auc).is_close_to(4.762203145560528e-06, TOL)


def test_confusion_matrix_write_read(tmp_path: Path):

    random = RandomState(2)

    ntrues = 100
    nfalses = 100

    true_samples = random.choice(ntrues + nfalses, ntrues, False)

    nsamples = 100
    samples = random.choice(ntrues + nfalses, nsamples, False)
    scores = random.randn(nsamples)
    idx = argsort(scores)

    pr_auc = 0.22942490919917213
    roc_auc = 0.1277

    cm = ConfusionMatrix(true_samples, nfalses, samples[idx])
    assert_that(cm.pr_curve.auc).is_close_to(pr_auc, TOL)
    assert_that(cm.roc_curve.auc).is_close_to(roc_auc, TOL)

    cm.write_pickle(tmp_path / "cm.pkl")
    cm = ConfusionMatrix.read_pickle(tmp_path / "cm.pkl")
    assert_that(cm.pr_curve.auc).is_close_to(pr_auc, TOL)
    assert_that(cm.roc_curve.auc).is_close_to(roc_auc, TOL)
