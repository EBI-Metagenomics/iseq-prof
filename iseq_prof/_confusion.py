from __future__ import annotations

import pickle
from pathlib import Path
from typing import Iterable, Optional, Union

import numba
from numpy import asarray, empty, linspace, searchsorted, trapz, zeros_like

__all__ = ["ConfusionMatrix", "ROCCurve", "PRCurve"]


class ConfusionMatrix:
    """
    Confusion matrix.

    Parameters
    ----------
    true_samples
        Set of all positive samples from the solution space.
    N
        Number of negative samples from the solution space.
    sorted_samples
        Samples sorted from the most to the least likely one to be considered positive.
        Those samples are usually scored by a method being assessed.
    sample_scores
        Score of each sample from `sorted_samples` in sequence. It therefore has to
        be a monotonically non-decreasing sequence of numbers.
    """

    def __init__(
        self,
        true_samples: Iterable[int],
        N: int,
        sorted_samples: Iterable[int],
        sample_scores: Optional[Iterable[float]] = None,
    ):
        if len(set(sorted_samples) - set(true_samples)) > N:
            raise ValueError("Invalid number of negative samples.")

        true_arr = asarray(true_samples, int)
        P = len(true_arr)

        sorted_arr = asarray(sorted_samples, int)
        self._num_sorted_samples = len(sorted_arr)

        self._TP = empty(len(sorted_arr) + 1, int)
        self._FP = empty(len(sorted_arr) + 1, int)

        self._N = N
        self._P = P
        self._set_tp_fp(true_arr, sorted_arr)
        if sample_scores is None:
            sample_scores = linspace(0, 1, len(sorted_arr))
        self._sample_scores = asarray(sample_scores, float)

    @property
    def sample_scores(self):
        return self._sample_scores

    def cutpoint(self, score: float) -> int:
        return searchsorted(self._sample_scores, score, side="right")

    def _set_tp_fp(self, true_samples, sorted_samples):
        true_samples.sort()
        ins_pos = searchsorted(true_samples, sorted_samples)
        set_tp_fp(self._TP, self._FP, ins_pos, true_samples, sorted_samples)

    def write_pickle(self, filepath: Union[str, Path]):
        with open(Path(filepath), "wb") as file:
            pickle.dump(self, file)

    @staticmethod
    def read_pickle(filepath: Union[str, Path]):
        with open(Path(filepath), "rb") as file:
            return pickle.load(file)

    @property
    def P(self) -> int:
        """
        Number of positive samples.
        """
        return self._P

    @property
    def N(self) -> int:
        """
        Number of negative samples.
        """
        return self._N

    @property
    def TP(self):
        """
        Number of true positive.
        """
        return self._TP

    @property
    def FP(self):
        """
        Number of false positive.
        """
        return self._FP

    @property
    def TN(self):
        """
        Number of true negative.
        """
        return self._N - self.FP

    @property
    def FN(self):
        """
        Number of false negative.
        """
        return self._P - self.TP

    @property
    def sensitivity(self):
        """
        Sensitivity.
        """
        if self._P == 0:
            return zeros_like(self.TP)
        return self.TP / self._P

    @property
    def tpr(self):
        """
        True positive rate.
        """
        return self.sensitivity

    @property
    def recall(self):
        """
        Recall.
        """
        return self.sensitivity

    @property
    def specificity(self):
        """
        Specificity.
        """
        return self.TN / self._N

    @property
    def selectivity(self):
        """
        Selectivity.
        """
        return self.specificity

    @property
    def tnr(self):
        """
        True negative rate.
        """
        return self.specificity

    @property
    def precision(self):
        """
        Precision.
        """
        from numpy import empty, nan

        r = empty(self._num_sorted_samples + 1)
        r[0] = nan
        r[1:] = self.TP[1:] / (self.TP[1:] + self.FP[1:])

        return r

    @property
    def ppv(self):
        """
        Positive predictive value.
        """
        return self.precision

    @property
    def npv(self):
        """
        Negative predictive value.
        """
        from numpy import empty, nan

        r = empty(self._num_sorted_samples + 1)
        r[-1] = nan
        r[:-1] = self.TN[:-1] / (self.TN[:-1] + self.FN[:-1])

        return r

    @property
    def fallout(self):
        """
        Fall-out.
        """
        return 1 - self.specificity

    @property
    def fpr(self):
        """
        False positive rate.
        """
        return self.fallout

    @property
    def fnr(self):
        """
        False negative rate.
        """
        return 1 - self.sensitivity

    @property
    def miss_rate(self):
        """
        Miss rate.
        """
        return self.fnr

    @property
    def for_(self):
        """
        False omission rate.
        """
        return 1 - self.npv

    @property
    def fdr(self):
        """
        False discovery rate.
        """
        return 1 - self.precision

    @property
    def accuracy(self):
        """
        Accuracy.
        """
        return (self.TP + self.TN) / (self._N + self._P)

    @property
    def f1score(self):
        """
        F1 score.
        """
        return 2 * self.TP / (2 * self.TP + self.FP + self.FN)

    @property
    def roc_curve(self) -> ROCCurve:
        from numpy import argsort

        if self._num_sorted_samples < 1:
            raise ValueError("Not enough sorted samples.")

        idx = argsort(self.fpr, kind="stable")
        return ROCCurve(self.fpr[idx], self.tpr[idx])

    @property
    def pr_curve(self) -> PRCurve:
        from numpy import argsort

        if self._num_sorted_samples < 1:
            raise ValueError("Not enough sorted samples.")

        idx = argsort(self.recall, kind="stable")
        return PRCurve(self.recall[idx], self.precision[idx])


T = numba.int_[:]


@numba.njit(numba.void(T, T, T, T, T), cache=True)
def set_tp_fp(TP, FP, ins_pos, true_samples, sorted_samples):
    TP[0] = 0
    FP[0] = 0
    i = 0
    while i < sorted_samples.shape[0]:
        FP[i + 1] = FP[i]
        TP[i + 1] = TP[i]

        j = ins_pos[i]
        if j == len(true_samples) or true_samples[j] != sorted_samples[i]:
            FP[i + 1] += 1
        else:
            TP[i + 1] += 1
        i += 1


class PRCurve:
    """
    Precision-Recall curve.
    """

    def __init__(self, recall: Iterable[float], precision: Iterable[float]):
        self._recall = asarray(recall, float)[1:]
        self._precision = asarray(precision, float)[1:]

    @property
    def recall(self):
        return self._recall

    @property
    def precision(self):
        return self._precision

    @property
    def auc(self) -> float:
        return trapz(self.precision, x=self.recall)


class ROCCurve:
    """
    ROC curve.
    """

    def __init__(self, fpr: Iterable[float], tpr: Iterable[float]):
        self._fpr = asarray(fpr, float)
        self._tpr = asarray(tpr, float)

    @property
    def fpr(self):
        return self._fpr

    @property
    def tpr(self):
        return self._tpr

    @property
    def auc(self) -> float:
        return trapz(self.tpr, x=self.fpr)
