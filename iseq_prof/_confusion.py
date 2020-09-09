from __future__ import annotations

import pickle
from pathlib import Path
from typing import Iterable, Optional

__all__ = ["ConfusionMatrix", "ROCCurve", "PRCurve"]


class ConfusionMatrix:
    """
    Confusion matrix.

    Parameters
    ----------
    true_samples
        Set of all positive samples from the solution space.
    N
        Number of negative samples.
    sorted_samples
        Samples sorted from the most to the least likely one to be considered positive.
    """

    def __init__(
        self,
        true_samples: Iterable[int],
        N: int,
        sorted_samples: Iterable[int],
        sample_scores: Optional[Iterable[float]] = None,
    ):
        from numpy import asarray, empty, linspace

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
        from numpy import searchsorted

        return searchsorted(self._sample_scores, score, side="right")

    def _set_tp_fp(self, true_samples, sorted_samples):
        from numpy import asarray, searchsorted

        true_arr = asarray(true_samples, int)
        true_arr.sort()

        ins_pos = searchsorted(true_arr, sorted_samples)

        self._TP[0] = 0
        self._FP[0] = 0
        i = 0
        while i < len(sorted_samples):
            self._FP[i + 1] = self._FP[i]
            self._TP[i + 1] = self._TP[i]

            j = ins_pos[i]
            if j == len(true_arr) or true_arr[j] != sorted_samples[i]:
                self._FP[i + 1] += 1
            else:
                self._TP[i + 1] += 1
            i += 1

    def write_pickle(self, filepath: Path):
        with open(filepath, "wb") as file:
            pickle.dump(self, file)

    @staticmethod
    def read_pickle(filepath: Path):
        with open(filepath, "rb") as file:
            return pickle.load(file)

    @property
    def P(self) -> int:
        return self._P

    @property
    def N(self) -> int:
        return self._N

    @property
    def TP(self):
        return self._TP

    @property
    def FP(self):
        return self._FP

    @property
    def TN(self):
        return self._N - self.FP

    @property
    def FN(self):
        return self._P - self.TP

    @property
    def sensitivity(self):
        """Sensitivity (also known as true positive rate.)"""
        return self.TP / self._P

    @property
    def tpr(self):
        return self.sensitivity

    @property
    def recall(self):
        return self.sensitivity

    @property
    def specifity(self):
        """Specifity (also known as true negative rate.)"""
        return self.TN / self._N

    @property
    def precision(self):
        from numpy import empty, nan

        r = empty(self._num_sorted_samples + 1)
        r[0] = nan
        r[1:] = self.TP[1:] / (self.TP[1:] + self.FP[1:])

        return r

    @property
    def npv(self):
        """Negative predictive value."""
        from numpy import empty, nan

        r = empty(self._num_sorted_samples + 1)
        r[-1] = nan
        r[:-1] = self.TN[:-1] / (self.TN[:-1] + self.FN[:-1])

        return r

    @property
    def fallout(self):
        """Fall-out (also known as false positive rate.)"""
        return 1 - self.specifity

    @property
    def fpr(self):
        return self.fallout

    @property
    def fnr(self):
        """False negative rate."""
        return 1 - self.sensitivity

    @property
    def fdr(self):
        """False discovery rate."""
        return 1 - self.precision

    @property
    def accuracy(self):
        """Accuracy."""
        return (self.TP + self.TN) / (self._N + self._P)

    @property
    def f1score(self):
        """F1 score (harmonic mean of precision and sensitivity)."""
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


class PRCurve:
    """
    Precision-Recall curve.
    """

    def __init__(self, recall: Iterable[float], precision: Iterable[float]):
        from numpy import asarray

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
        return auc(self.recall, self.precision)

    def plot(self, ax=None):
        xlabel = "recall (sensitivity)"
        ylabel = "precision"
        title = f"Precision-Recall curve (area={self.auc:6.4f})"
        ax = plot(self.recall, self.precision, xlabel, ylabel, title, ax)
        ax.plot([0, 1], [1, 0], linestyle="--")
        ax.legend(loc="lower left")
        return ax


class ROCCurve:
    """
    ROC curve.
    """

    def __init__(self, fpr: Iterable[float], tpr: Iterable[float]):
        from numpy import asarray

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
        return auc(self.fpr, self.tpr)

    def plot(self, ax=None):
        xlabel = "false positive rate"
        ylabel = "true positive rate"
        title = f"ROC curve (area={self.auc:6.4f})"
        ax = plot(self.fpr, self.tpr, xlabel, ylabel, title, ax)
        ax.plot([0, 1], [0, 1], linestyle="--")
        ax.legend(loc="lower right")
        return ax


def plot(x, y, xlabel, ylabel, title, ax):
    import seaborn as sns
    from matplotlib import pyplot as plt

    sns.set(color_codes=True)

    if ax is None:
        ax = plt.subplots()[1]
    ax.plot(x, y, label=title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax


def auc(x, y) -> float:
    left = x[0]
    area = 0.0
    for i in range(1, len(x)):
        width = x[i] - left
        area += width * y[i - 1]
        left = x[i]
    area += (1 - left) * y[-1]
    return area
