import logging
import numpy as np
import scipy.stats
from unittest import TestCase
import distributions

logger=logging.getLogger("test_distributions")


def fractional_error(a, b):
    return np.abs((a-b)/a)

def check_fractional_error(a, b, tolerance, message):
    if fractional_error(a, b)>tolerance:
        logger.error("Fractional error of {0} too large. Expected "+
            "{1} but found {2}".format(message, a, b))
        return False
    return True


class TestExponential(TestCase):
    def test_average(self):
        """
        This tests the theoretical mean and standard deviation.
        """
        lam=0.5
        rng=np.random.RandomState()
        ed=distributions.ExponentialDistribution(lam, 0.0)
        cnt=10000
        res=np.zeros(cnt)
        for i in range(cnt):
            res[i]=ed.sample(0, rng)
        lambda_estimator=1/np.average(res)
        too_low=lam < lambda_estimator*(1-1.96/np.sqrt(cnt))
        too_high=lam > lambda_estimator*(1+1.96/np.sqrt(cnt))
        self.assertTrue(not too_low)
        self.assertTrue(not too_high)
        variance=np.var(res)
        check_fractional_error(variance, np.power(lam, -2), 0.01, "variance")

    def test_integrals(self):
        """
        Are the hazard integral and its inverse really inverses?
        """
        tol=0.001
        ed=ExponentialDistribution(0.5, 3)
        for x in np.linspace(3, 7, num=5):
            xa=ed.hazard_integral(3, x)
            self.assertTrue(tol>abs(xa-ed.implicit_hazard_integral(xa, 3)))


class TestWeibull(TestCase):
    def test_samples(self):
        pass


