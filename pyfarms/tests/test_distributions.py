import logging
import numpy as np
import scipy.stats
from unittest import TestCase
import gspn.distributions as distributions

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
        te=0.3
        rng=np.random.RandomState()
        ed=distributions.ExponentialDistribution(lam, te)
        cnt=10000
        res=np.zeros(cnt)
        for i in range(cnt):
            res[i]=ed.sample(te, rng)
        lambda_estimator=1/(np.average(res)-te)
        logger.debug("Exponential estimator {0} lambda {1}".format(
            lambda_estimator, lam))
        too_low=lam < lambda_estimator*(1-1.96/np.sqrt(cnt))
        self.assertTrue(not too_low)
        too_high=lam > lambda_estimator*(1+1.96/np.sqrt(cnt))
        self.assertTrue(not too_high)
        variance=np.var(res-te)
        check_fractional_error(variance, np.power(lam, -2), 0.01, "variance")

    def test_integrals(self):
        """
        Are the hazard integral and its inverse really inverses?
        """
        tol=0.001
        ed=distributions.ExponentialDistribution(0.5, 3)
        for x in np.linspace(3, 7, num=5):
            xa=ed.hazard_integral(3, x)
            self.assertTrue(tol>abs(x-ed.implicit_hazard_integral(xa, 3)))

    def test_samples(self):
        """
        Sample from the scipy distribution and from ours. Compare.
        """
        cnt=10000
        rng=np.random.RandomState()
        samples=np.zeros(cnt)
        lam=2.0
        te=0.7
        now=1.1
        exp_dist=distributions.ExponentialDistribution(lam, te)
        for i in range(cnt):
            samples[i]=exp_dist.sample(now, rng)
        emp_dist=distributions.EmpiricalDistribution(samples)
        system=scipy.stats.expon.rvs(scale=1./lam, loc=now, size=cnt)
        system_dist=distributions.EmpiricalDistribution(system)
        ks_fit=emp_dist.compare_empirical(system_dist)
        logger.debug("Exponential test_samples ks {0}".format(ks_fit))
        self.assertTrue(ks_fit<1.63)


    def test_anderson_samples(self):
        """
        Sample from the scipy distribution and from ours using Anderson's
        method.
        """
        cnt=10000
        rng=np.random.RandomState()
        lam=2.0
        te=0.7
        now=1.1
        exp_dist=distributions.ExponentialDistribution(lam, te)
        samples=distributions.anderson_sample_tester(exp_dist, now, cnt, rng)
        emp_dist=distributions.EmpiricalDistribution(samples)
        system=scipy.stats.expon.rvs(scale=1./lam, loc=now, size=cnt)
        system_dist=distributions.EmpiricalDistribution(system)
        ks_fit=emp_dist.compare_empirical(system_dist)
        logger.debug("Exponential test_samples ks {0}".format(ks_fit))
        self.assertTrue(ks_fit<1.63)


class TestWeibull(TestCase):
    def test_samples(self):
        pass


