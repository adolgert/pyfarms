import logging
import numpy as np
import scipy.stats
import scipy
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
        have_a_test=False
        self.assertTrue(have_a_test)




class TestGamma(TestCase):
    def test_theoretical(self):
        """
        Test theoretical Gamma.
        """
        alpha=1.34
        beta=0.18
        theta=1.0/beta
        te=0.0
        now=0.0
        tol=0.0001
        rng=np.random.RandomState()
        ed=distributions.GammaDistribution(alpha, beta, te)
        cnt=10000
        res=np.zeros(cnt)
        for i in range(cnt):
            res[i]=ed.sample(now, rng)
        avg=np.average(res)
        # The calculation of the average varies a lot, but that makes
        # sense for a Gamma. What's a better statistic of the samples?
        logger.debug("Gamma avg {0} theory {1}".format(avg, alpha*theta))
        self.assertTrue(0.01>abs(avg-alpha*theta))
        variance=np.var(res)
        logger.debug("Gamma variance {0} theory {1}".format(
            var, alpha*theta**2))
        self.assert_true(tol>abs(var-alpha*theta**2))
        skew=scipy.stats.skew(res)
        logger.debug("Gamma skew {0} theory {1}".format(
            skew, 2/np.sqrt(alpha)))
        self.assert_true(tol>abs(skew-2/np.sqrt(alpha)))


    def test_integrals(self):
        """
        Are the Gamma hazard integral and its inverse really inverses?
        """
        alpha=1.34
        beta=0.18
        theta=1.0/beta
        te=0.2
        now=1.1
        tol=0.0001
        rng=np.random.RandomState()
        ed=distributions.GammaDistribution(alpha, beta, te)
        for x in np.linspace(now, 2*now, num=5):
            xa=ed.hazard_integral(now, x)
            self.assertTrue(tol>abs(x-ed.implicit_hazard_integral(xa, now)))

    def test_samples(self):
        """
        Sample from the Gamma scipy distribution and from ours. Compare.
        """
        cnt=10000
        rng=np.random.RandomState()
        samples=np.zeros(cnt)
        alpha=1.34
        beta=0.18
        theta=1.0/beta
        te=0.2
        now=1.1
        exp_dist=distributions.GammaDistribution(alpha, beta, te)
        for i in range(cnt):
            samples[i]=exp_dist.sample(now, rng)
        emp_dist=distributions.EmpiricalDistribution(samples)
        system=np.zeros(cnt)
        for i in range(cnt):
            v=now-1
            while v<now:
                v=scipy.stats.gamma.rvs(a=alpha, scale=1.0/beta, loc=0, size=1)
            system[i]=v
        system_dist=distributions.EmpiricalDistribution(system)
        ks_fit=emp_dist.compare_empirical(system_dist)
        logger.debug("Exponential test_samples ks {0}".format(ks_fit))
        self.assertTrue(ks_fit<1.63)


    def test_anderson_samples(self):
        """
        Sample from the Gamma scipy distribution and from ours using Anderson's
        method.
        """
        cnt=10000
        rng=np.random.RandomState()
        alpha=1.34
        beta=0.18
        theta=1.0/beta
        te=0.2
        now=1.1
        exp_dist=distributions.GammaDistribution(alpha, beta, te)
        samples=distributions.anderson_sample_tester(exp_dist, now, cnt, rng)
        emp_dist=distributions.EmpiricalDistribution(samples)
        system=np.zeros(cnt)
        for i in range(cnt):
            v=now-1
            while v<now:
                v=scipy.stats.gamma.rvs(a=alpha, scale=1.0/beta, loc=0, size=1)
            system[i]=v
        system_dist=distributions.EmpiricalDistribution(system)
        ks_fit=emp_dist.compare_empirical(system_dist)
        logger.debug("Exponential test_samples ks {0}".format(ks_fit))
        self.assertTrue(ks_fit<1.63)

