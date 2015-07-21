import logging
import itertools
#from enum import Enum
import numpy as np
import scipy.spatial.distance as distance
import distributions
import llcp
import sample
import runner
import point_process

logger=logging.getLogger(__file__)


##############################################################
# Within-farm disease model
##############################################################
class DiseaseState(object):
    susceptible=1
    latent=2
    subclinical=3
    clinical=4
    recovered=5

class DiseasePlace(object):
    def __init__(self, disease_model):
        self.disease_model=disease_model
        self.state=DiseaseState.susceptible


class DiseaseABTransition:
    def __init__(self, farm_place, a, b, distribution):
        self.place=farm_place
        self.a=a
        self.b=b
        self.dist=distribution

    def depends(self):
        return [self.place]

    def affected(self):
        return [self.place]

    def enabled(self, now):
        if self.place.state==self.a:
            return (True, self.dist(now))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.place.state=self.b


class InfectiousIntensity:
    """
    This is part of a transition which contributes to the hazard
    rate through a factor called an intensity.
    """
    def __init__(self, disease_model):
        self.disease_model=disease_model
        self.place=disease_model.place
    def depends(self):
        return [self.place]
    def intensity(self, now):
        if self.place.state in (DiseaseState.subclinical, DiseaseState.clinical):
            return 1
        return None

class DetectionIntensity(object):
    def __init__(self, disease_model):
        self.disease_model=disease_model
        self.place=disease_model.place
    def depends(self):
        return [self.place]
    def intensity(self, now):
        if self.place.state in (DiseaseState.clinical, DiseaseState.recovered):
            return 1
        return None


class InfectPartial:
    """
    This is the action-part of a transition. It isn't a whole transition.
    """
    def __init__(self, farm):
        self.farm=farm
    def depends(self):
        return [self.farm.place]
    def affected(self):
        return [self.farm.place]
    def enabled(self, now):
        if self.farm.place.state in (DiseaseState.susceptible,):
            return True
        else:
            return False
    def fire(self, now, rng):
        self.farm.place.state=DiseaseState.latent


class SendIntensity(object):
    """
    Whether and how much this farm can send.
    """
    def __init__(self, farm):
        self.farm=farm
    def depends(self):
        return [self.farm.place]
    def intensity(self, now):
        return None

class DiseaseModel(object):
    """
    This is a scenario model for disease state within a farm.
    """
    def __init__(self, farm):
        self.farm=farm
        self.place=DiseasePlace(self)

    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        t=DiseaseABTransition(self.place, DiseaseState.latent, DiseaseState.subclinical,
            lambda now: distributions.ExponentialDistribution(0.5, now))
        writer.add_transition(t)
        t=DiseaseABTransition(self.place, DiseaseState.subclinical, DiseaseState.clinical,
            lambda now: distributions.ExponentialDistribution(0.5, now))
        writer.add_transition(t)
        t=DiseaseABTransition(self.place, DiseaseState.clinical, DiseaseState.recovered,
            lambda now: distributions.ExponentialDistribution(0.5, now))
        writer.add_transition(t)

    def infectious_intensity(self):
        return InfectiousIntensity(self)

    def infection_partial(self):
        return InfectPartial(self)

    def detectable_intensity(self):
        return DetectionIntensity(self)

###############################################################
# Quarantine model
###############################################################
class QuarantinePlace(object):
    def __init__(self):
        self.state=False

class QuarantineTransition(object):
    def __init__(self, model):
        self.model=model
        self.detectable=model.farm.detectable_intensity()

    def depends(self):
        dep=[self.model.place]
        dep.extend(self.detectable.depends())
        return dep

    def affected(self):
        return [self.model.place]

    def enabled(self, now):
        if ((self.detectable.intensity(now) is not None)
                and (self.model.place.state is False)):
            return (True, distributions.ExponentialDistribution(0.1, now))
        else:
            return (None, None)
    def fire(self, now, rng):
        self.model.place.state=True

class QuarantineIntensity(object):
    def __init__(self, model):
        self.model=model
    def depends(self):
        return [self.model.place]
    def intensity(self, now):
        """
        Returns true if quarantine in effect.
        """
        return self.model.place.state

class QuarantineModel(object):
    def __init__(self, farm):
        self.farm=farm
        self.place=QuarantinePlace()

    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        writer.add_transition(QuarantineTransition(self))

    def quarantine(self):
        return QuarantineIntensity(self)



###############################################################
# Farm contains models.
###############################################################
class SendIntensity(object):
    def __init__(self, farm):
        self.farm=farm
        self.quarantine=farm.quarantine()

    def depends(self):
        return self.quarantine.depends()
    def intensity(self, now):
        return self.quarantine.intensity(now)==False

class Farm(object):
    def __init__(self, name, size=1000):
        self.name=name
        self.size=size
        self.disease=DiseaseModel(self)
        self.quarantine=QuarantineModel(self)

    def write_places(self, writer):
        self.disease.write_places(writer)
        self.quarantine.write_places(writer)

    def write_transitions(self, writer):
        self.disease.write_transitions(writer)
        self.quarantine.write_transitions(writer)

    def infectious_intensity(self):
        return self.disease.infectious_intensity()

    def infection_partial(self):
        return self.disease.infection_partial()

    def detectable_intensity(self):
        return self.disease.detectable_intensity()

    def send_shipments(self):
        return SendIntensity(self)

    def receive_shipments(self):
        return ReceiveIntensity(self)


##############################################################
# Kernel-based neighbor infection
##############################################################
class InfectTransition(object):
    """
    This transition brings together pieces from different models
    into a full transition.
    """
    def __init__(self, intensity, action):
        self.intensity=intensity
        self.action=action

    def depends(self):
        deps=self.intensity.depends()
        deps.extend(self.action.depends())
        return deps

    def affected(self):
        return self.action.affected()

    def enabled(self, now):
        intensity=self.intensity.intensity(now)
        if intensity is not None and self.action.enabled(now):
            rate=0.5*intensity
            return (True, distributions.ExponentialDistribution(rate, now))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.action.fire(now, rng)


class InfectNeighborModel(object):
    """
    This is a scenario model for infection of one farm by another.
    """
    def __init__(self, farma, farmb):
        self.farma=farma
        self.farmb=farmb

    def write_places(self, writer):
        pass

    def write_transitions(self, writer):
        writer.add_transition(InfectTransition(
            self.farma.infectious_intensity(), self.farmb.infection_partial()))
        writer.add_transition(InfectTransition(
            self.farmb.infectious_intensity(), self.farma.infection_partial()))


##############################################################
# Movement restrictions
##############################################################
class RestrictionPlace(object):
    def __init__(self):
        self.restricted_date=None

class RestrictionTransition(object):
    def __init__(self, farm, restriction_place):
        self.farm=farm
        self.place=restriction_place
        self.detectable=farm.detectable_intensity()
        self.te=None

    def depends(self):
        dep=self.detectable.depends()
        dep.extend([self.place])
        return dep

    def affected(self):
        return [self.place]

    def enabled(self, now):
        detectable=self.detectable.intensity(now) is not None
        unrestricted=self.place.restricted_date is None
        if detectable and unrestricted:
            if self.te is None:
                self.te=now
            return (True, distributions.ExponentialDistribution(1, self.te))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.place.restricted_date=now

class RestrictionIntensity(object):
    def __init__(self, place):
        self.place=place
    def depends(self):
        return [self.place]
    def intensity(self, now):
        if self.restriction_place is None:
            return None
        else:
            # some kind of increasing restriction
            return (now-self.place.restricted_date)


class MovementRestrictionsModel(object):
    def __init__(self, landscape):
        self.landscape=landscape
        self.restriction_place=RestrictionPlace()

    def write_places(self, writer):
        writer.add_place(self.restriction_place)

    def write_transitions(self, writer):
        for farm in self.landscape.farms:
            t=RestrictionTransition(farm, self.restriction_place)
            writer.add_transition(t)
    def restriction_intensity(self):
        return RestrictionIntensity(self.restriction_place)



##############################################################
# Direct Movement Model
##############################################################
class DirectTransition(object):
    """
    Represents direct contact where a farm makes a certain
    number of shipments a day.
    """
    def __init__(self, farm, model):
        self.farm=farm
        self.model=model
    def depends(self):
        pass

class DirectModel(object):
    def __init__(self, landscape):
        self.landscape=landscape

    def write_places(self, writer):
        pass

    def write_transitions(self, writer):
        pass


class Landscape(object):
    def __init__(self):
        self.farm_locations=point_process.thomas_point_process_2D(
            5, 0.1, 5, (0, 1, 0, 1))
        individual_cnt=self.farm_locations.shape[0]
        self.distances=distance.squareform(distance.pdist(self.farm_locations,
            "euclidean"))
        self.farms=list()
        for add_idx in range(individual_cnt):
            self.farms.append(Farm(add_idx))

class Scenario(object): pass

def Build():
    net=llcp.LLCP()
    landscape=Landscape()
    farms=landscape.farms

    for f in farms:
        f.write_places(net)
        f.write_transitions(net)

    for a, b in itertools.combinations(farms, 2):
        infect=InfectNeighborModel(a, b)
        infect.write_places(net)
        infect.write_transitions(net)

    restrictions=MovementRestrictionsModel(landscape)
    restrictions.write_places(net)
    restrictions.write_transitions(net)

    scenario=Scenario()
    scenario.landscape=landscape
    scenario.farms=farms
    return net, scenario


########################################
# This is the part that runs the SIR
def observer(transition, when):
    if isinstance(transition, DiseaseABTransition):
        print("AB {0} {1} {2} {3}".format(transition.place.disease_model.farm.name,
            transition.a, transition.b, when))
    elif isinstance(transition, InfectTransition):
        print("Infect {0} {1} {2}".format(transition.intensity.place.disease_model.farm.name,
            transition.action.farm.farm.name,
            when))
    elif isinstance(transition, QuarantineTransition):
        print("Quarantine {0}".format(when))
    else:
        print("Unknown transition {0}".format(transition))
    return True

def test_farm():
    rng=np.random.RandomState()
    rng.seed(33333)
    net, scenario=Build()
    initial_idx=0
    scenario.farms[initial_idx].disease.infection_partial().fire(0, rng)
    sampler=sample.NextReaction(net, rng)
    run=runner.RunnerFSM(sampler, observer)
    run.init()
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_farm()