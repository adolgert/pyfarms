import logging
import itertools
import collections
import copy
#from enum import Enum
import numpy as np
import numbers
import scipy.spatial.distance as distance
import gspn
import pyfarms.util as util

logger=logging.getLogger("pyfarms.farms")


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
        self.when=None


class DiseaseABTransition(object):
    def __init__(self, farm, farm_place, a, b, distribution):
        self.farm=(farm,)
        self.place=farm_place
        self.a=a
        self.b=b
        self.dist=distribution
        self.te=None

    def __str__(self):
        return "Disease({0} {1}:{2})".format(self.farm[0].name, self.a, self.b)

    def depends(self):
        return [self.place]

    def affected(self):
        return [self.place]

    def enabled(self, now):
        if self.place.state==self.a:
            if self.te is not None:
                now=self.te
            return (True, self.dist(now))
        else:
            self.te=None
            return (False, None)

    def fire(self, now, rng):
        self.te=None
        self.place.when=now
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
            return self.place.when
        return None


class InfectPartial:
    """
    This is the action-part of a transition. It isn't a whole transition.
    """
    def __init__(self, model):
        self.farm=(model.farm,)
        self.model=model
    def depends(self):
        return [self.model.place]
    def affected(self):
        return [self.model.place]
    def enabled(self, now):
        if self.model.place.state in (DiseaseState.susceptible,):
            return True
        else:
            return False
    def fire(self, now, rng):
        assert(self.model.place.state==DiseaseState.susceptible)
        self.model.place.state=DiseaseState.latent
        logger.debug("InfectPartial fire {0}".format(self.farm[0].name))


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


def read_naadsm_pdf(pdf_owner, ns):
    pdf_root=pdf_owner.find("probability-density-function", ns)
    stage_name=pdf_root.attrib["name"]
    dist=None
    for child in pdf_root:
        if child.tag=="gamma":
            alpha=float(child.find("alpha", ns).text)
            beta=float(child.find("beta", ns).text)
            dist=(stage_name, gspn.GammaDistribution, alpha, beta)
        elif child.tag=="point":
            days=int(child.text)
            if days>0:
                dist=(stage_name, gspn.UniformDistribution, days-0.5, days+0.5)
        elif child.tag=="uniform":
            a=float(child.find("a", ns).text)
            b=float(child.find("b", ns).text)
            dist=(stage_name, gspn.UniformDistribution, a, b)
        elif child.tag=="units":
            logger.debug("Ignoring the units tag.")
        else:
            logger.error("Unknown distribution {0}".format(child.tag))
    return dist


def read_naadsm_relational(obs, ns):
    obs_rel=obs.find("relational-function", ns)
    if obs_rel is None:
        logger.error("Could not find relational-function in {0}".format(
                obs))
        return tuple()
    x=list()
    y=list()
    for val in obs_rel:
        if val.tag=="value":
            x.append(float(val[0].text))
            y.append(float(val[1].text))
        elif val.tag=="x-units":
            pass
        elif val.tag=="y-units":
            pass
        else:
            logger.error("unknown tag in relational-function")
    return (x, y)


class DiseaseModel(object):
    """
    This is a scenario model for disease state within a farm.
    """
    def __init__(self):
        self.farm=None
        self.place=DiseasePlace(self)
        self.transitions=list()

    def clone(self, farm):
        dm=copy.copy(self)
        dm.farm=farm
        dm.place=DiseasePlace(dm)
        return dm

    def add_transition(self, name, start_state, end_state,
            distribution, *distargs):
        dist=[name, distribution, distargs]
        self.transitions.append([name, start_state, end_state, dist])

    def initial_infection(self):
        logger.debug("Initial infection at {0}".format(self.farm.name))
        self.place.state=DiseaseState.latent

    def from_naadsm_file(self, root, ns):
        name_to_state={"latent-period" : DiseaseState.latent,
            "infectious-subclinical-period" : DiseaseState.subclinical,
            "infectious-clinical-period" : DiseaseState.clinical,
            "immunity-period" : DiseaseState.recovered}
        transitions=list()
        for stage in root:
            stage_id=stage.tag
            logger.debug("DiseaseModel state {0}".format(stage_id))
            dist=read_naadsm_pdf(stage, ns)
            if dist is not None:
                transitions.append([stage_id, dist])
        self.transitions=list()
        for tidx in range(len(transitions)):
            t=transitions[tidx]
            start_state=name_to_state[t[0]]
            if tidx+1<len(transitions):
                end_state=name_to_state[transitions[tidx+1][0]]
            else:
                end_state=DiseaseState.susceptible
            logger.debug("DM {0} {1}:{2} {3}".format(t[0], start_state,
                end_state, t[1][2:]))
            self.add_transition(t[0], start_state, end_state, t[1][1],
                t[1][2:])

    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        for t in self.transitions:
            dist=t[3][1]
            args=t[3][2][0]
            def make_dist(dist, args):
                return lambda enable : dist(*args, te=enable)
            trans=DiseaseABTransition(self.farm, self.place, t[1], t[2],
                make_dist(dist, args))
            writer.add_transition(trans)

    def infectious_intensity(self):
        """
        Are these animals shedding virus?
        """
        return InfectiousIntensity(self)

    def infection_partial(self):
        return InfectPartial(self)

    def detectable_intensity(self):
        """
        Would a vet be able to see signs of disease?
        """
        return DetectionIntensity(self)

###############################################################
# DetectReport model
###############################################################

class DetectionPlace(object):
    def __init__(self):
        self.reported=False


class DetectionTransition(object):
    """
    This is detection and reporting, so it's the act of reporting.
    One farm sets its detection to true and tells the global one.
    """
    def __init__(self, model):
        self.te=None
        self.model=model
        self.farm=(model.farm,)
    def depends(self):
        dep=[self.model.place]
        dep.extend(self.detectable.depends())
        dep.extend(self.global_intensity.depends())
        for o in self.observers:
            dep.extend(o.depends())
        return dep
    def affected(self):
        a=[self.model.place]
        for o in self.observers:
            a.extend(o.affected())
        return a
    def enabled(self, now):
        # self.detectable is whether the disease is detectable.
        if ((self.detectable.intensity(now) is not None)
                and (not self.model.place.reported)):
            # If and when was the first observation?
            # This will be None or it will be a time.
            first_observation=self.global_intensity.intensity(now)

            if self.te is not None:
                now=self.te
            return (True, gspn.ExponentialDistribution(4, now))
        else:
            self.te=None
            return (None, None)
    def fire(self, now, rng):
        logger.debug("DetectionTransition fire {0}".format(id(self.model.place)))
        self.model.place.reported=True
        self.te=None
        for observer in self.observers:
            observer.fire(now, rng)

class DetectedIntensity(object):
    def __init__(self, model):
        self.model=model
        self.place=model.place
    def depends(self):
        return [self.place]
    def intensity(self, now):
        return self.place.reported


class DetectionModel(object):
    def __init__(self, global_model):
        self.farm=None
        self.global_model=global_model
    def __str__(self):
        if self.farm is None:
            return "DetectionModel({0}, {1}, {2})".format(
                id(self.global_model), self.detect, self.report)
        else:
            return "DetectionModel({0}, {1}, {2}, {3}, {4})".format(
                self.farm.name, id(self.place), id(self.global_model),
                self.detect, self.report)
    def from_naadsm_file(self, detect_model, ns):
        clinobs=detect_model.find("prob-report-vs-time-clinical", ns)
        clinobs_rel=clinobs.find("relational-function")
        self.detect=read_naadsm_relational(clinobs, ns)

        obs=detect_model.find("prob-report-vs-time-since-outbreak", ns)
        obs_rel=obs.find("relational-function")
        self.report=read_naadsm_relational(obs, ns)

    def clone(self, farm):
        dm=copy.copy(self)
        dm.farm=farm
        dm.place=DetectionPlace()
        dm.observers=list([self.global_model.on_detect()])
        return dm
    def is_detected(self):
        return DetectedIntensity(self)

    def on_detect(self, observer):
        """
        Observers of detection are partial transitions, which means
        they fire, have dependencies, and have affected places.
        """
        self.observers.append(observer)
    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        dt=DetectionTransition(self)
        dt.observers=self.observers
        dt.detectable=self.farm.detectable_intensity()
        dt.global_intensity=self.global_model.when_detected()
        dt.observers.append(self.global_model.on_detect())
        writer.add_transition(dt)



###############################################################
# Quarantine model
###############################################################
class QuarantinePlace(object):
    def __init__(self):
        self.state=False

class QuarantineTransition(object):
    def __init__(self, model):
        logger.debug("Quarantine transition create")
        self.model=model
        self.farm=(model.farm,)
        self.detectable=model.farm.detection.is_detected()
        self.te=None

    def depends(self):
        dep=[self.model.place]
        dep.extend(self.detectable.depends())
        return dep

    def affected(self):
        return [self.model.place]

    def enabled(self, now):
        """
        Quarantine happens "a day after" reporting.
        """
        # If it has been reported and not quarantined.
        if (self.detectable.intensity(now)
                and (not self.model.place.state)):
            if self.te is not None:
                now=self.te
            return (True, gspn.UniformDistribution(0.5, 1.5, now))
        else:
            self.te=None
            return (None, None)

    def fire(self, now, rng):
        self.te=None
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
    def __init__(self):
        self.farm=None
        self.place=QuarantinePlace()        

    def clone(self, farm):
        qm=copy.copy(self)
        qm.farm=farm
        qm.place=QuarantinePlace()
        return qm

    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        writer.add_transition(QuarantineTransition(self))

    def quarantine_intensity(self):
        """
        Is quarantine in effect?
        """
        return QuarantineIntensity(self)


class NoQuarantineIntensity(object):
    def __init__(self, model):
        self.model=model
    def depends(self):
        return []
    def affected(self):
        return []
    def intensity(self, now):
        return False

class NoQuarantineModel(object):
    def __init__(self):
        self.farm=None
    def clone(self, farm):
        qm=copy.copy(self)
        qm.farm=farm
        return qm
    def write_places(self, writer):
        pass
    def write_transitions(self, writer):
        pass
    def quarantine_intensity(self):
        return NoQuarantineIntensity(self)


###############################################################
# Farm contains models.
###############################################################
class SendIntensity(object):
    def __init__(self, farm):
        self.farm=(farm,)
        self.quarantine=farm.quarantine.quarantine_intensity()

    def depends(self):
        return self.quarantine.depends()
    def intensity(self, now):
        return not self.quarantine.intensity(now)

class ReceiveIntensity(object):
    def __init__(self, farm):
        self.farm=(farm,)
        self.quarantine=farm.quarantine.quarantine_intensity()
    def depends(self):
        return self.quarantine.depends()
    def intensity(self, now):
        return not self.quarantine.intensity(now)

class Farm(object):
    def __init__(self):
        self.name=None
        self.size=None
        self.disease=DiseaseModel()
        self.quarantine=NoQuarantineModel()

    def clone(self, name, size):
        fm=copy.copy(self)
        fm.name=name
        fm.size=size
        fm.disease=self.disease.clone(fm)
        fm.detection=self.detection.clone(fm)
        fm.quarantine=self.quarantine.clone(fm)
        return fm

    def write_places(self, writer):
        self.disease.write_places(writer)
        self.detection.write_places(writer)
        self.quarantine.write_places(writer)

    def write_transitions(self, writer):
        self.disease.write_transitions(writer)
        self.detection.write_transitions(writer)
        self.quarantine.write_transitions(writer)

    def infectious_intensity(self):
        return self.disease.infectious_intensity()

    def infection_partial(self):
        """
        A partial transition to infect a susceptible unit.
        """
        return self.disease.infection_partial()

    def detectable_intensity(self):
        """
        An intensity indicating whether disease at the farm
        is detectable.
        """
        return self.disease.detectable_intensity()

    def send_shipments(self):
        """
        Would this farm be allowed to send animals?
        """
        return SendIntensity(self)

    def receive_shipments(self):
        """
        Would this farm be allowed to receive animals?
        """
        return ReceiveIntensity(self)


##############################################################
# Model for whether the authorities have seen any reports.
##############################################################
class GlobalDetectPlace(object):
    def __init__(self):
        self.detected=False
        self.when=None

class GlobalDetectPartial(object):
    def __init__(self, model):
        self.model=model
    def depends(self):
        return [self.model.place]
    def affected(self):
        return [self.model.place]
    def intensity(self, now):
        return True
    def fire(self, now, rng):
        if not self.model.place.detected:
            logger.debug("First global detection")
            self.model.place.detected=True
            self.model.place.when=now

class GlobalDetectIntensity(object):
    def __init__(self, model):
        self.model=model
    def depends(self):
        return [self.model.place]
    def affected(self):
        return []
    def intensity(self, now):
        return self.model.place.when

class GlobalDetectionModel(object):
    def __init__(self):
        self.place=GlobalDetectPlace()

    def on_detect(self):
        return GlobalDetectPartial(self)

    def when_detected(self):
        return GlobalDetectIntensity(self)

    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        # This class doesn't own its own transitions.
        pass

##############################################################
# Kernel-based neighbor infection
##############################################################
class InfectTransition(object):
    """
    This transition brings together pieces from different models
    into a full transition.
    """
    def __init__(self, farms, intensity, action, rate):
        self.farm=farms
        self.intensity=intensity
        self.action=action
        self.rate=rate

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
            return (True, gspn.ExponentialDistribution(rate, now))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.action.fire(now, rng)


class InfectNeighborModel(object):
    """
    This is a scenario model for infection of one farm by another.
    Converts distance into a hazard rate.
    Model for probability is exponential P=exp(-r d) where we are told
    0.5 = exp(-r 1) so r=-ln(0.5)
    Given probability for infection in a day, the hazard rate l is:
    P=1 - exp(- l t) so l=-ln(1-P)
    where t is 1 day.
    """
    def __init__(self):
        self.farma=None
        self.farmb=None
        self.distance=None

    def set_special_factor(self, special):
        self.special=special

    def clone(self, farma, farmb, distance):
        inm=copy.copy(self)
        inm.farma=farma
        inm.farmb=farmb
        inm.distance=distance
        return inm

    def from_naadsm_file(self, root, ns):
        # This is for the exponenital model.
        self.p=float(root.find("prob-spread-1km", ns).text)
        wind=root.find("wind-direction-start", ns)
        self.wind_angle_begin=float(wind.find("value").text)
        wind=root.find("wind-direction-end", ns)
        self.wind_angle_end=float(wind.find("value").text)
        delay=root.find("delay", ns)
        self.pdf=read_naadsm_pdf(delay, ns)
        self.hazard=lambda dx: np.power(self.p, dx)

    def write_places(self, writer):
        pass

    def write_transitions(self, writer):
        base=self.hazard(self.distance)
        rate=base*self.special[self.farma.size]*self.special[self.farmb.size]
        writer.add_transition(InfectTransition((self.farma, self.farmb),
            self.farma.infectious_intensity(), self.farmb.infection_partial(),
            rate))

    def herd_factor(self, premises):
        """
        This herd factor is used for airborne exponential models,
        and maybe more. It's some way to account for size of herds
        affecting spread. It's a kind of rescaling.
        """
        special_factor=2 # Comes from NAADSM
        histogram=collections.defaultdict(int)
        for p in premises:
            histogram[p.size]+=1
        total=sum(histogram.values())
        running=0
        cumulants=list()
        for size, count in sorted(histogram.items()):
            running+=count
            cumulants.append([size, special_factor*running/total])
        factor_dict=dict()
        previous=0
        for size, val in cumulants:
            factor_dict[size]=0.5*(previous+val)
            previous=val
        logger.debug("Creating herd special factor {0}".format(factor_dict))
        return factor_dict

##############################################################
# Indirect Contact
##############################################################

class IndirectTransition(object):
    def __init__(self, landscape, farm_models,
            source_farm, source_idx, rate, dist_pdf):
        self.farm=source_farm
        self.source_idx=source_idx
        self.landscape=landscape
        self.rate=rate
        self.distance_pdf=dist_pdf
        # Distances is an array matrix. Take the row and delete
        # the self-to-self distance.
        self.source_intensity=self.farm.infectious_intensity()
        self.sending=self.farm.send_shipments()
        self.infectable=list()
        self.receiving=list()
        for target in farm_models:
            # Is this the right question? Want that they aren't quarantined.
            self.infectable.append(target.infection_partial())
            self.receiving.append(target.receive_shipments())
        self.affected_idx=source_idx

    def __str__(self):
        return "Indirect({0}-{1} {2})".format(self.farm.name,
                self.landscape.premises[affected_idx].name,
                self.landscape.production_type)

    def depends(self):
        d=self.source_intensity.depends()
        d.extend(self.sending.depends())
        for infectable in self.infectable:
            d.extend(infectable.depends())
        for receive in self.receiving:
            d.extend(receive.depends())
        return d

    def affected(self):
        return self.infectable[self.affected_idx].affected()

    def enabled(self, now):
        can_send=self.sending.intensity(now)
        am_hot=self.source_intensity.intensity(now) is not None
        if (not can_send) or (not am_hot):
            return (False, None)

        self.current_recipients=list()
        current_distances=list()
        for tidx in range(len(self.infectable)):
            uninfected=self.infectable[tidx].enabled(now)
            receiving=self.receiving[tidx].intensity(now)
            if uninfected and receiving:
                self.current_recipients.append(tidx)
                current_distances.append(
                    self.landscape.distances[self.farm_idx, tidx])
        current_distances=np.array(current_distances)
        if len(self.current_recipients) is 0:
            return (False, None)
        sort_idx=np.argsort(current_distances)
        prob_basket=np.zeros(len(current_distances), dtype=np.double)
        assert(self.distance_pdf[1]==gspn.UniformDistribution)
        uniform_max=self.distance_pdf[3]
        inner=self.distance_pdf[2]
        for didx in range(0, len(current_distances)-1):
            ptidx=sort_idx[didx]
            if inner<uniform_max:
                outer=0.5*(current_distances[ptidx]+current_distances[sort_idx[didx+1]])
                outer=min(outer, uniform_max)
                prob_basket[ptidx]=outer-inner
                inner=outer
        total_prob=np.sum(prob_basket)
        if not total_prob>0.0:
            return (False, None)
        self.prob_basket=prob_basket/total_prob
        self.overall_rate=self.rate
        return (True, gspn.ExponentialDistribution(self.overall_rate, now))

    def fire(self, now, rng):
        try:
            self.affected_idx=np.random.choice(self.current_recipients,
                p=self.prob_basket)
        except ValueError as ve:
            logger.error("Cannot choose from probabilities: {0}".format(self.prob_basket))
        logger.debug("Indirect fire {0} {1}".format(self.source_idx, self.affected_idx))
        self.infectable[self.affected_idx].fire(now, rng)


class IndirectModel(object):
    """
    The NAADSM model for indirect contact, which means trucks without animals.
    There is a different model for each source production type. There will
    be an insance of this model for each source farm and for each
    possible destination production type.
    """
    def __init__(self):
        pass

    def clone(self, farm_models, landscape, farm_idx):
        logger.debug("IndirectModel clone farm {0}".format(farm_idx))
        im=copy.copy(self)
        im.landscape=landscape.single_production(self.to_production, farm_idx)
        im.source_farm=farm_models[farm_idx]
        im.farm_models=[farm_models[i] for i in im.landscape.farm_indices]
        im.source_idx=farm_idx
        im.source_farm=farm_models[farm_idx]
        return im

    def from_naadsm_file(self, root, ns):
        self.from_production=root.attrib["from-production-type"]
        self.to_production=root.attrib["to-production-type"]
        self.contact_type=root.attrib["contact-type"] # direct or indirect
        # movement per day
        self.movement_rate=float(root.find("movement-rate/value", ns).text)
        # distance kilometers
        self.distance_pdf=read_naadsm_pdf(root.find("distance"), ns)
        self.delay=read_naadsm_pdf(root.find("delay", ns), ns)
        self.probability_infect=float(root.find("prob-infect").text)
        self.latent_can_infect=bool(
                root.find("latent-units-can-infect", ns).text)
        self.subclinical_can_infect=bool(
                root.find("subclinical-units-can-infect", ns).text)
        self.movement_control=read_naadsm_relational(
                root.find("movement-control", ns), ns)


    def write_places(self, writer):
        pass

    def write_transitions(self, writer):
        t=IndirectTransition(self.landscape, self.farm_models, self.source_farm,
                self.source_idx, self.movement_rate, self.distance_pdf)
        writer.add_transition(t)

    def __str__(self):
        s=("IndirectModel({from_production}, {to_production}) "+
            "movement={movement_rate} probability infect {probability_infect} "+
            "latent can infect {latent_can_infect} subclinical can infect "+
            "{subclinical_can_infect}")
        return s.format(**self.__dict__)



##############################################################
# Movement restrictions
##############################################################

class RestrictionPlace(object):
    def __init__(self):
        self.restricted_date=None

class RestrictionTransition(object):
    def __init__(self, farm, restriction_place):
        self.farm=(farm,)
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
            return (True, gspn.ExponentialDistribution(1, self.te))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.te=None
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
        self.farm=(farm,)
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


class Premises(object): pass

class Landscape(object):
    """
    The landscape is the world upon which the rules act.
    It knows where farms are, what types of production are at
    each farm. It has data about the world.

    The landscape contains sub-copies of itself for each
    production type.
    """
    def __init__(self):
        self.premises=list()
        self.farm_locations=np.zeros(0)
        self.distances=np.zeros((0,0))
        self.production_types=set()
        self.production_landscapes=dict()
        # self.farm_locations=gspn.thomas_point_process_2D(
        #     5, 0.1, 5, (0, 1, 0, 1))
        # individual_cnt=self.farm_locations.shape[0]
        # self.distances=distance.squareform(distance.pdist(self.farm_locations,
        #     "euclidean"))
        # self.premises=list()
        # for add_idx in range(individual_cnt):
        #     self.premises.append(Farm(add_idx))

    def add_premises(self, name, production_type, size, location):
        f=Premises()
        f.name=name
        f.production_type=production_type
        self.production_types.add(production_type)
        f.size=size
        assert(len(location)==2)
        assert(isinstance(location[0], numbers.Real))
        assert(isinstance(location[1], numbers.Real))
        f.latlon=location
        self.premises.append(f)

    def from_naadsm_file(self, root, ns):
        for unit in root.findall("herd", ns):
            unit_name=unit.find("id").text
            unit_type=unit.find("production-type").text
            self.production_types.add(unit_type)
            unit_size=int(unit.find("size").text)
            location=unit.find("location")
            lat=float(location.find("latitude").text)
            lon=float(location.find("longitude").text)
            latlon=np.array([lat, lon])
            self.add_premises(unit_name, unit_type, unit_size, latlon)
        self._build()

    def single_production(self, production_type, source_farm_idx):
        """
        Look at just the part of the landscape that's one production type.
        Distances matrix is no longer square. It's from any type A to
        B of only a single type.
        """
        l=Landscape()
        l.production_type=production_type
        farm_indices=[i for (i, x) in enumerate(self.premises)
            if x.production_type==production_type]
        if source_farm_idx in farm_indices:
            farm_indices.remove(source_farm_idx)
        l.farm_indices=np.array(farm_indices)
        l.premises=[self.premises[ii] for ii in l.farm_indices]
        l.farm_locations=self.farm_locations[l.farm_indices]
        l.distances=self.distances[source_farm_idx, l.farm_indices]
        return l

    def _build(self):
        self.farm_locations=np.array([x.latlon for x in self.premises])
        self.distances=distance.squareform(
            distance.pdist(self.farm_locations, util.distancekm))
        logger.debug("found {0} premises".format(len(self.premises)))


class Scenario(object):
    """
    The scenario is the set of rules for how the simulation world works.
    """
    def __init__(self):
        self.models_loaded=False

    def build_from_landscape(self, landscape):
        """
        Given a landscape, make instances from the models.
        """
        assert(self.models_loaded)

        self.farms=list()
        for p in landscape.premises:
            f=self.farm_models[p.production_type].clone(p.name, p.size)
            self.farms.append(f)

        #### airborne
        if len(self.spread_models)>0:
            a_spread_model=next (iter (self.spread_models.values()))
            special_factor=a_spread_model.herd_factor(landscape.premises)
            for v in self.spread_models.values():
                v.set_special_factor(special_factor)
            self.airborne=list()
            for aidx, bidx in itertools.permutations(range(len(self.farms)), 2):
                a=self.farms[aidx]
                b=self.farms[bidx]
                from_type=landscape.premises[aidx].production_type
                to_type=landscape.premises[bidx].production_type
                dx=landscape.distances[aidx, bidx]
                air_model=self.spread_models[(from_type, to_type)].clone(a, b, dx)
                self.airborne.append(air_model)
        else:
            self.airborne=list()

        ### Indirect and Direct Contact
        self.indirect=list()
        if len(self.contact_models)>0:
            for ind_idx in range(len(self.farms)):
                f=self.farms[ind_idx]
                p=landscape.premises[ind_idx]
                for cm in self.contact_models[p.production_type]:
                    im=cm.clone(self.farms, landscape, ind_idx)
                    self.indirect.append(im)


    def clone(self):
        """
        A scenario that has been built from the landscape can be cloned.
        This will clone all of its state.
        """
        sc=copy.copy(self)
        sc.farms=list()
        for f in self.farms:
            sc.farms.append(f.clone(f.name, f.size))
        sc.airborne=list()
        for a in self.airborne:
            sc.airborne.append(a.clone(a.farma, a.farmb, a.distance))
        return sc
 
    def write_gspn(self, net):
        """
        Given instances, write places and transitions into a net.
        """
        self.global_detection.write_places(net)
        self.global_detection.write_transitions(net)

        for f in self.farms:
            f.write_places(net)
            f.write_transitions(net)

        for airborne_instance in self.airborne:
            airborne_instance.write_places(net)
            airborne_instance.write_transitions(net)

        for indirect_instance in self.indirect:
            indirect_instance.write_places(net)
            indirect_instance.write_transitions(net)


    def from_naadsm_file(self, root, ns):
        """
        Read input specifications in order to produce a closet of models.
        """
        models=root.find("models")
        self.disease_by_type=dict()
        self.disease_by_id=dict()
        for disease_model in models.findall("disease-model", ns):
            production_type=disease_model.attrib["production-type"]
            production_id=disease_model.attrib["production-type-id"]
            dm=DiseaseModel()
            dm.from_naadsm_file(disease_model, ns)
            self.disease_by_type[production_type]=dm
            self.disease_by_id[production_id]=dm
        logger.debug("result of find quarantine {0}".format(models.find(
            "quarantine-model", ns)))
        if models.find("quarantine-model", ns) is not None:
            logger.debug("Using quarantine model")
            self.quarantine=QuarantineModel()
        else:
            self.quarantine=NoQuarantineModel()

        self.global_detection=GlobalDetectionModel()
        self.detect_models=dict()
        for detect_model in models.findall("detection-model", ns):
            production_type=detect_model.attrib["production-type"]
            production_id=detect_model.attrib["production-type-id"]
            dm=DetectionModel(self.global_detection)
            dm.from_naadsm_file(detect_model, ns)
            self.detect_models[production_type]=dm
            logger.debug(dm)

        self.spread_models=dict()
        for neighbor_model in models.findall(
                "airborne-spread-exponential-model", ns):
            from_production=neighbor_model.attrib["from-production-type"]
            to_production=neighbor_model.attrib["to-production-type"]
            im=InfectNeighborModel()
            im.from_naadsm_file(neighbor_model, ns)
            self.spread_models[(from_production, to_production)]=im

        self.contact_models=collections.defaultdict(list)
        for indirect_model in models.findall("contact-spread-model", ns):
            from_production=indirect_model.attrib["from-production-type"]
            to_production=indirect_model.attrib["to-production-type"]
            contact_type=indirect_model.attrib["contact-type"]
            if contact_type=="indirect":
                inm=IndirectModel()
                inm.from_naadsm_file(indirect_model, ns)
                logger.debug("from_naadsm_file: indirect {0}".format(inm))
                self.contact_models[from_production].append(inm)
            elif contact_type=="direct":
                logger.warn("Ignoring direct contact model")
            else:
                logger.warn("Unknown contact spread model {0}".format(
                        contact_type))

        self.farm_models=dict() # production_type => farm model
        for production_type in self.disease_by_type.keys():
            f=Farm()
            f.production_type=production_type
            f.disease=self.disease_by_type[production_type]
            f.quarantine=self.quarantine
            f.detection=self.detect_models[production_type]
            self.farm_models[production_type]=f

        self.models_loaded=True



def Build(scenario, landscape):
    net=gspn.LLCP()
    scenario.build_from_landscape(landscape)
    scenario.write_gspn(net)
    return net


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
    sampler=gspn.NextReaction(net, rng)
    run=gspn.RunnerFSM(sampler, observer)
    run.init()
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_farm()