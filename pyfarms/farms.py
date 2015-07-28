import logging
import itertools
import copy
#from enum import Enum
import numpy as np
import numbers
import scipy.spatial.distance as distance
import gspn
import pyfarms.util as util

logger=logging.getLogger("farms")


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


class DiseaseABTransition(object):
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
        else:
            logger.error("Unknown distribution {0}".format(child.tag))
    return dist


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
            logger.debug("DM {0} {1}:{2}".format(t[0], start_state,
                end_state))
            self.add_transition(t[0], start_state, end_state, t[1][1],
                t[1][2:])

    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        for t in self.transitions:
            dist=t[3][1]
            args=t[3][2][0]
            logger.debug("write transitions dist {0} args {1}".format(
                type(dist), args))
            trans=DiseaseABTransition(self.place, t[1], t[2],
                lambda now : dist(*args, te=now))
            writer.add_transition(trans)

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
            return (True, gspn.ExponentialDistribution(0.1, now))
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

    def quarantine(self):
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
    def quarantine(self):
        return NoQuarantineIntensity(self)


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
        fm.quarantine=self.quarantine.clone(fm)
        return fm

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
            return (True, gspn.ExponentialDistribution(rate, now))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.action.fire(now, rng)


class InfectNeighborModel(object):
    """
    This is a scenario model for infection of one farm by another.
    """
    def __init__(self):
        self.farma=None
        self.farmb=None
        self.distance=None

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
            return (True, gspn.ExponentialDistribution(1, self.te))
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


class Premises(object): pass

class Landscape(object):
    """
    The landscape is the world upon which the rules act.
    It knows where farms are, what types of production are at
    each farm. It has data about the world.
    """
    def __init__(self):
        self.premises=list()
        self.farm_locations=np.zeros(0)
        self.distances=np.zeros((0,0))
        self.production_types=set()
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
            unit_size=int(unit.find("size").text)
            location=unit.find("location")
            lat=float(location.find("latitude").text)
            lon=float(location.find("longitude").text)
            latlon=np.array([lat, lon])
            self.add_premises(unit_name, unit_type, unit_size, latlon)
        self._build()

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

        self.airborne=list()
        for aidx, bidx in itertools.combinations(range(len(self.farms)), 2):
            a=self.farms[aidx]
            b=self.farms[bidx]
            from_type=landscape.premises[aidx].production_type
            to_type=landscape.premises[bidx].production_type
            dx=landscape.distances[aidx, bidx]
            air_model=self.spread_models[(from_type, to_type)].clone(a, b, dx)
            self.airborne.append(air_model)

 
    def write_gspn(self, net):
        """
        Given instances, write places and transitions into a net.
        """
        for f in self.farms:
            f.write_places(net)
            f.write_transitions(net)
        for airborne_instance in self.airborne:
            airborne_instance.write_places(net)
            airborne_instance.write_transitions(net)


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
        if root.find("quarantine-model", ns) is not None:
            self.quarantine=QuarantineModel()
        else:
            self.quarantine=NoQuarantineModel()
        self.spread_models=dict()
        for neighbor_model in models.findall(
                "airborne-spread-exponential-model", ns):
            from_production=neighbor_model.attrib["from-production-type"]
            to_production=neighbor_model.attrib["to-production-type"]
            im=InfectNeighborModel()
            im.from_naadsm_file(neighbor_model, ns)
            self.spread_models[(from_production, to_production)]=im

        self.farm_models=dict() # production_type => farm model
        for production_type in self.disease_by_type.keys():
            f=Farm()
            f.disease=self.disease_by_type[production_type]
            f.quarantine=self.quarantine
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