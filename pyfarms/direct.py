"""
Compare a transliteration of NAADSM 3.2's direct contact spread
with an algorithm which produces a stochastic, continuous-time
direct spread. This runs both on the same scenario and produces
plots to compare.
"""
import logging
import numpy as np
import scipy.spatial as spatial
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.path as path
import matplotlib.patches as patches
import point_process

logger=logging.getLogger(__file__)

class Scenario(object): pass

def farm_locations(farm_cnt):
    # Put farms in space with some non-uniform distribution.
    kappa=5
    mu=farm_cnt/kappa
    pts=point_process.thomas_point_process_2D(kappa, 0.2, mu, (0, 1, 0, 1))
    cnt=len(pts)
    logger.debug("There are {0} points".format(cnt))
    zone_radius=0.3
    infect_radius=0.6

    # Create a zone somewhere in the bounds and zone all inside.
    zone_center=np.random.rand(2).reshape(1,2)
    zone_distance=spatial.distance.cdist(pts, zone_center).flatten()
    logger.debug("zone distance {0}".format(zone_distance))
    zoned=(zone_distance<zone_radius).flatten()

    # Quarantine everything near the zone center, including outside it.
    p=np.zeros(cnt)
    p[zone_distance<infect_radius]=0.4
    quarantined=np.random.binomial(1, p)==1
    logger.debug("zoned {0}".format(zoned))
    logger.debug("quarantined {0}".format(quarantined))
    possible_starts=np.where(np.logical_and(
            np.logical_not(zoned), np.logical_not(quarantined))==True)[0]
    start_idx=np.random.choice(possible_starts)
    scenario=Scenario()
    scenario.pts=pts
    scenario.zoned=zoned
    scenario.quarantined=quarantined
    scenario.start_idx=start_idx
    scenario.zone_center=zone_center
    scenario.zone_radius=zone_radius
    scenario.infect_radius=infect_radius
    scenario.uniform_kernel_max=0.4
    scenario.epsilon=0.01
    return scenario


def naadsm_direct(scenario, run_cnt=1):
    """
    This tries to replicate A4.1 Direct contact spread from
    NAADSM Model Specification 1.2.1.
    """
    pts=scenario.pts
    zoned=scenario.zoned
    quarantined=scenario.quarantined
    start_idx=scenario.start_idx
    start_loc=pts[start_idx].reshape(1,2)

    #3a Make a list of those that can be a recipient.
    recipients_idx=np.where(quarantined==False)[0]
    recipients_idx=recipients_idx[recipients_idx!=start_idx]
    recipients=pts[recipients_idx]
    zoned_recipients=zoned[recipients_idx]

    distances=spatial.distance.cdist(recipients, start_loc).flatten()

    hits=np.zeros(pts.shape[0], dtype=np.int)
    for run_idx in range(run_cnt):
        #3c Sample a number, distance, from the movement distance distribution.
        # Using uniform here.
        travel_distance=np.random.uniform(0, scenario.uniform_kernel_max)

        #3d Choose the closest
        closest_idx=np.argmin(np.abs(distances-travel_distance))
        closest_distance=distances[closest_idx]
        possibles=np.where(distances-closest_distance<scenario.epsilon)[0]
        logger.debug("possibles {0}".format(possibles))

        #3e If forbidden by zoning, drop it.
        unzoned=possibles[zoned_recipients[possibles]==False]
        if np.any(unzoned):
            hits[recipients_idx[np.random.choice(unzoned)]]+=1
        else:
            pass # nothing sent
    return hits


class Exposure(object):
    def __init__(self, desired_distance, difference):
        self.desired_distance=desired_distance
        self.difference=difference
        self.target_idx=None
        self.target_distance=None
        self.cumulative_size=0

def find_in_list(distances, sizes, exposure, within_radius):
    """
    Corresponds to the check_and_choose function.
    Given a unit, ask which exposure might accept it as a target.
    """
    for target_idx in within_radius:
        for exp in exposure:
            difference=abs(distances[target_idx]-
                    exp.exposure_distance)
            if exp.target_idx is None:
                exp.target_idx=target_idx
                exp.target_distance=distances[target_idx]
                exp.difference=difference
                exp.cumulative_size+=sizes[target_idx]
            elif abs(difference-exp.difference)<scenario.epsilon:
                chance=np.random()< scenario.size[target_idx]/exp.cumulative_size
                allowed=zoned[target_idx] is False and quarantined[target_idx] is False
                if allowed and chance:
                    exp.target_idx=target_idx
                    exp.target_distance=distances[target_idx]
                    exp.difference=difference
                    exp.cumulative_size+=sizes[target_idx]
            elif difference < exp.difference:
                exp.target_idx=target_idx
                exp.target_distance=distances[target_idx]
                exp.difference=difference
                exp.cumulative_size+=sizes[target_idx]
            else:
                pass


def naadsm_code(scenario, run_cnt=1):
    """
    This follows the code from NAADSM 3.2 in contact-spread-model.c.
    """
    pts=scenario.pts
    zoned=scenario.zoned
    quarantined=scenario.quarantined
    start_idx=scenario.start_idx
    start_loc=pts[start_idx].reshape(1,2)
    distances=spatial.distance.cdist(recipients, start_loc).flatten()

    hits=np.zeros(pts.shape[0], dtype=np.int)

    within_radius=list()
    for rad_idx in range(pts.shape[0]):
        if (distance_sq(pts[rad_idx], start_loc)<max_distance_sq and
                rad_idx!=start_idx):
            within_radius.append(rad_idx)

    for run_idx in range(run_cnt):
        exposure_cnt=np.random.poisson(scenario.infection_rate)
        exposure=list()
        for e_idx in range(exposure_cnt):
            exposure_distance=
                    np.random.uniform(0, scenario.uniform_kernel_max,
                    size=exposure_cnt)
            exposure.append(Exposure(exposure_distance))
        find_in_list(distances, scenario.sizes, exposure, within_radius)
        if any([x.target_idx is None for x in exposure]):
            find_in_list(distances, scenario.sizes, exposure, range(pts.shape[0]))

        for res in exposure:
            if not (quarantined[res.target_idx] or zoned[res.target_idx]):
                hits[target_idx]+=1


def faithful(scenario, run_cnt=1):
    """
    This is an attempt at an equivalent continuous-time model.
    """
    pts=scenario.pts
    zoned=scenario.zoned
    quarantined=scenario.quarantined
    start_idx=scenario.start_idx
    start_loc=pts[start_idx].reshape(1,2)

    #3a Make a list of those that can be a recipient.
    recipients_idx=np.where(quarantined==False)[0]
    recipients_idx=recipients_idx[recipients_idx!=start_idx]
    recipients=pts[recipients_idx]
    zoned_recipients=zoned[recipients_idx]

    #3d For each recipient, find whether it will really receive.
    distances=spatial.distance.cdist(recipients, start_loc).flatten()
    sort_idx=np.argsort(distances)
    logger.debug("faithful distances {0}".format(distances[sort_idx]))
    # What is the probability of sending a truck, ignoring zones?
    prob_basket=np.zeros(len(distances), dtype=np.double)
    inner=0.0
    # The -1 assumes the max of the kernel doesn't include at least one farm.
    for didx in range(0, len(distances)-1):
        ptidx=sort_idx[didx]
        if inner<scenario.uniform_kernel_max:
            outer=0.5*(distances[ptidx]+distances[sort_idx[didx+1]])
            outer=min(outer, scenario.uniform_kernel_max)
            prob_basket[ptidx]=outer-inner
            inner=outer
    prob_basket/=np.sum(prob_basket)
    logger.debug("faith prob_basket {0}".format(prob_basket[sort_idx]))

    # Now take trucks that would have gone to zones and either
    # don't send them or give them to the nearby available farms.
    adj_prob=np.zeros(len(distances), dtype=np.double)
    for ridx in range(0, len(distances)-1):
        possibles=np.where(distances-distances[ridx]<scenario.epsilon)[0]
        unzoned=possibles[zoned_recipients[possibles]==False]
        if len(unzoned)>0:
            adj_prob[unzoned]+=prob_basket[ridx]/len(unzoned)

    # Losing trucks makes the probability less than one. Readjust
    # and take that probability out by reducing the number of times
    # we draw a value.
    total_prob=np.sum(adj_prob)
    adj_prob /= total_prob
    logger.debug("faith adj_prob {0}".format(adj_prob[sort_idx]))
    indices=np.array(range(len(adj_prob)))

    hits=np.zeros(pts.shape[0], dtype=np.int)
    bins=np.array(range(len(adj_prob)+1))
    recipient_hits=np.histogram(np.random.choice(indices, p=adj_prob,
        size=run_cnt*total_prob), bins=bins)[0]
    hits[recipients_idx]=recipient_hits
    return hits    


def plot_comparison(scenario, hits0, hits1):
    plot_scatter(scenario, hits0, "scatter0.pdf")
    plot_scatter(scenario, hits1, "scatter1.pdf")
    plot_histogram(hits0, hits1)


def plot_histogram(hits0, hits1):
    nonzero=set(np.nonzero(hits0)[0])
    nonzero.update(set(np.nonzero(hits1)[0]))
    hits_idx=list(nonzero)
    h0=hits0[hits_idx]
    h1=hits1[hits_idx]
    N=len(h0)
    ind=np.arange(N)
    width=0.35

    fig, ax=plt.subplots()
    rects0=ax.bar(ind, h0, width, color="blue")
    rects1=ax.bar(ind+width, h1, width, color="tan")

    ax.set_ylabel("Counts")
    ax.set_title("Comparison of Infection Counts")
    #ax.set_xticks(ind+width)
    ax.set_xlabel("By Unit")

    ax.legend( (rects0[0], rects1[0]), ("NAADSM", "CT"))
    plt.savefig("hithisto.pdf")



def plot_scatter(scenario, hits0, filename):
    pts=scenario.pts
    zoned=scenario.zoned
    quarantined=scenario.quarantined
    start_idx=scenario.start_idx
    max_size=200
    min_size=10

    hitpts=pts[hits0>0]
    hithits=hits0[hits0>0]
    hithits=hithits*max_size/np.max(hithits) + min_size

    pts2=np.delete(pts, start_idx, 0)
    hits2=np.delete(hits0, start_idx)
    zoned2=np.delete(zoned, start_idx)
    quarantined2=np.delete(quarantined, start_idx)
    mispts=pts2[hits2==0]
    miszone=zoned2[hits2==0]
    misquar=quarantined2[hits2==0]

    colors=np.zeros(len(mispts))
    colors[:]=0.9
    colors[miszone]=0.3
    colors[misquar]=0.5
    #logger.debug("plotsizes {0}".format(sizes))

    fig, ax=plt.subplots()

    plt.scatter(hitpts[:,0], hitpts[:,1], s=hithits, c="blue", marker='o',
        alpha=0.5)
    plt.scatter(mispts[:,0], mispts[:,1], s=min_size, c=colors, marker='s',
        alpha=0.5)
    plt.scatter(pts[start_idx,0], pts[start_idx,1], s=min_size, c="orange", marker='*',
        alpha=1)
    zc=scenario.zone_center[0]
    plt.scatter(zc[0], zc[1], s=min_size, c="orange", marker='x',
        alpha=1)
    zone=path.Path.circle(center=scenario.zone_center[0],
        radius=scenario.zone_radius)
    patch=patches.PathPatch(zone, facecolor="yellow", edgecolor="gray",
        alpha=0.1)
    ax.add_patch(patch)
    plt.savefig(filename, format="pdf")


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    cnt=100
    scenario=farm_locations(cnt)
    run_cnt=10000
    hits0=naadsm_direct(scenario, run_cnt)
    hits1=faithful(scenario, run_cnt)
    print("hits1 {0}".format(hits1))
    plot_comparison(scenario, hits0, hits1)

