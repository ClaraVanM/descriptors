from aaindex import aaindex1
import numpy as np


def get_normalized_props(properties):
    aai_props = {prop: aaindex1[prop].values for prop in properties}
    # normalize
    for prop, dic in aai_props.items():
        for aa, value in dic.items():
            # "-":0 is not taken into account
            aai_props[prop][aa] = (value - np.mean(list(dic.values()))) / np.std(list(dic.values()))
    return aai_props


def moreaubroto_ac(sequence, lag=1, properties=["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102","CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]):
    """

    :param sequence: protein sequence
    :param lag: window length
    :param properties: the codes of the properties saved in the aaindex1 database
    :return: normalized moreau-broto autocorrelation
    """
    #check that lag has valid input
    if lag>=len(sequence) or lag <0:
        lag = 1
    sequence = sequence.upper()
    aai_props = get_normalized_props(properties)
    #calculate autocorrelation
    ac = {}
    for prop in properties:
        ac[prop] = 0
        for i in range(len(sequence)-lag):
                ac[prop] += aai_props[prop][sequence[i]]*aai_props[prop][sequence[i+lag]]
        ac[prop] = ac[prop]/(len(sequence)-lag)
    return ac


def moran_ac(sequence, lag=1, properties=["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102","CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]):
    sequence = sequence.upper()
    if lag >= len(sequence) or lag < 0:
        lag = 1
    aai_props = get_normalized_props(properties)
    moran = {}
    for prop in properties:
        numerator = 0
        for i in range(len(sequence)-lag):
                numerator += (aai_props[prop][sequence[i]]-np.mean(list(aai_props[prop].values())))*(aai_props[prop][sequence[i+lag]]-np.mean(list(aai_props[prop].values())))
        numerator = numerator/(len(sequence)-lag)
        devider = 0
        for i in range(len(sequence)):
            devider += (aai_props[prop][sequence[i]]-np.mean(list(aai_props[prop].values())))**2
        devider = devider/len(sequence)
        moran[prop] = numerator/devider
    return moran


def geary_ac(sequence, lag=1, properties=["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]):
    sequence = sequence.upper()
    if lag>=len(sequence) or lag<0:
        lag=1
    aai_props = get_normalized_props(properties)
    geary = {}
    for prop in properties:
        numerator = 0
        for i in range(len(sequence)-lag):
            numerator += (aai_props[prop][sequence[i]] - aai_props[prop][sequence[i+lag]])**2
        numerator = numerator/(2*(len(sequence)-lag))
        devider = 0
        for i in range(len(sequence)):
            devider += (aai_props[prop][sequence[i]]-np.mean(list(aai_props[prop].values())))**2
        devider = devider/(len(sequence)-1)
        geary[prop] = numerator/devider
    return geary

