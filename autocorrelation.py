from aaindex import aaindex1
import numpy as np


def get_normalized_props(properties):
    """

    :param properties: list of property identitynumbers from aaindex database
    :return: normalized propertie values
    """
    aai_props = {prop: aaindex1[prop].values for prop in properties}
    # normalize
    for prop, dic in aai_props.items():
        for aa, value in dic.items():
            # "-":0 is not taken into account
            aai_props[prop][aa] = (value - np.mean(list(dic.values()))) / np.std(list(dic.values()))
    return aai_props


def moreaubroto_ac(sequence, prop_values, k=1):
    """
    :param sequence: protein sequence
    :param prop_values: normalized amino acid property values
    :param k: window
    :return: normalized moreau-broto autocorrelation
    """
    count = 0
    ac = float(0)
    for i in range(len(sequence)-k):
        if not "X" in (sequence[i], sequence[i+k]):
            count +=1
            ac += prop_values[sequence[i]] * prop_values[sequence[i+k]]
    if count !=0:
        ac/count
    else:
        ac = 0
    return ac


def moran_ac(sequence, prop_values, k=1):
    """

    :param sequence: protein sequence
    :param prop_values: normalized amino acid property values
    :param k: window
    :return: moran autocorrelation value
    """
    numerator = float(0)
    count = 0
    mean = np.mean(list(prop_values.values()))
    for i in range(len(sequence)-k):
        if not 'X' in (sequence[i], sequence[i+k]):
            count += 1
            numerator += (prop_values[sequence[i]] - mean * (prop_values[sequence[i+k]] - mean))
    if count != 0:
          numerator /= count
    else:
        numerator = 0
    devider = 0
    count = 0
    for i in range(len(sequence)):
        if sequence[i] != "X":
            count +=1
            devider += (prop_values[sequence[i]] - mean)**2
    devider /= count
    return numerator/devider


def geary_ac(sequence, prop_values, k=1):
    """

    :param sequence: protein sequence
    :param prop_values: normalized amino acid property values
    :param k: window
    :return: geary autocorrelaton value
    """
    numerator = float(0)
    count = 0
    for i in range(len(sequence)-k):
        if not 'X' in (sequence[i], sequence[i+k]):
            count +=1
            numerator += (prop_values[sequence[i]] - prop_values[sequence[i+k]])**2
    if count !=0:
        numerator /= 2*count
    else: count = 0
    devider = 0
    count = 0
    mean = np.mean(list(prop_values.values()))
    for i in range(len(sequence)):
        if sequence[i] != "X":
            count += 1
            devider += (prop_values[sequence[i]] - mean)**2
    devider /= count
    return numerator/devider


def autocorrelation(sequence, function, lag = 30, properties=["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]):
    """

    :param sequence: protein sequence
    :param function: moreau_broto, moran or geary ac
    :param lag: max window
    :param properties: amino acid properties
    :return: list of values of autocorrelation from 1 to lag
    """
    aai_props = get_normalized_props(properties)
    ac_all_props = {}
    for prop, values in aai_props.items():
        ac = {}
        for i in range(lag):
            ac[i] = function(sequence, values, i)
        ac_all_props[prop] = ac
    return ac_all_props
