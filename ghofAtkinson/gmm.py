''' Implementation of the BC Hydro model'''
import numpy as np
import itertools
import model

def interpolateSpectra(spectra1 , spectra2 , period1, period2, period):
    slope = (spectra2 - spectra1)/(np.log10(period2) - np.log10(period1))
    return spectra1 + slope*(np.log10(period) - np.log10(period1))

def fixPeriods(augmentedSpectra, periods, augmentedPeriods):
    fixedSpectra = np.zeros((augmentedSpectra.shape[0], periods.shape[0]))
    sortedAugmented = sorted(augmentedPeriods)

    for i, per in enumerate(periods):
        if per in augmentedPeriods:
            fixedSpectra[:, i] = augmentedSpectra[:, 0]
            fixedSpectra[:, i] = augmentedSpectra[:, np.nonzero(augmentedPeriods == per)[0][0]]
        else:
            idxLow = np.nonzero(sortedAugmented < per)[0][-1]
            perLow = sortedAugmented[idxLow]
            idxLow = np.nonzero(augmentedPeriods == perLow)[0][0]

            idxHigh = np.nonzero(sortedAugmented > per)[0][0]
            perHigh = sortedAugmented[idxHigh]
            idxHigh = np.nonzero(augmentedPeriods == perHigh)[0][0]

            fixedSpectra[:, i] = interpolateSpectra(augmentedSpectra[:, idxLow] , augmentedSpectra[:, idxHigh] , augmentedPeriods[idxLow] , augmentedPeriods[idxHigh] , per)
    return fixedSpectra

def augmentPeriods(periods):
    availablePeriods = np.array([0.01, 0.07, 0.09, 0.11, 0.14, 0.18, 0.22, 0.27, 0.34, 0.42, 0.53, 0.65, 0.81, 1.01, 1.25, 1.56, 1.92, 2.44, 3.03, 3.7, 4.55, 5.88, 7.14, 9.09])
    modifiedPeriods = []
    for period in periods:
        if period in availablePeriods:
            modifiedPeriods.append(period)
        else:
            idxLow = np.nonzero(availablePeriods < period)[0][-1]
            modifiedPeriods.append(availablePeriods[idxLow])
            idxHigh = np.nonzero(availablePeriods > period)[0][0]
            modifiedPeriods.append(availablePeriods[idxHigh])
    return modifiedPeriods

def spectra(M, Rrup, Faba, Vs30, cascadia, epistemic, periods):
    periods = np.array(periods)
    idx = np.nonzero(periods <= 0.01)[0]
    periods[idx] = 0.01 # Periods less than eq to 0.01 are essentially pga. Setting this allows us to avoid errors while interpolating in log-log scale.
    augmentedPeriods = augmentPeriods(periods)

    R = Rrup
    nRow = len(M)
    nCol = len(augmentedPeriods)

    M = np.array([M] * nCol).transpose()
    R = np.array([R] * nCol).transpose()
    Faba = np.array([Faba] * nCol).transpose()
    Vs30 = np.array([Vs30] * nCol).transpose()
    cascadia = np.array([cascadia] * nCol).transpose()

    augmentedPeriods = np.array([augmentedPeriods] * nRow)

    augmentedSpectra = np.array([[model.computeSpectra(mag, r, faba, vs, cas, epistemic, per) for mag, r, faba, vs, cas, per in itertools.izip(mags, rs, fabas, vss, cass, pers)] for mags, rs, fabas, vss, cass, pers in itertools.izip(M, R, Faba, Vs30, cascadia, augmentedPeriods)])

    fixedSpectra = fixPeriods(augmentedSpectra, periods, augmentedPeriods[0])

    # convert to cm/s/s
    fixedSpectra = 10**fixedSpectra

    intraEventSigma = model.intraEventSigma(periods) * nRow
    interEventSigma = model.interEventSigma(periods) * nRow

    return fixedSpectra.tolist(), intraEventSigma, interEventSigma
