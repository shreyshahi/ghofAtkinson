''' Implementation of the BC Hydro model'''
import numpy as np
import itertools
import model

def interpolateSpectra(spectra1 , spectra2 , period1, period2, period):
    slope = (spectra2 - spectra1)/(np.log(period2) - np.log(period1))
    return spectra1 + slope*(np.log(period) - np.log(period1))

def fixPeriods(augmentedSpectra, periods, augmentedPeriods):
    fixedSpectra = np.zeros((augmentedSpectra.shape[0], periods.shape[0]))
    for i, per in enumerate(periods):
        if per in augmentedPeriods:
            fixedSpectra[:, i] = augmentedSpectra[:, 0]
            fixedSpectra[:, i] = augmentedSpectra[:, np.nonzero(augmentedPeriods == per)[0][0]]
        else:
            idxLow = np.nonzero(augmentedPeriods < per)[0][-1]
            idxHigh = np.nonzero(augmentedPeriods > per)[0][0]
            fixedSpectra[:, i] = interpolateSpectra(augmentedSpectra[:, idxLow] , augmentedSpectra[:, idxHigh] , augmentedPeriods[idxLow] , augmentedPeriods[idxHigh] , per)
    return fixedSpectra

def augmentPeriods(periods):
    availablePeriods = np.array([0.01, 0.02, 0.05, 0.075, 0.1, 0.15, 0.2, 0.250, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.5, 10.0])
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

def spectra(M, Rrup, Rhyp, eventType, Z, Faba, Vs30, periods):
    periods = np.array(periods)
    idx = np.nonzero(periods <= 0.01)[0]
    periods[idx] = 0.01 # Periods less than eq to 0.01 are essentially pga. Setting this allows us to avoid errors while interpolating in log-log scale.
    augmentedPeriods = augmentPeriods(periods)

    # compute R based on event type
    R = [rrup if ev == 0 else rhyp for rrup, rhyp, ev in itertools.izip(Rrup, Rhyp, eventType)]

    nRow = len(M)
    nCol = len(augmentedPeriods)

    M = np.array([M] * nCol).transpose()
    R = np.array([R] * nCol).transpose()
    eventType = np.array([eventType] * nCol).transpose()
    Z = np.array([Z] * nCol).transpose()
    Faba = np.array([Faba] * nCol).transpose()
    Vs30 = np.array([Vs30] * nCol).transpose()
    augmentedPeriods = np.array([augmentedPeriods] * nRow)

    augmentedSpectra = np.array([[model.computeSpectra(mag, r, evt, z, faba, vs, per) for mag, r, evt, z, faba, vs, per in itertools.izip(mags, rs, evts, zs, fabas, vss, pers)] for mags, rs, evts, zs, fabas, vss, pers in itertools.izip(M, R, eventType, Z, Faba, Vs30, augmentedPeriods)])

    fixedSpectra = fixPeriods(augmentedSpectra, periods, augmentedPeriods[0])

    intraEventSigma = np.ones(fixedSpectra.shape) * model.intraEventSigma()
    interEventSigma = np.ones(fixedSpectra.shape) * model.interEventSigma()

    return (np.exp(fixedSpectra)).tolist(), intraEventSigma.tolist(), interEventSigma.tolist()
