'''
BC hydro ground-motion model implementation
'''
import gmm

def spectra(M, Rrup, Rhyp, eventType, Z, Faba, Vs30, periods):
    '''
    The function takes the source, site, and distance parameters and returns the predicted response spectra, intra event sigma (:math:`\\phi`), and inter event sigma (:math:`\\tau`) at the requested periods.

    :param M: List of magnitudes for which the response spectra is needed.
    :type M: list(float)
    :param Rrup: List with closest distance to rupture for sites at which spectra is needed. Rrup is needed for Interface events. Filler values like -999 can be passed as Rrup for Intraslab events.
    :type Rrup: list(float)
    :param Rhyp: List with the hypocentral distances for sites at which spectra is needed. Rhyp is only needed for Intraslab events. Filler values like -999 can be passed as Rhyp for Interface events.
    :type Rhyp: list(float)
    :param eventType: List with 0 for Interface and 1 for Intraslab events.
    :type eventType: list(int)
    :param Z: List of hypocentral depths for each event.
    :type Z: list(float)
    :param Faba: List with 0 for forearcs and unknown sites, and 1 for backarc sites
    :type Faba: list(int) 0 for forearc or unknown, 1 for backarc
    :param Vs30: List of Vs30 at the sites where the spectra is needed.
    :type Vs30: list(float)
    :param periods: List of periods at which the spectra is needed.
    :type periods: list(float)
    '''

    argumentLengths = [len(M) , len(Rrup) , len(Rhyp) , len(eventType) , len(Z) , len(Faba) , len(Vs30)]

    if any([ (l != argumentLengths[0]) for l in argumentLengths]):
        raise Exception("Argument lengths not equal")

    if len(periods) == 0:
        raise Exception("Periods argument is empty")

    return gmm.spectra(M, Rrup, Rhyp, eventType, Z, Faba, Vs30, periods)
