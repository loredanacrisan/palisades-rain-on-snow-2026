"""SNOW-17 for one point at an hourly time step.

Python port of the SNOW-17 Fortran in NOAA-OWP/snow17 (https://github.com/NOAA-OWP/snow17,
src/snow19: PACK19.f, MELT19.f, ROUT19.f), Apache License 2.0; see LICENSE-snow17 in this folder.
Changes from the original: translated to Python; one zone with full snow cover (SI = 0, no areal
depletion), no ground melt, no snow depth, no forecaster adjustments, and rain/snow from the
single PXTEMP threshold. Model: Anderson (1973, 2006).
"""
import math

NEXLAG = 7                        # lag slots for an hourly step (5 // IDT + 2 in exsnow19)
JULDAY = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]


def melt19(idn, ta, s, p, idt):
    """Non-rain melt and heat-deficit change (MELT19, latitude < 54 N)."""
    mfmax, mfmin = p['MFMAX'] * idt / 6, p['MFMIN'] * idt / 6
    nmf, tipm = p['NMF'] * idt / 6, 1 - (1 - p['TIPM']) ** (idt / 6)
    mf = math.sin(idn * 2 * 3.1416 / 366) * (mfmax - mfmin) * 0.5 + (mfmax + mfmin) * 0.5
    cnhs = mf / mfmax * nmf * (s['TINDEX'] - min(ta, 0.0))
    s['TINDEX'] = min(s['TINDEX'] + tipm * (ta - s['TINDEX']), 0.0)
    return mf * max(ta - p['MBASE'], 0.0), cnhs


def rout19(excess, s, idt):
    """Lag and attenuate liquid water leaving the pack (ROUT19, full cover)."""
    fit = float(idt)
    packro = 0.0
    ex = s['EXLAG']
    if excess != 0.0:
        if excess < 0.1 or s['WE'] < 1.0:
            ex[0] += excess
        else:
            n = int((excess * 4.0) ** 0.3 + 0.5) or 1
            for i in range(1, n + 1):
                term = min(0.03 * fit / 6.0 * s['WE'] * n / (excess * (i - 0.5)), 150.0)
                flag = 5.33 * (1.0 - math.exp(-term))
                l2 = int((flag + fit) / fit + 1.0)
                por2 = (flag + fit - (l2 - 1) * idt) / fit
                ex[l2 - 1] += por2 * excess / n
                ex[l2 - 2] += (1.0 - por2) * excess / n
    if s['STORGE'] + ex[0] != 0.0:
        if s['STORGE'] + ex[0] < 0.1:
            packro = s['STORGE'] + ex[0]
            s['STORGE'] = 0.0
        else:
            el = ex[0] / fit
            term = min(500.0 * (el / 25.4) / (s['WE'] / 25.4) ** 1.3, 150.0)
            r1 = 1.0 / (5.0 * math.exp(-term) + 1.0)
            for _ in range(idt):
                out = (s['STORGE'] + el) * r1
                packro += out
                s['STORGE'] += el - out
            if s['STORGE'] <= 0.001:
                packro += s['STORGE']
                s['STORGE'] = 0.0
    s['EXLAG'] = ex[1:] + [0.0]
    return packro


def step(t_end, pxi, ta, s, p, pa, idt=1):
    """One PACK19 step for the interval ending at t_end (a datetime)."""
    kda = JULDAY[t_end.month - 1] + t_end.day
    i0 = JULDAY[2] + 21                               # March 21
    idn = kda - i0 if kda >= i0 else 365 - (i0 - kda)
    sfall = cnhspx = rain = rainm = 0.0
    if pxi > 0.0:
        if ta <= p['PXTEMP']:
            ts = min(ta, 0.0)
            sfall = pxi * p['SCF']
            s['WE'] += sfall
            cnhspx = -ts * sfall / 160.0
            if sfall > 1.5 * idt:
                s['TINDEX'] = ts
        else:
            rain = pxi
            rainm = 0.0125 * rain * max(ta, 0.0)
    d = dict(snowfall=sfall, rain=rain, melt=0.0, refreeze=0.0, outflow=0.0)
    pmelt, cnhs = melt19(idn, ta, s, p, idt)
    if rain > 0.25 * idt:                             # rain-on-snow energy balance
        ea = 0.90 * 2.7489e8 * math.exp(-4278.63 / (ta + 242.792))
        qn = 0.0612 * idt * (((ta + 273) * 0.01) ** 4 - 55.55)
        qe = 8.5 * (ea - 6.11) * p['UADJ'] * idt / 6
        qh = 7.5 * 0.000646 * pa * p['UADJ'] * idt / 6 * ta
        melt = max(qn + qe + qh + rainm, 0.0)
    else:
        melt = pmelt + rainm
    cnhs = max(cnhs, -s['NEGHS'])
    if melt >= s['WE'] and melt > 0:                  # the whole pack melts
        d['melt'] = s['WE'] + s['LIQW']
        d['outflow'] = d['melt'] + sum(s['EXLAG']) + s['STORGE'] + rain
        s.update(WE=0.0, LIQW=0.0, NEGHS=0.0, STORGE=0.0, EXLAG=[0.0] * NEXLAG)
        return d
    s['WE'] -= melt
    d['melt'] = melt
    water = melt + rain
    liqwmx = p['PLWHC'] * s['WE']
    s['NEGHS'] = min(max(s['NEGHS'] + cnhs + cnhspx, 0.0), 0.33 * s['WE'])
    excess = 0.0
    if water + s['LIQW'] >= liqwmx + (1 + p['PLWHC']) * s['NEGHS']:
        excess = water + s['LIQW'] - liqwmx - (1 + p['PLWHC']) * s['NEGHS']
        s['LIQW'] = liqwmx + p['PLWHC'] * s['NEGHS']
        d['refreeze'] = s['NEGHS']
    elif water >= s['NEGHS']:
        s['LIQW'] += water - s['NEGHS']
        d['refreeze'] = s['NEGHS']
    else:
        d['refreeze'] = water
    s['WE'] += d['refreeze']
    s['NEGHS'] -= d['refreeze']
    if s['NEGHS'] == 0.0:
        s['TINDEX'] = 0.0
    d['outflow'] = rout19(excess, s, idt)
    return d


def run(times_end, precip, temp, params, swe0, neghs0, tindex0, elev_m):
    """Run hourly SNOW-17. times_end: datetimes ending each hour; precip in mm, temp in C.
    Returns a list of dicts: SWE (ice + liquid), liquid water, heat deficit and the hour's fluxes, all mm."""
    e = elev_m / 100.0
    pa = 33.86 * (29.9 - 0.335 * e + 0.00022 * e ** 2.4)     # air pressure (mb) from elevation
    s = dict(WE=swe0, NEGHS=neghs0, LIQW=0.0, TINDEX=tindex0, STORGE=0.0, EXLAG=[0.0] * NEXLAG)
    out = []
    for t, px, ta in zip(times_end, precip, temp):
        d = step(t, px, ta, s, params, pa)
        d['liquid'] = s['LIQW'] + sum(s['EXLAG']) + s['STORGE']
        d['SWE'] = s['WE'] + d['liquid']
        d['heat_deficit'] = s['NEGHS']
        out.append(d)
    return out
