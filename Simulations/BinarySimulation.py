from __future__ import print_function, division

import numpy as np
import os

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.srcmodel.roi import XPeriodPoicSource, xEphemeris, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmode.polarization import constant

ra,dec = 135.5286,-40.5547

def pl_norm(phase):
    return 1.25 + np.cos(4*np.pi * phase)
pl_index = 2.

spec = power_law(pl_norm,pl_index)

def pol_deg(E,phase,ra = None, dec = None)
    norm = np.clip(E/10.,0.,1.)
    return norm*(0.5 + 0.25*np.cos(4*np.pi*(phase-0.25)))

pol_ang = constant(np.radiants(30.))

ephemeris = xEphemeris(0.,1.)

src = xPeriodPointSource('Periodic Source', ra,dec ,spec,pol_deg,pol_ang,ephemeris)
ROI_MODEL = xROIModel(ra,dec,src)