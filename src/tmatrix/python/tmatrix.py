"""
Sam Anderson
October 2017
"""

import warnings
import numpy as np
from pytmatrix.fortran_tm import pytmatrixf


class tmatrix(object):
    
    """T-Matrix scattering from nonspherical particles.

    Class for simulating scattering from nonspherical particles with the 
    T-Matrix method. Uses a wrapper to the Fortran code by M. Mishchenko.

    Usage instructions:

    First, the class should be be initialized. Any attributes (see below)
    can be passed as keyword arguments to the constructor. For example:
    sca = tmatrix.tmatrix(wavelength=2.0, m=complex(0,2))

    The properties of the scattering and the radiation should then be set 
    as attributes of this object. 

    Attributes: 
        radius: Equivalent radius.
        radius_type: If radius_type==tmatrix.RADIUS_EQUAL_VOLUME (default),
            radius is the equivalent volume radius.
            If radius_type==tmatrix.RADIUS_MAXIMUM, radius is the maximum 
            radius.
            If radius_type==tmatrix.RADIUS_EQUAL_AREA, 
            radius is the equivalent area radius.
        wavelength: The wavelength of incident light (same units as axi).
        m: The complex refractive index.
        axis_ratio: The horizontal-to-rotational axis ratio.
        shape: Particle shape. 
            tmatrix.SHAPE_SPHEROID: spheroid
            tmatrix.SHAPE_CYLINDER: cylinders; 
            tmatrix.SHAPE_CHEBYSHEV: Chebyshev particles (not yet 
                supported).
    """

    _attr_list = set(["radius_num","radius_max", "radius_type", "wavelength", 
        "m", "axial_ratio", "shape", "ddelt", "ndgs", "ndistr", "gpoints", 
        "ncoeff", "angles", "b_coeff", "gamma" ])

    _deprecated_aliases = {"axi": "radius",
        "lam": "wavelength",
        "eps": "axial_ratio",
        "rat": "radius_type",
        "np": "shape",
        "scatter": "orient"
    }

    RADIUS_EQUAL_VOLUME = 1.0
    RADIUS_EQUAL_AREA = 0.5
    RADIUS_MAXIMUM = 2.0

    SHAPE_SPHEROID = -1
    SHAPE_CYLINDER = -2
    SHAPE_CHEBYSHEV = 1

    DISTRIBUTION_MODIFIED_GAMMA = 1
    DISTRIBUTION_LOGNORMAL = 2
    DISTRIBUTION_POWERLAW = 3
    DISTRIBUTION_GAMMA = 4
    DISTRIBUTION_MODIFIED_POWERLAW = 5

    def __init__(self, **kwargs):
        self.radii = 1
        self.radius_max = 1.0
        self.radius_type = RADIUS_EQUAL_VOLUME
        self.wavelength = 1.0
        self.m = complex(1.53,0.008)
        self.axial_ratio = 1.0
        self.shape = SHAPE_SPHEROID
        self.ddelt = 1e-3
        self.distr = DISTRIBUTION_LOGNORMAL
        self.dvpoints = 2
        self.gqpoints = 5
        self.ncoeff = 60
        self.angles = 180
        self.bcoeff = 0.1
        self.gamma = 0.5
        
        self.reff = 0 
        self.veff = 0
        self.cext = 0
        self.csca = 0
        self.albedo = 0
        self.asym = 0
        self.f = 0

        for attr in self._attr_list:
            if attr in kwargs:
                self.__dict__[attr] = kwargs[attr]


    def generate_tables(self):
        """
        Generate aerosol tables using tmatrix.

        Args:
        """
        (self.reff, self.veff, self.cext, self.csca, self.albedo, self.asym, self.f) =\
        pytmatrixf.calcrand(self.radius_type, self.wavelength, self.m.real, \
        self.m.imag, self.axis_ratio, self.shape, self.ddelt, self.dvpoints, \
        self.radii, self.radius_max, self.bcoeff, self.gamma, \
        self.distr, self.gqpoints, self.angles, self.ncoeff)

        return (self.reff, self.veff, self.cext, self.csca, self.albedo, self.asym, self.f) 


    def mie(self, wl, m, sz, sig, nang ):
        """
        Generate Mie scattering values using lognormal distribution

        """
        B = 2 * sig**2
        
        (reff, veff, cext, csca, albedo, asym, f) =\
        pytmatrixf.calcrand(self.radius_type, wl, m.real, \
        m.imag, 1, SHAPE_SPHEROID, self.ddelt, self.dvpoints, \
        1, sz, B, 0, DISTRIBUTION_LOGNORMAL, self.gqpoints, nangs, self.ncoeff)

        qext = cext/(2*math.pi*reff**2)
        qsca = csca/(2*math.pi*reff**2)
        return (reff, veff, qext, qsca, albedo, asym, f) 


    def __getattr__(self, name):
        name = self._deprecated_aliases.get(name, name)  
        return object.__getattribute__(self, name)


    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)  


    def _init_tmat(self):
        """Initialize the scatter_random T-matrix.
        """

        if self.radius_type == Scatterer.RADIUS_MAXIMUM:
            # Maximum radius is not directly supported in the original
            # so we convert it to equal volume radius
            radius_type = Scatterer.RADIUS_EQUAL_VOLUME
            radius = self.equal_volume_from_maximum()
        else:
            radius_type = self.radius_type
            radius = self.radius

