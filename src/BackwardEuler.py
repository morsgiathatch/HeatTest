import numpy as np

class BackwardEuler:
    def __init__(self, rod_xl, rod_xr, tank_xl, tank_xr, t0, tf, beta, J, N, solar_flux, source=None):
        """
        BackwardEuler class to save heat state for analysis

        :param rod_xl: left end point of rod
        :type rod_xl: float
        :param rod_xr: right end point of rod
        :type rod_xr: float
        :param tank_xl: left end point of storage tank
        :type tank_xl: float
        :param tank_xr: right end point of storage tank
        :type tank_xr: float
        :param t0: initial time
        :type t0: float
        :param tf: final time
        :type tf: float
        :param beta: function to find value of beta given x. This must have derivative 0 everywhere and be discontinuous for a small amount of points.
        :type beta: function
        :param J: number of intervals in x direction
        :type J: int
        :param N: number of intervals in t direction
        :type N: int
        :param solar_flux: This represents the heat flux in at rod_xl. Must be non-negative during the day and non-positive at night.
        :type solar_flux: function
        :param source: This is the source term. For our purposes it is none but it could be useful in modeling lateral heat loss
        :type source: function
        """
        self.rod_xl = rod_xl
        self.rod_xr = rod_xr
        self.tank_xl = tank_xl
        self.tank_xr = tank_xr
        self.t0 = t0
        self.tf = tf
        self.beta = beta
        self.J = J
        self.N = N
        self.solar_flux = solar_flux
        self.source = source

        self.L = self.rod_xr - self.rod_xl                                  # Total system length
        self.H = self.tank_xr - self.tank_xl                                # Storage tank height
        self.T = self.tf - self.t0                                          # Total system time
        self.h = self.L / self.J                                            # Step size in x
        self.k = self.T / self.N                                            # Step size in t

        self.x = np.linspace(self.rod_xl, self.rod_xr, self.J + 1)
        self.t = np.linspace(self.t0, self.tf, self.N + 1)
        self.betas = np.array(self.beta(self.x))

    # def solve(self, status_bar=False):
