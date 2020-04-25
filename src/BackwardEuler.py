import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


class BackwardEuler:
    def __init__(self, rod_xl, rod_xr, tank_xl, tank_xr, t0, tf, beta, J, N, solar_flux, initial_conds, area=lambda x: np.ones(len(x)), solar_panel_area=1.0, pump_speed=lambda x: np.zeros(len(x)), source=None):
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
        :param solar_flux: this represents the heat flux in at rod_xl. Must be non-negative during the day and non-positive at night.
        :type solar_flux: function
        :param initial_conds: initial conditions of the system
        :type initial_conds: function
        :param area: optional cross-sectional area function of the rod. Default is 1
        :type area: function
        :param solar_panel_area: optional area of solar panel. Default is 1
        :type solar_panel_area: float
        :param source: this is the optional source term. For our purposes it is none but it could be useful in modeling lateral heat loss
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
        self.initial_conds = initial_conds
        self.area = area
        self.solar_panel_area = solar_panel_area
        self.pump_speed = pump_speed
        self.source = source

        self.L = self.rod_xr - self.rod_xl                                  # Total system length
        self.H = self.tank_xr - self.tank_xl                                # Storage tank height
        self.T = self.tf - self.t0                                          # Total system time
        self.h = self.L / self.J                                            # Step size in x
        self.k = self.T / self.N                                            # Step size in t

        self.x = np.linspace(self.rod_xl, self.rod_xr, self.J + 1)
        self.t = np.linspace(self.t0, self.tf, self.N + 1)
        self.betas = np.array(self.beta(self.x))
        self.areas = np.array(self.area(self.x))
        if self.source is None:
            self.source = lambda x: 0.0 * x
        self.lambd = self.betas * self.k / self.h**2 * self.areas
        self.pump_speeds = np.array(self.pump_speed(self.x))

        self.U = np.zeros((self.J + 1, self.N + 1))

    # There are undoubtedly good solvers out there but I wanted to showcase what I know
    def solve(self):
        """
        Solve the heat equation using the backward euler finite difference scheme described in the readme
        :return:
        """
        # Construct system
        maindiag = -2.0 * np.ones(self.J + 1) * self.lambd
        subdiag = np.ones(self.J + 1) * self.lambd
        # subdiag[self.J - 1] *= 2.0                                         # Homogeneous Neumann conditions
        supdiag = np.ones(self.J + 1)
        supdiag[1] = 2
        supdiag *= self.lambd
        # bottomdiag = np.ones(self.J + 1) * self.lambd                      # Continuity at endpoints
        bottomdiag = np.zeros(self.J + 1)
        data = np.array([bottomdiag.tolist(), subdiag.tolist(), maindiag.tolist(), supdiag.tolist()])
        M = sp.spdiags(data, np.array([-1 * self.J, -1, 0, 1]), self.J + 1, self.J + 1)
        Idiag = np.ones(self.J + 1) * self.areas
        A = sp.spdiags(Idiag, np.array([0]), self.J + 1, self.J + 1) - M

        # Set initial conditions
        self.U[:, 0] = self.initial_conds(self.x)

        # Solve system
        for i in range(1, self.N):
            b = np.zeros(self.J + 1)
            b[0] = 2 * self.k * self.lambd[0] / self.areas[0] * self.solar_flux(self.t[i])
            b += self.areas * self.source(self.x)

            # Add pump term
            vu = np.zeros(len(self.U[:, i]))
            vu[1:self.J] = self.U[2:self.J + 1, i] * self.areas[2: self.J + 1]*self.pump_speeds[2: self.J + 1] - self.U[0: self.J - 1, i] * self.areas[0: self.J - 1]*self.pump_speeds[0: self.J - 1]
            vu[0] = self.U[1, i]*self.areas[1]*self.pump_speeds[1] - self.U[self.J, i]*self.areas[self.J]*self.pump_speeds[self.J]
            vu[self.J] = self.U[0, i]*self.areas[0]*self.pump_speeds[0] - self.U[self.J - 1, i]*self.areas[self.J - 1]*self.pump_speeds[self.J - 1]

            self.U[:, i + 1] = (la.spsolve(A, self.areas * self.U[:, i] + b + self.k * vu / (2.0 * self.h)))
            self.U[self.J, i + 1] = 0.0                                    # Homogeneous Dirichlet condition
            self.U[:, i + 1] = self.U[:, i + 1].clip(min=0)                # Don't allow negative heat

    # adapted from https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
    def plot_solution(self, tank_only=False):
        """
        Plot a 3D surface of the evolution of heat energy in the system.

        :param tank_only: Plot the solution only for the portion of x that lie in the tank
        :type tank_only: bool
        :return:
        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Plot the surface.
        tank_xl_index = 0
        tank_xr_index = 0
        for i in range(0, len(self.x)):
            if self.x[i] >= self.tank_xl:
                tank_xl_index = i
                break

        for i in range(tank_xl_index, len(self.x)):
            if self.x[i] >= self.tank_xr:
                tank_xr_index = i
                break

        # Make data.
        if tank_only:
            X, T = np.meshgrid(self.x[tank_xl_index:tank_xr_index + 1], self.t)
            U = np.transpose(self.U[tank_xl_index:tank_xr_index + 1, :])
        else:
            X, T = np.meshgrid(self.x, self.t)
            U = np.transpose(self.U)

        surf = ax.plot_surface(X, T, U, cmap=cm.coolwarm, linewidth=0, antialiased=False)

        # Customize the z axis.
        ax.set_xlabel("x (meters)")
        ax.set_ylabel("t (days)")
        #ax.set_zlabel("Heat Energy (Joules)")

        # Add a color bar which maps values to colors.
        plt.show()
