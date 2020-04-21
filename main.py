from src import BackwardEuler as be
import numpy as np


# Main driver for program
def main():
    # Set parameters
    rod_xl = 0.0                                                            # left end point of system (meters)
    rod_xr = 10.0                                                           # right end point of system (meters)
    tank_xl = 4.0                                                           # left end point of tank (meters)
    tank_xr = 6.0                                                           # left end point of tank (meters)
    t0 = 0.0                                                                # initial time (days)
    tf = 3.                                                                 # final time (days)
    J = 100                                                                 # number of x intervals
    N = 1000                                                                # number of t intervals
    seconds_per_day = 86400
    beta_tank = 0.75e-6 * seconds_per_day                                   # concrete tank
    beta_rod = 22.8e-6 * seconds_per_day                                    # iron rod
    solar_energy_density = 1000.0 * seconds_per_day                         # joules / (sec * meter^2)

    # Set thermal diffusivities for system
    def beta(x):
        ret_beta = np.zeros(len(x))
        for i in range(0, len(x)):
            if tank_xl <= x[i] <= tank_xr:
                ret_beta[i] = beta_tank
            else:
                ret_beta[i] = beta_rod
        return ret_beta

    # simulate day and night
    def solar_flux(t):
        return solar_energy_density * np.sin(2.0 * np.pi * t)

    # Assume no thermal energy in system
    def initial_conds(x):
        return 0.0 * x

    solver = be.BackwardEuler(rod_xl=rod_xl, rod_xr=rod_xr, tank_xl=tank_xl, tank_xr=tank_xr, t0=t0, tf=tf, beta=beta,
                              J=J, N=N, solar_flux=solar_flux, initial_conds=initial_conds)
    solver.solve()
    solver.plotSolution(tank_only=False)
    solver.plotSolution(tank_only=True)

    return 0


if __name__ == '__main__':
    err_code = main()
    exit(err_code)

else:
    exit(1)
