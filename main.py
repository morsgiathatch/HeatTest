import numpy as np
from src import BackwardEuler as be


# Main driver for program
def main():
    # Set parameters
    rod_xl = 0.0                                                            # left end point of system (meters)
    rod_xr = 10.0                                                           # right end point of system (meters)
    tank_xl = 4.0                                                           # left end point of tank (meters)
    tank_xr = 6.0                                                           # right end point of tank (meters)
    t0 = 0.0                                                                # initial time (days)
    tf = 10.0                                                               # final time (days)
    J = 200                                                                 # number of x intervals
    N = 2000                                                                # number of t intervals
    seconds_per_day = 86400
    beta_tank = 0.5 * 0.75e-6 * seconds_per_day                             # concrete heat tank thermal diffusivity. Multiplying by 1/2 acts as a better store of heat
    beta_rod = 5 * 22.8e-6 * seconds_per_day                                # iron rod thermal diffusivity. Multiplying by 5 represents
    rod_spec_heat = 460.548                                                 # iron specific heat J / (kg * K)
    rod_dens = 7300                                                         # iron density kg / m^3
    solar_panel_width = 0.01                                                # solar panel width in meters. Only need solar panel width assuming a rectangular panel
    abs_rate = 0.215                                                        # absorption percentage of solar flux
    solar_energy_flux = (1000.0 * abs_rate * seconds_per_day /              # solar energy flux at earth surface at equator. Includes specific heat, density, and width for later convenience
                         (rod_spec_heat * rod_dens * solar_panel_width))
    rod_area = 0.025 ** 2 * np.pi                                           # cross sectional area of rod in meters^2
    tank_area = 0.5 ** 2 * np.pi                                            # cross sectional area of tank in meters^2

    # Set thermal diffusivities for system dependent on spatial variable
    def beta(x):
        ret_beta = np.zeros(len(x))
        for i in range(0, len(x)):
            if tank_xl <= x[i] <= tank_xr:
                ret_beta[i] = beta_tank
            else:
                ret_beta[i] = beta_rod
        return ret_beta

    # Set areas for system dependent on spatial variable
    def area(x):
        ret_area = np.zeros(len(x))
        for i in range(0, len(x)):
            if tank_xl <= x[i] <= tank_xr:
                ret_area[i] = tank_area
            else:
                ret_area[i] = rod_area
        return ret_area

    # Simulate day and night
    def solar_flux(t):
        flux = solar_energy_flux * np.sin(2.0 * np.pi * t)
        # Uncomment the below for insulation during night
        # if flux < 0.0:
        #     return 0.0
        return flux

    # Assume no thermal energy in system
    def initial_conds(x):
        return 0.0 * x

    # Assume no source
    def source(x):
        return np.zeros(len(x))

    # Solve system and plot solutions
    solver = be.BackwardEuler(rod_xl=rod_xl, rod_xr=rod_xr, tank_xl=tank_xl, tank_xr=tank_xr, t0=t0, tf=tf, beta=beta,
                              J=J, N=N, solar_flux=solar_flux, area=area, source=source, initial_conds=initial_conds)
    solver.solve()
    solver.plot_solution(tank_only=False)
    solver.plot_solution(tank_only=True)

    return 0


if __name__ == '__main__':
    err_code = main()
    exit(err_code)

else:
    exit(1)
