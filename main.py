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
    tf = 10.0                                                                # final time (days)
    J = 200                                                                 # number of x intervals
    N = 2000                                                                # number of t intervals
    seconds_per_day = 86400
    beta_tank = 0.75e-6 * seconds_per_day                                   # concrete heat tank. I choose a material with lower diffusivity to represent an insulated
    beta_rod = 22.8e-6 * seconds_per_day                                    # iron rod
    solar_energy_density = 1000.0 * seconds_per_day                         # joules / (day * meter^2)
    rod_area = 0.025 ** 2 * np.pi                                           # cross sectional area of rod in meters
    tank_area = 0.5 ** 2 * np.pi                                            # cross sectional area of tank in meters
    solar_panel_area = 1.0                                                  # in meter^2
    pump_speed_in_rod = 0.00001 * seconds_per_day                           # meters per day. Small number to represent slight increase in beta_rod. If it is too much, then heat flows to the end and accumulates
    heat_bleed_width = 0.05
    heat_bleed_amount = -80000                                              # Joules used per day. This appears touchy and dependent on J, N sizes

    # Set thermal diffusivities for system
    def beta(x):
        ret_beta = np.zeros(len(x))
        for i in range(0, len(x)):
            if tank_xl <= x[i] <= tank_xr:
                ret_beta[i] = beta_tank
            else:
                ret_beta[i] = beta_rod
        return ret_beta

    def area(x):
        ret_area = np.zeros(len(x))
        for i in range(0, len(x)):
            if tank_xl <= x[i] <= tank_xr:
                ret_area[i] = tank_area
            else:
                ret_area[i] = rod_area
        return ret_area

    # Take into account the contiuity of fluids
    def pump_speed(x):
        return np.zeros(len(x))
        # Model not working like I thought it would
        # ret_speed = np.zeros(len(x))
        # for i in range(0, len(x)):
        #     if tank_xl <= x[i] <= tank_xr:
        #         ret_speed[i] = pump_speed_in_rod * rod_area / tank_area
        #     else:
        #         ret_speed[i] = pump_speed_in_rod
        # return ret_speed

    # simulate day and night
    def solar_flux(t):
        flux = solar_energy_density * np.sin(2.0 * np.pi * t)
        # if flux < 0.0:
        #     return 0.0
        return flux

    # Assume no thermal energy in system
    def initial_conds(x):
        return 0.0 * x

    def source(x):
        ret_source = np.zeros(len(x))
        for i in range(0, len(x)):
            if tank_xl <= x[i] <= tank_xl + heat_bleed_width:
                ret_source[i] = heat_bleed_amount
        return ret_source

    solver = be.BackwardEuler(rod_xl=rod_xl, rod_xr=rod_xr, tank_xl=tank_xl, tank_xr=tank_xr, t0=t0, tf=tf, beta=beta,
                              J=J, N=N, solar_flux=solar_flux, area=area, pump_speed=pump_speed,
                              solar_panel_area=solar_panel_area, source=source, initial_conds=initial_conds)
    solver.solve()
    solver.plot_solution(tank_only=False)
    solver.plot_solution(tank_only=True)

    return 0


if __name__ == '__main__':
    err_code = main()
    exit(err_code)

else:
    exit(1)
