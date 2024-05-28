import numpy as np

def time_estimate(num_event, world_radius, activity = 10.87):  # activity in Bq/kg
    volume = 4/3 * np.pi * world_radius**3
    mass = 1.04 * 10**3 * volume
    frequency = activity * mass
    total_time = np.random.gamma(num_event, 1/frequency)
    return total_time

if __name__ == "__main__":
    print(time_estimate(20000000*500, 30))
