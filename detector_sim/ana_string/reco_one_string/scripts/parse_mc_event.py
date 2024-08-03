import pandas as pd
import numpy as np
import json
from tqdm import tqdm

# Define a function to calculate energy from momentum components and mass
def calculate_energy(px, py, pz, mass):
    return np.sqrt(px**2 + py**2 + pz**2 + mass**2)

# Constants
proton_mass = 0.938  # GeV/c^2
muon_mass = 0.105    # GeV/c^2

# Load JSON data from a file
input_json_path = "../data/mc_event.json"
with open(input_json_path, 'r') as file:
    data = json.load(file)

# Prepare the data for the DataFrame
records = []

for event in tqdm(data, desc="Processing events"):
    event_id = event["event_id"]
    proton = event["particles_in"][0]
    proton_energy = calculate_energy(proton["px"], proton["py"], proton["pz"], proton_mass)
    
    if event["particles_at_detector"]:
        muon = event["particles_at_detector"][0]
        muon_p = np.sqrt(muon["px"]**2 + muon["py"]**2 + muon["pz"]**2)
        muon_energy = calculate_energy(muon["px"], muon["py"], muon["pz"], muon_mass)
        cos_zenith = np.arccos(muon["pz"] / muon_p)
        azimuthal = np.arctan2(muon["py"], muon["px"])
    else:
        muon_energy = np.nan
        cos_zenith = np.nan
        azimuthal = np.nan
    
    spectrum_weight = event["weights"]["spectrum"]
    
    records.append({
        "eventID": event_id,
        "particles_in_energy": proton_energy,
        "particles_out_energy": muon_energy,
        "zenith_angle_out": cos_zenith,
        "azimuthal_angle_out": azimuthal,
        "weight_spectrum": spectrum_weight
    })

# Create DataFrame
df = pd.DataFrame(records)

# Save to parquet
output_parquet_path = "../data/mc_event.parquet"
df.to_parquet(output_parquet_path)
df.to_csv(output_parquet_path.replace("parquet","csv"))

print(f"Data saved to {output_parquet_path}")

