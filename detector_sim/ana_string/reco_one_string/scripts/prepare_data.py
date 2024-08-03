import pandas as pd
import numpy as np
import uproot
import awkward as ak
import awkward_pandas
from sklearn.model_selection import train_test_split

### Hardcoded values are marked with TODO
n_hDOMs_per_string = 21

def Prepare_one_string(truth_file, sim_file, out_file_train, out_file_test):

    # Open the parquet file, which contains truth angle information
    df = pd.read_parquet(truth_file)

    # Open the ROOT file, which contains photon hits information
    file = uproot.open(sim_file)
    # TODO this is complicated..., try to improve
    string_ID = int( sim_file.split("string_")[1].split(".")[0] )

    # Access the tree
    tree = file["NoisyPmtHit"]
    pmthits_ak = tree.arrays(library="ak")

    # Use only photon hits associated with signals for now
    signal = pmthits_ak.Type == 1 #TODO
    pmthits_signal = ak.to_dataframe(pmthits_ak[signal])

    # Calculate the number of unique DomId values for each entry
    # In other words, I am curious for each event how many DOM will have photons
    unique_domId_counts = pmthits_signal.groupby('EventId')['DomId'].nunique()

    # Filter the entries with at least 5 unique DomId
    filtered_event_ids = unique_domId_counts[unique_domId_counts >= 5].index #TODO
    filtered_df = pmthits_signal[pmthits_signal['EventId'].isin(filtered_event_ids)]

    # Make a copy of the filtered DataFrame to avoid SettingWithCopyWarning
    filtered_df = filtered_df.copy()

    # Create the unique_PMT_ID column
    filtered_df['unique_PMT_ID'] = filtered_df['DomId'] * n_hDOMs_per_string + filtered_df['PmtId']
    
    '''
    The training data structure is the following:
    Each row is one event, it has 21 columns, each correspond to first hit time among its PMT hits
    Note in this case, the geometry information is encoded in the column order
    '''

    # Step 1: Group by Dom_ID and get the minimal t0 value in each group
    min_t0_by_pmt = filtered_df.groupby(['EventId','DomId'])['t0'].min().reset_index()
    #min_t0_by_pmt.head(30)

    # Step 2: Extract unique EventIds and DomIds
    unique_event_ids = min_t0_by_pmt['EventId'].unique()
    dom_ids = np.arange(n_hDOMs_per_string*string_ID, n_hDOMs_per_string*(string_ID+1)) #TODO DOMID contains stringID..., this is confusing in coding

    # Step 3: Create a new dataframe with zeros
    columns = [f'DomId_{dom_id}' for dom_id in dom_ids]
    new_df = pd.DataFrame(0.0, index=unique_event_ids.astype('int'), columns=columns)

    # Step 4: Populate the new dataframe with t0 values
    for _, row in min_t0_by_pmt.iterrows():
        event_id = row['EventId'].astype('int')
        dom_id = row['DomId'].astype('int')
        t0_value = row['t0']
        if dom_id >= n_hDOMs_per_string*(string_ID+1): continue #TODO DOMID contains stringID..., this is confusing in coding
        #print (event_id, f'DomId_{dom_id}')
        new_df.at[event_id, f'DomId_{dom_id}'] = t0_value

    # TODO this is because again DomID contains stringID info... not convenient for merging
    new_column_names = {col: f'DomId_{int(col.split("_")[1]) % n_hDOMs_per_string}' for col in new_df.columns if col.startswith('DomId')}
    new_df = new_df.rename(columns=new_column_names)

    # Reset the index of new_df to make eventID a column and rename it to 'eventID'
    new_df = new_df.reset_index().rename(columns={'index': 'eventID'})

    # Perform the merge operation
    merged_df = new_df.merge(df[['eventID', 'zenith_angle_out', 'azimuthal_angle_out']], on='eventID', how='left')

    # Split the dataframe into 80% training and 20% testing sets
    train_df, test_df = train_test_split(merged_df, test_size=0.2, random_state=42)

    # Save the training set as a parquet file
    train_df.to_parquet(out_file_train)

    # Save the testing set as a parquet file
    test_df.to_parquet(out_file_test)

truth_file = "/Users/meihualin/Projects/trident/siMu_atm/detector_sim/ana_string/reco_one_string/data/mc_event.parquet"
sim_file_base = "/Users/meihualin/Projects/trident/siMu_atm/detector_sim/ana_string/reco_one_string/data/2024.08.02_muon_withnoise/"
out_file_base = "/Users/meihualin/Projects/trident/siMu_atm/detector_sim/ana_string/reco_one_string/data/2024.08.02_muon_withnoise_ML/"

from glob import glob
from tqdm import tqdm


for sim_file in tqdm( glob(sim_file_base + "*") ):

    out_file_train = sim_file.replace("withnoise", "withnoise_ML").replace(".root", "_train.parquet")
    out_file_test = sim_file.replace("withnoise", "withnoise_ML").replace(".root", "_test.parquet")
    Prepare_one_string(truth_file, sim_file, out_file_train, out_file_test)


# Merge train and test dataset
train_pieces = []
for train_data in glob( out_file_base + "*train*"):
    train_pieces.append(pd.read_parquet(train_data))

test_pieces = []
for test_data in glob( out_file_base + "*test*"):
    test_pieces.append(pd.read_parquet(test_data))

train_dataset = pd.concat(train_pieces, ignore_index = True)
test_dataset = pd.concat(test_pieces, ignore_index = True)

train_dataset.to_parquet(out_file_base + "train.parquet")
test_dataset.to_parquet(out_file_base + "test.parquet")
