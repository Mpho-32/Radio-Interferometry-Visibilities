## Mpho Manamela MNMMPH003 RI Assignment 1
## Inspecting Visibilities and Data Formats

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from casacore.tables import table, taql
from astropy.time import Time

def get_baseline_query(ant1, ant2):  ## TaQL query of getting the baseline from the antennas
    return f"ANTENNA1 == {ant1} && ANTENNA2 == {ant2}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate dynamic spectra from a CASA Measurement Set.")
    
    parser.add_argument("ms_path", help="Path to the Measurement Set directory")
    parser.add_argument("-f", "--field", help="Select field by name")
    parser.add_argument("-s", "--scans", nargs="+", help="Select specific scans")
    parser.add_argument("-a", "--antennas", help="Comma-separated list of antenna names")
    parser.add_argument("--data-column", default="DATA", help="Data column to plot")
    parser.add_argument("--spw", help="Spectral window ID to plot")

    args = parser.parse_args()

    ## Checks the path provided if there is a measurement set in it
    if not os.path.isdir(args.ms_path):
        parser.error(f"Measurement Set '{args.ms_path}' does not exist or is not a directory")

    ## Opening the measurement set
    ms = table(args.ms_path, readonly=True)
    ant1 = ms.getcol("ANTENNA1")
    ant2 = ms.getcol("ANTENNA2")
    baselines = list(set(zip(ant1, ant2)))

    ## Mapping antenna IDs to names
    ant_table = table(f"{args.ms_path}/ANTENNA")
    ant_names = ant_table.getcol("NAME")
    ant_table.close()

    ## Mapping field name to FIELD_ID
    field_id = None
    if args.field:
        field_table = table(f"{args.ms_path}/FIELD")
        field_names = field_table.getcol("NAME")
        for i, name in enumerate(field_names):
            if name == args.field:
                field_id = i  
                break
        field_table.close()
        if field_id is None:  
            parser.error(f"Field '{args.field}' not found in MS")

    ## Filtering baselines by antenna names 
    if args.antennas:
        selected_ants = args.antennas.split(",") ## Splitting the antennas with a comma
        ant_ids = [i for i, name in enumerate(ant_names) if name in selected_ants] ## Checking if names match 
        baselines = [(a1, a2) for a1, a2 in baselines if a1 in ant_ids and a2 in ant_ids] ## Includes only pairs where both antennas are in the selected ant_ids list

    ## Iterating over baselines and querying the 
    for ant1, ant2 in baselines:
        query = get_baseline_query(ant1, ant2)
        if args.field and field_id is not None: ## Checking if a field is specified and field_id was found
            query += f" && FIELD_ID == {field_id}"
        if args.scans: ## Checking if scan numbers were provided
            scans = [int(s) for s in args.scans]
            query += f" && SCAN_NUMBER IN {tuple(scans)}"
        if args.spw: ## Checking if a spectral window ID was provided
            query += f" && DATA_DESC_ID == {int(args.spw)}"
        data = ms.query(query) ## TaQL query on the Measurement Set to select the matching data
        print(f"Baseline {ant1}-{ant2} ({ant_names[ant1]}-{ant_names[ant2]}): {data.nrows()} rows selected")

        if data.nrows() == 0: ## Checking if the query returned no rows
            print(f"Skipping plotting for {ant_names[ant1]}-{ant_names[ant2]}: No data selected")
            continue

        # Extracting visibility data and metadata
        vis = data.getcol(args.data_column)  ## Getting the visibilities with the flags in them
        flags = data.getcol("FLAG")          ## Getting the flagged data
        vis[flags] = np.nan                  ## Masking flagged data
        amp = np.abs(vis)                    ## Amplitude is the visibilities without the flags
        time = data.getcol("TIME")           ## Time is in MJD seconds
        spw_table = table(f"{args.ms_path}/SPECTRAL_WINDOW")
        freq = spw_table.getcol("CHAN_FREQ")[int(args.spw or 0)] / 1e6  ## Frequency is in MHz
        spw_table.close()
        utc_time = Time(time / 86400, format="mjd", scale="utc").iso  ## Convert to UTC

        ## Plotting
        for corr in range(amp.shape[2]):
            plt.figure(figsize=(10, 6))
            plt.imshow(amp[:, :, corr].T, aspect="auto", origin="lower",
                    extent=[0, len(time) - 1, freq[0], freq[-1]], 
                    vmin=0, vmax=15, cmap='viridis')
            
            ## Set x-axis ticks to show time intervals (HH:MM:SS)
            num_ticks = 8  ## Number of ticks shown on the axes
            tick_indices = np.linspace(0, len(time) - 1, num_ticks, dtype=int)
            tick_labels = [t.split(" ")[1] for t in utc_time[tick_indices]]  ## Extract time part only without the date
            plt.xticks(tick_indices, tick_labels, rotation=45) ## Rotation prevents squashed ticks in the axes
            
            plt.xlabel("Time (UTC)")
            plt.ylabel("Frequency (MHz)")
            plt.title(f"Baseline: {ant_names[ant1]}-{ant_names[ant2]}, Corr {corr}")
            plt.colorbar(label="Amplitude")
            plt.tight_layout()  
            plt.savefig(f"/home/mpho/Documents/NASSP_MSc/RI/Assignment1/plots_multical/spectra_{ant_names[ant1]}_{ant_names[ant2]}_corr{corr}.png") ## Saving into a sub directory
            #plt.show()
            plt.close()

    ms.close()