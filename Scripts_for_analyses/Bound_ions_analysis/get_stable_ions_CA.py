import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances
from collections import defaultdict
import pandas as pd
import os

# ---------------- Parameters ----------------
cutoff_minor_groove = 5
cutoff_major_groove = 5
cutoff_backbone = 5
cutoff_histone = 5
cutoff_dna = 5
cutoff_gap = 10.0
ions_selection = "resname CA"

binding_threshold = 200


# ---------------- Helper functions ----------------
def update_binding_counts(binding_counts, ion_resid, binding_resid, binding_resname):
    if ion_resid not in binding_counts:
        binding_counts[ion_resid] = {}
    binding_counts[ion_resid][(binding_resid, binding_resname)] = binding_counts[ion_resid].get((binding_resid, binding_resname), 0) + 1

def update_ion_counts(ion_counts, ion_sites, current_ion_resids, current_binding_sites):
    for ion_resid, sites in zip(current_ion_resids, current_binding_sites):
        if ion_resid not in ion_counts:
            ion_counts[ion_resid] = 1
            ion_sites[ion_resid] = {site: 1 for site in sites}
        else:
            ion_counts[ion_resid] += 1
            for site in set(sites):
                ion_sites[ion_resid][site] = ion_sites[ion_resid].get(site, 0) + 1
                
def update_gap_counts(ion_counts, ion_resid):
    """Update counts for ions in the DNA gap region."""
    if ion_resid not in ion_counts:
        ion_counts[ion_resid] = 1
    else:
        ion_counts[ion_resid] += 1
        

def process_frame(CA_ions, upper_backbone, lower_backbone):
    """Process one trajectory frame for ion binding in the gap region."""
    gap_counts = {}
    binding_counts = {}

    dist_to_upper = distances.distance_array(CA_ions.positions, upper_backbone.positions)
    dist_to_lower = distances.distance_array(CA_ions.positions, lower_backbone.positions)

    for i, ion in enumerate(CA_ions):
        counted_sites = set()

        in_gap = (dist_to_upper[i].min() < cutoff_gap) and (dist_to_lower[i].min() < cutoff_gap)
        if in_gap:
            update_gap_counts(gap_counts, ion.resid)

            for j, dist in enumerate(dist_to_upper[i]):
                if dist < cutoff_gap:
                    site = (upper_backbone[j].resid, upper_backbone[j].resname)
                    if site not in counted_sites:
                        update_binding_counts(binding_counts, ion.resid, *site)
                        counted_sites.add(site)

            for k, dist in enumerate(dist_to_lower[i]):
                if dist < cutoff_gap:
                    site = (lower_backbone[k].resid, lower_backbone[k].resname)
                    if site not in counted_sites:
                        update_binding_counts(binding_counts, ion.resid, *site)
                        counted_sites.add(site)

    return gap_counts, binding_counts
            
def get_stable_ions(ion_counts, ion_sites, binding_threshold):
    """Return ions that meet the binding threshold in a given region."""
    stable_ions = {}
    for ion_resid, count in ion_counts.items():
        if count >= binding_threshold:
            stable_sites = {site: site_count for site, site_count in ion_sites[ion_resid].items() if site_count >= binding_threshold}
            if stable_sites:
                stable_ions[ion_resid] = stable_sites
    return stable_ions

def get_stable_joint_binding_sites(ion_sites, binding_threshold):
    """Return ions that are stably bound to both DNA and histone (interface)."""
    stable_ions = {}
    for ion_resid, site_dict in ion_sites.items():
        stable_sites = {site: count for site, count in site_dict.items() if count >= binding_threshold}
        if stable_sites:
            stable_ions[ion_resid] = stable_sites
    return stable_ions

def calculate_ion_indices_in_region(ions, region, cutoff):
    """Calculate ions within cutoff of a region and record binding sites."""
    if len(region) == 0:
        return np.array([]), []
    distances_array = distances.distance_array(ions.positions, region.positions)
    ion_indices = np.where(distances_array.min(axis=1) < cutoff)[0]
    all_binding_sites = []
    for i in ion_indices:
        close_sites_indices = np.where(distances_array[i] < cutoff)[0]
        binding_sites = [(region[idx].resid, region[idx].resname) for idx in close_sites_indices]
        all_binding_sites.append(binding_sites)
    return ions[ion_indices].resids, all_binding_sites

# ---------------- Main function ----------------
def process_system_and_run(system_name, topology_file, trajectory_file, run_index):
    print(f"Processing {system_name}, Run {run_index + 1}...")
    u = mda.Universe(topology_file, trajectory_file)

    # Define regions
    minor_groove = u.select_atoms("((resname DA) and (name N3 or name H2 )) or ((resname DT) and (name O2 )) or ((resname DG) and (name N3 or name N2)) or ((resname DC) and (name O2))")
    major_groove = u.select_atoms("((resname DA) and (name N7 or name N6 or name H8)) or ((resname DT ) and (name O4 or name C7)) or ((resname DG) and (name N7 or name O6 or name H8)) or ((resname DC) and (name N4 or name C5))")
    dna_backbone = u.select_atoms("(resname DA DT DG DC) and (name C3' or name C4' or name C5' or name O3' or name O5' or name P or name OP1 or name OP2)")
    histone = u.select_atoms("protein and not name H*")
    histone_tail = u.select_atoms('(resid 803:839 or resid 938:974 or resid 1073:1092 or resid 1175:1194 or resid 295:307 or resid 412:423 or resid 424:436 or resid 541:552 or resid 553:578 or resid 678:703) and not name H*')
    histone_core = u.select_atoms('(protein and not (resid 803:839 or resid 938:974 or resid 1073:1092 or resid 1175:1194 or resid 295:307 or resid 412:423 or resid 424:436 or resid 541:552 or resid 553:578 or resid 678:703)) and not name H*')
    acidic_patch = u.select_atoms("(protein and ((resid 349) or (resid 354) or (resid 357) or (resid 384) or (resid 385) or (resid 386) or (resid 479) or (resid 484) or (resid 487) or (resid 513) or (resid 514) or (resid 515) or (resid 657) or (resid 665) or (resid 782) or (resid 790))) and not name H*")
    acidic_residue = u.select_atoms("(protein and (resname GLU or resname ASP)) and not name H*")
    dna = u.select_atoms("nucleic and not name H*")
    upper_backbone = u.select_atoms("(resid 1:68 or resid 227:294) and (name C3' or name C4' or name C5' or name O3' or name O5' or name P or name OP1 or name OP2)")  
    lower_backbone = u.select_atoms("(resid 148:215 or resid 80:147) and (name C3' or name C4' or name C5' or name O3' or name O5' or name P or name OP1 or name OP2)")  
    CA_ions = u.select_atoms(ions_selection)

     # Initialize counters
    ions_in_minor_groove_counts, ions_in_minor_groove_sites = {}, {}
    ions_in_major_groove_counts, ions_in_major_groove_sites = {}, {}
    ions_in_backbone_counts, ions_in_backbone_sites = {}, {}
    ions_near_histone_counts, ions_near_histone_sites = {}, {}
    ions_near_dna_counts, ions_near_dna_sites = {}, {}
    ions_in_acidic_patch_counts, ions_in_acidic_patch_sites = {}, {}
    ions_in_acidic_residue_counts, ions_in_acidic_residue_sites = {}, {}
    ions_in_histone_tail_counts, ions_in_histone_tail_sites = {}, {}
    ions_in_histone_core_counts, ions_in_histone_core_sites = {}, {}
    ions_near_both_sites = {}
    total_gap_counts = {}
    total_binding_counts = {}

    # Loop over trajectory
    for ts in u.trajectory[2500:12500:25]:
        ions_in_minor_groove_current, minor_groove_sites_current = calculate_ion_indices_in_region(CA_ions, minor_groove, cutoff_minor_groove)
        ions_in_major_groove_current, major_groove_sites_current = calculate_ion_indices_in_region(CA_ions, major_groove, cutoff_major_groove)
        ions_in_backbone_current, backbone_sites_current = calculate_ion_indices_in_region(CA_ions, dna_backbone, cutoff_backbone)
        ions_near_histone_current, histone_sites_current = calculate_ion_indices_in_region(CA_ions, histone, cutoff_histone)
        ions_near_dna_current, dna_sites_current = calculate_ion_indices_in_region(CA_ions, dna, cutoff_dna)
        ions_in_acidic_patch_current, acidic_patch_sites_current = calculate_ion_indices_in_region(CA_ions, acidic_patch, cutoff_histone)
        ions_in_acidic_residue_current, acidic_residue_sites_current = calculate_ion_indices_in_region(CA_ions, acidic_residue, cutoff_histone)
        ions_in_histone_tail_current, histone_tail_sites_current = calculate_ion_indices_in_region(CA_ions, histone_tail, cutoff_histone)
        ions_in_histone_core_current, histone_core_sites_current = calculate_ion_indices_in_region(CA_ions, histone_core, cutoff_histone)

        update_ion_counts(ions_in_minor_groove_counts, ions_in_minor_groove_sites, ions_in_minor_groove_current, minor_groove_sites_current)
        update_ion_counts(ions_in_major_groove_counts, ions_in_major_groove_sites, ions_in_major_groove_current, major_groove_sites_current)
        update_ion_counts(ions_in_backbone_counts, ions_in_backbone_sites, ions_in_backbone_current, backbone_sites_current)
        update_ion_counts(ions_near_histone_counts, ions_near_histone_sites, ions_near_histone_current, histone_sites_current)
        update_ion_counts(ions_near_dna_counts, ions_near_dna_sites, ions_near_dna_current, dna_sites_current)
        update_ion_counts(ions_in_acidic_patch_counts, ions_in_acidic_patch_sites, ions_in_acidic_patch_current, acidic_patch_sites_current)
        update_ion_counts(ions_in_acidic_residue_counts, ions_in_acidic_residue_sites, ions_in_acidic_residue_current, acidic_residue_sites_current)
        update_ion_counts(ions_in_histone_tail_counts, ions_in_histone_tail_sites, ions_in_histone_tail_current, histone_tail_sites_current)
        update_ion_counts(ions_in_histone_core_counts, ions_in_histone_core_sites, ions_in_histone_core_current, histone_core_sites_current)
        gap_counts, binding_counts = process_frame(CA_ions, upper_backbone, lower_backbone)
       
        # Accumulate gap results
        for ion_index, count in gap_counts.items():
            if ion_index in total_gap_counts:
                total_gap_counts[ion_index] += count
            else:
                total_gap_counts[ion_index] = count

        for ion_index, sites in binding_counts.items():
            if ion_index not in total_binding_counts:
                total_binding_counts[ion_index] = {}
            for (resid, resname), site_count in sites.items():
                if (resid, resname) in total_binding_counts[ion_index]:
                    total_binding_counts[ion_index][(resid, resname)] += site_count
                else:
                    total_binding_counts[ion_index][(resid, resname)] = site_count
                    
        # Interface: ions bound to both DNA and histone
        common_ions = set(ions_near_histone_current).intersection(set(ions_near_dna_current))
        for ion_resid in common_ions:
            histone_index = np.where(ions_near_histone_current == ion_resid)[0][0]
            dna_index = np.where(ions_near_dna_current == ion_resid)[0][0]
            histone_binding_sites = histone_sites_current[histone_index]
            dna_binding_sites = dna_sites_current[dna_index]
            if ion_resid not in ions_near_both_sites:
                ions_near_both_sites[ion_resid] = {}
            frame_pairs = set((d_site, h_site) for d_site in dna_binding_sites for h_site in histone_binding_sites)
            for pair in frame_pairs:
                ions_near_both_sites[ion_resid][pair] = ions_near_both_sites[ion_resid].get(pair, 0) + 1

    # Stable binding results
    ions_in_dna_gap_stable = get_stable_ions(total_gap_counts, total_binding_counts, binding_threshold)
    ions_in_minor_groove_stable = get_stable_ions(ions_in_minor_groove_counts, ions_in_minor_groove_sites, binding_threshold)
    ions_in_major_groove_stable = get_stable_ions(ions_in_major_groove_counts, ions_in_major_groove_sites, binding_threshold)
    ions_in_backbone_stable = get_stable_ions(ions_in_backbone_counts, ions_in_backbone_sites, binding_threshold)
    ions_near_histone_stable = get_stable_ions(ions_near_histone_counts, ions_near_histone_sites, binding_threshold)
    ions_near_dna_stable = get_stable_ions(ions_near_dna_counts, ions_near_dna_sites, binding_threshold)
    ions_in_acidic_patch_stable = get_stable_ions(ions_in_acidic_patch_counts, ions_in_acidic_patch_sites, binding_threshold)
    ions_in_acidic_residue_stable = get_stable_ions(ions_in_acidic_residue_counts, ions_in_acidic_residue_sites, binding_threshold)
    ions_in_histone_tail_stable = get_stable_ions(ions_in_histone_tail_counts, ions_in_histone_tail_sites, binding_threshold)
    ions_in_histone_core_stable = get_stable_ions(ions_in_histone_core_counts, ions_in_histone_core_sites, binding_threshold)
    ions_near_both_stable = get_stable_joint_binding_sites(ions_near_both_sites, binding_threshold)

    rows = []
    for region, stable_dict in [
        ("Backbone", ions_in_backbone_stable),
        ("Major_Groove", ions_in_major_groove_stable),
        ("Minor_Groove", ions_in_minor_groove_stable),
        ("Histone_Tail", ions_in_histone_tail_stable),
        ("Histone_Core", ions_in_histone_core_stable),
        ("Acidic_Patch", ions_in_acidic_patch_stable),
        ("Gap_Region", ions_in_dna_gap_stable)
    ]:
        for ion, sites in stable_dict.items():
            for (resid, resname), count in sites.items():
                rows.append({
                    "Region": region,
                    "Ion_Resid": ion,
                    "Binding_Resid": resid,
                    "Binding_Resname": resname
                })

    for ion, site_pairs in ions_near_both_stable.items():
        for (dna_site, histone_site), count in site_pairs.items():
            rows.append({
                "Region": "Interface",
                "Ion_Resid": ion,
                "Binding_Resid": f"{dna_site[0]}:{dna_site[1]}",
                "Binding_Resname": f"{histone_site[0]}:{histone_site[1]}"
            })

    df = pd.DataFrame(rows)
    DNA_resnames = {'DA', 'DT', 'DG', 'DC'}
    
    
    # ---------------- Apply classification rules ----------------
    classifications = {
        'around_backbone': [],
        'around_gap': [],
        'around_major': [],
        'around_minor': [],
        'around_histone_tail': [],
        'around_histone_core': [],
        'around_acidic_patch': [],
        'around_Interface': []
    }

    for _, row in df.iterrows():
        region = row['Region']
        ion_index = row['Ion_Resid']

        if region != 'Interface':
            binding_info = (row['Binding_Resid'], row['Binding_Resname'])

            if region == 'Backbone':
                if all(ion_index not in df.loc[df['Region'] == r, 'Ion_Resid'].values 
                       for r in ['Major_Groove', 'Minor_Groove', 'Histone_Tail', 'Histone_Core', 'Interface']):
                    classifications['around_backbone'].append((ion_index, binding_info))

            elif region == 'Gap_Region': 
                if all(ion_index not in df.loc[df['Region'] == r, 'Ion_Resid'].values 
                       for r in ['Major_Groove', 'Minor_Groove']):
                    classifications['around_gap'].append((ion_index, binding_info))

            elif region == 'Major_Groove':
                if all(ion_index not in df.loc[df['Region'] == r, 'Ion_Resid'].values 
                       for r in ['Minor_Groove', 'Histone_Tail', 'Histone_Core', 'Interface']):
                    classifications['around_major'].append((ion_index, binding_info))

            elif region == 'Minor_Groove':
                if all(ion_index not in df.loc[df['Region'] == r, 'Ion_Resid'].values 
                       for r in ['Major_Groove', 'Histone_Tail', 'Histone_Core', 'Interface']):
                    classifications['around_minor'].append((ion_index, binding_info))

            elif region == 'Histone_Tail':
                if all(ion_index not in df.loc[df['Region'] == r, 'Ion_Resid'].values 
                       for r in ['Backbone', 'Major_Groove', 'Minor_Groove', 'Interface']):
                    classifications['around_histone_tail'].append((ion_index, binding_info))

            elif region == 'Histone_Core':
                if all(ion_index not in df.loc[df['Region'] == r, 'Ion_Resid'].values 
                       for r in ['Backbone', 'Major_Groove', 'Minor_Groove', 'Interface']):
                    classifications['around_histone_core'].append((ion_index, binding_info))

            elif region == 'Acidic_Patch':
                classifications['around_acidic_patch'].append((ion_index, binding_info))

        else:
 
            try:
                resid1_str, resname1 = str(row['Binding_Resid']).split(":")
                resid2_str, resname2 = str(row['Binding_Resname']).split(":")
                resid1 = int(resid1_str)
                resid2 = int(resid2_str)
            except Exception as e:
                continue

            dna_list, histone_list = [], []
            if resname1 in DNA_resnames:
                dna_list.append((resid1, resname1))
            else:
                histone_list.append((resid1, resname1))
            if resname2 in DNA_resnames:
                dna_list.append((resid2, resname2))
            else:
                histone_list.append((resid2, resname2))

            classifications['around_Interface'].append((ion_index, dna_list, histone_list))

    # ---------------- Output final results ----------------
    output_filename = f"result_{system_name}_run_{run_index + 1}.csv"
    rows_output = []

    for classification, items in classifications.items():
        if classification == 'around_Interface':
            ion_bindings = {}
            for ion_index, dna_list, histone_list in items:
                if ion_index not in ion_bindings:
                    ion_bindings[ion_index] = {'DNA': [], 'Histone': []}
                ion_bindings[ion_index]['DNA'].extend(dna_list)
                ion_bindings[ion_index]['Histone'].extend(histone_list)
            for ion_index, region_bindings in ion_bindings.items():
                dna_str = ", ".join([f"({resid}, {resname})" for resid, resname in region_bindings['DNA']]) if region_bindings['DNA'] else ""
                histone_str = ", ".join([f"({resid}, {resname})" for resid, resname in region_bindings['Histone']]) if region_bindings['Histone'] else ""
                combined_str = ", ".join([s for s in [dna_str, histone_str] if s])  # merge into one string
                rows_output.append({
                    "Classification": classification,
                    "Ion Index": ion_index,
                    "Bindings": combined_str
                })
        else:
            ion_dict = {}
            for ion_index, binding_info in items:
                ion_dict.setdefault(ion_index, []).append(binding_info)
            for ion_index, bindings in ion_dict.items():
                binding_str = ", ".join([f"({resid}, {resname})" for resid, resname in bindings])
                rows_output.append({
                    "Classification": classification,
                    "Ion Index": ion_index,
                    "Bindings": binding_str
                })

    df_output = pd.DataFrame(rows_output)
    df_output.to_csv(output_filename, index=False)

    print(f"Classification results have been saved to {output_filename}")


# ---------------- Batch processing entry ----------------
systems = [
    {
        "name": "CA_15mM_1264_Box_50",
        "topology": "/media/lenovo/Seagate Basic/remove_wat_nc_40ps/CA_15mM_1264_Box_50_rw.prmtop",
        "trajectory_files": [
            "/media/lenovo/Seagate Basic/remove_wat_nc_40ps/CA_15mM_1264_Box_50_run1_rw.nc",
            "/media/lenovo/Seagate Basic/remove_wat_nc_40ps/CA_15mM_1264_Box_50_run2_rw.nc",
            "/media/lenovo/Seagate Basic/remove_wat_nc_40ps/CA_15mM_1264_Box_50_run3_rw.nc"
        ]
    },
]

if __name__ == "__main__":
    for system in systems:
        for run_index, trajectory_file in enumerate(system["trajectory_files"]):
            process_system_and_run(system["name"], system["topology"], trajectory_file, run_index)
