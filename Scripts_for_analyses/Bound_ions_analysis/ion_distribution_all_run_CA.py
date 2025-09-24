import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# ==============================
# 1. Define simulation systems
# ==============================
# Each system includes one topology file (prmtop) and multiple trajectory files (nc)
systems = [
    {
        "name": "CA_15mM_1264_Box_50",
        "topology": "../CA_15mM_1264_Box_50_rw.prmtop",
        "trajectory_files": [
            "../CA_15mM_1264_Box_50_run1_rw.nc",
            "../CA_15mM_1264_Box_50_run2_rw.nc",
            "../CA_15mM_1264_Box_50_run3_rw.nc"
        ]
    },
]

# ==============================
# 2. Define atom selection strings (MDAnalysis syntax)
# ==============================
# DNA minor groove and major groove
minor_groove_selection = "((resname DA) and (name N3 or name H2 )) or ((resname DT) and (name O2 )) or ((resname DG) and (name N3 or name N2)) or ((resname DC) and (name O2))"
major_groove_selection = "((resname DA) and (name N7 or name N6 or name H8)) or ((resname DT ) and (name O4 or name C7)) or ((resname DG) and (name N7 or name O6 or name H8)) or ((resname DC) and (name N4 or name C5))"

# Whole histone & acidic patch
histone_selection = "protein and not name H*"
acidic_patch_selection = "(protein and ((resid 349) or (resid 354) or (resid 357) or (resid 384) or (resid 385) or (resid 386) or (resid 479) or (resid 484) or (resid 487) or (resid 513) or (resid 514) or (resid 515) or (resid 657) or (resid 665) or (resid 782) or (resid 790))) and not name H*"
acidic_residues_selection = "(resname GLU or resname ASP) and not name H* "

# Histone tail & histone core
histone_tail_selection = '(resid 803:839 or resid 938:974 or resid 1073:1092 or resid 1175:1194 or resid 295:307 or resid 412:423 or resid 424:436 or resid 541:552 or resid 553:578 or resid 678:703) and not name H*' 
histone_core_selection = '(protein and not (resid 803:839 or resid 938:974 or resid 1073:1092 or resid 1175:1194 or resid 295:307 or resid 412:423 or resid 424:436 or resid 541:552 or resid 553:578 or resid 678:703)) and not name H*'

# DNA and backbone
dna_selection = "nucleic and not name H*"
dna_backbone_selection = "(resname DA DT DG DC) and (name C3' or name C4' or name C5' or name O3' or name O5' or name P or name OP1 or name OP2)"
upper_backbone_selection = "(resid 1:68 or resid 227:294) and (name C3' or name C4' or name C5' or name O3' or name O5' or name P or name OP1 or name OP2)"  
lower_backbone_selection = "(resid 148:215 or resid 80:147) and (name C3' or name C4' or name C5' or name O3' or name O5' or name P or name OP1 or name OP2)"  

# ==============================
# 3. Define distance cutoffs (Å)
# ==============================
cutoff_minor_groove = 5
cutoff_major_groove = 5
cutoff_dna_backbone = 5
cutoff_histone = 5
cutoff_dna = 5
cutoff_gap = 10.0

# ==============================
# 4. Loop over all systems and trajectories
# ==============================
for system_index, system in enumerate(systems):
    system_name = system["name"]
    topology_file = system["topology"]
    trajectory_files = system["trajectory_files"]

    # Loop over different runs
    for run_index, traj_file in enumerate(trajectory_files):
        u = mda.Universe(topology_file, traj_file)

        # ---- Define ions and regions of interest ----
        ca_ions = u.select_atoms("resname CA")
        minor_groove = u.select_atoms(minor_groove_selection)
        major_groove = u.select_atoms(major_groove_selection)
        dna_backbone = u.select_atoms(dna_backbone_selection)
        dna = u.select_atoms(dna_selection)
        histone_tail = u.select_atoms(histone_tail_selection)
        histone_core = u.select_atoms(histone_core_selection)
        acidic_patch = u.select_atoms(acidic_patch_selection)
        histone = u.select_atoms(histone_selection)
        acidic_residues = u.select_atoms(acidic_residues_selection)
        dna_upper_backbone = u.select_atoms(upper_backbone_selection)
        dna_lower_backbone = u.select_atoms(lower_backbone_selection)

        # ==============================
        # 5. Define helper functions
        # ==============================

        # Compute minimum ion–region distance
        def calculate_distances(ions, region):
            if len(region) == 0:
                return np.array([])  
            distances_array = distances.distance_array(ions.positions, region.positions)
            return distances_array.min(axis=1)  

        # Count functions (exclusive binding to certain regions)
        def count_ions_in_only_minor_groove(distances_minor_groove, distances_histone):
            if distances_minor_groove.size == 0 or distances_histone.size == 0:
                return 0
            return np.sum((distances_minor_groove < cutoff_minor_groove) & (distances_histone >= cutoff_histone))

        def count_ions_in_only_major_groove(distances_major_groove, distances_histone):
            if distances_major_groove.size == 0 or distances_histone.size == 0:
                return 0
            return np.sum((distances_major_groove < cutoff_major_groove) & (distances_histone >= cutoff_histone))

        def count_ions_in_only_dna_backbone(distances_minor_groove, distances_major_groove, distances_dna_backbone, distances_histone):
            if distances_major_groove.size == 0 or distances_minor_groove.size == 0 or distances_dna_backbone.size == 0 or distances_histone.size == 0:
                return 0
            return np.sum((distances_major_groove >= cutoff_major_groove) & 
                          (distances_minor_groove >= cutoff_minor_groove) & 
                          (distances_dna_backbone < cutoff_dna_backbone) & 
                          (distances_histone >= cutoff_histone)) 

        def count_ions_in_only_dna(distances_dna, distances_histone):
            if distances_dna.size == 0 or distances_histone.size == 0:
                return 0
            return np.sum((distances_dna < cutoff_dna) & (distances_histone >= cutoff_histone))

        def count_ions_in_only_histone_tail(distances_histone_tail, distances_histone_core, distances_dna):
            if distances_dna.size == 0 or distances_histone_tail.size == 0:
                return 0
            return np.sum((distances_dna >= cutoff_dna) & (distances_histone_tail < cutoff_histone))

        def count_ions_in_only_histone_core(distances_histone_tail, distances_histone_core, distances_dna):
            if distances_dna.size == 0 or distances_histone_core.size == 0:
                return 0
            return np.sum((distances_dna >= cutoff_dna) & (distances_histone_core < cutoff_histone))

        def count_ions_in_only_acidic_patch(distances_acidic_patch, distances_dna):
            if distances_dna.size == 0 or distances_acidic_patch.size == 0:
                return 0
            return np.sum((distances_dna >= cutoff_dna) & (distances_acidic_patch < cutoff_histone)) 

        def count_ions_in_only_acidic_residues(distances_acidic_residues, distances_dna):
            if distances_dna.size == 0 or distances_acidic_residues.size == 0:
                return 0
            return np.sum((distances_dna >= cutoff_dna) & (distances_acidic_residues < cutoff_histone)) 

        def count_ions_in_only_histone(distances_histone, distances_dna):
            if distances_dna.size == 0 or distances_histone.size == 0:
                return 0
            return np.sum((distances_dna >= cutoff_dna) & (distances_histone < cutoff_histone))

        def count_ions_in_interface(distances_dna, distances_histone):
            if distances_dna.size == 0 or distances_histone.size == 0:
                return 0
            return np.sum((distances_dna < cutoff_dna) & (distances_histone < cutoff_histone))

        def count_ions_in_only_gap(distances_major_groove, distances_minor_groove, distances_upper_backbone, distances_lower_backbone, distances_histone):
            if distances_upper_backbone.size == 0 or distances_lower_backbone.size == 0:
                return 0
            return np.sum((distances_major_groove >= cutoff_major_groove) & 
                          (distances_minor_groove >= cutoff_minor_groove) & 
                          (distances_upper_backbone < cutoff_gap) & 
                          (distances_lower_backbone < cutoff_gap))

        def count_binding_ions(distances_dna, distances_histone):
            """Bound to at least DNA or histone"""
            if distances_dna.size == 0 and distances_histone.size == 0:
                return 0
            return np.sum((distances_dna < cutoff_dna) | (distances_histone < cutoff_histone))

        def count_unbinding_ions(distances_dna, distances_histone):
            """Not bound to either DNA or histone"""
            if distances_dna.size == 0 or distances_histone.size == 0:
                return 0
            return np.sum((distances_dna >= cutoff_dna) & (distances_histone >= cutoff_histone))

        # ==============================
        # 6. Write output files
        # ==============================
        output_file = f"ion_distribution_{system_name}_run_{run_index + 1}.csv"
        with open(output_file, "w") as file:
            file.write("Frame,only_minor_groove,only_major_groove,only_dna_backbone,only_dna,only_gap,only_histone_tail,only_histone_core,only_acidic_patch,only_acidic_residues,only_histone,interface,binding_ions,unbinding_ions\n")

            # Loop through trajectory frames (e.g., 2500–12500, step=25)
            for ts in u.trajectory[2500:12500:25]:
                print(f"Processing frame {ts.frame} of run {run_index + 1} in {system_name}...")

                # ---- Compute distances ----
                distances_minor_groove = calculate_distances(ca_ions, minor_groove)
                distances_major_groove = calculate_distances(ca_ions, major_groove)
                distances_dna_backbone = calculate_distances(ca_ions, dna_backbone)
                distances_dna = calculate_distances(ca_ions, dna)
                distances_histone_tail = calculate_distances(ca_ions, histone_tail)
                distances_histone_core = calculate_distances(ca_ions, histone_core)
                distances_acidic_patch = calculate_distances(ca_ions, acidic_patch)
                distances_acidic_residues = calculate_distances(ca_ions, acidic_residues)
                distances_histone = calculate_distances(ca_ions, histone)
                distances_upper_backbone = calculate_distances(ca_ions, dna_upper_backbone)
                distances_lower_backbone = calculate_distances(ca_ions, dna_lower_backbone)

                # ---- Count ions in different regions ----
                ions_in_only_minor_groove = count_ions_in_only_minor_groove(distances_minor_groove, distances_histone)
                ions_in_only_major_groove = count_ions_in_only_major_groove(distances_major_groove, distances_histone)
                ions_in_only_dna_backbone = count_ions_in_only_dna_backbone(distances_minor_groove, distances_major_groove, distances_dna_backbone, distances_histone)
                ions_in_only_dna = count_ions_in_only_dna(distances_dna, distances_histone)
                ions_in_only_gap = count_ions_in_only_gap(distances_major_groove, distances_minor_groove, distances_upper_backbone, distances_lower_backbone, distances_histone)
                ions_in_only_histone_tail = count_ions_in_only_histone_tail(distances_histone_tail, distances_histone_core, distances_dna)
                ions_in_only_histone_core = count_ions_in_only_histone_core(distances_histone_tail, distances_histone_core, distances_dna)
                ions_in_only_acidic_patch = count_ions_in_only_acidic_patch(distances_acidic_patch, distances_dna)
                ions_in_only_acidic_residues = count_ions_in_only_acidic_residues(distances_acidic_residues, distances_dna)
                ions_in_only_histone = count_ions_in_only_histone(distances_histone, distances_dna)
                ions_in_interface = count_ions_in_interface(distances_dna, distances_histone)
                binding_ions = count_binding_ions(distances_dna, distances_histone)
                unbinding_ions = count_unbinding_ions(distances_dna, distances_histone)
              
                # ---- Write results in CSV ----
                file.write(f"{ts.frame},{ions_in_only_minor_groove},{ions_in_only_major_groove},{ions_in_only_dna_backbone},{ions_in_only_dna},{ions_in_only_gap},{ions_in_only_histone_tail},{ions_in_only_histone_core},{ions_in_only_acidic_patch},{ions_in_only_acidic_residues},{ions_in_only_histone},{ions_in_interface},{binding_ions},{unbinding_ions}\n")

        print(f"Results for run {run_index + 1} in {system_name} written to {output_file}")

print("All runs processed.")


