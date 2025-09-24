#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 14:48:45 2019

@author: pengy10
"""

import pandas as pd
#"NaCl_12_6","CA_20mM_126_4","CA_26mM_126_4","MG_12mM_126_4","MG_17mM_126_4"
# Different sets of simulation runs
Model_list = ["all_NaCl_Box_30","all_CA_15mM_1264_Box_30","all_CA_15mM_1264_Box_50","all_NaCl_Box_50","all_MG_15mM_1264_Box_30","all_MG_15mM_1264_Box_50",]
for file_name in Model_list:
    # Unbound state is defined if the percentage of tail residues maintaining contacts with the DNA molecule is no more than the “cut-off”
    for cut_off_ratio in [0, 0.2]:
        tail_residence_time = {i: [] for i in ["H2A_N", "H2A_C", "H2B", "H3", "H4"]}

        # Read the measurements of tail-DNA contacts from simulations  
        H2A_contacts1 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h2a_C.dat", sep="\t", header=0)
        H2A_contacts2 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h2a_G.dat", sep="\t", header=0)
        H2B_contacts1 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h2b_D.dat", sep="\t", header=0)
        H2B_contacts2 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h2b_H.dat", sep="\t", header=0)

        H3_contacts1 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h3_A.dat", sep="\t", header=0)
        H3_contacts2 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h3_E.dat", sep="\t", header=0)
        H4_contacts1 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h4_B.dat", sep="\t", header=0)
        H4_contacts2 = pd.read_csv("../tail_dna/tail_DNA_mean_contacts_" + file_name + "_h4_F.dat", sep="\t", header=0)

        # Combine measurements from the long runs for calculations
        H2A_contacts = pd.concat([H2A_contacts1.iloc[0:1203, ], H2A_contacts2.iloc[0:1203, ]]).reset_index().iloc[:, 1:]
        H2B_contacts = pd.concat([H2B_contacts1.iloc[0:1203, ], H2B_contacts2.iloc[0:1203, ]]).reset_index().iloc[:, 1:]
        H3_contacts = pd.concat([H3_contacts1.iloc[0:1203, ], H3_contacts2.iloc[0:1203, ]]).reset_index().iloc[:, 1:]
        H4_contacts = pd.concat([H4_contacts1.iloc[0:1203, ], H4_contacts2.iloc[0:1203, ]]).reset_index().iloc[:, 1:]

        # Time step for frame intervals
        ts = 1  # ns

        # Define run size (e.g., 400 frames per run)
        run_size = 401
        num_runs = len(H2A_contacts) // run_size

        # Function to compute residence times within a single run
        def compute_residence_time(contacts, start_col, end_col, num_residues, tail_key):
            for run_index in range(num_runs):
                start_frame = run_index * run_size
                end_frame = start_frame + run_size

                count_bound_time = 0
                for i in range(start_frame, end_frame - 2):
                    contact0 = filter(lambda x: x >= 1, contacts.iloc[i, start_col:end_col])
                    contact1 = filter(lambda x: x >= 1, contacts.iloc[i + 1, start_col:end_col])

                    ratio0 = len(list(contact0)) / num_residues
                    ratio1 = len(list(contact1)) / num_residues

                    if ratio0 > cut_off_ratio:
                        count_bound_time += 1
                        if i == end_frame - 3:
                            tail_residence_time[tail_key].append(count_bound_time * ts)
                            count_bound_time = 0
                        elif ratio1 <= cut_off_ratio:
                            tail_residence_time[tail_key].append(count_bound_time * ts)
                            count_bound_time = 0
                    elif ratio0 <= cut_off_ratio:
                        count_bound_time = 0

        # Compute residence times for all tails
        compute_residence_time(H2A_contacts, 1, 14, 13, "H2A_N")
        compute_residence_time(H2A_contacts, 14, 26, 12, "H2A_C")
        compute_residence_time(H2B_contacts, 1, 27, 26, "H2B")
        compute_residence_time(H3_contacts, 1, 38, 37, "H3")
        compute_residence_time(H4_contacts, 1, 21, 20, "H4")

        # Convert results into a DataFrame
        df_residence_time = pd.DataFrame({"residence_time": tail_residence_time["H2A_N"],
                                          "tail_type": ["H2A_N"] * len(tail_residence_time["H2A_N"])})
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H2A_C"],
                                                                   "tail_type": ["H2A_C"] * len(tail_residence_time["H2A_C"])}))
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H2B"],
                                                                   "tail_type": ["H2B"] * len(tail_residence_time["H2B"])}))
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H3"],
                                                                   "tail_type": ["H3"] * len(tail_residence_time["H3"])}))
        df_residence_time = df_residence_time.append(pd.DataFrame({"residence_time": tail_residence_time["H4"],
                                                                   "tail_type": ["H4"] * len(tail_residence_time["H4"])}))

        # Filter residence times < 10 ns
        cut_off_time = 10
        df_residence_time_filter = df_residence_time.loc[df_residence_time['residence_time'] >= cut_off_time]

        # Write filtered residence times to a file
        with open(file_name + "_full_tail_residence_time" + str(cut_off_ratio) + ".csv", "w") as fwh:
            for i in range(0, len(df_residence_time_filter)):
                fwh.write(str(df_residence_time_filter.reset_index().iloc[i, 2]) + "," +
                          str(df_residence_time_filter.reset_index().iloc[i, 1]) + "\n")

        
                
        
        
                
        
