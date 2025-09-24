import MDAnalysis as mda
import numpy as np
from math import sqrt
from datetime import datetime
import os

# 开始时间记录
startTime = datetime.now()

# 定义体系和对应的 run 信息
systems = ["NaCl_1264_wo", "MG_15mM_1264_wo"]
runs = ["run1", "run2", "run3"]

# 定义不同SHL的残基范围
shl_ranges = [


    ((64, 66), (229, 231), (142, 144), (151, 153)),  # SHL -1 和 SHL +7
    ((54, 56), (239, 241), (132, 134), (161, 163)),# SHL -2 和 SHL +6
    ((44, 46), (249, 251), (122, 124), (171, 173)),  # SHL -3 和 SHL +5
    ((34, 36), (259, 261), (112, 114), (181, 183)),  # SHL -4 和 SHL +4  
    ((24, 26), (269, 271), (102, 104), (191, 193)),  # SHL -5 和 SHL +3
    ((14, 16), (279, 281), (92, 94), (201, 203)),  # SHL -6 和 SHL +2
    ((4, 6), (289, 291), (82, 84), (211, 213)),  # SHL -7 和 SHL +1
    
   
]

# 获取轨迹帧范围
firstframe, lastframe, step = 2000, 100000, 10

# 创建输出目录
output_dir = "distance_results"
os.makedirs(output_dir, exist_ok=True)

# 遍历每个体系和对应的 run
for system in systems:
    for run in runs:
        print(f"Processing {system}, {run}...")

        # 文件路径
        topology_file = f"../{system}_rw.prmtop"
        trajectory_file = f"../{system}_{run}_rw.nc"

        try:
            # 加载拓扑和轨迹
            u = mda.Universe(topology_file, trajectory_file)

            # 遍历每个 SHL 配对
            for idx, (res1, res2, res3, res4) in enumerate(shl_ranges):
                outfile = os.path.join(output_dir, f"{system}_{run}_shl{idx+1}_gap_along_time.dat")

                # 打开输出文件
                with open(outfile, 'w') as f:
                    for ts in u.trajectory[firstframe:lastframe:step]:  # 遍历指定范围的轨迹
                        # 定义选择
                        selSeg1 = u.select_atoms(
                            f" (resid {res1[0]}:{res1[1]} and not name H*) or "
                            f" (resid {res2[0]}:{res2[1]} and not name H*)"
                        )
                        selSeg2 = u.select_atoms(
                            f" (resid {res3[0]}:{res3[1]} and not name H*) or "
                            f" (resid {res4[0]}:{res4[1]} and not name H*)"
                        )

                        # 计算中心
                        if selSeg1.n_atoms > 0 and selSeg2.n_atoms > 0:
                            center1 = selSeg1.center_of_mass()
                            center2 = selSeg2.center_of_mass()
                            dist = sqrt((center1[0] - center2[0])**2 +
                                        (center1[1] - center2[1])**2 +
                                        (center1[2] - center2[2])**2)
                            f.write(f"{dist:.4f}\n")
                        else:
                            f.write("NaN\n")  # 如果选择为空，写入 NaN

                print(f"Output saved to {outfile}")
        except Exception as e:
            print(f"Error processing {system}, {run}: {e}")

# 输出完成信息
print(f"所有体系计算完成，总耗时: {datetime.now() - startTime}")



