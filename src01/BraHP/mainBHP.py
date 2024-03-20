import os 
import pandas as pd
from pyBhpProblem import BhpProblem  # 导入自定义问题接口

import geatpy as ea  # import geatpy
import numpy as np
import shutil

if __name__ == '__main__':

    flue_gas_composition=dict()
    flue_gas_composition["co2"]=0.1338
    flue_gas_composition["o2"]=0.0384
    flue_gas_composition["n2"]=0.6975
    #p_bray_M = 13e6
    p_bray_L = 7.5e6
    t_reaction = 525

    all_plant_results={}
    all_best_variables={}
    pinch_point_analysis_folder="./ppa_results"
    convergence_plots_folder="./convergence_plots"
    folders=[pinch_point_analysis_folder,convergence_plots_folder]
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)
        else:
            shutil.rmtree(folder)
            os.mkdir(folder)
    print(f"solving the case:  PL = {p_bray_L}, T = {t_reaction}")
    case_id=f"PL={p_bray_L}_T={t_reaction}"
    print(f"case_id: {case_id}")
    parameters={}
    parameters["flue_gas_composition"]=flue_gas_composition
    parameters["obj"] = "cop"
    problem = BhpProblem(parameters)
    algorithm = ea.soea_DE_currentToBest_1_bin_templet(problem,
                                    ea.Population(Encoding='RI', NIND=15),
                                    MAXGEN=120,  # 最大进化代数。
                                    logTras=1, #,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                    trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                    maxTrappedCount=100)  # # 进化停滞计数器最大上限值。
    algorithm.mutOper.F=0.95 #变异概率
    algorithm.recOper.XOVR = 0.95  # 重组概率

    # 求解
    res = ea.optimize(algorithm, seed=1, verbose=True, drawing=1, outputMsg=True, drawLog=True, saveFlag=False)
    print(res)
    best_vars=res["Vars"]
    print(best_vars)
    # copy the trace plot to folder
    shutil.copyfile("./Trace Plot.svg",f"{convergence_plots_folder}/{case_id}.svg")
    input={}
    input["p_bray_H"]=best_vars[0,0]
    input["p_bray_M"]=best_vars[0,1]
    input["Store_electrical_power"]=1e6
    plant_results=problem.solve(input)
    print(f"plant results: {plant_results}")
    