import os 
import pandas as pd
from pyCaesProblem import CaesProblem  # 导入自定义问题接口
from pyHENPinch import Hen_pinch_analyzer

import geatpy as ea  # import geatpy
import numpy as np
import shutil

if __name__ == '__main__':

    flue_gas_composition=dict()
    flue_gas_composition["co2"]=0.1338
    flue_gas_composition["o2"]=0.0384
    flue_gas_composition["n2"]=0.6975

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
    print(f"solving the case: PM = {1}, PL = {1}, T = {1}")
    case_id=f"PM={1}_PL={1}_T={1}"
    print(f"case_id: {case_id}")
    parameters={}
    parameters["flue_gas_composition"]=flue_gas_composition
    parameters["obj"] = "power"
    problem = CaesProblem(parameters)
    algorithm = ea.soea_DE_currentToBest_1_bin_templet(problem,
                                    ea.Population(Encoding='RI', NIND=8),
                                    MAXGEN=150,  # 最大进化代数。
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
    inputs={}
    inputs["p_bray_H"] = best_vars[0,0]#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = best_vars[0,1] #优化变量2，热泵循环中间压力
    inputs["Dehy_caoh2_in"]=best_vars[0,2]
    inputs["p_Dehy"] = best_vars[0,3] #变量4，反应器压力
    inputs["Dehy_overheating_temperature"] = best_vars[0,4] #变量2，脱水反应器过热温度

    inputs["p_bray_H_B"] = best_vars[0,5]
    inputs["p_bray_MH_B"] = best_vars[0,6]
    inputs["p_bray_ML_B"] = best_vars[0,7]
    inputs["p_Hydr"] = best_vars[0,8]
    inputs["Hydr_overheating_temperature"] = best_vars[0,9]

    inputs["Hydr_cao_in"]=best_vars[0,10]

    inputs["cn_m1"] = best_vars[0,11]
    inputs["m2"] = best_vars[0,12]#补充烟气
    inputs["T_X"] = 1
    inputs["m1"] = 0#供暖水

    plant_results=problem.solve(inputs)
    print(f"plant results: {plant_results}")
    print(plant_results["Dehy"]["pinch_analysis_text"])
    print(plant_results["Hydr"]["HEN"]["pinch_analysis_text"])