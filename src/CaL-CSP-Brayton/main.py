import os 
import pandas as pd
from pyCaLProblem import CaLProblem
from pyCaLPlant import Pinch_point_analyzer
import geatpy as ea
import numpy as np
import shutil

if __name__ == '__main__':
    # 实例化问题对象
    Tcarbs=[775,825,875]
    cao_conversions=[0.2,0.3,0.4,0.5]
    # Tcarbs=[775]
    # cao_conversions=[0.2,0.3]
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

    for Tcarb in Tcarbs:
        for X in cao_conversions:
            print(f"solving the case: Tcarb = {Tcarb}, X = {X}")
            case_id=f"Tcarb={Tcarb}_X={X}"
            print(f"case_id: {case_id}")
            parameters={}
            parameters["T_reaction"]=Tcarb
            parameters["m_cao_i"]=1.38
            parameters["p_co2_storage"]=75e5
            parameters["cao_conversion"]=X
            problem = CaLProblem(parameters)
            algorithm = ea.soea_DE_currentToBest_1_bin_templet(problem,
                                            ea.Population(Encoding='RI', NIND=60),
                                            MAXGEN=200,  # 最大进化代数。
                                            logTras=1, #,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                            trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                            maxTrappedCount=100)  # # 进化停滞计数器最大上限值。
            algorithm.mutOper.F=0.95 #变异概率
            algorithm.recOper.XOVR = 0.95  # 重组概率

            # 求解
            res = ea.optimize(algorithm, seed=1, verbose=True, drawing=1, outputMsg=True, drawLog=True, saveFlag=False)
            print(res)
            best_vars=res["Vars"]
            all_best_variables[case_id]=best_vars[0]

            # copy the trace plot to folder
            shutil.copyfile("./Trace Plot.svg",f"{convergence_plots_folder}/{case_id}.svg")
            
            # solve the plant status with the best choice   
            input={}
            input["m_cao_i"] = parameters["m_cao_i"]
            input["p_co2_storage"] = parameters["p_co2_storage"]
            input["T_reaction"] = parameters["T_reaction"]
            input["cao_conversion"]=parameters["cao_conversion"]
            input["p_reaction"]=best_vars[0,0]
            input["expansion_ratio"]=best_vars[0,1]
            input["T_co2_reactor_in"]=best_vars[0,2]
            input["T_cao_reactor_in"]=best_vars[0,3]
            input["T_co2_compressor_in"]=best_vars[0,4]
            input["T_co2_storage_turbine_in"]=best_vars[0,5]
            plant_results=problem._plant.solve(input)
            all_plant_results[case_id]=plant_results

            # save csv files for pinch point analysis
            ppa=Pinch_point_analyzer(plant_results)
            ppa.write_pyPinch_data_csv(f"{pinch_point_analysis_folder}/{case_id}.csv")
            
    df_plant_results=pd.DataFrame(all_plant_results).T
    df_plant_results.to_csv("plant_results.csv")

    df_best_variables=pd.DataFrame(all_best_variables,index=["p_reaction","expansion_ratio",
        "T_co2_reactor_in","T_cao_reactor_in","T_co2_compressor_in","T_co2_storage_turbine_in"]).T
    df_best_variables.to_csv("best_variables.csv")
    print("Succeed to output the results")

