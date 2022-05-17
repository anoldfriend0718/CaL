import os 
import pandas as pd
from pyCaLProblem import CaLProblem
from pyCaLPlant import Pinch_point_analyzer
import geatpy as ea
import numpy as np
import shutil

if __name__ == '__main__':
    # 实例化问题对象
    # Tcarbs=[600,630,660]
    # cao_conversions=[0.2,0.3,0.4,0.5]
    # T_water_reactor_out=[80,150]

    Tcarbs=[660]
    cao_conversions=[0.3]
    T_water_reactor_out=[80]
    flue_gas_composition=dict()
    flue_gas_composition["co2"]=0.1338
    flue_gas_composition["o2"]=0.0384
    flue_gas_composition["n2"]=0.6975
    vol_rate_flue_gas=6000/3600
    decarbonized_rate=0.9
    T_flue_gas=40
    
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
            for Tw in T_water_reactor_out:
                print(f"solving the case: Tcarb = {Tcarb}, Tw = {Tw}, X = {X}")
                case_id=f"Tcarb={Tcarb}_Tw={Tw}_X={X}"
                print(f"case_id: {case_id}")
                parameters={}
                parameters["T_carb"]=Tcarb
                parameters["cao_conversion"]=X
                parameters["flue_gas_composition"]=flue_gas_composition
                parameters["vol_rate_flue_gas"]=vol_rate_flue_gas
                parameters["decarbonized_rate"]=decarbonized_rate
                parameters["T_flue_gas"]=T_flue_gas
                parameters["T_water_reactor_out"]=Tw
                
                problem = CaLProblem(parameters)
                algorithm = ea.soea_DE_currentToBest_1_bin_templet(problem,
                                                ea.Population(Encoding='RI', NIND=30),
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
                input["T_flue_gas_reactor_in"]=best_vars[0,0]
                input["T_cao_reactor_in"]=best_vars[0,1]
                input["T_water_reactor_in"]=best_vars[0,2]
                plant_results=problem._plant.solve(input)
                print(f"plant results: {plant_results}")
                all_plant_results[case_id]=plant_results

                # save csv files for pinch point analysis
                ppa=Pinch_point_analyzer(plant_results)
                ppa.write_pyPinch_data_csv(f"{pinch_point_analysis_folder}/{case_id}.csv")
            
    df_plant_results=pd.DataFrame(all_plant_results).T
    df_plant_results.to_csv("plant_results.csv")

    df_best_variables = pd.DataFrame(all_best_variables, index=["T_flue_gas_reactor_in",
                                                                "T_cao_reactor_in", "T_water_reactor_in"]).T
    df_best_variables.to_csv("best_variables.csv")
    print("Succeed to output the results")

