import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")

import numpy as np
import shutil 
import json
import copy
from pyEconomicComparer import economic_comparer
from pyCaLAnalysis import CaLAnalyser

def save(results,folder,key,value):
    dst=os.path.join(folder,key)
    if not os.path.exists(dst):
        os.makedirs(dst)
    if isinstance(value,int):
        file_name=f"{value}.json"
    else:
        file_name=f"{round(value,3)}.json"
    file_path=os.path.join(dst,file_name)
    if os.path.exists(file_path):
            os.remove(file_path)
    with open(file_path,"w") as fp:
        json.dump(results,fp,indent=4)

def run(parameters,results_folder,key_parameter,value):
    flue_gas_composition={}
    flue_gas_composition["co2"]=0.1338
    flue_gas_composition["o2"]=0.0384
    flue_gas_composition["n2"]=0.6975 
    design_inputs={}
    design_inputs["flue_gas_composition"]=flue_gas_composition       
    design_inputs["user_heat_load"]=15.5e6 #15.5MW
    design_inputs["decarbonized_rate"]=0.9
    design_inputs["T_flue_gas"]=40
    design_inputs["T_carb"]=650
    design_inputs["T_calc"]=900
    design_inputs["calciner_thermal_loss"]=parameters["calciner_thermal_loss"] ### S10
    design_inputs["carbonator_thermal_loss"]=0.01
    design_inputs["cao_conversion"]=parameters["cao_conversion"] ### S1
    design_inputs["T_water_prod_out"]=parameters["T_water_prod_out"]  ### S2
    design_inputs["T_water_supply_in"]=design_inputs["T_water_prod_out"]-25
    design_inputs["storage_carbonator_distance"]=parameters["storage_carbonator_distance"] ### S3
    design_inputs["delta_T_pinch"]=parameters["delta_T_pinch"] ### S4
    design_inputs["calciner_capacity_factor"]=parameters["calciner_capacity_factor"]  ### S8
    design_inputs["CaCO3_storage_duration_hours"]=24*7
    design_inputs["HTCW"]=0
    design_inputs["HRCP"]=1
    design_inputs["obj"]="energy" 
    # design_inputs["obj"]="economic"  

    economic_inputs={}
    economic_inputs["limestone_price"]=70
    economic_inputs["calciner_cost_factor"]=parameters["calciner_cost_factor"]  ## S5
    economic_inputs["elec_price"]=parameters["elec_price"] #元/千瓦时  ##S6
    economic_inputs["operation_hours"]=parameters["operation_hours"]  # hours  ### S9
    economic_inputs["discount_ratio"]=6/100 #8%
    economic_inputs["operational_years"]=30
    economic_inputs["operation_labour_cost_indictor"]=0.025
    economic_inputs["maintain_cost_indictor"]=parameters["maintain_cost_indictor"]  ##S7

    ca=CaLAnalyser()
    results={}
    plant_results=ca.solve(design_inputs)
    results["plant"]=plant_results
    analysis_results=ca.analyze(plant_results,economic_inputs)
    results["metrics"]=analysis_results
    comparor=economic_comparer(economic_inputs)
    results["comparison"]=comparor.compare(results)

    save(results,results_folder,key_parameter,value)



if __name__=="__main__":
    base_parameters={
        "cao_conversion":0.3,
        "T_water_prod_out":85,
        "storage_carbonator_distance":100,
        "delta_T_pinch":15,
        "calciner_cost_factor":1,
        "elec_price":0.165,
        "maintain_cost_indictor":0.025,
        "calciner_capacity_factor":1,
        "operation_hours":3435,
        "calciner_thermal_loss":0.06,
    }
    results_folder="/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL/src/CaL-CC-HS/results/revised/SA"
    variables={
        "cao_conversion":[0.1,0.15,0.2,0.3,0.4,0.5],
        "T_water_prod_out":[80,85,90,100,110],
        "storage_carbonator_distance":[50,100,200,300],
        "delta_T_pinch":[10,15,20,25],
        "calciner_cost_factor":np.array([0.8,0.9,1,1.1,1.2])*base_parameters["calciner_cost_factor"],
        "elec_price":np.array([0.8,0.9,1,1.1,1.2])*base_parameters["elec_price"],
        # "maintain_cost_indictor":[0.015,0.025,0.035,0.045],
        "calciner_capacity_factor":[1,1.5,2,2.5,3],
        "operation_hours":[2000,3000,4000,5000,6000],
        "calciner_thermal_loss":[0.04,0.06,0.08,0.10]
    }

    for key_parameter in variables.keys():
        for value in variables[key_parameter]:
            print(f"solve parameter: {key_parameter} with value: {value} ")
            parameters=copy.deepcopy(base_parameters)
            parameters[key_parameter]=value
            run(parameters,results_folder,key_parameter,value)
    # save(base_parameters,results_folder,"tmp",3)
    
    






