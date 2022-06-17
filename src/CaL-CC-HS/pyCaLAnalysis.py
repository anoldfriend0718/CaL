import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")
from pyCostEstimator import Cost_Estimator

import math
import pandas as pd
import geatpy as ea
import shutil
from scipy.optimize import minimize_scalar
from pyCaLPlant import CalcinerSide,CarbonatorSide,Pinch_point_analyzer
from pyCaLProblem import CarbProblem,CalcProblem
import CoolProp.CoolProp as CP


class CaLAnalyser(object):
    def __init__(self,parameters):
        self._parameters=parameters
        self._calc=CalcProblem(parameters)
        self._carb=CarbProblem(parameters)
        self._cost_estimator=Cost_Estimator()

    def solve(self):
        # carbonator side 
        #optimize the carbonator side according to the energy efficiency
        algorithm = ea.soea_DE_currentToBest_1_bin_templet(self._carb,
                                                ea.Population(Encoding='RI', NIND=30),
                                                MAXGEN=200,  # 最大进化代数。
                                                logTras=1, #,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                                trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                                maxTrappedCount=100)  # # 进化停滞计数器最大上限值。
        algorithm.mutOper.F=0.95 #变异概率
        algorithm.recOper.XOVR = 0.95  # 重组概率
        res = ea.optimize(algorithm, seed=1, verbose=True, drawing=1, outputMsg=True, drawLog=True, saveFlag=False)
        print(f"carbonator optimation results: {res}")
        best_vars=res["Vars"]
        # solve the carbonator results with the best choice   
        plant_results={}
        carb_opt_results={}
        carb_opt_results["T_flue_gas_reactor_in"]=best_vars[0,0]
        carb_opt_results["T_cao_reactor_in"]=best_vars[0,1]
        carb_opt_results["T_water_reactor_in"]=best_vars[0,2]
        carb_results=self._carb.solve(carb_opt_results)
        plant_results["carb"]=carb_results

        #calciner side 
        mass_camix_calc=carb_results["m_cao_unr_out"]+carb_results["m_caco3_out"]
        calc_opt_results=minimize_scalar(self._calc.opt,args=(mass_camix_calc,),
            method='bounded',bounds=(0, 1),tol=1e-5)
        calc_mfrac=calc_opt_results.x
        calc_inputs={}
        calc_inputs["mass_camix_in"]=mass_camix_calc
        calc_inputs["mfrac"]=calc_mfrac
        calc_results=self._calc.solve(calc_inputs)
        plant_results["calc"]=calc_results

        return plant_results

    def analyze(self,plant_results):
        results={}
        results["economic"]=self.analyze_economic_metrics(plant_results)
        results["energy"]=self.analyze_energy_metrics(plant_results)
        return results

    def analyze_economic_metrics(self,plant_results):
        invCosts=self._cost_estimator.solve(plant_results)
        return invCosts
    
    def analyze_energy_metrics(self,plant_results):
        analysis_results={}
        # thermal storage rate
        total_auxiliary_power=plant_results["carb"]["carb_auxiliary_power"]+\
            plant_results["calc"]["calc_auxiliary_power"]
        total_power=plant_results["calc"]["We_calc"]-total_auxiliary_power
        analysis_results["total_power"]=total_power
        analysis_results["Q_hot_water"]=plant_results["carb"]["Q_hot_water"]
        analysis_results["heat_storage_eff"]=analysis_results["Q_hot_water"]/total_power
        
        #seperation effeciency
        vol_rate_flue_gas=plant_results["carb"]["vol_rate_flue_gas"]
        flue_gas_composition = plant_results["carb"]["flue_gas_composition"]
        co2_percent=flue_gas_composition["co2"]

        T_flue_gas_in = plant_results["carb"]["T_flue_gas_in"]
        mass_co2_i =vol_rate_flue_gas*co2_percent *\
            CP.PropsSI('D', 'T', T_flue_gas_in+273.15, 
                        'P', plant_results["carb"]["p_flue_gas_in"], "REFPROP::co2")
        M_CO2=44e-3
        mole_co2_in=mass_co2_i/M_CO2
        deconbonized_rate = plant_results["carb"]["deconbonized_rate"]
        specific_sep_power=self._specific_seperation_power(T_flue_gas_in,flue_gas_composition,deconbonized_rate)
        total_sep_power=mole_co2_in*deconbonized_rate*specific_sep_power
        analysis_results["total_sep_power"]=total_sep_power
        T_amb=plant_results["carb"]["T_amb"]+273.15
        T_hot_water=plant_results["carb"]["T_water_reactor_out"]+273.15
        sep_eff1=total_sep_power/(total_power-analysis_results["Q_hot_water"]*(1-T_amb/T_hot_water))
        #not include the compression power for co2 storage
        sep_eff2=total_sep_power/(total_power+plant_results["calc"]["CO2_compression_train"]["compressor_power"]
                -analysis_results["Q_hot_water"]*(1-T_amb/T_hot_water))
        analysis_results["sep_eff1"]=sep_eff1
        analysis_results["sep_eff2"]=sep_eff2
        
        
        return analysis_results

    # J/mol (CO2)
    def _specific_seperation_power(self,T,comps,deconbonized_rate):
        R=8.3145
        T=T+273.15
        gibbs=0
        ns=1
        for i in comps.keys():
            yi=comps[i]
            gibbs=gibbs+R*T*ns*yi*math.log(yi)
        ne=1-comps["co2"]*deconbonized_rate
        compe={}
        for i in comps.keys():
            if i=="co2":
                compe["co2"]=comps["co2"]*(1-deconbonized_rate)/ne
            else:
                compe[i]=comps[i]/ne
        gibbe=0
        for i in compe.keys():
            yi=compe[i]
            gibbe=gibbe+R*T*ne*yi*math.log(yi)
        sep_power=gibbe-gibbs
        specific_sep_power=sep_power/(deconbonized_rate*comps["co2"])
        return specific_sep_power
        
        
if __name__=="__main__":
        flue_gas_composition={}
        flue_gas_composition["co2"]=0.1338
        flue_gas_composition["o2"]=0.0384
        flue_gas_composition["n2"]=0.6975 
        parameters={}
        parameters["flue_gas_composition"]=flue_gas_composition
        parameters["vol_rate_flue_gas"]=6000/3600
        parameters["decarbonized_rate"]=0.9
        parameters["T_flue_gas"]=40
        parameters["T_carb"]=650
        parameters["T_calc"]=900
        parameters["cao_conversion"]=0.5
        parameters["T_water_reactor_out"]=80

        ca=CaLAnalyser(parameters)
        results={}
        plant_results=ca.solve()
        results["plant"]=plant_results
        analysis_results=ca.analyze(plant_results)
        results["metrics"]=analysis_results
        print(results)