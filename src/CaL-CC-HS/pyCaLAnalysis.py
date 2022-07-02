import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")


import math
import copy
import pandas as pd
from scipy.optimize import minimize_scalar,brentq
import geatpy as ea
import CoolProp.CoolProp as CP

from pyCaLProblem import CarbProblem,CalcProblem
from pyCostEstimator import Cost_Estimator
from pyEconomicComparer import economic_comparer


class CaLAnalyser(object):
    def __init__(self):
        self._tmp_carb_results={}
     
    def solve(self,parameters):
        # carbonator side
        flue_gas_rate_s=0.1
        flue_gas_rate_e=20 #TODO: hardcode, depending on the user hot load
        target_heat_load=parameters["user_heat_load"]
        # xtol=0.0001 #m3/s
        rtol=0.1/100 #0.1%    
        maxiter=20
        vol_flue_gas_rate=brentq(self._solve_carb_side,flue_gas_rate_s,flue_gas_rate_e,
            args=(target_heat_load,parameters),rtol=rtol,maxiter=maxiter)
        print(f"calculated flue gas volumetric rate: {vol_flue_gas_rate}")
        carb_results=self._tmp_carb_results

        #calciner side 
        calc_side=CalcProblem(parameters)
        mass_camix_calc=carb_results["m_cao_unr_out"]+carb_results["m_caco3_out"]
        calc_opt_results=minimize_scalar(calc_side.opt,args=(mass_camix_calc,),
            method='bounded',bounds=(0, 1),tol=1e-5)
        calc_mfrac=calc_opt_results.x
        calc_inputs={}
        calc_inputs["mass_camix_in"]=mass_camix_calc
        calc_inputs["mfrac"]=calc_mfrac
        calc_results=calc_side.solve(calc_inputs)

        #compose plant results
        plant_results={}
        plant_results["carb"]=carb_results
        plant_results["calc"]=calc_results
        return plant_results

    def _solve_carb_side(self,vol_rate_flue_gas,target_heat_load,other_parameters):
        # carbonator side 
        print(f"solving carb side with flue gas volumetric rate: {vol_rate_flue_gas}")
        carb_parameters = copy.deepcopy(other_parameters)
        carb_parameters["vol_rate_flue_gas"]=vol_rate_flue_gas
        carb_prb=CarbProblem(carb_parameters)
        
        #optimize the carbonator side according to the energy efficiency
        algorithm = ea.soea_DE_currentToBest_1_bin_templet(carb_prb,
                                                ea.Population(Encoding='RI', NIND=30), #它的取值范围一般在[5D，10D]之间（D为每个个体的维度）。
                                                MAXGEN=200,  # 最大进化代数。
                                                logTras=1, #,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                                trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                                maxTrappedCount=100)  # # 进化停滞计数器最大上限值。
         #变异概率 #[0.4，0.95]
         #增大F 可以加大算法的搜索空间，提高种群多样性，有利于算法搜索最优解，但会降低收敛速率。
         #减小F 可以增加算法的开发能力，提高算法的收敛速度，但同时增加陷入早熟收敛的风险。
        algorithm.mutOper.F=0.95
        #交叉概率因子(CR)起着平衡算法全局与局部搜索能力的作用。它的取值范围一般在[0.3, 0.9]之间。
        #增大CR 可以提高种群多样性，但可能会造成算法后期收敛速度变慢。
        #减小CR 有利于分析个体各维可分离问题
        algorithm.recOper.XOVR = 0.9  # 重组概率
        res = ea.optimize(algorithm, seed=1, verbose=True, drawing=1, outputMsg=True, drawLog=True, saveFlag=False)
        print(f"carbonator optimation results: {res}")
        print(f'lastPop:\n {res["lastPop"].Phen}')
        best_vars=res["Vars"]
        # solve the carbonator results with the best choice   
        
        carb_opt_results=self._compose_carb_opt_results(best_vars,carb_parameters)
        carb_results=carb_prb.solve(carb_opt_results)
        self._tmp_carb_results=carb_results
  
        Q_hot_water=carb_results["Q_hot_water"]
        error=Q_hot_water-target_heat_load
        print(f"the calculated user heat load: {Q_hot_water}, the error: {error}")
        return error

    def _compose_carb_opt_results(self,best_vars,carb_parameters):
        carb_opt_results={}
        if carb_parameters["HTCW"]==1 and carb_parameters["HRCP"]==1:
            carb_opt_results["T_flue_gas_reactor_in"]=best_vars[0,0]
            carb_opt_results["T_cao_reactor_in"]=best_vars[0,1]
            carb_opt_results["T_water_reactor_in"]=best_vars[0,2]
        elif carb_parameters["HTCW"]==1 and carb_parameters["HRCP"]==0:
            carb_opt_results["T_flue_gas_reactor_in"]=best_vars[0,0]
            carb_opt_results["T_cao_reactor_in"]=best_vars[0,1]
        elif carb_parameters["HTCW"]==0 and carb_parameters["HRCP"]==1:
            carb_opt_results["T_flue_gas_reactor_in"]=best_vars[0,0]
            carb_opt_results["m_water_in"]=best_vars[0,1]
        else:
            raise ValueError("Invalid HTCW or HRCP!")

        return carb_opt_results
    
    def analyze(self,plant_results,economic_inputs):
        results={}
        results["energy"]=self.analyze_energy_metrics(plant_results)
        results["economic"]=self.analyze_economic_metrics(plant_results,results["energy"],economic_inputs)
        
        return results

    def analyze_economic_metrics(self,plant_results,energy_analysis_results,economic_inputs):
        cost_estimator=Cost_Estimator()
        invCosts=cost_estimator.solve(plant_results,energy_analysis_results,economic_inputs)
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
        T_hot_water=plant_results["carb"]["T_water_prod_out"]+273.15
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
        design_inputs={}
        design_inputs["flue_gas_composition"]=flue_gas_composition
        # parameters["vol_rate_flue_gas"]=6000/3600
        design_inputs["user_heat_load"]=15.5e6 #15.5MW
        design_inputs["decarbonized_rate"]=0.9
        design_inputs["T_flue_gas"]=40
        design_inputs["T_carb"]=650
        design_inputs["T_calc"]=900
        design_inputs["cao_conversion"]=0.5
        design_inputs["T_water_supply_in"]=60
        design_inputs["T_water_prod_out"]=85
        design_inputs["HTCW"]=0
        design_inputs["HRCP"]=1
        # design_inputs["obj"]="energy" 
        design_inputs["obj"]="economic"  

        economic_inputs={}
        economic_inputs["elec_price"]=0.165 #元/千瓦时
        economic_inputs["operation_hours"]=3435 # hours
        economic_inputs["discount_ratio"]=6/100 #8%
        economic_inputs["operational_years"]=30
        economic_inputs["operation_labour_cost_indictor"]=0.025
        economic_inputs["maintain_labour_cost_indictor"]=0.025

        ca=CaLAnalyser()
        results={}
        plant_results=ca.solve(design_inputs)
        results["plant"]=plant_results
        analysis_results=ca.analyze(plant_results,economic_inputs)
        results["metrics"]=analysis_results
        comparor=economic_comparer(economic_inputs)
        results["comparison"]=comparor.compare(results)

        print(results)
