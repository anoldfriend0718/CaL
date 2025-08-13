import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")
import math
import json
import numpy as np
from pyHydrator import Hydrator
from scipy.optimize import fsolve
## construct cost indictor
construct_labour_cost_indictor=0.25
engineering_project_cost_indictor=0.175
TASC_multiplier=1.13
piping_integration_cost_indictor=0.05
#单次储能时间8小时
Single_run_time=8*3600

M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol

class Cost_Estimator(object):
    def __init__(self, parameters) -> None:
        self._hy = Hydrator(parameters)
        self._Eur2Dollars={2014:39.9116/30.3287,2018:35.3758/30.1162}
        self._Dollars2RMB={2024:7.2492}
        self._CEPCIs={2001:394.3, 2014:576.1,2017:567.5, 2018:603.1, 2020:596.2, 2022:816, 2024:798.8 }
        #https://toweringskills.com/financial-analysis/cost-indices/
        self._target_year=2024
        
    def solve(self,inputs,economic_inputs):
        self._res = self._hy.solve(inputs)
        #单次储能时间8小时
        Single_run_time=8*3600 #年循环次数（储能2560小时，释能2560小时，年运行时间5120小时 250天2000小时 320天2560小时
        #循环极限次数
        Limit_cycle_number = economic_inputs["caoh2_unit_price"]/ economic_inputs["caoh2_price/life_rate"]
        #安全存储比例
        cost={}
        cost["BH"]= self._res["BH"]["cost"]
        cost["Dehy"]= self._res["Dehy"]["cost"]
        cost["Hydr"]= self._res["Hydr"]["cost"]
        cost["BC"]= self._res["Brayton"]["cost"]
        CaOH2_storage_mass_flow = self._res["Dehy"]["dehydrator"]["in"]["m_camix"]
        CaOH2_storage_once = CaOH2_storage_mass_flow * Single_run_time / 1000#t
        total_power = self._res["Case_All"]["power_in"]
        total_power_out = self._res["Case_All"]["power"]
        total_hot_out = self._res["Case_All"]["hot_out"]
        N = economic_inputs["operational_years"]
        r = economic_inputs["discount_ratio"]

        #建设成本
        construct_costs={}
        equipment_costs=self.calculate_equipment_costs(cost)
        #设备成本
        construct_costs["equipment"]=equipment_costs
        #安装成本
        construct_costs["installation"]=piping_integration_cost_indictor*equipment_costs["total"]
        #劳动力成本
        construct_costs["labour"]=construct_labour_cost_indictor*(equipment_costs["total"]+construct_costs["installation"])
        #工程项目成本
        construct_costs["engineering&project"]=engineering_project_cost_indictor*(equipment_costs["total"]+construct_costs["installation"])
        #初始材料成本
        construct_costs["initial_material"]=CaOH2_storage_once*economic_inputs["caoh2_unit_price"]*1.3 #1.3是安全存储比例
        #实际总支出
        construct_costs["total as-spent"]=TASC_multiplier*(equipment_costs["total"]+construct_costs["installation"]+
            construct_costs["labour"]+construct_costs["engineering&project"])+construct_costs["initial_material"]

        ## 运营成本
        operation_costs={}
        #等效补充成本
        make_up_limestone_percentage=CaOH2_storage_once*(economic_inputs["Annual_cycle_count"]/Limit_cycle_number)
        operation_costs["make-up_limestone"]=make_up_limestone_percentage*economic_inputs["caoh2_unit_price"]
        #劳动力成本
        operation_costs["labour"]=economic_inputs["operation_labour_cost_indictor"]*construct_costs["total as-spent"]
        #维护成本
        operation_costs["maintain"]=economic_inputs["maintain_cost_indictor"]*construct_costs["total as-spent"]
        #电力成本
        charge_power=total_power/1000*economic_inputs["operation_hours"]
        operation_costs["electricity"]=economic_inputs["elec_price"]*(total_power)/1000* \
            economic_inputs["operation_hours"]
        #年运行成本
        operation_costs["total as-spent"]=operation_costs["labour"]+operation_costs["maintain"]+\
            operation_costs["electricity"]+operation_costs["make-up_limestone"]
        #年收益
        Annual_operating_income={}
        dischager_power = (total_power_out)/1000*economic_inputs["operation_hours"]
        Annual_operating_income["electricity"]=economic_inputs["elec_price_h"]*(total_power_out)/1000* \
            economic_inputs["operation_hours"]
        Annual_operating_income["hot"]=economic_inputs["hot_price_h"]*(total_hot_out)/1000* \
            economic_inputs["operation_hours"]
        #平准化度电成本
        LCOE = ((construct_costs["total as-spent"]+sum(((operation_costs["total as-spent"]-Annual_operating_income["hot"])/(1+r)**x) for x in range(1, N+1)))/
                sum((dischager_power/(1+r)**x) for x in range(1, N+1)))
        #内部收益率
        r_initial_guess=0.01
        solution = fsolve(self.calculate_irr, r_initial_guess, args=(construct_costs["total as-spent"],
                                 Annual_operating_income["electricity"]+Annual_operating_income["hot"]-operation_costs["total as-spent"],
                                 N))
        IRR = solution[0]
        #动态回收周期 先检查容量，没有收益不能运行这个
        dpp = self.discounted_payback_period(construct_costs["total as-spent"],
                                            r,
               Annual_operating_income["electricity"]+Annual_operating_income["hot"]-operation_costs["total as-spent"],
               N)
        investment_costs={}
        investment_costs["BH"]= self._res["BH"]
        investment_costs["Dehy"]= self._res["Dehy"]
        investment_costs["Hydr"]= self._res["Hydr"]
        investment_costs["BC"]= self._res["Brayton"]
        investment_costs["construction"]=construct_costs
        investment_costs["operation"]=operation_costs
        investment_costs["dischager_power"]=dischager_power
        investment_costs["Annual_operating_income"]=Annual_operating_income["electricity"]+Annual_operating_income["hot"]
        investment_costs["Annual_operating_profit"]=investment_costs["Annual_operating_income"]-operation_costs["total as-spent"]
        investment_costs["LCOE"]=LCOE
        investment_costs["IRR"]=IRR
        investment_costs["Round-trip"]=self._res["Case_All"]["Round-trip_eff"]
        investment_costs["DPP"]=dpp
    
        return investment_costs
    
    def calculate_material_costs(self,CaOH2_storage_mass_flow):
        CaOH2_storage_mass=CaOH2_storage_mass_flow*Single_run_time
        total_CaOH2_storage_mass_ton=(CaOH2_storage_mass)/1000 #循环量t
        initial_material_cost=caoh2_unit_price*total_CaOH2_storage_mass_ton/1e6
        return initial_material_cost
    
    def calculate_equipment_costs(self,design):
        equipment_costs={}
        equipment_costs.update(self.calculate_bretonHP_costs(design["BH"]))
        equipment_costs.update(self.calculate_dehydrator_costs(design["Dehy"]))
        equipment_costs.update(self.calculate_hydrator_costs(design["Hydr"]))
        equipment_costs.update(self.calculate_breton_costs(design["BC"]))

        equipment_costs["total"]=(equipment_costs["cost_BHall"]+
                                  equipment_costs["cost_Ball"]+
                                  equipment_costs["cost_deall"]+
                                  equipment_costs["cost_hyall"])

        #equipment_costs["total"]=np.sum(list(equipment_costs.values()))
        return equipment_costs
    
    def calculate_bretonHP_costs(self,bretonHP_design):
        invCosts={}  
        invCosts["cost_pc"]=self._cost_bretonHP_comp(bretonHP_design["cost_pc"])
        invCosts["cost_sc"]=self._cost_bretonHP_comp(bretonHP_design["cost_sc"])
        invCosts["cost_t"]=self._cost_bretonHP_turb(bretonHP_design["cost_t"])
        invCosts["cost_pe"] = self._cost_HEN(bretonHP_design["cost_pe"])
        invCosts["cost_se"] = self._cost_HEN(bretonHP_design["cost_se"])
        invCosts["cost_me"] = self._cost_HEN(bretonHP_design["cost_me"])
        invCosts["cost_hr"] = self._cost_HEN(bretonHP_design["cost_hr"])
        invCosts["cost_comp&turbs"] = self._convert_to_RMB_in_target_year(bretonHP_design["cost_comp&turb"],2020)
        invCosts["cost_exchangers"] = self._convert_to_RMB_in_target_year(bretonHP_design["cost_exchanger"],2020)
        invCosts["cost_BHall"] = invCosts["cost_comp&turbs"]+invCosts["cost_exchangers"]
        #invCosts["cost_"] = self._convert_to_RMB_in_target_year(bretonHP_design["cost_"],2020)
        # TODO: solid conveying system
        return invCosts
    
    def calculate_breton_costs(self,bretonHP_design):
        result={}  
        result["cost_Bpt"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_pt"],2020)
        result["cost_Bst"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_st"],2020)
        result["cost_Bpc"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_pc"],2020)
        result["cost_Bsc"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_sc"],2020)
        result["cost_Bpe"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_pe"],2020)
        result["cost_Bse"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_se"],2020)
        result["cost_Bhe"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_he"],2020)
        result["cost_Ble"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_le"],2020)
        result["cost_Bhr"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_hr"],2020)
        result["cost_Bic"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_ic"],2020)
        result["cost_Bct"]=self._convert_to_RMB_in_target_year(bretonHP_design["cost_ct"],2020)
        result["cost_Bgenerator"] =self._convert_to_RMB_in_target_year(bretonHP_design["cost_generator"],2020)

        result["cost_Bcomp&turb"]=result["cost_Bpc"]+result["cost_Bsc"]+result["cost_Bpt"]+result["cost_Bst"]
        result["cost_Bexchanger"]=result["cost_Bpe"]+result["cost_Bse"]+result["cost_Bhe"]+result["cost_Ble"]+result["cost_Bic"]+result["cost_Bct"]+result["cost_Bhr"]
        result["cost_Ball"]=result["cost_Bcomp&turb"]+result["cost_Bexchanger"]+result["cost_Bgenerator"]
        #invCosts["cost_"] = self._convert_to_RMB_in_target_year(bretonHP_design["cost_"],2020)
        # TODO: solid conveying system
        return result
    
    def calculate_dehydrator_costs(self,dehy_design):
        invCosts={}
        invCosts["cost_hen"] = self._convert_to_RMB_in_target_year(dehy_design["cost_hen"],2020)
        invCosts["cost_dehy"] = self._convert_to_RMB_in_target_year(dehy_design["cost_dehy"],2020)
        invCosts["cost_steam"] = self._convert_to_RMB_in_target_year(dehy_design["cost_steam_blow"],2020)
        invCosts["cost_fgf"] = self._convert_to_RMB_in_target_year(dehy_design["cost_flue_gas_fan"],2020)
        invCosts["cost_wp"] = self._convert_to_RMB_in_target_year(dehy_design["cost_water_pump"],2020)
        invCosts["cost_deall"] = invCosts["cost_hen"]+invCosts["cost_dehy"]+invCosts["cost_steam"]+invCosts["cost_fgf"]+invCosts["cost_wp"]
        # TODO: solid conveying system
        return invCosts
    def calculate_hydrator_costs(self,dehy_design):
        invCosts={}
        invCosts["cost_HEN"] = self._convert_to_RMB_in_target_year(dehy_design["cost_hen"],2020)
        invCosts["cost_hydr"] = self._convert_to_RMB_in_target_year(dehy_design["cost_hydr"],2020)
        invCosts["cost_Steam"] = self._convert_to_RMB_in_target_year(dehy_design["cost_steam_blow"],2020)
        invCosts["cost_Fgf"] = self._convert_to_RMB_in_target_year(dehy_design["cost_flue_gas_fan"],2020)
        invCosts["cost_hyall"] = invCosts["cost_HEN"]+invCosts["cost_hydr"]+invCosts["cost_Steam"]+invCosts["cost_Fgf"]
        # TODO: solid conveying system
        return invCosts

    def _cost_fluidized_bed_dehydrator(self,Qc):
        #reference:https://doi.org/10.1016/j.ijggc.2019.01.005
        #the referred year: 2014
        #the given Money unit: EUR
        #Qc unit: W
        year=2014
        FC_Dollars=(0.217*math.pow(Qc/1e6,0.65)+3.83)*self._Eur2Dollars[year]  #unit M$
        FC=self._convert_to_RMB_in_target_year(FC_Dollars,year)
        return FC
    
    def _cost_bretonHP_turb(self,TC_Dollars):
        year=2020
        #TC_Dollars=mass_co2*492.2(1-mius)*((Tin+273.15)/(Tout+273.15))*math.log((Tin+273.15)/(Tout+273.15))(1+math.exp(0.036*(Tin+273.15)-65.66))
        TC=self._convert_to_RMB_in_target_year(TC_Dollars,year)
        return TC 
    
    def _cost_bretonHP_comp(self,CC_Dollars):
        year=2020
        #CC_Dollars=mass_co2*59.1(1-mius)*(Pout/Pin)*math.log(Pout/Pin)
        CC=self._convert_to_RMB_in_target_year(CC_Dollars,year)
        return CC 
    
    def _cost_pump(self,W):
        #reference: Chemical process: design and integration
        #the referred year: 2001
        #the given Money unit: USD
        year=2001
        CPump_USD=9840*math.pow(W/1000/4,0.55)/1e6 #unit M$ 
        CPump=self._convert_to_RMB_in_target_year(CPump_USD,year)
        return CPump

    def _cost_HEN(self,CHEN_Dollars):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: EUR
        #Ai unit: m2,p: Pa
        year=2020
        #p=p/1e5
        #CHEN_Dollars=2546.9*math.pow(At,0.67)*math.pow(p,0.28)*self._Eur2Dollars[year]/1e6  #unit M$ 
        CHEN=self._convert_to_RMB_in_target_year(CHEN_Dollars,year)
        return CHEN
    
    def calculate_irr(self,r,a,b,N):
        return a-sum((b/(1+r)**x) for x in range(1, N+1))
    
    def discounted_payback_period(self,C, R, A, N):
        if C <= 0 or R < 0 or A <= 0:

             return 0  # 如果C小于等于0，或者R小于0，或者A小于等于0，返回0

        cumulative_discounted_cash_flow = 0.0
        year = 0
        max_years = N
        while cumulative_discounted_cash_flow < C:
            discounted_A = A / ((1 + R) ** (year + 1))
            cumulative_discounted_cash_flow += discounted_A
            year += 1.0
            if year > max_years:
                     return max_years+1
            if cumulative_discounted_cash_flow >= C and year > 0:
                
                previous_cumulative = cumulative_discounted_cash_flow - discounted_A
                interpolation_factor = (C - previous_cumulative) / (discounted_A - previous_cumulative)
                
                return year  # 直接返回year，因为A是定值，每年的现金流相同
    
    def _convert_to_RMB_in_target_year(self,cost_dollars,given_year):
        #a=b/596.2*798.8*7.2492
        cost_RMB=cost_dollars/self._CEPCIs[given_year]*self._CEPCIs[self._target_year]*self._Dollars2RMB[self._target_year]
        return cost_RMB

if __name__ == '__main__':
    parameters = dict()
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["isentropic_eff_mc"] = 0.88
    parameters["t_isentropic_eff_mc"] = 0.92
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["industrial_waste_heat_t"] =300 #℃
    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_amb"] = 20   #环境温度
    parameters["p_amb"] = 101325   #环境压力

    parameters["p_bray_L"] = 7.5e6
    parameters["Store_electrical_power"] = 50e6  # 机组容量
    parameters["p_water_supply_in"] = 2e5 
    parameters["water_pressure_drop_rate"] = 100 #100Pa/m
    parameters["water_pipe_length"] = 1000
    parameters["water_pump_hydraulic_efficiency"] = 0.75
    parameters["water_pump_mechanical_efficiency"] = 0.94

    parameters["cao_conversion"] = 0.95  #氧化钙转化率
    parameters["cao_purity"] = 0.98 #氢氧化钙含量
    parameters["dehydrator_eff"] = 0.95   #脱水器效率
    parameters["steam_pressure_loss_ratio"] = 0.01
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100

    parameters["hydrator_eff"] = 0.95   #水合器器效率
    parameters["p_bray_L_B"] = 7.5e6

    inputs={}
    inputs["p_bray_H"] =17318583#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12303975 #优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 1e5 #变量4，反应器压力
    inputs["Economic Model Selection"] = 2 #经济模型选择，1Tesio，2Nathan T
    inputs["Compressor power limit"] = 200e6#功率界限，影响齿轮离心和滚筒离心模型的选取，单位W，桶式离心需体积流量
    inputs["Turbine power limit"] = 35e6
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度
    inputs["min_temperature_exchange"] = 15
    inputs["min_temperature_HEN"] = 15
    
    inputs["p_bray_H_B"] = 30e6
    inputs["p_bray_MH_B"] = 16217752.142109105
    inputs["p_bray_ML_B"] = 12217752.142109105
    inputs["p_Hydr"] = 1e5
    inputs["Hydr_overheating_temperature"] = 40

    inputs["Hydr_cao_in"]=440
    inputs["T_X"] = 1


    economic_inputs={}
    #economic_inputs["limestone_price"]=70
    #economic_inputs["calciner_cost_factor"]=1
    economic_inputs["caoh2_unit_price"]=1500 #元/吨
    economic_inputs["caoh2_price/life_rate"] = 30
    economic_inputs["elec_price"]=0.26 #元/千瓦时  ##S6
    economic_inputs["hot_price_h"] = 0.1542 #每千瓦时0.1542元，42.84元/吉焦
    economic_inputs["elec_price_h"] = 1.5912#1.0827
    economic_inputs["Annual_cycle_count"]=320
    economic_inputs["operation_hours"]=economic_inputs["Annual_cycle_count"]*8 # hours 250天2000小时 320天2560小时
    economic_inputs["discount_ratio"]=6/100 #8%
    economic_inputs["operational_years"]=30
    economic_inputs["operation_labour_cost_indictor"]=0.025/2#劳动力比例
    #年运行时间为一般机组的一半
    economic_inputs["maintain_cost_indictor"]=0.025/2#维护比例 0.025,
    #年运行时间为一般机组的一半
    cost = Cost_Estimator(parameters)
    results = cost.solve(inputs,economic_inputs)
    data_for_json = {key: (value.item() if isinstance(value, np.floating) else value) for key, value in results.items()}
    json_output = json.dumps(data_for_json, indent=4)
    #print(json_output)
    print(results )
    #print(results["operation"]["total as-spent"])
    #print(results["IRR"])

