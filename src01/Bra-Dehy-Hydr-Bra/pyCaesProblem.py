import geatpy as ea
import numpy as np
import concurrent.futures
from pyHENPinch import Hen_pinch_analyzer
from Cp0massWrapper import Cp0mass_Wrapper
from pyDehydrator import Dehydrator
from pyBrayton import Brayton
from pyPinch import PyPinch
from pyHydrator import Hydrator
import pandas as pd
# global parameters
T_amb = 20
p_amb = 101325
P_amb = 101325
convey_consumption = 10e3/100  # reference: https://doi.org/10.1016/j.ecmx.2019.100025

flue_gas_composistion = dict()
flue_gas_composistion["co2"] = 0.1338
flue_gas_composistion["o2"] = 0.0384
flue_gas_composistion["n2"] = 0.6975
# storage_dehydrator_distance = 100  # reference: https://doi.org/10.1016/j.ecmx.2019.100025
dehydrator_eff = 0.95
Store_electrical_power = 1e6
cooling_eff = 0.8e-2  # reference: https://doi.org/10.1016/j.ecmx.2019.100025
isentropic_eff_mc = 0.88 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
t_isentropic_eff_mc = 0.92
industrial_waste_heat_t = 300
isentropic_eff_fan= 0.80 # reference: https://doi.org/10.1016/j.ecmx.2020.100038
mechanical_eff = 0.98 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
heat_transfer_loss_eff = 0.96
p_co2_storage=75e5 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
# delta_T_pinch=15 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
# dehydrator_thermal_loss=0.01  # reference: https://doi.org/10.1016/j.ecmx.2019.100025
flue_gas_pressure_loss_ratio = 0.01 #https://doi.org/10.1016/j.ecmx.2019.100025
steam_pressure_loss_ratio = 0.01
decarbon_flue_gas_pressure_loss_ratio  = 0.02 #estimated based on flue_gas_pressure_loss_ratio
dehydrator_pressure_loss=15e3 #https://doi.org/10.1016/j.jclepro.2019.02.049

t_reaction = 525
p_bray_H = 20e6
p_bray_M = 13e6
p_bray_L = 7.5e6

deltaTmin_SSHX=25 #20 #https://doi.org/10.1016/j.apenergy.2016.04.053
deltaTmin_SGHX=20 #15 #https://doi.org/10.1016/j.apenergy.2016.04.053
T_cooling_co2=20   #https://doi.org/10.1016/j.apenergy.2016.04.053
n_compression=6    #https://doi.org/10.1016/j.apenergy.2016.04.053
HTC_SSHX=200 ## TODO: too big ,150?
HTC_SGHX=300 #https://doi.org/10.1016/j.ecmx.2020.100039 ## TODO: too big ,150?

## feedwater 
water_pres_drop_rate=100 #Pa/m
hot_water_pipe_length=1000 #m
water_pump_hydraulic_efficiency=0.75
water_pump_mechanical_efficiency=0.94
class CaesProblem(ea.Problem):
    def __init__(self,Caes_parameters) -> None:
        name = 'CaesProblem'  # 初始化name（函数名称，可以随意设置）
        ## plant variables
        self._flue_gas_composistion=Caes_parameters["flue_gas_composition"]
        self._obj=Caes_parameters["obj"]

        Caes_parameters=self._compose_caes_parameters()
        self._Hydrator = Hydrator(Caes_parameters)
        # for plant constraints
        self._hot_util=10
        # for concurrent worker
        self._max_map_executor=12
        M, maxormins, Dim, varTypes, lb, ub, lbin, ubin = self._compose_generic_algorithm_paramters()
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def _compose_generic_algorithm_paramters(self):
        ## generic algorithm parameters 
        M = 1 # 优化目标个数
        maxormins = self._get_maxormins(M)  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 13  # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [16e6,9e6,300,1e5,20,20e6,13e6,12217752.142109105,1e5,40,350,8,0]  
        # 12热泵高压，中间压，3氢氧化钙入反应器温度,45储压力过热
        #678布雷顿高压，高低温中间压，910释压力过热,11氧化钙入反应器温度
        #12储能侧供暖水流量、13释能侧补烟气
        ub = [20e6,15e6,500,1e5,20,30e6,20e6,12217752.142109105,1e5,40,440,15,3] #  上界
        lbin = [1,1,1,1,1,1,1,1,1,1,1,1,1] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含） 
        ubin = [1,1,1,1,1,1,1,1,1,1,1,1,1]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）

        return M,maxormins,Dim,varTypes,lb,ub,lbin,ubin
    
    def _get_maxormins(self, M):
        # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        if self._obj=="energy":
            maxormins = [-1] * M
            return maxormins
        elif self._obj=="power":
            maxormins = [-1] * M
            return maxormins
        elif self._obj=="exergy":
            maxormins = [-1] * M
            return maxormins
        elif self._obj=="cop":
            maxormins = [-1] * M
            return maxormins
        else:
            raise ValueError("objective is invaid. It should be energy or economic.")

    def _compose_caes_parameters(self):
        parameters = dict()
        flue_gas_composistion = dict()
        flue_gas_composistion["co2"] = 0.1338
        flue_gas_composistion["o2"] = 0.0384
        flue_gas_composistion["n2"] = 0.6975
        parameters["flue_gas_composition"] = flue_gas_composistion
        parameters["isentropic_eff_mc"] = 0.88
        parameters["t_isentropic_eff_mc"] = 0.92
        parameters["mechanical_eff"] = 0.98   #机械效率
        parameters["min_temperature_exchange"] = 15 
        parameters["deltaTmin_SSHX"] = parameters["min_temperature_exchange"]+5   #固-固换热器最小温差
        parameters["deltaTmin_SGHX"] = parameters["min_temperature_exchange"]   #固-气换热器最小温差
        parameters["industrial_waste_heat_t"] =300 #℃
        parameters["heat_transfer_loss_eff"] = 0.96
        parameters["t_amb"] = 20   #环境温度
        parameters["p_amb"] = 101325   #环境压力

        parameters["p_bray_L"] = 7.5e6
        parameters["Store_electrical_power"] = 1e6

        parameters["cao_conversion"] = 0.95  #氧化钙转化率
        parameters["cao_purity"] = 0.98 #氢氧化钙含量
        parameters["dehydrator_eff"] = 0.95   #脱水器效率
        parameters["steam_pressure_loss_ratio"] = 0.01
        parameters["convey_consumption"] = 10e3/100
        parameters["storage_dehydrator_distance"] = 100

        parameters["hydrator_eff"] = 0.95   #水合器器效率
        parameters["p_bray_L_B"] = 7.5e6

        return parameters
    
    def evalVars(self, Vars):  # 目标函数
        # with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_map_executor) as map_executor_pool:
        with concurrent.futures.ProcessPoolExecutor(max_workers=self._max_map_executor) as map_executor_pool:
            futures=[]
            objs=[]
            constraints=[]
            case_num=Vars.shape[0]
            # map
            for i in np.arange(0,case_num,1):
                inputs={}
                inputs["p_bray_H"] = Vars[i,0]#优化变量1，热泵循环最高压力
                inputs["p_bray_M"] = Vars[i,1] #优化变量2，热泵循环中间压力
                inputs["Dehy_caoh2_in"]=Vars[i,2]
                inputs["p_Dehy"] = Vars[i,3] #变量4，反应器压力
                inputs["Dehy_overheating_temperature"] = Vars[i,4] #变量2，脱水反应器过热温度

                inputs["p_bray_H_B"] = Vars[i,5]
                inputs["p_bray_MH_B"] = Vars[i,6]
                inputs["p_bray_ML_B"] = Vars[i,7]
                inputs["p_Hydr"] = Vars[i,8]
                inputs["Hydr_overheating_temperature"] = Vars[i,9]

                inputs["Hydr_cao_in"]=Vars[i,10]

                inputs["cn_m1"] = Vars[i,11]
                inputs["m2"] = Vars[i,12]#补充烟气
                inputs["T_X"] = 1
                inputs["m1"] = 0#供暖水
                future=map_executor_pool.submit(self._Hydrator.solve,inputs)
                futures.append(future)
            # reduce
            for _, future in enumerate(futures):
                # print("I am here")
                results=future.result()
            
                    # objective function
                obj=self._get_obj(results)
                objs.append(obj)
                    #constrainst
                c0 = results["BH"]["primary_compressor"]["t_compressor_out"]-610
                c1 = results["BH"]["secondary_compressor"]["t_compressor_out"]-610
                c2 = results["Dehy"]["hot_utility"]-self._hot_util
                c3 = results["Hydr"]["HEN"]["hot_utility"]-self._hot_util
                c4 = 20-(results["Hydr"]["initialvalue"]["T_hydr"]-
                         results["Hydr"]["HEN"]["p"]["TTARGET"]["C_Steam"])
                c5 = results["Dehy"]["cold_utility"]-53000
                c6 = results["Hydr"]["HEN"]["cold_utility"] - 74000
                
                    # c3=T_amb+self._delta_T_pinch-result["T_cao_reactor_in"]
                constraints.append([c0,c1,c2,c3,c4,c5,c6])
               
        f=np.reshape(objs,(-1,1)) # 计算目标函数值矩阵
        CV=np.array(constraints) # 构建违反约束程度矩阵
        return f, CV
    def _get_obj(self,results):
        if self._obj=="energy":
            return results["Case_All"]["energy_eff"]
        elif self._obj=="power":
            a=results["Case_All"]["Round-trip_eff"]
            return a
        elif self._obj=="exergy":
            return results["Case_All"]["exergy_eff"]
        elif self._obj=="cop":
            if "invCosts" in results.keys() and "specific" in results["invCosts"].keys():
                return results["cop"]
            inf=1e6
            return inf
        else:
            raise ValueError("objective is invaid. It should be energy or economic.")

    def solve(self,inputs):
        results=self._Hydrator.solve(inputs)
        return results
    
