import geatpy as ea
import numpy as np
import concurrent.futures
from pyBrayton import Brayton
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

t_reaction = 465
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

class BProblem(ea.Problem):
    def __init__(self,BraytonHeatPump_parameters):
        name = 'BProblem'  # 初始化name（函数名称，可以随意设置）

         ##  plant variabels
        self._flue_gas_composition=BraytonHeatPump_parameters["flue_gas_composition"]
        #self._isentropic_eff_mc=BraytonHeatPump_parameters["isentropic_eff_mc"]
        #self._t_isentropic_eff_mc=BraytonHeatPump_parameters["t_isentropic_eff_mc"]
        #self._mechanical_eff=BraytonHeatPump_parameters["mechanical_eff"]
        #self._min_temperature_exchange=BraytonHeatPump_parameters["min_temperature_exchange"]
        #self._industrial_waste_heat_t=BraytonHeatPump_parameters["industrial_waste_heat_t"]
        #self._heat_transfer_loss_eff=BraytonHeatPump_parameters["heat_transfer_loss_eff"]
        #self._t_reaction=BraytonHeatPump_parameters["t_reaction"]
        #self._p_bray_M=BraytonHeatPump_parameters["p_bray_M"]
        #self._p_bray_L=BraytonHeatPump_parameters["p_bray_L"]
        self._obj=BraytonHeatPump_parameters["obj"]
        # initialize plant
        BraytonHeatPump_parameters=self._compose_BraytonHeatPump_parameters()
        self._Brayton = Brayton(BraytonHeatPump_parameters)
        # for plant constraints
        self._high_t = 615
        self._hot_util=10
        # for concurrent worker
        self._max_map_executor=20
        M, maxormins, Dim, varTypes, lb, ub, lbin, ubin = self._compose_generic_algorithm_paramters()
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def _compose_generic_algorithm_paramters(self):
        ## generic algorithm parameters 
        M = 1  # 优化目标个数
        maxormins = self._get_maxormins(M)  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 3  # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [300e5,130e5,75e5]  # 高压下界、中压下界
        ub = [300e5,300e5,300e5] # 高压上界、中压下界
        lbin = [1,1,1] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1,1,1]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）

        return M,maxormins,Dim,varTypes,lb,ub,lbin,ubin
    
    def _get_maxormins(self, M):
        # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        if self._obj=="energy":
            maxormins = [-1] * M
            return maxormins
        elif self._obj=="cop":
            maxormins = [-1] * M
            return maxormins
        else:
            raise ValueError("objective is invaid. It should be energy or economic.")
    def _compose_BraytonHeatPump_parameters(self):
        parameters={}
        parameters["flue_gas_composition"] = self._flue_gas_composition
        parameters["isentropic_eff_mc"] = isentropic_eff_mc  #等熵效率
        parameters["t_isentropic_eff_mc"] = t_isentropic_eff_mc
        parameters["mechanical_eff"] = mechanical_eff  #机械效率
        parameters["min_temperature_exchange"] = 15
        parameters["industrial_waste_heat_t"] =200 #℃
        parameters["heat_transfer_loss_eff"] = 0.96
        parameters["t_reaction_B"] = 465
        #parameters["p_bray_H"]=Vars[i,0]
        #parameters["p_bray_M"] = 13e6
        parameters["p_bray_L_B"] = 7.5e6
        parameters["T_amb"] = T_amb
        parameters["p_amb"] = p_amb
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
                inputs["p_bray_H"]=Vars[i,0]
                inputs["p_bray_M1"]=Vars[i,1]
                inputs["p_bray_M2"]=Vars[i,2]
                inputs["Hydrator_heat"]=1498957.696586326-30766.434959109065
                future=map_executor_pool.submit(self._Brayton.solve,inputs)
                futures.append(future)
            # reduce
            for _, future in enumerate(futures):
                # print("I am here")
                results=future.result()
            
                    # objective function
                obj=self._get_obj(results)
                objs.append(obj)
                    #constrainst
                c0=results["B_primary_compressor"]["t_compressor_out"]-self._high_t
                c1=results["B_primary_compressor"]["t_compressor_out"]-self._high_t #
                c2=results["B_secondary_compressor"]["t_compressor_out"]-self._high_t
                    # c3=T_amb+self._delta_T_pinch-result["T_cao_reactor_in"]
                constraints.append([c0,c1,c2])
               
        f=np.reshape(objs,(-1,1)) # 计算目标函数值矩阵
        CV=np.array(constraints) # 构建违反约束程度矩阵
        return f, CV
    def _get_obj(self,results):
        if self._obj=="energy":
            return results["evaluation_indicators"]["Energy efficiency"]
        elif self._obj=="cop":

            return results["evaluation_indicators"]["power"]


        else:
            raise ValueError("objective is invaid. It should be energy or economic.")

    def solve(self,inputs):
        results=self._Brayton.solve(inputs)
        return results