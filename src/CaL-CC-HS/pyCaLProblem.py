from pyCaLPlant import CalcinerSide, CarbonatorSide

import geatpy as ea
import numpy as np
import concurrent.futures

# global parameters
T_amb = 20
p_amb = 101325
convey_consumption = 10e3/100  # reference: https://doi.org/10.1016/j.ecmx.2019.100025
storage_carbonator_distance = 100  # reference: https://doi.org/10.1016/j.ecmx.2019.100025
cooling_eff = 0.8e-2  # reference: https://doi.org/10.1016/j.ecmx.2019.100025
isentropic_eff_mc = 0.87 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
isentropic_eff_fan= 0.80 # reference: https://doi.org/10.1016/j.ecmx.2020.100038
mechanical_eff = 0.97 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
p_co2_storage=75e5 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
delta_T_pinch=15 # reference: https://doi.org/10.1016/j.ecmx.2019.100025
carbonator_thermal_loss=0.01  # reference: https://doi.org/10.1016/j.ecmx.2019.100025
flue_gas_pressure_loss_ratio = 0.01 #https://doi.org/10.1016/j.ecmx.2019.100025
decarbon_flue_gas_pressure_loss_ratio  = 0.02 #estimated based on flue_gas_pressure_loss_ratio
carbonator_pressure_loss=15e3 #https://doi.org/10.1016/j.jclepro.2019.02.049


class CalcProblem(object):
    def __init__(self,parameters) -> None:
        ## plant variables
        self._T_calc=parameters["T_calc"]
        self._cao_conversion=parameters["cao_conversion"]
        self._flue_gas_composistion=parameters["flue_gas_composition"]
        self._decarbonized_rate=parameters["decarbonized_rate"]
        ## fixed parameters
        self._calciner_eff=0.99
        self._deltaTmin_SSHX=20
        self._deltaTmin_SGHX=15
        self._T_cooling_co2=20
        self._n_compression=6

        ## global parameters
        self._isentropic_eff_mc=isentropic_eff_mc
        self._mechanical_eff=mechanical_eff
        self._convey_consumption=convey_consumption
        self._storage_carbonator_distance=storage_carbonator_distance
        self._T_amb=T_amb
        self._p_amb=p_amb
        self._cooling_eff=cooling_eff
        self._p_co2_storage=p_co2_storage

        calc_para=self._compose_calc_parameters()
        self._calcs=CalcinerSide(calc_para)


    def _compose_calc_parameters(self):
        parameters={}
        parameters["flue_gas_composition"] = self._flue_gas_composistion
        parameters["isentropic_eff_mc"] = self._isentropic_eff_mc
        parameters["mechanical_eff"] = self._mechanical_eff
        parameters["decarbonized_rate"] = self._decarbonized_rate
        parameters["cao_conversion"] = self._cao_conversion
        parameters["calciner_eff"] = self._calciner_eff
        parameters["convey_consumption"] = self._convey_consumption
        parameters["storage_carbonator_distance"] = self._storage_carbonator_distance
        parameters["T_amb"] = self._T_amb
        parameters["p_amb"] = self._p_amb
        parameters["T_calc"]= self._T_calc
        parameters["deltaTmin_SSHX"]=self._deltaTmin_SSHX
        parameters["deltaTmin_SGHX"]=self._deltaTmin_SGHX
        parameters["T_cooling_co2"]=self._T_cooling_co2
        parameters["p_co2_storage"]=self._p_co2_storage
        parameters["n_compression"]=self._n_compression
        parameters["cooling_eff"]=self._cooling_eff
        return parameters

    def solve(self,inputs):
        results=self._calcs.solve(inputs)
        return results

    def opt(self, m,mfrac):
        results=self._calcs.calciner(mfrac,m)
        return results["We_calc"]




class CarbProblem(ea.Problem):  # 继承Problem父类
    def __init__(self,carbonator_parameters):
        name = 'CaLProblem'  # 初始化name（函数名称，可以随意设置）
        ##  plant variabels
        self._T_carb=carbonator_parameters["T_carb"]
        self._cao_conversion=carbonator_parameters["cao_conversion"]
        self._flue_gas_composistion=carbonator_parameters["flue_gas_composition"]
        self._vol_rate_flue_gas=carbonator_parameters["vol_rate_flue_gas"]
        self._decarbonized_rate=carbonator_parameters["decarbonized_rate"]
        self._T_flue_gas=carbonator_parameters["T_flue_gas"]
        self._T_water_reactor_out=carbonator_parameters["T_water_reactor_out"]

        # initialize plant
        carbonator_parameters=self._compose_carbonator_parameters()
        self._carbonatorSide = CarbonatorSide(carbonator_parameters)

        # for plant constraints
        self._hot_util=1e2

        # for concurrent worker
        self._max_map_executor=10

        ## generic algorithm parameters 
        M = 1  # 优化目标个数
        maxormins = [-1] * M  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 3  # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）

        lb = [T_amb+delta_T_pinch, #T_cao,in
              T_amb+delta_T_pinch, #T_flue_gas,in
              T_amb]  # T_cold water  决策变量下界
        ub = [self._T_carb-delta_T_pinch, #T_cao,in
              self._T_carb-delta_T_pinch, #T_flue_gas,in
              self._T_water_reactor_out]  #  T_hot water 决策变量上界
        lbin = [1,1,0] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1,1,0]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def _compose_carbonator_parameters(self):
        parameters={}

        parameters["flue_gas_composition"]=self._flue_gas_composistion
        parameters["isentropic_eff_mc"] = isentropic_eff_fan
        parameters["mechanical_eff"] = mechanical_eff
        parameters["flue_gas_pressure_loss_ratio"] = flue_gas_pressure_loss_ratio
        parameters["decarbon_flue_gas_pressure_loss_ratio"] = decarbon_flue_gas_pressure_loss_ratio
        parameters["carbonator_pressure_loss"]=carbonator_pressure_loss
        parameters["decarbonized_rate"]=self._decarbonized_rate
        parameters["vol_rate_flue_gas"]=self._vol_rate_flue_gas
        parameters["T_flue_gas"]=self._T_flue_gas
        parameters["T_water_reactor_out"]=self._T_water_reactor_out
        parameters["T_carb"]=self._T_carb
        parameters["cao_conversion"]=self._cao_conversion

        parameters["carbonator_eff"]=1-carbonator_thermal_loss
        parameters["convey_consumption"]=convey_consumption 
        parameters["storage_carbonator_distance"]=storage_carbonator_distance 
        parameters["cooling_eff "]=cooling_eff
        parameters["delta_T_pinch"]=delta_T_pinch 
        parameters["T_amb"]=T_amb 
        parameters["p_amb"] = p_amb
        parameters["p_water"] = parameters["p_amb"]

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
                input={}
                input["T_flue_gas_reactor_in"]=Vars[i,0]
                input["T_cao_reactor_in"]=Vars[i,1]
                input["T_water_reactor_in"]=Vars[i,2]
                future=map_executor_pool.submit(self._carbonatorSide.solve,input)
                futures.append(future)
            # reduce
            for _, future in enumerate(futures):
                result=future.result()
                # objective function
                obj=result["carb_heat_rec_eff"]
                objs.append(obj)
                #constrainst
                c1=result["hot_utility"]-self._hot_util #
                c2=-result["m_water_in"]
                constraints.append([c1,c2])
        f=np.reshape(objs,(-1,1)) # 计算目标函数值矩阵
        CV=np.array(constraints) # 构建违反约束程度矩阵
        return f, CV

    def solve(self,inputs):
        results=self._carbonatorSide.solve(inputs)
        return results
