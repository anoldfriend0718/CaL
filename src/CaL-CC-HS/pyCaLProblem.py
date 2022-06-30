import geatpy as ea
import numpy as np
import concurrent.futures
from pyCaLPlant import CalcinerSide, CarbonatorSide

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

calciner_thermal_loss=0.01
deltaTmin_SSHX=20  #https://doi.org/10.1016/j.apenergy.2016.04.053
deltaTmin_SGHX=15  #https://doi.org/10.1016/j.apenergy.2016.04.053
T_cooling_co2=20   #https://doi.org/10.1016/j.apenergy.2016.04.053
n_compression=6    #https://doi.org/10.1016/j.apenergy.2016.04.053
HTC_SSHX=200 ## TODO: too big ,150?
HTC_SGHX=300 #https://doi.org/10.1016/j.ecmx.2020.100039 ## TODO: too big ,150?

## feedwater 
water_pres_drop_rate=100 #Pa/m
hot_water_pipe_length=1000 #m
water_pump_hydraulic_efficiency=0.75
water_pump_mechanical_efficiency=0.94


class CalcProblem(object):
    def __init__(self,parameters) -> None:
        ## plant variables
        self._T_calc=parameters["T_calc"]
        self._cao_conversion=parameters["cao_conversion"]
        self._flue_gas_composistion=parameters["flue_gas_composition"]
        self._decarbonized_rate=parameters["decarbonized_rate"]
        calc_para=self._compose_calc_parameters()
        self._calcs=CalcinerSide(calc_para)

    def _compose_calc_parameters(self):
        parameters={}
        parameters["flue_gas_composition"] = self._flue_gas_composistion
        parameters["isentropic_eff_mc"] = isentropic_eff_mc
        parameters["mechanical_eff"] = mechanical_eff
        parameters["decarbonized_rate"] = self._decarbonized_rate
        parameters["cao_conversion"] = self._cao_conversion
        parameters["calciner_eff"] = 1-calciner_thermal_loss
        parameters["convey_consumption"] = convey_consumption
        parameters["storage_carbonator_distance"] = storage_carbonator_distance
        parameters["T_amb"] = T_amb
        parameters["p_amb"] = p_amb
        parameters["T_calc"]= self._T_calc
        parameters["deltaTmin_SSHX"]=deltaTmin_SSHX
        parameters["deltaTmin_SGHX"]=deltaTmin_SGHX
        parameters["T_cooling_co2"]=T_cooling_co2
        parameters["p_co2_storage"]=p_co2_storage
        parameters["n_compression"]=n_compression
        parameters["cooling_eff"]=cooling_eff
        parameters["HTC_SSHX"]=HTC_SSHX
        parameters["HTC_SGHX"]=HTC_SGHX
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
        self._HTCW=carbonator_parameters["HTCW"]
        self._HRCP=carbonator_parameters["HRCP"]
        self._obj=carbonator_parameters["obj"]

        # initialize plant
        carbonator_parameters=self._compose_carbonator_parameters()
        self._carbonatorSide = CarbonatorSide(carbonator_parameters)
        # for plant constraints
        self._hot_util=1e2
        # for concurrent worker
        self._max_map_executor=10
        M, maxormins, Dim, varTypes, lb, ub, lbin, ubin = self._compose_generic_algorithm_paramters()
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def _compose_generic_algorithm_paramters(self):
        ## generic algorithm parameters 
        M = 1  # 优化目标个数
        maxormins = self._get_maxormins(M)  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        if self._HTCW==1 and self._HRCP==1:
            Dim = 3  # 初始化Dim（决策变量维数）
            varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
            lb = [T_amb+delta_T_pinch, #T_cao,in
                T_amb+delta_T_pinch, #T_flue_gas,in
                T_amb]  # T_cold water  决策变量下界
            ub = [self._T_carb-delta_T_pinch, #T_cao,in
                self._T_carb-delta_T_pinch, #T_flue_gas,in
                self._T_water_reactor_out]  #  T_hot water 决策变量上界
            lbin = [1,0,0] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
            ubin = [1,1,0]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        elif self._HTCW==1 and self._HRCP==0:
            Dim = 2  # 初始化Dim（决策变量维数）
            varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
            lb = [T_amb+delta_T_pinch, #T_cao,in
                T_amb+delta_T_pinch] #T_flue_gas,in
            ub = [self._T_carb-delta_T_pinch, #T_cao,in
                self._T_carb-delta_T_pinch] #T_flue_gas,in
            lbin = [1,0] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
            ubin = [1,1]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        elif self._HTCW==0 and self._HRCP==1:
            Dim = 2  # 初始化Dim（决策变量维数）
            varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
            lb = [T_amb+delta_T_pinch, #T_flue_gas,in
                0]  # m_water,in
            ub = [self._T_carb-delta_T_pinch, #T_flue_gas,in
                10] # m_water,in
            lbin = [0,0] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
            ubin = [1,0]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        else:
            raise ValueError("Invalid HTCW or HRCP!")

        return M,maxormins,Dim,varTypes,lb,ub,lbin,ubin

    def _get_maxormins(self, M):
        # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        if self._obj=="energy":
            maxormins = [-1] * M
            return maxormins
        elif self._obj=="economic":
            maxormins = [1] * M
            return maxormins
        else:
            raise ValueError("objective is invaid. It should be energy or economic.")

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
        parameters["water_pressure_drop_rate"]=water_pres_drop_rate
        parameters["water_pipe_length"]=hot_water_pipe_length
        parameters["water_pump_hydraulic_efficiency"]=water_pump_hydraulic_efficiency
        parameters["water_pump_mechanical_efficiency"]=water_pump_mechanical_efficiency
        parameters["HTCW"]=self._HTCW
        parameters["HRCP"]=self._HRCP
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
                if self._HTCW==1 and self._HRCP==1:
                    inputs["T_flue_gas_reactor_in"]=Vars[i,0]
                    inputs["T_cao_reactor_in"]=Vars[i,1]
                    inputs["T_water_reactor_in"]=Vars[i,2]
                elif self._HTCW==1 and self._HRCP==0:
                    inputs["T_flue_gas_reactor_in"]=Vars[i,0]
                    inputs["T_cao_reactor_in"]=Vars[i,1]
                elif self._HTCW==0 and self._HRCP==1:
                    inputs["T_flue_gas_reactor_in"]=Vars[i,0]
                    inputs["m_water_in"]=Vars[i,1]
                else:
                    raise ValueError("Invalid HTCW or HRCP!")
                future=map_executor_pool.submit(self._carbonatorSide.solve,inputs)
                futures.append(future)
            # reduce
            for _, future in enumerate(futures):
                results=future.result()
                if self._HTCW==1:
                    # objective function
                    obj=self._get_obj(results)
                    objs.append(obj)
                    #constrainst
                    c0=1-results["is_succeed"]
                    c1=results["hot_utility"]-self._hot_util #
                    c2=-results["m_water_in"]
                    constraints.append([c0,c1,c2])
                elif self._HRCP==1:
                    # objective function
                    obj=self._get_obj(results)
                    objs.append(obj)
                    #constrainst
                    c0=1-results["is_succeed"]
                    c1=results["hot_utility"]-self._hot_util #
                    c2=results["T_cao_reactor_in"]-(self._T_carb-delta_T_pinch)
                    # c3=T_amb+delta_T_pinch-result["T_cao_reactor_in"]
                    constraints.append([c0,c1,c2])
                else:
                    raise ValueError("Invalid HTCW or HRCP!")
        f=np.reshape(objs,(-1,1)) # 计算目标函数值矩阵
        CV=np.array(constraints) # 构建违反约束程度矩阵
        return f, CV

    def _get_obj(self,results):
        if self._obj=="energy":
            return results["carb_heat_rec_eff"]
        elif self._obj=="economic":
            if "invCosts" in results.keys() and "specific" in results["invCosts"].keys():
                return results["invCosts"]["specific"]
            inf=1e6
            return inf
        else:
            raise ValueError("objective is invaid. It should be energy or economic.")

    def solve(self,inputs):
        results=self._carbonatorSide.solve(inputs)
        return results
