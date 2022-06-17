from pyCaLPlant import CalcinerSide, CarbonatorSide
import geatpy as ea
import numpy as np
import concurrent.futures

# global parameters

convey_consumption = 10e3/100
storage_carbonator_distance = 100
cooling_eff = 0.8e-2
T_amb = 20
p_amb = 101325
isentropic_eff_mc = 0.87
mechanical_eff = 0.97
p_co2_storage=75e5

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
    def __init__(self,parameters):
        name = 'CaLProblem'  # 初始化name（函数名称，可以随意设置）
        ##  plant variabels
        self._T_carb=parameters["T_carb"]
        self._cao_conversion=parameters["cao_conversion"]
        self._flue_gas_composistion=parameters["flue_gas_composition"]
        self._vol_rate_flue_gas=parameters["vol_rate_flue_gas"]
        self._decarbonized_rate=parameters["decarbonized_rate"]
        self._T_flue_gas=parameters["T_flue_gas"]
        self._T_water_reactor_out=parameters["T_water_reactor_out"]

        ## fixed plant parameters
        self._carbonator_eff = 0.99
        self._delta_T_pinch=15
        self._flue_gas_pressure_loss_1 = 0.04
        self._flue_gas_pressure_loss_2 = 0.04
       
       #global paramters
        self._convey_consumption = convey_consumption
        self._storage_carbonator_distance = storage_carbonator_distance
        self._cooling_eff = cooling_eff
        self._T_amb = T_amb
        self._p_amb = p_amb
        self._isentropic_eff_mc = isentropic_eff_mc
        self._mechanical_eff = mechanical_eff

        # initialize plant
        plant_parameters=self._compose_plant_parameters()
        self._plant = CarbonatorSide(plant_parameters)

        # for plant constraints
        self._hot_util=1e2

        # for concurrent worker
        self._max_map_executor=10

        ## generic algorithm parameters 
        M = 1  # 优化目标个数
        maxormins = [-1] * M  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 3  # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）

        lb = [self._T_amb+self._delta_T_pinch, #T_cao,in
              self._T_amb+self._delta_T_pinch, #T_flue_gas,in
              self._T_amb]  # T_cold water  决策变量下界
        ub = [self._T_carb-self._delta_T_pinch, #T_cao,in
              self._T_carb-self._delta_T_pinch, #T_flue_gas,in
              self._T_water_reactor_out]  #  T_hot water 决策变量上界
        lbin = [1,1,0] # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1,1,0]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def _compose_plant_parameters(self):
        parameters={}

        parameters["flue_gas_composition"]=self._flue_gas_composistion
        parameters["isentropic_eff_mc"] = self._isentropic_eff_mc
        parameters["mechanical_eff"] = self._mechanical_eff
        parameters["flue_gas_pressure_loss_1"] = self._flue_gas_pressure_loss_1
        parameters["flue_gas_pressure_loss_2"] = self._flue_gas_pressure_loss_2
        parameters["decarbonized_rate"]=self._decarbonized_rate
        parameters["vol_rate_flue_gas"]=self._vol_rate_flue_gas
        parameters["T_flue_gas"]=self._T_flue_gas
        parameters["T_water_reactor_out"]=self._T_water_reactor_out
        parameters["T_carb"]=self._T_carb
        parameters["cao_conversion"]=self._cao_conversion

        parameters["carbonator_eff"]=self._carbonator_eff 
        parameters["convey_consumption"]=self._convey_consumption 
        parameters["storage_carbonator_distance"]=self._storage_carbonator_distance 
        parameters["cooling_eff "]=self._cooling_eff 
        parameters["delta_T_pinch"]=self._delta_T_pinch 
        parameters["T_amb"]=self._T_amb 
        parameters["p_amb"] = self._p_amb
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
                future=map_executor_pool.submit(self._plant.solve,input)
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
        results=self._plant.solve(inputs)
        return results
