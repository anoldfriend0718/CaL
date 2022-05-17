from pyCaLPlant import Plant
import geatpy as ea
import numpy as np
import concurrent.futures

class CaLProblem(ea.Problem):  # 继承Problem父类
    def __init__(self,parameters):
        name = 'CaLProblem'  # 初始化name（函数名称，可以随意设置）
        ## plant variables 
        self._T_reaction=parameters["T_reaction"]
        self._m_cao_i=parameters["m_cao_i"]
        self._p_co2_storage=parameters["p_co2_storage"]
        self._cao_conversion=parameters["cao_conversion"]

        ## fixed plant parameters
        self._isentropic_eff_mt = 0.89
        self._isentropic_eff_mc = 0.87
        self._isentropic_eff_st = 0.75
        self._mechanical_eff = 0.97
        self._carbonator_eff = 0.99
        self._convey_consumption = 10e3/100
        self._storage_carbonator_distance = 100
        self._cooling_eff = 0.8e-2
        self._pres_loss_rec = 0.04
        self._pres_loss_stoic = 0.01
        self._pres_loss_mix = 0.06
        self._delta_T_pinch=15
        self._T_amb = 20

        # initialize plant
        parameters=self._compose_plant_parameters()
        self._plant = Plant(parameters)

        # for plant constraints
        self._p_amb=101325
        self._hot_util=1e3

        # for concurrent worker
        self._max_map_executor=10

        ## generic algorithm parameters 
        M = 1  # 优化目标个数
        maxormins = [-1] * M  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 6  # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）

        lb = [1.5e5, #pressure carbonator
              1.2, # expantion ratio
              self._T_amb+self._delta_T_pinch, #T_co2,in
              310, #T_cao,in
              self._T_amb+self._delta_T_pinch, # T_compressor,in
              250]  # T_main_turbine, in 决策变量下界
        ub = [15e5,#5e5 #pressure carbonator
              15, # expantion ratio
              self._T_reaction-self._delta_T_pinch, #T_co2,in
              self._T_reaction-self._delta_T_pinch, #T_cao,in
              250, # T_compressor,in
              650]  #  T_main_turbine, in 决策变量上界
        lbin = [1] * Dim  # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def _compose_plant_parameters(self):
        parameters={}
        parameters["isentropic_eff_mt"]=self._isentropic_eff_mt
        parameters["isentropic_eff_mc"]=self._isentropic_eff_mc
        parameters["isentropic_eff_st"]= self._isentropic_eff_st 
        parameters["mechanical_eff"]=self._mechanical_eff
        parameters["carbonator_eff"]=self._carbonator_eff 
        parameters["convey_consumption"]=self._convey_consumption 
        parameters["storage_carbonator_distance"]=self._storage_carbonator_distance 
        parameters["cooling_eff "]=self._cooling_eff 
        parameters["pres_loss_rec"]=self._pres_loss_rec 
        parameters["pres_loss_stoic"]=self._pres_loss_stoic 
        parameters["pres_loss_mix"]=self._pres_loss_mix 
        parameters["delta_T_pinch"]=self._delta_T_pinch 
        parameters["T_amb"]=self._T_amb 
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
                input["m_cao_i"] = self._m_cao_i
                input["p_co2_storage"] = self._p_co2_storage
                input["T_reaction"] = self._T_reaction
                input["cao_conversion"]=self._cao_conversion
                input["p_reaction"]=Vars[i,0]
                input["expansion_ratio"]=Vars[i,1]
                input["T_co2_reactor_in"]=Vars[i,2]
                input["T_cao_reactor_in"]=Vars[i,3]
                input["T_co2_compressor_in"]=Vars[i,4]
                input["T_co2_storage_turbine_in"]=Vars[i,5]
                future=map_executor_pool.submit(self._plant.solve,input)
                futures.append(future)
            # reduce
            for _, future in enumerate(futures):
                result=future.result()
                # objective function
                obj=result["plant_eff"]
                objs.append(obj)
                #constrainst
                c1=result["hot_utility"]-self._hot_util #
                c2=self._p_amb-result["p_co2_compressor_in"]
                constraints.append([c1,c2])
        f=np.reshape(objs,(-1,1)) # 计算目标函数值矩阵
        CV=np.array(constraints) # 构建违反约束程度矩阵
        return f, CV
