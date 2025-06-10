import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")
import math
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from Cp0massWrapper import Cp0mass_Wrapper

M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol
M_H20 = 18e-3  # kg/mol

class BraytonHeatPump(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
        self._flue_gas_composition = self._pw._norm_flue_gas_composition
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]#等熵效率
        self._t_isentropic_eff_mc = parameters["t_isentropic_eff_mc"]#透平等熵效率
        self._mechanical_eff = parameters["mechanical_eff"]#机械效率
        self._min_temperature_exchange=parameters["min_temperature_exchange"]#最小换热温差
        self._industrial_waste_heat_t=parameters["industrial_waste_heat_t"]#工业余热温度
        self._heat_transfer_loss_eff=parameters["heat_transfer_loss_eff"]#换热损失
        self._p_bray_L = parameters["p_bray_L"]
        self._p_amb = parameters["p_amb"]
        self._T_amb = parameters["t_amb"]
        self._Store_electrical_power=parameters["Store_electrical_power"]
        
    def solve(self,inputs,a):
        results = {}
        self._p_bray_H = inputs["p_bray_H"] 
        self._p_bray_M =inputs["p_bray_M"] 
        self._p_reaction =inputs["p_Dehy"]
        self._Dehy_ot = inputs["Dehy_overheating_temperature"] 
        #self._p_bray_H = inputs["p_bray_H"]
        #Basic input data
        initialvalue = self.initialvalue()
        results["initialvalue"] = initialvalue
        #t_equilibrium = self.equilibrium()
        self._t_reaction = a+self._Dehy_ot
        results["initialvalue"]["t_reaction"]= self._t_reaction
        #The primary compressor inlet is the starting point
        primary_compressor = self.compressor(self._t_reaction,
                                             self._p_bray_L,
                                             self._p_bray_M)
        results["primary_compressor"] = primary_compressor
        #Primary heat exchanger
        t_h_e1_in=results["primary_compressor"]["t_compressor_out"]
        primary_h_exchanger=self.h_exchanger(t_h_e1_in,
                                                     self._t_reaction+self._min_temperature_exchange,
                                                     self._p_bray_M)
        results["primary_h_exchanger"] = primary_h_exchanger
        #Secondary compressors
        secondary_compressor = self.compressor(self._t_reaction+self._min_temperature_exchange,
                                             self._p_bray_M,
                                             self._p_bray_H)
        results["secondary_compressor"] = secondary_compressor
        #Secondary heat exchanger
        t_h_e2_in=results["secondary_compressor"]["t_compressor_out"]
        secondary_h_exchanger=self.h_exchanger(t_h_e2_in,
                                                       self._t_reaction+self._min_temperature_exchange,
                                                       self._p_bray_H)
        results["secondary_h_exchanger"] = secondary_h_exchanger
        #heat_recovery (main heat exchanger)
        t_h_em_hh_in=results["secondary_h_exchanger"]["t_h_exchanger_out"]
        mian_h_exchanger=self.h_exchanger_main(t_h_em_hh_in)
        results["mian_h_exchanger"] = mian_h_exchanger
        #Turbine
        t_p_tur_in=results["mian_h_exchanger"]["t_h_exchangerm_hh_out"]
        turbine = self.turbine(t_p_tur_in,
                                                self._p_bray_H,
                                                self._p_bray_L)
        results["turbine"] = turbine
        #Industrial waste heat heating part
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        t_p_rec_in=results["turbine"]["t_turbine_out"]
        heat_recovery = self.heat_recovery(t_p_rec_in,
                                               self._p_bray_L,
                                               flue_gas_name)
        results["heat_recovery"] = heat_recovery
        #Data synthesis
        power_in=self._Store_electrical_power
        evaluation_indicators = self.evaluation_indicators(results,power_in)
        results["evaluation_indicators"] = evaluation_indicators
        results["energy_eff"] = results["evaluation_indicators"]["hot_output"]/(results["evaluation_indicators"]["power_cost"]+results["evaluation_indicators"]["hot_cost"])
        results["exergy_eff"] = results["evaluation_indicators"]["exergy_out"]/results["evaluation_indicators"]["exergy_in"]
        results["cop"] = results["evaluation_indicators"]["cop"]

        return results
    
    def initialvalue(self):
        results = {}
        results["isentropic_eff_mc"] = self._isentropic_eff_mc
        results["mechanical_eff"] = self._mechanical_eff
        results["min_temperature_exchange"]= self._min_temperature_exchange 
        results["industrial_waste_heat_t"] = self._industrial_waste_heat_t
        results["p_reaction"] = self._p_reaction
        results["p_bray_H"] = self._p_bray_H
        results["p_bray_M"] = self._p_bray_M
        results["p_bray_L"] = self._p_bray_L
        return results
    
    def equilibrium(self):
        p = self._p_reaction
        t = (-12845/((math.log(p/1e5))-16.508))-273.15
        return t

    def compressor(self, T_in, P_in, P_out):
        t_compressor_in = T_in
        p_compressor_in = P_in
        p_compressor_out = P_out
        h_compressor_in = CP.PropsSI('H', 'T', t_compressor_in+273.15, 'P', p_compressor_in, "REFPROP::co2")
        s_compressor_in = CP.PropsSI('S', 'T', t_compressor_in+273.15, 'P', p_compressor_in, "REFPROP::co2")
        s_hypothesis = s_compressor_in
        h_hypothesis = CP.PropsSI('H', 'S', s_hypothesis, 'P', p_compressor_out, "REFPROP::co2")
        h_compressor_out = h_compressor_in+(h_hypothesis-h_compressor_in)/self._isentropic_eff_mc
        power_compressor = ((h_hypothesis-h_compressor_in)/self._isentropic_eff_mc)/self._mechanical_eff
        e_lost_compressor = power_compressor*(1-self._mechanical_eff)
        t_compressor_out = CP.PropsSI('T', 'H', h_compressor_out, 'P', p_compressor_out, "REFPROP::co2")-273.15
        s_compressor_out = CP.PropsSI('S', 'H', h_compressor_out, 'P', p_compressor_out, "REFPROP::co2")
        results = {}
        results["t_compressor_in"] = t_compressor_in 
        results["p_compressor_in"] = p_compressor_in
        results["h_compressor_in"] = h_compressor_in
        results["s_compressor_in"] = s_compressor_in
        results["t_compressor_out"] = t_compressor_out 
        results["p_compressor_out"] = p_compressor_out
        results["h_compressor_out"] = h_compressor_out
        results["s_compressor_out"] = s_compressor_out
        results["power_compressor"] = power_compressor
        results["e_lost_compressor"] = e_lost_compressor
        results["exergy_lost_compressor"] = (self._T_amb+273.15)*(s_compressor_out-s_compressor_in)
        results["ex_co2"]={}
        results["ex_co2"]["in"]=self.E(t_compressor_in,p_compressor_in)
        results["ex_co2"]["out"]=self.E(t_compressor_out,p_compressor_out)
        return results
    
    def E(self,T,P):
        H1=CP.PropsSI('H', 'T', T+273.15, 'P', P, "REFPROP::co2")
        S1=CP.PropsSI('S', 'T', T+273.15, 'P', P, "REFPROP::co2")
        H0=CP.PropsSI('H', 'T', self._T_amb+273.15, 'P', self._p_amb, "REFPROP::co2")
        S0=CP.PropsSI('S', 'T', self._T_amb+273.15, 'P', self._p_amb, "REFPROP::co2")
        a=H1-H0-(self._T_amb+273.15)*(S1-S0)
        return a

    def h_exchanger(self,T_in,T_out,P):
        t_h_exchanger_in=T_in
        t_h_exchanger_out=T_out
        h_h_exchanger_in =CP.PropsSI('H', 'T', t_h_exchanger_in+273.15,  'P', P, "REFPROP::co2")
        h_h_exchanger_out=CP.PropsSI('H', 'T', t_h_exchanger_out+273.15, 'P', P, "REFPROP::co2")
        hot_out_h_exchanger=h_h_exchanger_in-h_h_exchanger_out
        results = {}
        results["t_h_exchanger_out"] = t_h_exchanger_out
        results["p_h_exchanger_out"] = P
        results["h_h_exchanger_out"] = h_h_exchanger_out
        results["s_h_exchanger_out"] = CP.PropsSI('S', 'T', t_h_exchanger_out+273.15, 'P', P, "REFPROP::co2")
        results["hot_out_h_exchanger"] = hot_out_h_exchanger
        a1=T_in+273.15
        a2=T_out+273.15
        results["Process_taste"]=(a1-a2-293.15*math.log(a1/a2))/(a1-a2)
        return results
        
    def h_exchanger_main(self,T_in):
        t_h_exchangerm_hh_in=T_in#高温
        t_h_exchangerm_ll_in=self._industrial_waste_heat_t-self._min_temperature_exchange
        t_h_exchangerm_ll_out=t_h_exchangerm_hh_in-self._min_temperature_exchange
        h_h_exchangerm_ll_in =CP.PropsSI('H', 'T', t_h_exchangerm_ll_in+273.15,  'P', self._p_bray_L, "REFPROP::co2")
        h_h_exchangerm_ll_out=CP.PropsSI('H', 'T', t_h_exchangerm_ll_out+273.15, 'P', self._p_bray_L, "REFPROP::co2")
        s_h_exchangerm_ll_out=CP.PropsSI('S', 'T', t_h_exchangerm_ll_out+273.15, 'P', self._p_bray_L, "REFPROP::co2")
        s_h_exchangerm_ll_in =CP.PropsSI('S', 'T', t_h_exchangerm_ll_in+273.15,  'P', self._p_bray_L, "REFPROP::co2")
        hot_ll_exchange=h_h_exchangerm_ll_out-h_h_exchangerm_ll_in
        b=hot_ll_exchange
        hot_hh_exchange= hot_ll_exchange+hot_ll_exchange*(1-self._heat_transfer_loss_eff)
        a=hot_hh_exchange
        h_h_exchangerm_hh_in =CP.PropsSI('H', 'T', t_h_exchangerm_hh_in+273.15,  'P', self._p_bray_H, "REFPROP::co2")
        s_h_exchangerm_hh_in =CP.PropsSI('S', 'T', t_h_exchangerm_hh_in+273.15,  'P', self._p_bray_H, "REFPROP::co2")
        h_h_exchangerm_hh_out=h_h_exchangerm_hh_in-hot_hh_exchange
        t_h_exchangerm_hh_out=CP.PropsSI('T', 'H', h_h_exchangerm_hh_out,  'P', self._p_bray_H, "REFPROP::co2")-273.15
        s_h_exchangerm_hh_out=CP.PropsSI('S', 'T', t_h_exchangerm_hh_out+273.15, 'P', self._p_bray_H, "REFPROP::co2")
        results = {}
        results["t_h_exchangerm_hh_out"] = t_h_exchangerm_hh_out
        results["t_h_exchangerm_hh_in"] = t_h_exchangerm_hh_in 
        results["p_h_exchangerm_hh_out"] = self._p_bray_H
        results["h_h_exchangerm_hh_out"] = h_h_exchangerm_hh_out
        results["s_h_exchangerm_hh_out"] = s_h_exchangerm_hh_out
        results["t_h_exchangerm_ll_in"] = t_h_exchangerm_ll_in
        results["t_h_exchangerm_ll_out"] = t_h_exchangerm_ll_out 
        results["p_h_exchangerm_ll_in"] = self._p_bray_L
        results["h_h_exchangerm_ll_in"] = h_h_exchangerm_ll_in
        results["s_h_exchangerm_ll_in"] = s_h_exchangerm_ll_in
        results["h_lost_exchangerm"] = hot_hh_exchange-hot_ll_exchange
        results["hot_exergy"]=h_h_exchangerm_hh_in-h_h_exchangerm_hh_out+(self._T_amb+273.15)*(
            s_h_exchangerm_hh_out-s_h_exchangerm_hh_in)
        results["cold_exergy"]=h_h_exchangerm_ll_out-h_h_exchangerm_ll_in+(self._T_amb+273.15)*(
            s_h_exchangerm_ll_in-s_h_exchangerm_ll_out)
        #results["hot_exergy"] = self.ex_calculations(T_in+273.15,
        #                                             t_h_exchangerm_hh_out+273.15,
        #                                             a) 
        #results["cold_exergy"] = self.ex_calculations(t_h_exchangerm_ll_out+273.15,
        #                                              t_h_exchangerm_ll_in+273.15,
        #                                              b)
        results["lost_exergy"] = results["hot_exergy"]-results["cold_exergy"]
        results["ex_co2"]={}
        results["ex_co2"]["h_in"]=self.E(T_in,self._p_bray_H)
        results["ex_co2"]["h_out"]=self.E(t_h_exchangerm_hh_out,self._p_bray_H)
        results["ex_co2"]["l_in"]=self.E(t_h_exchangerm_ll_in,self._p_bray_L)
        results["ex_co2"]["l_out"]=self.E(t_h_exchangerm_ll_out,self._p_bray_L)
        if t_h_exchangerm_hh_out-t_h_exchangerm_ll_in >= self._min_temperature_exchange:
            results["main_exchanger"]=1
        else:
            results["main_exchanger"]=0
        return results
    
    def ex_calculations(self,a1,a2,Q):
        gcpw=(a1-a2-293.15*math.log(a1/a2))/(a1-a2)
        exergy = Q*gcpw
        return exergy
    def ex_calculations1(self,a1,Q):
        gcpw=1-(293.15)/(a1+273.15)
        exergy = Q*gcpw
        return exergy
    def turbine(self, T_in,P_in,P_out):
        t_turbine_in = T_in
        p_turbine_in = P_in
        p_turbine_out = P_out
        h_turbine_in = CP.PropsSI('H', 'T', t_turbine_in+273.15, 'P', p_turbine_in, "REFPROP::co2")
        s_turbine_in = CP.PropsSI('S', 'T', t_turbine_in+273.15, 'P', p_turbine_in, "REFPROP::co2")
        s_hypothesis = s_turbine_in
        h_hypothesis = CP.PropsSI('H', 'S', s_hypothesis, 'P', p_turbine_out, "REFPROP::co2")
        h_turbine_out = h_turbine_in+(h_hypothesis-h_turbine_in)*self._t_isentropic_eff_mc
        power_turbine = (h_hypothesis-h_turbine_in)*self._t_isentropic_eff_mc*self._mechanical_eff
        e_lost_turbine = h_turbine_out-h_turbine_in-power_turbine
        t_turbine_out = CP.PropsSI('T', 'H', h_turbine_out, 'P', p_turbine_out, "REFPROP::co2")-273.15
        s_turbine_out = CP.PropsSI('S', 'H', h_turbine_out, 'P', p_turbine_out, "REFPROP::co2")
        results = {}
        results["t_turbine_out"] = t_turbine_out 
        results["p_turbine_out"] = p_turbine_out
        results["h_turbine_out"] = h_turbine_out
        results["s_turbine_out"] = s_turbine_out
        results["power_turbine"] = -power_turbine
        results["e_lost_turbine"] = -e_lost_turbine
        results["exergy_lost_turbine"] = (self._T_amb+273.15)*(s_turbine_out-s_turbine_in)
        results["ex_co2"]={}
        results["ex_co2"]["in"]=self.E(t_turbine_in,p_turbine_in)
        results["ex_co2"]["out"]=self.E(t_turbine_out,p_turbine_out)
        return results
    
    def heat_recovery(self, T_in,P,flue_gas_name):
        t_heat_recovery_in = T_in
        h_heat_recovery_in = CP.PropsSI('H', 'T', t_heat_recovery_in+273.15, 'P', P , "REFPROP::co2")
        s_heat_recovery_in = CP.PropsSI('S', 'T', t_heat_recovery_in+273.15, 'P', P , "REFPROP::co2")
        t_heat_recovery_out = self._industrial_waste_heat_t-self._min_temperature_exchange
        h_heat_recovery_out = CP.PropsSI('H', 'T', t_heat_recovery_out+273.15, 'P', P , "REFPROP::co2")
        s_heat_recovery_out = CP.PropsSI('S', 'T', t_heat_recovery_out+273.15, 'P', P , "REFPROP::co2")
        h_heat_recovery_supply = h_heat_recovery_out-h_heat_recovery_in
        h_heat_recovery_receive=h_heat_recovery_supply/self._heat_transfer_loss_eff
        h_lost_heat_recovery1=h_heat_recovery_receive-h_heat_recovery_supply

        fluid=flue_gas_name
        t_flue_gas_in = self._industrial_waste_heat_t
        if T_in+self._min_temperature_exchange>=160 :
            t_flue_gas_out = T_in+self._min_temperature_exchange
        else :
            t_flue_gas_out = 160
        h_flue_gas_in = CP.PropsSI('H', 'T', t_flue_gas_in+273.15, 'P', self._p_amb , fluid)
        h_flue_gas_out = CP.PropsSI('H', 'T', t_flue_gas_out+273.15, 'P', self._p_amb , fluid)
        n=h_heat_recovery_receive/(h_flue_gas_in-h_flue_gas_out)

        results = {}
        results["t_flue_gas_in"] = t_flue_gas_in
        results["t_flue_gas_out"] = t_flue_gas_out
        results["n"] = n

        results["t_heat_recovery_out"] = t_heat_recovery_out 
        results["p_heat_recovery_out"] = P
        results["h_heat_recovery_out"] = h_heat_recovery_out
        results["s_heat_recovery_out"] = CP.PropsSI('S', 'T', t_heat_recovery_out+273.15, 'P', P, "REFPROP::co2")
        results["hot_heat_recovery"] = h_heat_recovery_receive
        results["h_lost_heat_recovery"] = h_lost_heat_recovery1
        results["hot_exergy"] = self.ex_calculations(t_flue_gas_in+273.15,
                                                     t_flue_gas_out+273.15,
                                                     h_heat_recovery_receive) 
        results["cold_exergy"]=h_heat_recovery_out-h_heat_recovery_in+(self._T_amb+273.15)*(
            s_heat_recovery_in-s_heat_recovery_out)
        #results["cold_exergy"] = self.ex_calculations(T_in+273.15,
        #                                              t_heat_recovery_out+273.15,
        #                                              h_heat_recovery_supply)
        results["lost_exergy"] = results["hot_exergy"]-results["cold_exergy"]
        results["ex_co2"]={}
        results["ex_co2"]["in"]=self.E(T_in,P)
        results["ex_co2"]["out"]=self.E(t_heat_recovery_out,P)
        return results
    
    def evaluation_indicators(self,results,P_in):
        eva={}
        eva["h_lost"]=results["heat_recovery"]["h_lost_heat_recovery"]+results["mian_h_exchanger"]["h_lost_exchangerm"]
        eva["e_lost"]=results["primary_compressor"]["e_lost_compressor"]+results["secondary_compressor"]["e_lost_compressor"]+results["turbine"]["e_lost_turbine"]
        eva["lost_benchmark"]=eva["h_lost"]+eva["e_lost"]
        eva["power_cost_benchmark"]=results["primary_compressor"]["power_compressor"]+results["secondary_compressor"]["power_compressor"]-results["turbine"]["power_turbine"]
        eva["hot_cost_benchmark"]=results["heat_recovery"]["hot_heat_recovery"]
        eva["hot_output_benchmark"]=results["primary_h_exchanger"]["hot_out_h_exchanger"]+results["secondary_h_exchanger"]["hot_out_h_exchanger"]
        eva["cop"]=eva["hot_output_benchmark"]/eva["power_cost_benchmark"]
        eva["mass_flow"]=P_in/eva["power_cost_benchmark"]
        eva["power_cost"]=P_in
        eva["hot_cost"]=eva["mass_flow"]*eva["hot_cost_benchmark"]
        eva["hot_output"]=eva["mass_flow"]*eva["hot_output_benchmark"]
        eva["lost"]=eva["mass_flow"]*eva["lost_benchmark"]
        eva["flue_gas_mass_flow"] = eva["mass_flow"]*results["heat_recovery"]["n"]
        eva["exergy_in"] = (P_in*1) +results["heat_recovery"]["hot_exergy"]*eva["mass_flow"]
        b1=results["primary_h_exchanger"]["Process_taste"]
        b2=results["secondary_h_exchanger"]["Process_taste"]
        eva["exergy_out"] = eva["mass_flow"]*(results["primary_h_exchanger"]["hot_out_h_exchanger"]*b1+results["secondary_h_exchanger"]["hot_out_h_exchanger"]*b2)
        eva["exergy_lost"]={}
        eva["exergy_lost"]["e_lost"]=eva["e_lost"]*eva["mass_flow"]
        eva["exergy_lost"]["p_Comp_lost"] = results["primary_compressor"]["exergy_lost_compressor"]*eva["mass_flow"]
        eva["exergy_lost"]["s_Comp_lost"] = results["secondary_compressor"]["exergy_lost_compressor"]*eva["mass_flow"]
        eva["exergy_lost"]["Turb_lost"] = results["turbine"]["exergy_lost_turbine"]*eva["mass_flow"]
        eva["exergy_lost"]["mian_h_exchanger"] = results["mian_h_exchanger"]["lost_exergy"]*eva["mass_flow"]
        eva["exergy_lost"]["heat_recovery"] = results["heat_recovery"]["lost_exergy"]*eva["mass_flow"]
        eva["exergy_lost_all"]=(eva["exergy_lost"]["e_lost"]+eva["exergy_lost"]["p_Comp_lost"]
                                +eva["exergy_lost"]["s_Comp_lost"]+eva["exergy_lost"]["Turb_lost"]
                                +eva["exergy_lost"]["mian_h_exchanger"]+eva["exergy_lost"]["heat_recovery"])

        return eva
    
if __name__ == '__main__':
    
    parameters = dict() 
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["isentropic_eff_mc"] = 0.88  #压缩机等熵效率
    parameters["t_isentropic_eff_mc"] = 0.92 #透平等熵效率
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["min_temperature_exchange"] = 15 #变量1，换热器最小换热温差
    parameters["industrial_waste_heat_t"] =300 #℃工业余热温度，变量3
    parameters["heat_transfer_loss_eff"] = 0.96 #换热器热损失系数
    parameters["t_amb"] = 20#环境温度
    parameters["p_amb"] = 101325#环境压力
    parameters["p_bray_L"] = 7.5e6#热泵低压
    parameters["Store_electrical_power"]= 1e6


    inputs={}
    inputs["p_bray_H"] = 19447839.26865841#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12827110.4341202 #优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 1e5 #变量4，反应器压力
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度

    a=505
    BraytonHBs = BraytonHeatPump(parameters)
    results = BraytonHBs.solve(inputs,a)
    print(results)
    print(results["primary_compressor"]["t_compressor_out"],results["secondary_compressor"]["t_compressor_out"],results["evaluation_indicators"]["cop"])