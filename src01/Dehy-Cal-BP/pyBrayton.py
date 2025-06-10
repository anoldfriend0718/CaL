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

class Brayton(object):
    def __init__(self, parameters) -> None:
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]#等熵效率
        self._t_isentropic_eff_mc = parameters["t_isentropic_eff_mc"]#透平等熵效率
        self._mechanical_eff = parameters["mechanical_eff"]#机械效率
        self._min_temperature_exchange = parameters["min_temperature_exchange"]#最小换热温差
        self._heat_transfer_loss_eff = parameters["heat_transfer_loss_eff"]#换热损失

        self._p_bray_L_B = parameters["p_bray_L_B"]
        self._p_amb = parameters["p_amb"]
        self._T_amb = parameters["t_amb"]
        self._T_L = 32
        self._hydrator_eff = parameters["hydrator_eff"]
        
    def solve(self,a,inputs):
        self._p_bray_H_B = inputs["p_bray_H_B"] 
        self._p_bray_M_B = inputs["p_bray_MH_B"] 
        self._p_bray_m_B = inputs["p_bray_ML_B"] 
        self._p_reaction = inputs["p_Hydr"]
        self.R = inputs["R"]
        self.H_out= inputs["H_out"]
        results = {}
        #Basic input data
        self._t_reaction_B = inputs["HT"]
        initialvalue = self.initialvalue()
        results["B_initialvalue"] = initialvalue
        #The primary Turbine is the starting point
        primary_turbine = self.turbine(self._t_reaction_B-self._min_temperature_exchange,
                                                self._p_bray_H_B,
                                                self._p_bray_M_B)
        results["B_primary_turbine"] = primary_turbine
        #secondary heat exchanger
        t_sec_h_in=results["B_primary_turbine"]["t_turbine_out"]
        secondary_h_exchanger=self.h_exchanger(t_sec_h_in,
                                                     self._t_reaction_B-self._min_temperature_exchange,
                                                     self._p_bray_M_B)
        results["secondary_h_exchanger"] = secondary_h_exchanger
        #The secondary Turbine 
        secondary_turbine = self.turbine(self._t_reaction_B-self._min_temperature_exchange,
                                                self._p_bray_M_B,
                                                self._p_bray_L_B)
        results["B_secondary_turbine"] = secondary_turbine
        #High heat_recovery 
        t_h_re_in=results["B_secondary_turbine"]["t_turbine_out"]

        #转到分流部分
        primary_compressor = self.compressor(self._T_L,
                                             self._p_bray_L_B,
                                             self._p_bray_m_B)
        results["B_primary_compressor"] = primary_compressor
        #间冷
        t_pri_com_out=results["B_primary_compressor"]["t_compressor_out"]
        intercooler = self.cooling_tower(t_pri_com_out,
                                           self._T_L,
                                           self._p_bray_m_B)
        results["intercooler"] = intercooler
        #Secondary compressors
        secondary_compressor = self.compressor(self._T_L,
                                             self._p_bray_m_B,
                                             self._p_bray_H_B)
        results["B_secondary_compressor"] = secondary_compressor
        t_h_out=results["B_secondary_compressor"]["t_compressor_out"]+self._min_temperature_exchange
        #低温回热器
        Low_recovery=self.l_recovery_h(self.H_out,
                                                     t_h_out,
                                                     self._p_bray_L_B,
                                                     self._p_bray_H_B)
        results["Low_recovery"] = Low_recovery
        cooling_tower1 = self.cooling_tower(t_h_out,
                                           self._T_L,
                                           self._p_bray_L_B)
        results["cooling_tower1"] = cooling_tower1
        #再压缩机
        R_compressor = self.compressor(t_h_out,
                                             self._p_bray_L_B,
                                             self._p_bray_H_B)
        results["R_compressor"] = R_compressor
        #合流点
        H_R=(Low_recovery["h_h_recovery_lh_out"]*(1-self.R)+
             R_compressor["h_compressor_out"]*self.R)
        t_R=CP.PropsSI('T', 'H', H_R, 'P', self._p_bray_H_B, "REFPROP::co2")-273.15
        results["t_R"]=t_R
        #高温换热器
        results["t_C"]=results["R_compressor"]["t_compressor_out"]
        High_h_recovery=self.h_recovery_h(t_h_re_in,
                                                     self.H_out,
                                                     t_R,
                                                     self._p_bray_L_B,
                                                     self._p_bray_H_B)
        results["High_h_recovery"] = High_h_recovery
        #primary heat exchanger
        t_p_in=results["High_h_recovery"]["t_h_recovery_lh_out"]
        primary_h_exchanger=self.h_exchanger(t_p_in,
                                                     self._t_reaction_B-self._min_temperature_exchange,
                                                     self._p_bray_H_B)
        results["primary_h_exchanger"] = primary_h_exchanger

        evaluation_indicators = self.evaluation_indicators(results,a)
        results["evaluation_indicators"] = evaluation_indicators
        return results

    def initialvalue(self):
        results = {}
        results["isentropic_eff_mc"] = self._isentropic_eff_mc
        results["mechanical_eff"] = self._mechanical_eff
        results["min_temperature_exchange"]= self._min_temperature_exchange 
        results["t_reaction_B"] = self._t_reaction_B
        results["p_bray_H_B"] = self._p_bray_H_B
        results["p_bray_M_B"] = self._p_bray_M_B
        results["p_bray_L_B"] = self._p_bray_L_B
        return results
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
        results["hot_in_h_exchanger"] = -hot_out_h_exchanger
        a1=T_in+273.15
        a2=T_out+273.15
        results["Process_taste"]=(a1-a2-293.15*math.log(a1/a2))/(a1-a2)
        return results
    def h_recovery_h(self,T_in,T_out,t_in,Pl,Ph):
        t_h_exchangerm_hl_in=T_in#低压入口温度
        t_h_exchangerm_hl_out = T_out#低压出口温度
        t_h_exchangerm_lh_in=t_in#高压入口温度

        h_h_exchangerm_hl_in =CP.PropsSI('H', 'T', t_h_exchangerm_hl_in+273.15,  'P', Pl, "REFPROP::co2")
        s_h_exchangerm_hl_in =CP.PropsSI('S', 'T', t_h_exchangerm_hl_in+273.15,  'P', Pl, "REFPROP::co2")
        h_h_exchangerm_hl_out=CP.PropsSI('H', 'T', t_h_exchangerm_hl_out+273.15,  'P', Pl, "REFPROP::co2")
        s_h_exchangerm_hl_out=CP.PropsSI('S', 'T', t_h_exchangerm_hl_out+273.15, 'P', Pl, "REFPROP::co2")
        hot_hl_exchange=h_h_exchangerm_hl_in-h_h_exchangerm_hl_out
        hot_lh_exchange= hot_hl_exchange*self._heat_transfer_loss_eff

        h_h_exchangerm_lh_in =CP.PropsSI('H', 'T', t_h_exchangerm_lh_in+273.15,  'P', Ph, "REFPROP::co2")
        s_h_exchangerm_lh_in =CP.PropsSI('S', 'T', t_h_exchangerm_lh_in+273.15,  'P', Ph, "REFPROP::co2")
        h_h_exchangerm_lh_out=h_h_exchangerm_lh_in+hot_lh_exchange
        t_h_exchangerm_lh_out=CP.PropsSI('T', 'H', h_h_exchangerm_lh_out,  'P', Ph, "REFPROP::co2")-273.15
        s_h_exchangerm_lh_out=CP.PropsSI('S', 'T', t_h_exchangerm_lh_out+273.15, 'P', Ph, "REFPROP::co2")
        results = {}
        results["t_h_recovery_hl_out"] = t_h_exchangerm_hl_out 
        results["p_h_recovery_hl_out"] = self._p_bray_L_B
        results["h_h_recovery_hl_out"] = h_h_exchangerm_hl_out
        results["s_h_recovery_hl_out"] = s_h_exchangerm_hl_out
        results["t_h_recovery_lh_out"] = t_h_exchangerm_lh_out
        results["p_h_recovery_lh_out"] = self._p_bray_H_B
        results["h_h_recovery_lh_out"] = h_h_exchangerm_lh_out
        results["s_h_recovery_lh_out"] = s_h_exchangerm_lh_out
        results["h_lost_recovery"] = hot_hl_exchange-hot_lh_exchange
        results["hot_exergy"]=h_h_exchangerm_hl_in-h_h_exchangerm_hl_out+(self._T_amb+273.15)*(
            s_h_exchangerm_hl_out-s_h_exchangerm_hl_in)
        results["cold_exergy"]=h_h_exchangerm_lh_out-h_h_exchangerm_lh_in+(self._T_amb+273.15)*(
            s_h_exchangerm_lh_in-s_h_exchangerm_lh_out)
        results["lost_exergy"] = results["hot_exergy"]-results["cold_exergy"]
        results["ex_co2"]={}
        results["ex_co2"]["h_in"]=self.E(T_in,Pl)
        results["ex_co2"]["h_out"]=self.E(T_out,Pl)
        results["ex_co2"]["l_in"]=self.E(t_h_exchangerm_lh_in,Ph)
        results["ex_co2"]["l_out"]=self.E(t_h_exchangerm_lh_out,Ph)
        return results
    def l_recovery_h(self,T_in,T_out,Pl,Ph):
        t_h_exchangerm_hl_in=T_in#低压入口温度
        t_h_exchangerm_hl_out = T_out#低压出口温度
        t_h_exchangerm_lh_in=t_h_exchangerm_hl_out-self._min_temperature_exchange#高压入口温度

        h_h_exchangerm_hl_in =CP.PropsSI('H', 'T', t_h_exchangerm_hl_in+273.15,  'P', Pl, "REFPROP::co2")
        s_h_exchangerm_hl_in =CP.PropsSI('S', 'T', t_h_exchangerm_hl_in+273.15,  'P', Pl, "REFPROP::co2")
        h_h_exchangerm_hl_out=CP.PropsSI('H', 'T', t_h_exchangerm_hl_out+273.15,  'P', Pl, "REFPROP::co2")
        s_h_exchangerm_hl_out=CP.PropsSI('S', 'T', t_h_exchangerm_hl_out+273.15, 'P', Pl, "REFPROP::co2")
        hot_hl_exchange=h_h_exchangerm_hl_in-h_h_exchangerm_hl_out#热流
        hot_lh_exchange= hot_hl_exchange*self._heat_transfer_loss_eff

        h_h_exchangerm_lh_in =CP.PropsSI('H', 'T', t_h_exchangerm_lh_in+273.15,  'P', Ph, "REFPROP::co2")
        s_h_exchangerm_lh_in =CP.PropsSI('S', 'T', t_h_exchangerm_lh_in+273.15,  'P', Ph, "REFPROP::co2")
        h_h_exchangerm_lh_out=h_h_exchangerm_lh_in+hot_lh_exchange/(1-self.R)
        t_h_exchangerm_lh_out=CP.PropsSI('T', 'H', h_h_exchangerm_lh_out,  'P', Ph, "REFPROP::co2")-273.15
        s_h_exchangerm_lh_out=CP.PropsSI('S', 'T', t_h_exchangerm_lh_out+273.15, 'P', Ph, "REFPROP::co2")
        results = {}
        results["t_h_recovery_hl_out"] = t_h_exchangerm_hl_out 
        results["p_h_recovery_hl_out"] = self._p_bray_L_B
        results["h_h_recovery_hl_out"] = h_h_exchangerm_hl_out
        results["s_h_recovery_hl_out"] = s_h_exchangerm_hl_out
        results["t_h_recovery_lh_out"] = t_h_exchangerm_lh_out
        results["p_h_recovery_lh_out"] = self._p_bray_H_B
        results["h_h_recovery_lh_out"] = h_h_exchangerm_lh_out
        results["s_h_recovery_lh_out"] = s_h_exchangerm_lh_out
        results["h_lost_recovery"] = hot_hl_exchange-hot_lh_exchange
        results["hot_exergy"]=h_h_exchangerm_hl_in-h_h_exchangerm_hl_out+(self._T_amb+273.15)*(
            s_h_exchangerm_hl_out-s_h_exchangerm_hl_in)
        results["cold_exergy"]=h_h_exchangerm_lh_out-h_h_exchangerm_lh_in+(self._T_amb+273.15)*(
            s_h_exchangerm_lh_in-s_h_exchangerm_lh_out)*(1-self.R)
        results["lost_exergy"] = results["hot_exergy"]-results["cold_exergy"]
        results["ex_co2"]={}
        results["ex_co2"]["h_in"]=self.E(T_in,Pl)
        results["ex_co2"]["h_out"]=self.E(T_out,Pl)
        results["ex_co2"]["l_in"]=self.E(t_h_exchangerm_lh_in,Ph)*(1-self.R)
        results["ex_co2"]["l_out"]=self.E(t_h_exchangerm_lh_out,Ph)*(1-self.R)
        return results
     
    def cooling_tower(self,T_in,T_out,P):
        results = {}
        results["t_cooling_tower_in"] = T_in
        results["p_cooling_tower_in"] = P
        results["h_cooling_tower_in"] = CP.PropsSI('H', 'T', T_in+273.15, 'P', P, "REFPROP::co2")
        results["s_cooling_tower_in"] = CP.PropsSI('S', 'T', T_in+273.15, 'P', P, "REFPROP::co2")
        results["t_cooling_tower_out"] = T_out
        results["p_cooling_tower_out"] = P
        results["h_cooling_tower_out"] = CP.PropsSI('H', 'T', T_out+273.15, 'P', P, "REFPROP::co2")
        results["s_cooling_tower_out"] = CP.PropsSI('S', 'T', T_out+273.15, 'P', P, "REFPROP::co2")
        results["hot_cooling_tower"] = (results["h_cooling_tower_in"]
                                        -results["h_cooling_tower_out"])
        
        results["exergy_lost"] =((results["h_cooling_tower_in"]-results["h_cooling_tower_out"])
                                 -(self._T_amb+273.15)*(results["s_cooling_tower_in"]-
                                                        results["s_cooling_tower_out"]))
        results["ex_co2"]={}
        results["ex_co2"]["in"]=self.E(T_in,P)
        results["ex_co2"]["out"]=self.E(T_out,P)
        return results

    def E(self,T,P):
        H1=CP.PropsSI('H', 'T', T+273.15, 'P', P, "REFPROP::co2")
        S1=CP.PropsSI('S', 'T', T+273.15, 'P', P, "REFPROP::co2")
        H0=CP.PropsSI('H', 'T', self._T_amb+273.15, 'P', self._p_amb, "REFPROP::co2")
        S0=CP.PropsSI('S', 'T', self._T_amb+273.15, 'P', self._p_amb, "REFPROP::co2")
        a=H1-H0-(self._T_amb+273.15)*(S1-S0)
        return a
    def ex_calculations(self,a1,a2,Q):
        gcpw=(a1-a2-293.15*math.log(a1/a2))/(a1-a2)
        exergy = Q*gcpw
        return exergy
    def ex_calculations1(self,a1,Q):
        gcpw=1-(293.15)/(a1+273.15)
        exergy = Q*gcpw
        return exergy
    
    def evaluation_indicators(self,results,H_in):
        eva={}
        eva["hydrator_lost"] = H_in*(1-self._hydrator_eff)
 
        eva["power_benchmark"]=(-results["B_primary_compressor"]["power_compressor"]*(1-self.R)
                                -results["B_secondary_compressor"]["power_compressor"]*(1-self.R)
                                -results["R_compressor"]["power_compressor"]*self.R
                                +results["B_primary_turbine"]["power_turbine"]
                                +results["B_secondary_turbine"]["power_turbine"])
        eva["hot_in_benchmark"]=(results["primary_h_exchanger"]["hot_in_h_exchanger"]
                                 +results["secondary_h_exchanger"]["hot_in_h_exchanger"])
        
        eva["mass_flow"]=(H_in*self._hydrator_eff)/eva["hot_in_benchmark"]
        eva["re_Heat_in"]=H_in*self._hydrator_eff#循环用热
        eva["power"]=eva["mass_flow"]*eva["power_benchmark"]#循环产电
        eva["cooling_lost"]=results["cooling_tower1"]["hot_cooling_tower"]*eva["mass_flow"]


        eva["exergy"]={}
        eva["exergy"]["re_Heat_in"]=(results["primary_h_exchanger"]["hot_in_h_exchanger"]*results["primary_h_exchanger"]["Process_taste"]
                                 +results["secondary_h_exchanger"]["hot_in_h_exchanger"]*results["secondary_h_exchanger"]["Process_taste"])*eva["mass_flow"]
        eva["exergy"]["power"]=eva["power"]
        eva["exergy_efficiency"]=eva["exergy"]["power"]/eva["exergy"]["re_Heat_in"]
        eva["energy_efficiency"]=eva["power"]/eva["re_Heat_in"]
        return eva


if __name__ == '__main__':
    
    parameters = dict() 
    
    parameters["isentropic_eff_mc"] = 0.7  #等熵效率
    parameters["t_isentropic_eff_mc"] = 0.85
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["min_temperature_exchange"] = 15

    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_amb"] = 20
    parameters["p_amb"] = 101325

    parameters["hydrator_eff"] = 0.95   #水合器器效率
    parameters["p_bray_L_B"] = 7.885e6

    BraytonHBs = Brayton(parameters)
    Hydrator_heat=	1498957.696586326-30766.434959109065
    inputs={}
    inputs["p_bray_H_B"] = 30e6
    inputs["p_bray_MH_B"] = 16217752.142109105
    inputs["p_bray_ML_B"] = 12217752.142109105
    inputs["p_Hydr"] = 1e5
    inputs["R"]=0.37
    inputs["HT"] = 650
    inputs["H_out"]=214

    results = BraytonHBs.solve(Hydrator_heat,inputs)
    print(results)