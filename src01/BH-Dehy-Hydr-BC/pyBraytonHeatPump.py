import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")
import math
import numpy as np
import pandas as pd
import json
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from Cp0massWrapper import Cp0mass_Wrapper
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

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
        self._industrial_waste_heat_t=parameters["industrial_waste_heat_t"]#工业余热温度
        self._heat_transfer_loss_eff=parameters["heat_transfer_loss_eff"]#换热损失
        self._p_bray_L = parameters["p_bray_L"]
        self._p_amb = parameters["p_amb"]
        self._T_amb = parameters["t_amb"]
        self._Store_electrical_power=parameters["Store_electrical_power"]
        
    def solve(self,inputs):
        results = {}
        self._p_bray_H = inputs["p_bray_H"] 
        self._p_bray_M =inputs["p_bray_M"] 
        self._p_reaction =inputs["p_Dehy"]
        self._Dehy_ot = inputs["Dehy_overheating_temperature"] 
        self._min_temperature_exchange=inputs["min_temperature_exchange"]#最小换热温差
        #self._p_bray_H = inputs["p_bray_H"]
        #Basic input data
        initialvalue = self.initialvalue()
        results["initialvalue"] = initialvalue
        t_equilibrium = self.equilibrium()
        self._t_reaction = t_equilibrium+self._Dehy_ot
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
        results["cost"] = self.cost(results,inputs["Economic Model Selection"],
                                    inputs["Compressor power limit"],
                                    inputs["Turbine power limit"])
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
        results["f"] = ((t_compressor_out-550)**2)*5.4/100000+1
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
        t1=T_in-525
        t2=T_out-525
        t=(t1-t2)/(math.log(t1/t2))
        results["A"] = hot_out_h_exchanger/(t*300)
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
        t1=t_h_exchangerm_hh_in-t_h_exchangerm_ll_out
        t2=t_h_exchangerm_hh_out-t_h_exchangerm_ll_in
        t=(t1-t2)/(math.log(t1/t2))
        results["A"] = hot_hh_exchange/(t*300)
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
        results["f"] = 1
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
        t1=self._industrial_waste_heat_t-t_heat_recovery_out
        t2=t_flue_gas_out-t_heat_recovery_in
        t=(t1+t2)/2
        results["A"] = h_heat_recovery_receive/(t*300)
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
        
    def cost(self,results,a,b,c):
        a22=results["evaluation_indicators"]["mass_flow"]
        t_h_e1_in=results["primary_compressor"]["t_compressor_out"]
        t_h_e2_in=results["secondary_compressor"]["t_compressor_out"]
        t_h_em_hh_in=results["secondary_h_exchanger"]["t_h_exchanger_out"]
        t_p_tur_in=results["mian_h_exchanger"]["t_h_exchangerm_hh_out"]
        result={}
        a_pc=results["primary_compressor"]["power_compressor"]*a22/1000 #转kW
        a_sc=results["secondary_compressor"]["power_compressor"]*a22/1000 #转kW
        a_t=results["turbine"]["power_turbine"]*a22/1000 #转kW
        a_pe=results["primary_h_exchanger"]["hot_out_h_exchanger"]*a22
        t_pe=self.lmtd(t_h_e1_in,self._t_reaction+self._min_temperature_exchange,self._t_reaction,self._t_reaction)
        a_se=results["secondary_h_exchanger"]["hot_out_h_exchanger"]*a22
        t_se=self.lmtd(t_h_e2_in,self._t_reaction+self._min_temperature_exchange,self._t_reaction,self._t_reaction)
        a_me=results["mian_h_exchanger"]["h_lost_exchangerm"]/(1-self._heat_transfer_loss_eff)*a22
        t_me=self.lmtd(t_h_em_hh_in,t_p_tur_in,
                       self._industrial_waste_heat_t-self._min_temperature_exchange,
                       t_h_em_hh_in-self._min_temperature_exchange)
        a_hr=results["heat_recovery"]["hot_heat_recovery"]*a22
        t_hr=self.lmtd(results["heat_recovery"]["t_flue_gas_in"],
                       results["heat_recovery"]["t_flue_gas_out"],
                       results["turbine"]["t_turbine_out"],
                       results["heat_recovery"]["t_heat_recovery_out"])
        if a == 1:
            result["cost_pc"]=7331*a_pc**0.7865
            result["cost_sc"]=7331*a_sc**0.7865  
            result["cost_t"]=8279*a_t**0.6842 
            result["cost_pe"]=a_pe/t_pe* self.cHe(a_pe/t_pe)
            result["cost_se"]=a_se/t_se* self.cHe(a_se/t_se)
            result["cost_me"]=a_me/t_me* self.cRe(a_me/t_me)
            result["cost_hr"]=a_hr/t_hr* self.cHe(a_hr/t_hr)
        elif a == 2:
            pc_value = results["primary_compressor"]["power_compressor"]*a22
            pc_f=pc_value/1e6
            result["cost_pc"] = 1230000*pc_f**0.3992 if b > pc_value else 1230000*pc_f**0.3992*0.9814
            sc_value = results["secondary_compressor"]["power_compressor"]*a22
            sc_f=sc_value/1e6
            result["cost_sc"] = 1230000*sc_f**0.3992 if b > sc_value else 1230000*sc_f**0.3992*0.9814
            t_value=results["turbine"]["power_turbine"]*a22
            t_ff=results["turbine"]["f"]
            t_f=t_value/1e6
            result["cost_t"] = 406200*t_f**0.8*t_ff if c > t_value else 182600*t_f**0.5561*t_ff*0.9814
            result["cost_pe"]=a_pe/t_pe* self.cHe(a_pe/t_pe)
            result["cost_se"]=a_se/t_se* self.cHe(a_se/t_se)
            result["cost_me"]=351.81*(a_me/t_me)**0.7544/7.2492
            result["cost_hr"]=32.88*(a_hr/t_hr)**0.75*0.9814
        result["cost_comp&turb"]=result["cost_pc"]+result["cost_sc"]+result["cost_t"]
        result["cost_exchanger"]=result["cost_pe"]+result["cost_se"]+result["cost_me"]+result["cost_hr"]
        result["cost_all"]=result["cost_comp&turb"]+result["cost_exchanger"]
        return result

    def cHe(self, Xa):
        # 初始化输入和输出数据
        X = np.array([5000, 30000, 100000, 300000, 1000000,1000000000]).astype(float)  # 输入数据需转换为二维数组
        Y = np.array([1.9, 1.3, 1.1, 1, 1,1]).astype(float)  # 输出数据
        # 创建插值函数对象，这里使用线性插值
        f_linear = interp1d(X, Y, kind='linear')
        # 使用插值函数进行预测
        Y_pred = f_linear(Xa)
        return Y_pred
    
    def cRe(self, Xa):
        # 初始化输入和输出数据
        X = np.array([5000, 30000, 100000, 300000, 1000000,1000000000]).astype(float)  # 输入数据需转换为二维数组
        Y = np.array([5.89, 1.31, 1.22, 1.03, 0.94,0.94]).astype(float)  # 输出数据
        f_linear = interp1d(X, Y, kind='linear')
        Y_pred = f_linear(Xa)
        return Y_pred
    
    def cAir(self, Xa):
        # 初始化输入和输出数据
        X = np.array([5000, 30000, 100000, 300000, 1000000,1000000000]).astype(float)  # 输入数据需转换为二维数组
        Y = np.array([7.6, 2.4 , 1.3, 1.1, 1,1]).astype(float)  # 输出数据
        f_linear = interp1d(X, Y, kind='linear')
        Y_pred = f_linear(Xa)
        return Y_pred
 
    def lmtd(self,T1, T2, t1, t2):
        if (T1 - t2) - (T2 - t1)==0:
            LMTD = (T1 - t2)
        else:
            numerator = (T1 - t2) - (T2 - t1)
            denominator = math.log((T1 - t2) / (T2 - t1))
            LMTD = numerator / denominator
        return LMTD
    
    def cost2(self,results):
        a22=results["evaluation_indicators"]["mass_flow"]
        t_h_e1_in=results["primary_compressor"]["t_compressor_out"]
        t_h_e2_in=results["secondary_compressor"]["t_compressor_out"]
        t_h_em_hh_in=results["secondary_h_exchanger"]["t_h_exchanger_out"]
        t_p_tur_in=results["mian_h_exchanger"]["t_h_exchangerm_hh_out"]
        result={}
        a_pc=self.r_p(self._t_reaction,self._p_bray_L,t_h_e1_in,self._p_bray_M)
        result["1"] = a_pc**11
        b_pc=a22*self.r_m(self._t_reaction,self._p_bray_L)
        result["cost_pc"]=b_pc*47.1/(1-self._isentropic_eff_mc)*result["1"]*math.log(result["1"])*1.13
        a_sc=89
        b_sc=a22*self.r_m(self._t_reaction+self._min_temperature_exchange,self._p_bray_M)
        result["cost_sc"]=b_sc*47.1/(1-self._isentropic_eff_mc)*a_sc*math.log(a_sc)*1.13
        a_t=7
        b_t=a22*self.r_m(t_p_tur_in,self._p_bray_H)
        result["cost_t"]=b_t*392.2/(1-self._t_isentropic_eff_mc)*a_t*math.log(a_t)*(
            1+math.exp(0.036*(t_p_tur_in+273.15)-65.66))*1.13
        a_pe=results["primary_h_exchanger"]["hot_out_h_exchanger"]*a22/1000
        result["cost_pe"]=a_pe*32.3*1.13
        a_se=results["secondary_h_exchanger"]["hot_out_h_exchanger"]*a22/1000
        result["cost_se"]=a_se*32.3*1.13
        a_me=results["mian_h_exchanger"]["h_lost_exchangerm"]/(1-self._heat_transfer_loss_eff)*a22
        t_me=self.lmtd(t_h_em_hh_in,t_p_tur_in,
                       self._industrial_waste_heat_t-self._min_temperature_exchange,
                       t_h_em_hh_in-self._min_temperature_exchange)
        A_me = a_me / (1700 * t_me)
        result["cost_me"]=2546.9*(A_me**0.67)*((self._p_bray_H/10e5)**0.28)*1.13
        a_hr=results["heat_recovery"]["hot_heat_recovery"]*a22/1000
        result["cost_hr"]=a_hr*32.3*1.13
        result["cost_comp&turb"]=result["cost_pc"]+result["cost_sc"]+result["cost_t"]
        result["cost_exchanger"]=result["cost_pe"]+result["cost_se"]+result["cost_me"]+result["cost_hr"]
        result["cost_all"]=result["cost_comp&turb"]+result["cost_exchanger"]
        return result
    
    def r_m(self,t,p):
        # 定义温度和压力（单位：K和Pa）
        temperature = t+273.15  # 例如，300 K（约27°C）
        pressure = p  # 例如，101325 Pa（约1 atm）
        
        # 计算空气的密度（使用'Air'作为制冷剂名称）
        # 注意：CoolProp 中的'Air'模型是基于NIST的REFPROP数据，它可能不完全等同于实际的空气组成
        density_air = CP.PropsSI('D', 'T', temperature, 'P', pressure, 'Air')
        
        # 计算二氧化碳的密度（使用'CO2'作为制冷剂名称）
        density_co2 = CP.PropsSI('D', 'T', temperature, 'P', pressure, 'CO2')
        r_m = density_air/density_co2
        return r_m

    def r_p(self,t,p,T,P):
        r=self.r_m(t,p)
        fluid = 'Air'
        # 给定的密度（单位：kg/m^3）和温度（单位：K）
        target_density =r* CP.PropsSI('D', 'T', T+273.15, 'P', P, 'CO2')  # 示例密度，单位：kg/m^3
        given_temperature = T+273.15 # 示例温度，单位：K
        # 使用fsolve来找到使密度差异为零的压力值
        # 初始猜测的压力值（这里使用标准大气压作为示例）
        initial_guess_pressure = P  # 单位：Pa
        # 调用fsolve进行迭代求解
        solution = fsolve(self.density_difference, initial_guess_pressure, args=(given_temperature, target_density, fluid))
        # 获取求解得到的压力值
        calculated_pressure = solution[0]

        r_p=calculated_pressure/p
        return r_p

    def density_difference(self,pressure, temperature, target_density, fluid):
        # 计算在给定压力和温度下的密度
        calculated_density = CP.PropsSI('D', 'T', temperature, 'P', pressure, fluid)
        # 返回计算密度与目标密度之间的差异
        return calculated_density - target_density

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
    
    parameters["industrial_waste_heat_t"] =300 #℃工业余热温度，变量3
    parameters["heat_transfer_loss_eff"] = 0.96 #换热器热损失系数
    parameters["t_amb"] = 20#环境温度
    parameters["p_amb"] = 101325#环境压力
    parameters["p_bray_L"] = 7.5e6#热泵低压
    parameters["Store_electrical_power"]= 50e6


    inputs={}
    inputs["p_bray_H"] = 19447839.26865841#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12827110.4341202 #优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 1e5 #变量4，反应器压力
    inputs["Economic Model Selection"] = 1 #经济模型选择，1Tesio，2Nathan T
    inputs["Compressor power limit"] = 200e6#功率界限，影响齿轮离心和滚筒离心模型的选取，单位W，桶式离心需体积流量
    inputs["Turbine power limit"] = 35e6
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度
    inputs["min_temperature_exchange"] = 15 #变量1，换热器最小换热温差

    BraytonHBs = BraytonHeatPump(parameters)
    results = BraytonHBs.solve(inputs)
    data_for_json = {key: (value.item() if isinstance(value, np.floating) else value) for key, value in results.items()}
    json_output = json.dumps(data_for_json, indent=4)
    print(json_output)
    