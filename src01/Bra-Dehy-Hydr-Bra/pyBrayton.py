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
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
        self._flue_gas_composition = self._pw._norm_flue_gas_composition
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]#等熵效率
        self._t_isentropic_eff_mc = parameters["t_isentropic_eff_mc"]#透平等熵效率
        self._mechanical_eff = parameters["mechanical_eff"]#机械效率
        self._min_temperature_exchange = parameters["min_temperature_exchange"]#最小换热温差
        self._industrial_waste_heat_t = parameters["industrial_waste_heat_t"]#工业余热温度
        self._heat_transfer_loss_eff = parameters["heat_transfer_loss_eff"]#换热损失
        self._t_reaction = parameters["t_reaction"]
        self._p_bray_H = parameters["p_bray_H"]
        self._p_bray_M = parameters["p_bray_M"]
        self._p_bray_L = parameters["p_bray_L"]
        self._p_amb = parameters["p_amb"]
        self._T_amb = parameters["T_amb"]
        self._T_L = 32
        
    def solve(self,inputs):
        results = {}
        #Basic input data
        initialvalue = self.initialvalue()
        results["B_initialvalue"] = initialvalue
        #The primary Turbine is the starting point
        primary_turbine = self.turbine(self._t_reaction-self._min_temperature_exchange,
                                                self._p_bray_H,
                                                self._p_bray_M)
        results["B_primary_turbine"] = primary_turbine
        #secondary heat exchanger
        t_sec_h_in=results["B_primary_turbine"]["t_turbine_out"]
        secondary_h_exchanger=self.h_exchanger(t_sec_h_in,
                                                     self._t_reaction-self._min_temperature_exchange,
                                                     self._p_bray_M)
        results["secondary_h_exchanger"] = secondary_h_exchanger
        #The secondary Turbine 
        secondary_turbine = self.turbine(self._t_reaction-self._min_temperature_exchange,
                                                self._p_bray_M,
                                                self._p_bray_L)
        results["B_secondary_turbine"] = secondary_turbine
        #High heat_recovery 
        t_h_re_in=results["B_secondary_turbine"]["t_turbine_out"]
        High_h_recovery=self.h_recovery(t_h_re_in,
                                                     self._industrial_waste_heat_t,
                                                     self._p_bray_L,
                                                     self._p_bray_H)
        results["High_h_recovery"] = High_h_recovery

        #primary heat exchanger
        t_p_in=results["High_h_recovery"]["t_h_recovery_hh_out"]
        primary_h_exchanger=self.h_exchanger(t_p_in,
                                                     self._t_reaction-self._min_temperature_exchange,
                                                     self._p_bray_H)
        results["primary_h_exchanger"] = primary_h_exchanger
        #primary compressors
        primary_compressor = self.compressor(self._T_L,
                                             self._p_bray_L,
                                             self._p_bray_M)
        results["B_primary_compressor"] = primary_compressor
        #Secondary compressors
        t_pri_com_out=results["B_primary_compressor"]["t_compressor_out"]
        secondary_compressor = self.compressor(t_pri_com_out,
                                             self._p_bray_M,
                                             self._p_bray_H)
        results["B_secondary_compressor"] = secondary_compressor
        #Low heat_recovery 
        t_l_out=results["B_secondary_compressor"]["t_compressor_out"]+self._min_temperature_exchange
        Low_h_recovery=self.h_recovery(self._industrial_waste_heat_t,
                                                     t_l_out,
                                                     self._p_bray_L,
                                                     self._p_bray_H)
        results["Low_h_recovery"] = Low_h_recovery
        #Cooling towers
        cooling_tower = self.cooling_tower(t_l_out,
                                           self._T_L,
                                           self._p_bray_L)
        results["cooling_tower"] = cooling_tower

        #Industrial waste heat heating part
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        t_p_rec_in=results["Low_h_recovery"]["t_h_recovery_hh_out"]
        heat_recovery = self.heat_recovery(t_p_rec_in,
                                               self._p_bray_H,
                                               flue_gas_name)
        results["heat_recovery"] = heat_recovery
        #Data synthesis
        heat_in=inputs
        evaluation_indicators = self.evaluation_indicators(results,heat_in)
        results["evaluation_indicators"] = evaluation_indicators
        return results
    
    def initialvalue(self):
        results = {}
        results["isentropic_eff_mc"] = self._isentropic_eff_mc
        results["mechanical_eff"] = self._mechanical_eff
        results["min_temperature_exchange"]= self._min_temperature_exchange 
        results["industrial_waste_heat_t"] = self._industrial_waste_heat_t
        results["t_reaction"] = self._t_reaction
        results["p_bray_H"] = self._p_bray_H
        results["p_bray_M"] = self._p_bray_M
        results["p_bray_L"] = self._p_bray_L
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
        #results["t_compressor_in"] = t_compressor_in 
        #results["p_compressor_in"] = p_compressor_in
        #results["h_compressor_in"] = h_compressor_in
        #results["s_compressor_in"] = s_compressor_in
        results["t_compressor_out"] = t_compressor_out 
        results["p_compressor_out"] = p_compressor_out
        results["h_compressor_out"] = h_compressor_out
        results["s_compressor_out"] = s_compressor_out
        results["power_compressor"] = power_compressor
        results["e_lost_compressor"] = e_lost_compressor
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
        results["hot_out_h_exchanger"] = hot_out_h_exchanger
        return results
    def h_recovery(self,T_in,T_out,Pl,Ph):
        t_h_exchangerm_hh_in=T_out-self._min_temperature_exchange
        t_h_exchangerm_ll_in=T_in
        t_h_exchangerm_ll_out=T_out
        h_h_exchangerm_ll_in =CP.PropsSI('H', 'T', t_h_exchangerm_ll_in+273.15,  'P', Pl, "REFPROP::co2")
        h_h_exchangerm_ll_out=CP.PropsSI('H', 'T', t_h_exchangerm_ll_out+273.15, 'P', Pl, "REFPROP::co2")
        s_h_exchangerm_ll_in =CP.PropsSI('S', 'T', t_h_exchangerm_ll_in+273.15,  'P', Pl, "REFPROP::co2")
        hot_ll_exchange=h_h_exchangerm_ll_in-h_h_exchangerm_ll_out
        hot_hh_exchange= hot_ll_exchange*self._heat_transfer_loss_eff
        h_h_exchangerm_hh_in =CP.PropsSI('H', 'T', t_h_exchangerm_hh_in+273.15,  'P', Ph, "REFPROP::co2")
        h_h_exchangerm_hh_out=h_h_exchangerm_hh_in+hot_hh_exchange
        t_h_exchangerm_hh_out=CP.PropsSI('T', 'H', h_h_exchangerm_hh_out,  'P', Ph, "REFPROP::co2")-273.15
        s_h_exchangerm_hh_out=CP.PropsSI('S', 'T', t_h_exchangerm_hh_out+273.15, 'P', Ph, "REFPROP::co2")
        results = {}
        results["t_h_recovery_hh_out"] = t_h_exchangerm_hh_out 
        results["p_h_recovery_hh_out"] = self._p_bray_H
        results["h_h_recovery_hh_out"] = h_h_exchangerm_hh_out
        results["s_h_recovery_hh_out"] = s_h_exchangerm_hh_out
        results["t_h_recovery_ll_in"] = t_h_exchangerm_ll_in 
        results["p_h_recovery_ll_in"] = self._p_bray_L
        results["h_h_recovery_ll_in"] = h_h_exchangerm_ll_in
        results["s_h_recovery_ll_in"] = s_h_exchangerm_ll_in
        results["h_lost_recovery"] = hot_ll_exchange-hot_hh_exchange
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
        results["hot_cooling_tower"] = results["h_cooling_tower_in"]-results["h_cooling_tower_out"]
        return results
    def heat_recovery(self, T_in,P,flue_gas_name):
        t_heat_recovery_in = T_in
        h_heat_recovery_in = CP.PropsSI('H', 'T', t_heat_recovery_in+273.15, 'P', P , "REFPROP::co2")
        t_heat_recovery_out = self._industrial_waste_heat_t-self._min_temperature_exchange
        h_heat_recovery_out = CP.PropsSI('H', 'T', t_heat_recovery_out+273.15, 'P', P , "REFPROP::co2")
        h_heat_recovery_supply = h_heat_recovery_out-h_heat_recovery_in
        h_heat_recovery_receive=h_heat_recovery_supply/self._heat_transfer_loss_eff
        h_lost_heat_recovery1=h_heat_recovery_receive-h_heat_recovery_supply

        fluid=flue_gas_name
        t_flue_gas_in = self._industrial_waste_heat_t
        t_flue_gas_out = T_in+self._min_temperature_exchange
        h_flue_gas_in = CP.PropsSI('H', 'T', t_flue_gas_in+273.15, 'P', self._p_amb , fluid)
        h_flue_gas_out = CP.PropsSI('H', 'T', t_flue_gas_out+273.15, 'P', self._p_amb , fluid)
        h_heat_recovery_in = (h_heat_recovery_supply/0.985)/0.975
        n=h_heat_recovery_in/(h_flue_gas_in-h_flue_gas_out)

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
        return results
    
    def evaluation_indicators(self,results,H_in):
        eva={}
        eva["hydrator_lost"] = H_in*0.05
        eva["h_lost_benchmark"]=results["High_h_recovery"]["h_lost_recovery"]+results["Low_h_recovery"]["h_lost_recovery"]+results["heat_recovery"]["h_lost_heat_recovery"]
        eva["e_lost_benchmark"]=results["B_primary_turbine"]["e_lost_turbine"]+results["B_secondary_turbine"]["e_lost_turbine"]
        #
        eva["lost_benchmark"]=eva["h_lost_benchmark"]+eva["e_lost_benchmark"]
        eva["power_benchmark"]=-results["B_primary_compressor"]["power_compressor"]-results["B_secondary_compressor"]["power_compressor"]+results["B_primary_turbine"]["power_turbine"]+results["B_secondary_turbine"]["power_turbine"]
        eva["hot_cost_benchmark"]=results["heat_recovery"]["hot_heat_recovery"]
        eva["hot_in_benchmark"]=-results["primary_h_exchanger"]["hot_out_h_exchanger"]-results["secondary_h_exchanger"]["hot_out_h_exchanger"]
        eva["Energy efficiency"]=eva["power_benchmark"]/(eva["hot_in_benchmark"]/0.95+eva["hot_cost_benchmark"])
        eva["mass_flow"]=H_in*0.95/eva["hot_in_benchmark"]
        eva["re_Heat_in"]=H_in
        eva["hot_cost"]=eva["mass_flow"]*eva["hot_cost_benchmark"]
        eva["power"]=eva["mass_flow"]*eva["power_benchmark"]
        eva["lost"]=(eva["mass_flow"]*eva["lost_benchmark"])+eva["hydrator_lost"]
        eva["flue_gas_mass_flow"] = eva["mass_flow"]*results["heat_recovery"]["n"]
        return eva


if __name__ == '__main__':
    
    parameters = dict() 
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["isentropic_eff_mc"] = 0.88  #等熵效率
    parameters["t_isentropic_eff_mc"] = 0.92
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["min_temperature_exchange"] = 15
    parameters["industrial_waste_heat_t"] =300 #℃
    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_reaction"] = 465
    parameters["p_bray_H"] = 20e6
    parameters["p_bray_M"] = 14e6
    parameters["p_bray_L"] = 7.6e6
    parameters["T_amb"] = 20
    parameters["p_amb"] = 101325

    BraytonHBs = Brayton(parameters)
    Hydrator_heat=1376678.4 #1.38MW
    results = BraytonHBs.solve(Hydrator_heat)
    print(results)