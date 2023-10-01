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

class Braytonyz(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
        self._flue_gas_composition = self._pw._norm_flue_gas_composition
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]#压缩机等熵效率

        self._t_isentropic_eff_mc = parameters["t_isentropic_eff_mc"]#透平等熵效率
        self._mechanical_eff = parameters["mechanical_eff"]#机械效率
        self._min_temperature_exchange = parameters["min_temperature_exchange"]#最小换热温差
        self._industrial_waste_heat_t = parameters["industrial_waste_heat_t"]#工业余热温度
        self._heat_transfer_loss_eff = parameters["heat_transfer_loss_eff"]#换热损失
        self._t_reaction = parameters["t_reaction"]
        #self._p_bray_H = parameters["p_bray_H"]
        #self._p_bray_M = parameters["p_bray_M"]
        self._p_bray_L = parameters["p_bray_L"]
        self._p_amb = parameters["p_amb"]
        self._T_amb = parameters["T_amb"]
        self._T_L = 31.55
        
    def solve(self,inputs):
        results = {}
        self._p_bray_H = inputs["p_bray_H"]
        self._p_bray_M = inputs["p_bray_M"]
        #Basic input data
        initialvalue = self.initialvalue()
        results["B_initialvalue"] = initialvalue
        #The primary Turbine is the starting point
        t1=236.45
        primary_turbine = self.turbine(t1,
                                                (self._p_bray_H*0.99),
                                                (self._p_bray_L/0.99)/0.99)
        results["B_primary_turbine"] = primary_turbine

        #primary compressors
        primary_compressor = self.compressor(self._T_L,
                                             self._p_bray_L,
                                             self._p_bray_H)
        results["B_primary_compressor"] = primary_compressor
        self._t1 = results["B_primary_compressor"]["t_compressor_out"]
        #High heat_recovery 
        t_low = 45.65
        t_h_re_in=results["B_primary_turbine"]["t_turbine_out"]
        High_h_recovery=self.h_recovery(t_h_re_in,
                                                     t_low,
                                                     self._p_bray_L,
                                                     self._p_bray_H)
        results["High_h_recovery"] = High_h_recovery


        t_sec_h_in=results["High_h_recovery"]["t_h_recovery_ll_out"]
        secondary_h_exchanger=self.h_exchanger(t_sec_h_in,
                                                     t1,
                                                     self._p_bray_H*0.99)
        results["secondary_h_exchanger"] = secondary_h_exchanger

        #Cooling towers
        cooling_tower = self.cooling_tower(t_low,
                                           self._T_L,
                                           self._p_bray_L)
        results["cooling_tower"] = cooling_tower

        #Data synthesi
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
        t_h_exchangerm_hh_in=T_in
        t_h_exchangerm_hh_out=T_out

        t_h_exchangerm_ll_in=self._t1
        #t_h_exchangerm_ll_out
        h_h_exchangerm_hh_in =CP.PropsSI('H', 'T', t_h_exchangerm_hh_in+273.15,  'P', ((Pl/0.99)/0.99), "REFPROP::co2")
        h_h_exchangerm_hh_out=CP.PropsSI('H', 'T', t_h_exchangerm_hh_out+273.15, 'P', (Pl/0.99), "REFPROP::co2")
        s_h_exchangerm_hh_in =CP.PropsSI('S', 'T', t_h_exchangerm_hh_in+273.15,  'P', Pl, "REFPROP::co2")
        hot_hh_exchange=h_h_exchangerm_hh_in-h_h_exchangerm_hh_out
        hot_ll_exchange= hot_hh_exchange*self._heat_transfer_loss_eff
        h_h_exchangerm_ll_in =CP.PropsSI('H', 'T', t_h_exchangerm_ll_in+273.15,  'P', Ph, "REFPROP::co2")
        h_h_exchangerm_ll_out=h_h_exchangerm_ll_in+hot_ll_exchange
        t_h_exchangerm_ll_out=CP.PropsSI('T', 'H', h_h_exchangerm_ll_out,  'P', (Ph*0.99), "REFPROP::co2")-273.15
        s_h_exchangerm_ll_out=CP.PropsSI('S', 'T', t_h_exchangerm_ll_out+273.15, 'P', Ph, "REFPROP::co2")
        results = {}
        results["t_h_recovery_hh_out"] = t_h_exchangerm_hh_out 
        results["p_h_recovery_hh_out"] = self._p_bray_L/0.99
        results["h_h_recovery_hh_out"] = h_h_exchangerm_hh_out

        results["t_h_recovery_ll_out"] = t_h_exchangerm_ll_out 
        results["p_h_recovery_ll_out"] = self._p_bray_H*0.99
        results["h_h_recovery_ll_out"] = h_h_exchangerm_ll_out
        results["h_lost_recovery"] = hot_hh_exchange-hot_ll_exchange
        return results
    
    def cooling_tower(self,T_in,T_out,P):
        results = {}
        results["t_cooling_tower_in"] = T_in
        results["p_cooling_tower_in"] = P
        results["h_cooling_tower_in"] = CP.PropsSI('H', 'T', T_in+273.15, 'P', P/0.99, "REFPROP::co2")
        results["s_cooling_tower_in"] = CP.PropsSI('S', 'T', T_in+273.15, 'P', P/0.99, "REFPROP::co2")
        results["t_cooling_tower_out"] = T_out
        results["p_cooling_tower_out"] = P
        results["h_cooling_tower_out"] = CP.PropsSI('H', 'T', T_out+273.15, 'P', P, "REFPROP::co2")
        results["s_cooling_tower_out"] = CP.PropsSI('S', 'T', T_out+273.15, 'P', P, "REFPROP::co2")
        results["hot_cooling_tower"] = results["h_cooling_tower_in"]-results["h_cooling_tower_out"]
        return results



if __name__ == '__main__':
    
    parameters = dict() 
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["isentropic_eff_mc"] = 0.8  #压缩机等熵效率

    parameters["t_isentropic_eff_mc"] = 0.8#透平等熵效率
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["min_temperature_exchange"] = 15
    parameters["industrial_waste_heat_t"] =300 #℃
    parameters["heat_transfer_loss_eff"] = 1
    parameters["t_reaction"] = 465
    parameters["p_bray_H"] = 25e6
    parameters["p_bray_M"] = 14e6
    parameters["p_bray_L"] = 7.622e6
    parameters["T_amb"] = 20
    parameters["p_amb"] = 101325

    BraytonHBs = Braytonyz(parameters)

    inputs={}
    inputs["Hydrator_heat"] = 1609467.621084121
    inputs["p_bray_H"] = 9.432e6
    inputs["p_bray_M"] = 7.622e6
    results = BraytonHBs.solve(inputs)
    print(results)