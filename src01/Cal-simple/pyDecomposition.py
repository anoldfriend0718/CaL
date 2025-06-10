import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")

import math
import numpy as np
import CoolProp.CoolProp as CP
from Cp0massWrapper import Cp0mass_Wrapper
from pyBHP import BraytonHeatPump
from pyHENPinch import Hen_pinch_analyzer

#from pyMakeUpFlowBalaner import Make_Up_Flow_Balaner

M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol
M_H2O = 18e-3  # kg/mol
M_caco3 = 100e-3
M_co2 = 44e-3
C_cao = 112396 #J/mole
C_caoh2 = 54593 
C_h2o = 750

class Dehydrator(object):
    def __init__(self, parameters,) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
        self._bh = BraytonHeatPump(parameters)
        self._Store_electrical_power=parameters["Store_electrical_power"]

        self._cao_conversion = parameters["cao_conversion"]#氢氧化钙转化率
        self._cao_purity = parameters["cao_purity"]
        self._dehydrator_eff = parameters["dehydrator_eff"]#脱水器传热效率
        self._P_amb = parameters["p_amb"]#环境压力
        self._T_amb = parameters["t_amb"]#环境温度

        self._delta_H_Tref = -178e3  # J/mole CaO

        self._steam_pressure_loss_ratio = parameters["steam_pressure_loss_ratio"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mechanical_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_dehydrator_distance = parameters["storage_dehydrator_distance"]
        self.deltaTmin_SSHX = parameters["deltaTmin_SSHX"] = 25   #固-固换热器最小温差
        self.deltaTmin_SGHX = parameters["deltaTmin_SGHX"] = 20   #固-气换热器最小温差

    def solve(self,input):
        results = {}
        self._res = self._bh.solve(input,505)
        results["BH"] = self._res
        self._p_dehy = input["p_Dehy"]
        self._T_dehy = 900#反应器温度
   
        #self._BraytonHeatPump_cop = self._res["evaluation_indicators"]["cop"]
        #self._H_in = self._Store_electrical_power*self._BraytonHeatPump_cop*self._dehydrator_eff
        #Basic input data
        initialvalue = self.initialvalue()
        results["initialvalue"] = initialvalue     
        #High temperature section of the heat exchange network
        results["ht_HEN"],loss=self.high_tem_HEN(self._cao_conversion)
        #results["ht_HEN"] = self._Dehy_caoh2_in
       # dehydrator(self, Ti_flue_gas, Ti_cao, Ti_water,To_water,Tcarb, pcarb, X):
        dehydrator = self.dehydrator(results["ht_HEN"],
                                             self._T_dehy,
                                             self._P_amb,
                                             self._cao_conversion,
                                             self._cao_purity)
        results["dehydrator"] = dehydrator
        #Phase change part of the low temperature section of the heat exchange network
        m_steam = results["dehydrator"]["out"]["m_steam"]
        #pc_HEN,pc_lost = self.pc_exchange(m_steam)
        #results["pc_HEN"] = pc_HEN
        #results["pc_lost"]=pc_lost

        # steam blower
        steam_name = "co2"
        steam_pi = self._p_dehy
        steam_po = self._P_amb/(1-self._steam_pressure_loss_ratio)
        steam_blower = self.steam_blower_power(self._T_dehy,
                                         steam_pi,
                                         steam_po,
                                         m_steam,
                                         steam_name)
        results["steam_blower"] = steam_blower
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["dehydrator"]["in"]["m_caco3"],
            results["dehydrator"]["out"]["m_camix"])*(-1)
        evaluation_indicators = self.evaluation_indicators(results)
        results["evaluation_indicators"] = evaluation_indicators
        return results 


    def initialvalue(self):
        results = {}
        results["cao_conversion"] = self._cao_conversion
        results["dehydrator_eff"] = self._dehydrator_eff
        results["P_amb"] = self._P_amb
        results["T_amb"] = self._T_amb
        results["T_dehy"]=self._T_dehy
        results["delta_H_Tref"] = self._delta_H_Tref
        return results
    
    def high_tem_HEN(self,X):
        cp_cao_o = self._pw.cp0mass_mean("cao", self._T_dehy, self._T_amb+self.deltaTmin_SSHX)
        cp_caco3_o = self._pw.cp0mass_mean("caco3", self._T_dehy, self._T_amb+self.deltaTmin_SSHX)
        cp_caco3_i = self._pw.cp0mass_mean("caco3", self._T_dehy-80, self._T_amb)
        h_steam_in = CP.PropsSI('H', 'T', self._T_dehy+273.15, 'P', self._P_amb, "REFPROP::co2")
        h_steam_out = CP.PropsSI('H', 'T', self._T_amb+self.deltaTmin_SGHX+273.15, 'P', self._P_amb, "REFPROP::co2")

        a=cp_caco3_i*M_caco3
        b=(cp_caco3_o*M_caco3*(1-X))+(cp_cao_o*M_cao*X)
        t_out = self._T_amb+((self._T_dehy-self._T_amb-self.deltaTmin_SSHX)*b+
                             (h_steam_in-h_steam_out)*M_co2*X)/(a/0.96)
        a1=self._T_dehy-self.deltaTmin_SSHX
        if t_out> a1:
            t_out=a1
        else:
            t_out=t_out
        h_mole_lost=((self._T_dehy-self._T_amb-self.deltaTmin_SSHX)*b+
                             (h_steam_in-h_steam_out)*M_co2*X)-(t_out-self._T_amb)*a
        return t_out,h_mole_lost
    

    def dehydrator(self,T_solid_in, Tdehy, pdehy, X , Y):
        
        delta_H_Tr = -self._mole_reaction_heat(Tdehy, pdehy)
        heat_caoh2 = self._pw.cp0mass_mean("caco3", Tdehy, T_solid_in)*M_caco3*(Tdehy-T_solid_in)

        power = (delta_H_Tr*X+heat_caoh2)/0.97

        results ={}
        results["is_succeed"] = 1  
        results["delta_H_Tr"] = delta_H_Tr 
        results["heat_caco3"] = heat_caoh2
        results["mole_dehydrator_reactions"]=X
        results["mole_Dehydration_in"] = 1
        results["in"]={}
        results["in"]["m_caco3"] = results["mole_Dehydration_in"]*M_caco3
        results["in"]["mole_caco3"] = results["mole_Dehydration_in"]
        results["out"]={}
        results["out"]["mole_cao"]=results["mole_Dehydration_in"]*X
        results["out"]["m_cao"]=results["mole_Dehydration_in"]*X*M_cao
        results["out"]["m_caco3"] = results["mole_Dehydration_in"]*(1-X)*M_caco3
        results["out"]["mole_caco3"] = results["mole_Dehydration_in"]*(1-X)
        results["out"]["m_camix"] =results["out"]["m_cao"]+results["out"]["m_caco3"]
        results["out"]["m_steam"]= results["mole_Dehydration_in"]*X*M_co2
        results["out"]["mole_steam"] = results["mole_Dehydration_in"]*X

        results["exergy"] = {}
        results["exergy"]["power"]=power
        results["exergy"]["chemical_t"]=self._bh.ex_calculations1( Tdehy,delta_H_Tr*X)
        results["exergy"]["sensible_heat"] = self._bh.ex_calculations(T_solid_in, Tdehy,heat_caoh2)
        results["exergy"]["lost"] = (results["exergy"]["power"]-results["exergy"]["chemical_t"]-
                                     results["exergy"]["sensible_heat"])
    
        return results
    
    def _mole_reaction_heat(self, Tdehy, pdehy):
        Tref = 20
        cp_cao_mean_Tref_Tr = self._pw.cp0mass_mean("cao", Tref, Tdehy)
        cp_caco3_mean_Tref_Tr = self._pw.cp0mass_mean("caco3", Tref, Tdehy)

        delta_H_Tr = self._delta_H_Tref+(cp_caco3_mean_Tref_Tr*M_caco3
                                         - cp_cao_mean_Tref_Tr*M_cao)*(Tdehy-Tref)\
            - (CP.PropsSI('H', 'T', Tdehy+273.15, 'P', pdehy, "REFPROP::co2") -
               CP.PropsSI('H', 'T', Tref+273.15, 'P', pdehy, "REFPROP::co2"))*M_co2
               
        return delta_H_Tr
    
    def steam_blower_power(self, Ti, pi, po, mass_rate, fluid):
        hi = CP.PropsSI('H', 'T', Ti+273.15, 'P', pi, fluid)
        si = CP.PropsSI('S', 'T', Ti+273.15, 'P', pi, fluid)
        ho_s = CP.PropsSI('H', 'P', po, 'S', si, fluid)
        ho_c = (ho_s-hi)/self._isentropic_eff_mc+hi
        To = CP.PropsSI('T', 'P', po, 'H', ho_c, fluid)-273.15
        W = ((ho_c-hi)/self._mechanical_eff)*mass_rate
        results = {}
        results["p_out"] = po
        results["T_out"] = To
        results["power"] = W*(-1)
        return results
    
    def conveying_power(self, m_camix_in,m_camix_o):
        return self._convey_consumption*self._storage_dehydrator_distance * \
            (m_camix_in+m_camix_o)
    
    def hen_input(self,res,m1):
        inputs={}
        inputs["cao_conversion"]= self._cao_conversion
        inputs["cao_purity"] =self._cao_purity #氢氧化钙含量
        inputs["m_camix_out"] = res["dehydrator"]["out"]["m_camix"]*0.96
        inputs["m_camix_in"] = res["dehydrator"]["in"]["m_camix"]
        inputs["m_steam_out"] = res["dehydrator"]["out"]["m_steam"]*0.96
        inputs["m_flue_gas_in"] = res["BH"]["evaluation_indicators"]["flue_gas_mass_flow"]*0.96
        inputs["m_water_in"] = m1
        inputs["T_dehy"] = self._T_dehy
        inputs["T_steam_in"]= res["steam_blower"]["T_out"]
        inputs["p_amb"] = self._P_amb
        inputs["T_amb"] = self._T_amb
        inputs["T_flue_gas_bray_out"] = res["BH"]["heat_recovery"]["t_flue_gas_out"]
        inputs["p_flue_gas_bray_out"] = 101325
        inputs["T_flue_gas_dew"] =160 #(不一定)
        inputs["T_solid_in"] = self._Dehy_caoh2_in
        
        inputs["T_lsteam_in"] = 102
        inputs["T_lsteam_out"] = 98
        inputs["T_water_supply_in"] = 60
        inputs["T_water_reactor_in"] = 85
        inputs["T_delta_pinch"] = 20
        inputs["p_dehy_o"] = 101325
        inputs["p_dehy_i"] = 101325
        inputs["p_water_after_pump"] = 101325
        flue_gas_composistion = dict()
        flue_gas_composistion["co2"] = 0.1338
        flue_gas_composistion["o2"] = 0.0384
        flue_gas_composistion["n2"] = 0.6975
        inputs["flue_gas_composition"] = flue_gas_composistion
        return inputs
    
    def evaluation_indicators(self,results):
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        fluid=flue_gas_name
        
        a={}
        a["power_in"] = (-results["conveying_power"]-results["steam_blower"]["power"]+
                         results["dehydrator"]["exergy"]["power"])
        a["hot_stockpile"] =-self._delta_H_Tref*self._cao_conversion
        a["hot_lost"] = a["power_in"]-a["hot_stockpile"]
        a["energy_eff"] = (a["hot_stockpile"])/(a["power_in"])

        a["exergy_in"]=(a["power_in"])
        a["exergy_out"] = results["dehydrator"]["exergy"]["chemical_t"]
        a["exergy_eff"] = a["exergy_out"]/a["exergy_in"]
        return a 


    
if __name__ == '__main__':
    parameters = dict()
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["isentropic_eff_mc"] = 0.88
    parameters["t_isentropic_eff_mc"] = 0.92
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["min_temperature_exchange"] = 15 
    parameters["deltaTmin_SSHX"] = parameters["min_temperature_exchange"]+5   #固-固换热器最小温差
    parameters["deltaTmin_SGHX"] = parameters["min_temperature_exchange"]   #固-气换热器最小温差
    parameters["industrial_waste_heat_t"] =300 #℃
    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_amb"] = 20   #环境温度
    parameters["p_amb"] = 101325   #环境压力

    parameters["p_bray_L"] = 7.5e6
    parameters["Store_electrical_power"] = 1e6

    parameters["cao_conversion"] = 0.4  #氧化钙转化率
    parameters["cao_purity"] = 1 #氢氧化钙含量
    parameters["dehydrator_eff"] = 0.97   #脱水器效率
    parameters["steam_pressure_loss_ratio"] = 0.01
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100

    calcs = Dehydrator(parameters) 

    inputs={}
    inputs["p_bray_H"] = 19447839.26865841#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12827110.4341202 #优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 1e5 #变量4，反应器压力
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度
    results = calcs.solve(inputs)
    print(results)

