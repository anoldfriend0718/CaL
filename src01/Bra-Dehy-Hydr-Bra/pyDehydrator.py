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
from pyCostEstimator import Cost_Estimator
from pyBraytonHeatPump import BraytonHeatPump
from pyHENPinch import Hen_pinch_analyzer

#from pyMakeUpFlowBalaner import Make_Up_Flow_Balaner

M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol
M_H2O = 18e-3  # kg/mol

class Dehydrator(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
        self._bh = BraytonHeatPump(parameters)
        self._Store_electrical_power = parameters["Store_electrical_power"] #1MW
        self._res = self._bh.solve(self._Store_electrical_power)

        self._cao_conversion = parameters["cao_conversion"]#氢氧化钙转化率
        self._dehydrator_eff = parameters["dehydrator_eff"]#脱水器传热效率
        self._P_amb = parameters["P_amb"]#环境压力
        self._T_amb = parameters["T_amb"]#环境温度
        self._T_dehy = parameters["T_dehy"]#反应器温度
        self._delta_H_Tref = 104e3  # J/mole CaO504℃反应热
        self._BraytonHeatPump_cop = self._res["evaluation_indicators"]["cop"]
        self._H_in = self._Store_electrical_power*self._BraytonHeatPump_cop*self._dehydrator_eff
        self._steam_pressure_loss_ratio = parameters["steam_pressure_loss_ratio"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mcmechanical_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_dehydrator_distance = parameters["storage_dehydrator_distance"]
        self.deltaTmin_SSHX = parameters["deltaTmin_SSHX"] = 25   #固-固换热器最小温差
        self.deltaTmin_SGHX = parameters["deltaTmin_SGHX"] = 20   #固-气换热器最小温差

    def solve(self,input):
        results = {}
        results["BH"] = self._res
        #Basic input data
        initialvalue = self.initialvalue()
        results["initialvalue"] = initialvalue     
        #High temperature section of the heat exchange network
        ht_HEN = self.high_tem_HEN(self._res["heat_recovery"]["t_flue_gas_out"])
        results["ht_HEN"] = ht_HEN

        T_solid_in=ht_HEN
       # dehydrator(self, Ti_flue_gas, Ti_cao, Ti_water,To_water,Tcarb, pcarb, X):
        dehydrator = self.dehydrator(T_solid_in,
                                             self._T_dehy,
                                             self._P_amb,
                                             self._cao_conversion)
        results["dehydrator"] = dehydrator

        #Phase change part of the low temperature section of the heat exchange network
        m_steam = results["dehydrator"]["m_steam_out"]
        pc_HEN = self.pc_exchange(m_steam)
        results["pc_HEN"] = pc_HEN

        # steam blower
        steam_name = "water"
        steam_pi = self._P_amb
        steam_po = steam_pi/(1-self._steam_pressure_loss_ratio)
        steam_mass_rate = results["dehydrator"]["mass_steam_out"]
        steam_blower = self.steam_blower_power(self._T_dehy,
                                         steam_pi,
                                         steam_po,
                                         steam_mass_rate,
                                         steam_name)
        results["steam_blower"] = steam_blower
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["dehydrator"]["m_camix_in"],
            results["dehydrator"]["m_camix_out"])*(-1)
        hen_input = self.hen_input(results,
                                   input)
        he=Hen_pinch_analyzer(hen_input)
        he_text = he.write_pyPinch_data_text()
        hot_util, cold_util,total_HENA = he.solve(he_text)
        results["pinch_analysis_text"]=he_text
        results["total_HEN_area"]=total_HENA
        results["hot_utility"] = hot_util
        results["cold_utility"] = cold_util
        evaluation_indicators = self.evaluation_indicators(results,input)
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
        results["Store_electrical_power"] = self._Store_electrical_power
        results["BraytonHeatPump_cop"] = self._BraytonHeatPump_cop
        results["Q_dehydrator"] = self._H_in    
        return results
    
    def high_tem_HEN(self,t):
        cp_cao = self._pw.cp0mass_mean("cao", self._T_dehy, t-self.deltaTmin_SGHX+self.deltaTmin_SSHX)
        cp_caoh2 = self._pw.cp0mass_mean("caoh2", self._T_dehy-70, t-self.deltaTmin_SGHX)
        h_steam_in = CP.PropsSI('H', 'T', self._T_dehy+273.15, 'P', self._P_amb, "REFPROP::water")
        h_steam_out = CP.PropsSI('H', 'T', t+273.15, 'P', self._P_amb, "REFPROP::water")
        t_caoh2_out = t-self.deltaTmin_SGHX+((self._T_dehy-t+self.deltaTmin_SGHX-self.deltaTmin_SSHX)*cp_cao*M_cao+(h_steam_in-h_steam_out)*M_H2O)/(cp_caoh2*M_caoh2)     
        return t_caoh2_out
    
    def pc_exchange(self,m_h):
        h_hsteam_in =  CP.PropsSI('H', 'T', 105+273.15, 'P', self._P_amb, "REFPROP::water")
        h_hsteam_out =  CP.PropsSI('H', 'T', 95+273.15, 'P', self._P_amb, "REFPROP::water")
        h_lwater_in =  CP.PropsSI('H', 'T', 60+273.15, 'P', self._P_amb, "REFPROP::water")
        h_lwater_out =  CP.PropsSI('H', 'T', 85+273.15, 'P', self._P_amb, "REFPROP::water")
        m_l=m_h*(h_hsteam_in-h_hsteam_out)/(h_lwater_out-h_lwater_in)
        return m_l


    def dehydrator(self,T_solid_in, Tdehy, pdehy, X):
        
        delta_H_Tr = self._mole_reaction_heat(Tdehy, pdehy)
        heat_caoh2 = self._pw.cp0mass_mean("caoh2", Tdehy,T_solid_in)*(Tdehy-T_solid_in)*M_caoh2
        mole_Dehydration_reactions = self._H_in/(delta_H_Tr+heat_caoh2)
        mole_Dehydration_in = mole_Dehydration_reactions/X
        mole_caoh2_i = mole_Dehydration_in
        mass_caoh2_i = mole_caoh2_i*M_caoh2
        mole_cao_o = mole_Dehydration_reactions
        mass_cao_o = mole_cao_o*M_cao
        mole_caoh2_o = mole_Dehydration_in*(1-X)
        mass_caoh2_o = mole_caoh2_o*M_caoh2
        mass_camix_o=mass_caoh2_o+mass_cao_o
        mole_steam_o=mole_Dehydration_reactions
        mass_steam_o=mole_steam_o*M_H2O
        results ={}
        results["is_succeed"] = 1  
        results["delta_H_Tr"] = delta_H_Tr 
        results["heat_caoh2"] = heat_caoh2
        results["mole_dehydrator_reactions"]=mole_Dehydration_reactions
        results["mole_Dehydration_in"] = mole_Dehydration_in
        results["mass_caoh2_in"] = mass_caoh2_i
        results["m_camix_in"]= mass_caoh2_i
        results["mole_cao_o"] = mole_cao_o
        results["mass_cao_o"] = mass_cao_o
        results["mole_caoh2_o"] = mole_caoh2_o
        results["mass_caoh2_o"] = mass_caoh2_o
        results["mole_steam_out"] = mole_steam_o
        results["mass_steam_out"] = mass_steam_o
        results["m_steam_out"]= mass_steam_o
        results["m_camix_out"] = mass_camix_o
        return results
    
    def _mole_reaction_heat(self, Tdehy, pdehy):
        Tref = 504
        cp_cao_mean_Tref_Tr = self._pw.cp0mass_mean("cao", Tref, Tdehy)
        cp_caoh2_mean_Tref_Tr = self._pw.cp0mass_mean("caoh2", Tdehy, Tref)
        #热容：J/(Kg·K)
        delta_H_Tr = self._delta_H_Tref+((cp_caoh2_mean_Tref_Tr*M_caoh2)*(Tref-Tdehy)
                                         + (cp_cao_mean_Tref_Tr*M_cao)*(Tdehy-Tref))\
            + (CP.PropsSI('H', 'T', Tdehy+273.15, 'P', pdehy, "REFPROP::water") -
               CP.PropsSI('H', 'T', Tref+273.15, 'P', pdehy, "REFPROP::water"))*M_H2O
        return delta_H_Tr
    
    def steam_blower_power(self, Ti, pi, po, mass_rate, fluid):
        hi = CP.PropsSI('H', 'T', Ti+273.15, 'P', pi, fluid)
        si = CP.PropsSI('S', 'T', Ti+273.15, 'P', pi, fluid)
        ho_s = CP.PropsSI('H', 'P', po, 'S', si, fluid)
        ho_c = (ho_s-hi)/self._isentropic_eff_mc+hi
        To = CP.PropsSI('T', 'P', po, 'H', ho_c, fluid)-273.15
        W = ((ho_c-hi)/self._mechanical_eff)*mass_rate
        results = {}
        results["p_flue_gas_fan_out"] = po
        results["T_flue_gas_fan_out"] = To
        results["flue_gas_fan_power"] = W*(-1)
        return results
    
    def conveying_power(self, m_camix_in,m_camix_o):
        return self._convey_consumption*self._storage_dehydrator_distance * \
            (m_camix_in+m_camix_o)
    
    def hen_input(self,res,m1):
        inputs={}
        inputs["cao_conversion"]= 1
        inputs["m_camix_out"] = res["dehydrator"]["m_camix_out"]
        inputs["m_camix_in"] = res["dehydrator"]["m_camix_in"]
        inputs["m_steam_out"] = res["dehydrator"]["m_steam_out"]
        inputs["m_flue_gas_in"] = res["BH"]["evaluation_indicators"]["flue_gas_mass_flow"]
        inputs["m_water_in"] = m1
        inputs["T_dehy"] = self._T_dehy
        inputs["p_amb"] = self._P_amb
        inputs["T_amb"] = self._T_amb
        inputs["T_flue_gas_bray_out"] = res["BH"]["heat_recovery"]["t_flue_gas_out"]
        inputs["p_flue_gas_bray_out"] = 101325
        inputs["T_flue_gas_dew"] =160 #(不一定)
        inputs["T_solid_in"] = res["BH"]["heat_recovery"]["t_flue_gas_out"]-self.deltaTmin_SGHX#非定值-定值
        inputs["T_lsteam_in"] = 105
        inputs["T_lsteam_out"] = 95
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
    
    def evaluation_indicators(self,results,m2):
        a={}
        a["m_caoh2_in"] = results["dehydrator"]["mass_caoh2_in"]
        a["m_camix_out"] = results["dehydrator"]["m_camix_out"]
        a["m_steam_out"] = results["dehydrator"]["m_steam_out"]
        a["m_BHP_co2"] = results["BH"]["evaluation_indicators"]["mass_flow"]
        a["m_flue_gas"] = results["BH"]["evaluation_indicators"]["flue_gas_mass_flow"]
        a["m_heating_water"] = results["pc_HEN"]+m2
        a["power_in"] = -results["conveying_power"]-results["steam_blower"]["flue_gas_fan_power"]
        #a["hot_in"] = (self._Store_electrical_power*self._BraytonHeatPump_cop)
        a["hot_out"]= a["m_heating_water"]*(CP.PropsSI('H', 'T', 85+273.15, 'P', self._P_amb, "REFPROP::water")-CP.PropsSI('H', 'T', 60+273.15, 'P', self._P_amb, "REFPROP::water"))
        a["hot_stockpile"] = results["dehydrator"]["mole_dehydrator_reactions"]*68284.538#20℃下1mole反应热
        a["hot_lost"] = a["power_in"]+results["cold_utility"]+self._Store_electrical_power*self._BraytonHeatPump_cop*0.05
        a["energy_eff"] = (a["hot_out"]+a["hot_stockpile"])/(a["hot_out"]+a["hot_stockpile"]+a["hot_lost"])
        return a 


    
if __name__ == '__main__':
    parameters = dict()
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["cao_conversion"] = 1  #氧化钙转化率
    parameters["dehydrator_eff"] = 0.95   #脱水器效率
    parameters["T_amb"] = 20   #环境温度
    parameters["P_amb"] = 101325   #环境压力
    parameters["T_dehy"] = 525     #脱水器温度
    parameters["deltaTmin_SSHX"] = 25   #固-固换热器最小温差
    parameters["deltaTmin_SGHX"] = 20   #固-气换热器最小温差
    parameters["Store_electrical_power"] = 1e6
    parameters["steam_pressure_loss_ratio"] = 0.01
    parameters["isentropic_eff_mc"] = 0.88
    parameters["t_isentropic_eff_mc"] = 0.92
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["min_temperature_exchange"] = 15 
    parameters["industrial_waste_heat_t"] =300 #℃
    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_reaction"] = 525
    parameters["p_bray_H"] = 19447839.26865841
    parameters["p_bray_M"] = 12827110.4341202
    parameters["p_bray_L"] = 7.5e6
    parameters["p_amb"] = 101325

    parameters["mcmechanical_eff"] = 0.98
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100
    calcs = Dehydrator(parameters)
    m1=6
    results = calcs.solve(m1)
    print(results)
    print(results["pinch_analysis_text"])
