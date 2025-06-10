import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")
import json
import math
import numpy as np
import CoolProp.CoolProp as CP
from Cp0massWrapper import Cp0mass_Wrapper
from pyBraytonHeatPump import BraytonHeatPump
from pyHENPinch import Hen_pinch_analyzer

#from pyMakeUpFlowBalaner import Make_Up_Flow_Balaner

M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol
M_H2O = 18e-3  # kg/mol
M_caco3 = 100e-3
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
        self._heat_transfer_loss_eff=parameters["heat_transfer_loss_eff"]#换热损失
        self._P_amb = parameters["p_amb"]#环境压力
        self._T_amb = parameters["t_amb"]#环境温度
        self._T_flue_gas_i=parameters["industrial_waste_heat_t"]
        self._T_water_supply_in= 60
        self._p_water_supply_in = parameters["p_water_supply_in"]  
        self._T_water_prod_out = 85
        self._water_pressure_drop_rate=parameters["water_pressure_drop_rate"]
        self._water_pipe_length=parameters["water_pipe_length"]
        self._pump_hydraulic_eff=parameters["water_pump_hydraulic_efficiency"]
        self._pump_mechanical_eff=parameters["water_pump_mechanical_efficiency"]

        self._delta_H_Tref = 104e3  # J/mole CaO504℃反应热

        self._steam_pressure_loss_ratio = parameters["steam_pressure_loss_ratio"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mechanical_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_dehydrator_distance = parameters["storage_dehydrator_distance"]
        

    def solve(self,input):
        self.deltaTmin_SSHX = input["min_temperature_exchange"]+5  #固-固换热器最小温差
        self.deltaTmin_SGHX = input["min_temperature_exchange"]  #固-气换热器最小温差
        self.deltaTmin_HEN  = input["min_temperature_HEN"]
        results = {}
        self._res = self._bh.solve(input)
        results["BH"] = self._res
        self._p_dehy = input["p_Dehy"]
        self._T_dehy = self._res["initialvalue"]["t_reaction"]#反应器温度
        self._BraytonHeatPump_cop = self._res["evaluation_indicators"]["cop"]
        self._H_in = self._Store_electrical_power*self._BraytonHeatPump_cop*self._dehydrator_eff
        #Basic input data
        initialvalue = self.initialvalue()
        results["initialvalue"] = initialvalue     
        #High temperature section of the heat exchange network
        self._Dehy_caoh2_in,a,results=self.caoh2_in(results,1e-2)
        results["ht_HEN"] = self._Dehy_caoh2_in
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["dehydrator"]["in"]["m_camix"],
            results["dehydrator"]["out"]["m_camix"])*(-1)
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        results["flue_gas_fan_power"]=self.flue_gas_fan_power(self._T_flue_gas_i,
                                                            self._P_amb,
                                                            self._P_amb/0.99,
                                                            results["BH"]["evaluation_indicators"]["flue_gas_mass_flow"],
                                                            flue_gas_name)
        p_water_after_pump=self._p_water_supply_in+self._water_pipe_length*self._water_pressure_drop_rate
        results["water_pump_power"]=self._water_pump_power(self._p_water_supply_in,
            p_water_after_pump,a)
        evaluation_indicators,Caes = self.evaluation_indicators(results,a)
        results["evaluation_indicators"] = evaluation_indicators
        results["Caes"] = Caes
        results["cost"] = self.cost(results)
        return results 
    
    def caoh2_in(self,results,tolerance):
        # 定义搜索范围的上下界
        low = 350
        high = self._T_dehy
        mid = None
        a=12*self._Store_electrical_power/1e6
        
        while high - low > tolerance:
            mid = (low + high) / 2
            Dehy_caoh2_in = mid
            dehydrator = self.dehydrator(Dehy_caoh2_in,
                                             self._T_dehy,
                                             self._P_amb,
                                             self._cao_conversion,
                                             self._cao_purity)
            results["dehydrator"] = dehydrator
            steam_name = "water"
            steam_pi = self._p_dehy
            steam_po = self._P_amb/(1-self._steam_pressure_loss_ratio)
            steam_mass_rate = results["dehydrator"]["out"]["m_steam"]
            steam_blower = self.steam_blower_power(self._T_dehy,
                                            steam_pi,
                                            steam_po,
                                            steam_mass_rate,
                                            steam_name)
            results["steam_blower"] = steam_blower
            hen_input = self.hen_input(results,
                                    a,
                                    mid)
            he=Hen_pinch_analyzer(hen_input)
            he_text = he.write_pyPinch_data_text()
            hot_util, cold_util,total_HENA = he.solve(he_text)
            # 根据hot_util的值调整搜索范围
            if hot_util > 0:
                high = mid
            else:
                low = mid
        low = round(low, 3)
        step_size=1
        max_a=10000
        while a <= max_a:
            hen_input = self.hen_input(results,
                                    a,
                                    low)
            he=Hen_pinch_analyzer(hen_input)
            
            he_text = he.write_pyPinch_data_text()
            hot_util, cold_util,total_HENA = he.solve(he_text) 

            if hot_util <= 0 and step_size >= 0.1:
                a = a+step_size
            elif hot_util >= 0 and step_size >= 0.1:
                a = a-step_size*0.9
                step_size = step_size / 10
            else:
                results["HEN_exergy"],results["HEN_energy"]=he.exergy()
                results["pinch_analysis_text"]=he_text
                results["total_HEN_area"]=total_HENA/self._heat_transfer_loss_eff
                results["hot_utility"] = hot_util
                results["cold_utility"] = cold_util
                break 
        a = round(a, 3)
        return low, a,results

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
    
    def high_tem_HEN(self,t,X,Y):
        cp_cao_o = self._pw.cp0mass_mean("cao", self._T_dehy, t-self.deltaTmin_SGHX+self.deltaTmin_SSHX)
        cp_caoh2_o = self._pw.cp0mass_mean("caoh2", self._T_dehy, t-self.deltaTmin_SGHX+self.deltaTmin_SSHX)
        cp_caco3_o = self._pw.cp0mass_mean("caco3", self._T_dehy, t-self.deltaTmin_SGHX+self.deltaTmin_SSHX)
        cp_caoh2_i = self._pw.cp0mass_mean("caoh2", self._T_dehy-80, t-self.deltaTmin_SGHX)
        cp_caco3_i = self._pw.cp0mass_mean("caco3", self._T_dehy-80, t-self.deltaTmin_SGHX)
        h_steam_in = CP.PropsSI('H', 'T', self._T_dehy+273.15, 'P', self._P_amb, "REFPROP::water")
        h_steam_out = CP.PropsSI('H', 'T', t+273.15, 'P', self._P_amb, "REFPROP::water")
        a=(cp_caoh2_i*M_caoh2/X)+(cp_caco3_i*((M_caoh2/X/Y)-M_caoh2/X))
        b=(cp_caoh2_o*(M_caoh2/X-M_caoh2))+(cp_cao_o*M_cao)+(cp_caco3_o*((M_caoh2/X/Y)-M_caoh2/X))
        t_caoh2_out = t-self.deltaTmin_SGHX+((self._T_dehy-t+self.deltaTmin_SGHX-self.deltaTmin_SSHX)*b+(h_steam_in-h_steam_out)*M_H2O)/(a/0.96)
        h_mole_lost=((t_caoh2_out-t+self.deltaTmin_SGHX)*a/0.96) *0.04
        return t_caoh2_out,h_mole_lost
    
    def pc_exchange(self,m_h):
        h_hsteam_in =  CP.PropsSI('H', 'T', 105+273.15, 'P', self._P_amb, "REFPROP::water")
        h_hsteam_out =  CP.PropsSI('H', 'T', 95+273.15, 'P', self._P_amb, "REFPROP::water")
        h_lwater_in =  CP.PropsSI('H', 'T', 60+273.15, 'P', self._P_amb, "REFPROP::water")
        h_lwater_out =  CP.PropsSI('H', 'T', 85+273.15, 'P', self._P_amb, "REFPROP::water")
        
        m_l=m_h*(h_hsteam_in-h_hsteam_out)*0.96/(h_lwater_out-h_lwater_in)
        pc_lost = m_h*(h_hsteam_in-h_hsteam_out)*0.04

        return m_l,pc_lost

    def dehydrator(self,T_solid_in, Tdehy, pdehy, X , Y):
        
        delta_H_Tr = self._mole_reaction_heat(Tdehy, pdehy)
        heat_caoh2 = (self._pw.cp_camix_mean_Ci(T_solid_in,Tdehy,Y)*(Tdehy-T_solid_in)*M_caoh2)/X

        mole_Dehydration_reactions = self._H_in/(delta_H_Tr+heat_caoh2)
        mole_Dehydration_in = mole_Dehydration_reactions/X
        mole_caoh2_i = mole_Dehydration_in
        mass_caoh2_i = mole_caoh2_i*M_caoh2
        mass_in = mass_caoh2_i/Y
        mass_caco3_i = mass_in-mass_caoh2_i
        mole_caco3_i = mass_caco3_i/M_caco3

        mole_cao_o = mole_Dehydration_reactions
        mass_cao_o = mole_cao_o*M_cao
        mole_caoh2_o = mole_Dehydration_in*(1-X)
        mass_caoh2_o = mole_caoh2_o*M_caoh2
        mass_camix_o=mass_caoh2_o+mass_cao_o+mass_caco3_i

        mole_steam_o=mole_Dehydration_reactions
        mass_steam_o=mole_steam_o*M_H2O
        results ={}
        results["is_succeed"] = 1  
        results["delta_H_Tr"] = delta_H_Tr 
        results["heat_caoh2"] = heat_caoh2
        results["mole_dehydrator_reactions"]=mole_Dehydration_reactions
        results["mole_Dehydration_in"] = mole_Dehydration_in
        results["in"]={}
        results["in"]["m_camix"] = mass_in
        results["in"]["m_caoh2"] = mass_caoh2_i
        results["in"]["mole_caoh2"] = mole_caoh2_i
        results["in"]["m_caco3"] = mass_caco3_i
        results["in"]["mole_caco3"] = mole_caco3_i
        results["out"]={}
        results["out"]["m_camix"] = mass_camix_o
        results["out"]["m_cao"]=mass_cao_o
        results["out"]["mole_cao"]=mole_cao_o
        results["out"]["m_caoh2"]=mass_caoh2_o
        results["out"]["mole_caoh2"]=mole_caoh2_o
        results["out"]["m_caco3"] = mass_caco3_i
        results["out"]["mole_caco3"] = mole_caco3_i
        results["out"]["m_steam"]= mass_steam_o
        results["out"]["mole_steam"] = mole_cao_o
        results["exergy"] = {}
        results["exergy"]["hot_in"]=self._res["evaluation_indicators"]["exergy_out"]
        results["exergy"]["chemical_t"]=self._bh.ex_calculations1( Tdehy,delta_H_Tr*
                                                                      mole_Dehydration_reactions)
        results["exergy"]["chemical"] = (112.396-54.593)*1000*mole_Dehydration_reactions
        results["exergy"]["sensible_heat"] = self._bh.ex_calculations(T_solid_in, Tdehy,heat_caoh2*
                                                                      mole_Dehydration_reactions)
        results["exergy"]["lost"] = (results["exergy"]["hot_in"]-results["exergy"]["chemical_t"]-
                                     results["exergy"]["sensible_heat"])
        return results
    
    def _mole_reaction_heat(self, Tdehy, pdehy):
        Tref = 25
        cp_cao_mean_Tref_Tr = self._pw.cp0mass_mean("cao", Tref, Tdehy)
        cp_caoh2_mean_Tref_Tr = self._pw.cp0mass_mean("caoh2", Tdehy, Tref)
        #热容：J/(Kg·K)
        delta_H_Tr = self._delta_H_Tref+((cp_caoh2_mean_Tref_Tr*M_caoh2)*(Tref-Tdehy)
                                         + (cp_cao_mean_Tref_Tr*M_cao)*(Tdehy-Tref))\
            + (CP.PropsSI('H', 'T', Tdehy+273.15, 'P', pdehy, "REFPROP::water") -
               CP.PropsSI('H', 'T', 105+273.15, 'P', pdehy, "REFPROP::water") + 
               CP.PropsSI('H', 'T', 95+273.15, 'P', pdehy, "REFPROP::water") - 
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
        results["p_out"] = po
        results["T_out"] = To
        results["power"] = W*(-1)
        return results
    def flue_gas_fan_power(self, Ti, pi, po, mass_rate, fluid):
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
        return W*(-1)
    def _water_pump_power(self,pi,po,mass_rate):
        rho_water=1000
        vol_rate=mass_rate/rho_water
        dp=po-pi
        power=vol_rate*dp/self._pump_hydraulic_eff/self._pump_mechanical_eff*(-1)
        return power
    
    def conveying_power(self, m_camix_in,m_camix_o):
        return self._convey_consumption*self._storage_dehydrator_distance * \
            (m_camix_in+m_camix_o)
    
    def hen_input(self,res,m1,mid):
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
        inputs["T_solid_in"] = mid
        
        inputs["T_lsteam_in"] = 102
        inputs["T_lsteam_out"] = 98
        inputs["T_water_supply_in"] = 60
        inputs["T_water_reactor_in"] = 85
        inputs["T_delta_pinch"] = self.deltaTmin_HEN 
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
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        fluid=flue_gas_name
        t_flue_gas_in = results["BH"]["heat_recovery"]["t_flue_gas_out"]
        t_flue_gas_out = 160
        h_flue_gas_in = CP.PropsSI('H', 'T', t_flue_gas_in+273.15, 'P', self._P_amb , fluid)
        h_flue_gas_out = CP.PropsSI('H', 'T', t_flue_gas_out+273.15, 'P', self._P_amb , fluid)
        h_in = (h_flue_gas_in-h_flue_gas_out)*results["BH"]["evaluation_indicators"]["flue_gas_mass_flow"]

        delta_H_Tr = self._delta_H_Tref - (
               CP.PropsSI('H', 'T', 102+273.15, 'P', self._P_amb, "REFPROP::water") - 
               CP.PropsSI('H', 'T', 98+273.15, 'P', self._P_amb, "REFPROP::water"))*M_H2O

        a={}
        a["m_heating_water"] =m2
        a["power_in"] = -results["conveying_power"]-results["steam_blower"]["power"]-results["water_pump_power"]-results["flue_gas_fan_power"]
        a["hot_in_HP"] = self._Store_electrical_power*self._BraytonHeatPump_cop
        a["hot_in_gas"] = h_in
        a["hot_in"] = a["hot_in_HP"]+a["hot_in_gas"]
        a["hot_out"]= a["m_heating_water"]*(CP.PropsSI('H', 'T', 85+273.15, 'P', self._P_amb, "REFPROP::water")
                                            -CP.PropsSI('H', 'T', 60+273.15, 'P', self._P_amb, "REFPROP::water"))
        #a["hot_stockpile"] = results["dehydrator"]["mole_dehydrator_reactions"]*68284.538#20℃下1mole反应热
        
        a["hot_stockpile"] = delta_H_Tr*results["dehydrator"]["mole_Dehydration_in"]
        a["hot_lost"] = a["power_in"]+a["hot_in"]-a["hot_out"]-a["hot_stockpile"]
        a["hot_HEN_lost"] = a["hot_lost"] - (a["power_in"]+a["hot_in_HP"]*0.05)
        a["energy_eff"] = (a["hot_out"]+a["hot_stockpile"])/(a["power_in"]+a["hot_in"])

        a["exergy_in"]=(results["BH"]["evaluation_indicators"]["exergy_out"]
                        +results["HEN_exergy"]["C_flue_gas"]/0.96+a["power_in"])
        gcpw2=(25-293.15*math.log(358.15/333.15))/(25)
        a["exergy_out"] = results["HEN_exergy"]["C_water"]+results["dehydrator"]["exergy"]["chemical"]
        a["exergy_eff"] = a["exergy_out"]/a["exergy_in"]

        b={}
        b["power_in"]=self._Store_electrical_power+a["power_in"]
        b["hot_in"]=results["BH"]["evaluation_indicators"]["hot_cost"]+a["hot_in_gas"]
        b["hot_out"]=a["hot_out"]
        b["hot_stockpile"]=a["hot_stockpile"]
        b["h_lost"]=results["BH"]["evaluation_indicators"]["lost"]+a["hot_lost"]
        b["energy_eff"]=(b["hot_out"]+b["hot_stockpile"])/(b["power_in"]+b["hot_in"])
        b["exergy_eff"]=a["exergy_out"]/(results["BH"]["evaluation_indicators"]["exergy_in"]
                                         +results["HEN_exergy"]["C_flue_gas"]/0.96+a["power_in"])
        b["Cost"]={}
        b["Cost"]["rea"] = 533394*(a["hot_in_HP"]/293000)**0.48
        b["Cost"]["h_ex"] = 3197*results["total_HEN_area"]**0.67*1.01325**0.28
        return a ,b
    def cost(self,results):
        result={}
        result["cost_hen"]=3197*(results["total_HEN_area"]**0.67)
        result["cost_hen2"]=2957.636*(results["total_HEN_area"]**0.67)*10**0.28
        result["cost_dehy"]=19594*(self._Store_electrical_power*self._BraytonHeatPump_cop/1000)**0.5
        result["cost_steam_blow"]=750*((-results["steam_blower"]["power"]/1000)**0.71)*(1+0.2/(1-self._isentropic_eff_mc))
        result["cost_flue_gas_fan"] = 0.827*1000000/6.901*(-results["flue_gas_fan_power"]/445000)**0.67
        result["cost_water_pump"]= 0.103*1000000/6.901*(-results["water_pump_power"]/4000)**0.55
        
        result["cost_all"]=result["cost_hen"]+result["cost_dehy"]+result["cost_steam_blow"]+result["cost_flue_gas_fan"]+result["cost_water_pump"]
        return result

    
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
    parameters["industrial_waste_heat_t"] =300 #℃
    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_amb"] = 20   #环境温度
    parameters["p_amb"] = 101325   #环境压力

    parameters["p_bray_L"] = 7.5e6
    parameters["Store_electrical_power"] = 50e6
    parameters["p_water_supply_in"] = 2e5 
    parameters["water_pressure_drop_rate"] = 100 #100Pa/m
    parameters["water_pipe_length"] = 1000
    parameters["water_pump_hydraulic_efficiency"] = 0.75
    parameters["water_pump_mechanical_efficiency"] = 0.94

    parameters["cao_conversion"] = 0.95  #氧化钙转化率
    parameters["cao_purity"] = 0.98 #氢氧化钙含量
    parameters["dehydrator_eff"] = 0.95   #脱水器效率
    parameters["steam_pressure_loss_ratio"] = 0.01
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100

    calcs = Dehydrator(parameters) 

    inputs={}
    inputs["p_bray_H"] = 18716228#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12481150 #优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 1e5 #变量4，反应器压力
    inputs["Economic Model Selection"] = 2 #经济模型选择，1Tesio，2Nathan T
    inputs["Compressor power limit"] = 200e6#功率界限，影响齿轮离心和滚筒离心模型的选取，单位W，桶式离心需体积流量
    inputs["Turbine power limit"] = 35e6
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度
    inputs["min_temperature_exchange"] = 15 #变量1，换热器最小换热温差
    inputs["min_temperature_HEN"] = 15 #变量1，HEN最小换热温差
    results = calcs.solve(inputs)
    data_for_json = {key: (value.item() if isinstance(value, np.floating) else value) for key, value in results.items()}
    json_output = json.dumps(data_for_json, indent=4)
    print(json_output)


