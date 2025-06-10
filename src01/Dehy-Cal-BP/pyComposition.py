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
from pyDecomposition import Dehydrator
from pyPinch import PyPinch
from pyBrayton import Brayton
import pandas as pd
M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol
M_H2O = 18e-3  # kg/mol
M_co2 = 44e-3
M_caco3 = 100e-3
M_1= (44e-3*0.1336+32e-3*0.0384+28e-3*0.6975)/(0.1336+0.0384+0.6975)
M_2= (44e-3*0.01336+32e-3*0.0384+28e-3*0.6975)/(0.01336+0.0384+0.6975)
C_cao = 112396 #J/mole
C_caoh2 = 54593 
C_h2o = 750
C_co2 = 20033
C_caco3=1528

class Hydrator(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"],
                                   parameters["flue_gas_composition2"])
        self._de = Dehydrator(parameters)
        self._B = Brayton(parameters)

        self._cao_conversion = parameters["cao_conversion2"]#氢氧化钙转化率
        self._cao_purity = parameters["cao_purity"]#含量
        self._dehydrator_eff = parameters["dehydrator_eff"]#脱水器传热效率
        self._hydrator_eff = parameters["hydrator_eff"]#水合器传热效率
        self._P_amb = parameters["p_amb"]#环境压力
        self._T_amb = parameters["t_amb"]#环境温度    
        self._industrial_waste_heat_t = parameters["industrial_waste_heat_t2"]
        self._delta_H_Tref = -178e3  # J/mole 反应热

        self._steam_pressure_loss_ratio = parameters["steam_pressure_loss_ratio"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mechanical_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_dehydrator_distance = parameters["storage_dehydrator_distance"]
        self.deltaTmin_SSHX = parameters["deltaTmin_SSHX"]   #固-固换热器最小温差
        self.deltaTmin_SGHX = parameters["deltaTmin_SGHX"]   #固-气换热器最小温差


    def solve(self,inputs):
        self._res = self._de.solve(inputs)
        self._P_hydr = inputs["p_Hydr"] #反应器压力
        #self._Hydr_ot = inputs["Hydr_overheating_temperature"]
        self._T_hydr = 650#反应器温度
        results = {}
        results["BH"] = self._res["BH"]
        results["Dehy"] = self._res
        results["Caes"] = self._res["cost"]
        del results["Dehy"]["BH"]
        del results["Dehy"]["cost"]    
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        flue_gas_name2 = self._pw.get_flue_gas_refprop_name2()
        results["ht_HEN"],ht_loss,a_gas,results["ht_hen"]=self.high_tem_HEN(self._cao_conversion,
                                                               flue_gas_name,
                                                               flue_gas_name2)
        self.Hydr_cao_in  = results["ht_HEN"]
        self.Hydr_steam_in = self._T_hydr-self.deltaTmin_SGHX

        self._t= inputs["T_X"]
        self._reactants_ES_mole = results["Dehy"]["dehydrator"]["mole_dehydrator_reactions"]
        results["ht_loss"]=ht_loss*self._reactants_ES_mole
        results["energy_gas"]=a_gas*self._reactants_ES_mole
        #Basic input data
        results["Hydr"] = {}
        initialvalue = self.initialvalue(inputs)
        results["Hydr"]["initialvalue"] = initialvalue 
        #Hydrator
        hydrator = self.hydrator(results["Hydr"]["initialvalue"]["mole_hydrator_reactions"],
                                 flue_gas_name,
                                 flue_gas_name2)
        results["Hydr"]["hydrator"] = hydrator
        #brayton
        self._Hydrator_heat = results["Hydr"]["hydrator"]["Q_produce"]

        results["Bra"]=self._B.solve(self._Hydrator_heat,inputs)
        steam_name = flue_gas_name
        steam_pi = self._P_hydr
        steam_po = self._P_amb/(1-self._steam_pressure_loss_ratio)
        steam_blower = self.steam_blower_power(self._T_hydr,
                                         steam_pi,
                                         steam_po,
                                         results["Hydr"]["hydrator"]["mass_steam_in"],
                                         steam_name)
        results["steam_blower"] = steam_blower
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["Hydr"]["hydrator"]["m_camix_out"],
            results["Hydr"]["hydrator"]["m_camix_in"])*(-1)
        #results["Brayton"] = self._Bray.solve(self._Hydrator_heat,inputs)
        #Case all
        B_case ,a= self.B_case(results)
        results["Case_B"] = B_case
        results["Case_1"] = a
        All_case = self.All_case(results)
        results["Case_All"] = All_case
        return results  
    
    def equilibrium(self):
        p = self._P_hydr
        t = (-12845/((math.log(p/1e5))-16.508))-273.15
        return t
    def high_tem_HEN(self,X,flue_gas_name,flue_gas_name2):
        cp_cao_i = self._pw.cp0mass_mean("cao", self._T_hydr-self.deltaTmin_SSHX, self._T_amb)
        cp_caco3_o = self._pw.cp0mass_mean("caco3", self._T_hydr, self._T_amb+self.deltaTmin_SSHX)
        fluid=flue_gas_name
        fluid2=flue_gas_name2
        h_steam_in = CP.PropsSI('H', 'T', self._T_hydr-self.deltaTmin_SGHX+273.15, 'P', self._P_amb, fluid)
        h_steam_out = CP.PropsSI('H', 'T', self._industrial_waste_heat_t+273.15, 'P', self._P_amb, fluid)
        h_steam_in2 = CP.PropsSI('H', 'T', self._T_hydr+273.15, 'P', self._P_amb, fluid2)
        h_steam_out2 = CP.PropsSI('H', 'T', self._industrial_waste_heat_t+20+273.15, 'P', self._P_amb, fluid2)
        h_steam_in22 = CP.PropsSI('H', 'T', self._industrial_waste_heat_t+273.15, 'P', self._P_amb, fluid2)

        hen={}
        hen["out_caco3"]=cp_caco3_o*M_caco3*(self._T_hydr-self._T_amb-self.deltaTmin_SSHX)*X
        hen["out_cao"]=cp_cao_i*M_cao*(self._T_hydr-self._T_amb-self.deltaTmin_SSHX)*(1-X)
        hen["out_as"]=(X/(0.1338*0.9)-X)*(h_steam_in2-h_steam_out2)*M_2
        hen["in_bs"]=(h_steam_in-h_steam_out)*M_1*X/(0.1338*0.9)
        a=(cp_caco3_o*M_caco3*(self._T_hydr-self._T_amb-self.deltaTmin_SSHX)*X
           +cp_cao_i*M_cao*(self._T_hydr-self._T_amb-self.deltaTmin_SSHX)*(1-X))#re
        a_s =(X/(0.1338*0.9)-X)*(h_steam_in2-h_steam_out2)*M_2
        b_s=cp_cao_i*M_cao#冷
        b_g=(h_steam_in-h_steam_out)*M_1*X/(0.1338*0.9)
        t_s_out = self._T_amb+((a+a_s)*0.96-b_g)/b_s

        a_gas=(h_steam_out2-h_steam_in22)*(X/(0.1338*0.9)-X)*M_2
        a1=self._T_hydr-self.deltaTmin_SSHX
        if t_s_out> a1:
            t_s_out=a1
        else:
            t_s_out=t_s_out
        hen["in_cao"]=b_s*(t_s_out-self._T_amb)
        h_mole_lost=(a+a_s-(t_s_out-self._T_amb)*b_s-b_g)
        return t_s_out,h_mole_lost,a_gas,hen
    def initialvalue(self,input):
        results = {}
        results["cao_conversion"] = self._cao_conversion
        results["dehydrator_eff"] = self._dehydrator_eff
        results["P_amb"] = self._P_amb
        results["T_amb"] = self._T_amb
        results["T_hydr"]=self._T_hydr
        results["delta_H_Tref"] = self._delta_H_Tref
        results["Time_scale"] = input["T_X"]
        results["mole_hydrator_reactions"] = self._reactants_ES_mole*input["T_X"]
        return results 
    def hydrator(self,mole_re,flue_gas_name,flue_gas_name2):
        delta_H_Tr = -self._mole_reaction_heat(self._T_hydr, self._P_hydr)
        Q_produce = (delta_H_Tr*mole_re*self._cao_conversion-95e3*mole_re*0.05)

        results={}
        results["is_succeed"] = 1  
        results["delta_H_Tr"] = delta_H_Tr 
        results["Q_reactions"] = Q_produce
        results["mole_hydrator_reactions"]=mole_re

        results["mole_Hydration_in"] = mole_re
        results["mole_caco3_o"] = mole_re*self._cao_conversion+mole_re/0.9-mole_re
        results["mass_caco3_o"] = M_caco3*results["mole_caco3_o"]
        results["mole_cao_out"] = mole_re*(1-self._cao_conversion)
        results["mass_cao_out"] = mole_re*M_cao*(1-self._cao_conversion)  
        all =results["mole_caco3_o"]+results["mole_cao_out"]
        results["caco3"] = all*0.1*M_cao+all*0.9*M_caco3

        results["mole_cao_in"] = mole_re
        results["mass_cao_in"] = mole_re*M_cao
        results["mole_caco3_in"] = mole_re/0.9-mole_re
        results["mass_caco3_in"] = M_caco3*results["mole_caco3_in"]
        results["mole_steam_in"] = mole_re*self._cao_conversion/(0.1338*0.9)
        results["mass_steam_in"] = results["mole_steam_in"]*M_1
        results["mass_steam_out"] = (results["mole_steam_in"]-mole_re*self._cao_conversion)*M_2
        results["m_camix_in"]= results["mass_cao_in"]+results["mass_caco3_in"]
        results["m_camix_out"] = results["mass_caco3_o"]+results["mass_cao_out"]

        heat_cao =(self._pw.cp_camix_mean_Co(self.Hydr_cao_in,self._T_hydr,0.95,0.9)*
                   (self._T_hydr-self.Hydr_cao_in)*results["m_camix_in"])
        fluid=flue_gas_name
        fluid2=flue_gas_name2
        heat_steam=((CP.PropsSI('H', 'T', self._T_hydr+273.15, 'P', self._P_hydr,fluid)
                    -CP.PropsSI('H', 'T', self.Hydr_steam_in+273.15, 'P', self._P_hydr, fluid))
                    *results["mass_steam_in"])
        results["Q_produce"] = Q_produce -heat_cao- heat_steam
        results["exergy"] = {}       
        results["exergy"]["chemical_t"]=self._de.ex_calculations1( self._T_hydr,Q_produce)
        results["exergy"]["sensible_heat"] = (self._de.ex_calculations(self.Hydr_cao_in, self._T_hydr,heat_cao)
                                              +(self.E_steam(self._T_hydr,self._P_hydr)-
                                                self.E_steam(self.Hydr_steam_in,self._P_hydr))
                                                *results["mass_steam_in"])

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
    
    def E_steam(self,T,P):
        H1=CP.PropsSI('H', 'T', T+273.15, 'P', P, "REFPROP::co2")
        S1=CP.PropsSI('S', 'T', T+273.15, 'P', P, "REFPROP::co2")
        H0=CP.PropsSI('H', 'T', self._T_amb+273.15, 'P', self._P_amb, "REFPROP::co2")
        S0=CP.PropsSI('S', 'T', self._T_amb+273.15, 'P', self._P_amb, "REFPROP::co2")
        a=H1-H0-(self._T_amb+273.15)*(S1-S0)
        return a 

    def ex_calculations(self,b1,b2,Q):
        a1=b1+273.15
        a2=b2+273.15
        gcpw=(a1-a2-293.15*math.log(a1/a2))/(a1-a2)
        exergy = Q*gcpw
        return exergy
    def write_pyPinch_data_csv(self, path):
        data = {}
        data["CP"] = self._pinch_point_data["CP"]
        data["TSUPPLY"] = self._pinch_point_data["TSUPPLY"]
        data["TTARGET"] = self._pinch_point_data["TTARGET"]
        df = pd.DataFrame(data)
        df["CP"] = df["CP"]/1000
        df.to_csv(path, index=False)
        with open(path, "r+") as fp:
            lines = fp.readlines()
            lines.insert(0, f'Tmin, {self._T_delta_pinch},\n')
            fp.seek(0)
            fp.writelines(lines)
    def write_pyPinch_data_text(self):
        data = {}
        data["CP"] = self._pinch_point_data["CP"]
        data["TSUPPLY"] = self._pinch_point_data["TSUPPLY"]
        data["TTARGET"] = self._pinch_point_data["TTARGET"]
        df = pd.DataFrame(data)
        df["CP"] = df["CP"]/1000
        text = df.to_string(index=False)
        text = f'Tmin {self._T_delta_pinch} \n'+text
        return text
    def solvehen(self, input_text):
        pinch = PyPinch(input_text)
        # pinch.shiftTemperatures()
        # pinch.constructTemperatureInterval()
        # pinch.constructProblemTable()
        # pinch.constructHeatCascade()
        # pinch.constructShiftedCompositeDiagram()
        # pinch.constructCompositeDiagram()
        # pinch.constructGrandCompositeCurve()
        # HIntervalTable=pinch.constructEnthalpyIntervalTable()
        pinch.solve()
        hot_util = pinch.hotUtility*1e3  # W
        cold_util=pinch.coldUtility*1e3  # W
        HIntervalTable=pinch.EnthaphyIntervalTable
        streamPropertyTable=[]
        for i,stream in enumerate(pinch.streams):
            record={}
            record["type"]=stream["type"]
            record["MCP"]=stream["cp"]
            record["HTC"]=self._HTCs[self._materials[i]]
            streamPropertyTable.append(record)
        for record in HIntervalTable:
            hotStreams=record["hotStreams"]
            hotHTCs=[streamPropertyTable[i]["HTC"] for i in hotStreams]
            hotMCPs=[streamPropertyTable[i]["MCP"] for i in hotStreams]
            coldStreams=record["coldStreams"]
            coldHTCs=[streamPropertyTable[i]["HTC"] for i in coldStreams]
            coldMCPs=[streamPropertyTable[i]["MCP"] for i in coldStreams]
            record["hotQi"]=[self._myRound(cp*(record["hotshiftedTs"]-record["hotshiftedTt"])) for cp in hotMCPs]
            record["coldQi"]=[self._myRound(cp*(record["coldshiftedTt"]-record["coldshiftedTs"])) for cp in coldMCPs]
            hotHXA=[]
            coldHXA=[]
            for i,Qi in enumerate(record["hotQi"]):
                hotHXA.append(self._myRound(Qi/record["LMTD"]/hotHTCs[i]*1000))
            for i,Qi in enumerate(record["coldQi"]):
                coldHXA.append(self._myRound(Qi/record["LMTD"]/coldHTCs[i]*1000))
            record["hotHXA"]=hotHXA
            record["coldHXA"]=coldHXA

        totalHXA=0
        for record in HIntervalTable:
            for HXA in record["hotHXA"]:
                totalHXA=totalHXA+HXA
            for HXA in record["coldHXA"]:
                totalHXA=totalHXA+HXA

        return hot_util, cold_util,totalHXA
    def hot_exchange(self,mole):
        cp_cao = self._pw.cp0mass_mean("cao", self._T_hydr-self.deltaTmin_SSHX, self._T_amb)
        cp_caoh2 = self._pw.cp0mass_mean("caoh2", self._T_hydr, self._T_amb+self.deltaTmin_SSHX)
        h_steam_out = CP.PropsSI('H', 'T', self._T_hydr-self.deltaTmin_SGHX, 'P', self._P_amb, "REFPROP::water")
        h_steam_in = CP.PropsSI('H', 'T', 105+273.15, 'P', self._P_amb, "REFPROP::water")
        h_water_out = CP.PropsSI('H', 'T', 95+273.15, 'P', self._P_amb, "REFPROP::water")
        h_water_in = CP.PropsSI('H', 'T', self._T_amb+273.15, 'P', self._P_amb, "REFPROP::water")

        hot_h = cp_caoh2*mole*M_caoh2*(self._T_hydr-self._T_amb-self.deltaTmin_SSHX)
        hot_c = cp_cao*mole*M_cao*(self._T_hydr-self._T_amb-self.deltaTmin_SSHX)+((h_steam_out-h_steam_in)+(h_water_out-h_water_in))*mole*M_H2O
        results={}
        results["hot_h"] = hot_h
        results["hot_c"] = hot_c
        results["h_lost"]=(hot_c/0.96)*0.04
        results["heat_out"] =hot_h*0.96-hot_c
        return results
    def _myRound(self,a,ndigital=3):
        return round(a,ndigital) 
    def Heating_water(self,h1,h2,h3):
        h=h1+h2+h3

        h_1=CP.PropsSI('H', 'T', 85+273.15, 'P', self._P_amb, "REFPROP::water")
        h_2=CP.PropsSI('H', 'T', 60+273.15, 'P', self._P_amb, "REFPROP::water")
        m = h/(h_1-h_2)
        return m
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
        results["power"] = W*(-2)
        return results
    
    def conveying_power(self, m_camix_in,m_camix_o):
        return self._convey_consumption*self._storage_dehydrator_distance * \
            (m_camix_in+m_camix_o)
    def B_case(self,results):
        b={}
        b["power_lost"] = (-results["conveying_power"]-results["steam_blower"]["power"])
        b["power_out"]=results["Bra"]["evaluation_indicators"]["power"]-b["power_lost"]
        b["hot_stockpile"]=-self._delta_H_Tref*results["Hydr"]["initialvalue"]["mole_hydrator_reactions"]*self._cao_conversion
        b["hot_to_B"]=self._Hydrator_heat
        b["h_lost"]=b["hot_stockpile"]-b["power_out"]
        b["energy_eff"]=b["power_out"]/b["hot_stockpile"]
        b["exergy_eff"]=b["power_out"]/+results["Hydr"]["hydrator"]["exergy"]["chemical_t"]

        a={}
        a["power"]=-b["power_out"]/1000*8*1.025
        a["co2"]=-(results["Hydr"]["hydrator"]["mass_steam_in"]-
                    results["Hydr"]["hydrator"]["mass_steam_out"])*3.6*8*51.23
        a["caco3"] = 761.706*results["Hydr"]["hydrator"]["caco3"]**3.6*8
        a["cost"]=a["power"]+a["co2"]-a["caco3"]

        a["cha_power"] = results["Caes"]["power"]+results["Hydr"]["hydrator"]["caco3"]*258460/1000*8*0.425
        a["cha_hot"] = results["Caes"]["hot"]
        a["cha_caoh2"] =results["Caes"] ["caoh2"]
        a["cost_all"] = a["cost"]+a["cha_power"]+a["cha_hot"]+a["cha_caoh2"] 
        return b,a
    
    def All_case(self,results):
        a={}
        a["power_in"] =results["Dehy"]["evaluation_indicators"]["power_in1"]
        a["power_grind"] = results["Hydr"]["hydrator"]["caco3"]*258460
        a2=(-self._delta_H_Tref*self._reactants_ES_mole*self._cao_conversion+(
            CP.PropsSI('H', 'T',self._industrial_waste_heat_t+273.15, 'P', self._P_amb,"REFPROP::co2")-
            CP.PropsSI('H', 'T',self._T_amb+273.15, 'P', self._P_amb,"REFPROP::co2")
        )*self._reactants_ES_mole*self._cao_conversion*M_co2)
        a["c_energy_out"] = results["Dehy"]["evaluation_indicators"]["hot_stockpile"]*(1-self._cao_conversion)
        a["c_energy_in"] = a2
        a["c_energy_in1"]=a2-results["Dehy"]["evaluation_indicators"]["hot_stockpile"]*self._cao_conversion
        a["hot_in"] = results["BH"]["evaluation_indicators"]["hot_cost"]
        a["hot_out"] = results["Dehy"]["evaluation_indicators"]["ht_Q"]
        a["power"] = results["Case_B"]["power_out"]
        a["Round-trip_eff"] = a["power"]/(a["power_in"]+a["power_grind"])
        a["energy_eff"] = (a["hot_out"]+a["power"])/(a["power_in"]+a["hot_in"]+a["c_energy_in1"]+a["power_grind"])
        a3=(CP.PropsSI('H', 'T',self._industrial_waste_heat_t+273.15, 'P', self._P_amb,"REFPROP::co2")-
            CP.PropsSI('H', 'T',self._T_amb+273.15, 'P', self._P_amb,"REFPROP::co2"))*self._reactants_ES_mole*self._cao_conversion*M_co2
        a["c_exergy"]=((C_co2-C_caco3+C_caoh2)*self._cao_conversion*self._reactants_ES_mole+
                       self._de.ex_calculations(self._industrial_waste_heat_t, self._T_amb,a3))
        a["exergy_eff"] = (a["hot_out"]*0.15154341+a["power"])/(
            results["Dehy"]["evaluation_indicators"]["exergy_in"]+ a["c_exergy"]+a["power_grind"])
        a["C_BH"]=(results["Hydr"]["hydrator"]["mass_steam_in"]-
                   results["Hydr"]["hydrator"]["mass_steam_out"])/a["power"]*3.6e6
        return a



if __name__ == '__main__':
    parameters = dict()
    flue_gas_composistion = {}
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    flue_gas_composistion2 = {}
    flue_gas_composistion2["co2"] = 0.01338
    flue_gas_composistion2["o2"] = 0.0384
    flue_gas_composistion2["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["flue_gas_composition2"] = flue_gas_composistion2
    parameters["isentropic_eff_mc"] = 0.88
    parameters["t_isentropic_eff_mc"] = 0.92
    parameters["mechanical_eff"] = 0.98   #机械效率
    parameters["min_temperature_exchange"] = 15 
    parameters["deltaTmin_SSHX"] = parameters["min_temperature_exchange"]+5   #固-固换热器最小温差
    parameters["deltaTmin_SGHX"] = parameters["min_temperature_exchange"]   #固-气换热器最小温差
    parameters["industrial_waste_heat_t"] =350 #℃
    parameters["industrial_waste_heat_t2"] =180 #℃
    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_amb"] = 20   #环境温度
    parameters["p_amb"] = 101325   #环境压力
    parameters["p_amb"] = 101325   #环境压力
    parameters["p_bray_L"] = 7.5e6#热泵低压
    parameters["Store_electrical_power"]= 1e6


    parameters["cao_conversion"] = 0.95  #氧化钙转化率
    parameters["cao_purity"] = 0.9 #氢氧化钙含量
    parameters["dehydrator_eff"] = 0.97   #脱水器效率
    parameters["steam_pressure_loss_ratio"] = 0.01
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100

    parameters["hydrator_eff"] = 0.95   #水合器器效率
    parameters["cao_conversion2"] = 0.4  #氧化钙转化率
    parameters["p_bray_L_B"] = 7.5e6

    inputs={}
    inputs["p_bray_H"] = 17794921.326326743#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12258207.576473976#优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 1e5 #变量4，反应器压力
    inputs["t_HEN"]= 453.15530148433106
    #inputs["m_water"] = 5.708391319264835 
    inputs["m_water"] = 0
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度

    inputs["p_bray_H_B"] = 30e6
    inputs["p_bray_MH_B"] = 16217752.142109105
    inputs["p_bray_ML_B"] = 12217752.142109105
    
    inputs["p_Hydr"] = 1e5
    inputs["R"]=0.37
    inputs["HT"] = 650
    inputs["H_out"]=214

    inputs["Hydr_cao_in"]=440
    inputs["T_X"] = 1


    calre = Hydrator(parameters)
    results = calre.solve(inputs)
    print(results)



