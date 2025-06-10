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
from pyBraytonHeatPump import BraytonHeatPump
from pyPinch import PyPinch
import pandas as pd
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
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"],
                                   parameters["flue_gas_composition2"])
        self._bhp = BraytonHeatPump(parameters)
        self._cao_conversion = parameters["cao_conversion"]#氢氧化钙转化率
        self._X = self._cao_conversion
        self._cao_purity = parameters["cao_purity"]#含量
        self._Y = self._cao_purity
        self._dehydrator_eff = parameters["dehydrator_eff"]#脱水器传热效率
        self._P_amb = parameters["p_amb"]#环境压力
        self._T_amb = parameters["t_amb"]#环境温度

        self._delta_H_Tref = 104e3  # J/mole CaO

        self._steam_pressure_loss_ratio = parameters["steam_pressure_loss_ratio"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mechanical_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_dehydrator_distance = parameters["storage_dehydrator_distance"]
        self.deltaTmin_SSHX = parameters["deltaTmin_SSHX"] = 25   #固-固换热器最小温差
        self.deltaTmin_SGHX = parameters["deltaTmin_SGHX"] = 20   #固-气换热器最小温差

    def solve(self,input):
        results = {}
        self._res = self._bhp.solve(input)
        results["BH"] = self._res
        self._p_dehy = input["p_Dehy"]
        self._T_dehy = 525#反应器温度
        initialvalue = self.initialvalue()
        results["initialvalue"] = initialvalue   

        dehydrator = self.dehydrator(input["t_HEN"],
                                             self._T_dehy,
                                             self._P_amb,
                                             self._cao_conversion,
                                             self._cao_purity)
        results["dehydrator"] = dehydrator
        aaa=dehydrator["mole_Dehydration_in"]
        self._a=aaa
        self._m_camix_out_dehy=dehydrator["out"]["m_camix"]
        self._m_steam = dehydrator["out"]["m_steam"]
        self._t_HEN = input["t_HEN"]
        self._m_camix_in_dehy = dehydrator["in"]["m_camix"]
        #High temperature section of the heat exchange network
        self._T_delta_pinch = 15
        self._T_vaporization = CP.PropsSI('T', 'P', self._p_dehy, 'Q', 1, "REFPROP::water")-273.15
        self._pinch_point_data,hot_out,hot_hen,cold_hen,exergy_hen=self.pinch_point_data(input["m_water"])
        hen_text = self.write_pyPinch_data_text()
        hot_util, cold_util,total_HENA = self.solvehen(hen_text)
        results["HEN"]={}
        results["HEN"]["pinch_analysis_text"]=hen_text
        results["HEN"]["total_HEN_area"]=total_HENA
        results["HEN"]["hot_utility"] = hot_util
        results["HEN"]["cold_utility"] = cold_util
        results["HEN"]["p"]=self._pinch_point_data
        results["HEN"]["hot_out"] =hot_out
        results["HEN"]["hot"] =hot_hen
        results["HEN"]["cold"] =cold_hen
        results["HEN"]["exergy"] =exergy_hen

        #results["ht_HEN"],ht_loss,ht_m1,ht_Q=self.high_tem_HEN(self._cao_conversion)
       # dehydrator(self, Ti_flue_gas, Ti_cao, Ti_water,To_water,Tcarb, pcarb, X):
        
        results["ht_loss"]=hot_hen-cold_hen
        results["ht_m1"]=input["m_water"]
        results["ht_Q"]=hot_out

        #Phase change part of the low temperature section of the heat exchange network
        steam_name = "water"
        steam_pi = self._p_dehy
        steam_po = self._P_amb/(1-self._steam_pressure_loss_ratio)
        steam_blower = self.steam_blower_power(self._T_dehy,
                                         steam_pi,
                                         steam_po,
                                         self._m_steam,
                                         steam_name)
        results["steam_blower"] = steam_blower
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["dehydrator"]["in"]["m_caoh2"],
            results["dehydrator"]["out"]["m_camix"])*(-1)
        #results["hot_water"]=self.hot_water(self._res["heat_recovery"]["t_flue_gas_out"],
        #                                    self._res["evaluation_indicators"]["flue_gas_mass_flow"],
        #                                    results["ht_m1"],
        #                                    results["ht_Q"])
        evaluation_indicators ,b= self.evaluation_indicators(results)
        results["evaluation_indicators"] = evaluation_indicators
        results["cost"]=b
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
    
    def pinch_point_data(self,m):
        self._materials=["Ca","Gas","Gas","Ca","Water","Water","Water","Gas"]
        self._HTCs={"Ca":300,"Gas":600,"Water":2500}
        pinch_point_data = {}
        pinch_point_data["TSUPPLY"] = {}
        pinch_point_data["TTARGET"] = {}
        pinch_point_data["ENERGY"] = {}
        pinch_point_data["FLOWRATE"] = {}
        pinch_point_data["CP"] = {}
        # H1 cao
        pinch_point_data["TSUPPLY"]["H_CaM"] = self._T_dehy
        pinch_point_data["TTARGET"]["H_CaM"] = self._T_amb
        pinch_point_data["FLOWRATE"]["H_CaM"] = self._m_camix_out_dehy*0.96
        pinch_point_data["CP"]["H_CaM"] = self._m_camix_out_dehy*self._pw.cp_camix_mean_Co(
            self._T_amb,self._T_dehy,self._X,self._Y)
        pinch_point_data["ENERGY"]["H_CaM"] =pinch_point_data["CP"]["H_CaM"]*(
            self._T_dehy-self._T_amb)
        
        # H2: Gas out:Steam
        pinch_point_data["TSUPPLY"]["H_Steam"] = self._T_dehy
        pinch_point_data["TTARGET"]["H_Steam"] = self._T_vaporization+1
        pinch_point_data["FLOWRATE"]["H_Steam"] = self._m_steam*0.96
        pinch_point_data["ENERGY"]["H_Steam"] = self._m_steam * \
            (CP.PropsSI('H', 'T', self._T_dehy+273.15,
                        'P',self._p_dehy, "water") -
             CP.PropsSI('H', 'T', self._T_vaporization+1+ 273.15,
                        'P', self._p_dehy, "water"))
        pinch_point_data["CP"]["H_Steam"] = pinch_point_data["ENERGY"]["H_Steam"] / \
            (self._T_dehy-self._T_vaporization-1)

        # H3: water after phase change
        pinch_point_data["TSUPPLY"]["PC_water"] = self._T_vaporization-1
        pinch_point_data["TTARGET"]["PC_water"] = self._T_amb
        pinch_point_data["FLOWRATE"]["PC_water"] = self._m_steam*0.96
        pinch_point_data["ENERGY"]["PC_water"] = self._m_steam * \
            (CP.PropsSI('H', 'T', self._T_vaporization-1+273.15,
                        'P', self._p_dehy, "water") -
             CP.PropsSI('H', 'T', self._T_amb + 273.15,
                        'P', self._p_dehy, "water"))
        pinch_point_data["CP"]["PC_water"] = pinch_point_data["ENERGY"]["PC_water"] / \
            (self._T_vaporization-1-self._T_amb)
        
        # H4: water        phase change
        pinch_point_data["TSUPPLY"]["PC_w"] = self._T_vaporization+1
        pinch_point_data["TTARGET"]["PC_w"] = self._T_vaporization-1
        pinch_point_data["FLOWRATE"]["PC_w"] = self._m_steam*0.96
        pinch_point_data["ENERGY"]["PC_w"] = self._m_steam * \
            (CP.PropsSI('H', 'T', self._T_vaporization+1+273.15,
                        'P', self._p_dehy, "water") -
             CP.PropsSI('H', 'T', self._T_vaporization-1+273.15,
                        'P', self._p_dehy, "water"))
        pinch_point_data["CP"]["PC_w"] = pinch_point_data["ENERGY"]["PC_w"] / \
            (2)

        ## C1: Camix out: Ca(OH)2 
        pinch_point_data["TSUPPLY"]["C_Caoh2"] = self._T_amb
        pinch_point_data["TTARGET"]["C_Caoh2"] = self._t_HEN
        pinch_point_data["FLOWRATE"]["C_Caoh2"] = self._m_camix_in_dehy
        pinch_point_data["CP"]["C_Caoh2"] = (self._m_camix_in_dehy*
                                             self._pw.cp_camix_mean_Ci(
                                                self._T_amb,self._t_HEN,self._Y))
        pinch_point_data["ENERGY"]["C_Caoh2"] =pinch_point_data["CP"]["C_Caoh2"]*(self._t_HEN-self._T_amb)

        ## C2: water
        pinch_point_data["TSUPPLY"]["C_water"] = 60
        pinch_point_data["TTARGET"]["C_water"] = 85
        pinch_point_data["FLOWRATE"]["C_water"] = m
        pinch_point_data["ENERGY"]["C_water"] = m * \
            (CP.PropsSI('H', 'T', 85+273.15,
                        'P', self._P_amb, "REFPROP::water") -
            CP.PropsSI('H', 'T', 60+273.15,
                        'P', self._P_amb, "REFPROP::water"))
        pinch_point_data["CP"]["C_water"] = pinch_point_data["ENERGY"]["C_water"] / \
            (25)
       
        self._pinch_point_data = pinch_point_data
        a =(pinch_point_data["ENERGY"]["H_CaM"] + 
            pinch_point_data["ENERGY"]["H_Steam"] + 
            pinch_point_data["ENERGY"]["PC_water"]+
            pinch_point_data["ENERGY"]["PC_w"])/0.96
        b = (pinch_point_data["ENERGY"]["C_Caoh2"] + 
             pinch_point_data["ENERGY"]["C_water"])
        hot_out = pinch_point_data["ENERGY"]["C_water"]

        exergy={}
        exergy["H_CaM"]=self.ex_calculations(self._pinch_point_data["TSUPPLY"]["H_CaM"],
                                             self._pinch_point_data["TTARGET"]["H_CaM"],
                                             self._pinch_point_data["ENERGY"]["H_CaM"])
        exergy["H_Steam"]=self.ex_calculations(self._pinch_point_data["TSUPPLY"]["H_Steam"],
                                             self._pinch_point_data["TTARGET"]["H_Steam"],
                                             self._pinch_point_data["ENERGY"]["H_Steam"])
        exergy["PC_water"]=self.ex_calculations(self._pinch_point_data["TSUPPLY"]["PC_water"],
                                             self._pinch_point_data["TTARGET"]["PC_water"],
                                             self._pinch_point_data["ENERGY"]["PC_water"])
        exergy["PC_w"]=self.ex_calculations(self._pinch_point_data["TSUPPLY"]["PC_w"],
                                             self._pinch_point_data["TTARGET"]["PC_w"],
                                             self._pinch_point_data["ENERGY"]["PC_w"])
        exergy["C_Caoh2"]=self.ex_calculations(self._pinch_point_data["TSUPPLY"]["C_Caoh2"],
                                             self._pinch_point_data["TTARGET"]["C_Caoh2"],
                                             self._pinch_point_data["ENERGY"]["C_Caoh2"])
        exergy["C_water"]=self.ex_calculations(self._pinch_point_data["TSUPPLY"]["C_water"],
                                             self._pinch_point_data["TTARGET"]["C_water"],
                                             self._pinch_point_data["ENERGY"]["C_water"])
        exergy["hot_in"]=(exergy["H_CaM"]+exergy["H_Steam"]
                          +exergy["PC_water"]+exergy["PC_w"])/0.96
        exergy["cold_add"]=exergy["C_Caoh2"]+exergy["C_water"]
        exergy["lost"]=exergy["hot_in"]-exergy["cold_add"]

        return pinch_point_data,hot_out,a,b,exergy

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
    def _myRound(self,a,ndigital=3):
        return round(a,ndigital) 
    
    def high_tem_HEN(self,X):
        c=self._T_amb
        cp_cao_o = self._pw.cp0mass_mean("cao", self._T_dehy, 105)
        cp_caco3_o = self._pw.cp0mass_mean("caoh2", self._T_dehy, 105)
        cp_caco3_i = self._pw.cp0mass_mean("caoh2", self._T_dehy-80, 75)
        h_steam_in = CP.PropsSI('H', 'T', self._T_dehy+273.15, 'P', self._P_amb, "REFPROP::water")
        h_steam_out = CP.PropsSI('H', 'T', 105+273.15, 'P', self._P_amb, "REFPROP::water")
        h_steam_1 = CP.PropsSI('H', 'T', 95+273.15, 'P', self._P_amb, "REFPROP::water")
        h_steam_2 = CP.PropsSI('H', 'T', 25+273.15, 'P', self._P_amb, "REFPROP::water")

        a=cp_caco3_i*M_caoh2#cool
        b=((cp_caco3_o*M_caoh2*(1-X))+(cp_cao_o*M_cao*X))*(self._T_dehy-105)+(h_steam_in-h_steam_out)*M_H2O*X
        t_out = 75+b*0.96/a
        a1=self._T_dehy-self.deltaTmin_SSHX
        if t_out> a1:
            t_out=a1
        else:
            t_out=t_out

        hot=(self._pw.cp0mass_mean("cao", self._T_dehy, c)*X*(self._T_dehy-c)*M_cao+
             self._pw.cp0mass_mean("caoh2", self._T_dehy,c)*(1-X)*(self._T_dehy-c)*M_caoh2+
             (h_steam_in-h_steam_out+h_steam_1-h_steam_2)*M_H2O*X)
        cold=(self._pw.cp0mass_mean("caoh2", t_out, c)*(t_out-c)*M_caoh2)

        hotcp=(h_steam_out-h_steam_1)*M_H2O*X
        h_steam_3 = CP.PropsSI('H', 'T', 85+273.15, 'P', self._P_amb, "REFPROP::water")
        h_steam_4 = CP.PropsSI('H', 'T', 60+273.15, 'P', self._P_amb, "REFPROP::water")
        m1 = hotcp*0.96/(h_steam_3-h_steam_4)
        h_mole_lost=hot-cold+hotcp*0.04
        return t_out,h_mole_lost,m1,hotcp*0.96
    
    def dehydrator(self,T_solid_in, Tdehy, pdehy, X , Y):
        
        delta_H_Tr = self._mole_reaction_heat(Tdehy, pdehy)
        heat_caoh2 = ((self._pw.cp0mass_mean("caoh2", Tdehy, T_solid_in)*M_caoh2*(Tdehy-T_solid_in)*Y)+
                      (self._pw.cp0mass_mean("caco3", Tdehy, T_solid_in)*M_caco3*(Tdehy-T_solid_in)*(1-Y)))

        power = (delta_H_Tr*X+heat_caoh2/Y)/self._dehydrator_eff
        n=self._res["evaluation_indicators"]["hot_output"]/power

        results ={}
        results["is_succeed"] = 1  
        results["delta_H_Tr"] = delta_H_Tr 
        results["heat_caco3"] = heat_caoh2
        results["mole_dehydrator_reactions"]=X*n
        results["mole_Dehydration_in"] = 1*n
        results["in"]={}
        results["in"]["m_caoh2"] = results["mole_Dehydration_in"]*M_caoh2
        results["in"]["mole_caoh2"] = results["mole_Dehydration_in"]
        results["in"]["m_caco3"] = (results["mole_Dehydration_in"]/Y-results["mole_Dehydration_in"])*M_caco3
        results["in"]["mole_caco3"] = results["mole_Dehydration_in"]/Y-results["mole_Dehydration_in"]
        results["in"]["m_camix"] = results["in"]["m_caoh2"]+results["in"]["m_caco3"]
        results["out"]={}
        results["out"]["mole_cao"]=results["mole_Dehydration_in"]*X
        results["out"]["m_cao"]=results["mole_Dehydration_in"]*X*M_cao
        results["out"]["m_caoh2"] = results["mole_Dehydration_in"]*(1-X)*M_caoh2
        results["out"]["mole_caoh2"] = results["mole_Dehydration_in"]*(1-X)
        results["out"]["m_caco3"] =results["in"]["m_caco3"] 
        results["out"]["mole_caco3"] = results["in"]["mole_caco3"]

        results["out"]["m_camix"] =results["out"]["m_cao"]+results["out"]["m_caoh2"]+results["out"]["m_caco3"]
        results["out"]["m_steam"]= results["mole_Dehydration_in"]*X*M_H2O
        results["out"]["mole_steam"] = results["mole_Dehydration_in"]*X

        results["exergy"] = {}
        results["exergy"]["power"]=self._res["evaluation_indicators"]["exergy_out"]
        results["exergy"]["chemical_t"]=self.ex_calculations1( Tdehy,delta_H_Tr*X*n)
        results["exergy"]["sensible_heat"] = self.ex_calculations(T_solid_in, Tdehy,heat_caoh2*n)
        results["exergy"]["lost"] = (results["exergy"]["power"]-results["exergy"]["chemical_t"]-
                                     results["exergy"]["sensible_heat"])
    
        return results
    def ex_calculations(self,b1,b2,Q):
        a1=b1+273.15
        a2=b2+273.15
        gcpw=(a1-a2-293.15*math.log(a1/a2))/(a1-a2)
        exergy = Q*gcpw
        return exergy
    def ex_calculations1(self,a1,Q):
        gcpw=1-(293.15)/(a1+273.15)
        exergy = Q*gcpw
        return exergy
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
    
    def hot_water(self,t,m,m2,Q):
        h1=CP.PropsSI('H', 'T', t+273.15, 'P', self._P_amb, "REFPROP::water")
        h2=CP.PropsSI('H', 'T', 60+273.15, 'P', self._P_amb, "REFPROP::water")
        h11=CP.PropsSI('H', 'T', 85+273.15, 'P', self._P_amb, "REFPROP::water")
        a1=(h1-h2)*0.96*m
        m1=a1/(h11-h2)
        a={}
        a["m_hot"]=m1+m2
        a["hot_in"]=(h1-h2)*m
        a["hot_out"]=a1
        a["loss"]=(h1-h2)*(1-0.96)*m
        a["exergy_in"]=self.ex_calculations(t+273.15, 60+273.15,a["hot_in"])
        a["exergy_out"]=self.ex_calculations(85+273.15, 60+273.15,a1)
        a["exergy_loss"]= a["exergy_in"]- a["exergy_out"]
        a["Q"]=a1+Q
        return a

    def evaluation_indicators(self,results):
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        fluid=flue_gas_name
        a={}
        a["power_in1"] = (-results["conveying_power"]-results["steam_blower"]["power"]+1e6)
        a["energy_in"] = (a["power_in1"]+
                         self._res["evaluation_indicators"]["hot_cost"])
        a["hot_stockpile"] =self._delta_H_Tref*self._cao_conversion*self._a-598240/0.96
        a["ht_Q"]=results["ht_Q"]
        a["hot_lost"] = a["energy_in"]-a["hot_stockpile"]-a["ht_Q"]
        a["energy_eff"] = (a["hot_stockpile"]+a["ht_Q"])/(a["energy_in"])

        a["exergy_in"]=(-results["conveying_power"]-results["steam_blower"]["power"]+
                         self._res["evaluation_indicators"]["exergy_in"])
        a["exergy_out"] = (results["dehydrator"]["exergy"]["chemical_t"]+
                           self.ex_calculations(60+273.15, 85+273.15,results["ht_Q"]))
        a["exergy_eff"] = a["exergy_out"]/a["exergy_in"]

        b={}
        b["power"]=(-results["conveying_power"]-results["steam_blower"]["power"]+
                    self._res["evaluation_indicators"]["power_cost"])/1000*8*0.425
        #b["steam"]=self._res["evaluation_indicators"]["flue_gas_mass_flow"]*3.6*8*86
        b["hot"]=-(a["ht_Q"]-
                  self._res["evaluation_indicators"]["hot_cost"])*3600*8/1000000000*35.7
        b["caoh2"]=results["dehydrator"]["in"]["m_caoh2"]*3.6*8*660
        b["cost"]=b["power"]+b["caoh2"]-b["hot"]
        return a ,b


    
if __name__ == '__main__':
    parameters = dict()
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    flue_gas_composistion2 = dict()
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
    parameters["heat_transfer_loss_eff"] = 0.96
    parameters["t_amb"] = 20   #环境温度
    parameters["p_amb"] = 101325   #环境压力
    parameters["p_bray_L"] = 7.5e6#热泵低压
    parameters["Store_electrical_power"]= 1e6

    parameters["cao_conversion"] = 0.95  #氧化钙转化率
    parameters["cao_purity"] = 0.9 #氢氧化钙含量
    parameters["dehydrator_eff"] = 0.95   #脱水器效率
    parameters["steam_pressure_loss_ratio"] = 0.01
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100

    calcs = Dehydrator(parameters) 

    inputs={}
    inputs["p_bray_H"] = 17794921.326326743#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12258207.576473976#优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 101325  #变量4，反应器压力
    inputs["t_HEN"]=  453.15530148433106  
    inputs["m_water"] =5.708391319264835 
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度
    results = calcs.solve(inputs)
    print(results)

