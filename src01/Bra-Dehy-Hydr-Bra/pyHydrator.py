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
from pyDehydrator import Dehydrator
from pyBrayton import Brayton
from pyPinch import PyPinch
import pandas as pd
M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol
M_H2O = 18e-3  # kg/mol

class Hydrator(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
        self._de = Dehydrator(parameters)
        self._Bray = Brayton(parameters)
        self._Store_electrical_power = parameters["Store_electrical_power"] #1MW
        self._m1=6.542218380607665
        self._res = self._de.solve(self._m1)

        self._cao_conversion = parameters["cao_conversion"]#氢氧化钙转化率
        self._cao_purity = parameters["cao_purity"]
        self._dehydrator_eff = parameters["dehydrator_eff"]#脱水器传热效率
        self._P_amb = parameters["P_amb"]#环境压力
        self._T_amb = parameters["T_amb"]#环境温度
        self._T_hydr = parameters["T_hydr"]#反应器温度
        self._P_hydr = parameters["P_hydr"] #反应器压力
        self._industrial_waste_heat_t = parameters["industrial_waste_heat_t"]
        self._delta_H_Tref = 104e3  # J/mole CaO504℃反应热
        self._reactants_ES_mole = self._res["dehydrator"]["mole_dehydrator_reactions"]
        self._steam_pressure_loss_ratio = parameters["steam_pressure_loss_ratio"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mcmechanical_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_dehydrator_distance = parameters["storage_dehydrator_distance"]
        self.deltaTmin_SSHX = parameters["deltaTmin_SSHX"] = 25   #固-固换热器最小温差
        self.deltaTmin_SGHX = parameters["deltaTmin_SGHX"] = 20   #固-气换热器最小温差

    def solve(self,input):
        results = {}
        results["BH"] = self._res["BH"]
        results["Dehy"] = self._res
        results["Caes"] = self._res["Caes"]
        del results["Dehy"]["BH"]
        del results["Dehy"]["Caes"]
        
        #Basic input data
        results["Hydr"] = {}
        initialvalue = self.initialvalue(input)
        results["Hydr"]["initialvalue"] = initialvalue 

        #Hydrator
        hydrator = self.hydrator(results["Hydr"]["initialvalue"]["mole_hydrator_reactions"],
                                 results["Dehy"]["dehydrator"])
        results["Hydr"]["hydrator"] = hydrator

        #High temperature section of the heat exchange network
        h_HEN = self.high_tem_HEN(self._industrial_waste_heat_t,
                                                 results["Hydr"]["initialvalue"]["mole_hydrator_reactions"])
        results["Hydr"]["h_HEN"] = h_HEN
        #brayton
        Hydrator_heat = results["Hydr"]["hydrator"]["Q_produce"]-results["Hydr"]["h_HEN"]["heat_need"]
        results["Brayton"] = self._Bray.solve(Hydrator_heat)
        #pc_HEN
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        pc_HEN = self.pc_exchange(results["Brayton"]["heat_recovery"]["t_flue_gas_out"],
                                  results["Brayton"]["evaluation_indicators"]["flue_gas_mass_flow"],
                                  results["Hydr"]["hydrator"]["mass_steam_in"],
                                  flue_gas_name)
        results["Hydr"]["pc_HEN"] = pc_HEN
        #HEN 
        self._T_delta_pinch = 20
        self._pinch_point_data = self.pinch_point_data(results,input["m1"])
        hen_text = self.write_pyPinch_data_text()
        hot_util, cold_util,total_HENA = self.solvehen(hen_text)
        results["Hydr"]["HEN"]={}
        results["Hydr"]["HEN"]["pinch_analysis_text"]=hen_text
        results["Hydr"]["HEN"]["total_HEN_area"]=total_HENA
        results["Hydr"]["HEN"]["hot_utility"] = hot_util
        results["Hydr"]["HEN"]["cold_utility"] = cold_util

        hot_HEN = self.hot_exchange(results["Hydr"]["initialvalue"]["mole_hydrator_reactions"])
        results["Hydr"]["hot_HEN"] = hot_HEN
        #Case all
        m1 = self.Heating_water(results["Hydr"]["pc_HEN"]["heat_out"],
                                results["Hydr"]["hot_HEN"]["heat_out"],
                                results["Brayton"]["evaluation_indicators"]["hot_out"])
        B_case = self.B_case(results,m1)
        results["Case_B"] = B_case

        All_case = self.All_case(results)
        results["Case_All"] = All_case
        return results

    def initialvalue(self,input):
        results = {}
        results["cao_conversion"] = self._cao_conversion
        results["dehydrator_eff"] = self._dehydrator_eff
        results["P_amb"] = self._P_amb
        results["T_amb"] = self._T_amb
        results["T_hydr"]=self._T_hydr
        results["delta_H_Tref"] = self._delta_H_Tref
        results["Store_electrical_power"] = self._Store_electrical_power
        results["Time_scale"] = input["T_X"]
        results["mole_hydrator_reactions"] = self._reactants_ES_mole*input["T_X"]
        return results
    
    def hydrator(self,mole_re,input):
        delta_H_Tr = self._mole_reaction_heat(self._T_hydr, self._P_hydr)
        Q_produce = delta_H_Tr*mole_re

        results={}
        results["is_succeed"] = 1  
        results["delta_H_Tr"] = delta_H_Tr 
        results["Q_produce"] = Q_produce
        results["mole_hydrator_reactions"]=mole_re
        results["mole_Hydration_in"] = mole_re
        results["mole_caoh2_o"] = mole_re/self._cao_conversion
        results["mass_caoh2_o"] = mole_re*M_caoh2
        results["mole_cao_in"] = results["mole_Hydration_in"]
        results["mole_caoh2_in"] = (results["mole_Hydration_in"]/self._cao_conversion)-results["mole_Hydration_in"]
        results["mass_caoh2_in"] = results["mole_caoh2_in"]*M_caoh2
        results["mass_cao_in"] = results["mole_Hydration_in"]*M_cao
        results["mass_caco3_in"] = input["m_camix_in"]-input["mass_caoh2_in"]

        results["mole_steam_in"] = mole_re
        results["mass_steam_in"] = mole_re*M_H2O
        results["m_camix_in"]= results["mass_cao_in"]+results["mass_caoh2_in"]+results["mass_caco3_in"]
        results["m_camix_out"] = results["mass_caoh2_o"]/self._cao_purity
        return results

    def _mole_reaction_heat(self, Thydr, phydr):
        Tref = 25
        cp_cao_mean_Tref_Tr = self._pw.cp0mass_mean("cao", Tref, Thydr)
        cp_caoh2_mean_Tref_Tr = self._pw.cp0mass_mean("caoh2", Thydr, Tref)
        #热容：J/(Kg·K)
        delta_H_Tr = self._delta_H_Tref+((cp_caoh2_mean_Tref_Tr*M_caoh2)*(Tref-Thydr)
                                         + (cp_cao_mean_Tref_Tr*M_cao)*(Thydr-Tref))\
            + (CP.PropsSI('H', 'T', Thydr+273.15, 'P', phydr, "REFPROP::water") -
               CP.PropsSI('H', 'T', 105+273.15, 'P', phydr, "REFPROP::water")+ 
               CP.PropsSI('H', 'T', 95+273.15, 'P', phydr, "REFPROP::water") - 
               CP.PropsSI('H', 'T', Tref+273.15, 'P', phydr, "REFPROP::water"))*M_H2O
        return delta_H_Tr
    
    def high_tem_HEN(self,t,mole):
        cp_cao_cold = self._pw.cp0mass_mean("cao", self._T_hydr-self.deltaTmin_SSHX, t-self.deltaTmin_SSHX)
        cp_caoh2_cold = self._pw.cp0mass_mean("caoh2", self._T_hydr-self.deltaTmin_SSHX, t-self.deltaTmin_SSHX)
        cp_caco3_cold = self._pw.cp0mass_mean("caco3", self._T_hydr-self.deltaTmin_SSHX, t-self.deltaTmin_SSHX)
        cp_caoh2_hot = self._pw.cp0mass_mean("caoh2", self._T_hydr, t)
        cp_caoco3_hot = self._pw.cp0mass_mean("caco3", self._T_hydr, t)
        h_steam_hydr = CP.PropsSI('H', 'T', self._T_hydr+273.15, 'P', self._P_amb, "REFPROP::water")
        h_steam_out = CP.PropsSI('H', 'T', self._T_hydr+273.15-self.deltaTmin_SGHX, 'P', self._P_amb, "REFPROP::water")
        h_steam_in = CP.PropsSI('H', 'T', t+273.15-self.deltaTmin_SGHX, 'P', self._P_amb, "REFPROP::water")
        cp_cao_need = self._pw.cp0mass_mean("cao", self._T_hydr,self._T_hydr-self.deltaTmin_SSHX)
        cp_caoh2_need = self._pw.cp0mass_mean("caoh2", self._T_hydr, self._T_hydr-self.deltaTmin_SSHX)
        cp_caco3_need = self._pw.cp0mass_mean("caco3", self._T_hydr, self._T_hydr-self.deltaTmin_SSHX)

        mass_caoh2_hot = mole/self._cao_conversion*M_caoh2
        mass_caco3_hot = mass_caoh2_hot/self._cao_purity - mass_caoh2_hot
        mass_cao_cold = mole*M_cao
        mass_caoh2_cold = mass_caoh2_hot-mole*M_caoh2
        mass_caco3_cold = mass_caco3_hot

        hot_h = (cp_caoh2_hot*mass_caoh2_hot+cp_caoco3_hot*mass_caco3_hot)*(self._T_hydr-t)
        hot_c = (cp_cao_cold*mass_cao_cold+cp_caoh2_cold*mass_caoh2_cold+cp_caco3_cold*mass_caco3_cold)*(self._T_hydr-t)+(h_steam_out-h_steam_in)*mole*M_H2O
        heat_nes = (cp_cao_need*mass_cao_cold+cp_caoh2_need*mass_caoh2_cold+cp_caco3_need*mass_caco3_cold)*self.deltaTmin_SSHX + (h_steam_hydr-h_steam_out)*mole*M_H2O
        results={} 
        results["hot_h"] = hot_h
        results["hot_c"] = hot_c
        if hot_h*0.96 <= hot_c :
            results["heat_need"] =heat_nes+ (hot_c/0.96)-hot_h
            results["h_lost"]=hot_h*0.04
            a=hot_h*0.96 - (h_steam_out-h_steam_in)*mole*M_H2O
            b=cp_cao_cold*mass_cao_cold+cp_caoh2_cold*mass_caoh2_cold+cp_caco3_cold*mass_caco3_cold
            results["t_cao"] = t-self.deltaTmin_SSHX+a/b
            results["t_steam"] = self._T_hydr-self.deltaTmin_SGHX
        else:
            results["heat_need"] =heat_nes
            results["h_lost"] = 0
            results["t_cao"] = self._T_hydr-self.deltaTmin_SSHX
            results["t_steam"] = self._T_hydr-self.deltaTmin_SGHX
        return results
    
    def pc_exchange(self,t,m1,m2,flue_gas_name):
        h_hsteam_in =  CP.PropsSI('H', 'T', 95+273.15, 'P', self._P_amb, "REFPROP::water")
        h_hsteam_out =  CP.PropsSI('H', 'T', 105+273.15, 'P', self._P_amb, "REFPROP::water")

        fluid=flue_gas_name
        t_flue_gas_in = t
        t_flue_gas_out = 160
        h_flue_gas_inn = CP.PropsSI('H', 'T', self._industrial_waste_heat_t+273.15, 'P', self._P_amb , fluid)
        h_flue_gas_in = CP.PropsSI('H', 'T', t_flue_gas_in+273.15, 'P', self._P_amb , fluid)
        h_flue_gas_out = CP.PropsSI('H', 'T', t_flue_gas_out+273.15, 'P', self._P_amb , fluid)

        hot_h = (h_flue_gas_in-h_flue_gas_out)*m1
        hot_c = (h_hsteam_out-h_hsteam_in)*m2
        results={}
        results["hot_h"] = hot_h
        results["hot_c"] = hot_c
        if hot_h*0.96 >= hot_c :
            results["heat_out"] =hot_h*0.96-hot_c
            results["h_lost"]=hot_h*0.04
            results["flue_gas_mass_flow"] = m1
        else:
            results["heat_out"] = 0
            results["h_lost"] = (hot_c/0.96)*0.04
            results["flue_gas_mass_flow"] = m1+(hot_c/0.96-hot_h)/(h_flue_gas_inn-h_flue_gas_out)
        return results
    
    def pinch_point_data(self,res,m1):
        self._materials=["Ca","Gas","Gas","Ca","Water","Water","Water","Gas"]
        self._HTCs={"Ca":300,"Gas":600,"Water":2500}

        pinch_point_data = {}
        pinch_point_data["TSUPPLY"] = {}
        pinch_point_data["TTARGET"] = {}
        pinch_point_data["ENERGY"] = {}
        pinch_point_data["FLOWRATE"] = {}
        pinch_point_data["CP"] = {}

        ## H1: Camix out: CaOH2
        pinch_point_data["TSUPPLY"]["H_CaM"] = self._T_hydr
        pinch_point_data["TTARGET"]["H_CaM"] = self._T_amb
        pinch_point_data["FLOWRATE"]["H_CaM"] = res["Hydr"]["hydrator"]["m_camix_out"]*0.96
        pinch_point_data["CP"]["H_CaM"] = res["Hydr"]["hydrator"]["m_camix_out"]*0.96*self._pw.cp_camix_mean_Ci(self._T_amb,self._T_hydr,self._cao_purity)
        pinch_point_data["ENERGY"]["H_CaM"] =pinch_point_data["CP"]["H_CaM"]*(self._T_hydr-self._T_amb)
        # H2: flue gas
        pinch_point_data["TSUPPLY"]["C_flue_gas"] = res["Brayton"]["heat_recovery"]["t_flue_gas_out"]
        pinch_point_data["TTARGET"]["C_flue_gas"] = 160
        pinch_point_data["FLOWRATE"]["C_flue_gas"] = res["Brayton"]["evaluation_indicators"]["flue_gas_mass_flow"]*0.96
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["C_flue_gas"] = res["Brayton"]["evaluation_indicators"]["flue_gas_mass_flow"]*0.96 * \
            (CP.PropsSI('H', 'T', res["Brayton"]["heat_recovery"]["t_flue_gas_out"]+273.15,
                        'P', self._P_amb, flue_gas_name) -
             CP.PropsSI('H', 'T', 160+273.15,
                        'P', self._P_amb, flue_gas_name))
        pinch_point_data["CP"]["C_flue_gas"] = pinch_point_data["ENERGY"]["C_flue_gas"] / \
            (res["Brayton"]["heat_recovery"]["t_flue_gas_out"]-160)
        # H3: flue gas
        pinch_point_data["TSUPPLY"]["C2_flue_gas"] =300
        pinch_point_data["TTARGET"]["C2_flue_gas"] = 160
        pinch_point_data["FLOWRATE"]["C2_flue_gas"] =  (res["Hydr"]["pc_HEN"]["flue_gas_mass_flow"]-res["Brayton"]["evaluation_indicators"]["flue_gas_mass_flow"])*0.96
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["C2_flue_gas"] = pinch_point_data["FLOWRATE"]["C2_flue_gas"] * \
            (CP.PropsSI('H', 'T', 300+273.15,
                        'P', self._P_amb, flue_gas_name) -
             CP.PropsSI('H', 'T', 160+273.15,
                        'P', self._P_amb, flue_gas_name))
        pinch_point_data["CP"]["C2_flue_gas"] = pinch_point_data["ENERGY"]["C2_flue_gas"] / \
            (300-160)
        ## C1: Camix out: CaO
        pinch_point_data["TSUPPLY"]["C_CaoM"] = self._T_amb
        pinch_point_data["TTARGET"]["C_CaoM"] = res["Hydr"]["h_HEN"]["t_cao"]
        pinch_point_data["FLOWRATE"]["C_CaoM"] = res["Hydr"]["hydrator"]["m_camix_in"]
        pinch_point_data["CP"]["C_CaoM"] = res["Hydr"]["hydrator"]["m_camix_in"]*self._pw.cp_camix_mean_Co(self._T_amb,res["Hydr"]["h_HEN"]["t_cao"],self._cao_conversion,self._cao_purity)
        pinch_point_data["ENERGY"]["C_CaoM"] =pinch_point_data["CP"]["C_CaoM"]*(res["Hydr"]["h_HEN"]["t_cao"]-self._T_amb)
        ## C2: water
        pinch_point_data["TSUPPLY"]["C_water"] = 60
        pinch_point_data["TTARGET"]["C_water"] = 85
        pinch_point_data["FLOWRATE"]["C_water"] = m1
        pinch_point_data["ENERGY"]["C_water"] = m1 * \
            (CP.PropsSI('H', 'T', 85+273.15,
                        'P', self._P_amb, "REFPROP::water") -
            CP.PropsSI('H', 'T', 60+273.15,
                        'P', self._P_amb, "REFPROP::water"))
        pinch_point_data["CP"]["C_water"] = pinch_point_data["ENERGY"]["C_water"] / \
            (85-60)
        # C3: water before phase change
        pinch_point_data["TSUPPLY"]["be_PCwater"] = self._T_amb
        pinch_point_data["TTARGET"]["be_PCwater"] = 98
        pinch_point_data["FLOWRATE"]["be_PCwater"] = res["Hydr"]["hydrator"]["mass_steam_in"]
        pinch_point_data["ENERGY"]["be_PCwater"] =res["Hydr"]["hydrator"]["mass_steam_in"] * \
            (CP.PropsSI('H', 'T', 98+273.15,
                        'P', self._P_amb, "water") -
             CP.PropsSI('H', 'T', self._T_amb + 273.15,
                        'P', self._P_amb, "water"))
        pinch_point_data["CP"]["be_PCwater"] = pinch_point_data["ENERGY"]["be_PCwater"] / \
            (98-self._T_amb)
        # C4: water ———— phase change
        pinch_point_data["TSUPPLY"]["PC_water"] = 98
        pinch_point_data["TTARGET"]["PC_water"] = 102
        pinch_point_data["FLOWRATE"]["PC_water"] = res["Hydr"]["hydrator"]["mass_steam_in"]
        pinch_point_data["ENERGY"]["PC_water"] =res["Hydr"]["hydrator"]["mass_steam_in"] * \
            (CP.PropsSI('H', 'T', 102+273.15,
                        'P', self._P_amb, "water") -
             CP.PropsSI('H', 'T', 98 + 273.15,
                        'P', self._P_amb, "water"))
        pinch_point_data["CP"]["PC_water"] = pinch_point_data["ENERGY"]["PC_water"] / \
            (102-98)
        # C5: Gas in:Steam
        pinch_point_data["TSUPPLY"]["C_Steam"] = 102
        pinch_point_data["TTARGET"]["C_Steam"] = res["Hydr"]["h_HEN"]["t_steam"]
        pinch_point_data["FLOWRATE"]["C_Steam"] = res["Hydr"]["hydrator"]["mass_steam_in"]
        pinch_point_data["ENERGY"]["C_Steam"] = res["Hydr"]["hydrator"]["mass_steam_in"] * \
            (CP.PropsSI('H', 'T',  res["Hydr"]["h_HEN"]["t_steam"]+273.15,
                        'P', self._P_amb, "water") -
             CP.PropsSI('H', 'T', 102 + 273.15,
                        'P', self._P_amb, "water"))
        pinch_point_data["CP"]["C_Steam"] = pinch_point_data["ENERGY"]["C_Steam"] / \
            (res["Hydr"]["h_HEN"]["t_steam"]-102)
        return pinch_point_data
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
    def B_case(self,results,m1):
        b={}
        b["m_heating_water"] = m1
        b["power_in"]=results["Dehy"]["evaluation_indicators"]["power_in"]
        b["hot_in"]=results["Brayton"]["evaluation_indicators"]["hot_cost"]+results["Hydr"]["pc_HEN"]["hot_h"]
        b["hot_out"]=results["Hydr"]["pc_HEN"]["heat_out"]+results["Hydr"]["hot_HEN"]["heat_out"]+ results["Brayton"]["evaluation_indicators"]["hot_out"]
        b["h_lost"]=results["Brayton"]["evaluation_indicators"]["lost"]+results["Hydr"]["pc_HEN"]["h_lost"]+results["Hydr"]["hot_HEN"]["h_lost"]+b["power_in"]
        b["power"]=results["Brayton"]["evaluation_indicators"]["power"]
        b["cooling_tower"]=results["Brayton"]["cooling_tower"]["hot_cooling_tower"]*results["Brayton"]["evaluation_indicators"]["mass_flow"]
        b["energy_eff"]=(b["hot_out"]+b["power"])/(b["power_in"]+b["hot_in"]+b["cooling_tower"])
        return b
    
    def All_case(self,results):
        a={}
        a["power_in"] = results["Case_B"]["power_in"]+results["Caes"]["power_in"]
        a["hot_in"] = results["Case_B"]["hot_in"]+results["Caes"]["hot_in"]
        a["hot_out"] = results["Case_B"]["hot_out"]+results["Caes"]["hot_out"]
        a["power"] = results["Case_B"]["power"]
        a["energy_eff"] = (a["hot_out"]+a["power"])/(a["power_in"]+a["hot_in"])
        a["exergy_eff"] = (a["hot_out"]*0.15154341+a["power"])/(a["power_in"]+a["hot_in"]*0.414064)
        return a



if __name__ == '__main__':
    parameters = dict()
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["cao_conversion"] = 0.9  #氧化钙转化率
    parameters["cao_purity"] = 0.98 #氢氧化钙含量
    parameters["dehydrator_eff"] = 0.95   #脱水器效率
    parameters["T_amb"] = 20   #环境温度
    parameters["P_amb"] = 101325   #环境压力
    parameters["T_dehy"] = 525     #脱水器温度
    parameters["T_hydr"] = 465     #水合器温度
    parameters["P_hydr"] = 101325
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
    parameters["t_reaction_B"] = 465
    parameters["p_bray_H_B"] = 30000000
    parameters["p_bray_MH_B"] = 16217752.142109105
    parameters["p_bray_ML_B"] = 16217752.142109105
    parameters["p_bray_L_B"] = 7.5e6
    parameters["p_amb"] = 101325

    parameters["mcmechanical_eff"] = 0.98
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100
    calre = Hydrator(parameters)
    inputs={}
    inputs["T_X"] = 1
    inputs["m1"] = 0
    results = calre.solve(inputs)
    print(results)


