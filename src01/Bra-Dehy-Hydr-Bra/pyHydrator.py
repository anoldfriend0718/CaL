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
M_caco3 = 100e-3
C_cao = 112396 #J/mole
C_caoh2 = 54593 
C_h2o = 750

class Hydrator(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
        self._de = Dehydrator(parameters)
        self._Bray = Brayton(parameters)
        self._Store_electrical_power = parameters["Store_electrical_power"] #1MW

        self._cao_conversion = parameters["cao_conversion"]#氢氧化钙转化率
        self._cao_purity = parameters["cao_purity"]
        self._dehydrator_eff = parameters["dehydrator_eff"]#脱水器传热效率
        self._hydrator_eff = parameters["hydrator_eff"]#水合器传热效率
        self._P_amb = parameters["p_amb"]#环境压力
        self._T_amb = parameters["t_amb"]#环境温度    
        self._industrial_waste_heat_t = parameters["industrial_waste_heat_t"]
        self._delta_H_Tref = 104e3  # J/mole 反应热

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
        self._Hydr_ot = inputs["Hydr_overheating_temperature"]
        self.Hydr_cao_in  =inputs["Hydr_cao_in"]
        self.Hydr_steam_in = self.Hydr_cao_in+5
        self._m1=inputs["cn_m1"]   
        self._t= inputs["T_X"]
        t_equilibrium = self.equilibrium()
        self._T_hydr = t_equilibrium-self._Hydr_ot
        self._reactants_ES_mole = self._res["dehydrator"]["mole_dehydrator_reactions"]
        results = {}
        results["BH"] = self._res["BH"]
        results["Dehy"] = self._res
        results["Caes"] = self._res["Caes"]
        del results["Dehy"]["BH"]
        del results["Dehy"]["Caes"]    
        #Basic input data
        results["Hydr"] = {}
        initialvalue = self.initialvalue(inputs)
        results["Hydr"]["initialvalue"] = initialvalue 
        #Hydrator
        hydrator = self.hydrator(results["Hydr"]["initialvalue"]["mole_hydrator_reactions"],
                                 results["Dehy"]["dehydrator"])
        results["Hydr"]["hydrator"] = hydrator
        #brayton
        self._Hydrator_heat = results["Hydr"]["hydrator"]["Q_produce"]
        results["Brayton"] = self._Bray.solve(self._Hydrator_heat,inputs)
        #反应器补充
        results["Hydr"]["hydrator"]["exergy"]["hot_out"]=results["Brayton"]["evaluation_indicators"]["exergy"]["re_Heat_in"]
        results["Hydr"]["hydrator"]["exergy"]["lost"] = (results["Hydr"]["hydrator"]["exergy"]["chemical_t"]
                                                         -results["Hydr"]["hydrator"]["exergy"]["hot_out"]-
                                                         results["Hydr"]["hydrator"]["exergy"]["sensible_heat"])

        results["Hydr"]["flue_gas_mass_flow"]=results["Brayton"]["evaluation_indicators"]["flue_gas_mass_flow"]+inputs["m2"]
        #HEN 
        self._T_delta_pinch = 20
        self._pinch_point_data,hot_in,hot_out,hot_hen,cold_hen,exergy_hen=self.pinch_point_data(results,inputs["m1"],inputs["m2"])
        hen_text = self.write_pyPinch_data_text()
        hot_util, cold_util,total_HENA = self.solvehen(hen_text)
        results["Hydr"]["HEN"]={}
        results["Hydr"]["HEN"]["pinch_analysis_text"]=hen_text
        results["Hydr"]["HEN"]["total_HEN_area"]=total_HENA
        results["Hydr"]["HEN"]["hot_utility"] = hot_util
        results["Hydr"]["HEN"]["cold_utility"] = cold_util
        results["Hydr"]["HEN"]["p"]=self._pinch_point_data
        results["Hydr"]["HEN"]["hot_in"] =hot_in
        results["Hydr"]["HEN"]["hot_out"] =hot_out
        results["Hydr"]["HEN"]["hot"] =hot_hen
        results["Hydr"]["HEN"]["cold"] =cold_hen
        results["Hydr"]["HEN"]["exergy"] =exergy_hen
        #Case all
        B_case = self.B_case(results,inputs["m1"],hot_in , hot_out,hot_hen,cold_hen)
        results["Case_B"] = B_case
        All_case = self.All_case(results)
        results["Case_All"] = All_case
        return results  
    def equilibrium(self):
        p = self._P_hydr
        t = (-12845/((math.log(p/1e5))-16.508))-273.15
        return t
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
        results["m1"]=input["m1"]
        results["m2"]=input["m2"]
        return results 
    def hydrator(self,mole_re,input):
        delta_H_Tr = self._mole_reaction_heat(self._T_hydr, self._P_hydr)
        Q_produce = delta_H_Tr*mole_re

        results={}
        results["is_succeed"] = 1  
        results["delta_H_Tr"] = delta_H_Tr 
        results["Q_reactions"] = Q_produce
        results["mole_hydrator_reactions"]=mole_re
        results["mole_Hydration_in"] = mole_re
        results["mole_caoh2_o"] = mole_re/self._cao_conversion
        results["mass_caoh2_o"] = mole_re/self._cao_conversion*M_caoh2
        results["mole_cao_in"] = results["mole_Hydration_in"]
        results["mass_cao_in"] = results["mole_Hydration_in"]*M_cao
        results["mole_caoh2_in"] = (results["mole_Hydration_in"]/self._cao_conversion)-results["mole_Hydration_in"]
        results["mass_caoh2_in"] = results["mole_caoh2_in"]*M_caoh2       
        results["mass_caco3_in"] = results["mass_caoh2_o"]/self._cao_purity-results["mass_caoh2_o"]
        results["mole_caco3_in"] = results["mass_caco3_in"]/M_caco3
        results["mole_steam_in"] = mole_re
        results["mass_steam_in"] = mole_re*M_H2O
        results["m_camix_in"]= results["mass_cao_in"]+results["mass_caoh2_in"]+results["mass_caco3_in"]
        results["m_camix_out"] = results["mass_caoh2_o"]/self._cao_purity

        heat_cao = (self._pw.cp_camix_mean_Co(self.Hydr_cao_in,self._T_hydr,self._cao_conversion,
                                              self._cao_purity)*(self._T_hydr-self.Hydr_cao_in))*results["m_camix_in"]
        heat_steam=((CP.PropsSI('H', 'T', self._T_hydr+273.15, 'P', self._P_hydr, "REFPROP::water")
                    -CP.PropsSI('H', 'T', self.Hydr_steam_in+273.15, 'P', self._P_hydr, "REFPROP::water"))
                    *results["mass_steam_in"])
        results["Q_produce"] = Q_produce -heat_cao- heat_steam
        results["exergy"] = {}
        
        results["exergy"]["chemical_t"]=self._Bray.ex_calculations1( self._T_hydr,Q_produce)
        results["exergy"]["chemical"] = (112.396-54.593)*1000*mole_re
        results["exergy"]["sensible_heat"] = (self._Bray.ex_calculations(self.Hydr_cao_in, self._T_hydr,heat_cao)
                                              +(self.E_steam(self._T_hydr,self._P_hydr)-
                                                self.E_steam(self.Hydr_steam_in,self._P_hydr))
                                                *results["mass_steam_in"])

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
    def E_steam(self,T,P):
        H1=CP.PropsSI('H', 'T', T+273.15, 'P', P, "REFPROP::water")
        S1=CP.PropsSI('S', 'T', T+273.15, 'P', P, "REFPROP::water")
        H0=CP.PropsSI('H', 'T', self._T_amb+273.15, 'P', self._P_amb, "REFPROP::water")
        S0=CP.PropsSI('S', 'T', self._T_amb+273.15, 'P', self._P_amb, "REFPROP::water")
        a=H1-H0-(self._T_amb+273.15)*(S1-S0)
        return a 
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
    
    def pinch_point_data(self,res,m1,m2):
        self._materials=["Ca","Gas","Gas","Ca","Water","Water","Water","Gas"]
        self._HTCs={"Ca":300,"Gas":600,"Water":2500}

        pinch_point_data = {}
        pinch_point_data["TSUPPLY"] = {}
        pinch_point_data["TTARGET"] = {}
        pinch_point_data["ENERGY"] = {}
        pinch_point_data["FLOWRATE"] = {}
        pinch_point_data["CP"] = {}

        ## H1: Camix out: CaOH2
        pinch_point_data["TSUPPLY"]["C_Caoh2M"] = self._T_hydr
        pinch_point_data["TTARGET"]["C_Caoh2M"] = self._T_amb
        pinch_point_data["FLOWRATE"]["C_Caoh2M"] = res["Hydr"]["hydrator"]["m_camix_out"]*0.96
        pinch_point_data["CP"]["C_Caoh2M"] =(pinch_point_data["FLOWRATE"]["C_Caoh2M"]*
                                          self._pw.cp_camix_mean_Ci(self._T_amb,self._T_hydr,self._cao_purity))
        pinch_point_data["ENERGY"]["C_Caoh2M"] =pinch_point_data["CP"]["C_Caoh2M"]*(self._T_hydr-self._T_amb)
        # H2: flue gas
        pinch_point_data["TSUPPLY"]["C_flue_gas"] = res["Brayton"]["heat_recovery"]["t_flue_gas_out"]
        pinch_point_data["TTARGET"]["C_flue_gas"] = 160
        pinch_point_data["FLOWRATE"]["C_flue_gas"] = res["Brayton"]["evaluation_indicators"]["flue_gas_mass_flow"]*0.96
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["C_flue_gas"] = pinch_point_data["FLOWRATE"]["C_flue_gas"] * \
            (CP.PropsSI('H', 'T', res["Brayton"]["heat_recovery"]["t_flue_gas_out"]+273.15,
                        'P', self._P_amb, flue_gas_name) -
             CP.PropsSI('H', 'T', 160+273.15,
                        'P', self._P_amb, flue_gas_name))
        if res["Brayton"]["heat_recovery"]["t_flue_gas_out"]==160:
            pinch_point_data["CP"]["C_flue_gas"] = 0
        else:
            pinch_point_data["CP"]["C_flue_gas"] = pinch_point_data["ENERGY"]["C_flue_gas"] / \
            (res["Brayton"]["heat_recovery"]["t_flue_gas_out"]-160)
      
        # H3: flue gas
        pinch_point_data["TSUPPLY"]["C2_flue_gas"] = self._industrial_waste_heat_t
        pinch_point_data["TTARGET"]["C2_flue_gas"] = 160
        pinch_point_data["FLOWRATE"]["C2_flue_gas"] = m2*0.96
        pinch_point_data["ENERGY"]["C2_flue_gas"] = pinch_point_data["FLOWRATE"]["C2_flue_gas"] * \
            (CP.PropsSI('H', 'T', self._industrial_waste_heat_t+273.15,
                        'P', self._P_amb, flue_gas_name) -
             CP.PropsSI('H', 'T', 160+273.15,
                        'P', self._P_amb, flue_gas_name))
        pinch_point_data["CP"]["C2_flue_gas"] = pinch_point_data["ENERGY"]["C2_flue_gas"] / \
            (self._industrial_waste_heat_t-160)
        
        ## C1: Camix out: CaO
        pinch_point_data["TSUPPLY"]["C_CaoM"] = self._T_amb
        pinch_point_data["TTARGET"]["C_CaoM"] = self.Hydr_cao_in
        pinch_point_data["FLOWRATE"]["C_CaoM"] = res["Hydr"]["hydrator"]["m_camix_in"]
        pinch_point_data["CP"]["C_CaoM"] =(res["Hydr"]["hydrator"]["m_camix_in"]*
                                           self._pw.cp_camix_mean_Co(self._T_amb,self.Hydr_cao_in,self._cao_conversion,self._cao_purity))
        pinch_point_data["ENERGY"]["C_CaoM"] =pinch_point_data["CP"]["C_CaoM"]*(self.Hydr_cao_in-self._T_amb)
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
        pinch_point_data["TTARGET"]["C_Steam"] = self.Hydr_steam_in
        pinch_point_data["FLOWRATE"]["C_Steam"] = res["Hydr"]["hydrator"]["mass_steam_in"]
        pinch_point_data["ENERGY"]["C_Steam"] = res["Hydr"]["hydrator"]["mass_steam_in"] * \
            (CP.PropsSI('H', 'T', self.Hydr_steam_in+273.15,
                        'P', self._P_amb, "water") -
             CP.PropsSI('H', 'T', 102 + 273.15,
                        'P', self._P_amb, "water"))
        pinch_point_data["CP"]["C_Steam"] = pinch_point_data["ENERGY"]["C_Steam"] / \
            (self.Hydr_steam_in-102)
        
        hot = (pinch_point_data["ENERGY"]["C_flue_gas"] + pinch_point_data["ENERGY"]["C2_flue_gas"])/0.96
        a =(pinch_point_data["ENERGY"]["C_Caoh2M"] + 
            pinch_point_data["ENERGY"]["C_flue_gas"] + 
            pinch_point_data["ENERGY"]["C2_flue_gas"])/0.96
        b = (pinch_point_data["ENERGY"]["C_CaoM"] + 
             pinch_point_data["ENERGY"]["C_water"]+
             pinch_point_data["ENERGY"]["be_PCwater"]+
             pinch_point_data["ENERGY"]["PC_water"]+
             pinch_point_data["ENERGY"]["C_Steam"])
        hot_out = pinch_point_data["ENERGY"]["C_water"]

        exergy={}
        exergy["C_Caoh2M"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["C_Caoh2M"],
                                             pinch_point_data["TTARGET"]["C_Caoh2M"],
                                             pinch_point_data["ENERGY"]["C_Caoh2M"])
        exergy["C2_flue_gas"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["C2_flue_gas"],
                                             pinch_point_data["TTARGET"]["C2_flue_gas"],
                                             pinch_point_data["ENERGY"]["C2_flue_gas"])
        exergy["C_flue_gas"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["C_flue_gas"],
                                             pinch_point_data["TTARGET"]["C_flue_gas"],
                                             pinch_point_data["ENERGY"]["C_flue_gas"])
        
        exergy["PC_water"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["PC_water"],
                                             pinch_point_data["TTARGET"]["PC_water"],
                                             pinch_point_data["ENERGY"]["PC_water"])
        exergy["C_CaoM"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["C_CaoM"],
                                             pinch_point_data["TTARGET"]["C_CaoM"],
                                             pinch_point_data["ENERGY"]["C_CaoM"])
        exergy["be_PCwater"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["be_PCwater"],
                                             pinch_point_data["TTARGET"]["be_PCwater"],
                                             pinch_point_data["ENERGY"]["be_PCwater"])
        exergy["C_water"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["C_water"],
                                             pinch_point_data["TTARGET"]["C_water"],
                                             pinch_point_data["ENERGY"]["C_water"])
        exergy["C_Steam"]=self.ex_calculations(pinch_point_data["TSUPPLY"]["C_Steam"],
                                             pinch_point_data["TTARGET"]["C_Steam"],
                                             pinch_point_data["ENERGY"]["C_Steam"])
        exergy["hot_in"]=(exergy["C_Caoh2M"]+exergy["C2_flue_gas"]
                          +exergy["C_flue_gas"])/0.96
        exergy["cold_add"]=(exergy["PC_water"]+exergy["C_CaoM"]+exergy["C_water"]+
                            exergy["be_PCwater"]+exergy["C_Steam"])
        exergy["lost"]=exergy["hot_in"]-exergy["cold_add"]

        return pinch_point_data,hot,hot_out,a,b,exergy
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
    def B_case(self,results,m1,hot_in,hot_out,hot_hen,cold_hen):
        b={}
        b["m_heating_water"] = m1
        b["power_in"]=results["Dehy"]["evaluation_indicators"]["power_in"]*self._t
        b["hot_in"]=results["Brayton"]["evaluation_indicators"]["hot_cost_gas"]+hot_in
        b["hot_in_hydr"] = hot_in
        b["hot_out"]=hot_out
        b["h_lost_hydr"] =hot_in+results["Caes"]["hot_stockpile"]-self._Hydrator_heat
        b["h_lost"]=results["Brayton"]["evaluation_indicators"]["lost_all"]+b["h_lost_hydr"]
        b["power"]=results["Brayton"]["evaluation_indicators"]["power"]

        b["energy_eff_hydr"]=((b["hot_out"]+results["Brayton"]["evaluation_indicators"]["re_Heat_in"])
                              /(b["power_in"]+b["hot_in_hydr"]+results["Caes"]["hot_stockpile"]))
        b["energy_eff"]=(b["hot_out"]+b["power"]-b["power_in"])/(b["hot_in"]+results["Caes"]["hot_stockpile"])
        return b
    
    def All_case(self,results):
        a={}
        a["power_in"] =results["Caes"]["power_in"]
        a["hot_in"] = results["Case_B"]["hot_in"]+results["Caes"]["hot_in"]
        a["hot_out"] = results["Caes"]["hot_out"]
        a["power"] = results["Case_B"]["power"]-results["Case_B"]["power_in"]
        a["Round-trip_eff"] = a["power"]/a["power_in"]
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

    parameters["cao_conversion"] = 0.95  #氧化钙转化率
    parameters["cao_purity"] = 0.98 #氢氧化钙含量
    parameters["dehydrator_eff"] = 0.95   #脱水器效率
    parameters["steam_pressure_loss_ratio"] = 0.01
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_dehydrator_distance"] = 100

    parameters["hydrator_eff"] = 0.95   #水合器器效率
    parameters["p_bray_L_B"] = 7.5e6

    inputs={}
    inputs["p_bray_H"] = 19447839.26865841#优化变量1，热泵循环最高压力
    inputs["p_bray_M"] = 12827110.4341202 #优化变量2，热泵循环中间压力
    inputs["p_Dehy"] = 1e5 #变量4，反应器压力
    inputs["Dehy_caoh2_in"]=444.58
    inputs["Dehy_overheating_temperature"] = 20 #变量2，脱水反应器过热温度

    inputs["p_bray_H_B"] = 30e6
    inputs["p_bray_MH_B"] = 16217752.142109105
    inputs["p_bray_ML_B"] = 12217752.142109105
    inputs["p_Hydr"] = 1e5
    inputs["Hydr_overheating_temperature"] = 40

    inputs["Hydr_cao_in"]=440

    inputs["cn_m1"] = 12.34
    inputs["m2"] = 1.532#补充烟气
    inputs["T_X"] = 1
    inputs["m1"] = 0#供暖水
    calre = Hydrator(parameters)
    results = calre.solve(inputs)
    print(results)
    print(results["Hydr"]["HEN"]["pinch_analysis_text"])


