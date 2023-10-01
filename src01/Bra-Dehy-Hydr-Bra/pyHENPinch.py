import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")

from pyPinch import PyPinch
import pandas as pd
import CoolProp.CoolProp as CP
from Cp0massWrapper import Cp0mass_Wrapper


class Hen_pinch_analyzer(object):
    def __init__(self, inputs) -> None:

        self._X = inputs["cao_conversion"]
        self._Y = inputs["cao_purity"]
        self._materials=["Ca","Gas","Water","Ca","Gas","Water"]
        self._HTCs={"Ca":300,"Gas":600,"Water":2500}

        self._m_camix_out_dehy = inputs["m_camix_out"]
        self._m_camix_in_dehy = inputs["m_camix_in"]
        self._m_steam = inputs["m_steam_out"]
        self._m_flue_gas = inputs["m_flue_gas_in"]
        self._m_water_in = inputs["m_water_in"] #非定值

        self._T_dehy = inputs["T_dehy"]
        self._p_amb = inputs["p_amb"]
        self._T_amb = inputs["T_amb"]
        self._T_flue_gas_bray_out = inputs["T_flue_gas_bray_out"]  
        self._p_flue_gas_bray_out = inputs["p_flue_gas_bray_out"]
        self._T_flue_gas_dew = inputs["T_flue_gas_dew"]
        self._T_solid_in = inputs["T_solid_in"] #非定值
        self._T_water_supply_in= inputs["T_water_supply_in"]
        self._T_water_reactor_in = inputs["T_water_reactor_in"]
        self._p_steam_out = self._p_amb
        self._t_lsteam_in = inputs["T_lsteam_in"]
        self._t_lsteam_out = inputs["T_lsteam_out"]

        self._T_delta_pinch = inputs["T_delta_pinch"]

        self._p_dehy_o = inputs["p_dehy_o"]#
        self._p_dehy_i = inputs["p_dehy_i"]
        self._p_water = inputs["p_water_after_pump"]

        self._pw = Cp0mass_Wrapper(inputs["flue_gas_composition"])#
        self._flue_gas_composition = self._pw._norm_flue_gas_composition

        pinch_point_data = {}
        pinch_point_data["TSUPPLY"] = {}
        pinch_point_data["TTARGET"] = {}
        pinch_point_data["ENERGY"] = {}
        pinch_point_data["FLOWRATE"] = {}
        pinch_point_data["CP"] = {}

        ## H1: Camix out: CaO
        pinch_point_data["TSUPPLY"]["H_CaM"] = self._T_flue_gas_bray_out+5
        pinch_point_data["TTARGET"]["H_CaM"] = self._T_amb
        pinch_point_data["FLOWRATE"]["H_CaM"] = self._m_camix_out_dehy
        pinch_point_data["CP"]["H_CaM"] = self._m_camix_out_dehy*self._pw.cp_camix_mean_Co(self._T_amb,self._T_flue_gas_bray_out+5,self._X,self._Y)
        pinch_point_data["ENERGY"]["H_CaM"] =pinch_point_data["CP"]["H_CaM"]*(self._T_flue_gas_bray_out+5-self._T_amb)
        

        # H2: Gas out:Steam
        pinch_point_data["TSUPPLY"]["H_Steam"] = self._T_flue_gas_bray_out
        pinch_point_data["TTARGET"]["H_Steam"] = self._t_lsteam_in
        pinch_point_data["FLOWRATE"]["H_Steam"] = self._m_steam
        pinch_point_data["ENERGY"]["H_Steam"] = self._m_steam * \
            (CP.PropsSI('H', 'T', self._T_flue_gas_bray_out+273.15,
                        'P', self._p_dehy_o, "water") -
             CP.PropsSI('H', 'T', self._t_lsteam_in + 273.15,
                        'P', self._p_steam_out, "water"))
        pinch_point_data["CP"]["H_Steam"] = pinch_point_data["ENERGY"]["H_Steam"] / \
            (self._T_flue_gas_bray_out-self._t_lsteam_in)

        # H3: flue gas
        pinch_point_data["TSUPPLY"]["C_flue_gas"] = self._T_flue_gas_bray_out  #
        pinch_point_data["TTARGET"]["C_flue_gas"] = self._T_flue_gas_dew
        pinch_point_data["FLOWRATE"]["C_flue_gas"] = self._m_flue_gas
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["C_flue_gas"] = self._m_flue_gas * \
            (CP.PropsSI('H', 'T', self._T_flue_gas_bray_out+273.15,
                        'P', self._p_amb, flue_gas_name) -
             CP.PropsSI('H', 'T', self._T_flue_gas_dew+273.15,
                        'P', self._p_amb, flue_gas_name))
        pinch_point_data["CP"]["C_flue_gas"] = pinch_point_data["ENERGY"]["C_flue_gas"] / \
            (self._T_flue_gas_bray_out-self._T_flue_gas_dew)
        
        # H4: water after phase change
        pinch_point_data["TSUPPLY"]["PC_water"] = self._t_lsteam_out
        pinch_point_data["TTARGET"]["PC_water"] = self._T_amb
        pinch_point_data["FLOWRATE"]["PC_water"] = self._m_steam
        pinch_point_data["ENERGY"]["PC_water"] = self._m_steam * \
            (CP.PropsSI('H', 'T', self._t_lsteam_out+273.15,
                        'P', self._p_dehy_o, "water") -
             CP.PropsSI('H', 'T', self._T_amb + 273.15,
                        'P', self._p_steam_out, "water"))
        pinch_point_data["CP"]["PC_water"] = pinch_point_data["ENERGY"]["PC_water"] / \
            (self._t_lsteam_out-self._T_amb)

        ## C1: Camix out: Ca(OH)2 
        pinch_point_data["TSUPPLY"]["C_Caoh2"] = self._T_amb
        pinch_point_data["TTARGET"]["C_Caoh2"] = self._T_solid_in
        pinch_point_data["FLOWRATE"]["C_Caoh2"] = self._m_camix_in_dehy
        pinch_point_data["CP"]["C_Caoh2"] = self._m_camix_in_dehy*self._pw.cp_camix_mean_Ci(self._T_amb,self._T_solid_in,self._Y)
        pinch_point_data["ENERGY"]["C_Caoh2"] =pinch_point_data["CP"]["C_Caoh2"]*(self._T_solid_in-self._T_amb)

        ## C2: water
        pinch_point_data["TSUPPLY"]["C_water"] = self._T_water_supply_in
        pinch_point_data["TTARGET"]["C_water"] = self._T_water_reactor_in
        pinch_point_data["FLOWRATE"]["C_water"] = self._m_water_in
        pinch_point_data["ENERGY"]["C_water"] = self._m_water_in * \
            (CP.PropsSI('H', 'T', self._T_water_reactor_in+273.15,
                        'P', self._p_water, "REFPROP::water") -
            CP.PropsSI('H', 'T', self._T_water_supply_in+273.15,
                        'P', self._p_water, "REFPROP::water"))
        pinch_point_data["CP"]["C_water"] = pinch_point_data["ENERGY"]["C_water"] / \
            (self._T_water_reactor_in-self._T_water_supply_in)
       
        self._pinch_point_data = pinch_point_data

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

    def _myRound(self,a,ndigital=3):
        return round(a,ndigital)

    def solve(self, input_text):
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
