import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")


from pyPinch import PyPinch
import pandas as pd
import CoolProp.CoolProp as CP
from Cp0massWrapper import Cp0mass_Wrapper

class Pinch_point_analyzer(object):
    def __init__(self, inputs) -> None:
        self._materials=["Ca","Gas","Gas","Ca","Water"]
        self._HTCs={"Ca":300,"Water":2500,"Gas":600}

        self._m_caco3 = inputs["m_caco3_out"]
        self._m_cao_unr = inputs["m_cao_unr_out"]
        self._m_cao_i = inputs["m_cao_in"]
        self._m_decarbonized_flue_gas = inputs["m_deconbonized_flue_gas_out"]
        self._m_flue_gas_in = inputs["m_flue_gas_in"]
        self._m_water_in = inputs["m_water_in"]

        self._T_carb = inputs["T_carb"]
        self._p_amb = inputs["p_amb"]
        self._T_amb = inputs["T_amb"]
        self._T_flue_gas_fan_out = inputs["T_flue_gas_fan_out"]  
        self._p_flue_gas_fan_out = inputs["p_flue_gas_fan_out"]
        self._T_flue_gas_reactor_in = inputs["T_flue_gas_reactor_in"]
        self._T_cao_reactor_in = inputs["T_cao_reactor_in"]
        self._T_water_reactor_in = inputs["T_water_reactor_in"]
        self._T_decarbonized_flue_gas_out = inputs["T_decarbonized_flue_gas_out"]
        self._p_decarbonized_flue_gas_out = self._p_amb
        self._T_delta_pinch = inputs["T_delta_pinch"]

        self._p_carb_o = inputs["p_carb_o"]
        self._p_carb_i = inputs["p_carb_i"]
        self._p_water = inputs["p_water"]
        self._HTCW=inputs["HTCW"]
        self._HRCP=inputs["HRCP"]

        self._pw = Cp0mass_Wrapper(inputs["flue_gas_composition"], inputs["deconbonized_rate"])
        self._flue_gas_composition = self._pw._norm_flue_gas_composition_s
        self._deconbonized_rate = inputs["deconbonized_rate"]

        pinch_point_data = {}
        pinch_point_data["TSUPPLY"] = {}
        pinch_point_data["TTARGET"] = {}
        pinch_point_data["ENERGY"] = {}
        pinch_point_data["FLOWRATE"] = {}
        pinch_point_data["CP"] = {}

        ## H1: CaCO3+CaO_unr
        pinch_point_data["TSUPPLY"]["H_CaM"] = self._T_carb
        pinch_point_data["TTARGET"]["H_CaM"] = self._T_amb
        pinch_point_data["FLOWRATE"]["H_CaM"] = self._m_caco3+self._m_cao_unr
        pinch_point_data["ENERGY"]["H_CaM"] =\
            (self._m_caco3*self._pw.cp0mass_mean("caco3", self._T_amb, self._T_carb) +
             self._m_cao_unr*self._pw.cp0mass_mean("cao", self._T_amb, self._T_carb)) * \
            (self._T_carb-self._T_amb)
        pinch_point_data["CP"]["H_CaM"] = pinch_point_data["ENERGY"]["H_CaM"] / \
            (self._T_carb-self._T_amb)

        # H2: decarbonized flue gas
        pinch_point_data["TSUPPLY"]["H_flue_gas_decarb"] = self._T_carb
        pinch_point_data["TTARGET"]["H_flue_gas_decarb"] = self._T_decarbonized_flue_gas_out
        pinch_point_data["FLOWRATE"]["H_flue_gas_decarb"] = self._m_decarbonized_flue_gas
        decarb_flue_gas_name = self._pw.get_decarbonized_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["H_flue_gas_decarb"] = self._m_decarbonized_flue_gas * \
            (CP.PropsSI('H', 'T', self._T_carb+273.15,
                        'P', self._p_carb_o, decarb_flue_gas_name) -
             CP.PropsSI('H', 'T', self._T_decarbonized_flue_gas_out + 273.15,
                        'P', self._p_decarbonized_flue_gas_out, decarb_flue_gas_name))
        pinch_point_data["CP"]["H_flue_gas_decarb"] = pinch_point_data["ENERGY"]["H_flue_gas_decarb"] / \
            (self._T_carb-self._T_decarbonized_flue_gas_out)

        # C1: flue gas
        pinch_point_data["TSUPPLY"]["C_flue_gas"] = self._T_flue_gas_fan_out  # todo:T_flue_gas_fan_out
        pinch_point_data["TTARGET"]["C_flue_gas"] = self._T_flue_gas_reactor_in
        pinch_point_data["FLOWRATE"]["C_flue_gas"] = self._m_flue_gas_in
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["C_flue_gas"] = self._m_flue_gas_in * \
            (CP.PropsSI('H', 'T', self._T_flue_gas_reactor_in+273.15,
                        'P', self._p_carb_i, flue_gas_name) -
             CP.PropsSI('H', 'T', self._T_flue_gas_fan_out+273.15,
                        'P', self._p_flue_gas_fan_out, flue_gas_name))
        pinch_point_data["CP"]["C_flue_gas"] = pinch_point_data["ENERGY"]["C_flue_gas"] / \
            (self._T_flue_gas_reactor_in-self._T_flue_gas_fan_out)

        ## C2: CaO
        pinch_point_data["TSUPPLY"]["C_CaO"] = self._T_amb
        pinch_point_data["TTARGET"]["C_CaO"] = self._T_cao_reactor_in
        pinch_point_data["FLOWRATE"]["C_CaO"] = self._m_cao_i
        pinch_point_data["ENERGY"]["C_CaO"] = self._m_cao_i * \
            self._pw.cp0mass_mean("cao", self._T_amb, self._T_cao_reactor_in) * \
            (self._T_cao_reactor_in-self._T_amb)
        pinch_point_data["CP"]["C_CaO"] = pinch_point_data["ENERGY"]["C_CaO"] / \
            (self._T_cao_reactor_in-self._T_amb)

        ## C3: water
        if self._HRCP==1:
            pinch_point_data["TSUPPLY"]["C_water"] = self._T_amb
            pinch_point_data["TTARGET"]["C_water"] = self._T_water_reactor_in
            pinch_point_data["FLOWRATE"]["C_water"] = self._m_water_in
            pinch_point_data["ENERGY"]["C_water"] = self._m_water_in * \
                (CP.PropsSI('H', 'T', self._T_water_reactor_in+273.15,
                            'P', self._p_water, "REFPROP::water") -
                CP.PropsSI('H', 'T', self._T_amb+273.15,
                            'P', self._p_water, "REFPROP::water"))
            pinch_point_data["CP"]["C_water"] = pinch_point_data["ENERGY"]["C_water"] / \
                (self._T_water_reactor_in-self._T_amb)

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
