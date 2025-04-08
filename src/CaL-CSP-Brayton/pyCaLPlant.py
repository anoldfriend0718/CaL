import os
import sys
CaLRepo = os.environ.get("CaLRepo")
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")

import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
from pyPinch import PyPinch
from scipy import interpolate



class Cp0mass_Wrapper(object):
    def __init__(self) -> None:
        self.f_cao, self.f_caco3 = self._read_solid_properties()
        
    # T: oC
    # cp: J/kg/K
    def cp0mass(self,material,T,p=1e5):
        material = material.lower()
        if material=="co2":
            cp=self._cp0mass_co2(T,p)
        elif material=="cao":
            cp=self._cp0mass_cao(T)
        elif material=="caco3":
            cp=self._cp0mass_caco3(T)
        else:
            raise ValueError(f"Cp0mass of {material} is not supported")
        return cp

    # T: oC
    # cp: J/kg/K
    def cp0mass_mean(self, material, Ti, To, p=1e5,interval=10):
        Ts = self._arange(Ti, To, interval)
        cps = []
        for T in Ts:
            cp = self.cp0mass(material,T,p)
            cps.append(cp)
        cp_mean = np.mean(cps)
        return cp_mean
    
    def _read_solid_properties(self):
        data_csv = f"{CaLRepo}/data/cpmass_cao_caco3.csv"
        df = pd.read_csv(data_csv)
        f_cao = interpolate.interp1d(df["TEMP"], df["CAO"])
        f_caco3 = interpolate.interp1d(df["TEMP"], df["CACO3"])
        return f_cao, f_caco3

    def _cp0mass_cao(self,T):
        cp = float(self.f_cao(T))
        return cp

    def _cp0mass_caco3(self,T):
        cp = float(self.f_caco3(T))
        return cp

    def _cp0mass_co2(self,T, p):
        fluid = 'REFPROP::co2'
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp

    def _arange(self,start, stop, step=1, endpoint=True):
        arr = np.arange(start, stop, step)

        if endpoint and arr[-1] != stop:
            arr = np.concatenate([arr, [stop]])

        return arr
    

class Pinch_point_analyzer(object):
    def __init__(self, input) -> None:
        self._m_caco3 = input["m_caco3"]
        self._m_cao_unr = input["m_cao_unr"]
        self._m_co2_rec = input["m_co2_rec"]
        self._m_cao_i = input["m_cao_i"]
        self._m_co2_mix = input["m_co2_mix"]
        self._m_co2_stoic = input["m_co2_stoic"]
        self._T_reaction = input["T_reaction"]
        self._T_atm = input["T_atm"]
        self._T_co2_main_turbine_out = input["T_co2_main_turbine_out"]
        self._T_co2_compressor_in = input["T_co2_compressor_in"]
        self._T_co2_mix = input["T_co2_mix"]
        self._T_co2_reactor_in = input["T_co2_reactor_in"]
        self._T_co2_storage_turbine_in = input["T_co2_storage_turbine_in"]
        self._T_cao_reactor_in = input["T_cao_reactor_in"]
        self._p_reaction = input["p_reaction"]
        self._p_co2_main_turbine_out = input["p_co2_main_turbine_out"]
        self._p_co2_compressor_in = input["p_co2_compressor_in"]
        self._p_co2_compressor_out = input["p_co2_compressor_out"]
        self._p_co2_storage = input["p_co2_storage"]
        self._p_co2_storage_turbine_in = input["p_co2_storage_turbine_in"]
        self._T_min = input["T_min"]

        self._pw = Cp0mass_Wrapper()
        pinch_point_data = {}
        pinch_point_data["TSUPPLY"] = {}
        pinch_point_data["TTARGET"] = {}
        pinch_point_data["ENERGY"] = {}
        pinch_point_data["FLOWRATE"] = {}
        pinch_point_data["CP"] = {}
        ## H_CaM
        pinch_point_data["TSUPPLY"]["H_CaM"] = self._T_reaction
        pinch_point_data["TTARGET"]["H_CaM"] = self._T_atm
        pinch_point_data["FLOWRATE"]["H_CaM"] = self._m_caco3+self._m_cao_unr
        pinch_point_data["ENERGY"]["H_CaM"] =\
            (self._m_caco3*self._pw.cp0mass_mean("caco3", self._T_atm, self._T_reaction) +
             self._m_cao_unr*self._pw.cp0mass_mean("cao", self._T_atm, self._T_reaction)) * \
            (self._T_reaction-self._T_atm)
        pinch_point_data["CP"]["H_CaM"] = pinch_point_data["ENERGY"]["H_CaM"] / \
            (self._T_reaction-self._T_atm)
        ## H_CO2_reac
        pinch_point_data["TSUPPLY"]["H_CO2_reac"] = self._T_co2_main_turbine_out
        pinch_point_data["TTARGET"]["H_CO2_reac"] = self._T_co2_compressor_in
        pinch_point_data["FLOWRATE"]["H_CO2_reac"] = self._m_co2_rec
        pinch_point_data["ENERGY"]["H_CO2_reac"] = self._m_co2_rec * \
            (CP.PropsSI('H', 'T', self._T_co2_main_turbine_out+273.15,
                        'P', self._p_co2_main_turbine_out, "REFPROP::co2") -
             CP.PropsSI('H', 'T', self._T_co2_compressor_in + 273.15,
                        'P', self._p_co2_compressor_in, "REFPROP::co2"))
        pinch_point_data["CP"]["H_CO2_reac"] = pinch_point_data["ENERGY"]["H_CO2_reac"] / \
            (self._T_co2_main_turbine_out-self._T_co2_compressor_in)
        ## C_CaO
        pinch_point_data["TSUPPLY"]["C_CaO"] = self._T_atm
        pinch_point_data["TTARGET"]["C_CaO"] = self._T_cao_reactor_in
        pinch_point_data["FLOWRATE"]["C_CaO"] = self._m_cao_i
        pinch_point_data["ENERGY"]["C_CaO"] = self._m_cao_i * \
            self._pw.cp0mass_mean("cao", self._T_atm, self._T_cao_reactor_in) * \
            (self._T_cao_reactor_in-self._T_atm)
        pinch_point_data["CP"]["C_CaO"] = pinch_point_data["ENERGY"]["C_CaO"] / \
            (self._T_cao_reactor_in-self._T_atm)

        ## C_CO2_mix
        pinch_point_data["TSUPPLY"]["C_CO2_mix"] = self._T_co2_mix
        pinch_point_data["TTARGET"]["C_CO2_mix"] = self._T_co2_reactor_in
        pinch_point_data["FLOWRATE"]["C_CO2_mix"] = self._m_co2_mix
        pinch_point_data["ENERGY"]["C_CO2_mix"] = self._m_co2_mix * \
            (CP.PropsSI('H', 'T', self._T_co2_reactor_in+273.15,
                        'P', self._p_reaction, "REFPROP::co2") -
             CP.PropsSI('H', 'T', self._T_co2_mix+273.15,
                        'P', self._p_co2_compressor_out, "REFPROP::co2"))
        pinch_point_data["CP"]["C_CO2_mix"] = pinch_point_data["ENERGY"]["C_CO2_mix"] / \
            (self._T_co2_reactor_in-self._T_co2_mix)

        ## C_CO2_stoic
        pinch_point_data["TSUPPLY"]["C_CO2_stoic"] = self._T_atm
        pinch_point_data["TTARGET"]["C_CO2_stoic"] = self._T_co2_storage_turbine_in
        pinch_point_data["FLOWRATE"]["C_CO2_stoic"] = self._m_co2_stoic
        pinch_point_data["ENERGY"]["C_CO2_stoic"] = self._m_co2_stoic * \
            (CP.PropsSI('H', 'T', self._T_co2_storage_turbine_in+273.15,
                        'P', self._p_co2_storage_turbine_in, "REFPROP::co2") -
             CP.PropsSI('H', 'T', self._T_atm+273.15,
                        'P', self._p_co2_storage, "REFPROP::co2"))
        pinch_point_data["CP"]["C_CO2_stoic"] = pinch_point_data["ENERGY"]["C_CO2_stoic"] / \
            (self._T_co2_storage_turbine_in-self._T_atm)

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
            lines.insert(0, f'Tmin, {self._T_min},\n')
            fp.seek(0)
            fp.writelines(lines)

    def write_pyPinch_data_text(self):
        data = {}
        data["CP"] = self._pinch_point_data["CP"]
        data["TSUPPLY"] = self._pinch_point_data["TSUPPLY"]
        data["TTARGET"] = self._pinch_point_data["TTARGET"]
        df = pd.DataFrame(data)
        df["CP"] = df["CP"]/1000
        text=df.to_string(index=False)
        text=f'Tmin {self._T_min} \n'+text
        return text
    


    def solve(self, input_text):
        pinch = PyPinch(input_text)
        pinch.shiftTemperatures()
        pinch.constructTemperatureInterval()
        pinch.constructProblemTable()
        pinch.constructHeatCascade()
        hot_util = pinch.hotUtility*1e3  # W
        cold_uti = pinch.coldUtility*1e3  # W
        return hot_util, cold_uti



class Plant(object):
    def __init__(self,parameters) -> None:
        self._isentropic_eff_mt = parameters["isentropic_eff_mt"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._isentropic_eff_st = parameters["isentropic_eff_st"]
        self._mechanical_eff = parameters["mechanical_eff"]
        self._carbonator_eff = parameters["carbonator_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_carbonator_distance = parameters["storage_carbonator_distance"]
        self._cooling_eff = parameters["cooling_eff "]
        self._pres_loss_rec = parameters["pres_loss_rec"]
        self._pres_loss_stoic = parameters["pres_loss_stoic"]
        self._pres_loss_mix = parameters["pres_loss_mix"]
        self._delta_T_pinch = parameters["delta_T_pinch"]
        self._T_amb = parameters["T_amb"]
        self._delta_h=3178.6*1000 # J/kg CaO

    def solve(self,input):
        # unknown pressures
        results = self._calculate_unknown_pressure_by_loss(input)
        results = {**input, **results}
        results["T_min"] = self._delta_T_pinch
        results["T_atm"] = self._T_amb
        # carbonator
        carbonator_results = self.carbonator(input["T_co2_reactor_in"],
                                             input["T_cao_reactor_in"],
                                             input["m_cao_i"],
                                             input["T_reaction"],
                                             input["p_reaction"],
                                             input["cao_conversion"])
        results.update(carbonator_results)
        # main turbine
        T_co2_main_turbine_out, main_turbine_power = self.main_turbine_power(
            input["T_reaction"],
            input["p_reaction"],
            results["p_co2_main_turbine_out"],
            results["m_co2_mix"])
        results["T_co2_main_turbine_out"] = T_co2_main_turbine_out
        results["main_turbine_power"] = main_turbine_power
        # compressor
        T_co2_compressor_out, compressor_power = self.compressor_power(
            input["T_co2_compressor_in"],
            results["p_co2_compressor_in"],
            results["p_co2_compressor_out"],
            results["m_co2_rec"])
        results["T_co2_compressor_out"] = T_co2_compressor_out
        results["compressor_power"] = -compressor_power
        # storage turbine
        T_co2_storage_turbine_out, storage_turbine_power = \
            self.storage_turbine_power(
                input["T_co2_storage_turbine_in"],
                results["p_co2_storage_turbine_in"],
                results["p_co2_storage_turbine_out"],
                results["m_co2_stoic"])
        results["T_co2_storage_turbine_out"] = T_co2_storage_turbine_out
        results["storage_turbine_power"] = storage_turbine_power
        # CO2 mixture
        T_co2_mix = self._mix_co2(
            results["m_co2_stoic"],
            results["T_co2_storage_turbine_out"],
            results["m_co2_rec"],
            results["T_co2_compressor_out"],
            results["p_co2_mix"])
        results["T_co2_mix"] = T_co2_mix
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["m_cao_i"],
            results["m_cao_unr"],
            results["m_caco3"])*(-1)
        # pinch point analysis
        pa = Pinch_point_analyzer(results)

        pa_text=pa.write_pyPinch_data_text()
        hot_util, cold_util = pa.solve(pa_text)
        results["hot_utility"] = hot_util
        results["cold_utility"] = cold_util
        # cooling power
        results["cooling_power"] = self.cooling_power(cold_util)*(-1)
        # summary
        results["total_auxiliary_power"] = results["cooling_power"] + \
            results["conveying_power"]
        results["total_net_power"] = results["main_turbine_power"] + \
            results["storage_turbine_power"] + \
            results["compressor_power"] + \
            results["total_auxiliary_power"]
        results["plant_eff"]=results["total_net_power"]/(input["m_cao_i"]* \
                    input["cao_conversion"]*self._delta_h)

        return results

    def _mix_co2(self, m1, T1, m2, T2, p):
        H1 = m1*CP.PropsSI('H', 'T', T1+273.15, 'P', p, "REFPROP::co2")
        H2 = m2*CP.PropsSI('H', 'T', T2+273.15, 'P', p, "REFPROP::co2")
        h_mix = (H1+H2)/(m1+m2)
        T_mix = CP.PropsSI('T', 'H', h_mix, 'P', p, "REFPROP::co2")-273.15
        return T_mix

    def _calculate_unknown_pressure_by_loss(self,input):
        results = {}
        results["p_co2_main_turbine_out"] = input["p_reaction"]/input["expansion_ratio"]
        results["p_co2_compressor_in"] = results["p_co2_main_turbine_out"]*(1-self._pres_loss_rec)
        results["p_co2_mix"] = input["p_reaction"]/(1-self._pres_loss_mix)
        results["p_co2_compressor_out"] = results["p_co2_mix"]
        results["p_co2_storage_turbine_in"] = input["p_co2_storage"]*(1-self._pres_loss_stoic)
        results["p_co2_storage_turbine_out"] = results["p_co2_compressor_out"]
        return results

    def compressor_power(self, Ti, pi, po, mass_rate, fluid="REFPROP::co2"):
        hi = CP.PropsSI('H', 'T', Ti+273.15, 'P', pi, fluid)
        si = CP.PropsSI('S', 'T', Ti+273.15, 'P', pi, fluid)
        ho_s = CP.PropsSI('H', 'P', po, 'S', si, fluid)
        ho_c = (ho_s-hi)/self._isentropic_eff_mc+hi
        To = CP.PropsSI('T', 'P', po, 'H', ho_c, fluid)-273.15
        W = ((ho_c-hi)/self._mechanical_eff)*mass_rate
        return To, W

    def _turbine_power(self, Ti, pi, po, mass_rate, isentropic_eff, fluid="REFPROP::co2"):
        hi = CP.PropsSI('H', 'T', Ti+273.15, 'P', pi, fluid)
        si = CP.PropsSI('S', 'T', Ti+273.15, 'P', pi, fluid)
        ho_s = CP.PropsSI('H', 'P', po, 'S', si, fluid)
        ho_c = hi-(hi-ho_s)*isentropic_eff
        To = CP.PropsSI('T', 'P', po, 'H', ho_c, fluid)-273.15
        W = ((hi-ho_c)*self._mechanical_eff)*mass_rate
        return To, W

    def main_turbine_power(self, Ti, pi, po, mass_rate, fluid="REFPROP::co2"):
        return self._turbine_power(Ti, pi, po, mass_rate,
                                   self._isentropic_eff_mt, fluid)

    def storage_turbine_power(self, Ti, pi, po, mass_rate, fluid="REFPROP::co2"):
        return self._turbine_power(Ti, pi, po, mass_rate,
                                   self._isentropic_eff_st, fluid)

    def carbonator(self, Ti_co2, Ti_caO, mass_cao_i, Tr, pr, X):
        M_cao = 56e-3  # kg/mol
        M_caco3 = 100e-3  # kg/mol
        M_CO2 = 44e-3  # kg/mol

        delta_H_Tref = -178e3  # J/mole
        Tref = 20
        pw = Cp0mass_Wrapper()
        cp_cao_mean_Tref_Tr = pw.cp0mass_mean("cao", Tref, Tr)
        cp_caco3_mean_Tref_Tr = pw.cp0mass_mean("caco3", Tref, Tr)
        cp_co2_mean_Tref_Tr = pw.cp0mass_mean("co2", Tref, Tr, pr)
        delta_H_Tr = delta_H_Tref+(cp_caco3_mean_Tref_Tr*M_caco3
                                   - cp_cao_mean_Tref_Tr*M_cao
                                   - cp_co2_mean_Tref_Tr*M_CO2)*(Tr-Tref)
        # print(f"specific reaction heat released at Td=875oC is {delta_H_Tr/1000} KJ/mol")

        mole_cao_i = mass_cao_i/M_cao
        mole_cao_r = mole_cao_i*X
        mass_cao_r = mole_cao_r*M_cao
        mole_caco3_o = mole_cao_r
        mass_caco3_o = mole_caco3_o*M_caco3
        mass_co2_stoic = mole_cao_r*M_CO2

        Q_heat = mole_cao_r*delta_H_Tr*self._carbonator_eff
        # print(f"total reaction heat released at Td=875oC is {Q_heat/1000} kW")

        cp_cao_Tr = pw.cp0mass("cao", Tr)
        cp_cao_Ti = pw.cp0mass("cao", Ti_caO)
        hi_Tr = CP.PropsSI('H', 'T', Tr+273.15, 'P', pr, "REFPROP::co2")
        hi_Ti = CP.PropsSI('H', 'T', Ti_co2+273.15, 'P', pr, "REFPROP::co2")
        mass_co2_mix = (-Q_heat-mass_cao_i*(cp_cao_Tr*Tr-cp_cao_Ti*Ti_caO))/(hi_Tr-hi_Ti)

        results = {}
        results["m_co2_mix"] = mass_co2_mix
        results["m_co2_stoic"] = mass_co2_stoic
        results["m_co2_rec"] = mass_co2_mix-mass_co2_stoic
        results["m_cao_i"] = mass_cao_i
        results["m_cao_unr"] = mass_cao_i-mass_cao_r
        results["m_caco3"] = mass_caco3_o
        results["excess_index"] = mass_co2_mix/mass_co2_stoic

        return results

    def conveying_power(self, m_cao_i, m_cao_unr, m_caco3_o):
        return self._convey_consumption*self._storage_carbonator_distance * \
            (m_cao_i+m_cao_unr+m_caco3_o)

    def cooling_power(self, cold_utility):
        return self._cooling_eff*cold_utility




