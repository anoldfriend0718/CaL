import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")

from scipy import interpolate
from pyPinch import PyPinch
import pandas as pd
import numpy as np
import CoolProp.CoolProp as CP



class Cp0mass_Wrapper(object):
    def __init__(self, flue_gas_composistion, deconbonized_rate) -> None:
        self.f_cao, self.f_caco3 = self._read_solid_properties()
        total_mole_frac_s = flue_gas_composistion["co2"] + \
            flue_gas_composistion["n2"] +\
            flue_gas_composistion["o2"]
        self._co2_mole_frac_s = flue_gas_composistion["co2"]/total_mole_frac_s
        self._n2_mole_frac_s = flue_gas_composistion["n2"]/total_mole_frac_s
        self._o2_mole_frac_s = flue_gas_composistion["o2"]/total_mole_frac_s
        norm_flue_gas_composition={}
        norm_flue_gas_composition["co2"] = self._co2_mole_frac_s
        norm_flue_gas_composition["n2"] = self._n2_mole_frac_s
        norm_flue_gas_composition["o2"] = self._o2_mole_frac_s
        self._norm_flue_gas_composition_s=norm_flue_gas_composition

        total_mole_frac_e = total_mole_frac_s-flue_gas_composistion["co2"]*deconbonized_rate
        self._co2_mole_frac_e = flue_gas_composistion["co2"]*(1-deconbonized_rate)/total_mole_frac_e
        self._n2_mole_frac_e = flue_gas_composistion["n2"]/total_mole_frac_e
        self._o2_mole_frac_e = flue_gas_composistion["o2"]/total_mole_frac_e

    # T: oC
    # cp: J/kg/K
    def cp0mass(self, material, T, p=1e5):
        material = material.lower()
        if material == "flue_gas":
            cp = self._cp0mass_flue_gas(T, p)
        elif material == "decarbonized_flue_gas":
            cp = self._cp0mass_decarbonized_flue_gas(T, p)
        elif material == "co2":
            cp = self._cp0mass_co2(T, p)
        elif material == "cao":
            cp = self._cp0mass_cao(T)
        elif material == "caco3":
            cp = self._cp0mass_caco3(T)
        else:
            raise ValueError(f"Cp0mass of {material} is not supported")
        return cp

    # T: oC
    # cp: J/kg/K
    def cp0mass_mean(self, material, Ti, To, p=1e5, interval=10):
        if To < Ti:
            temp = Ti
            Ti = To
            To = temp
        Ts = self._arange(Ti, To, interval)
        cps = []
        for T in Ts:
            cp = self.cp0mass(material, T, p)
            cps.append(cp)
        cp_mean = np.mean(cps)
        return cp_mean

    def _read_solid_properties(self):
        data_csv = f"{CaLRepo}/data/cpmass_cao_caco3.csv"
        df = pd.read_csv(data_csv)
        f_cao = interpolate.interp1d(df["TEMP"], df["CAO"])
        f_caco3 = interpolate.interp1d(df["TEMP"], df["CACO3"])
        return f_cao, f_caco3

    def _cp0mass_cao(self, T):
        cp = float(self.f_cao(T))
        return cp

    def _cp0mass_caco3(self, T):
        cp = float(self.f_caco3(T))
        return cp

    def _cp0mass_co2(self, T, p):
        fluid = 'REFPROP::co2'
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp

    def _cp0mass_flue_gas(self, T, p):
        fluid = self.get_flue_gas_refprop_name()
        cp = self._cp0mass_gas(T, p, fluid)
        return cp

    def get_flue_gas_refprop_name(self):
        fluid = f'REFPROP::co2[{self._co2_mole_frac_s}]&\
                  nitrogen[{self._n2_mole_frac_s}]&\
                  oxygen[{self._o2_mole_frac_s}]'.replace(" ", "")

        return fluid

    def _cp0mass_decarbonized_flue_gas(self, T, p):
        fluid = self.get_decarbonized_flue_gas_refprop_name()
        cp = self._cp0mass_gas(T, p, fluid)
        return cp

    def get_decarbonized_flue_gas_refprop_name(self):
        fluid = f'REFPROP::co2[{self._co2_mole_frac_e}]&\
            nitrogen[{self._n2_mole_frac_e}]&\
            oxygen[{self._o2_mole_frac_e}]'.replace(" ", "")

        return fluid

    def _cp0mass_gas(self, T, p, fluid):
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp

    def _arange(self, start, stop, step=1, endpoint=True):
        arr = np.arange(start, stop, step)

        if endpoint and arr[-1] != stop:
            arr = np.concatenate([arr, [stop]])

        return arr


class Pinch_point_analyzer(object):
    def __init__(self, inputs) -> None:
        self._m_caco3 = inputs["m_caco3_out"]
        self._m_cao_unr = inputs["m_cao_unr_out"]
        self._m_cao_i = inputs["m_cao_in"]
        self._m_decarbonized_flue_gas = inputs["m_deconbonized_flue_gas_out"]
        self._m_flue_gas_in = inputs["m_flue_gas_in"]
        self._m_water_in = inputs["m_water_in"]

        self._T_carb = inputs["T_carb"]
        self._p_amb = inputs["p_amb"]
        self._T_amb = inputs["T_amb"]
        self._T_flue_gas_compressor_out = inputs["T_flue_gas_compressor_out"]  # todo:T_flue_gas_compressor_out
        self._p_flue_gas_compressor_out = inputs["p_flue_gas_compressor_out"]
        self._T_flue_gas_reactor_in = inputs["T_flue_gas_reactor_in"]
        self._T_cao_reactor_in = inputs["T_cao_reactor_in"]
        self._T_water_reactor_in = inputs["T_water_reactor_in"]
        self._T_decarbonized_flue_gas_out = inputs["T_decarbonized_flue_gas_out"]
        self._p_decarbonized_flue_gas_out = self._p_amb
        self._T_delta_pinch = inputs["T_delta_pinch"]

        self._p_carb = inputs["p_carb"]
        self._p_water = inputs["p_water"]

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
                        'P', self._p_carb, decarb_flue_gas_name) -
             CP.PropsSI('H', 'T', self._T_decarbonized_flue_gas_out + 273.15,
                        'P', self._p_decarbonized_flue_gas_out, decarb_flue_gas_name))
        pinch_point_data["CP"]["H_flue_gas_decarb"] = pinch_point_data["ENERGY"]["H_flue_gas_decarb"] / \
            (self._T_carb-self._T_decarbonized_flue_gas_out)

        # C1: flue gas
        pinch_point_data["TSUPPLY"]["C_flue_gas"] = self._T_flue_gas_compressor_out  # todo:T_flue_gas_compressor_out
        pinch_point_data["TTARGET"]["C_flue_gas"] = self._T_flue_gas_reactor_in
        pinch_point_data["FLOWRATE"]["C_flue_gas"] = self._m_flue_gas_in
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["C_flue_gas"] = self._m_flue_gas_in * \
            (CP.PropsSI('H', 'T', self._T_flue_gas_reactor_in+273.15,
                        'P', self._p_carb, flue_gas_name) -
             CP.PropsSI('H', 'T', self._T_flue_gas_compressor_out+273.15,
                        'P', self._p_flue_gas_compressor_out, flue_gas_name))
        pinch_point_data["CP"]["C_flue_gas"] = pinch_point_data["ENERGY"]["C_flue_gas"] / \
            (self._T_flue_gas_reactor_in-self._T_flue_gas_compressor_out)

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

    def solve(self, input_text):
        pinch = PyPinch(input_text)
        pinch.shiftTemperatures()
        pinch.constructTemperatureInterval()
        pinch.constructProblemTable()
        pinch.constructHeatCascade()
        hot_util = pinch.hotUtility*1e3  # W
        cold_uti = pinch.coldUtility*1e3  # W
        return hot_util, cold_uti


class CarbonatorSide(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"], parameters["decarbonized_rate"])
        self._flue_gas_composition = self._pw._norm_flue_gas_composition_s
        self._decarbonized_rate = parameters["decarbonized_rate"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mechanical_eff"]
        self._flue_gas_pressure_loss_1 = parameters["flue_gas_pressure_loss_1"]
        self._flue_gas_pressure_loss_2 = parameters["flue_gas_pressure_loss_2"]
        
        self._vol_rate_flue_gas = parameters["vol_rate_flue_gas"]  # m3/s
        self._T_flue_gas_i = parameters["T_flue_gas"]
        self._T_carb = parameters["T_carb"]
        # self._p_carb=parameters["p_carb"]
        self._cao_conversion = parameters["cao_conversion"]
        self._T_water_reactor_out = parameters["T_water_reactor_out"]
        self._p_water = parameters["p_water"]  # todo: change the water pressure

        self._carbonator_eff = parameters["carbonator_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_carbonator_distance = parameters["storage_carbonator_distance"]
        self._cooling_eff = parameters["cooling_eff "]
        self._delta_T_pinch = parameters["delta_T_pinch"]
        self._p_amb = parameters["p_amb"]
        self._T_amb = parameters["T_amb"]
        self._p_flue_gas_i = self._p_amb
        self._delta_H_Tref = -178e3  # J/mole CaO
        self._delta_h = 3178.6*1000  # J/kg CaO

    def solve(self, input):
        results = {}
        results["flue_gas_composition"] = self._flue_gas_composition
        results["deconbonized_rate"] = self._decarbonized_rate
        results["vol_rate_flue_gas"]=self._vol_rate_flue_gas 
        results["T_carb"] = self._T_carb
        results["p_amb"] = self._p_amb
        results["T_amb"] = self._T_amb
        results["T_delta_pinch"] = self._delta_T_pinch
        results["T_flue_gas_in"] = self._T_flue_gas_i
        results["p_flue_gas_in"] = self._p_flue_gas_i
        results["T_decarbonized_flue_gas_out"] = self._T_flue_gas_i  # keep same as the inlet flue gas
        results["T_water_reactor_out"] = self._T_water_reactor_out
        results["p_carb"] = self._p_amb/(1-self._flue_gas_pressure_loss_2)  # 压力略高于大气压，需扣减压力损失
        results["p_water"] = self._p_water
        results["cao_conversion"] = self._cao_conversion
        results = {**input, **results}

       # carbonator(self, Ti_flue_gas, Ti_cao, Ti_water,To_water,Tcarb, pcarb, X):
        carbonator_results = self.carbonator(input["T_flue_gas_reactor_in"],
                                             input["T_cao_reactor_in"],
                                             input["T_water_reactor_in"],
                                             self._T_water_reactor_out,
                                             self._T_carb,
                                             results["p_carb"],
                                             self._cao_conversion)
        results.update(carbonator_results)

        # flue gas compressor
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        flue_gas_pi = self._p_amb
        flue_gas_po = results["p_carb"]/(1-self._flue_gas_pressure_loss_1)
        flue_gas_mass_rate = results["m_flue_gas_in"]
        flue_gas_compressor_results = self.compressor_power(self._T_flue_gas_i,
                                                            flue_gas_pi,
                                                            flue_gas_po,
                                                            flue_gas_mass_rate,
                                                            flue_gas_name)
        results.update(flue_gas_compressor_results)

        # Q hot water
        h_To_water = CP.PropsSI('H', 'T', self._T_water_reactor_out+273.15,
                                'P', self._p_water, "REFPROP::water")
        h_Ti_water = CP.PropsSI('H', 'T', self._T_amb+273.15,
                                'P', self._p_water, "REFPROP::water")
        Q_hot_water = results["m_water_in"]*(h_To_water-h_Ti_water)
        results["Q_hot_water"] = Q_hot_water
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["m_cao_in"],
            results["m_cao_unr_out"],
            results["m_caco3_out"])*(-1)
        # pinch point analysis
        pa = Pinch_point_analyzer(results)
        pa_text = pa.write_pyPinch_data_text()
        hot_util, cold_util = pa.solve(pa_text)
        results["hot_utility"] = hot_util
        results["cold_utility"] = cold_util
        # cooling power
        results["cooling_power"] = self.cooling_power(cold_util)*(-1)
        # summary
        results["carb_auxiliary_power"] = results["cooling_power"] + \
            results["conveying_power"]+results["flue_gas_compressor_power"]

        results["carb_heat_rec_eff"] = results["Q_hot_water"]/((results["m_cao_in"] -
                                                        results["m_cao_unr_out"])*self._delta_h)

        return results

    def compressor_power(self, Ti, pi, po, mass_rate, fluid):
        hi = CP.PropsSI('H', 'T', Ti+273.15, 'P', pi, fluid)
        si = CP.PropsSI('S', 'T', Ti+273.15, 'P', pi, fluid)
        ho_s = CP.PropsSI('H', 'P', po, 'S', si, fluid)
        ho_c = (ho_s-hi)/self._isentropic_eff_mc+hi
        To = CP.PropsSI('T', 'P', po, 'H', ho_c, fluid)-273.15
        W = ((ho_c-hi)/self._mechanical_eff)*mass_rate
        results = {}
        results["p_flue_gas_compressor_out"] = po
        results["T_flue_gas_compressor_out"] = To
        results["flue_gas_compressor_power"] = W*(-1)
        return results

    def carbonator(self, Ti_flue_gas, Ti_cao, Ti_water, To_water, Tcarb, pcarb, X):
        M_cao = 56e-3  # kg/mol
        M_caco3 = 100e-3  # kg/mol
        M_CO2 = 44e-3  # kg/mol

        Tref = 20
        cp_cao_mean_Tref_Tr = self._pw.cp0mass_mean("cao", Tref, Tcarb)
        cp_caco3_mean_Tref_Tr = self._pw.cp0mass_mean("caco3", Tref, Tcarb)

        delta_H_Tr = self._delta_H_Tref+(cp_caco3_mean_Tref_Tr*M_caco3
                                         - cp_cao_mean_Tref_Tr*M_cao)*(Tcarb-Tref)\
            - (CP.PropsSI('H', 'T', Tcarb+273.15, 'P', pcarb, "REFPROP::co2") -
               CP.PropsSI('H', 'T', Tref+273.15, 'P', pcarb, "REFPROP::co2"))*M_CO2
        # print((cp_caco3_mean_Tref_Tr*M_caco3- cp_cao_mean_Tref_Tr*M_cao)*(Tcarb-Tref))
        # print((CP.PropsSI('H', 'T', Tcarb+273.15, 'P', pcarb, "REFPROP::co2") -
        #        CP.PropsSI('H', 'T', Tref+273.15, 'P', pcarb, "REFPROP::co2"))*M_CO2)
        # print(f"specific reaction heat released is {delta_H_Tr/1000} KJ/mol")

        mass_co2_i = self._vol_rate_flue_gas*self._flue_gas_composition["co2"] *\
            CP.PropsSI('D', 'T', self._T_flue_gas_i+273.15, 'P', self._p_flue_gas_i, "REFPROP::co2")
        mole_co2_i = mass_co2_i/M_CO2
        mole_co2_r = mole_co2_i*self._decarbonized_rate
        mass_co2_r = mole_co2_r*M_CO2

        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        mass_flue_gas = self._vol_rate_flue_gas *\
            CP.PropsSI('D', 'T', self._T_flue_gas_i+273.15, 'P', self._p_flue_gas_i, flue_gas_name)
        mass_deconbonized_flue_gas_out = mass_flue_gas-mole_co2_r*M_CO2

        mole_cao_r = mole_co2_r
        mass_cao_r = mole_cao_r*M_cao
        mole_cao_i = mole_cao_r/X
        mass_cao_i = mole_cao_i*M_cao

        mole_caco3_o = mole_cao_r
        mass_caco3_o = mole_caco3_o*M_caco3

        Q_heat = mole_cao_r*delta_H_Tr*self._carbonator_eff
        # print(f"total reaction heat released is {Q_heat/1000} kW")

        cp_cao_Tr = self._pw.cp0mass("cao", Tcarb)
        cp_cao_Ti = self._pw.cp0mass("cao", Ti_cao)

        h_Tr_flue_gas = CP.PropsSI('H', 'T', Tcarb+273.15, 'P', pcarb, flue_gas_name)
        h_Ti_flue_gas = CP.PropsSI('H', 'T', Ti_flue_gas+273.15, 'P', pcarb, flue_gas_name)  # 等压过程

        h_To_water = CP.PropsSI('H', 'T', To_water+273.15, 'P', self._p_water,
                                "REFPROP::water")  # todo: change water pressure
        h_Ti_water = CP.PropsSI('H', 'T', Ti_water+273.15, 'P', self._p_water, "REFPROP::water")
        mass_water = (-Q_heat-mass_cao_i*(cp_cao_Tr*Tcarb-cp_cao_Ti*Ti_cao) -
                      mass_flue_gas*(h_Tr_flue_gas-h_Ti_flue_gas)) /\
                     (h_To_water-h_Ti_water)
        # print(f"mass flow rate of water is {mass_water} kg/s ")

        results = {}
        deconbonized_flue_gas_composition = {}
        deconbonized_flue_gas_composition["co2"] = self._pw._co2_mole_frac_e
        deconbonized_flue_gas_composition["n2"] = self._pw._n2_mole_frac_e
        deconbonized_flue_gas_composition["o2"] = self._pw._o2_mole_frac_e
        results["deconbonized_flue_gas_composition"] = deconbonized_flue_gas_composition

        results["m_flue_gas_in"] = mass_flue_gas
        results["m_water_in"] = mass_water
        results["m_co2_capture"] = mass_co2_r
        results["m_cao_in"] = mass_cao_i
        results["m_deconbonized_flue_gas_out"] = mass_deconbonized_flue_gas_out
        results["m_cao_unr_out"] = mass_cao_i-mass_cao_r
        results["m_caco3_out"] = mass_caco3_o
        results["delta_H_Tcarb"]=delta_H_Tr
        return results

    def conveying_power(self, m_cao_i, m_cao_unr, m_caco3_o):
        return self._convey_consumption*self._storage_carbonator_distance * \
            (m_cao_i+m_cao_unr+m_caco3_o)

    def cooling_power(self, cold_utility):
        return self._cooling_eff*cold_utility


class CalcinerSide(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"], parameters["decarbonized_rate"])
        self._flue_gas_composition = self._pw._norm_flue_gas_composition_s
        self._decarbonized_rate = parameters["decarbonized_rate"]
        
        self._calciner_eff = parameters["calciner_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_carbonator_distance = parameters["storage_carbonator_distance"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mechanical_eff"]

        self._cao_conversion = parameters["cao_conversion"]
        self._p_amb = parameters["p_amb"]
        self._T_amb = parameters["T_amb"]
        self._T_calc = parameters["T_calc"]
        self._T_cooling_co2=parameters["T_cooling_co2"]
        self._p_co2_storage=parameters["p_co2_storage"]
        self._n_compression=parameters["n_compression"]
        self._p_calc = self._p_amb
        self._deltaTmin_SSHX = parameters["deltaTmin_SSHX"]
        self._deltaTmin_SGHX = parameters["deltaTmin_SGHX"]
        self._cooling_eff=parameters["cooling_eff"]
        self._delta_H_Tref = 178e3  # J/mole CaO

    def solve(self, input):
        mass_camix_in=input["mass_camix_in"]
        mfrac=input["mfrac"]
        results = {}
        calciner_results=self.calciner(mass_camix_in,mfrac)
        results.update(calciner_results)
        mass_co2_out=results["mass_co2_out"]
        T_co2_o=results["T_co2_o"]
        # print(mass_co2_out,T_co2_o,self._p_amb,self._p_co2_storage,self._T_cooling_co2,self._n_compression)
        compressor_results=self.co2_multiple_stage_compressor(mass_co2_out,
            T_co2_o,self._p_amb,self._p_co2_storage,self._T_cooling_co2,self._n_compression)
        results.update(compressor_results)

        results["conveying_power"]=self.conveying_power(mass_camix_in,results["mass_cao_out"])*(-1)
        results["cooling_power"]=self.cooling_power(results["cooling_energy"])*(-1)
        results["calc_auxiliary_power"]=results["conveying_power"]+results["cooling_power"]+\
            results["compressor_power"]

        return results

    def calciner(self, m, mfrac):
        # inflow mass rate is m, in which CaCO3 ocuppied X percentage
        M_cao = 56e-3  # kg/mol
        M_caco3 = 100e-3  # kg/mol
        M_CO2 = 44e-3  # kg/mol
        X = self._cao_conversion

        M_camix=X*M_caco3+(1-X)*M_cao
        mole_camix=m/M_camix
        mole_CO2_out = X*mole_camix
        mass_co2_out = mole_CO2_out*M_CO2
        mass_cao_out = mole_camix*M_cao
        mass_camix_1 = mfrac*m
        mass_camix_2 = (1-mfrac)*m

        # Hot CaO ~ cold CaCO3 heat exchanger
        T_camix_1, T_cao_o = self.cao_camix_heat_exchanger(mass_camix_1, mass_cao_out)
        # Hot CO2 ~ cold CaCO3 heat exchanger
        T_camix_2, T_co2_o = self.co2_camix_heat_exchanger(mass_camix_2, mass_co2_out)
        # CaCO3 mixer
        T_camix_reactor_in = self.caco3_mixer(mass_camix_1, T_camix_1, mass_camix_2, T_camix_2)
        # calciner
        delta_H_r_Tcalc = self.reaction_heat()
        delta_H=self.cp_camix(self._T_calc, X)*self._T_calc - \
                self.cp_camix(T_camix_reactor_in, X)*T_camix_reactor_in # based on mass
        Qe = X*mole_camix*delta_H_r_Tcalc+m*delta_H
        We = Qe/self._calciner_eff

        results = {}
        results["mfrac"]=mfrac
        results["mass_camix_1"]=mass_camix_1
        results["mass_camix_2"]=mass_camix_2
        results["mass_cao_out"]=mass_cao_out
        results["mass_co2_out"]=mass_co2_out
        results["T_cao_o"] = T_cao_o
        results["T_co2_o"] = T_co2_o
        results["T_camix_1"] = T_camix_1
        results["T_camix_2"] = T_camix_2
        results["T_camix_reactor_in"] = T_camix_reactor_in
        results["delta_H_Tcalc"] = delta_H_r_Tcalc
        results["Qe_calc"] = Qe
        results["We_calc"] = We
        return results




    def reaction_heat(self):
        Tref = self._T_amb
        Tcalc = self._T_calc
        pcalc = self._p_calc
        cp_caco3_mean = self._pw.cp0mass_mean("caco3", Tref, Tcalc)
        cp_cao_mean = self._pw.cp0mass_mean("cao", Tref, Tcalc)
        M_cao = 56e-3  # kg/mol
        M_caco3 = 100e-3  # kg/mol
        M_CO2 = 44e-3  # kg/mol

        # q1=(cp_cao_mean*M_cao-cp_caco3_mean*M_caco3)*(Tcalc-Tref)
        # q2=(CP.PropsSI('H', 'T', Tcalc+273.15, 'P', pcalc, "REFPROP::co2") -
        #       CP.PropsSI('H', 'T', Tref+273.15, 'P', pcalc, "REFPROP::co2"))*M_CO2
        # delta_H_Tr=self._delta_H_Tref+q1+q2
        delta_H_Tr = self._delta_H_Tref+(cp_cao_mean*M_cao-cp_caco3_mean*M_caco3)*(Tcalc-Tref)\
            + (CP.PropsSI('H', 'T', Tcalc+273.15, 'P', pcalc, "REFPROP::co2") -
                CP.PropsSI('H', 'T', Tref+273.15, 'P', pcalc, "REFPROP::co2"))*M_CO2
        return delta_H_Tr

    def cao_camix_heat_exchanger(self, m_camix, m_cao):
        T_camix_in = self._T_amb
        T_cao_in = self._T_calc
        cp_cao_mean_Tamb_Tr = self._pw.cp0mass_mean("cao", T_camix_in, T_cao_in)
        cp_camix_mean_Tamb_Tr = self.cp_camix_mean(T_camix_in, T_cao_in, self._cao_conversion)
        deltaTmin = self._deltaTmin_SSHX
        T_camix_out_n_1, T_cao_out_n_1 = self.SS_heat_exchanger(
            m_camix, cp_camix_mean_Tamb_Tr, m_cao, cp_cao_mean_Tamb_Tr,
            T_camix_in, T_cao_in, deltaTmin)
        T_camix_out_n = 0
        T_cao_out_n = 0
        eps = 1.5

        while abs(T_camix_out_n_1-T_camix_out_n) > eps or abs(T_cao_out_n_1-T_cao_out_n) > eps:
            T_camix_out_n = T_camix_out_n_1
            T_cao_out_n = T_cao_out_n_1
            cp_camix_mean = self.cp_camix_mean(T_camix_in, T_camix_out_n, self._cao_conversion)
            cp_cao_mean = self._pw.cp0mass_mean("cao", T_cao_in, T_cao_out_n)
            T_camix_out_n_1, T_cao_out_n_1 = self.SS_heat_exchanger(
                m_camix, cp_camix_mean, m_cao, cp_cao_mean,
                T_camix_in, T_cao_in, deltaTmin)
        return T_camix_out_n_1, T_cao_out_n_1

    def SS_heat_exchanger(self, mc, cpc, mh, cph, Tic, Tih, deltaTmin):
        if mc*cpc > mh*cph:
            Toh = Tic+deltaTmin
            Toc = mh*cph*(Tih-Toh)/(mc*cpc)+Tic
        else:
            Toc = Tih-deltaTmin
            Toh = Tih-mc*cpc*(Toc-Tic)/(mh*cph)
        return Toc, Toh

    def co2_camix_heat_exchanger(self, m_camin, m_co2):
        T_camix_in = self._T_amb
        T_co2_in = self._T_calc
        cp_camix_mean_Tamb_Tr = self.cp_camix_mean(T_camix_in, T_co2_in, self._cao_conversion)
        deltaTmin = self._deltaTmin_SGHX
        T_camix_out_n_1, T_co2_out_n_1 = self.SG_heat_exchanger(m_camin, cp_camix_mean_Tamb_Tr, m_co2,
                                                                T_camix_in, T_co2_in, self._p_calc, deltaTmin)
        T_camix_out_n = 0
        T_co2_out_n = 0
        eps = 1.5

        while abs(T_camix_out_n_1-T_camix_out_n) > eps or abs(T_co2_out_n_1-T_co2_out_n) > eps:
            T_camix_out_n = T_camix_out_n_1
            T_co2_out_n = T_co2_out_n_1
            cp_camix_mean = self.cp_camix_mean(T_camix_in, T_camix_out_n, self._cao_conversion)
            T_camix_out_n_1, T_co2_out_n_1 = self.SG_heat_exchanger(m_camin, cp_camix_mean, m_co2,
                                                                    T_camix_in, T_co2_in, self._p_calc, deltaTmin)
        return T_camix_out_n_1, T_co2_out_n_1

    def SG_heat_exchanger(self, mc, cpc, mh, Tic, Tih, pcalc, deltaTmin):
        # Hot fluid is Gas  CO2:
        Toh = Tic+deltaTmin
        Qh = mh*(CP.PropsSI('H', 'T', Tih+273.15, 'P', pcalc, "REFPROP::co2") -
                 CP.PropsSI('H', 'T', Toh+273.15, 'P', pcalc, "REFPROP::co2"))
        cph = Qh/(mh*(Tih-Toh))
        if mc*cpc > mh*cph:
            Toc = Qh/(mc*cpc)+Tic
        else:
            Toc = Tih-deltaTmin
            Qc = mc*cpc*(Toc-Tic)
            Hhi = mh*CP.PropsSI('H', 'T', Tih+273.15, 'P', pcalc, "REFPROP::co2")
            Hho = Hhi-Qc
            hhmo = Hho/mh
            Toh = CP.PropsSI('T', 'H', hhmo, 'P', pcalc, "REFPROP::co2")-273.15
            cph = Qc/(mh*(Tih-Toh))
            if mc*cpc > mh*cph:
                raise Exception("some error happen in SG_heat_exchanger")
        return Toc, Toh

    def caco3_mixer(self, m1, T1, m2, T2):
        Ton = 0
        Ton_1 = (m1*T1+m2*T2)/(m1+m2)
        eps = 1
        while abs(Ton_1-Ton) < eps:
            Ton = Ton_1
            cp1 = self.cp_camix_mean(T1, Ton_1, self._cao_conversion)
            cp2 = self.cp_camix_mean(T2, Ton_1, self._cao_conversion)
            Ton_1 = (m1*cp1*T1+m2*cp2*T2)/(m1*cp1+m2*cp2)
        return Ton_1

    def cp_camix_mean(self, T1, T2, X):
        cp_cao_mean = self._pw.cp0mass_mean("cao", T1, T2)
        cp_caco3_mean = self._pw.cp0mass_mean("caco3", T1, T2)
        mX=self.convert_X_to_mX(X)
        cp_camix_mean = cp_caco3_mean*mX+cp_cao_mean*(1-mX)
        return cp_camix_mean

    def cp_camix(self, T1, X):
        cp_cao = self._pw.cp0mass("cao", T1)
        cp_caco3 = self._pw.cp0mass("caco3", T1)
        mX=self.convert_X_to_mX(X)
        cp_camix = cp_caco3*mX+cp_cao*(1-mX)
        return cp_camix

    def convert_X_to_mX(self,X):
        M_cao = 56e-3  # kg/mol
        M_caco3 = 100e-3  # kg/mol
        molar_mass=X*M_caco3+(1-X)*M_cao
        mX=X*M_caco3/molar_mass
        return mX

    def co2_multiple_stage_compressor(self,m,Ti,pi,po,Tc,n):
        pr=pow(po/pi,1./n)
        i=1
        Ws=[]
        Tos=[]
        Qcs=[]
        pi_i=pi
        Ti_i=Ti
        while i<=n:
            po_i=pi_i*pr
            To_i,W_i=self.compressor_power(m,Ti_i, pi_i, po_i, "REFPROP::co2")
            hi_Ti= CP.PropsSI('H', 'T', To_i+273.15, 'P', po_i, "REFPROP::co2")
            hi_Tc = CP.PropsSI('H', 'T', Tc+273.15, 'P', po_i, "REFPROP::co2")
            Qc=m*(hi_Ti-hi_Tc)
            Ws.append(W_i)
            Tos.append(To_i)
            Qcs.append(Qc)
            i=i+1
            Ti_i=Tc
            pi_i=po_i
        Wt=np.sum(Ws)
        Qct=np.sum(Qcs)
        results={}
        results["compressor_power"]=Wt*(-1)
        results["cooling_energy"]=Qct
        return results

    def compressor_power(self, mass_rate,Ti, pi, po, fluid):
        hi = CP.PropsSI('H', 'T', Ti+273.15, 'P', pi, fluid)
        si = CP.PropsSI('S', 'T', Ti+273.15, 'P', pi, fluid)
        ho_s = CP.PropsSI('H', 'P', po, 'S', si, fluid)
        ho_c = (ho_s-hi)/self._isentropic_eff_mc+hi
        To = CP.PropsSI('T', 'P', po, 'H', ho_c, fluid)-273.15
        W = ((ho_c-hi)/self._mechanical_eff)*mass_rate
        return To,W

    def conveying_power(self, m_camix_i,m_cao_o):
        return self._convey_consumption*self._storage_carbonator_distance * \
            (m_camix_i+m_cao_o)

    def cooling_power(self, cold_utility):
        return self._cooling_eff*cold_utility

if __name__ == '__main__':
    parameters = dict()
    flue_gas_composistion = dict()
    flue_gas_composistion["co2"] = 0.1338
    flue_gas_composistion["o2"] = 0.0384
    flue_gas_composistion["n2"] = 0.6975
    parameters["flue_gas_composition"] = flue_gas_composistion
    parameters["isentropic_eff_mc"] = 0.87
    parameters["mechanical_eff"] = 0.97
    parameters["decarbonized_rate"] = 0.9
    parameters["cao_conversion"] = 0.64
    parameters["calciner_eff"] = 0.99
    parameters["convey_consumption"] = 10e3/100
    parameters["storage_carbonator_distance"] = 100
    parameters["T_amb"] = 20
    parameters["p_amb"] = 101325
    parameters["T_calc"] = 900
    parameters["deltaTmin_SSHX"] = 20
    parameters["deltaTmin_SGHX"] = 15
    parameters["T_cooling_co2"]=20
    parameters["p_co2_storage"]=75e5
    parameters["n_compression"]=6
    parameters["cooling_eff"]=0.01

    calcs = CalcinerSide(parameters)
    m = 57.23+29.21
    # mfrac = 57.23/m
    mfrac=0.6180339887498948
    results = calcs.calciner(m, mfrac)
    print(results)
