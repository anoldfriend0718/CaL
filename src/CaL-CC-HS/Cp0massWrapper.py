
import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")
from scipy import interpolate
import numpy as np
import pandas as pd
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
        f_cao = interpolate.interp1d(df["TEMP"], df["CAO"],fill_value="extrapolate")
        f_caco3 = interpolate.interp1d(df["TEMP"], df["CACO3"],fill_value="extrapolate")
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

    def cp_camix_mean(self, T1, T2, X):
        cp_cao_mean = self.cp0mass_mean("cao", T1, T2)
        cp_caco3_mean = self.cp0mass_mean("caco3", T1, T2)
        mX=self.convert_X_to_mX(X)
        cp_camix_mean = cp_caco3_mean*mX+cp_cao_mean*(1-mX)
        return cp_camix_mean

    def cp_camix(self, T1, X):
        cp_cao = self.cp0mass("cao", T1)
        cp_caco3 = self.cp0mass("caco3", T1)
        mX=self.convert_X_to_mX(X)
        cp_camix = cp_caco3*mX+cp_cao*(1-mX)
        return cp_camix

    def convert_X_to_mX(self,X):
        M_cao = 56e-3  # kg/mol
        M_caco3 = 100e-3  # kg/mol
        molar_mass=X*M_caco3+(1-X)*M_cao
        mX=X*M_caco3/molar_mass
        return mX

    def _cp0mass_gas(self, T, p, fluid):
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp

    def _arange(self, start, stop, step=1, endpoint=True):
        arr = np.arange(start, stop, step)

        if endpoint and arr[-1] != stop:
            arr = np.concatenate([arr, [stop]])

        return arr

