
import os
import sys
import math
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo = '/home/zyq0416/workspace/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")
from scipy import interpolate
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol
M_caco3 = 100e-3  # kg/mol
M_H2O = 18e-3  # kg/mol
class Cp0mass_Wrapper(object) :
    def __init__(self,flue_gas_composistion,flue_gas_composistion2) -> None:
        total_mole_frac = flue_gas_composistion["co2"] + \
            flue_gas_composistion["n2"] +\
            flue_gas_composistion["o2"]
        self._co2_mole_frac = flue_gas_composistion["co2"]/total_mole_frac
        self._n2_mole_frac = flue_gas_composistion["n2"]/total_mole_frac
        self._o2_mole_frac = flue_gas_composistion["o2"]/total_mole_frac
        norm_flue_gas_composition={}
        norm_flue_gas_composition["co2"] = self._co2_mole_frac
        norm_flue_gas_composition["n2"] = self._n2_mole_frac
        norm_flue_gas_composition["o2"] = self._o2_mole_frac
        self._norm_flue_gas_composition=norm_flue_gas_composition

        total_mole_frac2 = flue_gas_composistion2["co2"] + \
            flue_gas_composistion2["n2"] +\
            flue_gas_composistion2["o2"]
        self._co2_mole_frac2 = flue_gas_composistion2["co2"]/total_mole_frac2
        self._n2_mole_frac2 = flue_gas_composistion2["n2"]/total_mole_frac2
        self._o2_mole_frac2 = flue_gas_composistion2["o2"]/total_mole_frac2
        norm_flue_gas_composition2={}
        norm_flue_gas_composition2["co2"] = self._co2_mole_frac2
        norm_flue_gas_composition2["n2"] = self._n2_mole_frac2
        norm_flue_gas_composition2["o2"] = self._o2_mole_frac2
        self._norm_flue_gas_composition2=norm_flue_gas_composition2

        self.f_cao, self.f_caoh2,self.f_caco3= self._read_solid_properties()
    # T: oC
    # cp: J/kg/K
    def cp0mass(self, material, T, p=1e5):
        material = material.lower()
        if material == "water":
            cp = self._cp0mass_water(T, p)
        elif material == "flue_gas":
            cp = self._cp0mass_flue_gas(T, p)
        elif material == "flue_gas2":
            cp = self._cp0mass_flue_gas2(T, p)
        elif material == "co2":
            cp = self._cp0mass_co2(T, p)
        elif material == "cao":
            cp = self._cp0mass_cao(T)
        elif material == "caoh2":
            cp = self._cp0mass_caoh2(T)
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
#想改成基于焓值表计算的
    def _read_solid_properties(self):
        data_csv = f"{CaLRepo}/data/cpmass_cao_caoh2.csv"
        df = pd.read_csv(data_csv)
        f_cao = interpolate.interp1d(df["TEMP"], df["CAO"],fill_value="extrapolate")
        f_caoh2 = interpolate.interp1d(df["TEMP"], df["CAOH2"],fill_value="extrapolate")
        f_caco3 = interpolate.interp1d(df["TEMP"], df["CACO3"],fill_value="extrapolate")
        return f_cao, f_caoh2,f_caco3

    def _cp0mass_cao(self, T):
        cp = float(self.f_cao(T))
        return cp

    def _cp0mass_caoh2(self, T):
        cp = float(self.f_caoh2(T))
        return cp
    def _cp0mass_caco3(self, T):
        cp = float(self.f_caco3(T))
        return cp

    def _cp0mass_co2(self, T, p):
        fluid = 'REFPROP::co2'
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp
    def _cp0mass_water(self, T, p):
        fluid = 'REFPROP::water'
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp
    
    def _cp0mass_flue_gas(self, T, p):
        fluid = self.get_flue_gas_refprop_name()
        cp = self._cp0mass_gas(T, p, fluid)
        return cp
    def _cp0mass_flue_gas2(self, T, p):
        fluid = self.get_flue_gas_refprop_name2()
        cp = self._cp0mass_gas(T, p, fluid)
        return cp
    def get_flue_gas_refprop_name(self):
        fluid = f'REFPROP::co2[{self._co2_mole_frac}]&\
                  nitrogen[{self._n2_mole_frac}]&\
                  oxygen[{self._o2_mole_frac}]'.replace(" ", "")

        return fluid
    def get_flue_gas_refprop_name2(self):
        fluid = f'REFPROP::co2[{self._co2_mole_frac2}]&\
                  nitrogen[{self._n2_mole_frac2}]&\
                  oxygen[{self._o2_mole_frac2}]'.replace(" ", "")

        return fluid

    def cp_camix_mean(self, T1, T2, X,Y=1):
        cp_cao_mean = self.cp0mass_mean("cao", T1, T2)
        cp_caoh2_mean = self.cp0mass_mean("caoh2", T1, T2)
        cp_caco3_mean = self.cp0mass_mean("caco3", T1, T2)
        mX=self.convert_X_to_mX(X)
        cp_camix_mean = (cp_caoh2_mean*mX+cp_cao_mean*(1-mX))*Y+cp_caco3_mean*(1-Y)
        return cp_camix_mean
    def cp_camix_mean_Ci(self, T1, T2, Y):
        cp_cao_mean = self.cp0mass_mean("cao", T1, T2)
        cp_caoh2_mean = self.cp0mass_mean("caoh2", T1, T2)
        cp_caco3_mean = self.cp0mass_mean("caco3", T1, T2)
        a=M_caoh2
        b=M_caoh2/Y-M_caoh2

        cp_camix_mean = (cp_caoh2_mean*a+cp_caco3_mean*b)/(a+b)
        return cp_camix_mean
    def cp_camix_mean_Co(self, T1, T2, X,Y):
        cp_cao_mean = self.cp0mass_mean("cao", T1, T2)
        cp_caoh2_mean = self.cp0mass_mean("caoh2", T1, T2)
        cp_caco3_mean = self.cp0mass_mean("caco3", T1, T2)
        a=M_cao*X 
        b=M_caoh2*(1-X)
        c=(M_caoh2/Y)-M_caoh2
        cp_camix_mean = (cp_caoh2_mean*b+cp_cao_mean*a+cp_caco3_mean*c)/(a+b+c)
        return cp_camix_mean
    


    def cp_camix(self, T1, X,Y=1):
        cp_cao = self.cp0mass("cao", T1)
        cp_caoh2 = self.cp0mass("caoh2", T1)
        cp_caco3 = self.cp0mass("caco3",T1)
        mX=self.convert_X_to_mX(X)
        cp_camix = (cp_caoh2*mX+cp_cao*(1-mX))*Y+cp_caco3*(1-Y)
        return cp_camix

    def convert_X_to_mX(self,X):
        M_cao = 56e-3  # kg/mol
        M_caoh2 = 74e-3  # kg/mol
        molar_mass=X*M_caoh2+(1-X)*M_cao
        mX=X*M_caoh2/molar_mass
        return mX

    def _cp0mass_gas(self, T, p, fluid):
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp

    def _arange(self, start, stop, step=1, endpoint=True):
        arr = np.arange(start, stop, step)

        if endpoint and arr[-1] != stop:
            arr = np.concatenate([arr, [stop]])

        return arr
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
    pw = Cp0mass_Wrapper(parameters["flue_gas_composition"],parameters["flue_gas_composition2"])
    Q_cao = pw.cp0mass_mean("cao", 900, 20)*M_cao*880
    Q_caco3 = pw.cp0mass_mean("caco3", 900, 20)*M_caco3*880  #j/mol
    gcpw=(900-20-293.15*math.log((900+273.15)/293.15))/(900-20)
    exergy =Q_cao*gcpw
    print(Q_cao,exergy)

