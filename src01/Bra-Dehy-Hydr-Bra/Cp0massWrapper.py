
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
class Cp0mass_Wrapper(object):
    def __init__(self,flue_gas_composistion) -> None:
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

        self.f_cao, self.f_caoh2,self.f_caco3= self._read_solid_properties()
    # T: oC
    # cp: J/kg/K
    def cp0mass(self, material, T, p=1e5):
        material = material.lower()
        if material == "water":
            cp = self._cp0mass_water(T, p)
        elif material == "flue_gas":
            cp = self._cp0mass_flue_gas(T, p)
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
    def get_flue_gas_refprop_name(self):
        fluid = f'REFPROP::co2[{self._co2_mole_frac}]&\
                  nitrogen[{self._n2_mole_frac}]&\
                  oxygen[{self._o2_mole_frac}]'.replace(" ", "")

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
    parameters["flue_gas_composition"] = flue_gas_composistion
    pw = Cp0mass_Wrapper(parameters["flue_gas_composition"])
    cp_cao_o = pw.cp0mass_mean("cao", 525, 20)*M_cao*505
    cp_caoh2_o = pw.cp0mass_mean("caoh2", 525, 20)*M_caoh2*505
    h_steam_in = CP.PropsSI('H', 'T', 101+273.15, 'P', 101325, "REFPROP::water")
    h_steam_out = CP.PropsSI('H', 'T', 525+273.15, 'P', 101325, "REFPROP::water")
    b=(h_steam_out-h_steam_in)*M_H2O
    
    fluid=pw.get_flue_gas_refprop_name()
    h_flue_gas_in = CP.PropsSI('H', 'T', 400+273.15, 'P', 101325 , fluid)
    h_flue_gas_out = CP.PropsSI('H', 'T', 160+273.15, 'P', 101325 , fluid)#J/(kg)
    D1 = CP.PropsSI('D', 'T', 400+273.15, 'P', 101325 , fluid)
    D2 = CP.PropsSI('D', 'T', 160+273.15, 'P', 101325 , fluid)
    D3 = CP.PropsSI('D', 'T', 20+273.15, 'P', 101325 , fluid)#(kg/m3)
    A1 = (h_flue_gas_in-h_flue_gas_out)*D1*10000/1000/3600#kW
    A2 = (h_flue_gas_in-h_flue_gas_out)*D2*10000/1000/3600
    A3 = (h_flue_gas_in-h_flue_gas_out)*D3*10000/1000/3600

    #S_steam_in = CP.PropsSI('S', 'T', 20+273.15, 'P', 101325, "REFPROP::water")
    #S_steam_out = CP.PropsSI('S', 'T', 525+273.15, 'P', 101325, "REFPROP::water")

    #H_steam_in = CP.PropsSI('H', 'T', 20+273.15, 'P', 101325, "REFPROP::water")
    #H1 = CP.PropsSI('H', 'T', 98+273.15, 'P', 101325, "REFPROP::water")
    #H2 = CP.PropsSI('H', 'T', 102+273.15, 'P', 101325, "REFPROP::water")
    #H_steam_out = CP.PropsSI('H', 'T', 525+273.15, 'P', 101325, "REFPROP::water")
    #E_h201 = (H_steam_out-293.15*S_steam_out)-(H_steam_in-293.15*S_steam_in) #J/kg
    #F1 = (78-293.15*math.log((98+273.15)/(20+273.15)))/78
    #F2 = (423-293.15*math.log((525+273.15)/(102+273.15)))/423
    #F3 = 1-293.15/(273.15+100)
    #E_h202 = (H1-H_steam_in)*F1+(H_steam_out-H2)*F2+(H2-H1)*F3
   # print(cp_cao_o,cp_caoh2_o,b,A)
    #print(A1,A2,A3)
    C_cao = 112396 #J/mole
    C_caoh2 = 54593 
    F1 = (445-293.15*math.log((465+273.15)/(20+273.15)))/445
    E_cao_o = pw.cp0mass_mean("cao", 465, 20)*445*M_cao*15.140556*F1#J/(kg K) *K *kg/mole *mole/s=j/s
    E_caoh2_o = pw.cp0mass_mean("caoh2", 465, 20)*M_caoh2*445*15.140556*F1
    S_steam_in = CP.PropsSI('S', 'T', 20+273.15, 'P', 101325, "REFPROP::water")
    S_steam_out = CP.PropsSI('S', 'T', 465+273.15, 'P', 101325, "REFPROP::water")
    H_steam_in = CP.PropsSI('H', 'T', 20+273.15, 'P', 101325, "REFPROP::water")
    H_steam_out = CP.PropsSI('H', 'T', 465+273.15, 'P', 101325, "REFPROP::water")
    E_h201 = ((H_steam_out-293.15*S_steam_out)-(H_steam_in-293.15*S_steam_in))*M_H2O*15.140556 #J/kg 
    exergy_of_Reac = E_cao_o + C_cao*15.140556 +E_h201-C_caoh2*15.140556-E_caoh2_o
    #print(C_cao*15.140556,E_cao_o,C_caoh2*15.140556,E_caoh2_o,E_h201,exergy_of_Reac)

    T1 = 362.9
    P1 = 7.5e6
    s0 = CP.PropsSI('S', 'T', 20+273.15, 'P', 1e6, "REFPROP::co2")
    h0 = CP.PropsSI('H', 'T', 20+273.15, 'P', 1e6, "REFPROP::co2")
    s1 =  CP.PropsSI('S', 'T', T1+273.15, 'P', P1, "REFPROP::co2")
    h1 = CP.PropsSI('H', 'T', T1+273.15, 'P', P1, "REFPROP::co2")
    a = (h1-293.15*s1)-(h0-293.15*s0) #j/kg
    print(a)


    h_hsteam_in =  CP.PropsSI('H', 'T', 105+273.15, 'P', 101325, "REFPROP::water")
    h_hsteam_out =  CP.PropsSI('H', 'T', 95+273.15, 'P', 101325, "REFPROP::water")
    h_lwater_in =  CP.PropsSI('H', 'T', 60+273.15, 'P', 101325, "REFPROP::water")
    h_lwater_out =  CP.PropsSI('H', 'T', 85+273.15, 'P', 101325, "REFPROP::water")
        
    m_l=(h_hsteam_in-h_hsteam_out)*0.96/(h_lwater_out-h_lwater_in)
    pc_lost = (h_hsteam_in-h_hsteam_out)*0.04
    a= 1 * \
            (CP.PropsSI('H', 'T', 105+273.15,
                        'P', 101325, "water") -
             CP.PropsSI('H', 'T', 95+273.15,
                        'P', 101325, "water"))
    print(a,(h_hsteam_in-h_hsteam_out),h_lwater_out-h_lwater_in)

