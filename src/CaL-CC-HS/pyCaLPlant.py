import os
import sys
# CaLRepo = os.environ.get("CaLRepo")
CaLRepo ='/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL'
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")

import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
from pyPinch import PyPinch
from scipy import interpolate



class Cp0mass_Wrapper(object):
    def __init__(self,flue_gas_composistion,deconbonized_rate) -> None:
        self.f_cao, self.f_caco3 = self._read_solid_properties()
        total_mole_frac_s=flue_gas_composistion["co2"]+ \
            flue_gas_composistion["n2"]+\
            flue_gas_composistion["o2"]   
        self._co2_mole_frac_s=flue_gas_composistion["co2"]/total_mole_frac_s
        self._n2_mole_frac_s=flue_gas_composistion["n2"]/total_mole_frac_s
        self._o2_mole_frac_s=flue_gas_composistion["o2"]/total_mole_frac_s

        total_mole_frac_e=total_mole_frac_s-flue_gas_composistion["co2"]*deconbonized_rate
        self._co2_mole_frac_e=flue_gas_composistion["co2"]*(1-deconbonized_rate)/total_mole_frac_e
        self._n2_mole_frac_e=flue_gas_composistion["n2"]/total_mole_frac_e
        self._o2_mole_frac_e=flue_gas_composistion["o2"]/total_mole_frac_e

    # T: oC
    # cp: J/kg/K
    def cp0mass(self,material,T,p=1e5):
        material = material.lower()
        if material=="flue_gas":
            cp=self._cp0mass_flue_gas(T,p)
        elif material=="decarbonized_flue_gas":
            cp=self._cp0mass_decarbonized_flue_gas(T,p)
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

    def _cp0mass_flue_gas(self,T, p):
        fluid = self.get_flue_gas_refprop_name()
        cp=self._cp0mass_gas(T, p,fluid)
        return cp

    def get_flue_gas_refprop_name(self):
        fluid = f'REFPROP::co2[{self._co2_mole_frac_s}]&\
                  nitrogen[{self._n2_mole_frac_s}]&\
                  oxygen[{self._o2_mole_frac_s}]'.replace(" ","")
                  
        return fluid

    def _cp0mass_decarbonized_flue_gas(self,T, p):
        fluid = self.get_decarbonized_flue_gas_refprop_name()
        cp=self._cp0mass_gas(T, p,fluid)
        return cp

    def get_decarbonized_flue_gas_refprop_name(self):
        fluid = f'REFPROP::co2[{self._co2_mole_frac_e}]&\
            nitrogen[{self._n2_mole_frac_e}]&\
            oxygen[{self._o2_mole_frac_e}]'.replace(" ","")
            
        return fluid
    
    def _cp0mass_gas(self,T, p,fluid):
        cp = CP.PropsSI('C', 'T', T+273.15, 'P', p, fluid)  # J/kg/k
        return cp


    def _arange(self,start, stop, step=1, endpoint=True):
        arr = np.arange(start, stop, step)

        if endpoint and arr[-1] != stop:
            arr = np.concatenate([arr, [stop]])

        return arr
    

class Pinch_point_analyzer(object):
    def __init__(self, input) -> None:
        self._m_caco3 = input["m_caco3_out"]
        self._m_cao_unr = input["m_cao_unr_out"]
        self._m_cao_i = input["m_cao_in"]
        self._m_decarbonized_flue_gas=input["m_deconbonized_flue_gas_out"]
        self._m_flue_gas_in=input["m_flue_gas_in"]
        self._m_water_in=input["m_water_in"]

        self._T_carb = input["T_carb"]
        self._T_amb = input["T_amb"]
        self._T_flue_gas_in=input["T_flue_gas_in"]
        self._T_flue_gas_reactor_in = input["T_flue_gas_reactor_in"]
        self._T_cao_reactor_in = input["T_cao_reactor_in"]
        self._T_water_reactor_in = input["T_water_reactor_in"]
        self._T_decarbonized_flue_gas_out=input["T_decarbonized_flue_gas_out"]
        self._T_delta_pinch = input["T_delta_pinch"]
        
        self._p_carb = input["p_carb"]
        
        self._flue_gas_composition=input["flue_gas_composition"]
        self._deconbonized_rate=input["deconbonized_rate"]

        self._pw = Cp0mass_Wrapper(self._flue_gas_composition,self._deconbonized_rate)
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

        ## H2: decarbonized flue gas 
        pinch_point_data["TSUPPLY"]["H_flue_gas_decarb"] = self._T_carb
        pinch_point_data["TTARGET"]["H_flue_gas_decarb"] = self._T_decarbonized_flue_gas_out
        pinch_point_data["FLOWRATE"]["H_flue_gas_decarb"] = self._m_decarbonized_flue_gas
        decarb_flue_gas_name=self._pw.get_decarbonized_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["H_flue_gas_decarb"] = self._m_decarbonized_flue_gas * \
            (CP.PropsSI('H', 'T', self._T_carb+273.15,
                        'P', self._p_carb, decarb_flue_gas_name) -
             CP.PropsSI('H', 'T', self._T_decarbonized_flue_gas_out + 273.15,
                        'P', self._p_carb, decarb_flue_gas_name))
        pinch_point_data["CP"]["H_flue_gas_decarb"] = pinch_point_data["ENERGY"]["H_flue_gas_decarb"] / \
            (self._T_carb-self._T_decarbonized_flue_gas_out)

        ## C1: flue gas 
        pinch_point_data["TSUPPLY"]["C_flue_gas"] = self._T_flue_gas_in
        pinch_point_data["TTARGET"]["C_flue_gas"] = self._T_flue_gas_reactor_in
        pinch_point_data["FLOWRATE"]["C_flue_gas"] = self._m_flue_gas_in
        flue_gas_name=self._pw.get_flue_gas_refprop_name()
        pinch_point_data["ENERGY"]["C_flue_gas"] = self._m_flue_gas_in * \
            (CP.PropsSI('H', 'T', self._T_flue_gas_reactor_in+273.15,
                        'P', self._p_carb, flue_gas_name) -
             CP.PropsSI('H', 'T', self._T_flue_gas_in+273.15,
                        'P', self._p_carb, flue_gas_name))
        pinch_point_data["CP"]["C_flue_gas"] = pinch_point_data["ENERGY"]["C_flue_gas"] / \
            (self._T_flue_gas_reactor_in-self._T_flue_gas_in)

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
                        'P', self._p_carb, "REFPROP::water") -
             CP.PropsSI('H', 'T', self._T_amb+273.15,
                        'P', self._p_carb, "REFPROP::water"))
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
        text=df.to_string(index=False)
        text=f'Tmin {self._T_delta_pinch} \n'+text
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
        self._flue_gas_composition=parameters["flue_gas_composition"]
        self._decarbonized_rate=parameters["decarbonized_rate"]
        self._vol_rate_flue_gas=parameters["vol_rate_flue_gas"] # m3/s
        self._T_flue_gas=parameters["T_flue_gas"]
        self._T_carb=parameters["T_carb"]
        self._p_carb=parameters["p_carb"]
        self._cao_conversion=parameters["cao_conversion"]
        self._T_water_reactor_out=parameters["T_water_reactor_out"]

        self._carbonator_eff = parameters["carbonator_eff"]
        self._convey_consumption = parameters["convey_consumption"]
        self._storage_carbonator_distance = parameters["storage_carbonator_distance"]
        self._cooling_eff = parameters["cooling_eff "]
        self._delta_T_pinch = parameters["delta_T_pinch"]
        self._T_amb = parameters["T_amb"]
        self._delta_H_Tref = -178e3  # J/mole CaO
        self._delta_h=3178.6*1000 # J/kg CaO
        
        self._pw = Cp0mass_Wrapper(self._flue_gas_composition,self._decarbonized_rate)

    def solve(self,input):
        results={}
        results["flue_gas_composition"]=self._flue_gas_composition
        results["deconbonized_rate"]=self._decarbonized_rate
        results["T_carb"]=self._T_carb
        results["T_amb"]=self._T_amb
        results["T_delta_pinch"]=self._delta_T_pinch
        results["T_flue_gas_in"]=self._T_flue_gas
        results["T_decarbonized_flue_gas_out"]=self._T_flue_gas
        results["T_water_reactor_out"]=self._T_water_reactor_out
        results["p_carb"]=self._p_carb
        results["cao_conversion"]=self._cao_conversion
        results = {**input, **results}
        # carbonator(self, Ti_flue_gas, Ti_cao, Ti_water,To_water,Tcarb, pcarb, X):
        carbonator_results = self.carbonator(input["T_flue_gas_reactor_in"],
                                             input["T_cao_reactor_in"],
                                             input["T_water_reactor_in"],
                                             self._T_water_reactor_out,
                                             self._T_carb,
                                             self._p_carb,
                                             self._cao_conversion)
        results.update(carbonator_results) 

        # Q hot water 
        h_To_water = CP.PropsSI('H', 'T', self._T_water_reactor_out+273.15,
                                'P', self._p_carb, "REFPROP::water")
        h_Ti_water = CP.PropsSI('H', 'T', self._T_amb+273.15,
                                'P', self._p_carb, "REFPROP::water")
        Q_hot_water=results["m_water_in"]*(h_To_water-h_Ti_water)
        results["Q_hot_water"]=Q_hot_water
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["m_cao_in"],
            results["m_cao_unr_out"],
            results["m_caco3_out"])*(-1)
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

        results["plant_eff"] = results["Q_hot_water"]/((results["m_cao_in"] -
                            results["m_cao_unr_out"])*self._delta_h)

        return results

  
    def carbonator(self, Ti_flue_gas, Ti_cao, Ti_water,To_water,Tcarb, pcarb, X):
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
        # print(f"specific reaction heat released is {delta_H_Tr/1000} KJ/mol")

        mass_co2_i=self._vol_rate_flue_gas*self._flue_gas_composition["co2"]*\
                    CP.PropsSI('D', 'T', self._T_flue_gas+273.15, 'P', pcarb, "REFPROP::co2")
        mole_co2_i=mass_co2_i/M_CO2
        mole_co2_r=mole_co2_i*self._decarbonized_rate
        mass_co2_r=mole_co2_r*M_CO2

        flue_gas_name=self._pw.get_flue_gas_refprop_name()
        mass_flue_gas=self._vol_rate_flue_gas*\
            CP.PropsSI('D', 'T', self._T_flue_gas+273.15, 'P', pcarb, flue_gas_name)
        mass_deconbonized_flue_gas_out=mass_flue_gas-mole_co2_r*M_CO2

        mole_cao_r=mole_co2_r
        mass_cao_r = mole_cao_r*M_cao
        mole_cao_i=mole_cao_r/X
        mass_cao_i=mole_cao_i*M_cao

        mole_caco3_o = mole_cao_r
        mass_caco3_o = mole_caco3_o*M_caco3

        Q_heat = mole_cao_r*delta_H_Tr*self._carbonator_eff
        # print(f"total reaction heat released is {Q_heat/1000} kW")

        cp_cao_Tr = self._pw.cp0mass("cao", Tcarb)
        cp_cao_Ti = self._pw.cp0mass("cao", Ti_cao)

        h_Tr_flue_gas = CP.PropsSI('H', 'T', Tcarb+273.15, 'P', pcarb, flue_gas_name)
        h_Ti_flue_gas = CP.PropsSI('H', 'T', Ti_flue_gas+273.15, 'P', pcarb, flue_gas_name)
        
        h_To_water = CP.PropsSI('H', 'T', To_water+273.15, 'P', pcarb, "REFPROP::water")
        h_Ti_water = CP.PropsSI('H', 'T', Ti_water+273.15, 'P', pcarb, "REFPROP::water")
        mass_water = (-Q_heat-mass_cao_i*(cp_cao_Tr*Tcarb-cp_cao_Ti*Ti_cao) -
                      mass_flue_gas*(h_Tr_flue_gas-h_Ti_flue_gas)) /\
                     (h_To_water-h_Ti_water)
        # print(f"mass flow rate of water is {mass_water} kg/s ")
        
        results = {}
        deconbonized_flue_gas_composition={}
        deconbonized_flue_gas_composition["co2"]=self._pw._co2_mole_frac_e
        deconbonized_flue_gas_composition["n2"]=self._pw._n2_mole_frac_e
        deconbonized_flue_gas_composition["o2"]=self._pw._o2_mole_frac_e
        results["deconbonized_flue_gas_composition"]=deconbonized_flue_gas_composition

        results["m_flue_gas_in"]=mass_flue_gas
        results["m_water_in"] = mass_water
        results["m_co2_capture"]=mass_co2_r
        results["m_cao_in"] = mass_cao_i
        results["m_deconbonized_flue_gas_out"]=mass_deconbonized_flue_gas_out
        results["m_cao_unr_out"] = mass_cao_i-mass_cao_r
        results["m_caco3_out"] = mass_caco3_o
        return results

    def conveying_power(self, m_cao_i, m_cao_unr, m_caco3_o):
        return self._convey_consumption*self._storage_carbonator_distance * \
            (m_cao_i+m_cao_unr+m_caco3_o)

    def cooling_power(self, cold_utility):
        return self._cooling_eff*cold_utility




