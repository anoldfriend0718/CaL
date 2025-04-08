import os
import sys
CaLRepo = os.environ.get("CaLRepo")
# print(CaLRepo)
sys.path.append(f"{CaLRepo}/utilities/")

import math
import numpy as np
import CoolProp.CoolProp as CP
from Cp0massWrapper import Cp0mass_Wrapper
from pyPinchPointAnalyzer import Pinch_point_analyzer
from pyCostEstimator import Cost_Estimator
from pyMakeUpFlowBalaner import Make_Up_Flow_Balaner

M_cao = 56e-3  # kg/mol
M_caco3 = 100e-3  # kg/mol
M_CO2 = 44e-3  # kg/mol

class CarbonatorSide(object):
    def __init__(self, parameters) -> None:
        self._pw = Cp0mass_Wrapper(parameters["flue_gas_composition"], parameters["decarbonized_rate"])
        self._cost_estimator=Cost_Estimator()

        self._flue_gas_composition = self._pw._norm_flue_gas_composition_s
        self._decarbonized_rate = parameters["decarbonized_rate"]
        self._isentropic_eff_mc = parameters["isentropic_eff_mc"]
        self._mechanical_eff = parameters["mechanical_eff"]
        self._flue_gas_pressure_loss_ratio = parameters["flue_gas_pressure_loss_ratio"]
        self._decarbon_flue_gas_pressure_loss_ratio = parameters["decarbon_flue_gas_pressure_loss_ratio"]
        self._carbonator_pressure_loss=parameters["carbonator_pressure_loss"]

        self._vol_rate_flue_gas = parameters["vol_rate_flue_gas"]  # m3/s
        self._T_flue_gas_i = parameters["T_flue_gas"]
        self._T_carb = parameters["T_carb"]

        self._cao_conversion = parameters["cao_conversion"]
        
        self._T_water_supply_in= parameters["T_water_supply_in"]
        self._p_water_supply_in = parameters["p_water_supply_in"]  
        self._T_water_prod_out = parameters["T_water_prod_out"]
        self._water_pressure_drop_rate=parameters["water_pressure_drop_rate"]
        self._water_pipe_length=parameters["water_pipe_length"]
        self._pump_hydraulic_eff=parameters["water_pump_hydraulic_efficiency"]
        self._pump_mechanical_eff=parameters["water_pump_mechanical_efficiency"]

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
        self._HTCW=parameters["HTCW"]
        self._HRCP=parameters["HRCP"]

        make_up_balaner=Make_Up_Flow_Balaner()
        make_up_results=make_up_balaner.solve(self._cao_conversion)
        self._caco3_mole_fraction_carbonator_in=make_up_results["caco3_mole_fraction_calciner_outlet"] #Y



    def solve(self, inputs):
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
        results["T_water_supply_in"]=self._T_water_supply_in
        results["T_water_prod_out"] = self._T_water_prod_out
        results["p_carb_o"] = self._p_amb/(1-self._decarbon_flue_gas_pressure_loss_ratio)  # 压力略高于大气压，需扣减压力损失
        results["p_carb_i"]=self._carbonator_pressure_loss+results["p_carb_o"] 
        results["p_carb"]=(results["p_carb_o"]+results["p_carb_i"])/2
        results["p_water_supply_in"] = self._p_water_supply_in
        results["cao_conversion"] = self._cao_conversion
        results["HTCW"]=self._HTCW
        results["HRCP"]=self._HRCP
        results = {**inputs, **results}

       # carbonator(self, Ti_flue_gas, Ti_cao, Ti_water,To_water,Tcarb, pcarb, X):
        carbonator_results = self.carbonator(inputs,
                                             self._T_water_prod_out,
                                             self._T_carb,
                                             results["p_carb"],
                                             self._cao_conversion)
        results.update(carbonator_results)
        # flue gas fan
        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        flue_gas_pi = self._p_amb
        flue_gas_po = results["p_carb_i"]/(1-self._flue_gas_pressure_loss_ratio)
        flue_gas_mass_rate = results["m_flue_gas_in"]
        flue_gas_fan_results = self.flue_gas_fan_power(self._T_flue_gas_i,
                                                            flue_gas_pi,
                                                            flue_gas_po,
                                                            flue_gas_mass_rate,
                                                            flue_gas_name)
        results.update(flue_gas_fan_results)
        if results["T_flue_gas_reactor_in"]<results["T_flue_gas_fan_out"]:
            results["is_succeed"]=0

        if results["is_succeed"]==0:
            results["hot_utility"]=1e6
            results["carb_heat_rec_eff"]=-1
            return results

        # hot water
        p_water_after_pump=self._p_water_supply_in+self._water_pipe_length*self._water_pressure_drop_rate
        results["p_water_after_pump"]=p_water_after_pump
        results["water_pump_power"]=self._water_pump_power(self._p_water_supply_in,
            p_water_after_pump,results["m_water_in"])
        h_To_water = CP.PropsSI('H', 'T', self._T_water_prod_out+273.15,
                                'P', p_water_after_pump, "REFPROP::water")
        h_Ti_water = CP.PropsSI('H', 'T', self._T_water_supply_in+273.15,
                                'P', p_water_after_pump, "REFPROP::water")
        Q_hot_water = results["m_water_in"]*(h_To_water-h_Ti_water)
        results["Q_hot_water"] = Q_hot_water
        # conveying power
        results["conveying_power"] = self.conveying_power(
            results["m_camix_in"],
            results["m_camix_out"])*(-1)
        # pinch point analysis
        pa = Pinch_point_analyzer(results)
        pa_text = pa.write_pyPinch_data_text()
        hot_util, cold_util,total_HENA = pa.solve(pa_text)
        results["pinch_analysis_text"]=pa_text
        results["total_HEN_area"]=total_HENA
        results["hot_utility"] = hot_util
        results["cold_utility"] = cold_util

        ## The storage vessels are not insulated, so the residual sensible heat goes dispersed to the environment
        ## avoiding the addition of a dedicated air cooler. (https://doi.org/10.1016/j.ecmx.2020.100039)
        ## cooling power
        ## results["cooling_power"] = self.cooling_power(cold_util)*(-1)

        # energy metrics summary
        results["is_succeed"]=1
        results["carb_auxiliary_power"] = results["conveying_power"]+ \
            results["flue_gas_fan_power"]+results["water_pump_power"]
        results["carb_heat_rec_eff"] = results["Q_hot_water"]/(-results["Q_carbonation"]+results["hot_utility"])

        # calculate investment costs
        invCosts={}
        invCosts.update(self._cost_estimator.calculate_carbonator_costs(results))
        invCosts["total"]=np.sum(list(invCosts.values()))
        invCosts["specific"]=invCosts["total"]/results["m_water_in"]
        results["invCosts"]=invCosts

        return results

    def _water_pump_power(self,pi,po,mass_rate):
        rho_water=1000
        vol_rate=mass_rate/rho_water
        dp=po-pi
        power=vol_rate*dp/self._pump_hydraulic_eff/self._pump_mechanical_eff*(-1)
        return power
    
    def flue_gas_fan_power(self, Ti, pi, po, mass_rate, fluid):
        hi = CP.PropsSI('H', 'T', Ti+273.15, 'P', pi, fluid)
        si = CP.PropsSI('S', 'T', Ti+273.15, 'P', pi, fluid)
        ho_s = CP.PropsSI('H', 'P', po, 'S', si, fluid)
        ho_c = (ho_s-hi)/self._isentropic_eff_mc+hi
        To = CP.PropsSI('T', 'P', po, 'H', ho_c, fluid)-273.15
        W = ((ho_c-hi)/self._mechanical_eff)*mass_rate
        results = {}
        results["p_flue_gas_fan_out"] = po
        results["T_flue_gas_fan_out"] = To
        results["flue_gas_fan_power"] = W*(-1)
        return results

    def carbonator(self,inputs, T_water_reactor_out, Tcarb, pcarb, X):
        if self._HTCW==0 and self._HRCP==1:
            T_flue_gas_reactor_in=inputs["T_flue_gas_reactor_in"]
            mass_water=inputs["m_water_in"]
            T_water_reactor_in=self._T_water_prod_out
        else:
            raise ValueError("only support HRCP now!")

        delta_H_Tr = self._mole_reaction_heat(Tcarb, pcarb)

        mass_co2_i = self._vol_rate_flue_gas*self._flue_gas_composition["co2"] *\
            CP.PropsSI('D', 'T', self._T_flue_gas_i+273.15, 'P', self._p_flue_gas_i, "REFPROP::co2")
        mole_co2_i = mass_co2_i/M_CO2
        mole_co2_r = mole_co2_i*self._decarbonized_rate
        mass_co2_r = mole_co2_r*M_CO2

        flue_gas_name = self._pw.get_flue_gas_refprop_name()
        mass_flue_gas = self._vol_rate_flue_gas *\
            CP.PropsSI('D', 'T', self._T_flue_gas_i+273.15, 'P', self._p_flue_gas_i, flue_gas_name)
        mass_deconbonized_flue_gas_out = mass_flue_gas-mole_co2_r*M_CO2

        X = self._cao_conversion
        Y = self._caco3_mole_fraction_carbonator_in   
        mole_camix_in=mole_co2_r/(X-Y)
        M_camix_in=Y*M_caco3+(1-Y)*M_cao  
        mass_camix_in=mole_camix_in*M_camix_in
        mole_cao_i = mole_camix_in*(1-Y)
        mass_cao_i = mole_cao_i*M_cao
        mole_caco3_i = mole_camix_in*Y
        mass_caco3_i = mole_caco3_i*M_caco3

        mole_cao_o = mole_camix_in*(1-X)
        mass_cao_o = mole_cao_o*M_cao
        mole_caco3_o = mole_camix_in*X
        mass_caco3_o = mole_caco3_o*M_caco3
        mass_camix_o=mass_caco3_o+mass_cao_o

        Q_heat = mole_co2_r*delta_H_Tr*self._carbonator_eff
        # print(f"total reaction heat released is {Q_heat/1000} kW")
  
        h_Tr_flue_gas = CP.PropsSI('H', 'T', Tcarb+273.15, 'P', pcarb, flue_gas_name)
        h_Ti_flue_gas = CP.PropsSI('H', 'T', T_flue_gas_reactor_in+273.15, 'P', pcarb, flue_gas_name)  # 等压过程
        results = {}
        results["is_succeed"]=1
        if self._HRCP==1:
            Q_camix_i=-Q_heat-mass_flue_gas*(h_Tr_flue_gas-h_Ti_flue_gas)
            cp_camix_i_Tr = self._pw.cp_camix(Tcarb,Y)
            h_Ti_camix_i=cp_camix_i_Tr*Tcarb-Q_camix_i/mass_cao_i
            if h_Ti_camix_i<=0 or h_Ti_camix_i>cp_camix_i_Tr*Tcarb:
                results["is_succeed"]=0
            T_camix_carbonator_in=self._calculate_T_camix(cp_camix_i_Tr, h_Ti_camix_i,Y)
            results["T_cao_reactor_in"]=T_camix_carbonator_in
        else:
            raise ValueError("Only support HRCP now!")

        deconbonized_flue_gas_composition = {}
        deconbonized_flue_gas_composition["co2"] = self._pw._co2_mole_frac_e
        deconbonized_flue_gas_composition["n2"] = self._pw._n2_mole_frac_e
        deconbonized_flue_gas_composition["o2"] = self._pw._o2_mole_frac_e
        results["deconbonized_flue_gas_composition"] = deconbonized_flue_gas_composition
        results["T_water_reactor_in"]=T_water_reactor_in
        results["mole_camix_in"]=mole_camix_in
        results["mole_carbonation"]=mole_co2_r
        results["m_flue_gas_in"] = mass_flue_gas
        results["m_water_in"] = mass_water
        results["m_co2_capture"] = mass_co2_r
        results["m_cao_in"] = mass_cao_i
        results["m_caco3_in"] = mass_caco3_i
        results["m_camix_in"] = mass_camix_in
        results["m_cao_out"] = mass_cao_o
        results["m_caco3_out"] = mass_caco3_o
        results["m_camix_out"] = mass_camix_o
        results["m_deconbonized_flue_gas_out"] = mass_deconbonized_flue_gas_out

        results["Q_carbonation"]=mole_co2_r*delta_H_Tr
        results["delta_H_Tcarb"]=delta_H_Tr 
        return results

    def _mole_reaction_heat(self, Tcarb, pcarb):
        Tref = 20
        cp_cao_mean_Tref_Tr = self._pw.cp0mass_mean("cao", Tref, Tcarb)
        cp_caco3_mean_Tref_Tr = self._pw.cp0mass_mean("caco3", Tref, Tcarb)

        delta_H_Tr = self._delta_H_Tref+(cp_caco3_mean_Tref_Tr*M_caco3
                                         - cp_cao_mean_Tref_Tr*M_cao)*(Tcarb-Tref)\
            - (CP.PropsSI('H', 'T', Tcarb+273.15, 'P', pcarb, "REFPROP::co2") -
               CP.PropsSI('H', 'T', Tref+273.15, 'P', pcarb, "REFPROP::co2"))*M_CO2
               
        return delta_H_Tr

    def _calculate_T_camix(self, cp_camix_Tr, h_Ti_camix,Y):
        Ti_camix_n_1=h_Ti_camix/cp_camix_Tr
        Ti_camix_n=0
        eps=1.5
        while abs(Ti_camix_n_1-Ti_camix_n)>eps:
            Ti_camix_n=Ti_camix_n_1
            cp_camix_Ti=self._pw.cp_camix(Ti_camix_n_1,Y)
            Ti_camix_n_1=h_Ti_camix/cp_camix_Ti
        Ti_camix=Ti_camix_n_1
        return Ti_camix

    def conveying_power(self, m_camix_in,m_camix_o):
        return self._convey_consumption*self._storage_carbonator_distance * \
            (m_camix_in+m_camix_o)

    # def cooling_power(self, cold_utility):
    #     return self._cooling_eff*cold_utility

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
        self._HTC_SSHX = parameters["HTC_SSHX"]
        self._HTC_SGHX = parameters["HTC_SGHX"]
        self._cooling_eff=parameters["cooling_eff"]
        self._delta_H_Tref = 178e3  # J/mole CaO

        make_up_balaner=Make_Up_Flow_Balaner()
        make_up_results=make_up_balaner.solve(self._cao_conversion)
        self.make_up_ratio=make_up_results["make_up_ratio"]
        self.calcination_fraction=make_up_results["calcination_fraction"] #fcalc
        self.caco3_mole_fraction_calciner_out=make_up_results["caco3_mole_fraction_calciner_outlet"] #Y

    def solve(self, mass_camix,mass_caco3_make_up,mfrac):
        results = {}
        calciner_results=self.calciner(mass_camix, mass_caco3_make_up,mfrac)
        results.update(calciner_results)
 
        # CO2 compression
        mass_co2_out=results["mass_co2_out"]
        T_co2_o=results["T_co2_o"]
        # print(mass_co2_out,T_co2_o,self._p_amb,self._p_co2_storage,self._T_cooling_co2,self._n_compression)
        compressor_results=self.co2_multiple_stage_compressor(mass_co2_out,
            T_co2_o,self._p_amb,self._p_co2_storage,self._T_cooling_co2,self._n_compression)
        compressor_results["cooling_power"]=self.cooling_power(compressor_results["cooling_energy"])*(-1)
        results["CO2_compression_train"]=compressor_results

        results["conveying_power"]=self.conveying_power(results["mass_camix_circulated"]+mass_caco3_make_up,results["mass_cao_out"])*(-1)

        results["calc_auxiliary_power"]=results["conveying_power"]+results["CO2_compression_train"]["cooling_power"]+\
            results["CO2_compression_train"]["compressor_power"]

        return results

    def calciner(self, mass_camix, mass_caco3_make_up,mfrac):
        # recirculated inflow mass rate is mass_Camix, in which CaCO3 ocuppied X percentage
        # make_up flow is mass_make_up which is full of fresh CaCO3

        X = self._cao_conversion
        M_camix=X*M_caco3+(1-X)*M_cao
        mole_camix=mass_camix/M_camix
        mole_camix_circulated=(1-self.make_up_ratio)*mole_camix #dispose r Camix for balance
        mass_camix_circulated=mole_camix_circulated*M_camix

        mole_caco3_make_up=mass_caco3_make_up/M_caco3
        mole_caco3_in=X*mole_camix_circulated+mole_caco3_make_up
        mole_CO2_out = mole_caco3_in*self.calcination_fraction
        mass_co2_out = mole_CO2_out*M_CO2

        mole_CaCO3_out=mole_camix*self.caco3_mole_fraction_calciner_out
        mass_caco3_out=mole_CaCO3_out*M_caco3

        mole_cao_out = mole_camix*(1-self.caco3_mole_fraction_calciner_out)
        mass_cao_out = mole_cao_out*M_cao

        mass_camix_out=mass_cao_out+mass_caco3_out

        mass_calciner_in = mass_camix_circulated+mass_caco3_make_up
        mass_camix_in_1 = mfrac*mass_calciner_in
        mass_camix_in_2 = (1-mfrac)*mass_calciner_in
        self.caco3_mole_fraction_calciner_in=mole_CO2_out/mole_cao_out

        # Hot camix outlet ~ cold camix inlet heat exchanger
        T_camix_i_1, T_camix_o,HXA_cao_camix = self.two_camix_heat_exchanger(mass_camix_in_1, mass_camix_out)
        # Hot CO2 ~ cold ca mixture heat exchanger
        T_camix_i_2, T_co2_o,HXA_co2_camix = self.co2_camix_heat_exchanger(mass_camix_in_2, mass_co2_out)
        # ca mixture mixer
        T_camix_calc_in = self.caco3_mixer(mass_camix_in_1, T_camix_i_1, mass_camix_in_2, T_camix_i_2)
        # calciner
        delta_H_r_Tcalc = self.reaction_heat()
        delta_H=self._pw.cp_camix(self._T_calc, self.caco3_mole_fraction_calciner_in)*self._T_calc - \
                self._pw.cp_camix(T_camix_calc_in, self.caco3_mole_fraction_calciner_in)*T_camix_calc_in # based on mass
        Qe = mole_CO2_out*delta_H_r_Tcalc+(mass_camix_circulated+mass_caco3_make_up)*delta_H
        Qe_loss=mole_caco3_make_up*(1-X)*delta_H_r_Tcalc
        We = Qe/self._calciner_eff
        We_loss=Qe_loss/self._calciner_eff

        results = {}
        results["p_calc"]=self._p_amb
        results["mfrac"]=mfrac
        results["calcination_fraction"]=self.calcination_fraction
        results["caco3_mole_fraction_calciner_inlet"]=self.caco3_mole_fraction_calciner_in
        results["caco3_mole_fraction_calciner_outlet"]=self.caco3_mole_fraction_calciner_out
        results["mole_camix_in"]=mole_camix
        results["mole_caco3_in"]=mole_caco3_in
        results["mole_calcination"]=mole_CO2_out
        results["mole_calcination_noused"]=mole_caco3_make_up*(1-X)
        results["mass_camix_circulated"]=mass_camix_circulated
        results["mass_caco3_make_up"]=mass_caco3_make_up
        results["mass_camix_1"]=mass_camix_in_1
        results["mass_camix_2"]=mass_camix_in_2
        results["mass_cao_out"]=mass_cao_out
        results["mass_co2_out"]=mass_co2_out
        results["mass_caco3_out"]=mass_caco3_out
        results["T_camix_o"] = T_camix_o
        results["T_co2_o"] = T_co2_o
        results["T_camix_i_1"] = T_camix_i_1
        results["T_camix_i_2"] = T_camix_i_2
        results["T_camix_reactor_in"] = T_camix_calc_in
        results["delta_H_Tcalc"] = delta_H_r_Tcalc
        results["Qe_calc"] = Qe
        results["We_calc"] = We
        results["Qe_calc_loss"] = Qe_loss
        results["We_calc_loss"] = We_loss
        results["HEN_area_cao_camix"]=HXA_cao_camix
        results["HEN_area_co2_camix"]=HXA_co2_camix
        results["total_HEN_area"]=HXA_cao_camix+HXA_co2_camix
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

    def two_camix_heat_exchanger(self, m_camix_i, m_camix_o):
        T_camix_i_hi = self._T_amb
        T_camix_o_hi = self._T_calc
        cp_camix_o_mean_Tamb_Tr = self._pw.cp_camix_mean(T_camix_i_hi, T_camix_o_hi,self.caco3_mole_fraction_calciner_out)
        cp_camix_i_mean_Tamb_Tr = self._pw.cp_camix_mean(T_camix_i_hi, T_camix_o_hi, self.caco3_mole_fraction_calciner_in)
        deltaTmin = self._deltaTmin_SSHX
        results= self.SS_heat_exchanger(m_camix_i, cp_camix_i_mean_Tamb_Tr, m_camix_o, 
            cp_camix_o_mean_Tamb_Tr,T_camix_i_hi, T_camix_o_hi, deltaTmin)
        T_camix_i_ho_n_1=results["Toc"]
        T_camix_o_ho_n_1=results["Toh"]
        T_camix_i_ho_n = 0
        T_camix_o_ho_n = 0
        eps = 2

        while abs(T_camix_i_ho_n_1-T_camix_i_ho_n) > eps or abs(T_camix_o_ho_n_1-T_camix_o_ho_n) > eps:
            T_camix_i_ho_n = T_camix_i_ho_n_1
            T_camix_o_ho_n = T_camix_o_ho_n_1
            cp_camix_i_mean = self._pw.cp_camix_mean(T_camix_i_hi, T_camix_i_ho_n, self.caco3_mole_fraction_calciner_in)
            cp_camix_o_mean = self._pw.cp_camix_mean(T_camix_o_hi, T_camix_o_ho_n,self.caco3_mole_fraction_calciner_out)
            results = self.SS_heat_exchanger(m_camix_i, cp_camix_i_mean, m_camix_o, cp_camix_o_mean,
                T_camix_i_hi, T_camix_o_hi, deltaTmin)
            T_camix_i_ho_n_1=results["Toc"]
            T_camix_o_ho_n_1=results["Toh"]
        HXA=results["HXA"]
        return T_camix_i_ho_n_1, T_camix_o_ho_n_1,HXA

    def SS_heat_exchanger(self, mc, cpc, mh, cph, Tic, Tih, deltaTmin):
        if mc*cpc > mh*cph:
            Toh = Tic+deltaTmin
            Toc = mh*cph*(Tih-Toh)/(mc*cpc)+Tic
        else:
            Toc = Tih-deltaTmin
            Toh = Tih-mc*cpc*(Toc-Tic)/(mh*cph)
        Q=mc*cpc*(Toc-Tic)
        deltaT1=(Tih-Toc)/2
        deltaT2=(Toh-Tic)/2
        lmtd=(deltaT1-deltaT2)/(math.log(deltaT1)-math.log(deltaT2))
        HXA=Q/(self._HTC_SSHX*lmtd)
        results={}
        results["Toc"]=Toc
        results["Toh"]=Toh
        results["Q"]=Q
        results["LMTD"]=lmtd
        results["HXA"]=HXA
        return results

    def co2_camix_heat_exchanger(self, m_camin, m_co2):
        T_camix_in = self._T_amb
        T_co2_in = self._T_calc
        cp_camix_mean_Tamb_Tr = self._pw.cp_camix_mean(T_camix_in, T_co2_in, self.caco3_mole_fraction_calciner_in)
        deltaTmin = self._deltaTmin_SGHX
        results = self.SG_heat_exchanger(m_camin, cp_camix_mean_Tamb_Tr, m_co2,
                        T_camix_in, T_co2_in, self._p_calc, deltaTmin)
        T_camix_out_n_1=results["Toc"]
        T_co2_out_n_1=results["Toh"]

        T_camix_out_n = 0
        T_co2_out_n = 0
        eps = 1.5

        while abs(T_camix_out_n_1-T_camix_out_n) > eps or abs(T_co2_out_n_1-T_co2_out_n) > eps:
            T_camix_out_n = T_camix_out_n_1
            T_co2_out_n = T_co2_out_n_1
            cp_camix_mean = self._pw.cp_camix_mean(T_camix_in, T_camix_out_n, self.caco3_mole_fraction_calciner_in)
            results = self.SG_heat_exchanger(m_camin, cp_camix_mean, m_co2,
                              T_camix_in, T_co2_in, self._p_calc, deltaTmin)
            T_camix_out_n_1=results["Toc"]
            T_co2_out_n_1=results["Toh"]
        HXA=results["HXA"]
        return T_camix_out_n_1, T_co2_out_n_1,HXA

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
        
        Q=mc*cpc*(Toc-Tic)
        deltaT1=Tih-Toc
        deltaT2=Toh-Tic
        lmtd=(deltaT1-deltaT2)/(math.log(deltaT1)-math.log(deltaT2))
        HXA=Q/(self._HTC_SGHX*lmtd)
        results={}
        results["Toc"]=Toc
        results["Toh"]=Toh
        results["Q"]=Q
        results["LMTD"]=lmtd
        results["HXA"]=HXA
        return results

    def caco3_mixer(self, m1, T1, m2, T2):
        Ton = 0
        Ton_1 = (m1*T1+m2*T2)/(m1+m2)
        eps = 1
        while abs(Ton_1-Ton) < eps:
            Ton = Ton_1
            cp1 = self._pw.cp_camix_mean(T1, Ton_1, self.caco3_mole_fraction_calciner_in)
            cp2 = self._pw.cp_camix_mean(T2, Ton_1, self.caco3_mole_fraction_calciner_in)
            Ton_1 = (m1*cp1*T1+m2*cp2*T2)/(m1*cp1+m2*cp2)
        return Ton_1

    def co2_multiple_stage_compressor(self,m,Ti,pi,po,Tc,n):
        pr=pow(po/pi,1./n)
        i=1
        pis=[]
        pos=[]
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
            pis.append(pi_i)
            pos.append(po_i)
            Ws.append(W_i)
            Tos.append(To_i)
            Qcs.append(Qc)
            i=i+1
            Ti_i=Tc
            pi_i=po_i
        Wt=np.sum(Ws)
        Qct=np.sum(Qcs)
        results={}
        results["compressor_Ti"]=Ti
        results["compressor_To_stages"]=Tos
        results["compressor_Tc"]=Tc
        results["compressor_pi_stages"]=pis
        results["compressor_po_stages"]=pos
        results["compressor_power_stages"]=Ws
        results["compressor_cooling_energy_stages"]=Qcs
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
