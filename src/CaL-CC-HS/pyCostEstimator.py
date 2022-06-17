
import math
import numpy as np

class Cost_Estimator(object):
    def __init__(self) -> None:
        self._Eur2Dollars={2014:39.9116/30.3287,2018:35.3758/30.1162}
        

    def solve(self,design):
        invCosts={}
        # Calciner side
        # Calciner
        Qcalc=design["calc"]["We_calc"]
        invCosts["Ccalc"]=self._cost_external_heated_calciner(Qcalc)
        # CO2 compression train
        PCCT=-design["calc"]["CO2_compression_train"]["compressor_power"]
        invCosts["CCO2ct"]=self._cost_CO2_compression_train(PCCT)
        # TODO: HEN
        # TODO: CO2 fan
        # TODO: solid conveying system
        
        # # Carbonator side 
        # # carbonator
        # # if the reaction heat is recovered at the carbonator wall
        # if design["calc"]["HTCW"]==1:
        #     Qcarb=-design["carb"]["Q_carbonation"]
        #     invCosts["Ccarb"]=self._cost_fluidized_bed_carbonator(Qcarb)
        # else:
        #     rho_CaO=3297.63  # data from Aspen
        #     V_CaO=design["carb"]["m_cao_in"]/rho_CaO
        #     invCosts["Ccarb"]=self._cost_entrained_flow_carbonator(V_CaO)
        
        Qcarb=-design["carb"]["Q_carbonation"]
        invCosts["Ccarb"]=self._cost_fluidized_bed_carbonator(Qcarb)
        # Flue gas fan 
        PFF=-design["carb"]["flue_gas_fan_power"]
        invCosts["CFF"]=self._cost_flue_gas_fan(PFF)
        # Heat tranfer network
        AHENCarb=design["carb"]["total_HEN_area"]
        pHEN=design["carb"]["p_carb"]
        invCosts["CHENCarb"]=self._cost_HEN(AHENCarb,pHEN)
        # TODO: solid conveying system

        # summary 
        invCosts["total"]=np.sum(list(invCosts.values()))

        return invCosts




        

        

    def _cost_entrained_flow_carbonator(self,Vi):
        #reference: https://doi.org/10.1016/j.ijggc.2019.01.005
        #the referred Year: 2014
        #the given Money unit: EUR 
        #Vi unit: m3/s
        year=2014
        EU=85.9*math.pow(Vi,0.5)*self._Eur2Dollars[year] #unit M$
        return EU

    def _cost_fluidized_bed_carbonator(self,Qc):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: EUR
        #Qc unit: W
        year=2018
        FC=16591*math.pow(Qc/1000,0.67)*self._Eur2Dollars[year]/1e6  #unit M$
        return FC
    
    def _cost_external_heated_calciner(self,Qc):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: EUR
        #Qc unit: W
        year=2018
        HC=13140*math.pow(Qc/1000,0.67)*self._Eur2Dollars[year]/1e6  #unit M$
        return HC

    def _cost_CO2_compression_train(self,Pt):
        #reference:https://doi.org/10.1016/j.ecmx.2020.100038
        #the referred year: 2018
        #the given Money unit: USD
        #Pt unit: W
        year=2018
        CCT=7331*math.pow(Pt/1e3,0.7865)/1e6  #unit M$
        return CCT

    def _cost_flue_gas_fan(self,Pi):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: EUR
        #Pi unit: W
        year=2018
        CF=103193*math.pow(Pi/445/1000,0.67)*self._Eur2Dollars[year]/1e6  #unit M$ 
        return CF

    def _cost_HEN(self,At,p=1e5):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: EUR
        #Ai unit: m2,p: Pa
        year=2018
        p=p/1e5
        CHEN=2546.9*math.pow(At,0.67)*math.pow(p,0.28)*self._Eur2Dollars[year]/1e6  #unit M$ 
        return CHEN





