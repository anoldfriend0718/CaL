
import math
import numpy as np

## construct cost indictor
construct_labour_cost_indictor=0.5
engineering_project_cost_indictor=0.35
TASC_multiplier=1.13
piping_integration_cost_indictor=0.05
M_cao = 56e-3  # kg/mol
M_caco3 = 100e-3  # kg/mol

class Cost_Estimator(object):
    def __init__(self) -> None:
        self._Eur2Dollars={2014:39.9116/30.3287,2018:35.3758/30.1162}
        self._Dollars2RMB={2020:6.901}
        self._CEPCIs={2001:394.3, 2014:576.1,2018:603.1,2020:596.2}
        self._target_year=2020
        
    def solve(self,design,energy_analysis_results,economic_inputs):
        ## construct costs
        construct_costs={}
        equipment_costs=self.calculate_equipment_costs(design,economic_inputs)
        construct_costs["equipment"]=equipment_costs
        construct_costs["installation"]=piping_integration_cost_indictor*equipment_costs["total"]
        construct_costs["labour"]=construct_labour_cost_indictor*(equipment_costs["total"]+construct_costs["installation"])
        construct_costs["engineering&project"]=engineering_project_cost_indictor*(equipment_costs["total"]+construct_costs["installation"])
        construct_costs["initial_material"]=self.calculate_material_costs(design,economic_inputs)
        construct_costs["total as-spent"]=TASC_multiplier*(equipment_costs["total"]+construct_costs["installation"]+
            construct_costs["labour"]+construct_costs["engineering&project"])+construct_costs["initial_material"]

        ## operational costs
        operation_costs={}
        make_up_limestone_percentage=design["calc"]["make-up_limestone_percentage"]
        operation_costs["make-up_limestone"]=design["carb"]["mole_camix_in"]*make_up_limestone_percentage \
            *M_caco3*economic_inputs["operation_hours"]*3.6*economic_inputs["limestone_price"]/1e6
        operation_costs["labour"]=economic_inputs["operation_labour_cost_indictor"]*construct_costs["total as-spent"]
        operation_costs["maintain"]=economic_inputs["maintain_cost_indictor"]*construct_costs["total as-spent"]
        operation_costs["electricity"]=economic_inputs["elec_price"]*(energy_analysis_results["total_power"])/1000* \
            economic_inputs["operation_hours"]/1e6
        operation_costs["total as-spent"]=operation_costs["labour"]+operation_costs["maintain"]+\
            operation_costs["electricity"]+operation_costs["make-up_limestone"]
  
        ## compose results
        investment_costs={}
        investment_costs["construction"]=construct_costs
        investment_costs["operation"]=operation_costs
    
        return investment_costs

    def calculate_material_costs(self,design,economic_inputs):
        CaCO3_storage_mass_circulated=design["calc"]["CaCO3_storage_mass_circulated"]
        CaCO3_storage_mass_make_up=design["calc"]["CaCO3_storage_mass_make_up"]
        total_CaCO3_storage_mass_ton=(CaCO3_storage_mass_circulated+CaCO3_storage_mass_make_up)/1000
        initial_material_cost=economic_inputs["limestone_price"]*total_CaCO3_storage_mass_ton/1e6
        return initial_material_cost

    def calculate_equipment_costs(self,design,economic_inputs):
        equipment_costs={}
        equipment_costs.update(self.calculate_calciner_costs(design["calc"],economic_inputs))
        equipment_costs.update(self.calculate_carbonator_costs(design["carb"]))
        equipment_costs["total"]=np.sum(list(equipment_costs.values()))
        return equipment_costs

    def calculate_carbonator_costs(self,carb_design):
        invCosts={}
        # # carbonator
        # # if the reaction heat is recovered at the carbonator wall
        if carb_design["HTCW"]==1:
            Qcarb=-carb_design["Q_carbonation"]
            invCosts["Ccarb"]=self._cost_fluidized_bed_carbonator(Qcarb)
        else:
            print
            v_flue_gas=carb_design["vol_rate_flue_gas"]
            invCosts["Ccarb"]=self._cost_entrained_flow_carbonator(v_flue_gas)
        
        # Flue gas fan 
        PFF=-carb_design["flue_gas_fan_power"]
        invCosts["CFF"]=self._cost_flue_gas_fan(PFF)
        # Heat tranfer network
        AHENCarb=carb_design["total_HEN_area"]
        pHENCarb=carb_design["p_carb"]
        invCosts["CHENCarb"]=self._cost_HEN(AHENCarb,pHENCarb)
        #water pump
        WPump=-carb_design["water_pump_power"]
        invCosts["CWaterPump"]=self._cost_pump(WPump)

        # TODO: solid conveying system
        return invCosts


    def calculate_calciner_costs(self,calc_design,economic_inputs):
        invCosts={}
        # Calciner
        Qcalc=calc_design["We_calc"]
        invCosts["Ccalc"]=self._cost_rotary_kiln_calciner(Qcalc)*economic_inputs["calciner_cost_factor"]
        # CO2 compression train
        PCCT=-calc_design["CO2_compression_train"]["compressor_power"]
        invCosts["CCO2ct"]=self._cost_CO2_compression_train(PCCT)
 
        # Heat tranfer network
        AHENCalc=calc_design["total_HEN_area"]
        pHENCalc=calc_design["p_calc"]
        invCosts["CHENCalc"]=self._cost_HEN(AHENCalc,pHENCalc)

        # TODO: solid conveying system
        return invCosts
    
    def _cost_entrained_flow_carbonator(self,Vi):
        #reference: https://doi.org/10.1016/j.ijggc.2019.01.005, https://doi.org/10.1016/j.rser.2019.109252
        #the referred Year: 2014
        #the given Money unit: EUR 
        #Vi unit: m3/s
        year=2014
        EU_Dollars=85.9*1e-3*math.pow(Vi,0.5)*self._Eur2Dollars[year] #unit Dollars
        EU=self._convert_to_RMB_in_target_year(EU_Dollars,year)
        return EU

    def _cost_fluidized_bed_carbonator(self,Qc):
        #reference:https://doi.org/10.1016/j.ijggc.2019.01.005
        #the referred year: 2014
        #the given Money unit: EUR
        #Qc unit: W
        year=2014
        FC_Dollars=(0.217*math.pow(Qc/1e6,0.65)+3.83)*self._Eur2Dollars[year]  #unit M$
        FC=self._convert_to_RMB_in_target_year(FC_Dollars,year)
        return FC
    
    def _cost_oxy_cfb_calciner(self,Qc):
        #reference:https://doi.org/10.1016/j.ecmx.2020.100039
        #the referred year: 2018
        #the given Money unit: EUR
        #Qc unit: W
        year=2018
        HC_Dollars=13140*math.pow(Qc/1000,0.67)*self._Eur2Dollars[year]/1e6  #unit M$
        HC=self._convert_to_RMB_in_target_year(HC_Dollars,year)
        return HC

    def _cost_rotary_kiln_calciner(self,Qc):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: Dollar
        #Qc unit: W
        year=2018
        HC_Dollars=533394*math.pow(Qc/1000/293,0.48)/1e6  #unit M$
        HC=self._convert_to_RMB_in_target_year(HC_Dollars,year)
        return HC

    def _cost_CO2_compression_train(self,Pt):
        #reference:https://doi.org/10.1016/j.ecmx.2020.100039
        #the referred year: 2018
        #the given Money unit: USD
        #Pt unit: W
        year=2018
        CCT_Dollars=7331*math.pow(Pt/1e3,0.7865)/1e6  #unit M$
        CCT=self._convert_to_RMB_in_target_year(CCT_Dollars,year)
        return CCT

    def _cost_flue_gas_fan(self,Pi):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: EUR
        #Pi unit: W
        year=2018
        CF_Dollars=103193*math.pow(Pi/445/1000,0.67)*self._Eur2Dollars[year]/1e6  #unit M$ 
        CF=self._convert_to_RMB_in_target_year(CF_Dollars,year)
        return CF

    def _cost_HEN(self,At,p=1e5):
        #reference:https://doi.org/10.1016/j.jclepro.2019.02.049
        #the referred year: 2018
        #the given Money unit: EUR
        #Ai unit: m2,p: Pa
        year=2018
        p=p/1e5
        CHEN_Dollars=2546.9*math.pow(At,0.67)*math.pow(p,0.28)*self._Eur2Dollars[year]/1e6  #unit M$ 
        CHEN=self._convert_to_RMB_in_target_year(CHEN_Dollars,year)
        return CHEN

    
    def _cost_pump(self,W):
        #reference: Chemical process: design and integration
        #the referred year: 2001
        #the given Money unit: USD
        year=2001
        CPump_USD=9840*math.pow(W/1000/4,0.55)/1e6 #unit M$ 
        CPump=self._convert_to_RMB_in_target_year(CPump_USD,year)
        return CPump

    
    def _convert_to_RMB_in_target_year(self,cost_dollars,given_year):
        cost_RMB=cost_dollars/self._CEPCIs[given_year]*self._CEPCIs[self._target_year]*self._Dollars2RMB[self._target_year]
        return cost_RMB 

    






