
import math
import numpy as np

## construct cost indictor
construct_labour_cost_indictor=0.5
engineering_project_cost_indictor=0.35
TASC_multiplier=1.13
piping_integration_cost_indictor=0.05
M_cao = 56e-3  # kg/mol
M_caoh2 = 74e-3  # kg/mol

class Cost_Estimator(object):
    def __init__(self) -> None:
        self._Eur2Dollars={2014:39.9116/30.3287,2018:35.3758/30.1162}
        self._Dollars2RMB={2020:6.901}
        self._CEPCIs={2001:394.3, 2014:576.1,2018:603.1,2020:596.2,2022:816}
        #https://toweringskills.com/financial-analysis/cost-indices/
        self._target_year=2022
        
    def solve(self,design,energy_analysis_results,economic_inputs):
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
        make_up_limestone_percentage=design["dehy"]["make-up_limestone_percentage"]
        operation_costs["make-up_limestone"]=design["dehy"]["mole_camix_in"]*make_up_limestone_percentage \
            *M_caoh2*economic_inputs["operation_hours"]*3.6*economic_inputs["limestone_price"]/1e6
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
        CaOH2_storage_mass_circulated=design["dehy"]["CaOH2_storage_mass_circulated"]
        CaOH2_storage_mass_make_up=design["dehy"]["CaOH2_storage_mass_make_up"]
        total_CaOH2_storage_mass_ton=(CaOH2_storage_mass_circulated+CaOH2_storage_mass_make_up)/1000 #循环量和补充量
        initial_material_cost=economic_inputs["limestone_price"]*total_CaOH2_storage_mass_ton/1e6
        return initial_material_cost
    
    def calculate_equipment_costs(self,design,economic_inputs):
        equipment_costs={}
        equipment_costs.update(self.calculate_bretonHP_costs(design["BtHP"],economic_inputs))
        equipment_costs.update(self.calculate_dehydrator_costs(design["dehy"]))
        equipment_costs.update(self.calculate_hydrator_costs(design["hydr"]))
        equipment_costs.update(self.calculate_breton_costs(design["bret"]))
        equipment_costs["total"]=np.sum(list(equipment_costs.values()))
        return equipment_costs
    
    def calculate_bretonHP_costs(self,bretonHP_design):
        mass_co2=bretonHP_design["mass_co2"]
        T6=bretonHP_design["T6"]
        T7=bretonHP_design["T7"]
        T8=bretonHP_design["T8"]
        T9=bretonHP_design["T9"]
        P1=bretonHP_design["P1"]
        P2=bretonHP_design["P2"]
        P3=bretonHP_design["P3"]
        P4=bretonHP_design["P4"]
        invCosts={}
        
        invCosts["bHP_turb1"]=self._cost_bretonHP_turb(mass_co2,T6,T7,mius=0.88)
        invCosts["bHP_turb2"]=self._cost_bretonHP_turb(mass_co2,T8,T9,mius=0.88)
        invCosts["bHP_comp1"]=self._cost_bretonHP_comp(mass_co2,P1,P2,mius=0.88)
        invCosts["bHP_comp2"]=self._cost_bretonHP_comp(mass_co2,P3,P4,mius=0.88)
        
        
        # Heat tranfer network
        AHENbretonHP=bretonHP_design["total_HEN_area"]
        pHENbretonHP=bretonHP_design["p_bretonHP"]
        invCosts["CHENbretonHP"]=self._cost_HEN(AHENbretonHP,pHENbretonHP)
        

        # TODO: solid conveying system
        return invCosts
    def calculate_dehydrator_costs(self,dehy_design):
        invCosts={}
        Qdehy=-dehy_design["Q_dehydrator"]
        invCosts["Cdehy"]=self._cost_fluidized_bed_dehydrator(Qdehy)
        # Heat tranfer network
        AHENdehy=dehy_design["total_HEN_area"]
        pHENdehy=dehy_design["p_dehy"]
        invCosts["CHENdehy"]=self._cost_HEN(AHENdehy,pHENdehy)
        #water pump
        WPump=dehy_design["water_pump_power"]
        invCosts["CWaterPump"]=self._cost_pump(WPump)

        # TODO: solid conveying system
        return invCosts



    def _cost_fluidized_bed_dehydrator(self,Qc):
        #reference:https://doi.org/10.1016/j.ijggc.2019.01.005
        #the referred year: 2014
        #the given Money unit: EUR
        #Qc unit: W
        year=2014
        FC_Dollars=(0.217*math.pow(Qc/1e6,0.65)+3.83)*self._Eur2Dollars[year]  #unit M$
        FC=self._convert_to_RMB_in_target_year(FC_Dollars,year)
        return FC
    
    def _cost_bretonHP_turb(self,mass_co2,Tin,Tout,mius):
        year=2019
        TC_Dollars=mass_co2*492.2(1-mius)*((Tin+273.15)/(Tout+273.15))*math.log((Tin+273.15)/(Tout+273.15))(1+math.exp(0.036*(Tin+273.15)-65.66))
        TC=self._convert_to_RMB_in_target_year(TC_Dollars,year)
        return TC 
    
    def _cost_bretonHP_comp(self,mass_co2,Pin,Pout,mius):
        year=2019
        CC_Dollars=mass_co2*59.1(1-mius)*(Pout/Pin)*math.log(Pout/Pin)
        CC=self._convert_to_RMB_in_target_year(CC_Dollars,year)
        return CC 
    
    def _cost_pump(self,W):
        #reference: Chemical process: design and integration
        #the referred year: 2001
        #the given Money unit: USD
        year=2001
        CPump_USD=9840*math.pow(W/1000/4,0.55)/1e6 #unit M$ 
        CPump=self._convert_to_RMB_in_target_year(CPump_USD,year)
        return CPump
    

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

    def _convert_to_RMB_in_target_year(self,cost_dollars,given_year):
        cost_RMB=cost_dollars/self._CEPCIs[given_year]*self._CEPCIs[self._target_year]*self._Dollars2RMB[self._target_year]
        return cost_RMB 

