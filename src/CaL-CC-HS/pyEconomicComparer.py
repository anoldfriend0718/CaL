import math
import json
import numpy as np

class Wulumuqi_project(object):
    def __init__(self,economic_inputs) -> None:
        self.investment_cost=120*596.2/541.7 #百万人民币，2020
        self.operation_labour_cost_indictor=economic_inputs["operation_labour_cost_indictor"] 
        # self.maintain_cost_indictor=economic_inputs["maintain_cost_indictor"]
        self.maintain_cost_indictor=0.025 ## hardcode, not different from the CaL project
        self._electric_price=economic_inputs["elec_price"] #RMB/kWh
        # self._operation_hours=economic_inputs["operation_hours"]
        self._operation_hours=3435
        self.electricity=5489*1e4 #kWh
        self.electric_cost=self.electricity*self._electric_price/1e6 
        self.operation_cost=(self.operation_labour_cost_indictor+self.maintain_cost_indictor)*self.investment_cost\
            +self.electric_cost
        self.energy_efficiency=0.97

class economic_analyzer(object):
    def __init__(self) -> None:
        pass 

    #investment_cost 百万元
    #operation_cost 百万元
    #annular_heat J
    def LCOH(self,investment_cost,operation_cost,discount_rate,operation_year,annular_heat):
        total_cost=investment_cost*1e6
        total_heat=0
        for n in np.arange(1,operation_year+1):
            total_cost=total_cost+operation_cost*1e6/math.pow(1+discount_rate,n)
            total_heat=total_heat+annular_heat/1e9/math.pow(1+discount_rate,n)
        lcoh=total_cost/total_heat #yuan/GJ
        return lcoh 

class economic_comparer(object):
    def __init__(self,economic_inputs) -> None:
        self._discount_ratio=economic_inputs["discount_ratio"]
        self._annular_operation_hours=economic_inputs["operation_hours"]
        self._operation_years=economic_inputs["operational_years"]
        
        self._p2=Wulumuqi_project(economic_inputs)
        self._analyzer=economic_analyzer()

    def compare(self,plant_results):
        results={}
        electricity1=plant_results["metrics"]["energy"]["total_power"] #W
        electricity2=self._p2.electricity*1e3/self._p2._operation_hours #W
        mass_rate_CO2=plant_results["plant"]["carb"]["m_co2_capture"]
        specific_capture_energy=(electricity1-electricity2)/mass_rate_CO2 #J/kg
        results["specific_capture_energy"]=specific_capture_energy/1e6 #GJ/t

        investment_cost1=plant_results["metrics"]["economic"]["construction"]["total as-spent"]
        operation_cost1=plant_results["metrics"]["economic"]["operation"]["total as-spent"]
        annular_heat1=plant_results["metrics"]["energy"]["Q_hot_water"]*self._annular_operation_hours*3600 #J
        LCOH1=self._analyzer.LCOH(investment_cost1,operation_cost1,self._discount_ratio,
                self._operation_years,annular_heat1)
        results["LCOH1"]=LCOH1 #RMB/GJ

        investment_cost2=self._p2.investment_cost
        operation_cost2=self._p2.operation_cost
        annular_heat2=plant_results["metrics"]["energy"]["Q_hot_water"]*self._p2._operation_hours*3600 #J
        LCOH2=self._analyzer.LCOH(investment_cost2,operation_cost2,self._discount_ratio,
                self._operation_years,annular_heat2)
        results["LCOH2"]=LCOH2 #RMB/GJ

        mass_rate_CO2_per_hot_GJ=mass_rate_CO2/(plant_results["metrics"]["energy"]["Q_hot_water"]/1e9) #kg/GJ
        specific_capture_cost=(LCOH1-LCOH2)/mass_rate_CO2_per_hot_GJ*1000 #RMB/t
        results["specific_capture_cost"]=specific_capture_cost
        return results


if __name__=="__main__":
    with open("/home/anoldfriend/Workspace/MyRepo/thermodynamics/CaL/src/CaL-CC-HS/tmp_results4.json","r") as fp:
        plant_results=json.load(fp)
    economic_inputs={}
    economic_inputs["elec_price"]=0.165 #元/千瓦时
    economic_inputs["operation_hours"]=3435 # hours
    economic_inputs["discount_ratio"]=6/100 #8%
    economic_inputs["operational_years"]=30
    economic_inputs["operation_labour_cost_indictor"]=0.025
    economic_inputs["maintain_cost_indictor"]=0.025

    comparer=economic_comparer(economic_inputs)
    results=comparer.compare(plant_results)
    print(results)
        


        
