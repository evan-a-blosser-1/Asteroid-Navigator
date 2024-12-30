"""

"""
import sys
from astroquery.jplhorizons import Horizons
import numpy as np
################## Personal Packages
sys.dont_write_bytecode = True

############################################
#%% Ephemris Lookup Funciton
def ephemeris(body, start_date, end_date, step): 

  obj = Horizons(id=body,
                 location='500@10',
                 epochs={'start':start_date, 'stop':end_date, 'step':step})
  
  
  oe = obj.elements(refplane='ecliptic')
  
  return oe


#%%
###############################
sys.dont_write_bytecode = True
import constants as C
const   = C.constants()
sun     = C.sun()
mercury = C.mercury()
venus   = C.venus()
earth   = C.earth()
moon    = C.moon()
mars    = C.mars()
jupiter = C.jupiter()
saturn  = C.saturn()
uranus  = C.uranus()
neptune = C.neptune()
pluto   = C.pluto()
##########################################################
#####################  Horizens Target Selection Function
def target_selection(Target='Earth'):
      if Target == 'Mercury':
        ID = '1'
        mu_p = mercury.mass
        r_p  = mercury.Re
      elif Target == 'Venus':
        ID = '2'
        mu_p = venus.mass
        r_p  = venus.Re
      elif Target == 'Earth':
        ID = '3'
        mu_p = earth.mass
        r_p  = earth.Re
      elif Target == 'Mars':
        ID = '4'
        mu_p = mars.mass
        r_p  = mars.Re
      elif Target =='Jupiter':
        ID = '5'
        mu_p = jupiter.mass
        r_p  = jupiter.Re
      elif Target == 'Saturn':
        ID = '6'
        mu_p = saturn.mass
        r_p  = saturn.Re
      elif Target == 'Uranus':
        ID = '7'
        mu_p = uranus.mass
        r_p  = uranus.Re
      elif Target == 'Neptune':
        ID = '8'
        mu_p = neptune.mass
        r_p  = neptune.Re
      elif Target == 'Pluto':
        ID = '9'
        mu_p = pluto.mass
        r_p  = pluto.Re
      else:
        print('| ERROR!: {Target} is not a valid target. Please select a valid target from the list below: |')
        print('| Mercury | Venus | Earth | Mars | Jupiter | Saturn | Uranus | Neptune | Pluto |')
        
      return ID, mu_p, r_p
  
  
#%% Departure Delta-V Calculation
def DelV_Hom(Location='Earth',
             Destination='Saturn',
             R1_park_alt=600,
             R2_park_alt=600,
             Start_Date='2021-01-01',
             End_Date='2021-01-02',
             Step='1d'):
    #########################################
    # G*M of the Sun
    mu_sun = 6.67428e-20*1.9884e30
    ##########################################
    ##### Astroquery, both targets
    ID_dep, mu_dep, r_dep = target_selection(Location)
    ID_arr, mu_arr, r_arr = target_selection(Destination)
    oe_dep = ephemeris(ID_dep, Start_Date, End_Date, Step)
    oe_arr = ephemeris('606', Start_Date, End_Date, Step)    
    # Semi-Major Axis of depature planet orbit
    R1 = oe_dep['a'][0]
    # Semi-Major Axis of arrival planet orbit
    R2 = oe_arr['a'][0] 
    ##########################################################
    ############################################## Departure
    # Heliocentric velocity of craft on departure from earth           
    v_helio_craft = np.sqrt(mu_sun*((2/R1)-(2/(R1+R2))))             
    # Departure planet's Heliocentric Velocity                                    
    v_helio_planet = np.sqrt(mu_sun/R1)        
    # The velocity of the departure hyperbola                          
    v_inf_dep =  v_helio_craft - v_helio_planet                 
    # Geocentric Spacecraft velocity                                   
    v_geo = np.sqrt(mu_dep/(R1_park_alt + r_dep))                          
    # Geocentric Space craft Velocity at perigee of Dep. Hyp.          
    v_geo_dep = np.sqrt(v_inf_dep**2 + (2*mu_dep)/(R1_park_alt+r_dep)) 
    # Departure Delta-V                                                
    Delta_V_dep = v_geo_dep - v_geo    
    #######################################################################
    ################# Arrival Delta-V #####################################
    #######################################################################
    # Heliocentric velocity of craft on departure from earth              
    v_helio_craft_arrive = np.sqrt(mu_sun*((2/R2)-(2/(R1+R2))))           
    # Earth's Heliocentric Velocity                                       
    v_helio_arrive_planet = np.sqrt(mu_sun/R2)                            
    # The velocity of the departure hyperbola                             
    v_inf_ari = v_helio_craft_arrive - v_helio_arrive_planet              
    # Geocentric Spacecraft velocity                                      
    v_geo_Rende = np.sqrt(mu_arr/(R2_park_alt+r_arr))                     
    # Geocentric Space craft Velocity at perigee of Dep. Hyp.             
    v_geo_ari = np.sqrt(v_inf_ari**2 + (2*mu_arr)/(R2_park_alt+r_arr))    
    # Departure Delta-V                                                   
    Delta_V_ari = v_geo_ari - v_geo_Rende                                 
    #####################################################################
    ################# Time of Flight ####################################
    # Time of Flight from Earth to Saturn for Hohman           
    TOF_Sec = (np.pi)/(np.sqrt(mu_sun))*((R1 + R2)/2)**(3/2)   
    # Time of Flight in years                                  
    TOF_days = TOF_Sec/(60*60*24)                              
    TOF_years = TOF_Sec/(60*60*24*365.256)                     
    ####################
    # TOTAL DELTA-V    #
    # Before Maneuvers #
    ######################
    # Total Delta-V Req. ######################################
    Delta_V_To_Destination = abs(Delta_V_dep) + abs(Delta_V_ari) 
    ###########################################################
    Mission_Time_Data_Message = f"""
{'-'*42}
| Mission TOF: {TOF_days} Days
|            : {TOF_years} Years
{'-'*42}
| Mission Start Time: {Mission_Start_Time_Horizons} Days
{'-'*42}
| Mission End Time: {Mission_Time} Days
{'-'*42}
"""                                    
    print(Mission_Time_Data_Message)  
    
          
    return  Delta_V_To_Destination 
    
    
    

    

def Fuel_Consumption(Craft_mass=2000,I_sp=300,Delta_V_To=10):
    # Accel. of Gravity (km/s^2) 
    g_0 = 9.81e-3    
    # Propellent Req. for the departure ##################################
    Fuel_Mass_Ratio_0 = 1 - np.exp(-Delta_V_To/(I_sp*g_0))      
    # Spacecraft Initial mass Equation                                   
    Mass_Prop_0 = (Fuel_Mass_Ratio_0*Craft_mass)/(1 - Fuel_Mass_Ratio_0) 
    # Return fuel needed for transfer
    return Mass_Prop_0







#%%






