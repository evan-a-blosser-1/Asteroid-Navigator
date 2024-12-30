"""
    WORK NOTES:  Double check all planetrary values!!!!
                  - cite sources for all values
    
    
    
    
"""


# This file contains parameters and constants used in the simulation.
class constants:
  def __init__(self):
    self.au = 149597870.7  # km (Astronomical Almanac 2024)
    self.G = 6.67428e-20   # km^3/kg/s^2 (Astronomical Almanac 2024)
    self.SRP0 = 4.56e-6    # N/m^2 Solar radiation pressure at 1 au
    self.day = 86400       # seconds
    self.mu_sun = 6.67428e-20*1.9884e30 # mass = 1.9884e30 kg (NASA NSSDC/GSFC 2024)


const = constants()
G = const.G # km^3/kg/s^2

# Major bodies physical parameters
class sun:
  def __init__(self):
    self.mass = 1.9884e30 # kg (NASA NSSDC/GSFC 2024)
    self.GM = G*self.mass
    self.Re = 695700 # km (volumetric mean radius - JPL Horizons rev. 2013)
    self.spin = 2.86533e-6 # rad/s (NASA NSSDC/GSFC 2024)
    self.eps = 7.25 # deg (obliquity of the ecliptic - NASA NSSDC/GSFC 2024)
    
    
class mercury:
  def __init__(self):
    self.mass = 3.3011e23 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 2439.7 # km (equatorial radius - Astronomical Almanac 2024)
        
class venus:
  def __init__(self):
    self.mass = 4.8675e24 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 6051.8 # km (equatorial radius - Astronomical Almanac 2024)
    
############################################################################
class earth:
    def __init__(self):
      self.mass = 5.9722e24 # kg
      self.GM = G*self.mass # km^3/s^2
      self.Re = 6378.1366 # km (equatorial radius - Astronomical Almanac 2024)
      self.spin = 7.292115e-5 # rad/s (rotational speed - Astronomical Almanac 2024)
      self.eps = 23.4392911 # deg (obliquity of the ecliptic - Astronomical Almanac 2024)
      self.J2 = 1.0826359e-3 # (Astronomical Almanac 2024)

class moon:
  def __init__(self):
    self.mass = 7.345828157e22 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass # km^3/s^2
    self.Re = 1737.4 # km (equatorial radius - Astronomical Almanac 2024)
    self.spin = 2.6617e-6 # rad/s (rotational speed - JPL Horizons rev. 2013)
    self.eps = 6.67 # deg (obliquity of the ecliptic - JPL Horizons rev. 2013)
#############################################################################
    
class mars:
  def __init__(self):
    self.mass = 6.4171e23 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 3389.5 # km (equatorial radius - Astronomical Almanac 2024)
    

class jupiter:
  def __init__(self):
    self.mass = 1.8982e27 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 69911 # km (equatorial radius - Astronomical Almanac 2024)
    
class saturn: 
  def __init__(self):
    self.mass = 5.6834e26 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 58232 # km (equatorial radius - Astronomical Almanac 2024)
    
class uranus:
  def __init__(self):
    self.mass = 8.6810e25 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 25362 # km (equatorial radius - Astronomical Almanac 2024)
    
class neptune:
  def __init__(self):
    self.mass = 1.02413e26 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 24622 # km (equatorial radius - Astronomical Almanac 2024)  
    
class pluto:
  def __init__(self):
    self.mass = 1.303e22 # kg (Astronomical Almanac 2024)
    self.GM = G*self.mass
    self.Re = 1188.3 # km (equatorial radius - Astronomical Almanac 2024) 
    





# # Asteroid physical parameters
# class Apophis:
#   def __init__(self):
#     self.mass = 5.31e10               # kg (Diogo 2021)
#     self.GM = G*self.mass
#     self.spinRate = 30.4              # rot/hr
#     self.Re = 0.1935                  # km (volumetric mean radius - Diogo 2021)
#     self.spin =  30.4                 # (30.4h per rotation Diogo 2021)
#     self.Poly_x = -0.002150           # volInt.c output
#     self.Poly_y = -0.001070           # Polyhedron Center of Mass (CM)
#     self.Poly_z = -0.000308           # -x,y,z (km)
#     self.gamma  = 0.2848196900000026  # Asteroid Scale to Km
# class DA1950_Prograde:
#     def __init__(self):
#         self.mass = 1.0e10                # kg
#         self.GM = G * self.mass
#         self.Re = 0.58                    # km (volumetric mean radius)
#         self.spin = 2.1216                # (2.1216 rot/hr)
#         self.Poly_x = 0.016325            # volInt.c output
#         self.Poly_y = 0.026771            # Polyhedron Center of Mass (CM)
#         self.Poly_z = 0.016897            # -x,y,z (km)
#         self.gamma = 1.0004426659700023   # Asteroid Scale to Km

# # Dictionary to map class names to class objects
# class_map = {
#     'Apophis': Apophis,
#     'DA1950_Prograde': DA1950_Prograde,
# }



# def create_new_model(class_name, mass, GM, Re, spin, Poly_x, Poly_y, Poly_z, gamma):
#     # Define a new class dynamically
#     new_class = type(class_name, (object,), {
#         '__init__': lambda self: None,
#         'mass': mass,
#         'GM': GM,
#         'Re': Re,
#         'spin': spin,
#         'Poly_x': Poly_x,
#         'Poly_y': Poly_y,
#         'Poly_z': Poly_z,
#         'gamma': gamma
#     })
    
#     # Add the new class to the class_map dictionary
#     class_map[class_name] = new_class

# # Example usage:
# # create_new_model('NewModel', 1.0e10, G * 1.0e10, 0.58, 2.1216, 0.016325, 0.026771, 0.016897, 1.0004426659700023)