# This file contains parameters and constants used in the simulation.
class constants:
  def __init__(self):
    self.au = 149597870.7 # km (Astronomical Almanac 2024)
    self.G = 6.67428e-20 # km^3/kg/s^2 (Astronomical Almanac 2024)
    self.SRP0 = 4.56e-6 # N/m^2 Solar radiation pressure at 1 au
    self.day = 86400 # seconds

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