import math
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt

#Functions related to time
  
def julian(dt):
    """Finds the julian date based on the SAM Photovoltaic Model Techinical Reference Update, Paul Gilman, Aron Dobos, Nicholas DiOrio, Janine Freeman, Steven Janzou, and David Ryberg, of National Renewable Energy Laboratory. March 2018. After Michalsky 1988."""
    
    days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour + dt.minute/60
    timezone = dt.utcoffset()
    
    k=0
    if (year%4 == 0):
        k = 1
    
    jdoy = day + sum(days[0:(month-1)])
    if month > 2:
        jdoy = jdoy + k
    
    try:
        tutc = hour - timezone
        if tutc<0:
            tutc = tutc+24
        if tutc>24:
            tutc = tutc-24
            
    except:
        tutc = hour
        
    return (32916.5 + 365*(year-1949)+(year-1949)//4+jdoy+tutc/24-51545)

def Greenwich_mean_siderial_time(julian, dt):
    """Calculates the siderial time from the julian and datetime (requires timezone (pytz)) (see NREL, Michalsky, Walraven)."""
    a = 6.697375+0.0657098242*julian+(dt.hour + dt.minute/60-dt.utcoffset().seconds/3600)
    a = a - 24*(a//24)
    if a < 0:
        a = a + 24
    return(a)

#Functions relating to radiation

def mean_extraterrestrial_radiation(datetime):
    """Extraterrestrial radiation incident on a plane normal to the radiation. Eq. 1.4.1b Duffie and Beckman."""
    
    doy = datetime.dt.dayofyear+datetime.dt.hour/24+datetime.dt.minute/1440-datetime.dt.utcoffset().seconds/86400

    if doy>365:
        doy = doy-365
    elif doy<0:
        doy = doy+365
            
    B = (doy-1)*2*math.pi/365
    return(1367*(1.000110 + 0.034221 * math.cos(B) + 0.001280 * math.sin(B)+ 0.000719 * math.cos(2*B)+0.000077*math.sin(2*B)))


def incident_radiation(radiation, angle_of_incidence):
    """Returns the radiation incident to a collector (this is actually just taking the incident component of the radiation and therefore will work for extraterrestrial or measured radiation)."""
    if (0<angle_of_incidence<math.pi/2):
        return radiation*math.cos(angle_of_incidence)
    elif (angle_of_incidence == 0):
        return radiation
    else:
        return 0


    
    
#Functions relating to the angle of the sun

def refraction_corrected_elevation(elevation_angle):
    """Returns the refraction corrected elevation angle in radians. From Eq. 4.17 NREL"""
    
    alpha_0d = math.degrees(elevation_angle)
    if alpha_0d>-0.56:
        r = 3.51561*(0.1594+0.0196*alpha_0d +0.00002*alpha_0d**2)/(1+0.505*alpha_0d +0.0845*alpha_0d**2)
    elif alpha_0d<=-0.56:
        r = 0.56
    if (alpha_0d+r)>90:
        return(math.pi/2)
    elif (alpha_0d+r<=90):
        return(math.radians(alpha_0d+r))   




class EclipticCoordinate(object):
    """Ecliptic coordinate of the sun apparent to the earth.
        Mean Longitude - mnlong (deg)
        Mean Anomaly - mnanom (rad)
        Ecliptic Longitude (rad)
        Obliquity (rad)
        Right Ascension (rad)
        Declination Angle (rad)
        """
    
    def __init__(self, dt):
        self.dt = dt
        self.julian = julian(dt)
        self.mnlong = self.set_mean_longitude()
        self.mnanom = self.set_mean_anomaly()
        self.eclong = self.set_ecliptic_longitude()
        self.obleq = self.set_oblequity()
        self.ra = self.set_right_ascension()
        self.decrad = self.set_declination()
        
    def __repr__(self):
        return "dt = {7},julian = {0}, mnlong = {1}, mnanom = {2}, eclong = {3}, obliq = {4}, ra = {5}, decrad = {6}".format(self.julian, self.mnlong, self.mnanom, self.eclong, self.obleq, self.ra, self.decrad, self.dt)
    
    def set_mean_longitude (self):
        """Returns the mean solar longitude in degrees (see NREL, Michalsky, Walraven)."""
        mnlong = 280.46+0.9856474*self.julian
        mnlong = mnlong - 360*(mnlong//360)
        if mnlong < 0:
            mnlong = mnlong + 360
        return(mnlong)

    def set_mean_anomaly (self):
        """Returns the mean solar anomaly in radians (see NREL, Michalsky, Walraven)"""
        mnanom = (357.528 + 0.9856003*self.julian)
        mnanom = mnanom - 360*(mnanom//360)
        if mnanom < 0:
            mnanom = mnanom + 360
            
        return(math.radians(mnanom))

    def set_ecliptic_longitude(self):
        """Returns the ecliptic longitude in radians (see NREL, Michalsky, Walraven)"""
        eclong = (self.mnlong+1.915*math.sin(self.mnanom)+0.02*math.sin(2*self.mnanom))
        eclong = eclong-360*(eclong//360)
        
        if eclong <0:
            eclong = eclong + 360
        eclong = math.radians(eclong)

        return eclong

    def set_oblequity(self):
        """Returns the obliquity angle in radians(see NREL, Michalsky, Walraven)."""

        obleq = 23.439-0.0000004*self.julian

        return math.radians(obleq)

    def set_right_ascension(self):
        """Returns the right ascension"""

        if math.cos(self.eclong)<0:
            ra = math.atan(math.cos(self.obleq)*math.sin(self.eclong)/math.cos(self.eclong))+math.pi
        elif math.cos(self.obleq)*math.sin(self.eclong)<0:
            ra = math.atan(math.cos(self.obleq)*math.sin(self.eclong)/math.cos(self.eclong))+2*math.pi
        else:
            ra = math.atan(math.cos(self.obleq)*math.sin(self.eclong)/math.cos(self.eclong))
        return ra

    
    def set_declination(self):
        """Returns declination angle in radians."""
        return(math.asin(math.sin(self.obleq)*math.sin(self.eclong)))


    def equation_of_time(self):
        """Returns equation of time"""
        a = 1/15*(self.mnlong-180/math.pi*self.ra)
        
        if (-0.33<=a<=0.33):
            return a
        
        elif a<-0.33:
            return(a+24)
        
        elif a>0.33:
            return (a-24)
        
        else:
            print("Issue calculating equation of time.")
            return None
    
    

 
    
class Position(object):
    """Coordinate as decimal latitude and longitude."""
    def __init__(self, latitude, longitude):
        self.lat = latitude
        self.lon = longitude
        
    def __repr__(self):
        return "{0},{1}".format(self.lat, self.lon)
    
    #Time related methods
    
    def local_mean_siderial_time(self, gmst):
        """Calculates the siderial time at the local position from Greenwich Mean Siderial Time, and the longitude (see NREL, Michalsky, Walraven)."""
        a = gmst + self.lon/15
        a = a - 24*(a//24)
        if a < 0:
            a = a + 24
        return(a)

    def hour_angle(self, lmst, right_ascension):
        """Determines hour angle based on the local mean siderial time and the right ascension angle."""

        b = 15*math.pi/180*lmst-right_ascension

        if (b < -math.pi):
            b = b+ 2*math.pi
            
        elif b > math.pi:
            b = b - (2*math.pi)
        
        return(b)
    
    #Methods for solar position relative to the point lat/long.
    
    def elevation_angle(self, declination, hour_angle):
        """Determines the elevation angle based on inputs. Radians"""
        x = math.sin(declination)*math.sin(math.radians(self.lat)) + math.cos(declination)*math.cos(math.radians(self.lat))*math.cos(hour_angle)

        if (x>=-1) and (x<=1):
            return math.asin(x)
        else:
            return np.sign(x)*math.pi/2
        
        
    def zenith_angle(self, declination, hour_angle):
        zenith = math.pi/2-self.elevation_angle(declination, hour_angle)
        return zenith
        

    def azimuth_angle(self, elevation_angle, declination, hour_angle):
        """Returns the azimuth angle in radians based on
        Elevation Angle (radians) = elrad
        Latitude (degrees) = lat
        Declination Angle (radians) = decrad"""

        a = (math.sin(elevation_angle)*math.sin(math.radians(self.lat))-math.sin(declination))/(math.cos(elevation_angle)*math.cos(math.radians(self.lat)))
        try:
            b = math.acos(a)
        except:
            if a > 1:
                b = 0
            elif (math.cos(elevation_angle)==0) or a<-1:
                b = math.pi

        if (hour_angle < -math.pi):
            return (b)
        elif ((-math.pi<=hour_angle) and (hour_angle<=0)) or (hour_angle>=math.pi):
            return(math.pi-b)
        elif (0<hour_angle) and (hour_angle<math.pi):
            return (math.pi+b)
        
    #Methods for sunrise/sunset    
    
    def sunrise_hour_angle(self, declination):
        """Returns the hour angle at sunrise/sunset. If there is no sunrise/sunset. Returns 0 if the sun is down all day, pi if the sun is up all day."""
        a = -math.tan(math.radians(self.lat))*math.tan(declination)
        if a>=1:
            return(0)
        
        elif a<=-1:
            return math.pi
        
        elif (-1<a<1):
            return(math.acos(a))
    
    def sunrise_hour(self, declination, timezone, EOT):
        """Determines the hour at sunrise on a given day of the year (declination)."""
        a = 12-1/15*180/math.pi*self.sunrise_hour_angle(declination)-(self.lon/15-timezone)-EOT
        return (a)
    
    def sunrise_hour(self, declination, timezone, EOT):
        """Determines the hour at sunset on a given day of the year (declination)."""
        a = 12+1/15*180/math.pi*self.sunrise_hour_angle(declination)-(self.lon/15-timezone)-EOT
        return (a)
    
    
class SolarConfig(Position):
    """SolarConfig is a subclass of Position. It contains orientation attributes in addition to latitude and longitude. Methods will determine the angle of incidence etc.
    Azimuth Angle between north and the base of the panel with 0 facing West."""
    
    def __init__(self, latitude, longitude, slope, surface_azimuth):
        Position.__init__(self,latitude, longitude)
        self.slop = slope
        self.surface_azimuth = surface_azimuth
    
    def angle_of_incidence(self, zenith_angle, azimuth_angle):
        a = math.sin(zenith_angle)*math.cos(math.radians(self.surface_azimuth)- azimuth_angle)*math.sin(radians(self.slope))+math.cos(zenith_angle)*math.cos(math.radians(self.slope))
        
        if a<-1:
            return (math.pi)
        elif a>1:
            return 0
        elif (-1<=a<=1):
            return math.acos(a)
        
    def HDKR_diffuse_radiation(self, angle_of_incidence, zenith_angle, extraterrestiral_irradiance, Eb, Ed, albedo = 0):
        """Calculates the diffuse radiation based on the Hay, Davies, Klucher and Reindl model. Including the reflectance which is initialised at zero (albedo). Parameters use are named after those in NREL.
        Ibh : Incident beam irradiadiance on "Plane of Array".
        Igh : Total irradiance on "Plane of Array".
        Eb : Beam horizontal irradiance (measured).
        Ed : Diffuse horizontal irradiance (measured).
        Rb : Ratio of beam incident to "Plane of Array" to horizontal beam.
        Ai : The anisotropy index for forward scattering circumsolar diffuse irradiance.
        f : Modulating factor for horizontal brightening correction.
        s : Horison brightening correction factor.
        
        Note that the isohor term contains both isotropic and horizonal brightening iso*(1+x)."""
        
        Ibh = Eb*math.cos(zenith_angle)
        
        Igh = Ibh + Ed
        
        Rb = math.cos(angle_of_incidence)/math.cos(zenith_angle)
        
        Ai = Ibh/extraterrestrial_irradiance
        
        f = math.sqrt(Ibh/Igh)
        
        s = (math.sin(self.slope/2))**3
        
        cir = Ed*Ai*Rb
        iso = Ed*(1-Ai)*(1+math.cos(self.slope))/2
        isohor = iso*(1+f*s)
        
        #Ground reflected diffuse (albedo is ground reflectance, if not provided = 0)
        Ir = albedo*(Igh)*(1-math.cos(self.slope))/2
        
        return (isohor+cir+Ir)
    
    
    
     
    
    


    

    

        
