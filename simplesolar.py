import math
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt

#global variables
daysinmonth = {"jan" : 31, "feb" : 28, "mar" : 31, "apr" : 30, "may" : 31, "jun" : 30, "jul" : 31, "aug" : 31, "sep" : 30, "oct" : 31, "nov" : 30, "dec" : 31}

#Functions (typically not related to the position)

  
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

def gmst(julian, dt):
    a = 6.697375+0.0657098242*julian+dt.hour + dt.minute/60+dt.utcoffset().seconds/3600
    a = a - 24*(a//24)
    if a < 0:
        a = a + 24
    return(a)

def lmst(gmst, lon):
    a = gmst + lon/15
    a = a - 24*(a//24)
    if a < 0:
        a = a + 24
    return(a)

def elev_angle(lat, decrad, hrrad):
    """Determines the elevation angle based on inputs. Radians"""
    x = math.sin(decrad)*math.sin(math.radians(lat)) + math.cos(decrad)*math.cos(math.radians(lat))*math.cos(hrrad)
    
    if (x>=-1) and (x<=1):
        return math.asin(x)
    else:
        return np.sign(x)*math.pi/2

def refraction_corrected_elevation_angle(elrad):
    """takes elevation angle in radians and corrects for refraction."""
    alpha_0d = math.degrees(elrad)
    if alpha_0d>-0.56:
        r = 3.51561*(0.1594+0.0196*alpha_0d +0.00002*alpha_0d**2)/(1+0.505*alpha_0d +0.0845*alpha_0d**2)
    elif alpha_0d<=-0.56:
        r = 0.56
    if (alpha_0d+r)>90:
        return(math.pi/2)
    elif (alpha_0d+r<=90):
        return(math.radians(alpha_0d+r))

def azimuth_angle(elrad, lat, decrad, hrrad):
    """Returns the azimuth angle in radians based on
    Elevation Angle (radians) = elrad
    Latitude (degrees) = lat
    Declination Angle (radians) = decrad"""
    
    a = (math.sin(elrad)*math.sin(math.radians(lat))-math.sin(decrad))/(math.cos(elrad)*math.cos(math.radians(lat)))
    try:
        b = math.acos(a)
    except:
        if a > 1:
            b = 0
        elif (math.cos(elrad)==0) or a<-1:
            b = math.pi

    if (hrrad < -math.pi):
        return (b)
    elif ((-math.pi<=hrrad) and (hrrad<=0)) or (hrrad>=math.pi):
        return(math.pi-b)
    elif (0<hrrad) and hrrad<math.pi:
        return (math.pi+b)



class EclipticCoordinate(object):
    """Ecliptic coordinate of the sun apparent to the earth.
        Mean Longitude - mnlong (deg)
        Mean Anomaly - mnanom (rad)
        Ecliptic Longitude (rad)
        Obliquity (rad)
        Right Ascension (rad)
        Declination Angle (rad)"""
    
    def __init__(self, dt):
        self.dt = dt
        self.julian = julian(dt)
        self.mnlong = self.set_mnlong()
        self.mnanom = self.set_mnanom()
        self.eclong = self.set_eclong()
        self.obleq = self.set_obleq()
        self.ra = self.set_ra()
        self.decrad = self.set_decrad()
        
    def __repr__(self):
        return "dt = {7},julian = {0}, mnlong = {1}, mnanom = {2}, eclong = {3}, obliq = {4}, ra = {5}, decrad = {6}".format(self.julian, self.mnlong, self.mnanom, self.eclong, self.obleq, self.ra, self.decrad, self.dt)
    
    def set_mnlong (self):
        """Returns the mean solar longitude in degrees (see NREL, Michalsky, Walraven)"""
        mnlong = 280.46+0.9856474*self.julian
        mnlong = mnlong - 360*(mnlong//360)
        if mnlong < 0:
            mnlong = mnlong + 360
        return(mnlong)

    def set_mnanom (self):
        """Returns the mean solar anomaly in radians (see NREL, Michalsky, Walraven)"""
        mnanom = (357.528 + 0.9856003*self.julian)
        mnanom = mnanom - 360*(mnanom//360)
        if mnanom < 0:
            mnanom = mnanom + 360
            
        return(math.radians(mnanom))

    def set_eclong(self):
        """Returns the ecliptic longitude in radians (see NREL, Michalsky, Walraven)"""
        eclong = (self.mnlong+1.915*math.sin(self.mnanom)+0.02*math.sin(2*self.mnanom))
        eclong = eclong-360*(eclong//360)
        
        if eclong <0:
            eclong = eclong + 360
        eclong = math.radians(eclong)

        return eclong

    def set_obleq(self):
        """Returns the obliquity angle in radians(see NREL, Michalsky, Walraven)."""

        obleq = 23.439-0.0000004*self.julian

        return math.radians(obleq)

    def set_ra(self):
        """Returns the right ascension"""

        if math.cos(self.eclong)<0:
            ra = math.atan(math.cos(self.obleq)*math.sin(self.eclong)/math.cos(self.eclong))+math.pi
        elif math.cos(self.obleq)*math.sin(self.eclong)<0:
            ra = math.atan(math.cos(self.obleq)*math.sin(self.eclong)/math.cos(self.eclong))+2*math.pi
        else:
            ra = math.atan(math.cos(self.obleq)*math.sin(self.eclong)/math.cos(self.eclong))
        return ra
    
    def set_decrad(self):
        """Returns declination."""
        return(math.asin(math.sin(self.obleq)*math.sin(self.eclong)))

class Position(object):
    """Coordinate as decimal latitude and longitude."""
    def __init__(self, latitude, longitude):
        self.lat = latitude
        self.lon = longitude
        
    def __repr__(self):
        return "{0}_{1}".format(self.lat, self.lon)

class SolarPosition(object):
    """Has angluar infomation based on EclipticCoordinate and Position.
        Greenwich Mean Siderial Time = gmst
        Local Mean Siderial Time = lmst
        Hour Angle (rad) = hrrad
        Elevation Angle (rad) = elrad
        Refraction Corrected Elevation Angle (rad) = cerad
        Azimuth Angle (rad) = azrad
        Zenith Angle (rad) = znrad"""
    
    def __init__(self, Position, EclCoord):
        self.gmst = gmst(EclCoord.julian, EclCoord.dt)
        self.lmst = lmst(self.gmst, Position.lon) + Position.lon/15
        self.hrrad = self.set_hrrad(EclCoord.ra)
        self.elrad = elev_angle(Position.lat, EclCoord.decrad, self.hrrad)
        self.cerad = refraction_corrected_elevation_angle(self.elrad)
        self.azrad = azimuth_angle(self.elrad, Position.lat, EclCoord.decrad, self.hrrad)
        self.znrad = math.pi/2-self.cerad
        
    def __repr__(self):
        return "gmst = {0}, lmst = {1}, hrrad = {2}, elrad = {3}, cerad = {4}, azrad = {5}, znrad = {6}".format(self.gmst, self.lmst, self.hrrad, self.elrad, self.cerad, self.azrad, self.znrad)

    def set_hrrad(self, ra):
        """sets hour angle"""

        b = 15*math.pi/180*self.lmst-ra

        if (b < -math.pi):
            b = b+ 2*math.pi
            
        elif b > math.pi:
            b = b - (2*math.pi)
        
        return(b)


# def declination_angle(dt):
#     """Angular position of the sun at solar noon by more accurate calculation. Eq. 1.6.1b Duffie and Beckman
#     The angle between the rays of the Sun and the plane of the Earth's equator."""
    
#     n = (dt - datetime.datetime(dt.year//4*4,1,1)).total_seconds()
#     B = (n-1)*2*math.pi/31556926
#     b = 0.006918 - 0.399912*math.cos(B) + 0.070257*math.sin(B) - 0.006758*math.cos(2*B)+0.000907*math.sin(2*B)-0.002697*math.cos(3*B)+0.00148*math.sin(3*B)
    
#     return (math.degrees(b))

# """Classes are:
#         Positon - has lat long and all methods relating only to those attributes
#         Orientation - now redundant has only attributes relating to the slope and surface azimuth
#         SolarConfig - has position and orientation attributes and methods that require both (angle of incidence)""" 


    
    
# class Position(object):
#     """Coordinate as decimal latitude and longitude."""
#     def __init__(self, latitude, longitude):
#         self.lat = latitude
#         self.lon = longitude
        
#     def __repr__(self):
#         return "{0}_{1}".format(self.latitude, self.longitude)
    
#     #Methods relating to time:
    
#     def solar_time(self, UTS_datetime):
#         """Returns solar time for input clock time. Eq. 1.5.2 Duffie and Beckman."""
        
#         n = UTS_datetime.dayofyear+UTS_datetime.hour/24+UTS_datetime.minute/1440

#         B = (n-1)*2*math.pi/365
#         E = 229.2*(0.000075 + 0.001868*math.cos(B) - 0.032077*math.sin(B) - 0.014615*math.cos(2*B) - 0.04089*math.sin (2*B))

#         return(UTS_datetime + pd.Timedelta(minutes = 4*(self.longitude)+E))
    
#     def hour_angle(self, UTS_datetime = None, solarTime = None):
#         """Takes a pd.datetime value for solar time and returns the solar hour angle to the second."""
#         try:
#             decimal_hour_time = solarTime.hour + solarTime.minute/60 +solarTime.second/3600
            
#         except:
#             solar_time = self.solar_time(UTS_datetime)
#             decimal_hour_time = solar_time.hour+solar_time.minute/60+solar_time.second/3600
            
#         return(360*(decimal_hour_time-12)/24)

#     def sun_rise(self, UTS_datetime = None, declinationAngle = None):
#         """Returns the hour angle at sunrise/sunset."""
#         try:
#             delta = math.radians(declinationAngle)
#             phi = math.radians(self.latitude)
            
        
#         except:
#             delta = math.radians(declination_angle(UTS_datetime))
#             phi = math.radians(self.latitude)

#         return (math.degrees(math.acos(-math.tan(phi)*math.tan(delta))))                                     
                                     
                                     
                                     
#     #Methods for Position:

#     def zenith_angle(self, UTS_datetime = None, declinationAngle = None, hourAngle = None):
#         """Determines the zenith angle based on inputs. """

#         phi = math.radians(self.latitude)
#         try:
#             delta = math.radians(declinationAngle)
#         except:
#             delta = math.radians(declination_angle(UTS_datetime))
#         try:
#             omega = math.radians(hourAngle)
#         except:
#             omega = math.radians(self.hour_angle(UTS_datetime)) 
          
#         x = math.cos(phi)*math.cos(delta)*math.cos(omega)+math.sin(phi)*math.sin(delta)

#         if abs(math.degrees(math.acos(x))) > 90:
#             return(None)
#         else:
#             return (math.degrees(math.acos(x)))
        
#     def elevation_angle(self, UTS_datetime = None, declinationAngle = None, hourAngle = None):
#         """Determines the elevation angle based on inputs. """

#         phi = math.radians(self.latitude)
#         try:
#             delta = math.radians(declinationAngle)
#         except:
#             delta = math.radians(declination_angle(UTS_datetime))
#         try:
#             omega = math.radians(hourAngle)
#         except:
#             omega = math.radians(self.hour_angle(UTS_datetime)) 

#         x = np.arcsin(np.cos(delta)*np.cos(phi)*np.cos(omega)+np.sin(delta)*np.sin(phi))

#         if x >= 0:
#             return (np.degrees(x))

#         else:
#             return (None)        
        
#     def solar_azimuth(self, UTS_datetime = None, declinationAngle = None, hourAngle = None, zenithAngle = 'Not Provided'):

#         """The solar azimuth. Eq 1.6.6 Duffie and Beckman."""
#         phi = math.radians(self.latitude)
#         try:
#             delta = math.radians(declinationAngle)
#         except:
#             delta = math.radians(declination_angle(UTS_datetime))
#         try:
#             omega = math.radians(hourAngle)
#         except:
#             omega = math.radians(self.hour_angle(UTS_datetime))
            
#         try:
#             theta_z = np.radians(zenithAngle)
            
#         except:
#             try:
#                 theta_z = np.radians(self.zenith_angle(UTS_datetime))
#             except:
#                 return None
            
#         a = (np.cos(theta_z)*np.sin(phi)-np.sin(delta))/(np.sin(theta_z)*np.cos(phi))
        
#         if abs(a)>1:
#             a=np.sign(a)*1
            
#         x = np.degrees(np.sign(omega)*abs(np.arccos(a)))
        
#         if phi >= 0:
#             return(x)
#         elif x>0:
#             return 180-x
#         elif x<0:
#             return (-x-180)
#         else:
#             return 0

#     def solar_azimuth_topocentric(self, UTS_datetime = None, declinationAngle = None, hourAngle = None):
#         """Returns the solar azimuth at a time/latitude. Eq. 7 page section 12.6 from 'Fundimentals of Renewable
#                 Energy Processes' by Aldo Vieira da Rosa."""

#         phi = math.radians(self.latitude)
#         try:
#             delta = math.radians(declinationAngle)
#         except:
#             delta = math.radians(declination_angle(UTS_datetime))
#         try:
#             omega = math.radians(hourAngle)
#         except:
#             omega = math.radians(self.hour_angle(UTS_datetime))

#         x = np.sin(omega)/(np.sin(phi)*np.cos(omega)-np.cos(phi)*np.tan(delta))

#         if  (np.sign(omega) == 1) and (np.sign(x) == 1):
#             return (180 + np.degrees(np.arctan(x)))

#         elif (np.sign(omega) == 1) and (np.sign(x) == -1):
#             return (360 + np.degrees(np.arctan(x)))

#         elif (np.sign(omega) == -1) and (np.sign(x) == 1):
#             return (x)

#         elif (np.sign(omega) == -1) and (np.sign(x) == -1):
#             return (180 + np.degrees(np.arctan(x)))

#         else:
#             print("Something went wrong calculating solar azimuth.")
#             return(None)        
        

    #Methods that produce dataframes containing data for plotting.

#     def getSolarTime(self, frequency):
#         """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding solar time at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
#         df = pd.DataFrame(date_rng, columns=['date'])
#         df['SolarTime'] = df.apply(lambda x: self.solar_time(x['date']), axis = 1)
#         return df

#     def getHourAngle(self, frequency):
#         """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding hour angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
#         df = pd.DataFrame(date_rng, columns=['date'])
#         df['HourAngle'] = df.apply(lambda x: self.hour_angle(x['date']), axis = 1)
#         return df    

#     def getSunRise(self):
#         """Returns a dataframe with approximate hour angle at sunrise/sunset at the Position self for each day of the year."""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= 'D', closed = 'left')
#         df = pd.DataFrame(date_rng, columns=['date'])
#         df['SunRise'] = df.apply(lambda x: self.zenith_angle(x['date']), axis = 1)
#         return df
    
#     def getZenith(self, frequency):
#         """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding zenith angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
#         df = pd.DataFrame(date_rng, columns=['date'])
#         df['Zenith'] = df.apply(lambda x: self.zenith_angle(x['date']), axis = 1)
#         return df

#     def getElevation(self, frequency):
#         """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding elevation angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
#         df = pd.DataFrame(date_rng, columns=['date'])
#         df['Elevation'] = df.apply(lambda x: self.elevation_angle(x['date']), axis = 1)
#         return df            
        
#     def getAzimuth(self, frequency):
#         """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding azimuth angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
#         df = pd.DataFrame(date_rng, columns=['date'])
#         df['Azimuth'] = df.apply(lambda x: self.solar_azimuth(x['date']), axis = 1)
#         return df
    
    
#     def getTopoAzimuth(self, frequency):
#         """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding topocentric azimuth angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
#         df = pd.DataFrame(date_rng, columns=['date'])
#         df['Topocentric_Azimuth'] = df.apply(lambda x: self.solar_azimuth_topocentric(x['date']), axis = 1)
#         return df
    
#     def getAllInfo(self, frequency):
#         """Returns a dataframe with specified information for the current year for timestamps at 'frequency' at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
#         date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        
#         df = pd.DataFrame(date_rng, columns=['date'])
        
#         df['SolarTime'] = df.apply(lambda x: self.solar_time(x['date']), axis = 1)
        
#         df['HourAngle'] = df.apply(lambda x: self.hour_angle(x['date'], x['SolarTime']), axis = 1)

#         df['DeclinationAngle'] = df.date.apply(declination_angle)
        
#         df['Zenith'] = df.apply(lambda x: self.zenith_angle(x['date'], x['DeclinationAngle'], x['HourAngle']), axis = 1)
        
#         df['Elevation'] = df.apply(lambda x: self.elevation_angle(x['date'], x['DeclinationAngle'], x['HourAngle']), axis = 1)
        
#         df['Azimuth'] = df.apply(lambda x: self.solar_azimuth(x['date'], x['DeclinationAngle'], x['HourAngle'], x['Zenith']), axis = 1)
        
#         df['Topocentric_Azimuth'] = df.apply(lambda x: self.solar_azimuth_topocentric(x['date'], x['DeclinationAngle'], x['HourAngle']), axis = 1)
        
#         return df

# class SolarOrientation(object):
#     """Orientation and slope of the site."""

#     def __init__(self, orientation, slope):
#         self.orientation = orientation
#         self.slope = slope
    
#     def __repr__(self):
#         return "{0}_{1}".format(self.orientation, self.slope)


# class SolarConfig(object):
#     """Solar configuration to be considered."""

    
#     def __init__(self, position, orientation, slope):
#         """Create a solar site.
#             name        - the site id.
#             latitude    - the latitude of the site in decimal degrees.
#             longitude   - the longitude of the site as a float in decimal degrees.
#             slope       - the slope of the solar collector in decimal degrees.
#             orientation - the orientation of the collector in decimal degrees. """
        
#         self.position = position
#         self.orientation = orientation
#         self.slope = slope
        
            
#     def getCoordinates(self):
#         """Returns the coordinates of the site."""
#         return self.position


#     def getOrientation(self):
#         """Returns the orientation of the site."""
#         return [self.orientation, self.slope]
    
#     def angle_of_incidence(self, UTS_datetime):
#         """Returns the angle of incidence. Eq. 1.6.2 Duffie and Beckman."""

#         delta = math.radians(declination_angle(UTS_datetime))
#         phi = math.radians(self.position.latitude)
#         omega = math.radians(self.position.hour_angle(UTS_datetime))
#         beta = math.radians(self.slope)
#         gamma = math.radians(self.orientation)

#         a = math.sin(delta)*math.sin(phi)*math.cos(beta)
#         b = math.sin(delta)*math.cos(phi)*math.sin(beta)*math.cos(gamma)
#         c = math.cos(delta)*math.cos(phi)*math.cos(beta)*math.cos(omega)
#         d = math.cos(delta)*math.sin(phi)*math.sin(beta)*math.cos(gamma)*math.cos(omega)
#         e = math.cos(delta)*math.sin(beta)*math.sin(gamma)*math.sin(omega)

#         return (min(90,math.degrees(math.acos((a-b+c+d+e)))))





