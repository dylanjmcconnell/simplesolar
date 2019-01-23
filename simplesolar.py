import math
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt

#Functions (typically not related to the position)

def declination_angle(dt):
    """Angular position of the sun at solar noon by more accurate calculation. Eq. 1.6.1b Duffie and Beckman
    The angle between the rays of the Sun and the plane of the Earth's equator."""
    n = (dt - datetime.datetime(dt.year//4*4,1,1)).total_seconds()
    B = (n-1)*2*math.pi/31556926
    b = 0.006918 - 0.399912*math.cos(B) + 0.070257*math.sin(B) - 0.006758*math.cos(2*B)+0.000907*math.sin(2*B)-0.002697*math.cos(3*B)+0.00148*math.sin(3*B)
    return (math.degrees(b))

"""Classes are:
        Positon - has lat long and all methods relating only to those attributes
        Orientation - now redundant has only attributes relating to the slope and surface azimuth
        SolarConfig - has position and orientation attributes and methods that require both (angle of incidence)""" 

    
class Position(object):
    """Coordinate as decimal latitude and longitude."""
    def __init__(self, latitude, longitude):
        self.latitude = latitude
        self.longitude = longitude
        
    def __repr__(self):
        return "{0}_{1}".format(self.latitude, self.longitude)
    
    #Methods relating to time:
    
    def solar_time(self, UTS_datetime):
        """Returns solar time for input clock time. Eq. 1.5.2 Duffie and Beckman."""
        
        n = UTS_datetime.dayofyear+UTS_datetime.hour/24+UTS_datetime.minute/1440

        B = (n-1)*2*math.pi/365
        E = 229.2*(0.000075 + 0.001868*math.cos(B) - 0.032077*math.sin(B) - 0.014615*math.cos(2*B) - 0.04089*math.sin (2*B))

        return(UTS_datetime + pd.Timedelta(minutes = 4*(self.longitude)+E))
    
    def hour_angle(self, UTS_datetime = None, solarTime = None):
        """Takes a pd.datetime value for solar time and returns the solar hour angle to the second."""
        try:
            decimal_hour_time = solarTime.hour + solarTime.minute/60 +solarTime.second/3600
            
        except:
            solar_time = self.solar_time(UTS_datetime)
            decimal_hour_time = solar_time.hour+solar_time.minute/60+solar_time.second/3600
            
        return(360*(decimal_hour_time-12)/24)

    def sun_rise(self, UTS_datetime = None, declinationAngle = None):
        """Returns the hour angle at sunrise/sunset."""
        try:
            delta = math.radians(declinationAngle)
            phi = math.radians(self.latitude)
            
        
        except:
            delta = math.radians(declination_angle(UTS_datetime))
            phi = math.radians(self.latitude)

        return (math.degrees(math.acos(-math.tan(phi)*math.tan(delta))))                                     
                                     
                                     
                                     
    #Methods for Position:

    def zenith_angle(self, UTS_datetime = None, declinationAngle = None, hourAngle = None):
        """Determines the zenith angle based on inputs. """

        phi = math.radians(self.latitude)
        try:
            delta = math.radians(declinationAngle)
        except:
            delta = math.radians(declination_angle(UTS_datetime))
        try:
            omega = math.radians(hourAngle)
        except:
            omega = math.radians(self.hour_angle(UTS_datetime)) 
          
        x = math.cos(phi)*math.cos(delta)*math.cos(omega)+math.sin(phi)*math.sin(delta)

        if abs(math.degrees(math.acos(x))) > 90:
            return(None)
        else:
            return (math.degrees(math.acos(x)))
        
    def elevation_angle(self, UTS_datetime = None, declinationAngle = None, hourAngle = None):
        """Determines the elevation angle based on inputs. """

        phi = math.radians(self.latitude)
        try:
            delta = math.radians(declinationAngle)
        except:
            delta = math.radians(declination_angle(UTS_datetime))
        try:
            omega = math.radians(hourAngle)
        except:
            omega = math.radians(self.hour_angle(UTS_datetime)) 

        x = np.arcsin(np.cos(delta)*np.cos(phi)*np.cos(omega)+np.sin(delta)*np.sin(phi))

        if x >= 0:
            return (np.degrees(x))

        else:
            return (None)        
        
    def solar_azimuth(self, UTS_datetime = None, declinationAngle = None, hourAngle = None, zenithAngle = 'Not Provided'):

        """The solar azimuth. Eq 1.6.6 Duffie and Beckman."""
        phi = math.radians(self.latitude)
        try:
            delta = math.radians(declinationAngle)
        except:
            delta = math.radians(declination_angle(UTS_datetime))
        try:
            omega = math.radians(hourAngle)
        except:
            omega = math.radians(self.hour_angle(UTS_datetime))
            
        try:
            theta_z = np.radians(zenithAngle)
            
        except:
            try:
                theta_z = np.radians(self.zenith_angle(UTS_datetime))
            except:
                return None
            
        a = (np.cos(theta_z)*np.sin(phi)-np.sin(delta))/(np.sin(theta_z)*np.cos(phi))
        
        if abs(a)>1:
            a=np.sign(a)*1
            
        x = np.degrees(np.sign(omega)*abs(np.arccos(a)))
        
        if phi >= 0:
            return(x)
        elif x>0:
            return 180-x
        elif x<0:
            return (-x-180)
        else:
            return 0

    def solar_azimuth_topocentric(self, UTS_datetime = None, declinationAngle = None, hourAngle = None):
        """Returns the solar azimuth at a time/latitude. Eq. 7 page section 12.6 from 'Fundimentals of Renewable
                Energy Processes' by Aldo Vieira da Rosa."""

        phi = math.radians(self.latitude)
        try:
            delta = math.radians(declinationAngle)
        except:
            delta = math.radians(declination_angle(UTS_datetime))
        try:
            omega = math.radians(hourAngle)
        except:
            omega = math.radians(self.hour_angle(UTS_datetime))

        x = np.sin(omega)/(np.sin(phi)*np.cos(omega)-np.cos(phi)*np.tan(delta))

        if  (np.sign(omega) == 1) and (np.sign(x) == 1):
            return (180 + np.degrees(np.arctan(x)))

        elif (np.sign(omega) == 1) and (np.sign(x) == -1):
            return (360 + np.degrees(np.arctan(x)))

        elif (np.sign(omega) == -1) and (np.sign(x) == 1):
            return (x)

        elif (np.sign(omega) == -1) and (np.sign(x) == -1):
            return (180 + np.degrees(np.arctan(x)))

        else:
            print("Something went wrong calculating solar azimuth.")
            return(None)        
        

    #Methods that produce dataframes containing data for plotting.

    def getSolarTime(self, frequency):
        """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding solar time at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        df = pd.DataFrame(date_rng, columns=['date'])
        df['SolarTime'] = df.apply(lambda x: self.solar_time(x['date']), axis = 1)
        return df

    def getHourAngle(self, frequency):
        """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding hour angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        df = pd.DataFrame(date_rng, columns=['date'])
        df['HourAngle'] = df.apply(lambda x: self.hour_angle(x['date']), axis = 1)
        return df    

    def getSunRise(self):
        """Returns a dataframe with approximate hour angle at sunrise/sunset at the Position self for each day of the year."""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= 'D', closed = 'left')
        df = pd.DataFrame(date_rng, columns=['date'])
        df['SunRise'] = df.apply(lambda x: self.zenith_angle(x['date']), axis = 1)
        return df
    
    def getZenith(self, frequency):
        """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding zenith angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        df = pd.DataFrame(date_rng, columns=['date'])
        df['Zenith'] = df.apply(lambda x: self.zenith_angle(x['date']), axis = 1)
        return df

    def getElevation(self, frequency):
        """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding elevation angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        df = pd.DataFrame(date_rng, columns=['date'])
        df['Elevation'] = df.apply(lambda x: self.elevation_angle(x['date']), axis = 1)
        return df            
        
    def getAzimuth(self, frequency):
        """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding azimuth angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        df = pd.DataFrame(date_rng, columns=['date'])
        df['Azimuth'] = df.apply(lambda x: self.solar_azimuth(x['date']), axis = 1)
        return df
    
    
    def getTopoAzimuth(self, frequency):
        """Returns a dataframe with timestamps at 'frequency' for the current year, and the corresponding topocentric azimuth angle at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        df = pd.DataFrame(date_rng, columns=['date'])
        df['Topocentric_Azimuth'] = df.apply(lambda x: self.solar_azimuth_topocentric(x['date']), axis = 1)
        return df
    
    def getAllInfo(self, frequency):
        """Returns a dataframe with specified information for the current year for timestamps at 'frequency' at Position self. Frequency as per Pandas date_range eg 'H' for hour, '10min' for ten minute) see https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timeseries-offset-aliases"""
        
        date_rng = pd.date_range(start='1/1/%d' %pd.Timestamp.utcnow().year, end='1/1/%d' %(pd.Timestamp.utcnow().year+1), freq= frequency, closed = 'left')
        
        df = pd.DataFrame(date_rng, columns=['date'])
        
        df['SolarTime'] = df.apply(lambda x: self.solar_time(x['date']), axis = 1)
        
        df['HourAngle'] = df.apply(lambda x: self.hour_angle(x['date'], x['SolarTime']), axis = 1)

        df['DeclinationAngle'] = df.date.apply(declination_angle)
        
        df['Zenith'] = df.apply(lambda x: self.zenith_angle(x['date'], x['DeclinationAngle'], x['HourAngle']), axis = 1)
        
        df['Elevation'] = df.apply(lambda x: self.elevation_angle(x['date'], x['DeclinationAngle'], x['HourAngle']), axis = 1)
        
        df['Azimuth'] = df.apply(lambda x: self.solar_azimuth(x['date'], x['DeclinationAngle'], x['HourAngle'], x['Zenith']), axis = 1)
        
        df['Topocentric_Azimuth'] = df.apply(lambda x: self.solar_azimuth_topocentric(x['date'], x['DeclinationAngle'], x['HourAngle']), axis = 1)
        
        return df

class SolarOrientation(object):
    """Orientation and slope of the site."""

    def __init__(self, orientation, slope):
        self.orientation = orientation
        self.slope = slope
    
    def __repr__(self):
        return "{0}_{1}".format(self.orientation, self.slope)


class SolarConfig(object):
    """Solar configuration to be considered."""

    
    def __init__(self, position, orientation, slope):
        """Create a solar site.
            name        - the site id.
            latitude    - the latitude of the site in decimal degrees.
            longitude   - the longitude of the site as a float in decimal degrees.
            slope       - the slope of the solar collector in decimal degrees.
            orientation - the orientation of the collector in decimal degrees. """
        
        self.position = position
        self.orientation = orientation
        self.slope = slope
        
            
    def getCoordinates(self):
        """Returns the coordinates of the site."""
        return self.position


    def getOrientation(self):
        """Returns the orientation of the site."""
        return [self.orientation, self.slope]
    
    def angle_of_incidence(self, UTS_datetime):
        """Returns the angle of incidence. Eq. 1.6.2 Duffie and Beckman."""

        delta = math.radians(declination_angle(UTS_datetime))
        phi = math.radians(self.position.latitude)
        omega = math.radians(self.position.hour_angle(UTS_datetime))
        beta = math.radians(self.slope)
        gamma = math.radians(self.orientation)

        a = math.sin(delta)*math.sin(phi)*math.cos(beta)
        b = math.sin(delta)*math.cos(phi)*math.sin(beta)*math.cos(gamma)
        c = math.cos(delta)*math.cos(phi)*math.cos(beta)*math.cos(omega)
        d = math.cos(delta)*math.sin(phi)*math.sin(beta)*math.cos(gamma)*math.cos(omega)
        e = math.cos(delta)*math.sin(beta)*math.sin(gamma)*math.sin(omega)

        return (min(90,math.degrees(math.acos((a-b+c+d+e)))))





