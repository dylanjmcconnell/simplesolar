import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

#Functions for handling time.

def to_solar_time(UTS_datetime, longitude):
    """Returns solar time for input clock time. Eq. 1.5.2 Duffie and Beckman."""
    n = UTS_datetime.dayofyear+UTS_datetime.hour/24+UTS_datetime.minute/1440
    
    B = (n-1)*2*math.pi/365
    E = 229.2*(0.000075 + 0.001868*math.cos(B) - 0.032077*math.sin(B) - 0.014615*math.cos(2*B) - 0.04089*math.sin (2*B))
    
    return(UTS_datetime + pd.Timedelta(minutes = 4*(longitude)+E))


def hour_angle(UTS_datetime, longitude):
    """Takes a pd.datetime value for solar time and returns the solar hour angle to the second."""
    
    solar_time = to_solar_time(UTS_datetime, longitude)
    decimal_hour_time = solar_time.hour+solar_time.minute/60+solar_time.second/3600
    
    return(360*(decimal_hour_time-12)/24)


#Functions for angles

def declination_angle(UTS_datetime):
    """Angular position of the sun at solar noon by more accurate calculation. Eq. 1.6.1b Duffie and Beckman."""
    n = UTS_datetime.dayofyear+UTS_datetime.hour/24+UTS_datetime.minute/1440
    B = (n-1)*2*math.pi/365
    b = 0.006918 - 0.399912*math.cos(B) + 0.070257*math.sin(B) - 0.006758*math.cos(2*B)+0.000907*math.sin(2*B)-0.002697*math.cos(3*B)+0.00148*math.sin(3*B)
    return (math.degrees(b))

def angle_of_incidence(UTS_datetime, surface_azimuth, slope, latitude, longitude):
    """Returns the angle of incidence. Eq. 1.6.2 Duffie and Beckman."""
    
    delta = math.radians(declination_angle(UTS_datetime))
    phi = math.radians(latitude)
    omega = math.radians(hour_angle(UTS_datetime, longitude))
    beta = math.radians(slope)
    gamma = math.radians(surface_azimuth)
    
    a = math.sin(delta)*math.sin(phi)*math.cos(beta)
    b = math.sin(delta)*math.cos(phi)*math.sin(beta)*math.cos(gamma)
    c = math.cos(delta)*math.cos(phi)*math.cos(beta)*math.cos(omega)
    d = math.cos(delta)*math.sin(phi)*math.sin(beta)*math.cos(gamma)*math.cos(omega)
    e = math.cos(delta)*math.sin(beta)*math.sin(gamma)*math.sin(omega)
    
    return (min(90,math.degrees(math.acos((a-b+c+d+e)))))


def sun_rise(latitude, UTS_datetime):
    """Returns the hour angle at sunrise/sunset."""
    
    delta = math.radians(declination_angle(UTS_datetime))
    phi = math.radians(latitude)
    
    return (math.degrees(math.acos(-math.tan(phi)*math.tan(delta))))


def zenith_angle(latitude, longitude, UTS_datetime):
    
    """Returns the angle between the sun and a vertical line at that location. The angle of
            incidence between the sun and a flat plate collector. Eq. 1.6.5 Duffie and Beckman."""
    delta = math.radians(declination_angle(UTS_datetime))
    phi = math.radians(latitude)
    omega = math.radians(hour_angle(UTS_datetime, longitude)) 
    
    
    a = math.cos(phi)*math.cos(delta)*math.cos(omega)+math.sin(phi)*math.sin(delta)

    if abs(math.degrees(math.acos(a))) > 90:
        return(None)
    else:
        return (math.degrees(math.acos(a)))



def solar_azimuth(latitude, longitude, UTS_datetime):

    """The solar azimuth. Eq 1.6.6 Duffie and Beckman."""
    
    delta = math.radians(declination_angle(UTS_datetime))
    phi = math.radians(latitude)
    omega = math.radians(hour_angle(UTS_datetime, longitude)) 

    if zenith_angle(latitude, longitude, UTS_datetime) is None:
        return None
    
    else:
        theta_z = np.radians(zenith_angle(latitude, longitude, UTS_datetime))
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

def solar_azimuth_topocentric(latitude, longitude, UTS_datetime):
    """Returns the solar azimuth at a time/latitude. Eq. 7 page section 12.6 from 'Fundimentals of Renewable
            Energy Processes' by Aldo Vieira da Rosa."""
    
    delta = math.radians(declination_angle(UTS_datetime))
    phi = math.radians(latitude)
    omega = math.radians(hour_angle(UTS_datetime, longitude))    
    
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

    
def elevation_angle(latitude, longitude, UTS_datetime):
    """Determines the elevation angle based on inputs."""
    
    delta = math.radians(declination_angle(UTS_datetime))
    phi = math.radians(latitude)
    omega = math.radians(hour_angle(UTS_datetime, longitude))   
    
    x = np.arcsin(np.cos(delta)*np.cos(phi)*np.cos(omega)+np.sin(delta)*np.sin(phi))
    
    if x >= 0:
        return (np.degrees(x))
    
    else:
        return (None)
    
    
    
    
    
