import math
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import xarray as xr
import pytz
import types
import progressbar
import os


def get_files(mypath):
    """Returns a list of files that have filenames starting with "firstnine" as list of lists [fullpath, filename]."""
    fileinfo = []
    for root, dirs, files in os.walk(mypath, topdown = False):
        for name in files:
            if name[-2:] == 'nc':
                fileinfo.append([os.path.join(root,name), name])
    return (fileinfo)

#Functions relating to radiation

def test_melb():
    loc = Position(-37, 144)
    con = SolarConfig(-37, 144, 25, 0)
    df = loc.get_radiation_data('Users/felixsilberstein/Desktop/sandbox/*.nc')
    return (df)


def mean_extraterrestrial_radiation(datetime):
    """Extraterrestrial radiation incident on a plane normal to the radiation. Eq. 1.4.1b Duffie and Beckman."""
    
    doy = datetime.timetuple().tm_yday+datetime.timetuple().tm_hour/24+datetime.timetuple().tm_min/1440

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
        self.julian = self.julian()
        self.gmst = self.greenwich_mean_siderial_time()
        self.mnlong = self.set_mean_longitude()
        self.mnanom = self.set_mean_anomaly()
        self.eclong = self.set_ecliptic_longitude()
        self.obleq = self.set_oblequity()
        self.ra = self.set_right_ascension()
        self.decrad = self.set_declination()
        self.eot = self.set_equation_of_time()
        
    def __str__(self):
        return "dt = {}, julian = {}, decrad = {}".format(self.dt, self.julian, self.decrad)
    
    def julian(self):
        """Finds the julian date based on the SAM Photovoltaic Model Techinical Reference Update, Paul Gilman, Aron Dobos, Nicholas DiOrio, Janine Freeman, Steven Janzou, and David Ryberg, of National Renewable Energy Laboratory. March 2018. After Michalsky 1988."""

        utctimetuple = self.dt.utctimetuple()

        year = utctimetuple.tm_year
        jdoy = utctimetuple.tm_yday
        hour = utctimetuple.tm_hour + utctimetuple.tm_min/60

        return (32916.5 + 365*(year-1949)+(year-1949)//4+jdoy+hour/24-51545)

    def greenwich_mean_siderial_time(self):
        """Calculates the siderial time from the julian and datetime (requires timezone (pytz)) (see NREL, Michalsky, Walraven)."""
        utctimetuple = self.dt.utctimetuple()
        hour = utctimetuple.tm_hour+utctimetuple.tm_min/60
        
        a = 6.697375+0.0657098242*self.julian+(hour)
        a = a - 24*(a//24)
        if a < 0:
            a = a + 24
        return(a)
    
    
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


    def set_equation_of_time(self):
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
        
    def __str__(self):
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
    
    def get_radiation_data(self, path = '/data/marble/sandbox/jsilberstein/yearstransposed/*.nc'):
        """Gets all avaliable solar data at the location path for the location from the nearest location. Returns it as dataframe."""

        radiation_data = xr.open_mfdataset(path)
        dsloc = radiation_data.sel(longitude=self.lon, latitude=self.lat, method='nearest')
        print(datetime.datetime.now())
        df = dsloc.to_dataframe()
        radiation_data.close()
        print(datetime.datetime.now())
        return (df)

    def adjust_time(self, index_dt):
        """Adds column to data frame refering to a timestamp that is adjusted for the satelite delay. Refer to BoM meta data document"""

        if (-44<= self.lat <=-10):
            latitude = -self.lat
        else:
            print("latitude outside range of data.")
            return(None)

        dates_dict = {'1990-01-01' : datetime.datetime(1989, 12, 1, tzinfo = datetime.timezone.utc ),
                        '1993-01-01' : datetime.datetime(1993, 1, 1, tzinfo = datetime.timezone.utc ),
                        '1994-07-01' : datetime.datetime(1994, 7, 1, tzinfo = datetime.timezone.utc ),
                        '1995-06-11' : datetime.datetime(1995, 6, 11, tzinfo = datetime.timezone.utc ),
                        '2003-05-21' : datetime.datetime(2003, 5, 21, tzinfo = datetime.timezone.utc ),
                        '2005-11-01' : datetime.datetime(2005, 11, 1, tzinfo = datetime.timezone.utc ),
                        '2010-07-01' : datetime.datetime(2010, 7, 1, tzinfo = datetime.timezone.utc ),
                        '2016-03-22' : datetime.datetime(2016, 3, 22, tzinfo = datetime.timezone.utc )}

        delay = {'1990-01-01': [[45.7, 46.7, 47.7, 48.7, 49.6, 50.5, 51.2, 51.8], [38.7, 39.7, 40.7, 41.7, 42.6, 43.5, 44.2, 44.8]],
                 '1993-01-01': [[47.2, 48.2, 49.3, 50.2, 51.1, 52.0, 52.7, 53.3], [40.7, 41.7, 42.8, 43.7, 44.6, 45.5, 46.2, 46.8]],
                 '1994-07-01': [[46.7, 47.7, 48.8, 49.7, 50.6, 51.5, 52.2, 52.8], [40.5, 41.5, 42.6, 43.5, 44.4, 45.3, 46.0, 46.6]],
                 '1995-06-11': [[46.7, 47.7, 48.8, 49.7, 50.6, 51.5, 52.2, 52.8], [39.7, 40.7, 41.8, 42.7, 43.6, 44.5, 45.2, 45.8]],
                 '2003-05-21': [[39.9, 41.0, 42.0, 43.0, 43.9, 44.7, 45.5, 46.0], [27.9, 29.0, 30.0, 31.0, 31.9, 32.7, 33.5, 34.0]],
                 '2005-11-01': [[46.2, 47.2, 48.3, 49.2, 50.1, 51.0, 51.7, 52.3], [46.2, 47.2, 48.3, 49.2, 50.1, 51.0, 51.7, 52.3]],
                 '2010-07-01': [[44.7, 45.7, 46.8, 47.7, 48.6, 49.5, 50.2, 50.8], [44.7, 45.7, 46.8, 47.7, 48.6, 49.5, 50.2, 50.8]],
                 '2016-03-22': [[36.0, 36.9, 37.0, 37.9, 38.4, 38.6, 38.9, 39.1], [36.0, 36.9, 37.0, 37.9, 38.4, 38.6, 38.9, 39.1]]}

        lats = [10, 15, 20, 25, 30, 35, 40, 44]

        key = max(k for k in dates_dict if dates_dict[k] <= index_dt)

        #identifies whether the time is in column a or b
        if index_dt.hour in [18, 19, 20, 21, 23, 0, 1, 2, 3, 5, 6, 7, 8, 9, 11]:
            a = 0
        elif index_dt.hour in [22, 4, 10]:
            a = 1
        else:
            return(index_dt)

        offset = np.interp(latitude, lats, delay[key][a])

        return (index_dt + datetime.timedelta(minutes = offset))

    def adjust_time_df(self, df):
        """Applies the time adjustment to all rows in dataframe (inplace)"""
        df['time_adjusted'] = df.index.map(lambda x: self.adjust_time(pytz.timezone('UTC').localize(x.to_pydatetime())))


    def get_clearsky_df(self, df):
        """Fills a given pd.DataFrame with solar position/radiation data based on a clear sky model. DataFrame must have datetime index. Note that the isinstance checks have played up in the past."""
        df['time_adjusted'] = df.index.map(lambda x: self.adjust_time(pytz.timezone('UTC').localize(x.to_pydatetime())))
        df['ecliptic'] = df.apply(lambda x: EclipticCoordinate(x.time_adjusted), axis = 1)
        df['lmst'] = df.apply(lambda x: self.local_mean_siderial_time(x.ecliptic.gmst), axis = 1)
        df['hour_angle'] = df.apply(lambda x: self.hour_angle(x.lmst, x.ecliptic.ra), axis = 1)
        df['elevation_angle'] = df.apply(lambda x: self.elevation_angle(x.ecliptic.decrad, x.hour_angle), axis = 1)
        df['zenith'] = df.apply(lambda x: self.zenith_angle(x.ecliptic.decrad, x.hour_angle), axis = 1)
        df['azimuth'] = df.apply(lambda x: self.azimuth_angle(x.elevation_angle, x.ecliptic.decrad, x.hour_angle), axis = 1)
        df['mean_et_rad'] = df.apply(lambda x: mean_extraterrestrial_radiation(x.ecliptic.dt), axis = 1)
        if isinstance(self, SolarConfig):
            df['angle_of_incidence'] = df.apply(lambda x: self.angle_of_incidence(x.zenith, x.azimuth, x.ecliptic.decrad, x.hour_angle), axis = 1)
            df['et_radiation'] = df.apply(lambda x: incident_radiation(x.mean_et_rad, x.angle_of_incidence), axis = 1)                    
        elif isinstance(self, Position):
            df['et_radiation'] = df.apply(lambda x: incident_radiation(x.mean_et_rad, x.zenith), axis = 1)
        else:
            print("Not a SolarConfig or Position variable, please ensure class is appropriate.")


    def get_clearsky(self, date):
        """Get the clearsky inforamtion for a single point return as a dictionary input datetime.datetime."""
        ecliptic = EclipticCoordinate(date)
        lmst = self.local_mean_siderial_time(ecliptic.gmst)
        hour_angle = self.hour_angle(lmst, ecliptic.ra)
        elevation_angle = self.elevation_angle(ecliptic.decrad, hour_angle)
        zenith_angle = self.zenith_angle(ecliptic.decrad, hour_angle)
        azimuth =  self.azimuth_angle(elevation_angle, ecliptic.decrad, hour_angle)
        mean_et_rad = mean_extraterrestrial_radiation(ecliptic.dt)
        if isinstance(self, SolarConfig):
            angle_of_incidence = self.angle_of_incidence(zenith_angle, azimuth, ecliptic.decrad, hour_angle)
            et_radiation = incident_radiation(mean_et_rad, angle_of_incidence)
            info = {'ecliptic': ecliptic, 'lmst': lmst, 'hour_angle': math.degrees(hour_angle), 'elevation_angle': math.degrees(elevation_angle), 'zenith_angle': math.degrees(zenith_angle), 'azimuth_angle': math.degrees(azimuth), 'angle_of_incidence': math.degrees(angle_of_incidence), 'et_radiation': et_radiation}

        elif isinstance(self, Position):
            et_radiation = incident_radiation(mean_extraterrestrial_radiation(ecliptic.dt), zenith_angle)
            info = {'ecliptic': ecliptic, 'lmst': lmst, 'hour_angle': math.degrees(hour_angle), 'elevation_angle': math.degrees(elevation_angle), 'zenith_angle': math.degrees(zenith_angle), 'azimuth_angle': math.degrees(azimuth), 'et_radiation': et_radiation}
        else:
            print("Not a SolarConfig or Position variable, please ensure class is appropriate.")

        return(info)

    
    def checktype(self):
        """This function prints the object class. It is useful to check if the Position/SolarConfig object is behaving properly. If it stops working you can usually fix the problem by reinstallising the Position/SolarConfig object."""
        print('auto')
        if isinstance(self, SolarConfig):
            print("SolarConfig", type(self))
            
        elif isinstance(self, Position):
            print("Position", type(self))
            
        else:
            print("Not a Position/SolarConfig class")
            
        
    
    
class SolarConfig(Position):
    """SolarConfig is a subclass of Position. It contains orientation attributes in addition to latitude and longitude. Methods will determine the angle of incidence etc.
    Surface Azimuth of 0 for north facing and 180 for south facing.."""
    
    def __init__(self, latitude, longitude, slope, surface_azimuth):
        Position.__init__(self,latitude, longitude)
        self.slope = slope
        self.surface_azimuth = surface_azimuth
    
    def angle_of_incidence(self, zenith_angle, azimuth_angle, declination, hour_angle):
        """Calculates angle of incidence based to the collector."""
        sunrise = self.sunrise_hour_angle(declination)

        #if the sun is down then 0 incident angle (shadow of earth).
        if abs(hour_angle)>abs(sunrise):
            return(math.pi)


        a = math.sin(zenith_angle)*math.cos(math.radians(self.surface_azimuth)- azimuth_angle)*math.sin(math.radians(self.slope))+math.cos(zenith_angle)*math.cos(math.radians(self.slope))
        
        if a<-1:
            return (math.pi)
        elif a>1:
            return 0
        elif (-1<=a<=1):
            return math.acos(a)
        
    def hdkr_planeofarray_radiation(self, angle_of_incidence, zenith_angle, mean_et_rad, Eb, Eg, albedo = 0):
        """Calculates the diffuse radiation based on the Hay, Davies, Klucher and Reindl model. Including the reflectance which is initialised at zero (albedo). Parameters use are named after those in NREL.
        Ibh : Incident beam irradiadiance on "Plane of Array".
        Igh : Total irradiance on "Plane of Array".
        Eg  : Total horizontal irradiance (ghi)
        Eb : Beam horizontal irradiance (measured).
        Ed : Diffuse horizontal irradiance (measured).
        Rb : Ratio of beam incident to "Plane of Array" to horizontal beam.
        Ai : The anisotropy index for forward scattering circumsolar diffuse irradiance.
        f : Modulating factor for horizontal brightening correction.
        s : Horison brightening correction factor.
        
        Note that the isohor term contains both isotropic and horizonal brightening iso*(1+x)."""
        if (Eb <= 0) or (Eg <= 0):
            return None


        H = max(mean_et_rad*math.cos(zenith_angle),0)

        Ed = max(Eg - Eb*math.cos(zenith_angle), 0)

        Ibh = Eb*math.cos(zenith_angle)
        if Ibh > mean_et_rad:
            print('Beam irradiance (Ibh) > mean extraterrestrial radiation (H)')
            return (None)


        Igh = Ibh + Ed

        Rb = math.cos(angle_of_incidence)/math.cos(zenith_angle)
        
        Ai = Ibh/H
        
        if (Igh <= 0) or (Ibh <= 0):
            return None

        f = math.sqrt(Ibh/Igh)
        
        s = (math.sin(math.radians(self.slope)/2))**3
        
        cir = max(Ed*Ai*Rb, 0)
        iso = Ed*(1-Ai)*(1+math.cos(math.radians(self.slope)))/2
        isohor = iso*(1+f*s)
        
        #Ground reflected diffuse (albedo is ground reflectance, if not provided = 0)
        Ir = albedo*(Igh)*(1-math.cos(math.radians(self.slope)))/2
        
        #POA beam radiation
        Ib = max(Eb*math.cos(angle_of_incidence), 0)
        hdkr = {'Beam': Ib,
                'Isentropic' : iso,
                'Horisonal' : isohor,
                'Circumsolar' : cir,
                'Ground_reflected': Ir,
                'Total_POA' : (Ib + isohor + cir + Ir)}
        return (hdkr['Total_POA'])

    def hdkr_df(self, df):
        df['hdkr_radiation'] = df.apply(lambda x: self.hdkr_planeofarray_radiation(x.angle_of_incidence, x.zenith, x.mean_et_rad, x.dni, x.ghi), axis = 1)


