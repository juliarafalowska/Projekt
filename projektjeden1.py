# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:28:53 2022

@author: admin
"""

from math import sqrt, atan, sin, cos, degrees, radians, tan
import numpy as np
import statistics as st
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2
    
    def xyz2blh(self, X, Y, Z):
        r = sqrt(X**2 + Y**2)
        lat_prev = atan(Z / (r * (1 - self.ecc2)))
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat_prev)) ** 2)
        h = r / cos(lat_prev) - N
        lat_next = atan((Z / r) * (((1 - self.ecc2 * N / (N + h))**(-1))))
        epsilon = 0.0000001 / 206265
        while abs(lat_prev - lat_next) < epsilon:
            lat_prev = lat_next
            N = self.a / sqrt(1 - self.ecc2 * (sin(lat_prev)) ** 2)
            h = r / cos(lat_prev) - N
            lat_next = atan((Z / r) * (((1 - self.ecc2 * N / (N + h))**(-1))))
        lat = lat_prev
        lon = atan(Y / X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat)) ** 2)
        h = r / cos(lat) - N
        return degrees(lat), degrees(lon), h
    
    def blh2xyz(self, fi, lam, h):
        N = self.a / sqrt(1 - self.ecc2 * (sin(fi)) ** 2)
        X = (N + h) * cos(fi) * cos(lam)
        Y = (N + h) * cos(fi)*sin(lam)
        Z = (N * (1 - self.ecc2) + h) * sin(fi)
        return X, Y, Z
    
    def neu(fi, l, X_sr, Y_sr, Z_sr, data):
        X=[]
        Y=[]
        Z=[]
        for line in data:
            X.append(line[1])
            Y.append(line[2])
            Z.append(line[3])
        X_sr = float(X_sr)
        Y_sr = float(Y_sr)
        Z_sr = float(Z_sr)
        dX = []
        dY = []
        dZ = []
        NEU = []
        dN = []
        dE = []
        dU = []
        for x, y, z in zip(X, Y, Z):
            delta_X = x - X_sr
            delta_Y = y - Y_sr
            delta_Z = z - Z_sr
            dX.append(delta_X)
            dY.append(delta_Y)
            dZ.append(delta_Z)
        Rt = np.matrix([((-sin(fi) * cos(l)), (-sin(fi) * sin(l)), (cos(fi))),
                ((-sin(l)), (cos(l)), (0)),
                ((cos(fi) * cos(l)), (cos(fi) * sin(l)), (sin(fi)))])
        for x, y, z in zip(dX, dY, dZ):
            d = np.matrix([x, y, z])
            d = d.T
            neu = Rt * d
            NEU.append(neu)
            dN.append(float(neu[0]))
            dE.append(float(neu[1]))
            dU.append(float(neu[2]))

        return (NEU, dN, dE, dU)
    
    def u2000(self, fi, lam, m_0, h):
        N = self.a /(sqrt(1-self.ecc2 * sin(fi)**2))
        t = tan(fi)
        n2 = self.ecc2 * cos(lam)**2
 #       lam = degrees(lam)
    
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam_0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam_0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam_0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam_0 = 24
        
        lam = radians(lam)
        lam_0 = radians(lam_0)
        l = lam - lam_0
    
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
    
    
        sigma = self.a* ((A_0*fi) - (A_2*sin(2*fi)) + (A_4*sin(4*fi)) - (A_6*sin(6*fi)))
    
        x = sigma + ((l**2)/2) * (N*sin(fi)*cos(fi)) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*cos(fi)) * (1 + ((((l**2)/6) * (cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))

        x00 = round(x * m_0, 3)
        y00 = round(y * m_0 + (s*1000000) + 500000, 3)   
    
        return x00, y00, h 
    
    def u92(self, fi, lam, m_0):
        e_2 = self.ecc2/(1-self.ecc2)
        N = self.a/(sqrt(1-self.ecc2 * sin(fi)**2))
        t = tan(fi)
        n2 = e_2 * cos(lam)**2
        lam_0 = radians(19)
        l = lam - lam_0
    
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
    
        sigma = self.a* ((A_0*fi) - (A_2*sin(2*fi)) + (A_4*sin(4*fi)) - (A_6*sin(6*fi)))
    
        x = sigma + ((l**2)/2) * (N*sin(fi)*cos(fi)) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*cos(fi)) * (1 + ((((l**2)/6) * (cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))

        x92 = round(x * m_0 - 5300000, 3)
        y92 = round(y * m_0 + 500000, 3)   
    
        return x92, y92 