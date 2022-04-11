# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:06:05 2022

@author: admin
"""


import numpy as np
import projektjeden1 as p

el_grs84 = p.Transformacje(model = "wgs84")
plik = "wsp_inp.txt"

tablica = np.genfromtxt(plik, delimiter = ",", skip_header = 4)
rows,cols = np.shape(tablica)

blh = np.zeros((rows,cols))
xy2000 = np.zeros((rows,3))
xy92 = np.zeros((rows,2))
neu = np.zeros((rows,cols))
az_elev_dis = np.zeros((rows,4))

tablica_ze_wsp = np.zeros((rows,8))
tablica_ze_wsp2 = np.zeros((rows,7))

for i in range(rows):
    blh[i] = el_grs84.xyz2blh(tablica[i,0], tablica[i,1], tablica[i,2])
    xy2000[i] = el_grs84.u2000(blh[i,0], blh[i,1], 0.999923, 21)
    xy92[i] = el_grs84.u92(blh[i,0], blh[i,1], 0.9993)
    
   # neu[i] = el_grs84.neu(tablica[i,0], tablica[i,1], tablica[i,2], tablica[i,0]+1, tablica[i,1]+1)
    tablica_ze_wsp[i,0:3] = blh[i]
    tablica_ze_wsp[i,3:6] = xy2000[i]
    tablica_ze_wsp[i,6:8] = xy92[i]
    
np.savetxt("wsp_BLH.txt", tablica_ze_wsp, delimiter=',')