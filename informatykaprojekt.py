# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:06:05 2022

@author: Julia
"""


import numpy as np
import main as m

el_grs84 = m.Transformacje(model = "wgs84")
plik = "wsp_inp.txt"

tablica = np.genfromtxt(plik, delimiter = ",", skip_header = 4)
rows,cols = np.shape(tablica)

blh = np.zeros((rows,cols))
xy2000 = np.zeros((rows,3))
xy92 = np.zeros((rows,2))
neu = np.zeros((rows,9))
az_elev_dis = np.zeros((rows,6))

tablica_ze_wsp = np.zeros((rows,13))

for i in range(rows):
    blh[i] = el_grs84.xyz2blh(tablica[i,0], tablica[i,1], tablica[i,2])
    xy2000[i] = el_grs84.u2000(blh[i,0], blh[i,1], 0.999923, 21)
    xy92[i] = el_grs84.u92(blh[i,0], blh[i,1], 0.9993)
    neu[i] = neu(tablica[i,0], tablica[i,1], tablica[i,2], tablica[0,0], tablica[0,1], tablica[0,2], blh[i,0], blh[i,1], 21)
    az_elev_dis[i] = el_grs84.azym_elew(neu[i,0], neu[i,1], neu[i,2],tablica[i,0], tablica[i,1], tablica[i,2])
    tablica_ze_wsp[i,0:3] = blh[i]
    tablica_ze_wsp[i,3:6] = xy2000[i]
    tablica_ze_wsp[i,6:8] = xy92[i]
    tablica_ze_wsp[i,8:11] = neu[i]
    tablica_ze_wsp[i,11:13] = az_elev_dis[i]
    
np.savetxt("wsp_BLH.txt", tablica_ze_wsp, delimiter=',')