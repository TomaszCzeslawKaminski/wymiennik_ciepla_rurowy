#!/usr/bin/python
# -*- coding: utf-8 -*-

import wzory_fizyczne

m1 = 0.3 # kg/s
cp1 = 2500. # J/kgK

m2 = 0.1 # kg/s
cp2 = 4190. # J/kgK

k = 300. # W/(m^2 * K)
A = 10. # m^2

Tin1 = 400. # K
Tin2 = 300. # K

R = wzory_fizyczne.obliczenia_R(m1,cp1,m2,cp2)
print( f"R = {R: .4f}" )

NTU = wzory_fizyczne.NTU(k,A,m1,cp1,m2,cp2)
print( f"NTU = {NTU: .4f}" )

e = wzory_fizyczne.efektywnosc_wyniennika_wspol_pradowego(NTU, R)
print( f"e = {e: .4f}" )

# calc:
Tout2 = wzory_fizyczne.wymiennik_tout2(e,Tin1,Tin2)
Tout1 = wzory_fizyczne.wymiennik_tout1(e,Tin1,Tin2,m1,cp1,m2,cp2)

print( f"Tout2 = {Tout2: .1f}*C" )
print( f"Tout1 = {Tout1: .1f}*C" )
