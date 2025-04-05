import scipy.optimize as so
import math


###################
# Stale           #
###################

# wartosci opalowe gazow [kJ/um3]

Q_w_CH4_kJ_mu3 = 35818
Q_w_C2H6_kJ_mu3 = 63748
Q_w_C3H8_kJ_mu3 = 91251

Q_w_H2_kJ_mu3 = 10790
Q_w_CO_kJ_mu3 = 12690

# cieplo wlasciwe gazow w warunkach normalnych 0 C [kJ/um3K]

C_p_CH4_kJ_um3K = 1.5499
C_p_C2H6_kJ_um3K = 2.2098
C_p_C3H8_kJ_um3K = 3.0484
C_p_CO2_kJ_um3K = 1.5998
C_p_N2_kJ_um3K = 1.2946
C_p_O2_kJ_um3K = 1.3095

C_p_H2_kJ_um3K = 1.28

C_p_CO_kJ_um3K = 1.3038
C_p_CO_kJ_kgK = 1.043

# gestosc gazow w warunkach normalnych 0 C [kg/um3]

g_CH4_kg_um3K = 0.717
g_C2H6_kg_um3K = 1.357
g_C3H8_kg_um3K = 2.004
g_CO2_kg_um3K = 1.977
g_N2_kg_um3K = 1.251
g_O2_kg_um3K = 1.429
g_SO2_kg_um3K = 2.9263
g_H2O_kg_um3K = 0.8040
g_H2_kg_um3K = 0.0899
#g_CO_kg_um3K = 0.0


######################
# Funkcje pomocnicze #
######################

def propocjonalnosc(x1, y1, x2, y2, dana):
	return y1+(dana-x1)*(y2-y1)/(x2-x1)

def moc_kotla(wartosc_opalowa_paliwa, ilosc_podawanego_paliwa):
# wartosc_opalowa_paliwa, MJ/kg
# ilosc_podawanego_paliwa, kg/h
	return (wartosc_opalowa_paliwa * ilosc_podawanego_paliwa)/3600. * 1000. # kW

def ilosc_potrzebnego_paliwa(wartosc_opalowa_paliwa, moc):
# wartosc_opalowa_paliwa, MJ/kg
# moc, kW
	return ((moc/1000.)*3600.)/wartosc_opalowa_paliwa # kg/h

###################
# Wlasciwosci     #
###################

def lepkosc_dynamiczna(lepkosc_kinematyczna, gestosc):
        return lepkosc_kinematyczna * gestosc # Pa*s

#*****************************************************************************
# powietrze
#*****************************************************************************


def gestosc__powietrze(T):

# T, [C]
# gestosc, [kg/m3]

# approximation range: -50 - 1200 C
	a6 = 3.51339e-18	     #  3.356e-19    (9.553%)
	a5 = -1.49119e-14	     #  1.176e-15    (7.887%)
	a4 = 2.55675e-11	     #  1.547e-12    (6.052%)
	a3 = -2.31452e-08	     #  9.456e-10    (4.085%)
	a2 = 1.24366e-05	     #  2.683e-07    (2.157%)
	a1 = -0.00447181	     #  3.088e-05    (0.6906%)
	a0 = 1.28977		     #  0.001049     (0.08134%)

	return a6*T**6+a5*T**5+a4*T**4+a3*T**3+a2*T**2+a1*T+a0


def cieplo_wlasciwe__powietrze(T):

# T, [C]
# cieplo_wlasciwe (heat_capacity), J/(um3*K)

# approximation range: 120 - 1200 C
	a4 = 7.19829e-14	#  3.759e-14    (52.22%)
	a3 = -2.59827e-10	#  9.809e-11    (37.75%)
	a2 = 2.4375e-07		#  9.809e-11    (37.75%)
	a1 = 0.000140147	#  2.928e-05    (20.89%)
	a0 = 0.989717		#  0.003059     (0.309%)


	if -50 <= T < -20:
		cieplo_wlasciwe = 1.013
	elif -25 <= T < -1:
		cieplo_wlasciwe = 1.009;
	elif -1 <= T < 69:
		cieplo_wlasciwe = 1.005;
	elif 69 <= T < 120:
		cieplo_wlasciwe = 1.009;
	else:
		cieplo_wlasciwe = a4*T**4+a3*T**3+a2*T**2+a1*T+a0
	

	return(cieplo_wlasciwe * 1000)



def cieplo_wlasciwe_objetosciowe__powietrze(cieplo_wlasciwe_powietrza):
# volumetric_heat_capacity

	k = 1.4 # Heat capacity ratio for powietrze, c_p/c_v=k

	return cieplo_wlasciwe_powietrza/k



def lepkosc_kinematyczna__powietrze(T):
# kinematic_viscosity
	return (3.90204e-09*T**3 + 6.61374e-05*T**2 + 0.0976732 * T+ 13.3477) * 10e-7

def lepkosc_dynamiczna__powietrze(T):
# dynamic_viscosity
	return lepkosc_kinematyczna__powietrze(T) * gestosc__powietrze(T)


def przewodnosc_cieplna__powietrze(T):
# thermal_conductivity
	a0=4.13737e-3
	a1=8.03174e-5
	a2=-1.53e-8

	return a2*T**2 + a1*T + a0

#*****************************************************************************
# N2
#*****************************************************************************

def gestosc__N2(T):
# T - temperatura, C
        if T < 700:
                return 1.24732880698351-0.00399758605838068*T+8.3403151702e-6*T**2-8.83188189921944e-9*T**3+3.51012077964067e-12*T**4 # kg/m3
        else:
                return  -2.799999999999999e-04*T+5.479999999999998e-01 # kg/m3

def cieplo_wlasciwe__N2(T):
# T - temperatura, C
        return (1.0303666343356-1.77911759272111e-5*T+4.83875975705508e-7*T**2-3.19947115606659e-10*T**3+2.54549443443642e-14*T**4)*1000. # J/(kg*K)

def przewodnosc_cieplna__N2(T):
# T - temperatura, C
        return (24.2406559740494+0.0785523217396528*T-3.07499849822791e-05*T**2)*1e-3 # W/(m*K)

def lepkosc_kinematyczna__N2(T):
# T - temperatura, C
        return (12.9242386015498+0.0911714122664745*T+7.31026611401453e-05*T**2)*1e-6 # m2/s

#*****************************************************************************
# O2
#*****************************************************************************

def gestosc__O2(T):
# T - temperatura, C
        if T < 700:
                return  1.42616491974594-0.00452967951342485*T+9.30278217425703e-06*T**2-9.71630536323288e-09*T**3+3.82335279157188e-12*T**4 # kg/m3
        else:
                 return -3.200000000000001e-04*T+6.260000000000001e-01 # kg/m3

def cieplo_wlasciwe__O2(T):
# T - temperatura, C
        return (0.913645145646257+0.000143750122247202*T+6.9131594830604e-07*T**2-1.15257809741335e-09*T**3+5.25872976142044e-13*T**4)*1000. # J/(kg*K)

def przewodnosc_cieplna__O2(T):
# T - temperatura, C
        return (24.1976851684988+0.0872903646302637*T-2.52668799183036e-05*T**2)*1e-3 # W/(m*K)

def lepkosc_kinematyczna__O2(T):
# T - temperatura, C
        return (13.2302126509281+0.0936598185859315*T+7.72843455277227e-05*T**2)*1e-6 # m2/s

#*****************************************************************************
# CO2
#*****************************************************************************

def gestosc__CO2(T):
# T - temperatura, C
        if T < 600.:
                return  1.968-0.00655845238095236*T+1.52245634920634e-05*T**2-1.81992063492062e-08*T**3+7.98809523809518e-12*T**4 # kg/m3
        else:
                return  -4.875000000000000e-04*T+9.104999999999999e-01 # kg/m3

def cieplo_wlasciwe__CO2(T):
# T - temperatura, C
        return (0.816+0.00109185714285715*T-1.20674603174609e-06*T**2+9.14285714285822e-10*T**3-3.25396825396883e-13*T**4)*1000. # J/(kg*K)

def przewodnosc_cieplna__CO2(T):
# T - temperatura, C
        return (14.3508158508158+0.0889172494172495*T-1.68414918414919e-05*T**2)*1e-3 # W/(m*K)

def lepkosc_kinematyczna__CO2(T):
# T - temperatura, C
        return (7.05874125874122+0.0498752913752914*T+5.9079254079254e-05*T**2)*1e-6 # m2/s

#*****************************************************************************
# SO2
#*****************************************************************************

def cielpo_wlasciwe__SO2(T):

	a2 = -3.82576e-07	#  7.501e-09    (1.961%)
	a1 = 0.000868621	#  7.012e-06    (0.8073%)
	a0 = 1.73785		#  0.001355     (0.07798%)

	return (a2*pow(T,2) + a1*T + a0)*1000 # Cp_CO2[J/um3K]

#*****************************************************************************
# Shell Thermia Oil B
#*****************************************************************************

def gestosc__Shell_Thermia_Oil_B(T):
        return -0.650*(T+273.15) + 1053.589 # kg/m3

def cieplo_wlasciwe__Shell_Thermia_Oil_B(T):
        return 3.594*(T+273.15) + 830.934 # [J/(kg*K)]

def przewodnosc_cieplna__Shell_Thermia_Oil_B(T):
        return -0.0001*(T+273.15) + 0.1555 

def lepkosc_kinematyczna__Shell_Thermia_Oil_B(T):
# od 0 do 300 C
        if T>0.: lepkosc = -5.125*T+230
        if T>40.: lepkosc = -0.338333*T+38.5333
        if T>100.: lepkosc = -0.035*T+8.2
        if T>200.: lepkosc = -0.007*T+2.6
        return lepkosc # mm2/s

#*****************************************************************************
# Woda
#*****************************************************************************

def woda_gestosc(T):
# T - temperatura, C
	return -0.0035*T**2 - 0.0435*T + 999.9 # kg/m3

def woda_cieplo_wlasciwe(T):
# T - temperatura, C
	return 0.021875*T**2 - 2.05*T + 4226 # J/(kg*K)

def woda_przewodnosc_ciepla(T):
# T - temperatura, C
	return -7.083e6*T**2 + 0.001952*T + 0.558

def woda_lepkosc_dynamiczna(T):
# T - temperatura, C
	return 1793.636e-6 * 0.968**T # Pa*s

#*****************************************************************************
# Wlasciwosci przy zalozeniu addytywnosci
#*****************************************************************************

def gestosc(uN2, uO2, uCO2, T):
# T - temperatura, C
        return uN2*gestosc__N2(T) + uO2*gestosc__O2(T) + uCO2*gestosc__CO2(T) # kg/m3

def cieplo_wlasciwe(uN2, uO2, uCO2, T):
# T - temperatura, C
        return uN2*cieplo_wlasciwe__N2(T) + uO2*cieplo_wlasciwe__O2(T) + uCO2*cieplo_wlasciwe__CO2(T) # J/(kg*K)

def lepkosc_kinematyczna(uN2, uO2, uCO2, T):
# T - temperatura, C
        return uN2*lepkosc_kinematyczna__N2(T) + uO2*lepkosc_kinematyczna__O2(T) + uCO2*lepkosc_kinematyczna__CO2(T) # m2/s

def lepkosc_dynamiczna(uN2, uO2, uCO2, T):
# T - temperatura, C
        return lepkosc_kinematyczna(uN2, uO2, uCO2, T) * gestosc(uN2, uO2, uCO2, T) # Pa*s

def przewodnosc_cieplna(uN2, uO2, uCO2, T):
# T - temperatura, C
        return uN2*przewodnosc_cieplna__N2(T) + uO2*przewodnosc_cieplna__O2(T) + uCO2*przewodnosc_cieplna__CO2(T) # W/(m*K)

def Liczba_Prandtla(uN2, uO2, uCO2, T):
# T - temperatura, C
        return cieplo_wlasciwe(uN2, uO2, uCO2, T) * lepkosc_dynamiczna(uN2, uO2, uCO2, T) / przewodnosc_cieplna(uN2, uO2, uCO2, T)

###################
# Geometria       #
###################

def obliczanie_powierzchni_kola(srednica):
        return math.pi * (srednica/2.)**2.

def obliczanie_obwodu_kola(srednica):
        return math.pi * srednica

def obliczanie_powierzchni_zwojnicy(dlugosc, srednica):
        return dlugosc * obliczanie_obwodu_kola(srednica)

###################
# Przeplyw        #
###################

def liczba_Reynoldsa(w, d, v):
# w - predkosc charakterystyczna plynu, [m/s]
# d - wymiar charakterystyczny zagadnienia, [m]
# v - lepkosc kinematyczna plynu, [m2/s]
        return w * d / v

def obliczanie_predkosci_plynu_z_przeplywu_masowego(przeplyw_masowy, powierzchnia, gestosc):
        return przeplyw_masowy/(gestosc*powierzchnia)

def obliczanie_predkosci_plynu_z_przeplywu_masowego_dla_rury(przeplyw_masowy, srednica, gestosc):
        powierzchnia = math.pi * (srednica/2.)**2.
        return przeplyw_masowy/(powierzchnia*gestosc) # m/s

def Liczba_Nusselta_dla_przeeplywu_laminarnego(Re, lepkosc_f, lepkosc_w, Pr, d, L):
        return 1.86*(lepkosc_f/lepkosc_w)**0.14*Re**(1/3)*Pr**(1/3)*(d/L)**(1/3)

def Liczba_Nusselta_dla_przeeplywu_przejsciowego(Re, Pr_f, Pr_w):
# liczba Nusselta dla liczby Reynoldsa w zakresie 2000 < Re_f < 10000 - Gogol W. Wymiana ciepla, tablice i wykresy, 1974
# Pr_f - liczba Prandtla dla medium
# Pr_w - liczba Prandtla dla medium w temperaturze scianki

# aproksymacja liczby K0 od Re z danych na stronie 104 Gogol W. Wymiana ciepla, tablice i wykresy, 1974
        Re *= 1e-3
        K0 = -48.1348540844692+49.3113660131855*Re-18.9841921613281*Re**2+4.14258767842297*Re**3-0.492461201943884*Re**4+0.0300764683151742*Re**5-0.00073918520187297*Re**6

        return K0 * Pr_f**0.43 * (Pr_f/Pr_w)**0.25

def wspolczynnik_przejmowania_ciepla(Liczba_Nusselta, przewodnosc_cieplna, wymiar_charakterystyczny):
        return Liczba_Nusselta*przewodnosc_cieplna/wymiar_charakterystyczny

def obliczanie_predkosci_maksymalnej_z_Re_dla_przeplywu_turbulentnego(Re, w_sr):
	if Re > 4000:
		st = 0.791
	
	if Re > 10000:
		st = 0.8

	if Re > 40000:
		st = 0.811

	if Re > 100000:
		st = 0.818

	if Re > 400000:
		st = 0.827

	if Re > 1000000:
		st = 0.841

	if Re > 4000000:
		st = 0.866

	return w_sr/st
	
	

###################
# Wymiana ciepla  #
###################

def strumien_ciepla(strumien_masy, cieplo_wlasciwe, temperatura):
# strumien_masy, [kg/s]
# cieplo_wlasciwe, [J/(kg*K)]
# temperatura, [K]
	return strumien_masy*cieplo_wlasciwe*temperatura

def strumien_ciepla_odebrany(strumien_masy, cieplo_wlasciwe, temperatura_na_wejsciu, temperatura_na_wyjsciu):
# strumien_masy, [kg/s]
# cieplo_wlasciwe, [J/(kg*K)]
# temperatura, [K]
	return strumien_masy*cieplo_wlasciwe*(temperatura_na_wyjsciu-temperatura_na_wejsciu)

def cieplo_odebrane_z_temperatur(strumien_masy, cieplo_wlasciwe, temperatura_na_wejsciu, temperatura_na_wyjsciu):
# strumien_masy, [kg/s]
# cieplo_wlasciwe, [J/(kg*K)]
# temperatura, [K]
	return abs(strumien_masy*cieplo_wlasciwe*(temperatura_na_wyjsciu-temperatura_na_wejsciu))

def potrzebna_ilosc_plynu_dla_odebrania_strumienia_ciepla(strumien_siepla, c_p, zmiana_tyemperatury):
# c_p = 4189.9 - cieplo wlasciwe wody J/(kg*K)
        return strumien_siepla/(c_p*zmiana_tyemperatury) #  strumien masy, kg/s

def strumien_ciepla_emitowany_przez_bryle_gazu_na_skutek_radiacji(A,e_w,o,e_g,T_g,a_g,T_w):
# 
	return A*ei_w*o*(e_g*T_g**4-a_g*T_w**4)

def efektywna_zdolnosc_emisyjna_sciany(e_w):
	return (e_w+1)/2

def strumien_ciepla_konwekcja(alpha,A,T_g,T_w):
	return alpha*A*(T_g-T_w)

def temperatura_srednia(T1,T2):
	return (T1+T2)/2

def temperatura_scianki(T_c,z,T_g):
	return T_c+z(T_g-T_c)

def temperatura_na_wyjsciu_z_rury(t_w, t_f0, alpha, d, m, c_pf, x):
# t_w - temperatura scianki, C
# t_fo - temperatura medium na wlocie, C
# d - srednica, m
# m - strumien masy, kg/s
# c_pf - cieplo wlasciwe
# x - odleglosc, m
        return t_w + (t_f0 - t_w)* math.exp(-(alpha*math.pi*d*x)/(m*c_pf))


# empiryczne

def wspolczunnik_przenikania_ciepla_dla_rury_w_budynku(temperatura_powierzchni=20., temperatura_powietrza=20.):
	return 9.4+0.052*(temperatura_powierzchni-temperatura_powietrza) # W/(m2*K)

###################
# Radiacja        #
###################

# emisyjnosci z ksiazki teoria wymiany ciepla, Jan Madejski

def emisyjnosc_CO2(p, L, T):
# T - temperatura, C - aproksymacja 800 - 1400 C
# p - cisnienie parcjalne, kPa
# L - srednia droga promieni, m

        if p*L > 0.1:
                n = 0.614
                a = 0.08697
                b = -0.04108

        if p*L > 0.93:
                n = 0.391
                a = 0.07814
                b = -0.03321

        if p*L > 4.:
                n = 0.374
                a = 0.07613
                b = -0.03038

        if p*L > 10.:
                n = 0.314
                a = 0.07791
                b = -0.02573

        if p*L > 70.:
                n = 0.310
                a = 0.07350
                b = -0.02081

        k = a + b*(T/1000.)

        return 1.-math.exp(-k*(p*L)**n)

def emisyjnosc_H2O(p, L, T):
# T - temperatura, C - aproksymacja 800 - 1400 C
# p - cisnienie parcjalne, kPa
# L - srednia droga promieni, m

        if p*L > 0.1:
                n = 0.945
                a = 0.04433
                b = -0.02552

        if p*L > 0.93:
                n = 0.814
                a = 0.03892
                b = -0.02027

        if p*L > 4.:
                n = 0.692
                a = 0.04210
                b = -0.01979

        if p*L > 10.:
                n = 0.530
                a = 0.05729
                b = -0.02375

        if p*L > 70:
                n = 0.395
                a = 0.09700
                b = -0.03809

        k = a + b*(T/1000.)

        return 1.-math.exp(-k*(p*L)**n)

def emisyjnosc_spalin(L, p_CO2, p_H2O, T):
# T - temperatura, C - aproksymacja 800 - 1400 C
# p - cisnienie parcjalne, kPa
# L - srednia droga promieni, m

        e_CO2 = 0.
        e_H2O = 0.
        beta = 0.

        if p_CO2 != 0:
                e_CO2 = emisyjnosc_CO2(p_CO2, L, T)

        if p_H2O != 0:
                e_H2O = emisyjnosc_H2O(p_H2O, L, T)
                beta = 1.+(0.6225-0.1346*math.log10(p_H2O*L))*(p_H2O/100.)**0.86
        
        return e_CO2 + beta*e_H2O - e_CO2*beta*e_H2O


def absorbcyjnosc_spalin(L, p_CO2, p_H2O, T_f, T_w):

        a_CO2 = 0.
        a_H2O = 0.
        beta = 0.

        if p_CO2 != 0:
                a_CO2 = emisyjnosc_CO2(p_CO2*T_w/T_f, L, T_w) * (T_f/T_w)**0.65

        if p_H2O != 0:
                a_H2O = emisyjnosc_H2O(p_H2O*T_w/T_f, L, T_w) * (T_f/T_w)**0.45
                beta = 1.+(0.6225-0.1346*math.log10(p_H2O*L))*(p_H2O/100.)**0.86

        return a_CO2 + beta*a_H2O - a_CO2*beta*a_H2O


def gestosc_radiacyjnego_strumienia_ciepla(emisyjnosc_scianki, T_g, T_w, L, p_CO2, p_H2O):

        e_e = (emisyjnosc_scianki + 1.)/2.

        SSB = 5.67037321e-8 # stala_promieniowanie_ciala_doskonale_czarnego, W/(m2*K2)

        e_s = emisyjnosc_spalin(L, p_CO2, p_H2O, T_g)
        a_s = absorbcyjnosc_spalin(L, p_CO2, p_H2O, T_g, T_w)

        return e_e * SSB * (e_s*(T_g+273.15)**4.-a_s*(T_w+273.15)**4.)

def wspolczynnik_przenikania_ciepla_radiacyjny(q, T_w, T_g):
        if q==0.:
                return 0

        return q/(abs(T_g-T_w))


###################
# Spalanie        #
###################

#-----------------------------------------
# zapotrzebowanie powietrza

def teoretyczne_zapotrzebowanie_powietrza__wegiel_kamienny(C=0.,H=0.,S=0.,O=0.):
# sklad w %masowych
# wegiel kaminenny np: C=58,5%, H=3,9%, O=9,6%, N=1%, S=1,3%, W=10,3%
        return (8.89*C+26.7*H+3.33*(S-O))/100.

def teoretyczne_zapotrzebowanie_powietrza__gaz(CH4=0,C2H6=0,C3H8=0,CO2=0,N2=0,O2=0,H2=0,CO=0,H2S=0): # %mol
        return 0.0476*(0.5*CO+0.5*H2+2.*C2H6+1.5*H2S+2.*CH4+3.5*C2H6+5.*C3H8-O2) # [um3/kg]

def teoretyczne_zapotrzebowanie_powietrza__gaz_ziemny__wzor_mniej_dokladny(wartosc_opalowa):
        return 1.1*wartosc_opalowa/4186.8 # [um3/kg]

#-----------------------------------------
# wartosc opalowa

def wartosc_opalowa_gazu(CH4=0,C2H6=0,C3H8=0,CO=0,H2=0): # udzialy objetosciowe skladnikow gazu
        return (35818.*CH4+63748.*C2H6+91251.*C3H8+10790.*H2+12690.*CO)*1000. # J/um3

#-----------------------------------------

def cieplo_wlasciwe_gazu_w_20C(CH4=0,C2H6=0,C3H8=0,CO=0,H2=0,CO2=0,N2=0,O2=0): # udzialy objetosciowe skladnikow gazu
        return (1.5499*CH4+2.2098*C2H6+3.0484*C3H8+1.5998*CO2+1.2946*N2+1.3059*O2+1.28*H2)*1000. # J/(um3*K)
        
#-----------------------------------------

def entalpia_gazu_w_temperaturze_20C(CH4=0,C2H6=0,C3H8=0,CO=0,H2=0,CO2=0,N2=0,O2=0): # udzialy objetosciowe skladnikow gazu
        return cieplo_wlasciwe_gazu_w_20C(CH4,C2H6,C3H8,CO,H2,CO2,N2,O2)*20.

#-----------------------------------------

def teoretyczna_objetosc_spalin_suchych__gaz_ziemny(CO2=0,CO=0,H2S=0,N2=0,CH4=0,C2H6=0,C3H8=0,V_pow=0):
        return 0.01*(CO2+CO+H2S+N2+CH4+2.*C2H6+3.*C3H8)+0.79*V_pow

def teoretyczna_objetosc_spalin_mokrych__gaz_ziemny(V_s_t_sp=0,H2S=0,H2=0,CH4=0,C2H6=0,C3H8=0, d = 0.124):
        return V_s_t_sp + 0.01*(H2S+H2+2.*CH4+3.*C2H6+4.*C3H8)+d*0.01

def rzeczywista_objetosc_spalin_suchych__gaz_ziemny(CO2=0, CO=0, H2S=0, N2=0, O=0, V_pow=0, h=0):
        return 0.01*(CO2+CO+H2S+N2+2.3*CH4+2.*C2H6+3.*C3H8)+0.79*h*V_pow

def rzeczywista_objetosc_spalin_mokrych__gaz_ziemny(CO2=0, CO=0, H2S=0, N2=0, O=0, V_pow=0, W=0, d=0, h=0):
        return rzeczywista_objetosc_spalin_suchych__gaz_ziemny(C, H, S, O, N, V_pow,h)+0.01*(2.*CH4+3.*C2H6+4.*C3H8)+0.0161*h*V_pow

#-----------------------------------------

def teoretyczna_objetosc_spalin_suchych__wegiel_kamienny(C=0, H=0, S=0, O=0, N=0, V_pow=0):
# sklad w %masowych
# wegiel kaminenny np: C=58,5%, H=3,9%, O=9,6%, N=1%, S=1,3%, W=10,3%
        return 1.867*C/100.+0.7*S/100.+0.8*N/100.+0.79*V_pow # [um3/kg]

def teoretyczna_objetosc_spalin_mokrych__wegiel_kamienny(C=0, H=0, S=0, O=0, N=0, V_pow=0, W=0, d=10.):
# d=10g/um3 - zawartosc wilgoci w powietrzu
        return teoretyczna_objetosc_spalin_suchych__wegiel_kamienny(C, H, S, O, N, V_pow) + 11.2*H/100.+1.24*W/100.+0.00124*d*V_pow # [um3/kg]

def rzeczywista_objetosc_spalin_suchych__wegiel_kamienny(C=0, H=0, S=0, O=0, N=0, V_pow=0, h=0):
        return 1.867*C/100.+0.7*S/100.+0.8*N/100.+0.79*h*V_pow+0.21*(h-1.)*V_pow # [um3/kg]

def rzeczywista_objetosc_spalin_mokrych__wegiel_kamienny(C=0, H=0, S=0, O=0, N=0, V_pow=0, W=0, d=0, h=0):
        return rzeczywista_objetosc_spalin_suchych__wegiel_kamienny(C, H, S, O, N, V_pow,h) + 11.2*H/100.+1.24*W/100.+0.00124*d*h*V_pow # [um3/kg]

#-----------------------------------------

def temperatura_teoretyczna_spalin(I,V_rz_sp,c_p):
# I - teoretyczna entalpia spalin, kJ/um3,
# V_rz_sp - rzeczywista objetosc spalin mokrych, um3/um3,
# c_p - cieplow walsciwe spalin, kJ/(um3*K)
        return  I / (V_rz_sp * c_p) # [C]


##########################
# Przeliczanie jednostek #
##########################

def ppm_do_procent_molowy(ppm): # ppm
        return ppm*10**-4 # [%mol]

def cal_do_mm(cal):
        return 25.4*cal

def cal_do_m(cal):
        return 0.0254*cal

##########################
#   Rurka spietrzajaca   #
##########################

def liczba_Re_dla_rurki(srednica_odbioru_cisnienia, pretkosc, lepkosc_kinetyczna):
        Re = srednica_odbioru_cisnienia*pretkosc/lepkosc_kinetyczna
        if Re>=200:
                return Re
        else:
                print( "Re=", Re, "- zamala liczba Re, rurka nie spelnia wymagan, musi byc >= 200" )
                return 0

def obliczanie_predkosci_dla_rurki(alpha, liczba_ekspansji, roznica_cisnien, gestosc):
# alpha - wspolczynnik wzorcowania dla rorek znormalizowanych praktycznie wynosi 1
# roznica_cisnien - miedzy cisnieniem calkowitym a statycznym = cisnienie dynamiczne
        return alpha*(1.-liczba_ekspansji)*(2.*roznica_cisnien/gestosc)**0.5

def obliczanie_liczby_ekspansji(wykladnik_izentropy, cisnienie_statyczne, roznica_cisnien):
# wykladnik_izentropy = 1,4 - dla powietrza
        p1 = roznica_cisnien/cisnienie_statyczne
        k1 = 1./wykladnik_izentropy
        k2 = (wykladnik_izentropy-1.)/(6.*wykladnik_izentropy**2.)

        return 1.-(1.-k1*p1+k2*p1**2.)**0.5


##########################
#   Spadek cisnienia     #
##########################

def wspolczynnik_strat_liniowych__przeplyw_laminarny(Re):
        return 64./Re

def wspolczynnik_strat_liniowych__przeplyw_turbulentny__scianki_hydraulicznie_gladkie(Re):
# wzor dla scianki hydraulicznie gladkiej - 3000 < Re < 80000
        return 0.316/(Re**0.25)

def wzor_Colebrooka_Whitea(Re, chropowatosc_zwgledna, wartosc_startowa=0.089):
        return so.fsolve(lambda h: 1.+(2.*math.log10(2.51/(Re*math.sqrt(h))+chropowatosc_zwgledna/3.71))*math.sqrt(h), wartosc_startowa)
   

def spadek_cisnien(wspolczynnik_strat_liniowych, dlugosc_rury, srednica, gestosc, predkosc):
        return wspolczynnik_strat_liniowych*dlugosc_rury*gestosc*(predkosc**2.)/(2.*srednica)

##########################
#   chemia               #
##########################

def objetosc_jednego_mola_gazu_doskonalego(T=0., p=101325.):
# T - temperatura w C
	return (8.314472*(T+273.15))/p # m3

def obliczanie_masy_jednego_mola(wzor):
	if wzor == 'CO2': return 12.0107+2*15.9994 # g
	if wzor == 'CO': return 12.0107+15.9994 # g
	if wzor == 'O2': return 2*15.9994 # g
	if wzor == 'N2': return 2*14.0067 # g
	if wzor == 'H2': return 2*1.00794 # g
	if wzor == 'SO2': return 32.065+2*15.9994 # g
	if wzor == 'CH4': return 12.0107+4*1.00794 # g

def obliczanie_gestosci_gazu(uCO2=0., uCO=0., uO2=0., uN2=0., uH2=0.,uSO2=0., uCH4=0., T=0., p=101325.):
	V = objetosc_jednego_mola_gazu_doskonalego(T, p)
	qCO2 = (obliczanie_masy_jednego_mola('CO2')/1000.)/V # kg/m3 - nie mylic z um3 bo dla V moze byc inna temperatura i cisnienie niz normalne
	qCO = (obliczanie_masy_jednego_mola('CO')/1000.)/V # kg/m3 - nie mylic z um3 bo dla V moze byc inna temperatura i cisnienie niz normalne
	qO2 = (obliczanie_masy_jednego_mola('O2')/1000.)/V # kg/m3 - nie mylic z um3 bo dla V moze byc inna temperatura i cisnienie niz normalne
	qN2 = (obliczanie_masy_jednego_mola('N2')/1000.)/V # kg/m3 - nie mylic z um3 bo dla V moze byc inna temperatura i cisnienie niz normalne
	qH2 = (obliczanie_masy_jednego_mola('H2')/1000.)/V # kg/m3 - nie mylic z um3 bo dla V moze byc inna temperatura i cisnienie niz normalne
	qSO2 = (obliczanie_masy_jednego_mola('SO2')/1000.)/V # kg/m3 - nie mylic z um3 bo dla V moze byc inna temperatura i cisnienie niz normalne
	qCH4 = (obliczanie_masy_jednego_mola('CH4')/1000.)/V # kg/m3 - nie mylic z um3 bo dla V moze byc inna temperatura i cisnienie niz normalne

	return uCO2*qCO2 + uCO*qCO + uO2*qO2 + uN2*qN2 + uH2*qH2 + uSO2*qSO2 + uCH4*qCH4

###################################
# wymienniki                      #
###################################

def obliczenia_R(m1,cp1,m2,cp2):
        if cp1*m1 < cp2*m2:
               return cp1*m1/(cp2*m2)
        if cp1*m1 > cp2*m2:
               return cp2*m2/(cp1*m1)

def NTU(k,A,m1,cp1,m2,cp2):
        if cp1*m1 > cp2*m2:
               return k*A/(cp2*m2)
        if cp1*m1 < cp2*m2:
               return k*A/(cp1*m1)

def efektywnosc_wyniennika_wspol_pradowego(NTU, R):
        return (1.-math.exp(-NTU*(1.+R)))/(1.+R)

def efektywnosc_wyniennika_przeciw_pradowego(NTU, R):
        return (1.-math.exp(-NTU*(1.-R)))/(1.-R*math.exp(-NTU*(1.-R)))

def wymiennik_tout2(e,Tin1,Tin2):
        return Tin2+e*(Tin1-Tin2)

def wymiennik_tout1(e,Tin1,Tin2,m1,cp1,m2,cp2):
        return Tin1-e*m2*cp2*(Tin1-Tin2)/(m1*cp1)


