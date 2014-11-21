'''  
  Compute the correlation of epsilon at two different periods for shallow
  crustal earthquakes in the PEER NGA database
 
  Reference: Baker, J.W. and Jayaram, N. (2008), "Correlation of spectral 
    acceleration values from NGA ground motion models," Earthquake 
    Spectra, 24 (1), 299-317. 
 
  INPUT
  T1, T2    = the two periods of interest. The periods may be equal, and 
              there is no restriction on which is larger.
 
  OUTPUT
  rho       = predicted correlation coefficient
 
'''
import numpy as np

def baker_jayaram_correlation(T1, T2):
    T_min = min(T1, T2)
    T_max = max(T1, T2)
    
    C1 = (1-np.cos(np.pi/2 - np.log(T_max/max(T_min, 0.109)) * 0.366 ))
    if T_max < 0.2:
        C2 = 1 - 0.105*(1 - 1./(1+np.exp(100*T_max-5)))*(T_max-T_min)/(T_max-0.0099)
    
    if T_max < 0.109:
        C3 = C2
    else:
        C3 = C1

    C4 = C1 + 0.5 * (np.sqrt(C3) - C3) * (1 + np.cos(np.pi*(T_min)/(0.109)))

    if T_max <= 0.109:
        rho = C2
    elif T_min > 0.109:
        rho = C1
    elif T_max < 0.2:
        rho = min(C2, C4)
    else:
        rho = C4
    
    return rho



print(baker_jayaram_correlation(1,0.3))