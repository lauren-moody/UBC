# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 19:28:58 2020

@author: lauren
"""

import numpy as np
import matplotlib.pyplot as plt

#mean mosquito development rate parameter
rh = 0.15460
ha = 33255.57
hh = 50543.49
th = 301.67

#probability mean and standard deviation
meanp = 0.8
stdevp = 0.1

#human recovery rate mean and standard deviation
hmeanrec = 0.125
hstdevrec = 0.05 
    
days = 210 #30 weeks
#-----------------------------------
#Set up arrays for storing the indiviudals in each human population
q=np.zeros(days*1)
X1=np.zeros(days*1)
X2=np.zeros(days*1)
X3=np.zeros(days*1)
X4=np.zeros(days*1)

#set up parameters with initial values for model ---------------       
X2[0] = 140 #infected
X3[0] = 198 #infectious
X4[0] = 0 #removed/recovered
N = 8275165 #population of veracruz
X1[0] = N-X2[0]-X3[0] #susceptible

#Set up arrays for storing the individuals in each mosquito population
q=np.zeros(days*1)
Z1=np.zeros(days*1)
Z2=np.zeros(days*1)
Z3=np.zeros(days*1)

Z2[0] = 1400 #infected
Z3[0] = 1980 #infectious
# X4[0] = 0 #removed/recovered
P = 10000000 #mosquito population
Z1[0] = P - Z2[0] - Z3[0] #susceptible

emerge = 16.05 #mosquito emergence rate - fitting parameter

#parameters/lists for multiple runs
run = 0
totalruns = 300
increase = 0
cases = []

#connect to spreadsheet
from xlwt import Workbook 

wb=Workbook()
sheet1 = wb.add_sheet('Sheet 1')

#main function
def tempdengue():
    run=0
    while run<=totalruns: #set up number of runs
        for i in range(1,days): #this loop for each day
            #temperature array based on real temperature data from Veracruz
            temp = np.array([29.25,28.5,29.125,28.875,29.25,29.625,29.25,27.75,28.5,26.75,26,26,27.125,28.25,29.75,29.75,30.5,30.375,30.375,29.125,28.125,28.25,26.5,27.625,28.25,28.875,28.875,29.625,29.875,29.25,30.125,29.25,30,29.25,28.125,25.125,25.75,25.25,26.125,27.625,28.125,29,28.625,28.5,27.625,27.5,28.875,27.5,28.25,28.5,29.5,28.875,28.625,29.125,29.5,27,27.625,28.125,27.875,27.625,27.875,27.5,28.25,26.375,26.75,27.375,27.625,27.75,28.375,28.25,28.375,29.125,28.625,28.75,28.25,28.625,29,30.25,29.125,29.25,28.625,29,29.625,28.5,29,29.125,28.875,28.75,26.75,27.5,28,29,29.625,29.75,29.75,29.75,29.5,28.875,28.875,28.5,29.625,28.375,28.75,29.25,28.75,29.625,28.375,27.5,27.5,27.875,28.375,27,25.875,27,27.625,27.875,27.5,27.5,28.25,28.75,28.625,28,27.625,28.375,27.875,28.625,28.75,28.125,27.5,26.625,28.25,28.25,28.125,27.5,27.625,26.875,27.375,26.625,27,27,27.25,28,25,24.75,26.375,26.875,24.75,25.625,26.625,28.125,27.875,28.125,27.125,25.75,26.25,25.5,23.625,23.25,24.625,26.375,27.125,22.5,20.875,22.375,23.875,24.25,25,25.125,25,23.625,23.75,24.125,24.25,21.875,19.5,21.75,21.875,21.5,21.875,22.625,23.5,22.75,22.25,23,24.125,23.5,23.25,24.875,25.625,25,25.75,25.75,25.66666667,24.125,22.875,23.625,23.25,23.875,24.125,22.75,22.75,24.125,22.5,21.625,21.5,22,23,23.625,21.875,17.625,19.83333333,18.375,20.875])
            t = temp[i] #daily temperature is average temperature for given day determined from data
            t = t+increase #for increasing temperatures under emission scenarios

            #bite rate
            b = (0.0043*t) - 0.0645
            bite = b
    #       
            #intrinsic incubation period
            incubatehuman = (np.exp(np.exp((0.56))) * (np.exp(1/(2*13.7))))  # incubation rate / day
            
            #probability mosquito to human
            bmh = np.zeros((1,1))
            for y in range (bmh.shape[0]):
                for z in range (bmh.shape[1]):
                    bmh[y,z] = meanp + np.random.normal()*stdevp
                    if bmh[y,z]<0:
                        bmh[y,z] = 0
                    elif bmh[y,z]>1:
                        bmh[y,z] = 1
                    else:
                        bmh[y,z] = bmh[y,z]

            #probability human to mosquito
            bhm = ((0.001044*(t))*((41-t)**0.5))*5.2
            bitehm = bhm
            
            #mosquito mortality rate
            m = 0.8692 - 0.1590*t + 0.01095*t**2 - 0.0003408*t**3 + 0.00000405*t**4 + 0.1 
            mo = m
            
            #extrinsic incubation period
            E = (np.exp(np.exp((2.9-(0.07765*t))))) * (np.exp(1/(2*4.9)))
            eip = E
            
            #mosquito development rate
            d = ((rh*((t+273.15)/298.15)*np.exp((ha/1.987)*((1/298.15)-(1/(t+273.15)))))/(1+np.exp((hh/1.987)*((1/th)-(1/(t+273.15))))))*1.5
            dev = d
            
            #human recovery rate
            hrecover = np.zeros((1,1))
            for f in range (hrecover.shape[0]):
                for g in range (hrecover.shape[1]):
                    hrecover[f,g] = hmeanrec + np.random.normal()*hstdevrec
                    if hrecover[f,g] <0:
                        hrecover[f,g] = 0
                    elif hrecover[f,g]>1:
                        hrecover[f,g] = 1
                    else:
                        hrecover[f,g] = hrecover[f,g]
            
            #reproduction number
            R = ((bite**2)*(bmh[y,z])*(bitehm)*(np.exp(mo*(1/eip)))*dev)/(hrecover[f,g]*(mo**3))
        
            #days     
            j=i-1
                
            #flow rates between human compartments
            F12 = ((Z3[j]*bite*bmh[y,z])*X1[j])/N
            F23 = ((1/incubatehuman)*X2[j])
            F34 = (hrecover[f,g]*X3[j])
            
            #calculate derivatives
            dX1 = -F12
            dX2 = F12-F23
            dX3 = F23-F34
            dX4 = F34
        
        	#increment compartments
            X1[i]=X1[j]+dX1
            X2[i]=X2[j]+dX2
            X3[i]=X3[j]+dX3
            X4[i]=X4[j]+dX4
        
        	# save time
            q[i]=i
            
            #flow rates between mosquito compartments
            G01 = (dev*emerge*P)
            G12 = ((X3[j]*bite*bitehm)*Z1[j])/P 
            G23 = ((1/eip)*Z2[j])
            G10 = mo*Z1[j]
            G20 = mo*Z2[j]
            G30 = mo*Z3[j]
            
            #calculate derivatives
            dZ1 = G01-G12-G10
            dZ2 = G12-G23-G20
            dZ3 = G23-G30
        
        	#increment compartments
            Z1[i]=Z1[j]+dZ1
            Z2[i]=Z2[j]+dZ2
            Z3[i]=Z3[j]+dZ3
            
        	# save time
            q[i]=i
            
            #print reproduction number
            print(R)
            
            #save cases per day
            cases.append(X3[i]) 
            casenumbers=X3[i]
            
            #export number of cases per day to spreadsheet
            print(casenumbers)
            print(run)
            sheet1.write(i,run,casenumbers)
            wb.save('printvaluesveracruz.xls')

        run=run+1 #next run

tempdengue() #run main function
    
#additional features:

print("The fraction of the human population that is susceptible is " + str((X1[days-1]/N)*100) + "%.")
print("The fraction of the human population that is infected is " + str((X2[days-1]/N)*100) + "%.")
print("The fraction of the human population that is infectious is " + str((X3[days-1]/N)*100) + "%.")
print("The fraction of the human population that is recovered or dead is " + str((X4[days-1]/N)*100) + "%.")
print("The number of cases is " + str(X2[days]))

#plot the data
q=q #convert time vector to months
plt.figure(1) #initialize plot
plt.title('Spread of Disease in Humans')
# if  producers:plt.plot(t,X1,'g-',label='Prod') #plot
plt.plot(q,X1/N,'r-',label='Susceptible')
plt.plot(q,X2/N,'b-',label='Infected')
plt.plot(q,X3/N,'c-',label='Infectious')
plt.plot(q,X4/N,'k-',label='Recovered')

plt.legend(loc=1)
plt.xlabel('Time (days)')
plt.ylabel('Fraction of Population')
plt.show()