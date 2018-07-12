#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 15:48:49 2017

ASTR 3830 Project - Cratering Simulation
Adam Cummings
Version 1.1
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import time

tStart = time.time()

# Create a 500x500 km field, broken into 100 x 100 m pixels
# 500 km / 100 m = 5000, so field is 5000 x 5000 pixels
kmFieldLength = 500
meterPixLength = 100
numPixels = int(kmFieldLength * 1000 / meterPixLength) 

# create a square matrix of zeros
field = [ [0] * numPixels for _ in range(numPixels)]
print("Field Created.. T = {}s".format(round(time.time()-tStart, 3)))

def percentCovered():
    # Returns the decimal percent of the field that is covered
    fieldPixels = np.ravel(field) # unravel field into 1D array
    coveredPixels = np.extract(fieldPixels == 1, fieldPixels) # extract all pixels that have value 1
    return float(len(coveredPixels))/len(fieldPixels)

def getDistance(x1, y1, x2, y2):
    # Return distance in pixels between (x1, y1) and (x2, y2) using pythagorean theorem
    #   x1, y1, x2, y2 are matrix indices representing pixel coordinates  
    pixDistance = int(math.sqrt(math.pow((x1 - x2), 2) +  math.pow((y1 - y2), 2)))
    return pixDistance

def impact(x, y, r):
    # Add an impact crater with center (x, y) and radius r on the field
    #   x, y, r are integer # pixels. 20 pix = 1 km. 
    #   x and y must be between 0 and numPixels
    #   Look at an rxr sub-matrix, simulate an impact at its center, then add it to the field 
    for row in range(r*2):
        for col in range(r*2): #check a square area with side length 2 * r
            if(getDistance(r, r, col, row) <= r): #check to see if a particular pixel is within circular radius of center
                if (row + y - r) >= 0 and (col + x - r) >= 0:
                    if (row + y - r) < numPixels and (col + x - r) < numPixels: #check to see if that pixel is within field, ignore if it bleeds off
                        field[row + y - r][col + x - r] = 1

powerLaw = (np.array([np.random.power(3) for _ in range(100000)]) * -1 + 1.1) * 14.3
#This power law generates a nice-looking distribution of radii
plt.hist(powerLaw, 15, normed = True, histtype = 'bar', rwidth = 0.9)
plt.title("Crater radii generated randomly according to power law\n" r'$P = 2x^2$' ", sampled 100k times")
plt.xlabel("Radius [km]")
plt.xlim([1, 17])
plt.figtext(x = 0.19, y = -0.05, s = "min {}, max {}, mean {}, median {}".format(round(np.min(powerLaw), 2), round(np.max(powerLaw), 2), round(np.average(powerLaw), 3), round(np.median(powerLaw), 3)))
plt.show()

tSaturated = 0
fracCovered = 0
craters = []
#time increments by 10 years every step. Maximum of one impact every 10 years.
inputTime = 800000 # This value has a large effect on program runtime
inputTime -= inputTime%3 #Force inputTime to be divisible by 3
plotIncrement = int(inputTime/3)
for t in range(1, inputTime+1):
    if np.random.random() < 0.01: #Only simulate an impact 1% of the time. Avg one impact (with avg diameter 10 km) every 1000 years
        
        # Generate random x and y positions to place the center of asteroid impact
        xPos = int(np.random.random() * numPixels) # np.random.random() generates a random float between 0 and 1
        yPos = int(np.random.random() * numPixels)
        
        # Generate a crater radius with average value of 5.0 km
        kmRad = (np.random.power(3) * - 1 + 1.1) * 14.3
        pixRad = int(kmRad * 1000 / meterPixLength)
        
        impact(xPos, yPos, pixRad)
        craters.append((xPos, yPos, pixRad, t))
        
        #Considered saturated when the field is ~90% filled, which happens around the ~5500th impact
        if len(craters) == 5500:
            tSaturated = t
            fracCovered = percentCovered()
    if(t%plotIncrement == 0):
        print("Plotting graph... T = {}s".format(round(time.time()-tStart, 2)))
        plt.imshow(field, cmap = 'Blues', origin = 'lower')
        plt.title("Craters randomly placed in 500x500 km field\nt = {} million years".format(round(t*10**-5, 1)))
        plt.xlabel("Pixels along x-axis")
        plt.ylabel("Pixels along y-axis")
        plt.figtext(x = 0.26, y = -0.05, s = "Field is {}% covered by {} total impacts".format(round(percentCovered()*100, 1), len(craters)))
        plt.show()


print("Impacts Simulated....... T = {}s".format(round(time.time()-tStart, 2)))
        

#After simulating all impacts, remove the craters that have been covered up
uncoveredCraters = list(reversed(craters))
for parent in reversed(craters):
    for crater in uncoveredCraters:
        if getDistance(crater[0], crater[1], parent[0], parent[1]) < parent[2] and parent[2] > crater[2] and parent[3] > crater[3]:
            #if parent is at a later time than crater, and parent has a larger radius, and crater centerpoint is inside parent's radius...
            uncoveredCraters.remove(crater)
uncoveredCraters.reverse()


print("Found uncovered craters. T = {}s".format(round(time.time()-tStart, 2)))


times = [ele[3] for ele in uncoveredCraters]
fitCoeff = curve_fit(lambda t,a,b: a * np.log(b*t+1),  range(len(times)),  times, p0 = (5*10**5, 3*10**-3))

xrange = np.arange(len(times)*2) #plot regression from 0 to numUncoveredCraters * 2
yrange = fitCoeff[0][0] * np.log(fitCoeff[0][1]* xrange + 1)
plt.plot(xrange, yrange, 'r', alpha = 0.5, lw = 5, ls = '--')
plt.plot(range(len(times)), times, 'r', lw = 3)
plt.hlines(tSaturated, 0, len(times)*2, 'k', linestyles='dashed',alpha = 0.5)
plt.xlabel("Number Uncovered Craters")
plt.ylabel("Time in tens of years")
plt.title("Plot of uncovered craters vs time, with curve fit")
plt.grid()
plt.figtext(x = 0.23, y = -0.07, s = "Time = {}E5 log({}E-3 * numCraters + 1). Saturation\ntime is {} million years, where field is {}% covered".format(round(fitCoeff[0][0]*10**-5, 1), round(fitCoeff[0][1]*10**3, 1), round(tSaturated*10**-5, 1), round(fracCovered*100, 1)))
plt.show()
            

print("Curve fit and plotted time vs uncovered craters. T = {}s".format(round(time.time()-tStart, 3)))


#=========================================================================================
#=========================================================================================

field = [ [0] * numPixels for _ in range(numPixels)]
print("Field Reset.... T = {}s".format(round(time.time()-tStart, 3)))

#this is the second power law to use for changed assumptions
powerLaw2 = (np.array([np.random.power(8) for _ in range(100000)]) * -1 + 1.05) * 31.1
plt.hist(powerLaw2, 25, normed = True, histtype = 'bar', rwidth = 0.9)
plt.title("Crater radii generated randomly according to power law\n" r'$P = 7x^7$' ", sampled 100k times")
plt.xlabel("Radius [km]")
plt.xlim([1, 17])
plt.figtext(x = 0.19, y = -0.05, s = "min {}, max {}, mean {}, median {}".format(round(np.min(powerLaw2), 2), round(np.max(powerLaw2), 2), round(np.average(powerLaw2), 3), round(np.median(powerLaw2), 3)))
plt.show()

craters = []
#time increments by 10 years every step. Maximum of one impact every 10 years.
for t in range(1, inputTime+1):
    if np.random.random() < 0.01: #Only simulate an impact 1% of the time. Avg one impact (with avg diameter 10 km) every 1000 years
        
        # Generate random x and y positions to place the center of asteroid impact
        xPos = int(np.random.random() * numPixels) # np.random.random() generates a random float between 0 and 1
        yPos = int(np.random.random() * numPixels)
        
        # Generate a crater radius with average value of 5.0 km
        kmRad = (np.random.power(8) * - 1 + 1.05) * 31.1
        pixRad = int(kmRad * 1000 / meterPixLength)
        
        impact(xPos, yPos, pixRad)
        craters.append((xPos, yPos, pixRad, t))
        
        #Considered saturated when the field is ~90% filled, which happens around the ~5500th impact
        if len(craters) == 5500:
            tSaturated = t
            fracCovered = percentCovered()
    if(t%plotIncrement == 0):
        print("Plotting graph... T = {}s".format(round(time.time()-tStart, 2)))
        plt.imshow(field, cmap = 'Blues', origin = 'lower')
        plt.title("Craters randomly placed in 500x500 km field\nt = {} million years".format(round(t*10**-5, 1)))
        plt.xlabel("Pixels along x-axis")
        plt.ylabel("Pixels along y-axis")
        plt.figtext(x = 0.26, y = -0.05, s = "Field is {}% covered by {} total impacts".format(round(percentCovered()*100, 1), len(craters)))
        plt.show()


print("Impacts Simulated....... T = {}s".format(round(time.time()-tStart, 2)))
        

uncoveredCraters = list(reversed(craters))
for parent in reversed(craters):
    for crater in uncoveredCraters:
        if getDistance(crater[0], crater[1], parent[0], parent[1]) < parent[2] and parent[2] > crater[2] and parent[3] > crater[3]:
            uncoveredCraters.remove(crater)
uncoveredCraters.reverse()


print("Found uncovered craters. T = {}s".format(round(time.time()-tStart, 2)))


times = [ele[3] for ele in uncoveredCraters]
fitCoeff = curve_fit(lambda t,a,b: a * np.log(b*t+1),  range(len(times)),  times, p0 = (5*10**5, 3*10**-3))

xrange = np.arange(len(times)*2) #plot regression from 0 to numUncoveredCraters * 2
yrange = fitCoeff[0][0] * np.log(fitCoeff[0][1]* xrange + 1)
plt.plot(xrange, yrange, 'r', alpha = 0.5, lw = 5, ls = '--')
plt.plot(range(len(times)), times, 'r', lw = 3)
plt.hlines(tSaturated, 0, len(times)*2, 'k', linestyles='dashed',alpha = 0.5)
plt.xlabel("Number Uncovered Craters")
plt.ylabel("Time in tens of years")
plt.title("Plot of uncovered craters vs time, with curve fit")
plt.grid()
plt.figtext(x = 0.23, y = -0.07, s = "Time = {}E5 log({}E-3 * numCraters + 1). Saturation\ntime is {} million years, where field is {}% covered".format(round(fitCoeff[0][0]*10**-5, 1), round(fitCoeff[0][1]*10**3, 1), round(tSaturated*10**-5, 1), round(fracCovered*100, 1)))
plt.show()
            

print("Curve fit and plotted time vs uncovered craters. T = {}s".format(round(time.time()-tStart, 3)))
