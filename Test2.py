import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d
import matplotlib.colors
import math
import os
import subprocess
#Variables
rAntennas = 0.05
nAntennas = 30
Max_Time = 5*(10**(-6))
Steps = 10000.0
time = list(np.arange(0., Max_Time, Max_Time/Steps))
c = 2.99792458 * (10**(8))
def simulation(x_value, y_value, filenumber):
	StringFile = str(filenumber)
	XMLfile = open("Project8Phase3_WithRoot_Template.xml","r")
	contents = XMLfile.readlines()
	XMLfile.close()
	if x_value == 0:
		x_value = 0.00000000000000001
	if y_value == 0:
		y_value = 0.00000000000000001
	x = str(x_value)
	y = str(y_value)
	NEWXMLfile = open("Project8Phase3_WithRoot_Template.xml", "w")
	contents[210] = "\t    <x_uniform value_min=\"" + x + "\" value_max=\"" + x + "\"/>\n"
	contents[211] = "\t    <y_uniform value_min=\"" + y + "\" value_min=\"" + y + "\"/>\n"
	contents = "".join(contents)
	NEWXMLfile.write(contents)
	NEWXMLfile.close()
	del contents
	# Run the simulation for a specific position adn change name of .txt file
	os.system("LocustSim config=LocustPhase3Template.json")
	os.system("mv Newfile.txt" + " " + "Newfile" + StringFile + ".txt")
	return
def get_antenna_posy(nAntennas, rAntennas, i):
	antenna_possy = (rAntennas*np.sin(2*np.pi/nAntennas*i))
	return antenna_possy
def get_antenna_posx(nAntennas, rAntennas, i):
	antenna_possx = (rAntennas*np.cos(2*np.pi/nAntennas*i))
	return antenna_possx
def get_antenna_dist(x0, y0, i):
	return np.sqrt((x0-get_antenna_posx(x0,y0,i))**2 + (y0-get_antenna_posy(x0, y0, i))**2)
def get_cyclotron_frequency(filenumber, i):
	StringFile = str(filenumber)
	FileDFT = open("Newfile" + StringFile + ".txt", "r")
	lines = FileDFT.readlines()
	FileDFT.close()
	yaxis = []
	for n in range(int(Steps)):
		List = lines[n].split()
		Voltage = float(List[i])
		yaxis.append(Voltage)	
	DFT = np.fft.fft(yaxis)
	Magnitude = np.absolute(DFT)
	Max_Index = np.argmax(Magnitude)
	frequencies = np.fft.fftfreq(n=10000, d = 5*(10**(-6))/10000)
	if Max_Index > 5000:
		Max_Index = 10000 - Max_Index
	return frequencies[Max_Index]
def get_phase(x0, y0, i):
	d = get_antenna_dist(x0, y0, i)
	return phi = -2. * np.pi * get_cyclotron_frequency(filenumber, i) * d / c




#Functions I haven't finished yet



def shift_by_phase(s, phi,i):
#	return phi shifted signal for ith channel
def shifted_sum(s, x0, y0):
#Did this while really tired so it might be very messed up :)
#	Total = []
#	Sum = 0
#	for n in range (nAntennas):
#		if n == 0:
#			DFT0 = np.fft.fft(shift_by_phase(s, get_phase(x0, y0, n), n ))
#			ABS = np.absolute(DFT0)
#			MAXI = np.argmax(ABS)
#			Angles1 = math.atan2(DFT0.imag[MAXI], DFT0.real[MAX_I})
#			for d in range(len(DFT0)):
#				if d != MAXI:
#					DFT0[d] = 0
#				else:
#					DFT0[d] = DFT0[d]
#			yaxis0 = np.fft.ifft[DFT0]
#			Total = yaxis0
#		else:
#			DFT = np.fft.fft(shift_by_phase(s, get_phase(x0, y0, n), n ))
#			ABS = np.absolute(DFT)
#			MAXI = np.argmax(ABS)
#			Angles2 = math.atan2(DFT2.imag[MAXI], DFT2.real[MAX_I])
#			if Angles2 - Angles1 < 0:
#				shift = (((Angles2 - Angles1) + 2*np.pi)/MAXI/2/(np.pi)*Steps)
#			else:
#				shift = ((Angles2 - Angles1)/ MAXI/2/(np.pi)*Steps)
#			for d in range(len(DFT0)):
#				if d != MAXI:
#					DFT[d] = 0
#				else:
#					DFT[d] = DFT[d]/(np.exp(comples(0,1)*2*np.pi*shift*d/100000))
#			yaxis= np.fft.ifft[DFT]
#		for x in range(len(yaxis)):
#			Total[x] = Total[x] * yaxis[x]
#	Sum = sum(Total)
#	return Sum 

def Get_Grid():
