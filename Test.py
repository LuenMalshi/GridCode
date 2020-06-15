import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d 
import matplotlib.colors 
import math
import os
import subprocess
#Creating the Grid
fig, ax = plt.subplots(1)
phi = np.linspace(0, 2*np.pi, 100)
r = np.sqrt(0.05**2)
x1 = r * np.cos(phi)
x2 = r * np.sin(phi)
ax.plot (x1,x2)
#Variables
Max_Time = 5*(10**(-6)) 
yaxis =[]               #Will be a list of voltages                                     
Steps = 10000.0		#The number of points	
time = list(np.arange (0., Max_Time, Max_Time/Steps))
ShiftValues = []	#The result will be a list of the sum of k for each position
#Starting values for the radius and filenumber
r_value = (0.05)/3.0    
filenumber = 1
#Changing the position of the electron in the .xml file
for iteration_r in range(1,3):				#Outside loop changes the value of r
	r_value = r_value*iteration_r
	for iteration_theta in range(0,8):		#Inside loop changes the value of theta for a given r
		StringFile = str(filenumber)
		XMLfile = open ("Project8Phase3_WithRoot_Template.xml", "r")
		contents = XMLfile.readlines()
		XMLfile.close()
		theta = (2*np.pi/8.0)*iteration_theta
		x_value = r_value*np.cos(theta)
		y_value = r_value*np.sin(theta)
		if x_value == 0:			#When I tried x and y = 0, the simulation wouldn't run
			x_value = 0.00000000000000001
		if y_value == 0:
			y_value = 0.00000000000000001
		#Changing x and y position on the xml file.
		x = str(x_value)
		y = str(y_value)
		NEWXMLfile = open("Project8Phase3_WithRoot_Template.xml", "w")
		contents[210] = "\t    <x_uniform value_min=\"" + x + "\" value_max=\"" + x  + "\"/>\n"
		contents[211] = "\t    <y_uniform value_min=\"" + y + "\" value_max=\"" + y  + "\"/>\n"
		contents = "".join(contents)
		NEWXMLfile.write(contents)
		NEWXMLfile.close()
		del contents
		# Run the simulation for a specific position and change name of .txt file
		os.system ("LocustSim config=LocustPhase3Template.json")
		os.system("mv Newfile.txt"+ " " + "Newfile" + StringFile + ".txt")
		# Compute the sum of k values that each channel has been shifted from the first
		# for a specific electron position
		FileDFT = open ("Newfile" + StringFile + ".txt", "r")
		lines = FileDFT.readlines()
		FileDFT.close()
		#Compute the k values for each channel and sum
		for Channels in range(30):
			for i in range(int(Steps)):
				List = lines[i].split()
				Voltage = float(List [Channels])
				yaxis.append(Voltage)
			DFT = np.fft.fft(yaxis)
			Magnitude = np.absolute(DFT)
			Max_Index = np.argmax(Magnitude)
			if Max_Index > 5000:
				Max_Index = 10000 - Max_Index
			if Channels == 0:
				#atan2 gives angles between -pi and pi
				Angles1 = math.atan2(DFT.imag[Max_Index], DFT.real[Max_Index])
				Total = 0
			else:
				Angles2 = math.atan2(DFT.imag[Max_Index], DFT.real[Max_Index])
				if Angles2 - Angles1 < 0:   
					#Second angle should be positive since you're adding a phase shift?     
					shift = (((Angles2-Angles1)+2*np.pi)/Max_Index/2/(np.pi)*Steps)
				else:
					shift = ((Angles2-Angles1)/Max_Index/2/(np.pi)*Steps)
				Total += shift
			DFT = []
			yaxis = []
		ShiftValues.append(Total) 
		lines = []
		#Plotting the points with the sum of k values
		ax.plot(x_value, y_value, 'bo')
		strv = str(ShiftValues [filenumber-1])
		ax.annotate(strv, xy = (x_value, y_value))
		#Change the file number/name
		filenumber +=1 
print ShiftValues
ax.set_aspect(1)
plt.show()

#Testing that the graphs line up
#I had this code on a different file and it was mixed with the code above
#as I was testing whether the graphs lined up or not. 
#
#filenumber = 1 #The number of the file you want to open
#stringfile = str(filenumber)
#File = open("Newfile" + stringfile + ".txt", "r")
#lines = File.readlines()
#File.close()
#for channels in range (30):
#	for i in range(int(Steps)):
#		list = lines[i].split()
#		voltage = float(list[channels])
#		yaxis.append(voltage)
#	example = np.fft.fft(yaxis)
#	magnitude = np.absolute(example)
#	m = np.argmax(example)
#	if channels == 0:
# 		yaxis0 = yaxis
#	else:
#		plt.plot(time, yaxis0)
#		for d in range(len(example)):	
#			if d != m:
#				example[d] = 0
#			else:
#				example[d] = example[d]/(np.exp(complex(0,1)*2*np.pi*shift*d/10000))
#		yaxis = []
#		yaxis = np.fft.ifft(example)
#		plt.plot(time, yaxis)
#		plt.show() 	
#	example = []
#	yaxis = []
