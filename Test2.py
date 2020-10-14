import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d
import matplotlib.colors as mcolors
import math
import os
import subprocess

#Variables
rAntennas = 0.05
nAntennas = 30
dt = 5e-10
nSteps = 10
time = np.arange(0, nSteps*dt, dt)

np.random.seed(4)

c = 2.99792458e8

def simulation(x_value, y_value, filenumber):
	StringFile = str(filenumber)
	XMLfile = open("Project8Phase3_WithRoot_Template.xml","r")
	contents = XMLfile.readlines()
	XMLfile.close()
	x = str(x_value)
	y = str(y_value)
	NEWXMLfile = open("Project8Phase3_WithRoot_Template.xml", "w")
	contents[210] = "\t    <x_fix value=\"" + x + "\"/>\n"
	contents[211] = "\t    <y_fix value=\"" + y + "\"/>\n"
	contents = "".join(contents)
	NEWXMLfile.write(contents)
	NEWXMLfile.close()
	del contents
	# Run the simulation for a specific position and change name of .txt file
	os.system("LocustSim config=LocustPhase3Template.json")
	os.system("mv Newfile.txt" + " " + "Newfile" + StringFile + ".txt")
	return


def fake_data(x0, y0, t):
    s = np.zeros(shape=(nAntennas,nSteps),dtype=np.complex_)
    for i in range(nAntennas):
        r = get_antenna_dist(x0,y0,i)
        #s[i] = np.exp(-1j * 2. * np.pi * get_cyclotron_frequency() * ( t - r / c))
        s[i] = 1. / r * np.exp(-1j * 2. * np.pi * get_cyclotron_frequency() * ( - r / c))
    return s


def get_antenna_pos(nAntennas, rAntennas, i):
	phase = 2*np.pi / nAntennas * i
	return [rAntennas * np.cos(phase), rAntennas * np.sin(phase)]


def get_antenna_dist(x0, y0, i):
    r = get_antenna_pos(nAntennas,rAntennas,i)
    return np.sqrt((x0-r[0])**2 + (y0-r[1])**2)

def get_cyclotron_frequency():
    return 25.6e9

def get_phase(x0, y0, i):
	d = get_antenna_dist(x0, y0, i)
	return -2. * np.pi * get_cyclotron_frequency() * d / c


def add_noise(s, sigma):
    return s+np.random.normal(scale=sigma, size=s.shape)


#return phi shifted signal for ith channel 
def shift_by_phase(s, phi,i):
    return s * np.exp( 1j * phi)


def shifted_sum(s, x0, y0, weighted=False):
    total_signal = 0
    normalization = 0 #s^2
    for i in range(nAntennas):
        phi_shift = get_phase(x0, y0, i)
        d = get_antenna_dist(x0, y0, i)
        shifted_sig = shift_by_phase(s[i], phi_shift, i) #/ d
        normalization += 1. / d**2
        total_signal +=  shifted_sig

    #total_signal = total_signal / np.sqrt(normalization)

    return total_signal

def get_grid(a,b,c,d):
    x = np.linspace(a, b, 200)
    y = np.linspace(c, d, 200)

    X,Y = np.meshgrid(x,y,indexing="ij")

    return X, Y

#Checking around the x axis first
if __name__ == "__main__":
    #Generate data
	y0 = 0
#	x0 = -0.015
	distlist = []
	fulllist = []
	for x0 in np.arange (-0.03, 0.03, 0.005):
		sData = fake_data(x0,y0,time) 
		kNoise = 50
		distlist.append(np.abs(x0))
#	sData = add_noise(sData, kNoise) + 1j * add_noise(sData, kNoise)
#	print(sData)
#	color map
#		colors1 = plt.cm.binary(np.linspace(0., 1, 1))
#		colors2 = plt.cm.viridis(np.linspace(0, 1, 255))
#		colors = np.vstack((colors1, colors2))
#		mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
		X,Y = get_grid(-rAntennas, rAntennas, -rAntennas, rAntennas)
		Z = 0 * X	    
#	    print X, Y
	    #compute sum for each grid point
		for i in range(len(X)):
 		       for j in range(len(X[0])):
        		    if X[i][j]**2 + Y[i][j]**2 < 0.99*rAntennas**2:
       	     	 		   Z[i][j] = np.abs(sum(shifted_sum(sData , X[i][j], Y[i][j])))
	    #1
	    #Creating a list D with the Z values only at x axis
	        n = 99
	    	D =[]
	    	for i in range(len(X)):
			D.append(Z[i][n])
#			print X[i][n], Y[i][n]                   #Check X and Y values
#		print D   					 #Z values at X axis
		ml = np.where(max(D)) 				 #Index of Max D  value
	    	half = max(D)/2.0				 #Half Max Value
	    	HD1 = []			
	    	HD2 = []
	    	for i in range(len(X)):
			if X[i][n] < x0:
				HD1.append(np.abs(D[i]-half))    #Creating a list with the difference of Z values from Half max to the left.
				HD2.append(100000)
			else:
				HD1.append(100000)
				HD2.append(np.abs(D[i]-half))	 #Creating a list with the difference of Z values from Half max to the right
		X1, X2, Y1, Y2, D1, D2, n1, n2 = X, X, Y, Y, D, D, n, n
		xleft = np.abs(X1[np.argmin(HD1)][n1] - x0)
		xright = np.abs(X2[np.argmin(HD2)][n2]- x0)
	    #2
            #Finding only the closer value. I was getting 2 values corresponding to half max. 
	    #This part of the code chooses the x value closer to the electron location. 
		if  xleft > 0.01:
			HD1 = []
			for i in range(len(X1)):
				if X1[i][n1] < x0 :
					if np.abs(X1[i][n1] - x0) > 0.01:
						HD1.append(1000000)
					else:		
						HD1.append(np.abs(D1[i]-half))
				else:
					HD1.append(10000)
#		print HD1
		if  xright > 0.01:
			HD2 = []
			for i in range(len(X2)):
				if X2[i][n2] > x0 :
					if np.abs(X2[i][n2] - x0) > 0.01:
						HD2.append(1000000)
					else:		
						HD2.append(np.abs(D2[i]-half))
				else:
					HD2.append(10000)
#This part of the code deals with precision. For an electron position, The code above determines Value of half max and X coordinate where that half max happens. In this part, if the potential corresponding to this X value is more than 100 eV away from actual half max, the code zooms in to find a more accurate X value, for both left and right sides. 

		if HD2[np.argmin(HD2)] > 100:
			X2,Y2 = get_grid(x0 - 0.005, X2[np.argmin(HD2)][n2] + 0.005, -0.005, 0.005)
			
			HD2= []
			Z2 = 0 * X2
			for i in range(len(X2)):
				for j in range(len(X2[0])):
					if X2[i][j]**2 + Y2[i][j]**2 < 0.99*rAntennas**2:
						Z2[i][j] = np.abs(sum(shifted_sum(sData, X2[i][j], Y2[i][j])))	
			ylist = Y2[0]
			n2 = np.argmin(np.abs(ylist))
#			print Y2[0]
			D2 = []
	    		for i in range(len(X2)):
				D2.append(np.abs(Z2[i][n2]))
			half = max(D2)/2.0
			print max(D2)		
			for i in range(len(X2)):
				if X2[i][n2] > x0 :
					if np.abs(X2[i][n2] - x0) > 0.01:
						HD2.append(1000000)
					else:		
						HD2.append(np.abs(D2[i]-half))
				else:
					HD2.append(10000)
		if HD1[np.argmin(HD1)] > 100:
			X1,Y1 = get_grid( X1[np.argmin(HD1)][n2] - 0.005,x0 + 0.005, -0.005, 0.005)
			
			HD1= []
			Z1 = 0 * X1
			for i in range(len(X1)):
				for j in range(len(X1[0])):
					if X1[i][j]**2 + Y1[i][j]**2 < 0.99*rAntennas**2:
						Z1[i][j] = np.abs(sum(shifted_sum(sData, X1[i][j], Y1[i][j])))	
			ylist = Y1[0]
			n1 = np.argmin(np.abs(ylist))
#			print Y2[0]
			D1 = []
	    		for i in range(len(X1)):
				D1.append(np.abs(Z1[i][n1]))
			half = max(D1)/2.0
#			print max(D1)		
			for i in range(len(X1)):
				if X1[i][n1] < x0 :
					if np.abs(X1[i][n1] - x0) > 0.01:
						HD1.append(1000000)
					else:		
						HD1.append(np.abs(D1[i]-half))
				else:
					HD1.append(1000000)

	#Printing values just to check. 
#		print np.argmin(HD1), np.argmin(HD2), HD1[np.argmin(HD1)], HD2[np.argmin(HD2)], X1[np.argmin(HD1)][n1], X2[np.argmin(HD2)][n2]
#		print max(D), D1[np.argmin(HD1)], D2[np.argmin(HD2)]
		fulllist.append(np.abs(X1[np.argmin(HD1)][n1]-X2[np.argmin(HD2)][n2]))
	print distlist, fulllist
	#Plotting
	plt.plot(distlist, fulllist)
	plt.show()		
#       fig = plt.figure(figsize=[7,7])
#       ax = fig.gca()

#	plotting
#       plt.plot(x0,y0,".",color="r")
#       s = plt.pcolor(X, Y, Z, cmap=mymap)
#	plt.xlabel("X (m)")
#	plt.ylabel("Y (m)")
#	cbar = fig.colorbar(cs)
#        plt.show()


