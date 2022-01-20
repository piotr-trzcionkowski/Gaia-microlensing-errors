#!/usr/bin/env python3.6

import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

	
def symulacja(t, dt, xa, ya, pa, vda, xb, yb, pb, vdb, pc, vdc, fs, f):
	if (pc > 0):
		xc = -pc
		yc = 0.0
	for i in range (800):
		ia = 45 #[deg] kat inklinacji
		xa1 = xa + pa * np.cos(i * 2 * np.pi / 730)
		ya1 = ya + pa * np.sin(i * 2 * np.pi / 730) * np.cos(ia * np.pi / 180)
		xb1 = xb + pb * np.cos(i * 2 * np.pi / 730)
		yb1 = yb + pb * np.sin(i * 2 * np.pi / 730) * np.cos(ia * np.pi / 180)
		xs = fs * xa1 + xb1 * (1 - fs)
		ys = fs * ya1 + yb1 * (1 - fs)
		if (pc>0):
			xc1=xc+pc*np.cos(i*2*np.pi/730)
			yc1=yc+pc*np.sin(i*2*np.pi/730)*np.cos(ia*np.pi/180)
			xc+=dt*vdc[0]
			yc+=dt*vdc[1]
			f.write("%i %f %f %f %f %f %f %f %f %f\n" % (t, xa1, ya1, xb1, yb1, xs, ys, xc1, yc1, fs))
		else:	
			f.write("%i %f %f %f %f %f %f %f\n" % (t, xa1, ya1, xb1, yb1, xs, ys, fs))
		xa+=dt*vda[0]
		ya+=dt*vda[1]
		xb+=dt*vdb[0]
		yb+=dt*vdb[1]
		t+=dt

def generujplik(t, xa, ya, pa, vda, xb, yb, pb, vdb, pc, vdc, fs):
	#Obliczanie pozycji gwiazd i wypadkowej na niebie w czasie od 0 do t i od 0 do -t (w tej kolejnosci)
	f = open("gwiazdapresort.dat", "w")
	#f.write("#T RA1 Dec1 RA2 Dec2 RA_wyp Dec_wyp\n")
	symulacja(t, +1, xa, ya, pa, vda, xb, yb, pb, vdb, pc, vdc, fs, f)
	symulacja(t, -1, xa, ya, pa, vda, xb, yb, pb, vdb, pc, vdc, fs, f)
	f.close()

	# Sortowanie powyzszych obliczen po czasie, zeby wykresy nie mialy brzydkiego przeskoku
	fn = 'gwiazdapresort.dat'
	sorted_fn = 'gwiazda.dat'

	with open(fn,'r') as first_file:
		rows = first_file.readlines()
		sorted_rows = sorted(rows, key=lambda x: int(x.split()[0]), reverse=False)
		with open(sorted_fn, 'w') as second_file:
			second_file.write("#T RA1 Dec1 RA2 Dec2 RA_wyp Dec_wyp fs\n")
			for row in sorted_rows:
				second_file.write(row)

def obliczparalakse():
	# Obliczanie paralaksy dla podanej F
	#paralaksa = 'gwiazda.dat'
	pozycje=np.loadtxt('gwiazda.dat', comments='#', delimiter=' ', unpack=False)
	ruchX=(pozycje[730][5]-pozycje[0][5])/730
	ruchY=(pozycje[730][6]-pozycje[0][6])/730
	paralaksaG=0
	pozycjepar=np.zeros((730,2))
	pozycjepar[0][0]=pozycje[0][5]
	pozycjepar[0][1]=pozycje[0][6]
	odleglosci=np.zeros((365,2))
	odleglosci[0][0]=pozycjepar[365][0]-pozycjepar[0][0]
	odleglosci[0][1]=pozycjepar[365][1]-pozycjepar[0][1]

	for i in range(730):
		pozycjepar[i][0]=pozycje[i][5]-i*ruchX
		pozycjepar[i][1]=pozycje[i][6]-i*ruchY

	for i in range(365):
		odleglosci[i][0]=(pozycjepar[i+365][0]-pozycjepar[i][0])/2
		odleglosci[i][1]=(pozycjepar[i+365][1]-pozycjepar[i][1])/2
		if paralaksaG < (odleglosci[i][0] ** 2 + odleglosci[i][1] ** 2) ** 0.5:
			paralaksaG = (odleglosci[i][0] ** 2 + odleglosci[i][1] ** 2) ** 0.5
	return paralaksaG, ruchX, ruchY

def rysujwykres(pc, FS, va, vb):
	#Rysowanie wykresu
	plt.subplots(figsize=(8,8))
	if (pc>0):
		Tp, RA1p, Dec1p, RA2p, Dec2p, RA_wypp, Dec_wypp, RA_c, Dec_c, fs = np.loadtxt('gwiazda.dat', delimiter=' ', unpack=True)
		plt.plot(RA1p, Dec1p, label='Source')
		plt.plot(RA2p, Dec2p, label='Blend')
		plt.plot(RA_wypp, Dec_wypp, label='Detected star', linewidth=3)
		plt.plot(RA_c, Dec_c, label='simulated')
	else:
		Tp, RA1p, Dec1p, RA2p, Dec2p, RA_wypp, Dec_wypp, fs = np.loadtxt('gwiazda.dat', delimiter=' ', unpack=True)
		plt.plot(RA1p, Dec1p, label='Source')
		plt.arrow(RA1p[0] + 0.2, Dec1p[0] + 0.2, va[0]/5, va[1]/5, width=0.05, color='blue')
		plt.plot(RA2p, Dec2p, label='Blend')
		plt.arrow(RA2p[0] - 0.1, Dec2p[0] + 0.3, vb[0]/5, vb[1]/5, width=0.05, color='orange')
		plt.plot(RA_wypp, Dec_wypp, label='Detected star', linewidth=3)
		plt.arrow(RA_wypp[0] + 0.5, Dec_wypp[0], (FS * va[0] + vb[0] * (1 - FS))/5, (FS * va[1] + vb[1] * (1 - FS))/5, width=0.05, color='green')
	plt.xlabel('RA [mas]')
	plt.ylabel('Dec [mas]')
	# plt.title('Zmiana polozenia gwiazd rzeczywistch i obserwowanej wzgledem polozenia soczewki')
	plt.legend()
	plt.text(2, -4, 'fs=%.2f' %(FS))
	plt.savefig('fs=%.2f.png' %(FS))
	plt.show()

def zapisywaniedanych():
	#pytanie o zapisanie danych
	pyt = 'n'
	# pyt=input('Czy chcesz zapisac symulowane pozycje gwiazd i obrazu?(y/n)')
	if pyt=='y':
		dane = open('fs=%.2f.dat' %(fs), 'w')
		print('Plik zostal osobno zapisany')
	else:
		print('Plik nie zostal osobno zapisany')


print("Co chcesz wykonac?\n\
	(1) Narysowac pozycje soczewki i blendy\n\
	(2) Obliczyc zaleznosc paralaksy\n\
	(3) Wyznaczyc ruch zblendowanego ukladu\n\
	(4) Wyznaczyc zaleznosc paralaks od fs\n\
	")
pyt=input()
if pyt == '1':
	
	#Zdefiniowanie parametrow poczatkowych
	t = 0 							#czas t=0
	pa=1./8. 						#[mas] - paralaksa zrodla
	xa=-pa 							#polozenie zrodla w momencie soczewkowania skorygowane o paralakse
	ya=0.0
	va=[-2.0, 2.0] 					#[mas/yr] ruch wlasny zrodla
	vda=[va[0]/730.0, va[1]/730.0] 	#[mas/d]

	pb=0.5							#[mas] - paralaksa blendy
	xb=-pb							#polozenie blendy w momencie soczewkowania skorygowane o jej paralakse
	yb=0.0
	vb=[4.0, 4.0] 					#[mas/yr] ruch wlasny blendy 
	vdb=[vb[0]/730.0, vb[1]/730.0] 	#[mas/d]

	fs=float(input("Podaj blending parameter fs = "))
	FS=fs	
	if (fs >= 1) or (fs < 0):
		print("Podałeś niefizyczny parametr fs. Następnym razem podaj 0 < fs < 1.")
	else:
		generujplik(t, xa, ya, pa, vda, xb, yb, pb, vdb, 0, 0, fs)
		rysujwykres(0, FS, va, vb)
		zapisywaniedanych()

elif pyt == '2':
	
	#Zdefiniowanie parametrow poczatkowych
	t = 0 							#czas t=0
	pa=1./8. 						#[mas] - paralaksa zrodla
	xa=-pa 							#polozenie zrodla w momencie soczewkowania skorygowane o paralakse
	ya=0.0
	va=[-2.0, 2.0] 					#[mas/yr] ruch wlasny zrodla
	vda=[va[0]/730.0, va[1]/730.0] 	#[mas/d]


	file= open('funkcja paralaks.dat', 'w')
	file.write('#x\tpi_S\tpi_B\tpi_G\tpi_G/pi_S\tfs\n') #x=pi_S/pi_L
	for j in np.linspace(2, 10, num=10):

		pb=j*pa
		xb=-pb							#polozenie blendy w momencie soczewkowania skorygowane o jej paralakse
		yb=0.0
		vb=[4.0, 4.0] 					#[mas/yr] ruch wlasny blendy 
		vdb=[vb[0]/730.0, vb[1]/730.0] 	#[mas/d]


		for i in range(20):
	#		fs = i / 10
	#		F = fs / (1 - fs)
			F = ((i + 1) / 10)
			fs = F / (F + 1)
			generujplik(t, xa, ya, pa, vda, xb, yb, pb, vdb, 0, 0,fs)
			file.write('%f\t%f\t%f\t%f\t%f\t%f\n' % (pa/pb, pa, pb, obliczparalakse()[0], obliczparalakse()[0]/pa, fs))
	file.close()
	
	#Rysowanie wykresu paralaks
	x, pi_S, pi_L, pi_G, GS, fs = np.loadtxt('funkcja paralaks.dat', delimiter='\t', unpack=True)
	
	#fig, axs = plt.subplots()
	#oints = np.array([fs, GS]).T.reshape(-1,1,2)
	#segments= np.concatenate([points[:-1], points[1:]], axis=1)
	#norm = plt.Normalize(x.min(), x.max())
	#lc = LineCollection(segments, cmap='nipy_spectral', norm=norm)
	#lc.set_array(F)
	#lc.set_linewidth(2)
	#line = axs.add_collection(lc)
	#fig.colorbar(line, ax=axs)
	#axs.set_xlim(fs.min(), fs.max())
	#axs.set_ylim(GS.min(), GS.max())
	#plt.set_color_cycle
	plt.figure(figsize=(14,12))
	plt.plot(fs, GS, '.', c='red', label='pi_G/pi_S (fs)')
	plt.plot(pi_S/pi_L, GS, '.', c="blue", label='pi_G/pi_S (pi_S/pi_L)')
	plt.xlabel('')
	plt.ylabel('pi_G/pi_S')
	plt.title('Zaleznosc pi_G/pi_S w funkcjach (pi_S/pi_L) i (fs) dla różnych stosunków jasności')
	plt.savefig('GS(SB).png')
	plt.legend()
	plt.tight_layout()
	plt.show()

elif pyt == '3':
	
	#Zdefiniowanie parametrow poczatkowych
	t = 0 							#czas t=0
	pa=1./8. 						#[mas] - paralaksa zrodla
	xa=-pa 							#polozenie zrodla w momencie soczewkowania skorygowane o paralakse
	ya=0.0
	va=[-2.0, -4.0] 				#[mas/yr] ruch wlasny zrodla
	vda=[va[0]/730.0, va[1]/730.0] 	#[mas/d]

	pb=1./4.						#[mas] - paralaksa blendy
	xb=-pb							#polozenie blendy w momencie soczewkowania skorygowane o jej paralakse
	yb=0.0
	vb=[-5.0, -5.0] 				#[mas/yr] ruch wlasny blendy 
	vdb=[vb[0]/730.0, vb[1]/730.0] 	#[mas/d]

	fs=0.5
	generujplik(t, xa, ya, pa, vda, xb, yb, pb, vdb, 0, 0, fs)
	vdc = [0 , 0]
	pc, vdc[0], vdc[1] = (obliczparalakse())
	generujplik(t, xa, ya, pa, vda, xb, yb, pb, vdb, pc, vdc, fs)
	rysujwykres(pc)
	zapisywaniedanych()

elif pyt == '4':
	
	file= open('piB(piG, piS, fs).dat', 'w')
	file.write('#    piS      piB       fs      piG    piB_f    piB_e       vSx      vSy      vBx      vBy      vGx      vGy     vBfx      vBfy    vBx_e    vBy_e\n')

	
	#Zdefiniowanie parametrow poczatkowych
	t = 0 							#czas t=0
	pa=1./8. 						#[mas] - paralaksa zrodla
	xa=-pa 							#polozenie zrodla w momencie soczewkowania skorygowane o paralakse
	ya=0.0
	va=[-2.0, 2.0] 					#[mas/yr] ruch wlasny zrodla
	vda=[va[0]/730.0, va[1]/730.0] 	#[mas/d]

	for j in np.linspace(1.25, 8, num=5):

		pb=j*pa
		xb=-pb							#polozenie blendy w momencie soczewkowania skorygowane o jej paralakse
		yb=0.0
		vb=[4.0, 4.0] 					#[mas/yr] ruch wlasny blendy 
		vdb=[vb[0]/730.0, vb[1]/730.0] 	#[mas/d]

		for i in range(70):
			fs = (1 + i) / 100
			if (fs == 1):
				continue
	#		F = ((i + 1) / 10)
	#		fs = F / (F + 1)
			generujplik(t, xa, ya, pa, vda, xb, yb, pb, vdb, 0, 0, fs)
			guess = (obliczparalakse()[0] - pa * fs) / (1 - fs) #parallax
			guessvx = (obliczparalakse()[1]*730 - va[0] * fs) / (1 - fs)
			guessvy = (obliczparalakse()[2]*730 - va[1] * fs) / (1 - fs)
	#		file.write('%f\t%f\t%f\t%f\t%f\t%f\n' % (pa, pb, fs, obliczparalakse()[0], guess, np.abs(guess - pb) / pb))
			file.write("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n" % (pa, pb, fs, obliczparalakse()[0], guess, np.abs(guess - pb) / pb, va[0], va[1], vb[0], vb[1], obliczparalakse()[1], obliczparalakse()[2], guessvx, guessvy, np.abs(guessvx-vb[0]), np.abs(guessvy-vb[1])))
	#		file.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (pa, pb, fs, obliczparalakse()[0], guess, np.abs(guess - pb) / pb, vb[0], vb[1], guessvx, guessvy, np.abs(guessvx-vb[0]), np.abs(guessvy-vb[1])))
	file.close()

else:
	print('Podaj wlasciwy numer')