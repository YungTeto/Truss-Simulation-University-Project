# #############################################################################
# Institut fuer Mechanik, TU Dortmund
# Kurs: Mechanik IV
# SS 2021
# -----------------------------------------------------------------------------
# FEM fuer 2D-Fachwerke
# ---> plot routine
#
# #############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def show_results(old_coord, u, conn, nel, nen, ndf, nnp, E):
	coord = old_coord + u  # updated node coordinates
	
	coorde = np.zeros((nel, nen * ndf))
	old_coorde = np.zeros((nel, nen * ndf))
	xe = np.zeros((nel, ndf))
	ye = np.zeros((nel, ndf))
	old_xe = np.zeros((nel, ndf))
	old_ye = np.zeros((nel, ndf))
	sig = np.zeros(nel)
	
	
	#--------------------------------------------------------------------------
	# determine coordinates and stresses of truss elements
	#--------------------------------------------------------------------------
	
	# loop over elements
	for e in range(nel):
		
		# Boolean Matrices
		L = np.zeros((nen*ndf, ndf*nnp))
		for P in range(nen):
			for j in range(ndf):
				L[2*P+j, (conn[e, P]-1)*ndf+j] = 1
		
		
		coorde[e] = np.transpose(np.dot(L, coord))  # save xe
		old_coorde[e] = np.transpose(np.dot(L, old_coord)) # save old xe
		
		
		# save the x and y coordinates of element e
		for i in range(ndf):
			xe[e, i] = coorde[e, 2 * i]
			ye[e, i] = coorde[e, 2 * i + 1]
			old_xe[e, i] = old_coorde[e, 2 * i]
			old_ye[e, i] = old_coorde[e, 2 * i + 1]
		
		# Calculate the old/new length of truss element e
		l0 = np.sqrt(np.power((old_xe[e, 1]-old_xe[e, 0]), 2) + np.power((old_ye[e, 1]-old_ye[e, 0]), 2))
		l1 = np.sqrt(np.power((xe[e, 1]-xe[e, 0]), 2) + np.power((ye[e, 1]-ye[e, 0]), 2))
		
		# Calculate strain of truss element e
		eps = (l1-l0)/l0
		
		# Calculate stress of truss element e
		sig[e] = E * eps
	
	
	#--------------------------------------------------------------------------
	# plot the system, colored in regard to stresses
	#--------------------------------------------------------------------------
	
	# define colormap for plotting stresses
	CMap = mpl.cm.get_cmap('jet')
	normalizedMap = mpl.colors.Normalize(vmin=min(sig), vmax=max(sig))
	
	# plot colorbar as legend
	sm = plt.cm.ScalarMappable(cmap=CMap, norm=normalizedMap)
	cbar = plt.colorbar(sm)
	cbar.set_label('Normalspannungen [MPa]')
	
	# plot elements with color according to stress value
	for e in range(nel):
		color_e = CMap(normalizedMap(sig[e]))
		plt.plot(xe[e], ye[e], color=color_e)
	
	
	
	# Save the nodal coordinates in x and y
	x = np.zeros(int(len(coord) / 2))
	y = np.zeros(int(len(coord) / 2))
	for i in range(int(len(coord)/2)):
		x[i] = coord[2 * i]
		y[i] = coord[2 * i + 1]
	
	# Plot the nodes
	plt.plot(x, y, "ko")
	
	
	
	# finalize plot
	plt.grid()
	plt.show()
