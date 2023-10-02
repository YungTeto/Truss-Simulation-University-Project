import numpy as np
from Input import get_data
import Graph


# shape funtion for 1d quadratic truss elements
def shape1d_lin(xi_local):
    N1 =  0.5*(1.0-xi_local)
    N2 =  0.5*(1.0+xi_local)
    
    
    N = np.array([N1, N2])
    
    gamma1 = -0.5
    gamma2 = 0.5
    
    gamma = np.array([gamma1, gamma2])
    
    return N, gamma 
    
# Jacobian matrix
def Jacobian(H_local, x_local, gamma_local):
    
    J = gamma_local[0]*np.dot(H_local[0,:], x_local)+gamma_local[1]*np.dot(H_local[1,:], x_local)
    
    return J

#get_data()
nnp, ndf, ndm , nel, nen, nqp, drltDoFs, freeDoFs, ud, x, conn, fpre, xi, w8, E, A = get_data()

u = np.zeros(nnp*ndf)
fsur = np.zeros(nnp*ndf)
K = np.zeros((nnp*ndf, nnp*ndf))

#element loop
for e in range(nel):
    
    
    
    #Boolean Matrix nen*ndf x nnp*ndf
    L = np.zeros((nen*ndf, ndf*nnp))
    for node in range(nen):
        for dof in range(ndf):
            L[ndf*node+dof, (conn[e, node]-1)*ndf+dof] = 1
                        

    #local coordinates of each element
    xe = np.transpose(np.dot(L, x))
    
    
    #auxiliary matrix (Hilfsmatrix)
    phi = np.arctan2((xe[3] - xe[1]), (xe[2] - xe[0]))
    
    H = np.array([(np.cos(phi), np.sin(phi), 0, 0),
                  (0, 0, np.cos(phi),np.sin(phi))])
    
    
    # initialize element stiffness matrix Ke
    Ke = np.zeros((nen*ndf, nen*ndf))
    
    # Integrationsschleife / integration loop
    for q in range(nqp):
        
        # Ansatzfunktion und dessen Ableitung / shape function an derivate
        N, gamma = shape1d_lin(xi[q])
        
        
        # Jacobi Matrix
        J = Jacobian(H, xe, gamma)
        
        # element stiffness Matrix / Element Steifigkeitsmatrix
        Z = gamma[0]*H[0,:]+gamma[1]*H[1,:]
        K_quer = E*A*np.outer(Z, Z)       
        Ke = Ke + K_quer*w8[q]/J
        
        
    # assamble Ke into global stiness matrix K
    K = K + np.dot(np.transpose(L), np.dot(Ke, L))
        
    
    
# initialize partial stiffness matrices
KFF = np.zeros((len(freeDoFs), len(freeDoFs)))
KFD = np.zeros((len(freeDoFs), len(drltDoFs)))
KDD = np.zeros((len(drltDoFs), len(drltDoFs)))


# rearrange K
for i, dof in enumerate(freeDoFs):
    for j, dof_j in enumerate(freeDoFs):
        KFF[i, j] = K[dof-1, dof_j-1]
        
for i, dof in enumerate(freeDoFs):
    for j, dof_j in enumerate(drltDoFs):
        KFD[i, j] = K[freeDoFs[i]-1, drltDoFs[j]-1]
        
for i, dof in enumerate(drltDoFs):
    for j, dof_j in enumerate(drltDoFs):
        KDD[i, j] = K[drltDoFs[i]-1, drltDoFs[j]-1]
        
KDF = np.transpose(KFD)
              
        
        
# solve K*u=f
uf = np.linalg.solve(KFF, fpre - np.dot(KFD, ud))
freac = np.dot(KDF, uf)+np.dot(KDD, ud)



# assemble global vector u
for i, dof in enumerate(freeDoFs):
    u[dof-1] = uf[i]
    
for i, dof in enumerate(drltDoFs):
    u[dof-1] = ud[i]

for i, dof in enumerate(freeDoFs):
    fsur[dof-1] = fpre[i]

for i, dof in enumerate(drltDoFs):
    fsur[dof-1] = freac[i]
    
    
    
    
# generate stress plots
Graph.show_results(x, u, conn, nel, nen, ndf, nnp, E)


np.savez("STUsolution_", x=x, xi=xi, w8=w8, drltDoFs=drltDoFs, ud=ud, freeDoFs=freeDoFs, fpre=fpre, nnp=nnp, ndf=ndf, ndm=ndm, nel=nel, nen=nen, nqp=nqp, E=E, A=A, K=K, uf=uf, freac=freac, u=u, fsur=fsur)











