import numpy as np

def get_data():
    #l=9.25
    x = np.array([749,     5*749,         
                  0,       4*749,
                  749,     4*749,
                  0,       3*749,
                  749,     3*749,
                  2*749,   2*749,
                  5*749,   2*749,
                  0,       749,
                  749,     749,
                  3*749,   749,
                  4*749,   749,
                  5*749,   749,
                  0,       0,
                  2*749,   0,
                  4*749,   0,
                  ])
    #x = x*l
    
    conn = np.array([(1, 2),
                      (1 ,3),
                      (2, 3),
                      (2, 4),
                      (3, 4),
                      (3, 5),
                      (4, 5),
                      (4, 8),
                      (4, 9),
                      (5, 9),
                      (5, 6),
                      (8, 9),
                      (6, 9),
                      (9, 10),
                      (6, 10),
                      (10, 11),
                      (7, 11),
                      (7, 12),
                      (11, 12),
                      (8, 13),
                      (9, 13),
                      (9, 14),
                      (10, 14),
                      (10, 15),
                      (11, 15),
                      (12, 15),
                      (13, 14),
                      (14, 15)]) 
             
    
#integration points
    
    xi = np.array((-1.0/np.sqrt(3), 1.0/np.sqrt(3)))
    w8 = np.array((1, 1))
    
    drltDoFs = np.array([1, 2, 3, 4, 7, 8])
    ud       = np.array([0, 0, 0, 0, 0, 0])



    freeDoFs = np.array([5, 6, 9, 10,11,12, 13,                          14,                            15, 16, 17, 18, 19, 20, 21,   22, 23, 24, 25, 26, 27, 28, 29, 30])
    fpre     = np.array([0, 0, 0, 0, 0, 0,  np.cos(96*(np.pi/180))*7000 ,np.sin(96*(np.pi/180))*7000  ,  0,  0,  0,  0,  0,  0, -6538, 0, 0,  0,  0,  0,  0,  0,   0,  0, ])
    


    nnp = 15 #number of node points
    ndf = 2 #number of degrees of freedom
    ndm = 2 #number of dimensions
    nel = 28 #number of elements
    nen = 2 #number of element notes
    nqp = 2 #AnzahlIntegrationspunkte
    
    
    # Material parameters
    
    E = 210000 #Emodul
    A = 150 #Querschnittsfläche Stäbe
    
    
    
    
    
    return nnp, ndf, ndm , nel, nen, nqp,  drltDoFs, freeDoFs, ud, x, conn, fpre, xi, w8, E, A