import numpy as np
from bellpolytope import BellPolytope





if __name__ == '__main__':
    N=4
    K=3
    poly = BellPolytope(N,K)  
    
    vertices=poly.getVertices()
    distributions=np.matrix(vertices)
    NumberOfCoefficients=(N*K)**2
    
    with open('randomfunctionals.txt','w') as f:
        f.write('Random functionals \n')
            
    for j in range (0,10):
        functional=np.random.uniform(-63,1,size=NumberOfCoefficients)
    
        values=np.dot(distributions,np.transpose(functional))
        c=np.amax(values)
    
        BellNormalisedFunctional=functional/c
        with open('randomfunctionals.txt','a') as f:
            np.savetxt(f, BellNormalisedFunctional, fmt='%.2f')
            f.write('\n')        