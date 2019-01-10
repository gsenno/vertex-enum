import numpy as np
from bellpolytope import BellPolytope
  

if __name__ == '__main__':

    N=3
    K=2
    poly = BellPolytope(N,K)    

    ineffResistInequalities = poly.getInefficiencyResistantInequalities()
    vertices33 = BellPolytope(N,K+1).getVertices();

    relevantIneqs = list(filter(lambda ineq : min(
        map(lambda vertex : np.dot(ineq[1:],vertex),vertices33))<-1,
        ineffResistInequalities))
   
    numberOfCoefficients=N**2*(K+1)**2
    Ineq=np.zeros((len(relevantIneqs),numberOfCoefficients))
    for i in range (0,len(relevantIneqs)):
        for j in range (0,numberOfCoefficients):
            Ineq[i][j]=relevantIneqs[i][1+j]
    with open('Ineqs.txt','wb') as f:
            np.savetxt(f, Ineq, fmt='%.2f')

    