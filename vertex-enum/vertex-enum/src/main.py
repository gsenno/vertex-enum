import numpy as np
from bellpolytope import BellPolytope
  

if __name__ == '__main__':


    poly = BellPolytope(3,2)    

    ineffResistInequalities = poly.getInefficiencyResistantInequalities()
    vertices33 = BellPolytope(3,3).getVertices();

    relevantIneqs = list(filter(lambda ineq : min(
        map(lambda vertex : np.dot(ineq[1:],vertex),vertices33))<-1,
        ineffResistInequalities))
   

    Ineq=np.zeros((20,81))
    for i in range (0,20):
        for j in range (0,81):
            Ineq[i][j]=relevantIneqs[i][1+j]
    with open('Ineqs.txt','wb') as f:
            np.savetxt(f, Ineq, fmt='%.2f')

    