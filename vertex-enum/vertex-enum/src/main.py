import numpy as np
from bellpolytope import BellPolytope


    
if __name__ == '__main__':
 
    poly = BellPolytope(3,2)    
    ineffResistInequalities = poly.getInefficiencyResistantInequalities()
    
    vertices33 = BellPolytope(3,3).getVertices();
    relevantIneqs = filter(lambda ineq : min(
        map(lambda vertex : np.dot(ineq[1:],vertex),vertices33))<=-2,
        ineffResistInequalities)
    
    print(len(relevantIneqs))
    print(relevantIneqs[11][1:])
    