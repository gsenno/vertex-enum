import numpy as np
from bellpolytope import BellPolytope
  

def Alg(functional[i][j],j,H[j-1]):
    
    return H.tolist()


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
    functional=np.zeros((len(relevantIneqs),numberOfCoefficients))
    functional[i]=Ineq[i][:]
    
    with open('Ineqs.txt','w') as f:
                f.write('Minimal representative functional \n')
                
    minimalfunctionals =[]
    permutedfunctional=[]
    for i in range (0,len(relevantIneqs)):
        H=[[] for x in xrange(numberOfCoefficients)]
        for k in range (1, numberOfCoefficients):
            H.append(k) #IDENTITY    

        for j in range (1, numberOfCoefficients+1):
            H=Alg(functional[i][:], j, H)
    
        for l in range (0, len(H)):
            for s in range (0, numberOfCoefficients):
                permutedfunctional[s]=functional[i][H[l][s]]
            minimalfunctionals.append(permutedfunctional)
    
    output = set()
    for x in minimalfunctionals:
        output.add(x)
    print output
    

    