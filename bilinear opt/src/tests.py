'''
Created on 30 dic. 2018

@author: gsenno
'''
def printMatrices(prob):
    print('X00 = ', prob.matsolX(0))
    print('X01 = ', prob.matsolX(1))
    print('X02 = ', prob.matsolX(2))
    print('X00+X01+X02 = ',prob.matsolX(0)+prob.matsolX(1)+prob.matsolX(2))

    print('-------------------------------------------------------------')    
    
    print('X10 = ', prob.matsolX(3))
    print('X11 = ', prob.matsolX(4))
    print('X12 = ', prob.matsolX(5))
    print('X10+X11+X12 = ',prob.matsolX(3)+prob.matsolX(4)+prob.matsolX(5))

    print('-------------------------------------------------------------')    
    
    print('X20 = ', prob.matsolX(6))
    print('X21 = ', prob.matsolX(7))
    print('X22 = ', prob.matsolX(8))
    print('X20+X21+X22 = ',prob.matsolX(6)+prob.matsolX(7)+prob.matsolX(8))
    
    print('-------------------------------------------------------------')    

    print('Y00 = ', prob.matsolY(0))
    print('Y01 = ', prob.matsolY(1))
    print('Y02 = ', prob.matsolY(2))
    print('Y00+Y01+Y02 = ',prob.matsolY(0)+prob.matsolY(1)+prob.matsolY(2))

    print('-------------------------------------------------------------')    
    
    print('Y10 = ', prob.matsolY(3))
    print('Y11 = ', prob.matsolY(4))
    print('Y12 = ', prob.matsolY(5))
    print('Y10+Y11+Y12 = ',prob.matsolY(3)+prob.matsolY(4)+prob.matsolY(5))

    print('-------------------------------------------------------------')    
    
    print('Y20 = ', prob.matsolY(6))
    print('Y21 = ', prob.matsolY(7))
    print('Y22 = ', prob.matsolY(8))
    print('Y20+Y21+Y22 = ',prob.matsolY(6)+prob.matsolY(7)+prob.matsolY(8))