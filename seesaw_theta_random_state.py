from __future__ import division
from ncpol2sdpa import *
from qutip import basis, tensor, ket2dm, expect, qeye
from numpy import pi as π
from numpy import e
from numpy.linalg import eigh, multi_dot, svd
import numpy as np
import qutip as qt
import itertools
import pickle
from time import time, strftime
from seesaw_utils import *

def get_maximum_guessing_probability(state, measA, measB, first=False, flag='angles'):
    probabilities = joint_probabilities(state, measA, measB, flag)

    moments = []
    symbols = []
    k = 0
    #Marginals first
    for i in range(4):
        moments.append(P_0_A[i][0]+P_1_A[i][0]+P_2_A[i][0]+P_3_A[i][0]-probabilities[k])
        k += 1
    
    for j in range(4):
        moments.append(P_0_B[j][0]+P_1_B[j][0]+P_2_B[j][0]+P_3_B[j][0]-probabilities[k])
        k += 1

    #Joint probabilities
    for i in range(4):
        for j in range(4):
            moments.append(P_0_A[i][0]*P_0_B[j][0]+P_1_A[i][0]*P_1_B[j][0]+P_2_A[i][0]*P_2_B[j][0]+P_3_A[i][0]*P_3_B[j][0]-probabilities[k])
            k += 1
    # We need to normalize the top-left elements of the moment matrices
    moments.append("0[0,0]+1[0,0]+2[0,0]+3[0,0]-1")
    if first:
        sdpRelaxation.get_relaxation(level, objective=guessing_probability,
                                     momentequalities=moments,
                                     substitutions=substitutions,
                                     extraobjexpr=extraobjexpr)
    else:
        sdpRelaxation.process_constraints(momentequalities=moments,)

    sdpRelaxation.solve(solver='mosek')
    return sdpRelaxation, moments, probabilities
    

level = 2             # Level of the hierarchy you ar going to use

A_configuration = [2, 2, 2, 2]
B_configuration = [2, 2, 2, 2]

# We do four copies à la Nieto-Silleras
P_0_A = generate_measurements(A_configuration, 'P_0_A')
P_0_B = generate_measurements(B_configuration, 'P_0_B')
substitutions1 = projective_measurement_constraints(P_0_A, P_0_B)
P_1_A = generate_measurements(A_configuration, 'P_1_A')
P_1_B = generate_measurements(B_configuration, 'P_1_B')
substitutions2 = projective_measurement_constraints(P_1_A, P_1_B)
P_2_A = generate_measurements(A_configuration, 'P_2_A')
P_2_B = generate_measurements(B_configuration, 'P_2_B')
substitutions3 = projective_measurement_constraints(P_2_A, P_2_B)
P_3_A = generate_measurements(A_configuration, 'P_3_A')
P_3_B = generate_measurements(B_configuration, 'P_3_B')
substitutions4 = projective_measurement_constraints(P_3_A, P_3_B)
substitutions = {**substitutions1, **substitutions2, **substitutions3, **substitutions4}

# Which settings randomness is certified from
x, y = 0, 0

# Parameters for having a random state in a fixed range
init = 0.2
width = 0.1

# Introduce state and measurements
Θ = (π / 4) * (init + width * np.random.uniform())
psi = np.cos(Θ) * tensor([basis(2, 0), basis(2, 0)]) + np.sin(Θ) * tensor([basis(2, 1), basis(2, 1)])

# Angles, coordinate corresponds to measurement setting
θA = np.arccos(1 - 2 * np.random.uniform(0, 1, (4, 1)))
φA = 2 * π * np.random.uniform(0, 1, (4, 1))
θB = np.arccos(1 - 2 * np.random.uniform(0, 1, (4, 1)))
φB = 2 * π * np.random.uniform(0, 1, (4, 1))
measA = flatten([[th, ph] for th, ph in zip(θA, φA)])
measB = flatten([[th, ph] for th, ph in zip(θB, φB)])

guessing_probability = -(P_0_A[x][0]*P_0_B[y][0]
                         + (P_1_B[y][0] - P_1_A[x][0]*P_1_B[y][0])
                         + (P_2_A[x][0] - P_2_A[x][0]*P_2_B[y][0]) 
                         + (-P_3_A[x][0] - P_3_B[y][0] + P_3_A[x][0]*P_3_B[y][0]))
extraobjexpr = "-3[0,0]"

sdpRelaxation = SdpRelaxation([flatten([P_0_A, P_0_B]),
                               flatten([P_1_A, P_1_B]),
                               flatten([P_2_A, P_2_B]),
                               flatten([P_3_A, P_3_B])],
                               normalized=False, verbose=1)
                               
sdp, moments, _ = get_maximum_guessing_probability(psi, measA, measB, first=True, flag='angles')
print('Pguess = ' + str(-sdp.primal))

info_to_RBM = []
anglesA = np.array(measA).T
anglesB = np.array(measB).T
info_to_RBM.append([Θ, anglesA, anglesB, -sdp.primal])

ineq = 0
dual_values = []
for i, moment in enumerate(moments[:-1]):
    value = sdp.get_dual(moment)[0][0,0]-sdp.get_dual(moment)[1][0,0]
    dual_values.append(value)
    
# Ineq = np.block([
                 # [0, np.array(dual_values[:4])],
                 # [np.array(dual_values[4:8]).reshape((4,1)), np.array([[dual_values[8+4*j+i]for j in range(4)] for i in range(4)])]
                 # ])
                 
RowBlock1 = np.hstack(([0], np.array(dual_values[:4])))
RowBlock2 = np.hstack((np.array(dual_values[4:8]).reshape((4,1)), np.array([[dual_values[8+4*j+i]for j in range(4)] for i in range(4)])))
Ineq = np.vstack((RowBlock1, RowBlock2))

θA = [measA[2*i] for i in range(4)]
φA = [measA[2*i + 1] for i in range(4)]
θB = [measB[2*i] for i in range(4)]
φB = [measB[2*i + 1] for i in range(4)]

A_init = [[ket2dm(projplus(θ, φ)).full()] for θ, φ in zip(θA, φA)]
B_init = [[ket2dm(projplus(θ, φ)).full()] for θ, φ in zip(θB, φB)]

psi_matrix = psi.full().reshape((2, 2))

bell, A_optim, B_optim, counts = seesaw(Ineq, psi_matrix)

measA = [meas for meas in A_optim[1:]]
measB = [meas for meas in B_optim[1:]]

sdp2, moments2, _ = get_maximum_guessing_probability(psi, measA, measB, first=False, flag='projectors')
print('Pguess = ' + str(-sdp2.primal))
anglesA = np.array([get_angles_from_projector(meas[0]) for meas in measA])
anglesB = np.array([get_angles_from_projector(meas[0]) for meas in measB])
info_to_RBM.append([Θ, anglesA, anglesB, -sdp2.primal])

count = 0
loop_count = 0
best_pguess = np.array([sdp.primal, sdp2.primal]).max()
sdpRelaxation.verbose = 0    # Set for avoiding all the information about the solution of the SDP
while (count < 50) & (loop_count < 1e5):
    sdp = sdp2
    moments = moments2
    dual_values = []
    for i, moment in enumerate(moments[:-1]):
        value = sdp.get_dual(moment)[0][0,0]-sdp.get_dual(moment)[1][0,0]
        dual_values.append(value)
    RowBlock1 = np.hstack(([0], np.array(dual_values[:4])))
    RowBlock2 = np.hstack((np.array(dual_values[4:8]).reshape((4,1)), np.array([[dual_values[8+4*j+i]for j in range(4)] for i in range(4)])))
    Ineq = np.vstack((RowBlock1, RowBlock2))
    # Ineq = np.block([
                 # [0, np.array(dual_values[:4])],
                 # [np.array(dual_values[4:8]).reshape((4,1)), np.array([[dual_values[8+4*j+i]for j in range(4)] for i in range(4)])]
                 # ])
    A_init = measA
    B_init = measB
    _, A_optim, B_optim, _ = seesaw(Ineq, psi_matrix)
    
    measA = [meas for meas in A_optim[1:]]
    measB = [meas for meas in B_optim[1:]]
    sdp2, moments2, _ = get_maximum_guessing_probability(psi, measA, measB, first=False, flag='projectors')
    if abs(sdp2.primal) < abs(best_pguess):
        best_pguess = sdp2.primal
        count = 0
    else:
        count += 1
    anglesA = flatten([get_angles_from_projector(meas[0]) for meas in measA])
    anglesB = flatten([get_angles_from_projector(meas[0]) for meas in measB])
    info_to_RBM.append([Θ, anglesA, anglesB, -sdp2.primal])
    loop_count += 1
    if (loop_count % 500) == 0:
        with open('seesaw_' + strftime('%d%H%M%S') + '.pickle', 'wb') as f:
            pickle.dump(info_to_RBM, f)
    print(loop_count, count, 'Pguess = ' + str(-sdp2.primal), 'Best = ' + str(-best_pguess))
    
with open('seesaw_' + strftime('%d%H%M%S') + '.pickle', 'wb') as f:
        pickle.dump(info_to_RBM, f)