# seesaw_utils.py
# (c) Alejandro Pozas-Kerstjens, 2018
from ncpol2sdpa import *
from qutip import basis, tensor, ket2dm, expect, qeye
from numpy import pi as p
from numpy import e
from numpy.linalg import eigh, multi_dot, svd
import numpy as np
import qutip as qt
import itertools
import pickle
import os
from time import time


def current_time():
    return int(round(time() * 1000))
    
    
def projplus(θ, φ):
    '''Gives the projector on the positive outcome of a qubit observable
    characterized by the Bloch-sphere direction (θ, φ)
    '''
    return qt.Qobj(np.array([np.cos(θ / 2), e **(1j * φ) * np.sin(θ / 2)]))


def get_angles_from_projector(projector):
    '''Gives the Bloch-sphere direction (θ, φ) of a positive-outcome projector.
    It is the inverse of the function projplus(θ, φ)
    '''
    θ = 2 * np.arccos(np.sqrt(projector[0][0]))
    φ = -np.imag(np.log(projector[0][1] / (np.cos(θ / 2) * np.sin(θ / 2))))
    if np.isnan(φ):
        φ = 0
    return θ, φ


def joint_probabilities(state, measA, measB, flag='angles'):
    '''Computes the probability vector p(ab|xy) for a given state and measurements.
    :format state: Qutip two-qubit quantum state, in vector (no density matrix) form
    :format measA: list of measurements. Each element is itself a list of projectors/POVM elements.
                   For two-outcome measurements, the structure is [[m1], [m2]...], where mi
                   characterizes the projector onto the positive outcome of measurement i.
    :format measB: list of measurements. Each element is itself a list of projectors/POVM elements.
                   For two-outcome measurements, the structure is [[m1], [m2]...], where mi
                   characterizes the projector onto the positive outcome of measurement i.
    :format flag: 'angles' if the measurements are given by its Bloch angles (θ, φ), 'projectors'
                  if the measurements are given by the projector matrices themselves.
    '''
    rho = ket2dm(state)
    
    if flag == 'angles':
        θA = [measA[2*i] for i in range(4)]
        φA = [measA[2*i + 1] for i in range(4)]
        θB = [measB[2*i] for i in range(4)]
        φB = [measB[2*i + 1] for i in range(4)]
        A = [[ket2dm(projplus(θ, φ))] for θ, φ in zip(θA, φA)]
        B = [[ket2dm(projplus(θ, φ))] for θ, φ in zip(θB, φB)]
    elif flag == 'projectors':
        A = [[qt.Qobj(a[0])] for a in measA]
        B = [[qt.Qobj(b[0])] for b in measB]
    else:
        raise Exception('The parameter flag should be either \'angles\' or \'projectors\'')

    probabilities = flatten([[expect(tensor(a[0], qeye(2)), rho) for a in A],
                             [expect(tensor(qeye(2), b[0]), rho) for b in B],
                             [expect(tensor(a[0], b[0]), rho) for a, b in itertools.product(A, B)]])

    return probabilities
    
    
def randU(N):
    '''Generates a random NxN unitary matrix, distributed uniformly
    according to the Haar measure.'''
    X = (np.random.randn(N, N) + 1j * np.random.randn(N, N)) / np.sqrt(2)
    Q, R = np.linalg.qr(X)
    R = np.diag(np.diag(R) / abs(np.diag(R)))
    U = np.dot(Q, R)
    return U


def rndmeas(N):
    '''Generates a random binary measurement operator'''
    random_number = np.random.rand()
    r = np.diag([1, 0] if random_number > 0.5 else [0, 1])
    for _ in range(4):
        U = randU(N)
        r = np.dot(U, np.dot(r, U.conj().T))
    return 2 * r - np.eye(N)


def updateA(ineq, stateSch, B):
    '''Updates the measurements of party A according to www.arxiv.org/abs/1006.3032
    :format ineq: matrix where the ij element corresponds to the inequality
                  coefficient corresponding to measurement A_i B_j
    :format stateSch: vector of the Schmidt form of the state. This is computed
                      in seesaw(), but for our case it essentially is
                      [cos(Θ), sin(Θ)]
    :format B: as given in previous functions, [[proj_m1], [proj_m2]...], where now
               we only have projectors.
    '''
    dims = ineq.shape
    X = np.array([[[np.sum([ineq[μ][ν] * stateSch[i] * B[ν][0][i][j] * stateSch[j] for ν in range(dims[1])])
                                                                                   for i in range(2)]
                                                                                   for j in range(2)]
                                                                                   for μ in range(dims[0])])
    for μ in range(len(X)):
        X[μ] = (X[μ] + X[μ].conj().T) / 2    # Fix small rounding errors, the matrix should be Hermitian
    Xdiag = [eigh(X[μ])[0] for μ in range(dims[0])]
    mats = [eigh(X[μ])[1] for μ in range(dims[0])]
    A = [np.diag([1 if a >= 1e-10 else 0 for a in Xdiag[μ+1]]) for μ in range(dims[0] - 1)]
    for i, a in enumerate(A):
        if (np.sum(a) == 2) | (np.sum(a) == 0):
            diagonal = np.zeros((2,))
            diagonal[np.argmax(Xdiag[i + 1])] = 1
            A[i] = np.diag(diagonal)
    A = np.concatenate(([np.eye(2)], A))
    A = np.array([[multi_dot([mats[μ], A[μ], mats[μ].conj().T])] for μ in range(dims[0])])
    return A


def updateB(ineq, stateSch, A):
    '''Updates the measurements of party B according to www.arxiv.org/abs/1006.3032
    :format ineq: matrix where the ij element corresponds to the inequality
                  coefficient corresponding to measurement A_i B_j. The first row
                  and column are for identity operators.
    :format stateSch: vector of the Schmidt form of the state. This is computed
                      in seesaw(), but for our case it essentially is
                      [cos(Θ), sin(Θ)]
    :format A: as given in previous functions, [[proj_m1], [proj_m2]...], where now
               we only have projectors.
    '''
    dims = ineq.shape
    Y = np.array([[[np.sum([ineq[μ][ν] * stateSch[i] * A[μ][0][i][j] * stateSch[j] for μ in range(dims[0])])
                                                                                   for i in range(2)]
                                                                                   for j in range(2)]
                                                                                   for ν in range(dims[1])])
    for ν in range(len(Y)):
        Y[ν] = (Y[ν] + Y[ν].conj().T) / 2    # Fix small rounding errors, the matrix should be Hermitian
    Ydiag = [eigh(Y[ν])[0] for ν in range(dims[1])]
    mats = [eigh(Y[ν])[1] for ν in range(dims[1])]
    B = [np.diag([1 if a >= 1e-10 else 0 for a in Ydiag[ν+1]]) for ν in range(dims[1] - 1)]
    for i, b in enumerate(B):
        if (np.sum(b) == 2) | (np.sum(b) == 0):
            diagonal = np.zeros((2,))
            diagonal[np.argmax(Ydiag[i + 1])] = 1
            B[i] = np.diag(diagonal)
    B = np.concatenate(([np.eye(2)], B))
    B = np.array([[multi_dot([mats[ν], B[ν], mats[ν].conj().T])] for ν in range(dims[1])])
    return B


def seesaw(ineq, state, pre_A=None, pre_B=None):
    '''Performs a seesaw either one or a number of times, and outputs the measurements that maximize the inequality given.
    If measurements are given, the optimization is run once (all the process is deterministic). If not, measurements are
    initialized randomly and the seesaw is run a number of times with different initial measurements so as to try to avoid
    falling into local minima.
    :format state: the state should be in the form of a 2-D matrix such that the coefficient s[i][j] corresponds to the ket
                   |i>|j>. In our case, this matrix is [[cos(Θ), 0], [0, sin(Θ)]]-
    :format ineq: same as before
    '''
    
    # Write the Schmidt decomposition of the state
    (_, d, _) = svd(state)
    dims = ineq.shape
    tries = 1     # Number of times we will run the seesaw
    
    # If no priors are given, initialize random measurements
    if (pre_A == None):
        pre_A = np.array([[rndmeas(2)] for _ in range(dims[0] - 1)])
        tries = 1e4
    if (pre_B == None):
        pre_B = np.array([[rndmeas(2)] for _ in range(dims[0] - 1)])
        tries = 1e4
    
    # Add identity operators
    pre_A = np.concatenate((np.array([[np.eye(2)]]), pre_A))
    pre_B = np.concatenate((np.array([[np.eye(2)]]), pre_B))
    
    # Compute initial Bell operator
    pre_bell = np.real(np.sum([[[[ineq[μ][ν] * d[i] * pre_A[μ][0][i][j] * pre_B[ν][0][i][j] * d[j] for ν in range(dims[1])]
                                                                                                   for μ in range(dims[0])]
                                                                                                   for i in range(2)]
                                                                                                   for j in range(2)]))
    # Do first iteration for having conditions to begin the loop
    A = updateA(ineq, d, pre_B)
    B = updateB(ineq, d, A)
    bell = np.real(np.sum([[[[ineq[μ][ν] * d[i] * A[μ][0][i][j] * B[ν][0][i][j] * d[j] for ν in range(dims[1])]
                                                                                       for μ in range(dims[0])]
                                                                                       for i in range(2)]
                                                                                       for j in range(2)]))
    
    tol = 1e-6 if d[0] < 0.99 else 1e-7
    solutions = []
    bells = []
    for idx in range(int(tries)):
        # If we are running more than one try, use the measurements obtained before
        # for the first one but initialize new in each run
        if ((tries > 1) and (idx != 1)):
            pre_A = np.array([[rndmeas(2)] for _ in range(dims[0] - 1)])
            pre_A = np.concatenate((np.array([[np.eye(2)]]), pre_A))
            pre_B = np.array([[rndmeas(2)] for _ in range(dims[0] - 1)])
            pre_B = np.concatenate((np.array([[np.eye(2)]]), pre_B))
        count = 0
        # The seesaw itself
        while (np.abs(pre_bell - bell) > tol) & (count < 1e6):
            pre_A = A
            pre_B = B
            pre_bell = bell
            A = updateA(ineq, d, B)
            B = updateB(ineq, d, A)
            bell = np.real(np.sum([[[[ineq[μ][ν] * d[i] * A[μ][0][i][j] * B[ν][0][i][j] * d[j] for ν in range(dims[1])]
                                                                                               for μ in range(dims[0])]
                                                                                               for i in range(2)]
                                                                                               for j in range(2)]))
            count +=1
        solutions.append([A, B])
        bells.append(bell)
    bells = max(bells)
    A, B = solutions[np.argmax(bells)]
    
    return bell, A, B, count