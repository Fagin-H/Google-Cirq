import cirq
from cirq import Simulator
from cirq.circuits import InsertStrategy
from circuitgenerators import order_find
import matplotlib.pyplot as plt
import numpy as np

'''
Creates functions for figuring out s
'''
def ntoc ( r, n ):

    a = []

    if ( r == 0.0 ):
        return a

    r2 = r
    a.append(int ( r2 ))

    for i in range ( 1, n + 1 ):
        if r2 - float(a[i-1]) == 0:
            return a
        else:
            r2 = 1.0 / ( r2 - float ( a[i-1] ) )
            a.append(int ( r2 ))

    return a



def ctond ( n, a ):

    p = np.zeros ( n + 1)
    q = np.zeros ( n + 1)

    for i in range ( n + 1 ):
        if p[i] or q[i] > 10000:
            break
        else:
            if ( i == 0 ):
                p[i] = a[i] * 1 + 0
                q[i] = a[i] * 0 + 1
            elif ( i == 1 ):
                p[i] = a[i] * p[i-1] + 1
                q[i] = a[i] * q[i-1] + 0
            else:
                p[i] = a[i] * p[i-1] + p[i-2]
                q[i] = a[i] * q[i-1] + q[i-2]

    return p[:-1], q[:-1]

def ntond(r, n): # Given a number r and a number of itterations n will return possible numerator and denomiator values.
    values = ntoc ( r, n )
    return ctond ( len(values) - 1, values )
       
def egcd(a, b): # Finds the gcd of 2 numbers.
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def find_y(N, a, q, reps = 5): # Uses the order finding algorithm to give possible values for y
    
    simulator = Simulator() # Sets up simulator
    circuit = cirq.Circuit()

    number_qubits = q # Sets up all qubits to be used.
    m_qubits = [cirq.GridQubit(i,0) for i in range(2*number_qubits)]
    x_qubits = [cirq.GridQubit(i,0) for i in range(2*number_qubits, 3*number_qubits)]
    b_qubits = [cirq.GridQubit(i,0) for i in range(3*number_qubits, 4*number_qubits + 1)]
    zero_qubit = cirq.GridQubit(4*number_qubits + 1,0)
    print('Number of qubits: {}'.format(4*number_qubits + 5))
    
    circuit.append(order_find(a, N, b_qubits, x_qubits, zero_qubit, m_qubits), strategy=InsertStrategy.NEW) # Runs the simulation.
    
    return simulator.run(circuit, repetitions = reps).histogram(key = 'q') # Returns a dictionary of the measurment results.


def find_r(Q, N, y): # Finds possible values for r given a y value.
    rs = []
    if y != 0: # Makes sure y is not equal to 0
        test_yQ = y/Q # Finds y/Q to find d/s.
        
        d, s = ntond(test_yQ, 20) # Finds possible values for d and s.
        
        for i in range(len(s)):
            gcd = egcd(d[i], s[i])[0] # Reduces d/s to its smallest form.
            s[i] = s[i]/gcd
            
            if np.absolute(y/Q-d[i]/s[i]) < 1/(2*Q) and s[i] < N: # Checks conditions on s before adding to the list of possible values.
                rs.append(s[i])
    return rs

def factor_N(N, a, reps = 5):
        
    for i in range(N): # Finds Q and q such that N^2 < Q < 2N^2.
        q = i
        Q = 2**q
        if N**2 <= Q:
            break
        
    possible_ys = find_y(N, a, q, reps) # Finds possible values for y using the quantum order finder
    
    possible_rs = []
    
    for y in possible_ys: # Finds possible values for r with the given values of y.
        rs = find_r(Q, N, y)
        for r in rs:
            possible_rs.append(r)
            
    possible_rs = list(dict.fromkeys(possible_rs)) # Deletes repeated values in possible_rs.
    
    rs = []
    
    for r in possible_rs: # Finds values of r that can be used to find factors of N
        if r%2 == 0 and (a**(r/2))%N != N-1:
            rs.append(r)
    
    if len(rs) == 0:
        return False, rs
    else:
        return True, rs
    
    
def shors(N, a = 0, reps =  5, attemps = 1):

    if a == 0: # Picks random value for a that doesnt alreay share a factor.
        a = np.random.randint(2,N)
        while egcd(a, N)[0] != 1:
            a = np.random.randint(2,N)
    
    for i in range(attemps): # Runs the factoring algorithm untill it succseds or runs out of attempts.
        found, rs = factor_N(N, a, reps)
        if found:
            break
    
    if found:
        
        return True, (egcd(a**(rs[0]/2) + 1, N), egcd(a**(rs[0]/2) - 1, N))
    else:
        return False, []
    
    
    
def doall(N, a = 0, reps = 5, attemps = 1):
    found, factors = shors(N, a, reps, attemps)

    if found:
        print('{}/{} = {}'.format(N,factors[0],N/factors[0]))
        print('{}/{} = {}'.format(N,factors[1],N/factors[1]))
        print()
    else:
        print('Failed to find factor')
        print()






















