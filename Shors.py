import cirq
from cirq.ops import CZ, H, Z, X, ZPowGate as ZP, CZPowGate as CZP
from cirq import Simulator
from cirq.circuits import InsertStrategy
from circuitgenerators import qft, psiadd, modadd, modmult, c_ua, order_find
import matplotlib.pyplot as plt
import numpy as np

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

def ntond(n, r):
    values = ntoc ( r, n )
    return ctond ( len(values) - 1, values )
       
def egcd(a, b): # Finds the gcd of 2 numbers.
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

simulator = Simulator()
circuit = cirq.Circuit()


N = 35

checks = True

while checks:
    a = np.random.randint(2,N)
    if egcd(a,N)[0] == 1:
        checks = False
        
for q in range(N):
    Q = 2**q
    if N**2 <= Q and Q < 2*N**2:
        break


number_qubits = a.bit_length()

m_qubits = [cirq.GridQubit(i,0) for i in range(2*number_qubits)]
x_qubits = [cirq.GridQubit(i,0) for i in range(2*number_qubits, 3*number_qubits)]
b_qubits = [cirq.GridQubit(i,0) for i in range(3*number_qubits, 4*number_qubits + 1)]
zero_qubit = cirq.GridQubit(4*number_qubits + 1,0)


circuit.append(order_find(a, N, b_qubits, x_qubits, zero_qubit, m_qubits), strategy=InsertStrategy.NEW)
print('Made circuit')

results = simulator.run(circuit, repetitions=20)
#print(results.histogram(key = 'q'))

top_results = [i for i in results.histogram(key = 'q')]

rs = []

for y in top_results:
    if y != 0:
        test_yQ = y/Q
        
        values = ntoc(test_yQ, 20)
        d, s = ctond(len(values)-1, values)
        
        for i in range(len(s)):
            gcd = egcd(d[i], s[i])[0]
            s[i] = s[i]/gcd
            if a**s[i]%N == 1 and s[i] < N:
                rs.append(s[i])
            
rs = list(dict.fromkeys(rs))

true_results = []

for test_r in rs:
    if egcd(a**(test_r/2) + 1, N)[0] != N and test_r%2 == 0 and egcd(a**(test_r/2) - 1, N)[0] != N:
        true_results.append((test_r,egcd(a**(test_r/2) + 1, N)[0],egcd(a**(test_r/2) - 1, N)[0]))

if len(true_results) > 0:
    print('r = {}'.format(true_results[0][0]))
    print('a = {}'.format(a))
    print('{}/{} = {}'.format(N,true_results[0][1],N/true_results[0][1]))
    print('{}/{} = {}'.format(N,true_results[0][2],N/true_results[0][2]))
    print()
else:
    print('Failed to find factor')
    print()






















