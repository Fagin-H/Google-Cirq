from cirq.ops import CZPowGate as CZP, H, X, SWAP, measure, ZPowGate as ZP
from numpy import zeros
import cirq

# Creates a generator for a quantum fourier transform circuit.
# qubts is a list of qubits in a linear array.
# tol is the maximum number of qubits the controlled gates will be applied, lower tol will lead to less accurate results.
# Set mea to True if you want it to add measurements at the end.
# initial takes an array of 0's and 1's relating to the starting state of the qubits, leave blank to set to all 0's.
# control is a list of qubits to be used as control bits.
def qft(qubits, tol = None, mea = False, initial = None, control = ()):
    size = len(qubits)
    if tol == None: # Sets maximum distance of controlled gates to cover entire circuit if not given.
        tol = size
    if initial == None: # Sets array for qubit starting state to all 0's if not given.
        initial = zeros(size)
        
    for i in range(size): # Applies an X gate to all qubits starting in the state 1.
        if initial[i] == 1:
            yield cirq.ControlledGate(X, num_controls = len(control))(*control, qubits[i])
    
    for i in range(size): # Applies Hadamard and CPhase gates to each block of qubits.
        yield cirq.ControlledGate(H, num_controls = len(control))(*control, qubits[i])
        for j in range(1,min(size - i, tol)):
            yield cirq.ControlledGate(CZP(exponent = 1/(2**(j))), num_controls = len(control))(*control, qubits[i], qubits[i+j])
            
    for i in range(size//2): # Applies swap gates to switch order of qubits at the end of the circuit.
        yield cirq.ControlledGate(SWAP, num_controls = len(control))(*control, qubits[i],qubits[size-i-1])
    
    for i in range(size): # Adds measurement operators at the end of the circuit if mea is set to true.
        if mea:
            yield measure(qubits[i], key='q{}'.format(str(i)))

# Creates a generator for a quantum adder, adding classical input a to a quantum fourier transformed number b.
# a is an integer number to be added to b.
# b_qubits is are the qubits that hold the quantum fourier transformed number b.
# control is a list of qubits to be used as control bits.
def psiadd(a, b_qubits, sin = 1, control = ()):
    quno = len(b_qubits) # This section of code turns the number a into a binary number then into a list, adding extra zeros if it is too small.
    a_ = [int(i) for i in bin(a)[2:]]
    a_.reverse()
    if len(a_) < quno:
        for i in range(quno-len(a_)):
            a_.append(0)
    
    
    for i in range(quno): # This applies classically calculated Z rotations to b_qubits based on a.
        w = 0
        for j in range(quno-i):
            if a_[j] == 1:
                w += 1/2**(quno - j - i - 1)
        yield cirq.ControlledGate(ZP(exponent = sin * w), num_controls = len(control))(*control, b_qubits[quno-i-1])