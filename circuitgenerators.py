from cirq.ops import CZPowGate as CZP, H, X, SWAP, measure, ZPowGate as ZP, CNOT as CX
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


# Creates a generator for a quantum adder modulo N, adding classical input a to a quantum fourier transformed number b.
# a is an integer number to be added to b.
# N is the modulo number
# b_qubits is are the qubits that hold the quantum fourier transformed number b.
# anc is a list of qubits, anc[0] should be an ancilla qubit set to 0, anc[1] and anc[2] are control qubits that can be used with other circuits.
def modadd(a, N, b_qubits, anc):
    
    def reversecir(cir): # Creates a function to reverse a quantum circuit, used to calculate the inverse qft.
        cir = [i for i in cir]
        cir.reverse()
        for gate in cir:
           yield gate**-1
    
    for gate in psiadd(a, b_qubits, 1, [anc[1], anc[2]]): # Adds a to b in the fouier space conditioned on the 2 qubits in anc.
        yield gate
    for gate in psiadd(N, b_qubits, -1): # Subtracts N from b in the fouier space.
        yield gate
    
    for gate in reversecir(qft(b_qubits)): # Performs the inverse qft, applies turns anc[0] to 1 if the msf is 1, then applies the qft.
        yield gate
    yield CX(b_qubits[0], anc[0])
    for gate in qft(b_qubits):
        yield gate
        
    for gate in psiadd(N, b_qubits, 1, [anc[0]]): # Adds N to b in the fouier space conditioned on anc[0].
        yield gate
    for gate in psiadd(a, b_qubits, -1, [anc[1], anc[2]]): # Subtracts a from b in the fouier space conditioned on the 2 qubits in anc.
        yield gate
    
    for gate in reversecir(qft(b_qubits)): # Performs the inverse qft, turns anc[0] to 0 ensuring it will end in the same state, then applies the qft.
        yield gate
    yield X(b_qubits[0])
    yield CX(b_qubits[0], anc[0])
    yield X(b_qubits[0])
    for gate in qft(b_qubits):
        yield gate
        
    for gate in psiadd(a, b_qubits, 1, [anc[1], anc[2]]): # Adds a to b in the fouier space conditioned on the 2 qubits in anc.
        yield gate
























