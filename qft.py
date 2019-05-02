from cirq.ops import CZPowGate as CZP, H, X, SWAP, measure
from numpy import zeros

# Creates a generator for a quantum fourier transform circuit
# qubts is a list of qubits in a linear array
# tol is the maximum number of qubits the controlled gates will be applied, lower tol will lead to less accurate results
# Set mea to True if you want it to add measurements at the end
# initial takes an array of 0's and 1's relating to the starting state of the qubits, leave blank to set to all 0's
def qft(qubits, tol = None, mea = False, initial = None):
    size = len(qubits)
    if tol == None: # Sets maximum distance of controlled gates to cover entire circuit if not given
        tol = size
    if initial == None: # Sets array for qubit starting state to all 0's if not given
        initial = zeros(size)
        
    for i in range(size): # Applies an X gate to all qubits starting in the state 1
        if initial[i] == 1:
            yield X(qubits[i])
    
    for i in range(size): # Applies Hadamard and CPhase gates to each block of qubits
        yield H(qubits[i])
        for j in range(1,min(size - i, tol)):
            yield CZP(exponent = 1/(2**(j)))(qubits[i],qubits[i+j])
            
    for i in range(size//2): # Applies swap gates to switch order of qubits at the end of the circuit
        yield SWAP(qubits[i],qubits[size-i-1])
    
    for i in range(size): # Adds measurement operators at the end of the circuit if mea is set to true
        if mea:
            yield measure(qubits[i], key='q{}'.format(str(i)))