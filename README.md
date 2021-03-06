# Google-Cirq
Testing out Google Cirq

qft:
    To run the quantum fourier transform code, import qft, running cirq.Circuit.from_ops(qft(qubits)) will create a quantum fourier transform circuit on a list of qubits given.
    Including an array of 0's and 1's relating to the starting states of qubits as initial = [] will also set the starting states of the qubits to the desired state.
    Adding the variable tol = t will restrict all controlled gates to only act on qubits distance t away.
    Setting mea = True will add measurement at the end of the circuit.
    Including control = [] for some list of qubits will add control qubits to the entire circuit.
    
    Example:
        cirq.Circuit.from_ops(qft([cirq.GridQubit(i,0) for i in range(10)], initial = [0,0,0,0,0,0,0,0,1,0], tol = 6, mea = False)) 
        
        This gives a circuit that you can see below, each block only contains at most 6 qubits. Plotting the result of this will give a sine wave of 2 cycles but will be jagged due to the approximation of tol = 6.
        
        (0, 0): ───────H───@───────@────────@─────────@──────────@───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────×───────────────────
                           │       │        │         │          │                                                                                                                                                                                                                                                                                                                                       │
        (1, 0): ───────────@^0.5───┼────────┼─────────┼──────────┼─────────H───@───────@────────@─────────@──────────@───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───×───────────────
                                   │        │         │          │             │       │        │         │          │                                                                                                                                                                                                                                                                                   │   │
        (2, 0): ───────────────────@^0.25───┼─────────┼──────────┼─────────────@^0.5───┼────────┼─────────┼──────────┼─────────H───@───────@────────@─────────@──────────@───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───┼───×───────────
                                            │         │          │                     │        │         │          │             │       │        │         │          │                                                                                                                                                                                                                               │   │   │
        (3, 0): ────────────────────────────@^(1/8)───┼──────────┼─────────────────────@^0.25───┼─────────┼──────────┼─────────────@^0.5───┼────────┼─────────┼──────────┼─────────H───@───────@────────@─────────@──────────@───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───┼───┼───×───────
                                                      │          │                              │         │          │                     │        │         │          │             │       │        │         │          │                                                                                                                                                                           │   │   │   │
        (4, 0): ──────────────────────────────────────@^(1/16)───┼──────────────────────────────@^(1/8)───┼──────────┼─────────────────────@^0.25───┼─────────┼──────────┼─────────────@^0.5───┼────────┼─────────┼──────────┼─────────H───@───────@────────@─────────@──────────@───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼───┼───┼───┼───×───
                                                                 │                                        │          │                              │         │          │                     │        │         │          │             │       │        │         │          │                                                                                                                       │   │   │   │   │
        (5, 0): ─────────────────────────────────────────────────@^0.031──────────────────────────────────@^(1/16)───┼──────────────────────────────@^(1/8)───┼──────────┼─────────────────────@^0.25───┼─────────┼──────────┼─────────────@^0.5───┼────────┼─────────┼──────────┼─────────H───@───────@────────@─────────@──────────────────────────────────────────────────────────────────────────────┼───┼───┼───┼───×───
                                                                                                                     │                                        │          │                              │         │          │                     │        │         │          │             │       │        │         │                                                                              │   │   │   │
        (6, 0): ─────────────────────────────────────────────────────────────────────────────────────────────────────@^0.031──────────────────────────────────@^(1/16)───┼──────────────────────────────@^(1/8)───┼──────────┼─────────────────────@^0.25───┼─────────┼──────────┼─────────────@^0.5───┼────────┼─────────┼──────────H───@───────@────────@──────────────────────────────────────────────┼───┼───┼───×───────
                                                                                                                                                                         │                                        │          │                              │         │          │                     │        │         │              │       │        │                                              │   │   │
        (7, 0): ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────@^0.031──────────────────────────────────@^(1/16)───┼──────────────────────────────@^(1/8)───┼──────────┼─────────────────────@^0.25───┼─────────┼──────────────@^0.5───┼────────┼─────────H───@───────@────────────────────────┼───┼───×───────────
                                                                                                                                                                                                                             │                                        │          │                              │         │                      │        │             │       │                        │   │
        (8, 0): ───X─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────@^0.031──────────────────────────────────@^(1/16)───┼──────────────────────────────@^(1/8)───┼──────────────────────@^0.25───┼─────────────@^0.5───┼────────H───@───────────┼───×───────────────
                                                                                                                                                                                                                                                                                 │                                        │                               │                     │            │           │
        (9, 0): ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────@^0.031──────────────────────────────────@^(1/16)────────────────────────@^(1/8)───────────────@^0.25───────@^0.5───H───×───────────────────

psiadd:
    This function creates a generator that will apply rotational Z gates to a list of given qubits b_qubits such that if b_qubits is the quantum fourier transform of a classical number b, the output will be a + b in the fourier space.
    Setting sin = -1 will instead subtract a from b.
    Including control = [] for some list of qubits will add control qubits to the entire circuit.
    For example importing psiadd then running circuit.append(psiadd(3,qubits), strategy=InsertStrategy.NEW) where qubits is a list of qubits to act on will append n rotational Z gates to the circuit in this case adding 3 to fourier transformed number in qubits.

modadd:
    Applies a + b (MOD N) conditioned on 2 given qubits.
    a is an integer number to be added to b.
    N is the modulo number
    b_qubits is are the qubits that hold the quantum fourier transformed number b.
    anc is a list of qubits, anc[0] should be an ancilla qubit set to 0, anc[1] and anc[2] are control qubits that can be used with other circuits.
    Make sure N > b and that b always has a buffer such that the msf is 0 at the start.
    For example importing modadd then running Circuit.append(modadd(1, 6, qubits, [cirq.GridQubit(3,i) for i in range(3)]),strategy=InsertStrategy.NEW) where qubits are the qubits that have the qft number b stored will return a circuit that gives 1 + b (MOD 6) if anc[1] and anc[2] are both in state 1.
    
modmult:
    Send |x>|b> --> |x>|b + (ax)MOD N> conditioned on a control qubit.
    Calling modmult(a, N, b_qubits, x_qubits, anc) will return a generator where a and N are the values above, b_qubits are the qubits holding the value b, x_qubits hold the value x, anc is a list where anc[0] is an ancilla qubit set to 0 and anc[1] is the control qubit.
    Make sure N > b, x and that b and x always have a buffer such that the msf is 0 at the start.

c_ua:
    Creates a generator for a quantum multiplyer that takes |x>|0> -> |(ax)MOD N>|0> controlled by a qubit c.
    zeros_qubits is a list of ancilla qubits all set to 0.
    x_qubits is a list of qubits where x is stored.
    zeros_qubits should have the same size as x_qubits.
    anc is a list where anc[0] is an ancilla qubit that should be set to 0, anc[1] is the control qubit.
    c_ua(a, N, zeros_qubits, x_qubits, anc) will take |x>_n|0>_n to |(ax)MOD N>_n|0>_n
	
order_find:
	Uses the example in https://arxiv.org/pdf/quant-ph/0205095.pdf to find the order given N and a with n the bit length of N. 
	zeros_qubits is a list of n ancilla qubits all set to 0.
    x_qubits is a list of n qubits set to 0.
    m_qubits is a list of 2n qubits set to 0.
    anc0 is an ancilla qubit that should be set to 0.