# Google-Cirq
Testing out Google Cirq

qft:
    To run the quantum fourier transform code, import qft, running cirq.Circuit.from_ops(qft(qubits)) will create a quantum fourier transform circuit on a list of qubits given.
    Including an array of 0's and 1's relating to the starting states of qubits as initial = [] will also set the starting states of the qubits to the desired state.
    Adding the variable tol = t will restrict all controlled gates to only act on qubits distance t away.
    Setting mea = True will add measurement at the end of the circuit.
    
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
    For example importing psiadd, then running circuit.append(psiadd(3,qubits), strategy=InsertStrategy.NEW) where qubits is a list of qubits to act on will append n rotational Z gates to the circuit in this case adding 3 to fourier transformed number in qubits.