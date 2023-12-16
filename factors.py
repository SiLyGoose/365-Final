from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.circuit.library import QFT
from qiskit.primitives import Sampler
from fractions import Fraction
from math import gcd

def c_amod15(a, N):
    """
    Controlled multiplication by a mod 15.
    This is hard-coded for simplicity.
    """
    # for i in range(1, N):
    #     if gcd(i, N) != 1:
    #         print(i)

    if a not in [2, 4, 7, 8, 11, 13]:
        raise ValueError("'a' must not have common factors with 15")
    U = QuantumCircuit(4)
    if a in [2, 13]:
        U.swap(2, 3)
        U.swap(1, 2)
        U.swap(0, 1)
    if a in [7, 8]:
        U.swap(0, 1)
        U.swap(1, 2)
        U.swap(2, 3)
    if a in [4, 11]:
        U.swap(1, 3)
        U.swap(0, 2)
    if a in [7, 11, 13]:
        for q in range(4):
            U.x(q)
    
    U = U.to_gate()
    U.name = f"{a} mod {N}"
    c_U = U.control()
    return c_U

def phase_estimation(
        controlled_operation: QuantumCircuit,
        psi_prep: QuantumCircuit,
        precision: int
    ):
    """
    Carry out phase estimation on a simulator.
    Args:
        controlled_operation: The operation to perform phase estimation on,
                              controlled by one qubit.
        psi_prep: Circuit to prepare |ψ>
        precision: Number of counting qubits to use
    Returns:
        float: Best guess for phase of U|ψ>
    """
    control_register = QuantumRegister(precision)
    output_register = ClassicalRegister(precision)

    target_register = QuantumRegister(psi_prep.num_qubits)
    qc = QuantumCircuit(control_register, target_register, output_register)

    # Prepare |ψ>
    qc.compose(psi_prep,
               qubits=target_register,
               inplace=True)

    # Do phase estimation
    for index, qubit in enumerate(control_register):
        qc.h(qubit)
        for _ in range(2**index):
            qc.compose(
                controlled_operation,
                qubits=[qubit] + list(target_register),
                inplace=True,
            )

    qc.compose(
        QFT(precision, inverse=True),
        qubits=control_register,
        inplace=True
    )

    qc.measure(control_register, output_register)

    measurement = Sampler().run(qc, shots=1).result().quasi_dists[0].popitem()[0]
    return measurement / 2**precision


if __name__ == "__main__":
    """
    Shor's algorithm to find the factor of a number consists of two parts.
    1. Reduction of the factoring problem to the problem of order-finding. (Classical Step)
    2. Quantum algorithm to solve the order-finding problem. (Quantum Step)
    """
    psi_prep = QuantumCircuit(4) # quantum circuit with 4 qubits
    """
    X gate is applied to the first cubit
    """
    psi_prep.x(0)

    """
    [Classical Step]
    1. Pick a pseudo-random number such that a < N and a >= 2.
    2. Compute gcd(a, N). Can be done using Euclidean algorithm.
    3. If gcd(a, N) != 1, then there is a non-trivial factor of N. Move to step 5.
        | in the case of 15, [3, 5, 6, 9, 10, 12] do not satisfy this requirement
    4. If gdc(a, N) == 1, use the period-finding subroutine to find r, the period of the following function
        | f(x) = a^x mod N
        | the smallest integer r for which f(x + r) = f(x)
        | phase estimation is used to find period
    5. If r is odd, go back to step 1.
    6. If a^(r/2) == -1 (mod N), go back to step 1.
    7. The factors of N are gcd(a^(r/2) +- 1, N).

    [Quantum Step]
    1. Start with a pair of input and output qubit registers with log2(N) qubits each.
        | Initialized to N^(-1/2) ∑{subscript}x|x>|0> where x runs from 0 to N-1
    2. Construct f(x) as a quantum function and apply it to the above state.
        | N^(-1/2) ∑{subscript}x|x>|f(x)>
    3. Apply the quantum Fourier transform on the input register. 
        | 
    """
    a = 8
    N = 15

    FACTOR_FOUND = False
    ATTEMPT = 0
    """
    Shor's algorithm is probabilistic and thus runs until a non-trivial factor for N is found.
    """
    while not FACTOR_FOUND:
        ATTEMPT += 1
        print(f"\nAttempt {ATTEMPT}")

        phase = phase_estimation(
            c_amod15(a, N),
            psi_prep,
            precision=8
        )
        frac = Fraction(phase).limit_denominator(N)
        r = frac.denominator

        print(f"r = {r}")

        if phase != 0:
            # Guess for a factor is gcd(x^{r/2} +- 1 , 15)
            guesses = [gcd(a ** (r // 2) - 1, N), gcd(a ** (r // 2) + 1, N)]
            print(f"Guessed factors: {guesses[0]} and {guesses[1]}")
            for guess in guesses:
                if guess not in [1, N] and (N % guess) == 0:
                    # Guess is a factor!
                    print(f"Non-trivial factor found: {guess}")
                    FACTOR_FOUND = True