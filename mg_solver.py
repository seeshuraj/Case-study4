import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve


def f_source(x, y):
    return 2 * (np.pi ** 2) * np.sin(np.pi * x) * np.sin(np.pi * y)


def initialize_problem(N):
    h = 1.0 / (N + 1)
    x = np.linspace(h, 1 - h, N)
    y = np.linspace(h, 1 - h, N)
    X, Y = np.meshgrid(x, y)
    b = f_source(X, Y).reshape(-1)
    return h, b


def generate_laplacian_matrix(N, h):
    A = lil_matrix((N * N, N * N))

    for i in range(N):
        for j in range(N):
            idx = i * N + j
            A[idx, idx] = 4.0
            if i > 0:
                A[idx, idx - N] = -1.0
            if i < N - 1:
                A[idx, idx + N] = -1.0
            if j > 0:
                A[idx, idx - 1] = -1.0
            if j < N - 1:
                A[idx, idx + 1] = -1.0

    return (A / (h * h)).tocsr()


def weighted_jacobi(A, x, b, omega, num_iter):
    D_inv = 1.0 / A.diagonal()
    for _ in range(num_iter):
        r = b - A.dot(x)
        x = x + omega * (D_inv * r)
    return x


def restrict(residual, N):
    Nc = N // 2
    res_2d = residual.reshape(N, N)
    restricted = np.zeros((Nc, Nc))

    for i in range(Nc):
        for j in range(Nc):
            i2, j2 = 2 * i, 2 * j
            restricted[i, j] = (res_2d[i2, j2] +
                                0.5 * (res_2d[i2+1, j2] + res_2d[i2, j2+1]) +
                                0.25 * res_2d[i2+1, j2+1]) / 4.0
    return restricted.reshape(-1)


def prolong(correction, N):
    Nc = N // 2
    correction_2d = correction.reshape(Nc, Nc)
    fine = np.zeros((N, N))

    for i in range(Nc):
        for j in range(Nc):
            fine[2*i, 2*j] = correction_2d[i, j]

    for i in range(1, N-1, 2):
        for j in range(N):
            fine[i, j] = 0.5 * (fine[i-1, j] + fine[i+1, j])
    for i in range(N):
        for j in range(1, N-1, 2):
            fine[i, j] = 0.5 * (fine[i, j-1] + fine[i, j+1])

    return fine.reshape(-1)


def v_cycle(A_list, x, b, level, lmax, omega, nu):
    A = A_list[level]
    x = weighted_jacobi(A, x, b, omega, nu)
    r = b - A.dot(x)

    if level == lmax:
        x = spsolve(A, b)
        return x

    N = int(np.sqrt(len(r)))
    r_coarse = restrict(r, N)

    x_coarse = np.zeros_like(r_coarse)
    x_coarse = v_cycle(A_list, x_coarse, r_coarse, level + 1, lmax, omega, nu)

    correction = prolong(x_coarse, N)
    x += correction

    x = weighted_jacobi(A, x, b, omega, nu)
    return x


def build_multigrid_levels(N0, lmax):
    A_levels = []
    N = N0
    for level in range(lmax + 1):
        h = 1.0 / (N + 1)
        A = generate_laplacian_matrix(N, h)
        A_levels.append(A)
        N = N // 2
    return A_levels


def run_multigrid_solver(N, lmax, max_cycles=50, tol=1e-7, omega=2/3, nu=3):
    h, b = initialize_problem(N)
    A_levels = build_multigrid_levels(N, lmax)
    A0 = A_levels[0]

    x = np.zeros_like(b)
    residuals = []

    start_time = time.time()

    for cycle in range(max_cycles):
        x = v_cycle(A_levels, x, b, 0, lmax, omega, nu)
        res = np.linalg.norm(b - A0.dot(x))
        residuals.append(res)
        print(f"Cycle {cycle+1}: Residual = {res:.2e}")

        if res < tol:
            print("Converged.")
            break
        if res > 1e10:
            print("Diverging! Stopping...")
            break

    total_time = time.time() - start_time
    print(f"Total time: {total_time:.2f} sec, Cycles: {len(residuals)}")

    return x, residuals, total_time


def benchmark_mg_convergence():
    Ns = [16, 32, 64, 128, 256]
    results = []

    for N in Ns:
        print(f"\nRunning for N = {N} with 2-level MG...")
        lmax_2 = 1
        _, res_2, time_2 = run_multigrid_solver(N=N, lmax=lmax_2)
        min_res_2 = res_2[-1]
        cycles_2 = len(res_2)

        print(f"Running for N = {N} with max-level MG...")
        lmax_max = int(np.log2(N // 8))
        _, res_m, time_m = run_multigrid_solver(N=N, lmax=lmax_max)
        min_res_m = res_m[-1]
        cycles_m = len(res_m)

        results.append({
            "N": N,
            "2L_Cycles": cycles_2, "2L_MinRes": min_res_2, "2L_Time": time_2,
            "MaxL_Cycles": cycles_m, "MaxL_MinRes": min_res_m, "MaxL_Time": time_m
        })

    import csv
    with open("mg_comparison_results.csv", "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print("\nBenchmarking complete. Results saved to mg_comparison_results.csv.")


if __name__ == "__main__":
    N = 128
    lmax = 3
    solution, res_history, runtime = run_multigrid_solver(N=N, lmax=lmax)

    plt.figure()
    plt.semilogy(range(1, len(res_history)+1), res_history, marker='o')
    plt.xlabel("V-cycle")
    plt.ylabel("Residual norm (log scale)")
    plt.title(f"Multigrid Convergence (N={N}, lmax={lmax})")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("mg_convergence.png")
    plt.show()

    # Uncomment to run full benchmark
    benchmark_mg_convergence()
