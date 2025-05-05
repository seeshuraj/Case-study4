# Multigrid Solver for 2D Poisson Equation

This project implements a **Multigrid (MG) method** to solve the 2D Poisson equation using a recursive V-cycle with weighted Jacobi smoothing. The convergence behavior is visualized and benchmarked across different grid sizes.

---

## 🔧 Features

- V-cycle multigrid implementation with adjustable levels
- Weighted Jacobi smoother
- Full-weighting restriction and bilinear prolongation
- Convergence plots and benchmark comparisons
- CSV output of cycles, residuals, and timings
- Ready-to-run and extensible Python code

---

## 📁 Files

| File | Description |
|------|-------------|
| `mg_solver.py` | Main Python script to run the multigrid solver |
| `mg_convergence.png` | Log plot showing convergence over V-cycles |
| `mg_comparison_results.csv` | Benchmark results for various grid sizes |
| `multigrid_assignment_report.md` | Detailed write-up and explanation |
| `requirements.txt` | Required Python packages |

---

## ▶️ How to Run

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Run the solver:
   ```bash
   python3 mg_solver.py
   ```

   It will output residuals per cycle and generate `mg_convergence.png`.

---

## 🧪 Benchmarking

The script runs V-cycle multigrid for various grid sizes (`N = 16, 32, 64, 128, 256`) using:
- 2-Level MG
- Max-level MG (deepest possible)

Results are saved in:
```
mg_comparison_results.csv
```

---

## 📊 Sample Output

![Convergence Plot](mg_convergence.png)

---

## 📄 Report

See `multigrid_assignment_report.md` for the full write-up.
---

## 👤 Author

**Seeshuraj Bhoopalan**

---
