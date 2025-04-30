from scipy.optimize import minimize
import numpy as np

# Define the inequality constraint function
def ineq_constraint(x):
    r1, r2, r3, b1, b2, b3, g1, g2, g3 = x  # g3 is fixed to 1

    lhs = r1 * (b2 * g3 + b3 * g2)
    rhs = r2 * (b1 * g3 + b3 * g1)
    return rhs - lhs  # want lhs < rhs → rhs - lhs > 0

# Define the equality constraint functions
def eq1(x): return x[0] + x[3] + x[6] - 1  # r1 + b1 + g1 = 1
def eq2(x): return x[1] + x[4] + x[7] - 1  # r2 + b2 + g2 = 1
def eq3(x): return x[2] + x[5] + x[8] - 1     # r3 + b3 + g3 = 4
def order1(x): return x[0] - x[1]          # r1 >= r2 → r1 - r2 ≥ 0
def order2(x): return x[1] - x[2]          # r2 >= r3 → r2 - r3 ≥ 0

# Initial guess
x0 = np.array([0.5, 0.4, 0.3, 0.2, 0.2, 2.7, 0.3, 0.4, 0.0])

# Bounds (if needed)
bounds = [(0, None)] * 9  # r1 to g2 ≥ 0, unbounded above for now

# Constraints list
constraints = [
    {'type': 'ineq', 'fun': ineq_constraint},
    {'type': 'eq', 'fun': eq1},
    {'type': 'eq', 'fun': eq2},
    {'type': 'eq', 'fun': eq3},
    {'type': 'ineq', 'fun': order1},
    {'type': 'ineq', 'fun': order2}
]

# Dummy objective function (we just want feasibility)
result = minimize(lambda x: 0, x0, method='SLSQP', bounds=bounds, constraints=constraints)

if result.success:
    r1, r2, r3, b1, b2, b3, g1, g2, g3 = result.x
    print("Feasible solution found:")
    print(f"r1={r1:.4f}, r2={r2:.4f}, r3={r3:.4f}")
    print(f"b1={b1:.4f}, b2={b2:.4f}, b3={b3:.4f}")
    print(f"g1={g1:.4f}, g2={g2:.4f}, g3={g3:.4f}")
else:
    print("No feasible solution found.")