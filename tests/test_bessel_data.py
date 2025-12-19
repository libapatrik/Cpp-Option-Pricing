import numpy as np
import scipy.special as sp
import pandas as pd

"""
Goal: Ensure that the modifiedBesselI (C++, in Utils) produces the correct numbers
Check: Compare the values from modifiedBesselI against the SciPy's Bessel function

"""

# ! Needs to be fixed properly
def generate_test_data(filename="bessel_test_data.csv"):
    # 1. Define Test Scenarios
    # Real orders (nu)
    nus = [0.0, 0.5, 1.0, 2.5, 5.0]

    # Complex arguments (z)
    # Grid covering small, transition, and large magnitudes
    real_parts = np.concatenate([
        np.linspace(-5, 5, 11),       # Small
        np.linspace(-35, 35, 15)      # Large/Transition
    ])
    imag_parts = np.concatenate([
        np.linspace(-5, 5, 11),
        np.linspace(-35, 35, 15)
    ])

    data = []

    for nu in nus:
        for re in real_parts:
            for im in imag_parts:
                z = complex(re, im)

                # Calculate reference value using Scipy
                # scipy.special.iv handles complex arguments robustly
                expected = sp.iv(nu, z)

                # Check for NaNs or Inf (skip unstable points if any)
                if np.isfinite(expected):
                    data.append({
                        "nu": nu,
                        "z_real": re,
                        "z_imag": im,
                        "expected_real": expected.real,
                        "expected_imag": expected.imag
                    })

    # 2. Save to CSV
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False, float_format='%.18e')
    print(f"Generated {len(df)} test cases in {filename}")

if __name__ == "__main__":
    generate_test_data()