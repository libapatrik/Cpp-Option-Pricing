import numpy as np
import scipy.special as sp
import pandas as pd
import os


"""
Goal: Ensure that the modifiedBesselI (C++, in Utils) produces the correct numbers
Check: Compare the values from modifiedBesselI against the SciPy's Bessel function

"""

def generate_test_data(filename="bessel_test_data.csv"):
    # 1. Define Test Scenarios
    # Real orders (nu)
    nus = [0.0, 0.5, 1.0, 2.5, 5.0]

    # Complex arguments (z)
    # Use unique points covering small, medium, and larger magnitudes
    # Avoid extreme values that cause overflow in Bessel functions
    real_parts = np.unique(np.concatenate([
        np.linspace(-5, 5, 11),        # Small range
        np.linspace(-20, 20, 9),       # Medium range (avoid ±35 overflow)
    ]))
    imag_parts = np.unique(np.concatenate([
        np.linspace(-5, 5, 11),
        np.linspace(-20, 20, 9),
    ]))

    data = []

    for nu in nus:
        for re in real_parts:
            for im in imag_parts:
                z = complex(re, im)

                # Calculate reference value using Scipy
                # scipy.special.iv handles complex arguments robustly
                expected = sp.iv(nu, z)

                # Check for NaNs or Inf (skip unstable points if any)
                # np.isfinite works on complex by checking both real and imag parts
                if np.isfinite(expected.real) and np.isfinite(expected.imag):
                    data.append({
                        "nu": nu,
                        "z_real": re,
                        "z_imag": im,
                        "expected_real": expected.real,
                        "expected_imag": expected.imag
                    })

    # 2. Save to CSV in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, filename)
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False, float_format='%.18e')
    print(f"Generated {len(df)} test cases in {filepath}")

if __name__ == "__main__":
    generate_test_data()
