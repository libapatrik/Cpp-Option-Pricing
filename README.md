# Option Pricing Library

## Table of Contents

1. [Black-Scholes](#black-scholes)
   - [Analytical Pricing](#bs-analytical-pricing)
   - [Greeks](#bs-greeks)

2. [Dupire Local Volatility](#dupire-local-volatility)
   - [Local Volatility Extraction](#local-volatility-extraction)
   - [Finite Difference Pricing](#fd-pricing)
   

3. [Heston Stochastic Volatility](#heston-stochastic-volatility)
   - [Monte Carlo Schemes](#heston-mc-schemes)

4. [Heston Stochastic Local Volatility](#heston-stochastic-local-volatility)

5. [Market Data](#market-data)
   - [Discount Curves](#discount-curves)
   - [Volatility Surfaces](#volatility-surfaces)
   - [Interpolation Schemes](#interpolation-schemes)

6. [PDE Framework](#pde-framework)
   - [Grid Construction](#grid-construction)
   - [Theta-Method Solvers](#theta-method-solvers)
   - [Boundary Conditions](#boundary-conditions)
   - [Greeks via Finite Differences](#fd-greeks)

7. [Utilities](#utilities)
   - [Thomas Algorithm](#thomas-algorithm)
   - [Numerical Derivatives](#numerical-derivatives)
   - [Newton Method](#newton)
   - [COS Method](#cos)
   - [CIR Sampling](#cir-sampling)
   - [Characteristic Function of Integrated Variance Process](#chf-int-var)

8. [Building and Testing](#building-and-testing)
   - [Test Suite](#test-suite)