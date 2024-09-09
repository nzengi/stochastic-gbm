# Geometric Brownian Motion (GBM) in Rust

This Rust crate provides an implementation of the **Geometric Brownian Motion (GBM)** model, a common stochastic process used to simulate the price movement of financial assets over time.

## Overview

The Geometric Brownian Motion (GBM) is defined by the stochastic differential equation:

dS = μ * S * dt + σ * S * dW

```swfit
Where:
- `S` is the asset price,
- `μ` is the drift (mean or trend),
- `σ` is the volatility (standard deviation of returns),
- `dW` is the Wiener process increment (Brownian motion),
- `dt` is the time increment.
```

GBM is commonly used in finance for modeling asset prices as it enforces the non-negativity of the asset price.

## Installation

Add the following to your `Cargo.toml` to include this crate in your project:

```toml
[dependencies]
stochastic_gbm = "0.1.0"
```

## Usage
Below is a simple example of how to use the Geometric Brownian Motion (GBM) in Rust:

use stochastic_gbm::GeometricBrownianMotion;

```rust
fn main() {
    // Create a new GBM model with the following parameters:
    // - Drift (mu): 0.2
    // - Volatility (sigma): 0.4
    // - Number of paths: 50
    // - Number of steps: 200
    // - Time horizon: 1.0
    // - Initial asset value (S_0): 500.0
    let gbm = GeometricBrownianMotion::new(0.2, 0.4, 50, 200, 1.0, 500.0);

    // Simulate the paths
    let paths = gbm.simulate();

    // Print the simulated paths
    for path in paths {
        println!("{:?}", path);
    }
}
```

## Parameters
When creating an instance of GeometricBrownianMotion, the following parameters are required:

mu (f64): The drift, representing the expected return or trend of the asset's price.
sigma (f64): The volatility, representing the variability or standard deviation of the returns.
n_paths (usize): The number of independent simulated paths.
n_steps (usize): The number of time steps in each path.
t_end (f64): The total time duration of the simulation.
s_0 (f64): The initial value of the asset at time t = 0.

## Output
The simulate function returns a 2D vector (Vec<Vec<f64>>) where each inner vector represents a simulated price path. Each path contains n_steps + 1 values, including the initial value S_0.

## Example Simulation
In this example, we simulate 50 paths of an asset's price over a time horizon of 1 year, with 200 time steps. The drift (expected return) is set to 20% and the volatility is set to 40%.



