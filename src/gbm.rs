use rand::Rng;

/// The Geometric Brownian Motion (GBM) model simulates the price movement
/// of an asset over time using the following formula:
///
/// dS = μ * S * dt + σ * S * dW
///
/// Where:
/// - `S` is the asset price,
/// - `μ` is the drift (expected return),
/// - `σ` is the volatility (standard deviation of returns),
/// - `dW` is the Wiener process increment (Brownian motion),
/// - `dt` is the time increment.
pub struct GeometricBrownianMotion {
    pub mu: f64,
    pub sigma: f64,
    pub n_paths: usize,
    pub n_steps: usize,
    pub t_end: f64,
    pub s_0: f64,
}

impl GeometricBrownianMotion {
    /// Creates a new instance of the Geometric Brownian Motion model.
    ///
    /// # Arguments
    ///
    /// * `mu` - The drift (mean or expected return) of the asset's price.
    /// * `sigma` - The volatility (standard deviation of returns) of the asset.
    /// * `n_paths` - Number of simulated paths.
    /// * `n_steps` - Number of steps in each path.
    /// * `t_end` - Total time duration of the simulation.
    /// * `s_0` - Initial value of the asset (price at t=0).
    ///
    /// # Returns
    ///
    /// A new instance of `GeometricBrownianMotion`.
    pub fn new(mu: f64, sigma: f64, n_paths: usize, n_steps: usize, t_end: f64, s_0: f64) -> Self {
        Self {
            mu,
            sigma,
            n_paths,
            n_steps,
            t_end,
            s_0,
        }
    }

    /// Simulates the asset price paths using the Euler-Maruyama method for Geometric Brownian Motion.
    ///
    /// # Returns
    ///
    /// A 2D vector where each inner vector represents a simulated path of asset prices.
    ///
    /// Each path has `n_steps + 1` values, including the initial value `s_0`.
    pub fn simulate(&self) -> Vec<Vec<f64>> {
        let dt = self.t_end / self.n_steps as f64; // Time step size
        let mut rng = rand::thread_rng(); // Random number generator
        let mut paths = vec![vec![self.s_0; self.n_steps + 1]; self.n_paths]; // Initialize paths

        // Simulate each path
        for i in 0..self.n_paths {
            for j in 1..=self.n_steps {
                let d_w = rng.gen::<f64>() * dt.sqrt(); // Wiener process increment
                paths[i][j] = paths[i][j - 1] * (1.0 + self.mu * dt + self.sigma * d_w).exp(); // Euler-Maruyama update
            }
        }

        paths
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gbm_simulation() {
        let gbm = GeometricBrownianMotion::new(0.2, 0.4, 50, 200, 1.0, 500.0);
        let paths = gbm.simulate();
        assert_eq!(paths.len(), 50);
        assert_eq!(paths[0].len(), 201); // n_steps + 1
    }
}
