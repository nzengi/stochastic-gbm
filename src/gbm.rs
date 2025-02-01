use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Normal};
use thiserror::Error;

/// Errors that can occur during GBM simulation
#[derive(Error, Debug)]
pub enum GBMError {
    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),
    #[error("Distribution error: {0}")]
    DistributionError(#[from] rand_distr::NormalError),
}

/// Configuration parameters for the GBM simulation
#[derive(Debug, Clone)]
pub struct GBMParams {
    /// Drift (expected return) of the asset
    pub mu: f64,
    /// Volatility (standard deviation of returns)
    pub sigma: f64,
    /// Initial price
    pub s_0: f64,
    /// Time horizon in years
    pub t_end: f64,
}

/// Simulation parameters
#[derive(Debug, Clone)]
pub struct SimulationParams {
    /// Number of paths to simulate
    pub n_paths: usize,
    /// Number of time steps per path
    pub n_steps: usize,
}

/// Result of a GBM simulation
#[derive(Debug)]
pub struct SimulationResult {
    /// Price paths [path_index][time_step]
    pub paths: Vec<Vec<f64>>,
    /// Time points corresponding to each step
    pub time_points: Vec<f64>,
    /// Parameters used for the simulation
    pub params: GBMParams,
}

/// Geometric Brownian Motion simulator
#[derive(Debug)]
pub struct GeometricBrownianMotion {
    params: GBMParams,
    sim_params: SimulationParams,
    rng: ThreadRng,
}

impl GeometricBrownianMotion {
    /// Creates a new GBM simulator with the given parameters
    ///
    /// # Arguments
    ///
    /// * `mu` - Drift (expected return)
    /// * `sigma` - Volatility
    /// * `n_paths` - Number of paths to simulate
    /// * `n_steps` - Number of time steps
    /// * `t_end` - Time horizon in years
    /// * `s_0` - Initial price
    ///
    /// # Returns
    ///
    /// A Result containing either the GBM simulator or an error
    pub fn new(
        mu: f64,
        sigma: f64,
        n_paths: usize,
        n_steps: usize,
        t_end: f64,
        s_0: f64,
    ) -> Result<Self, GBMError> {
        // Validate parameters
        if sigma <= 0.0 {
            return Err(GBMError::InvalidParameter(
                "Volatility must be positive".to_string(),
            ));
        }
        if t_end <= 0.0 {
            return Err(GBMError::InvalidParameter(
                "Time horizon must be positive".to_string(),
            ));
        }
        if s_0 <= 0.0 {
            return Err(GBMError::InvalidParameter(
                "Initial price must be positive".to_string(),
            ));
        }
        if n_paths == 0 || n_steps == 0 {
            return Err(GBMError::InvalidParameter(
                "Number of paths and steps must be positive".to_string(),
            ));
        }

        Ok(Self {
            params: GBMParams {
                mu,
                sigma,
                s_0,
                t_end,
            },
            sim_params: SimulationParams { n_paths, n_steps },
            rng: rand::thread_rng(),
        })
    }

    /// Simulates price paths using the Euler-Maruyama method
    ///
    /// # Returns
    ///
    /// A SimulationResult containing the price paths and simulation metadata
    pub fn simulate(&mut self) -> Result<SimulationResult, GBMError> {
        let dt = self.params.t_end / self.sim_params.n_steps as f64;
        let normal = Normal::new(0.0, 1.0)?;
        
        let mut paths = vec![vec![self.params.s_0; self.sim_params.n_steps + 1]; self.sim_params.n_paths];
        let time_points: Vec<f64> = (0..=self.sim_params.n_steps)
            .map(|i| i as f64 * dt)
            .collect();

        let drift = (self.params.mu - 0.5 * self.params.sigma.powi(2)) * dt;
        let vol_sqrt_dt = self.params.sigma * dt.sqrt();

        for i in 0..self.sim_params.n_paths {
            for j in 1..=self.sim_params.n_steps {
                let z = normal.sample(&mut self.rng);
                paths[i][j] = paths[i][j - 1] * (drift + vol_sqrt_dt * z).exp();
            }
        }

        Ok(SimulationResult {
            paths,
            time_points,
            params: self.params.clone(),
        })
    }

    /// Returns statistical properties of the simulated paths at a given time index
    ///
    /// # Arguments
    ///
    /// * `result` - The simulation result
    /// * `time_idx` - The time index to analyze
    ///
    /// # Returns
    ///
    /// A tuple containing (mean, standard_deviation, min, max)
    pub fn get_statistics(result: &SimulationResult, time_idx: usize) -> (f64, f64, f64, f64) {
        let prices: Vec<f64> = result.paths.iter().map(|path| path[time_idx]).collect();
        
        let mean = prices.iter().sum::<f64>() / prices.len() as f64;
        
        let variance = prices.iter()
            .map(|&price| (price - mean).powi(2))
            .sum::<f64>() / (prices.len() - 1) as f64;
        
        let std_dev = variance.sqrt();
        let min = prices.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max = prices.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        (mean, std_dev, min, max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gbm_initialization() {
        let result = GeometricBrownianMotion::new(0.05, 0.3, 1000, 252, 1.0, 100.0);
        assert!(result.is_ok());

        let result = GeometricBrownianMotion::new(0.05, -0.3, 1000, 252, 1.0, 100.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_simulation_basic_properties() {
        let mut gbm = GeometricBrownianMotion::new(0.05, 0.3, 1000, 252, 1.0, 100.0).unwrap();
        let result = gbm.simulate().unwrap();

        assert_eq!(result.paths.len(), 1000);
        assert_eq!(result.paths[0].len(), 253);
        assert_eq!(result.time_points.len(), 253);

        for path in &result.paths {
            assert_relative_eq!(path[0], 100.0);
        }
    }

    #[test]
    fn test_statistics() {
        let mut gbm = GeometricBrownianMotion::new(0.05, 0.3, 1000, 252, 1.0, 100.0).unwrap();
        let result = gbm.simulate().unwrap();
        
        let (mean, std_dev, min, max) = GeometricBrownianMotion::get_statistics(&result, 252);
        
        assert!(mean > 0.0);
        assert!(std_dev > 0.0);
        assert!(min < mean);
        assert!(max > mean);
    }
}
