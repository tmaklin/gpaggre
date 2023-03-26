// BSD 3-Clause License
// 
// Copyright (c) 2023, Tommi Mäklin (tommi `at' maklin.fi)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// Multiple output Gaussian process model for compositional data
// observed over several time points.

functions {
    // Generalized inverse Gaussian prior for the kernel lengthscale,
    // this penalizes both near-zero and very large values. See
    // section 10.3.1 in the Stan documentation for details:
    // https://mc-stan.org/docs/stan-users-guide/fit-gp.html
    real generalized_inverse_gaussian_lpdf(real x, int p,
					   real a, real b) {
	// Generalized inverse Gaussian lpdf
	return p * 0.5 * log(a / b)
	    - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
	    + (p - 1) * log(x)
	    - (a * x + b / x) * 0.5;
    }
}

data {
    // Gallup data
    int<lower=1> N_obs; // N observed polls.
    int<lower=1> P_obs; // P political parties.
    vector<lower=0>[N_obs] time_from_start_obs; // Time from first observation.
    matrix<lower=0, upper=100>[N_obs, P_obs] party_support; // Party support in percentages (0.0 to 100.0).

    // Pollster info
    int<lower=1> N_pollsters; // Total number of polling companies.
    int<lower=1, upper=N_pollsters> pollsters[N_obs]; // Integer vector indicating the pollling company.

    // Predicted points
    int<lower=0> N_pred; // How many points to predict.
    vector<lower=0>[N_pred] time_from_start_pred; // Which time points to predict values for.
}

transformed data {
    // Merge the observed and predicted values
    int<lower=0> N = N_obs + N_pred; // Covariance matrix should have this size
    vector<lower=0>[N] time_from_start; // 
    
    // Noise to keep the covariance matrices positive (semi) definite.
    real delta = 1e-9;

    // Data is constrained to sum up to 100%, unbound it to (-inf, inf)
    // via inverse additive logistic transformations
    int<lower=0> P = P_obs - 1; // Constraint means we have P-1 unbounded data points
    matrix[N_obs, P] unbounded_support;
    for (n in 1:N_obs) {
	// Unbound constrained support by taking elementwise log of P_obs - 1 elements
	// and subtracting log of the last element from each value.
	for (p in 1:P) {
	    unbounded_support[n, p] = log(party_support[n, p]/100.0) - log(party_support[n, P_obs]/100.0);
	}
    }

    // Merge observed and predicted values
    for (n1 in 1:N_obs) time_from_start[n1] = time_from_start_obs[n1];
    for (n2 in 1:N_pred) time_from_start[N_obs + n2] = time_from_start_pred[n2];
}

parameters {
    // GP covariance function parameters
    real<lower=0> time_magnitude; // Controls how large the covariances are in general.
    real<lower=0> time_lengthscale; // Controls how the covariance between observations decays in distance.

    // GP covariance function form
    // Inferring `ratio` means we also infer the form of the covariance function.
    // With ratio -> inf the covariance function approaches the squared exponential covariance.
    real<lower=0> ratio;

    // Polling error (each pollster is different) (residual noise).
    real<lower=0> sigma[N_pollsters];

    // GP should have mean zero and any deviations in the data are inferred by this parameter.
    real bias;

    // Used to sample from the GP, see https://mc-stan.org/docs/stan-users-guide/fit-gp.html
    cholesky_factor_corr[P] L_omega;
    matrix[N, P] eta;
}

transformed parameters {
    // Explicitly code for the latent representation of the GP.
    matrix[N, P] f; // GP values at all evaluated points

    {
	// Create the covariance matrix
	matrix[N, N] K_f;
	matrix[N, N] L_k; // Cholesky factorization of K_f

	for (n1 in 1:N) {
	    for (n2 in (n1 + 1):N) {
		// Use the rational quadratic covariance function. In
		// a Bayesian context this has the benefit that we can
		// infer a wide range of different covariance function forms.
		K_f[n1, n2] = exp(-ratio*log1p(((time_from_start[n1] - time_from_start[n2])^2/(2*ratio*time_lengthscale^2))));
		K_f[n2, n1] = K_f[n1, n2]; // Covariance matrices are symmetric
	    }
	    // Add `delta` to the diagonal to keep the matrix positive (semi) definite
	    K_f[n1, n1] = 1.0 + delta;
	}

	// Sample the latent GP `f` with an affine transformation of `eta` and `L_omega¬
	// see the multi output GP documentation in the Stan manual for details:
	// https://mc-stan.org/docs/stan-users-guide/fit-gp.html
	L_k = cholesky_decompose(K_f);
	f = L_k * eta * diag_pre_multiply(rep_vector(time_magnitude, P), L_omega);
    }
}

model {
    // Generalized inverse Gaussian prior for lengthscale and the ratio.
    //
    // This prior distribution allows for reasonably large values while
    // penalizing both near-zero and extremely large lengthscales.
    target += generalized_inverse_gaussian_lpdf(time_lengthscale | -1, 1.0/10.0, 10.0);
    target += generalized_inverse_gaussian_lpdf(ratio | -1, 1.0/10.0, 10.0);

    // Half Student-t prior for magnitudes
    time_magnitude ~ student_t(4, 0, 1);

    // Inverse gamma distribution on the polling errors
    // This format assumes they are reasonably small but not zero.
    for (i in 1:N_pollsters) {
	sigma[i] ~ inv_gamma(1, 1);
    }

    // lkj_corr_cholesky is a distribution on the Cholesky factorization of a covariance matrix
    L_omega ~ lkj_corr_cholesky(3);

    // Draw standard normal distributed variables for sampling `f`
    to_vector(eta) ~ std_normal();

    // Assume the classes are balanced by default by placing a standard normal prior
    bias ~ std_normal();

    // Assume the actual observations follow a multivariate logistic-normal distribution.
    for (i in 1:N_obs) {
	unbounded_support[i, ] ~ multi_normal(bias + f[i, ], diag_matrix(rep_vector(sigma[pollsters[i]], P)));
    }
}

generated quantities {
    // Predict unbounded support on the requested time points
    matrix[N_pred, P] pred;

    // Assume variance in unobserved points is the sum of the pollster variances
    real sigma_sum = 0.0;
    for (i in 1:N_pollsters) {
	sigma_sum += sigma[i];
    }

    // Predict points by sampling from multi_normal distribution with mean `bias + f` and diagonal covariance.
    for (n2 in 1:N_pred) {
	pred[n2, ] = to_row_vector(multi_normal_rng(bias + f[N_obs + n2, ], diag_matrix(rep_vector(sigma_sum, P))));
    }
}
