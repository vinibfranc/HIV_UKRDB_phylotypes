data {
 int<lower = 1> num_groups;
 int<lower = num_groups> num_y;
 int<lower = 1, upper = num_groups> group[num_y];
 real y[num_y];
}

parameters {
	real<lower=2, upper=7> y_mean_pop;
 real<lower=0> y_sd;
 real<lower=0> y_sd_group;
 vector[num_groups] group_effects_unscaled;
}

transformed parameters {
 vector[num_groups] y_by_group = y_mean_pop +
 group_effects_unscaled * y_sd_group;
 real y_mean[num_y];
 for (i in 1:num_y) {
   y_mean[i] = y_by_group[group[i]];
 }
}

model {
	// Update priors
	y_mean_pop ~ normal(4.5, 1) T[2, 7]; // truncated normal because no values below 2 and above 7
	y_sd ~ normal(0.85, 0.05) T[0, ]; // based on histogram of posterior
	y_sd_group ~ normal(0.05, 0.3) T[0, ]; // based on histogram of posterior

 group_effects_unscaled ~ normal(0, 1);
 // Disable likelihood
 // y ~ normal(y_mean, y_sd);
}
