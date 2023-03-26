# Gaussian process poll aggregator for Finnish 2023 parlamentary elections
Stan code for fitting a Gaussian process regression model to Finnish
parlamentary election polling data from after the 2019 elections until
just before the 2023 elections.

## Dependencies
- R (tested with v4.1.3).
- rstan (tested with v2.21.8).

## Structure
- `data/`: Polling data from 2019 to 2023.
- `src/`: Stan code for the model and R scripts for data processing & plotting.

## Quickstart
Clone the repository and run
```
Rscript fit_model.R
Rscript plot_fit.R
```

This will create the `gp_poll_aggregator_fit.pdf` file which will show the
results of the Gaussian process model defined in
`src/gp_poll_aggregator.stan` fit to
`data/polling_data_2019-2023.tsv`.

## License
Licensed under the [BSD-3-Clause license](https://opensource.org/licenses/BSD-3-Clause). A copy of the license is supplied with the project, or can alternatively be obtained from [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
