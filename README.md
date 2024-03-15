
<!-- README.md is generated from README.Rmd. Please edit that file -->

# twdtwr

The `twdtwr` package provides an implementation of the Time-Weighted
Dynamic Time Warping (TWDTW) algorithm in R.

> **Note**: This is a project I created to learn more about the twdtw
> algorithm.

## Installation

To install the `twdtwr` package, you can use the `devtools` package to
install directly from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("m3nin0-labs/twdtwr")
```

## Usage

The `twdtwr` package simplifies the application of TWDTW analysis to
pairs of time-series data. The main function provided by the package is
`twdtw`, which compares two time-series datasets based on their temporal
dynamics.

``` r
library(twdtwr)

# Assuming time_series_a, time_series_b, dates_a, and dates_b are predefined
result <- twdtwr::twdtw(time_series_a, time_series_b, dates_a, dates_b)
```

## Reference

The implementation of TWDTW in `twdtwr` is based on the following
publication:

Maus, Victor, Gilberto Camara, Ricardo Cartaxo, Alber Sanchez, Fernando
M. Ramos, and Gilberto R. de Queiroz. 2016. “A Time-Weighted Dynamic
Time Warping Method for Land-Use and Land-Cover Mapping.” IEEE Journal
of Selected Topics in Applied Earth Observations and Remote Sensing 9
(8): 3729–39. <https://doi.org/10.1109/JSTARS.2016.2517118>.

## Contributing

Contributions to `twdtwr` are welcome, including bug reports, feature
requests, and pull requests. Please see the GitHub repository for more
details on how to contribute.

## License

The `twdtwr` package is open source and licensed under the MIT License.
