# High-dimensional linear models

The files are structured in the same way as the linear model experiments (found in `../lm/`), but with an added screening step using [Sure Independent Screening](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2008.00674.x).

We separate the code for MIO because of its extra dependencies (`get_outputs_hd_mio.R`, `lmhdselect-mio.R`). For the rest of the methods, to obtain outputs for a specific experimental setting (Setting 1, 2, or 3), run the corresponding `lmhdelect-s*.R` file, with options for parameters such as $n, p, \Sigma$. The usage options for these files are similar, see the options for `lmhdselect-s1.R` below.

```
Usage: lmhdselect-s1.R [options]


Options:
        -n DOUBLE, --numsamp=DOUBLE
                sample size [default = 100]

        -p DOUBLE, --params=DOUBLE
                params size [default = 500]

        -b DOUBLE, --bootsamp=DOUBLE
                resample size [default = 1000]

        -s DOUBLE, --sigma=DOUBLE
                error sd [default = 1]

        -e DOUBLE, --equal=DOUBLE
                equal coef indicator [default = 1]

        -h, --help
                Show this help message and exit
```

Outputs from the `lmhdselect` files give performance metrics for the e-values and other methods. The file `get_outputs_hd.R` gives actual fitting functions to obtain these metrics, while `lmhdselect-s*.R` provide parametrized wrappers around them.