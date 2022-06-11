# Linear models

We separate the code for MIO because of its extra dependencies (`get_outputs_mio.R`, `lmselect-mio.R`). For the rest of the methods, to obtain outputs for a specific experimental setting (Setting 1, 2, or 3), run the corresponding `lmselect-s*.R` file, with options for parameters such as $n, p, \Sigma$. The usage options for these files are similar, see the options for `lmselect-s1.R` below.

```
Usage: lmselect-s1.R [options]


Options:
        -n DOUBLE, --numsamp=DOUBLE
                sample size [default = 500]   

        -p DOUBLE, --params=DOUBLE
                params size [default = 100]   

        -b DOUBLE, --bootsamp=DOUBLE
                resample size [default = 1000]

        -s DOUBLE, --sigma=DOUBLE
                error sd [default = 1]        

        -e DOUBLE, --equal=DOUBLE
                equal coef indicator [default = 1]

        -h, --help
                Show this help message and exit
```

Outputs from the `lmselect` files give performance metrics for the e-values and other methods. The file `get_outputs.R` gives actual fitting functions to obtain these metrics, while `lmselect-s*.R` provide parametrized wrappers around them.