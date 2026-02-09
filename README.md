# auxfinder
Finding auxiliary variables for imputation models

`auxfinder` finds auxiliary variables for imputation models. For details and more examples, see the
working paper here: [Zenodo Working Paper](https://doi.org/10.5281/zenodo.18536805)

To install `auxfinder`, type

    . net install auxfinder, replace from(https://raw.githubusercontent.com/fbittmann/auxfinder/main/)

Minimal example

    . webuse nhanes2, clear
    . replace hlthstat = . if hlthstat == 8
    . auxfinder hlthstat albumin vitaminc copper hdresult, testvars(sampl - lead) seed(123)

---

Main changes:

    16jan2026
    - first released on Github


