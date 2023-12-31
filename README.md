# QSEA

There are two programs, **qsea_ON_finiteT** for the $O(N)$ symmetric model at finite temperature, and **qsea_N2_finiteT** for a $N=2$ model with two fields.

The programs calculate the potential using QSEA for different temperatures, and is multithreaded.

## How to run
For both, edit the **qsea-params.json** file to change the parameters of the run.

### O(N) symmetric model
An example of a valid parameter file is:

```json
{
    "N": 1, #number of fields
    "T": [0, 0.4, 0.8, 1.2, 1.6, 2.0], #a list of temperatures to calculate the potential for
    "v": 1, #parameter for the potential (explained below)
    "Gamma": 0.5, #parameter for the potential (explained below)
    "params":{
        "kStart": 4, #starting k for the flow
        "phiSteps": 10000, #number of steps in the field
        "phiRange": [-10, 10], #range of the field
        "logAbsErr": -8, #log_10 of the permitted absolute error for the RK45 algorithm
        "logRelErr": -6 #log_10 of the permitted relative error for the RK45 algorithm
    },
    "outputFileName": "data" #name of the output file
}
```

Where the potenial is defined as

$$
V(\phi) = \Gamma^4 \left(\left(\frac{\phi}{v}\right)^2 -\frac{1}{2}\left(\frac{\phi}{v}\right)^4\right)
$$

### General N=2 model
An example of a valid parameter file is:

```json
{
    "N": 2,
    "T": [0, 0.1, 0.275, 0.45, 0.625, 0.8],
    "m12": 0.014, #mass squared
    "m22": 0.014, #mass squared
    "lambda1": 0.1775,
    "lambda2": 0.185,
    "lambda12": 0.1,
    "alpha": -0.09,
    "params":{
        "kStart": 2,
        "phi1Steps": 1000,
        "phi2Steps": 100,
        "phi1Range": [-15, 15],
        "phi2Range": [-10, 10],
        "logAbsErr": -7,
        "logRelErr": -6
    },
    "outputFileName": "data"
}
```

Where the potenial is defined as

$$
V(\phi,\varphi) = \frac{1}{2}\left(m_1^2 \phi^2 + m_2^2 \varphi^2\right) + \frac{1}{3!}\alpha \phi^3 + \frac{1}{4!}\left(\lambda_1^2 \phi^2 + \lambda_2^2 \varphi^2\right) + \frac{1}{4} \phi^2 \varphi^2
$$

## Results
Result are outputted as follows:

```json
{
    "phi": [], #list of phi (i.e. x axis)
    "runs": [ #can be a list of any length
        {
            "T": 0.0, #each run has a defined temperature
            "pot": [] #and a list for the potential at that temperature (i.e. y axis)
        }
    ],
    "treeLevel": [] #list of the tree level potential (i.e. y axis)
}
```

The results save the data for the potential V(\phi, \varphi = 0) for the multiple fields case.

## Plotting
The python files **plot_O1.py** and **plot_O2.py** plot specificially for the $O(1)$ and $O(2)$ cases, but can be generalised for $N \geq 3$.

The file **plot.py** in the **qsea_N2_finiteT** plotting folder plots only for the $N=2$ case.

They require the [cosmoTransitions](https://github.com/clwainwright/CosmoTransitions) package to calculate the perturbative potentials for comparison.

## Recompiling the C++ files
If you wish to recompile the files, you can use the _make_ command. It requires two dependencies - [openMP](https://www.openmp.org/) and [nlohmann/JSON](https://github.com/nlohmann/json).
