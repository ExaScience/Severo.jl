# Severo.jl

A software package for analysis and exploration of single-cell RNA-seq datasets

## Introduction

Severo.jl was created for the scalable analysis of single-cell RNA-seq datasets: the implementation is designed
to minimize hardware resources such as memory while maximizing performance through parallelism and hardware optimizations.

The package provides a toolbox of different algorithms and statistical methods from which the user can pick and choose.

Please see the [documentation](https://exascience.github.io/Severo.jl/dev/) for instructions and examples.

## Prerequisites

The package depends on `Severo_jll` which is available from https://github.com/Exascience/Severo_jll.jl.

To install `Severo_jll`, you'll need to run
```julia
import Pkg; Pkg.add("https://github.com/Exascience/Severo_jll.jl")
```

## Contributing

Contributions are encouraged. If there are additional features you would like to use, please open an [issue](https://github.com/Exascience/Severo.jl/issues) or [pull request](https://github.com/Exascience/Severo.jl/pulls).

Additional examples and documentation improvements are also very welcome.
