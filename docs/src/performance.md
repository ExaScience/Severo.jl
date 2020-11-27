# Performance

|    Function   | `Cell.jl` | Seurat | Speedup |
|:--------------|:----------|:------:|:-------:|
| filter_counts |  0.066542 |  0.434 | 6.71x   |
| normalize     |  0.06     |  0.212 | 3.53x   |
| highly_var..  |  0.362585 |  1.114 | 3.07x   |
