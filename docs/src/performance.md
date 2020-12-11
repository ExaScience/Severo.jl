# Performance

## PBMC 3K

|    Function   | `Cell.jl` | Seurat | Speedup |
|:--------------|:----------|:------:|:-------:|
| filter_counts |  0.066542 |  0.434 | 6.71x   |
| normalize     |  0.06     |  0.212 | 3.53x   |
| highly_var..  |  0.362585 |  1.114 | 3.07x   |
| scale(hvf)    |  0.020628 |  2.66  | 128x   |
| pca           |  0.123401 |  6.54  | 52.9x   |

## 1M nn cut

|    Function   | `Cell.jl` | Seurat | Speedup |
|:--------------|:----------|:-------:|:-----:|
| filter_counts |  74.04964 | 1490.6  | 20x   |
| normalize     |  72.24885 | 1844.9  | 25x   |
| highly_var..  |  36.32924 | 2716.39 | 75x   |
| scale(hvf)    |  6.894031 | 19086.6 | 2800x |
| pca           |  70.62987 | 27290.0 | 390x  |
| shared_nn     | 947.86579 | 2872.46 | 3x    |
| cluster       | 819.51106 | 9976.62 | 12x   |
| de            | 3835.9645 | 6 days  | x     |
