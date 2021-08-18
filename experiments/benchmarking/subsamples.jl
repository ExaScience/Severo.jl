using Cell
using ZipFile
using CSV, DataFrames

X = read_data("/data/thaber/1.3M_nn/1M_neurons_full.h5")
outputfolder = "/data/thaber/1.3_subsamples"
mkpath(outputfolder)

r = ZipFile.Reader("/data/thaber/subsamples.zip")
try
  for f in r.files
    println("processing $f")
    df = CSV.read(f, DataFrame)
    if ! isempty(df)
      cellnames = df.cell
      Y = X[cellnames,:]

      name, _ = splitext(basename(f.name))
      fname = joinpath(outputfolder, string(name, ".h5ad"))

      write_h5ad(fname, Y)
    end
  end
finally
  close(r)
end
