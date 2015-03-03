### Data formatting
library(data.table)
chemdat <- fread("data_trimmed.csv", header=TRUE)
variables <- c("pH", "DOC", "Ca", "Total.Al")

## Leafed vs. Leafless
chemdat[Month <= 10 & Month >=5, canopystate := "leafed"]
chemdat[Month > 10 | Month < 5, canopystate := "leafless"]

## Get rid of negative concentrations
chemdat[Ca < 0, Ca := NA]
chemdat[Total.Al < 0, Total.Al := NA]
chemdat[DOC < 0, DOC := NA]

## Convert Ca from ueq to umol
chemdat[, Ca := Ca / 2]

## Correct names
chemdat[Sample %in% c("Be ", "Be", "Bulk"), Sample := "BE"]
chemdat[Sample == "Yp", Sample := "YP"]
chemdat[Flux == "BP", Sample := "BP"]

## Compute Ca/Al ratio
chemdat[Total.Al > 0, Ca.Al := Ca / Total.Al]
variables <- c(variables, "Ca.Al")

## Create flux-sample combination
chemdat[, Sample.Flux := sprintf("%s_%s", Sample, Flux)]

## Select only necessary fluxes and variables
idvars <- c("ID", "Year", "Month", "Day", "BulkSeq",
            "Flux", "Sample", "Sample.Flux", "canopystate")
fluxes <- c("BP", "LL", "TF", "SF", "SSF")
chem.sub <- chemdat[Flux %in% fluxes & Sample != "SSF",
                    c(idvars,variables),
                    with=FALSE]

