### Load and format data -- updated data is "sub.dat"

## Load data
in.table <- read.csv("compiled_data.txt")
Variables <- c('pH','Ca','Total.Al','DOC')
full <- in.table[,c('Year','Month','Day',
                    'Bulk.vs..Seq','Sample.1','Study.Site',
                    Variables)]

## Data frame corrections
growing <- c(5, 10)
full <- within(full, {
        ## Calculate canopy state
        canopystate <- NA
        canopystate[Month >= growing[1] & Month <= growing[2]] <- "leafed"
        canopystate[is.na(canopystate)] <- "leafless"

        ## Fix tree labels
        Study.Site[Study.Site == 'Be' | Study.Site == 'Be '] <- 'BE'
        Study.Site[Study.Site == 'Bulk'] <- 'BP'
        Study.Site[Study.Site == 'Yp'] <- 'YP'

        ## Combine site and species labels
        samptype <- factor(sprintf("%s_%s", Study.Site, Sample.1))

        ## Remove negative Al values
        Total.Al[Total.Al < 0] <- NA

        ## Calculate Ca_Al ratios
        Ca.Al <- 0.5 * Ca / Total.Al
        Ca.Al[Ca.Al == Inf] <- NA
})

Variables <- c(Variables, "Ca.Al")

## Keep only meaningful variables
sample.set <- c("BP_BP", "LL_LL", "BE_TF", "BE_SF", "YP_TF", "YP_SF")
sub.dat <- subset(full, samptype %in% sample.set &
                          Bulk.vs..Seq == "Bulk")
sub.dat <- within(sub.dat, {
        Sample.1 <- factor(Sample.1)
        Study.Site <- factor(Study.Site)
        samptype <- factor(samptype)
})
canopy.set <- c("leafed", "leafless")
