### Reanalysis

source("format_data.R")

# ### Non-parametric bootstrap function
# bootstrap <- function(data, iterations, size = length(data)){
#         
#         s <- sapply(1:iterations, function(x) sample(data,
#                                                      size = length(data),
#                                                      replace = TRUE))
#         means <- colMeans(s, na.rm=TRUE)
#         ci <- quantile(means, c(0.025, 0.5, 0.975))
#         return(ci)
# }
# 
# n.bs <- 5000
# 
# sample.set <- unique(chem.sub[Flux.Sample != "SSF_SSF", Flux.Sample])
# 
# ci.list <- lapply(sample.set, function(s) {
#         a <- sapply(variables, function(y) bootstrap(
#                 chem.sub[Flux.Sample == s, y, with=FALSE][[y]], n.bs))
#         l <- sapply(variables, function(y) bootstrap(
#                 chem.sub[Flux.Sample == s & canopystate == "leafed", 
#                          y, with=FALSE][[y]],
#                 n.bs))
#         d <- sapply(variables, function(y) bootstrap(
#                 chem.sub[Flux.Sample == s & canopystate == "leafless",
#                          y, with=FALSE][[y]],
#                 n.bs))
#         out <- data.frame(rbind(a, l, d))
#         out$quantile <- rep(c("2.5%", "50%", "97.5%"), 3)
#         out$canopystate <- rep(c("All", "Leafed", "Leafless"), 1, each = 3)
#         return(out)
# })
# 
# ci.table <- do.call(rbind, ci.list)
# ci.table$sample <- rep(sample.set, 1, each = nrow(ci.list[[1]]))
# ci.table <- ci.table[,c("sample", "canopystate", "quantile", variables)]

ci.table.cs <- chem.sub[,lapply(.SD, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE),
                     by=c("Sample.Flux", "canopystate"),
                     .SDcols=c("pH", "DOC", "Ca", "Total.Al", "Ca.Al")]
ci.table.all <- chem.sub[,lapply(.SD, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE),
                         by=Sample.Flux,
                         .SDcols=c("pH", "DOC", "Ca", "Total.Al", "Ca.Al")]
ci.table.all[,canopystate := "all"]
setcolorder(ci.table.all, colnames(ci.table.cs))
ci.table <- rbind(ci.table.all, ci.table.cs)
ci.table[, Sample := gsub("^(.*)_.*", "\\1", Sample.Flux)]
ci.table[, Flux := gsub("^.*_(.*)", "\\1", Sample.Flux)]
#write.csv(ci.table, "confidence_intervals_raw.csv")


