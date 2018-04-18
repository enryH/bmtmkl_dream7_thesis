#load results of grid searches obtained from 
# "ch4_grid_search_parallel.R
###################################################################################################
mygrid <- function (alpha, beta){
  grid <- expand.grid(alpha_lambda  = alpha, beta_lambda  = beta,
                      alpha_upsilon = alpha, beta_upsilon = beta, 
                      alpha_gamma   = alpha, beta_gamma   = beta,
                      alpha_omega   = alpha, beta_omega   = beta, 
                      alpha_epsilon = alpha, beta_epsilon = beta  )
  return(grid)
}
###################################################################################################
load("data/results_171224_00-43_df.RData")
#
select <- c(39967, 41425,40696, 42154, 42883, 43612, 44341, 45070, 45799)
table_max <- matrix(0, length(result_i), length(select))
i=1
for (k in select) {
  cat("grid entry: ", k , "\n")
  table_max[,i]<-  unlist(runmodel(setparameters_from_grid(k)))
  i = i + 1
}
colnames(table_max) <- select
rownames(table_max) <- names(result_i)
rm("i", "k", "select")
#
table_max
currenttable <- xtable(table_max)
print(currenttable, type="latex", file=("paper/Tables/tab_min_test_error_comp"))
# grid 2 
grid2 <- mygrid(alpha=c(10^(-10),  1, 10), beta=c(10^(-10), 0.01, 1) )
# 3^10 = 59049 combinations 

load("data/results_180108_02-53_df.RData")
results_aroundone <- df
View(results_aroundone)
load("data/results_180107_21-03_df.RData")
View(df)
#########################################################################################
View(results_default)
load("data/results_180107_15-13_df.RData")
cat("alpha: ", unique(df[,1]), "beta: ", unique(df[,2]))
par_max_ELBO <- df[39303,]
par_min_ELBO <- df[73,]
par_min_MSE_ins <- df[52561,]
par_min_MSE_ofs <- df[52561,]
#########################################################################################
View(df)
results_default <- df

load("data/results_180104_19-20_df.RData")
View(df)


results_grid2 <- 