## Using R for an example to show to load the Tiled Data
## and perform a simple GLM with LASSO on the data. Note: using R because GLMnet
## is much faster than scikitlearn for performing LASSO and finding optimial regularization
## This model is simplified (i.e. does not include PCA including 1000 Genomes Data, Age, Gender or use Adaptive LASSO) and does not preform bootstrapping to find which variants are stable
# This demo assumes you are using tiled data after the basic filtering step which has a different format
# than the full tiled data or the subset of the tileddata.

# Inorder to use Python
reticulate::use_python('/usr/local/bin/python3')
reticulate::py_discover_config()

# R libraries
library(Matrix)
library(foreach)
suppressMessages(library(glmnet))
library(reticulate)
library(methods)

# Python libraries
scipy <- import("scipy")
np <- import("numpy")

# Loading in Filtered Tiled Data
# See 2xpu4-4zz18-bmvaczs8gw7di41/README_Filtered_2021_05.md for detailed descriptions of inputs

# X.npy, Xr.npy, Xc.npy: compoments used to generate a sparse matrix of the filtered data
Xdata_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/X.npy'
Xrdata_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/Xr.npy'
Xcdata_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/Xc.npy'

# XPCA.npy: matrix of top 20 PCA components of the 1-hot representation of the tiled data
XPCA_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/XPCA.npy'

# y.npy: corresponding AD phenotype value for each row of the filtered tile data, 1 for case and 0 for control
ydata_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/y.npy'

# varvals.npy: vector indicating the tile variant represented by each column in the filtered tiled data matrix
tilevariant_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/varvals.npy'

# tiletag.npy: vector of tile tag for each column in the filtered tiled matrix
tiletagnumber_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/tiletag.npy'

# zygosity.npy: vector indicating zygosity for each column in the filtered tiled data matrix, 1 if tile variant present in one allele, 2 if present in both alleles
zygosity_file = '../keep/by_id/2xpu4-4zz18-bmvaczs8gw7di41/zygosity.npy'

# annotations.csv:  mapping from tile variants to HGVS annotations relative to hg38
annotation_file = '../keep/by_id/2xpu4-4zz18-1kvsfedea77ban7/annotations.csv'

# Loading in X array and make into vector
Xdata <- as.vector(np$load(Xdata_file))
Xr <- as.integer(as.vector(np$load(Xrdata_file))) + 1
Xc <- as.integer(as.vector(np$load(Xcdata_file))) + 1

# Create a new sparse matrix out of X array
Xmat <- sparseMatrix(Xr,Xc,x = Xdata)
rm(Xr,Xc,Xdata) # removing to use less memory

# Loading in PCA (this PCA is for just the AD WGS genomes)
XPCA <- as.vector(np$load(XPCA_file))
XPCA <- as(matrix(XPCA,nco=20,nrow=length(XPCA)/20),"sparseMatrix")

# Combining X sparse matrix with PCA components
Xmat = cbind(Xmat,XPCA)

# Loading the y array and make into vector
ynump <- np$load(ydata_file)
y <- as.vector(ynump)
rm(ynump) # removing to use less memory

# Loading in tile variant numbers, tile tag numbers, and zygosity
varvals <- as.integer(as.vector(np$load(tilevariant_file)))
tiletag <- as.integer(as.vector(np$load(tiletagnumber_file)))
zygosity <- as.integer(as.vector(np$load(zygosity_file)))

# Setting up training and testing data for cross validation (80:20 training:testing split)
dt = sort(sample(nrow(Xmat), nrow(Xmat)*0.8))
Xtrain <- Xmat[dt,]
Xtest <- Xmat[-dt,]
ytrain <- y[dt]
ytest <- y[-dt]

# Caluclating weights for Weighted Linear Regression since our case/control ratio is not 50:50
fraction_0_train <- rep(1 - sum(ytrain == 0) / length(ytrain), sum(ytrain == 0))
fraction_1_train <- rep(1 - sum(ytrain == 1) / length(ytrain), sum(ytrain == 1))
fraction_0_test <- rep(1 - sum(ytest == 0) / length(ytest), sum(ytest == 0))
fraction_1_test <- rep(1 - sum(ytest == 1) / length(ytest), sum(ytest == 1))

wtrain <- rep(1,length(ytrain))
wtest <- rep(1,length(ytest))

wtrain[ytrain == 0] <- fraction_0_train
wtrain[ytrain == 1] <- fraction_1_train
wtest[ytest == 0] <- fraction_0_test
wtest[ytest == 1] <- fraction_1_test

# Finding Best Regularization (LASSO)
cv.lasso.class <- cv.glmnet(Xtrain, ytrain, family='binomial', alpha=1, nfolds=10, parallel=FALSE, keep=FALSE,standardize=FALSE, weights=wtrain, type.measure='class', trace.it=1)

# Plot and Save Figure
png('regularization_lasso.png')
plot(cv.lasso.class)
dev.off()

# Assessing model
preds = predict(cv.lasso.class, newx = Xtest, s = "lambda.min", weights = wtest)
outcome = assess.glmnet(preds, newy = ytest, family = "binomial", weight = wtest)
print(outcome$class)

# Confusion Matrix
cm = confusion.glmnet(cv.lasso.class, newx = Xtest, newy = ytest, weights = wtest, s = "lambda.min")
print(cm)

# Find regularization and nonzero coefficents corresponding to minimum "error"
coefVec <- coef(cv.lasso.class, s = "lambda.min") # minimum error
coefVec <- coefVec[-1]
idxnzmin <- which(coefVec != 0)
nzcoefVal <- coefVec[idxnzmin]

# Adding to varvals labels to represent PC components
PClist <- seq(1, 20)
PClist <- paste0("PC", PClist)
varvals = append(varvals, PClist)

# Sorting to find largest coefficent values, collecting tiledata in dataframe
varvals <- varvals[idxnzmin]
zygosity <- zygosity[idxnzmin]
tiletag <- tiletag[idxnzmin]

# Displaying tile information for non-zero coefficents
tiledata <- data.frame("nonnzerocoefs" = nzcoefVal, "tiletag" = tiletag, "varvals" = varvals, "zygosity" = zygosity)
idxsort <- order(abs(tiledata$nonnzerocoefs), decreasing = TRUE)
tiledata <- tiledata[idxsort,]
print(head(tiledata,20))

# Finding HGVS annotations for top 20 tile variants in model
# Note the following HGVS annotations may be of importance
# Tile 9553646, variant 2 --> chr19:g.44908684T>C, https://www.ncbi.nlm.nih.gov/snp/rs429358

# Helper function to look up tile variants in the annotation file
print_annotation <- function(tiletag, varvals) {
   pattern = paste("", tiletag, as.character(varvals), "", sep=",")
   command = paste("grep", pattern, annotation_file)
   returncode <- system(command)
   # return code is 1 if and only if there is no matching result, i.e., if the tile variant is the same as ref
   return(returncode)
}

print("Finding annotations to top 20 tile variants")

for (i in 1:20) {
   tiletag = tiledata$tiletag[i]
   varvals = as.numeric(tiledata$varvals[i])-1 # -1 offset required for filtered data
   rc = print_annotation(tiletag, varvals)
   if (varvals != 1 & rc == 1) { # if that tile variant is ref, look up the most common tile variant in that position instead
     for (i in 1:varvals-1) {
       rc = print_annotation(tiletag, i)
       if (rc == 0) {
         break
       }
     }
   }
}
