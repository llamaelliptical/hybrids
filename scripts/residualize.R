residualize <- function(x, y) {
  
  # Twisting and turning the residuals matrices to make an array of residual shapes without the mean centering
  b <- matrix(t(x), ncol = 3, byrow = FALSE)
  tb <- t(b)
  mat <- matrix(tb, ncol = 3, byrow = TRUE)  # Renamed from `matrix` to `mat`
  vector <- as.vector(c(mat[,1], mat[,2], mat[,3]))
  array <- array(vector, dim = c(3, dim(y$coords)[1], dim(y$coords)[3]))

  # Create an empty array with the same dimensions
  d <- array(dim = c(dim(y$coords)[1], 3, dim(y$coords)[3])) 
  
  # Loop through each slice and transpose
  for (i in 1:dim(array)[3]) {
    d[,,i] <- t(array[,,i])
  }
  
  # Assign dimnames only if they exist
  dimnames(d)[[3]] <- dimnames(y$coords)[[3]]
  
  # Add the consensus shape to these residuals
  a <- array(dim = c(dim(y$coords)[[1]], 3, dim(y$coords)[[3]]))  # Fixed hardcoded 213
  for (i in 1:dim(y$coords)[[3]]) {
    a[,,i] <- d[,,i] + y$consensus
  }
  
  # Assign dimnames safely
  dimnames(a)[[3]] <- dimnames(y$coords)[[3]]
  
  
  # Perform GPA alignment using gpagen
  if (!requireNamespace("geomorph", quietly = TRUE)) {
    stop("Package 'geomorph' is required but not installed.")
  }
  a <- geomorph::gpagen(a)
  
  return(list(coords= a$coords, residuals= d, csize=a$Csize))
  # Ensure function returns the processed array
}

# I don't know if I trust the csize here