.generatePseudoGUID = function(size) {
  chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
  ret = c()
  indices = sample(1:nchar(chars), size, replace=TRUE)

  for (i in 1:size) {
    ret = c(ret, substr(chars, indices[i], indices[i]))
  }

  return(paste(ret, collapse=""))
}

SelectionType = list(NONE=0, LEAVES=1, NODE=2)
