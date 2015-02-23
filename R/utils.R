.generatePseudoGUID = function(size) {
  chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
  ret = c()
  indices = sample(1:nchar(chars), size, replace=T)

  for (i in 1:size) {
    ret = c(ret, substr(chars, indices[i], indices[i]))
  }

  return(paste(ret, collapse=""))
}
