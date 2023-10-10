# line tool 
bez2poly <- function(x){
  x <- strsplit( gsub("[ ].*[ ]","l",x) ,split = "[a-zA-Z]")[[1]][-1]
  x <- as.data.frame(strsplit(x,split = ","))
  x <- apply(x,1,as.numeric)
  x <- cbind(cumsum(x[,1]), cumsum(x[,2]) )
  x <- rbind(x,x[1,])
  colnames(x) <- c("x", "y")
  return(x)
}

tmat <- function(x,n=3){
  tmp <- matrix(0,3,3)
  x <- as.numeric(strsplit(gsub("[a-z()]","",x),split = ",")[[1]])
  tmp[c(1,2),c(1,1)] <- x[1:2]
  tmp[c(1,2),c(2,2)] <- x[3:4]
  tmp[c(1,2),c(3,3)] <- x[5:6]
  tmp[3,3] <- 1
  return(tmp)
}

rec2poly <- function(x){
  x <- attributes(x)
  x <- setNames( as.numeric(unlist(x)[c("x","y","width","height")]) , c("x","y","width","height") )
  df <- data.frame(
    x = x["x"] + x["width"] * c(0,0,1,1,0),
    y = x["y"] + x["height"] * c(0,1,1,0,0))
  return(df)
}

get_shape <- function(x){
  xx <- attributes(x)
  if( "d" %in% names(xx) ){
    return( bez2poly( attr(x,"d") ) ) # path annot
  } else if( "width"  %in% names(xx) ){
    return(rec2poly( x )) # rect annot
  } 
}

poly_area <- function(df){
  res <- t(df)
  x <- abs( sum(  ((res[1,] * c(res[2,-1],res[2,1])) - (res[2,] * c(res[1,-1],res[1,1]) )) )/2 )
  return(x)
}
