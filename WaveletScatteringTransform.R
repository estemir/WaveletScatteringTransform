###### Scattering Transform with depth Depth and maximal Scale 2^(-MaxScale+1) ######

## Input: An matrix Img which is the image to be transformed, a constant Sigma1, corresponding to the Gauss Window,
## a constant Sigma2 and a constant Xi in R^2, both corresponding to the Morlet wavelet,
## a number Depth, a number MaxScale and a number Rot corresponding to the number of rotations,
## Output: A vector containing the scattering coefficients
## Required package: plyr

ScatNetwork <- function(Img,Sigma1,Sigma2,Xi,Size,Depth,MaxScale,Rot){

##### Constants - these are just so that the filters (Morlet and Gauss) have the "right" size #####

Size <- 13
scale <- -1

##### Preliminary functions #####

### 2d convolution using fourier trafo ###

conv2d_fft <- function(img,filter){
  ## get dimensions
  dim_filter <- dim(filter)
  dim_img <- dim(img)
  
  ## padding
  img_pad <- cbind(img,matrix(0,nrow=dim_img[1],ncol=dim_filter[2]-1))
  img_pad <- rbind(img_pad,matrix(0,nrow=dim_filter[1]-1,ncol=dim(img_pad)[2]))
  filter_pad <- cbind(filter,matrix(0,nrow=dim_filter[1],
                                    ncol=dim(img_pad)[2]-dim_filter[2]))
  filter_pad <- rbind(filter_pad,matrix(0,nrow=dim(img_pad)[1]-dim_filter[1],
                                        ncol=dim(filter_pad)[2]))
  
  n <- dim(img_pad)[1]
  m <- dim(img_pad)[2]
  
  ## fft of img
  fft_img <- mvfft(t(mvfft(t(img_pad))))
  
  ## fft of filter
  fft_filter <- mvfft(t(mvfft(t(filter_pad))))
  
  ## convolution via convolution theorem
  product <- fft_filter * fft_img
  
  ## inverse fourier transform
  ifft_img <- mvfft(t(mvfft(t(product),inverse=TRUE)),inverse=TRUE)/(n*m)
  
  ## true convolution is just submatrix with dimension of original image
  convolved <- ifft_img[floor(dim_filter[1]/2):(dim_img[1]+floor(dim_filter[1]/2)-1),
                        floor(dim_filter[2]/2):(dim_img[2]+floor(dim_filter[2]/2)-1)]
  
  ## adjust labels
  # axis labels
  x_axis <- c(1:dim_img[1])
  y_axis <- c(1:dim_img[2])
  # make matrix out of convolved
  convolved <- matrix((convolved), nrow = dim_img[1], ncol = dim_img[2], 
                      dimnames = list(x_axis, y_axis))
  
  return(convolved)
}


### calculate a gaussian window ###

## this function calculates the gaussian with sd sigma at the point u

gauss2d <- function(sigma,x,y){
  g <- 1/(2*pi*sigma^2)*exp(-0.5*(x^2+y^2)/sigma)
  return(g)
}

## this function calculates a window of dimensions 2*size+1 and at a determined scale

gauss2d_window <- function(sigma1,scale,size){
  n <- 2^(-scale)*seq(-size,size,by=1)
  z <- 2^(-2*scale)*outer(n,n,gauss2d,sigma=sigma1)
  return(z)
}


### 2d Morlet wavelet ###

## this function just calculates the morlet wavelet at a point u in R^2

mor2d <- function(xi,sigma,u){
  g <- exp(-0.5*(sum(u^2))/sigma^2)
  t <- exp(1i * t(u) %*% xi) * g
  return(as.vector(t))
}

## this function calculates the morlet window for filtering
# rotation is done with rot (note the change of angle to get the right direction)
# the scale is how much we want to zoom in
# the size lets us manipulate the dimension of the filter: dim = 2*size+1

mor2d_window <- function(xi1,sigma1,scale,rot,size){
  n <- seq(-size,size,by=1)
  z <- NULL
  A <- array(cbind(rep(n,length(n)),rep(n,each=length(n))),dim=c(length(n),length(n),2))
  B <- A
  A[,,1] <- apply(A,c(1,2),function(u) 2^scale*(cos(rot)*u[1] + sin(rot)*u[2]))
  A[,,2] <- apply(B,c(1,2),function(u) 2^scale*(-sin(rot)*u[1] + cos(rot)*u[2]))
  B <- 2^(2*scale)*apply(A,c(1,2),mor2d,xi=xi1,sigma=sigma1)
  B <- B-mean(B)
  return(B)
}

### Downsampling ###

DownSample <- function(Img,SampleRate){
  dir <- ceiling(sqrt(SampleRate))
  nrx <- floor(dim(Img)[1]/dir)
  nry <- floor(dim(Img)[2]/dir)
  ind <- as.matrix(expand.grid(seq(dir,nrx*dir,by=dir),seq(dir,nry*dir,by=dir)))
  return(matrix(Img[ind],nrow=nrx))
}


### Cartesian product ###

CartProduct = function(CurrentMatrix, NewElement)
{
  
  if (length(dim(NewElement)) != 0 )
  {
    warning("New vector has more than one dimension.")
    return (NULL)
  }
  
  if (length(dim(CurrentMatrix)) == 0)
  {
    CurrentRows = length(CurrentMatrix)
    CurrentMatrix = as.matrix(CurrentMatrix, nrow = CurrentRows, ncol = 1)
  } else {
    CurrentRows = nrow(CurrentMatrix)
  }
  
  var1 = replicate(length(NewElement), CurrentMatrix, simplify=FALSE)
  var1 = do.call("rbind", var1)
  
  var2 = rep(NewElement, CurrentRows)
  var2 = matrix(var2[order(var2)], nrow = length(var2), ncol = 1)
  
  CartProduct = cbind(var1, var2)
  return (CartProduct)
}


### Scattering network

  ## Create the two windows
  GaussWindow <- gauss2d_window(Sigma1,-MaxScale,Size)
  m <- CartProduct(1:(MaxScale-1),1:Rot)
  MorletWindow <- aaply(m,1,function(x){mor2d_window(Xi,Sigma2,-x[1],x[2]/2,Size)}) 
  ## MorletWindow[s,r,,] will be the window with scale s and rotation r
  MorletWindow <- array(MorletWindow,c(MaxScale-1,Rot,2*Size+1,2*Size+1))
  
  ## in Scat we will save the intermediate scattering results
  Scat <- list()
  ## in blur we will save the blurred results
  Blur <- list()
  ## Saving the good indices
  Indi <- list()
  
  ### the scattering function
  Scattering <- function(Step,Indices){
    
    ## generate a matrix containing all imaginable (though maybe not possible) scale/rotation pairs
    A <- array(NA,c(rep(c(MaxScale-1,Rot),Step),dim(Img)))
    ## here we apply the convolution with the corresponding Morlet wavelet to the
    ## possible indices (i.e. Indices from input)
    Result <- aaply(Indices,1,
                    function(x){
                      x <- matrix(x,nrow=1)
                      ## create all image indices
                      AllIndices <- CartProduct(x[,1:(2*(Step-1)),drop=FALSE],1:dim(Img)[1])
                      AllIndices <- CartProduct(AllIndices,1:dim(Img)[2])
                      ## get 'image' from last step
                      ToConvImg <- matrix(Scat[[Step-1]][AllIndices],nrow=dim(Img)[1])
                      ## convolve with morlet window
                      return(Mod(conv2d_fft(ToConvImg,
                                            MorletWindow[x[2*Step-1],x[2*Step],,])))
                    })
    ## fill the big matrix with the scattered images at the possible indices
    ImgIndices <- CartProduct(Indices,1:dim(Img)[1])
    ImgIndices <- CartProduct(ImgIndices,1:dim(Img)[2])
    A[ImgIndices] <- Result
    return(A)
  }
  
  Blurring <- function(Step,Indices){
    
    Result <- aaply(Indices,1,function(x){
      ## create all image indices
      AllIndices <- CartProduct(1:(dim(Img)[1]),1:(dim(Img)[2]))
      AllIndices <- cbind(t(matrix(rep(x,nrow(AllIndices)),ncol=nrow(AllIndices))),AllIndices)
      ## blurring with Gauss window and downsampling
      return(DownSample(Re(conv2d_fft(matrix(Scat[[Step]][AllIndices],nrow=dim(Img)[1]),GaussWindow)),
                        2^(2*(MaxScale-1))))
    }
    )
    ## how big is the blurred image?
    dir <- ceiling(sqrt(2^(2*(MaxScale-1))))
    nrx <- floor(dim(Img)[1]/dir)
    nry <- floor(dim(Img)[2]/dir)
    ## create indices
    ImgIndices <- CartProduct(Indices,1:nrx)
    ImgIndices <- CartProduct(ImgIndices,1:nry)
    ## generate a matrix containing all imaginable (though maybe not possible) scale/rotation pairs
    A <- array(NA,c(rep(c(MaxScale-1,Rot),Step),c(nrx,nry)))
    ## save result
    A[ImgIndices] <- Result
    ## return result as vector
    return(as.vector(Result))
  }
  
  ### do the first step seperately
  
  ## blur the original image
  Blur <- as.vector(DownSample(Re(conv2d_fft(Img,GaussWindow)),2^(2*(MaxScale-1))))
  
  ## scattering
  Indices <- CartProduct(1:(MaxScale-1),1:Rot)
  Scat[[1]] <- array(aaply(Indices,1,function(x){Mod(conv2d_fft(Img,MorletWindow[x[1],x[2],,]))}),
                     c(MaxScale-1,Rot,dim(Img)))
  ## we will save all the blurred results as one long vector. this is the second blurring step
  Blur <- c(Blur,as.vector(aaply(Scat[[1]],1:2,
                                 function(X){DownSample(Re(conv2d_fft(X,GaussWindow)),2^(2*(MaxScale-1)))})))
  #  print(length(as.vector(aaply(Scat[[1]],1:2,
  #                               function(X){DownSample(conv2d_fft(X,GaussWindow),2^(2*(MaxScale-1)))}))))
  ## now the 'good' indices
  Indi[[1]] <- Indices
  colna <- c("Scale1","Rot1")
  colnames(Indi[[1]]) <- colna
  
  ## if Depth==1 we are done
  if(Depth==1){
    return(Blur)
  }
  
  ## else continue
  for(Step in 2:Depth){
    ## get the indices
    # first generate lists corresponding to scales and rotations
    s <- replicate(Step,1:(MaxScale-1),simplify = FALSE)
    r <- replicate(Step,1:Rot,simplify = FALSE)
    # c(rbind(s,r)) interleaves them
    Indices <- as.matrix(expand.grid(c(rbind(s,r))))
    
    ## only keep those with ascending scales
    ## Here, we pick out all the indices with ascending scales
    ## aaply(Indices[,c(TRUE,FALSE)],1,diff) gives the matrix of differences along rows
    ## aaply(Indices[,c(TRUE,FALSE)],1,diff) > 0 checks if the difference is positive
    ## finally we apply "all" along rows to get just the ascending rows
    ## big matrix with all indices
    
    Indices <- Indices[aaply(aaply(Indices[,c(TRUE,FALSE)],1,diff) > 0,1,all),]
    colna <- c(colna,paste(c("Scale",as.character(Step)),collapse=""),
               paste(c("Rot",as.character(Step)),collapse=""))
    
    ## save progress
    Indi[[Step]] <- Indices 
    colnames(Indi[[Step]]) <- colna
    Scat[[Step]] <- Scattering(Step,Indices)
    Blur <- c(Blur,Blurring(Step,Indices))
  }
  ## return the large Blur vector
  return(Blur)
}
