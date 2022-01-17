# LOAD FUNCTION FOR PLOTTING COLOUR BAR OF CHOICE
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), cex.axis0, log=FALSE, title='', labels=TRUE, add0=FALSE) {
  scale = (length(lut))/(max-min)
  
  #dev.new(width=0.5, height=5)
  if(log==FALSE){
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', 
         ylab='', main=title)
  } else if(log==TRUE) {
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', 
         ylab='', main=title)
  }
  axis(2, ticks, las=1, labels=labels,cex.axis=cex.axis0)
  for (i in 1:(length(lut))) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
