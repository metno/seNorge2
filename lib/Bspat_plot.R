#!
pointspat<-function(x=NULL,y=NULL,yvar=NULL,yvar1=NULL,ydqc=NULL,xvar=NULL,brk=NULL,col=NULL,
                    namefileout=NULL,mtxt=NULL,mtxt1=NULL,orog=NULL,pos="bottomright",
                    xl=NULL,yl=NULL,bnd=NULL,cx=NULL,colext=NULL,
                    yvartext=NULL,yvartextcex=NULL,xvar.orog=NULL,
                    colpoints=NULL,legcex=NULL) {
# ydqc==0 good; ydqc==-1 missing; ydqc==1 erroneous
#----------------------------------------------------------------------------------
  if (is.null(legcex)) legcex<-0.8
  if (is.null(colext)) {
    col1<-c("#FFFFFF",col,"#000000")
  } else {
    col1<-c(colext[1],col,colext[2])
  }
  if (is.null(xvar)){
    png(file=namefileout,width=1200,height=1200)
    image(xvar.orog,breaks=seq(0,2500,length=10),col=gray.colors(9),main=mtxt,
          xlab="",ylab="",xlim=xl,ylim=yl,cex.main=1.6)
    contour(xvar.orog,levels=c(500,1500),drawlabels=F,col="black",lwd=1,add=T)
    plot(bnd,add=T)
    legtxt<-vector()
# good
    for (b in 1:length(col1)) {
      cond<-NULL
      if (b==1) cond<-which(yvar<brk[b])
      if (b==length(col1)) cond<-which(yvar>=brk[b-1])
      if (b>1 & b<length(col1)) cond<-which( (yvar>=brk[b-1]) & (yvar<brk[b]))
      if (b==1) legtxt[b]<-paste("<",brk[b],sep="")
      if (b==length(col1)) legtxt[b]<-paste(">=",brk[b-1],sep="")
      if (b>1 & b<length(col1)) legtxt[b]<-paste(brk[b-1]," ",brk[b],sep="")
      points(x[cond],y[cond],pch=19,cex=cx[b],col=col1[b])
      points(x[cond],y[cond],cex=(cx[b]+0.1),col=colpoints)
    }
    legend(x=pos,legend=legtxt,cex=legcex,fill=col1)
    if (length(mtxt1)>0) mtext(mtxt1,side=3,cex=1.6)
    if (!is.null(xvar.orog)) contour(xvar.orog,levels=c(0,500,1500),drawlabels=F,col="black",lwd=0.8,add=T)
    dev.off()
    return()
  }
  if (is.null(x)){
    png(file=namefileout,width=1200,height=1200)
    image(xvar,breaks=brk,col=col,main=mtxt,
          xlab="",ylab="",xlim=xl,ylim=yl,cex.main=1.6)
    contour(xvar,levels=c(0,-10,-20),drawlabels=F,col=c("red","darkgreen","darkblue"),lwd=2,add=T)
    plot(bnd,add=T)
    if (length(colext)==2) {
      mx<-cellStats(xvar,max)
      mn<-cellStats(xvar,min)
      image(xvar,breaks=c(brk[length(brk)],mx),add=T,col=colext[2])
      image(xvar,breaks=c(mn,brk[1]),add=T,col=colext[1])
    }
#  points(x,y)
    legtxt<-vector()
# good
    for (b in 1:length(col1)) {
      if (b==1) legtxt[b]<-paste("<",brk[b],sep="")
      if (b==length(col1)) legtxt[b]<-paste(">=",brk[b-1],sep="")
      if (b>1 & b<length(col1)) legtxt[b]<-paste(brk[b-1]," ",brk[b],sep="")
    }
    if (!is.null(xvar.orog)) contour(xvar.orog,levels=c(0,100,250,500,1500),drawlabels=F,col="black",lwd=0.8,add=T)
    legend(x=pos,legend=legtxt,cex=legcex,fill=col1)
    if (length(mtxt1)>0) mtext(mtxt1,side=3,cex=1.6)
    dev.off()
    return()
  }
  if (length(ydqc)==0) ydqc<-rep(0,length=length(yvar))
  if (length(colpoints)==0) colpoints<-"black"
  png(file=namefileout,width=1200,height=1200)
  plot(x[ydqc==0],y[ydqc==0],main=mtxt,xlab="",ylab="",xlim=xl,ylim=yl,cex.main=1.6)
  plot(bnd,add=T)
#  if (!is.null(orog)) contour(orog,levels=c(1,200,700,1500),add=T,col="black",lwd=1)
  image(xvar,breaks=brk,add=T,col=col)
  if (!is.null(xvar.orog)) contour(xvar.orog,levels=c(0,500,1500),drawlabels=F,col="black",lwd=0.8,add=T)
  contour(xvar,levels=c(0,-10,-20),drawlabels=F,col=c("red","darkgreen","darkblue"),lwd=2,add=T)
  if (length(colext)==2) {
    mx<-cellStats(xvar,max)
    mn<-cellStats(xvar,min)
    image(xvar,breaks=c((brk[length(brk)]+0.001),mx),add=T,col=colext[2])
    image(xvar,breaks=c(mn,(brk[1]-0.001)),add=T,col=colext[1])
  }
#  points(x,y)
  legtxt<-vector()
# good
  for (b in 1:length(col1)) {
    cond<-NULL
    if (b==1) cond<-which(yvar<brk[b] & ydqc>=0)
    if (b==length(col1)) cond<-which(yvar>=brk[b-1] & ydqc>=0)
    if (b>1 & b<length(col1)) cond<-which( (yvar>=brk[b-1]) & (yvar<brk[b]) & ydqc>=0)
    if (b==1) legtxt[b]<-paste("<",brk[b],sep="")
    if (b==length(col1)) legtxt[b]<-paste(">=",brk[b-1],sep="")
    if (b>1 & b<length(col1)) legtxt[b]<-paste(brk[b-1]," ",brk[b],sep="")
    points(x[cond],y[cond],pch=19,cex=cx[b],col=col1[b])
    points(x[cond],y[cond],cex=(cx[b]+0.1),col=colpoints)
  }
# erroneous 
  points(x[ydqc==1],y[ydqc==1],col="black",pch=4,lwd=3)
# missing
  if (length(yvar1)>0) {
    for (b in 1:length(col1)) {
      cond<-NULL
      if (b==1) cond<-which(yvar1<brk[b] & ydqc==-1)
      if (b==length(col1)) cond<-which(yvar1>=brk[b-1] & ydqc==-1)
      if (b>1 & b<length(col1)) cond<-which( (yvar1>=brk[b-1]) & (yvar1<brk[b]) & ydqc==-1)
      points(x[cond],y[cond],pch=19,cex=0.5,col=col1[b])
      points(x[cond],y[cond],cex=0.6,col=colpoints)
    }
  }
#
  if (!is.null(yvartext)) {
    text(x,y,labels=yvartext,cex=yvartextcex)
  }
#
  legend(x=pos,legend=legtxt,cex=legcex,fill=col1)
  if (length(mtxt1)>0) mtext(mtxt1,side=3,cex=1.6)
  if (!is.null(xvar.orog)) contour(xvar.orog,levels=c(0,500,1500),drawlabels=F,col="black",lwd=0.8,add=T)
  dev.off()
  return()
}

#!
rainspatplot<-function(x=NULL,y=NULL,yvar=NULL,yvar1=NULL,ydqc=NULL,xvar=NULL,brk=NULL,col=NULL,
                       namefileout=NULL,mtxt=NULL,mtxt1=NULL,pos="bottomright",
                       xl=NULL,yl=NULL,bnd=NULL,cx=NULL,colext=NULL,
                       yvartext=NULL,yvartextcex=NULL,xvar.orog=NULL,
                       colpoints=NULL,legcex=NULL) {
# ydqc==0 good; ydqc==-1 missing; ydqc==1 erroneous
#----------------------------------------------------------------------------------
  if (is.null(legcex)) legcex<-0.8
  col1<-c(colext[1],col,colext[2])
  if (length(ydqc)==0) ydqc<-rep(0,length=length(yvar))
  if (length(colpoints)==0) colpoints<-"black"
  y.miss<-which(ydqc==-1)
  png(file=namefileout,width=1200,height=1200)
  plot(x[ydqc==0],y[ydqc==0],main=mtxt,xlab="",ylab="",xlim=xl,ylim=yl,cex.main=1.6)
  plot(bnd,add=T)
#  if (!is.null(orog)) contour(orog,levels=c(1,200,700,1500),add=T,col="black",lwd=1)
  image(xvar.orog,breaks=c(0,500,1000,1500,2000,2500),col=gray(seq(0.7,1,length=5)),add=T)
  xvar1<-xvar
  xvar1[xvar==0]<-NA
  image(xvar1,breaks=brk,add=T,col=col)
#  contour(xvar,levels=c(0,-10,-20),drawlabels=F,col=c("red","darkgreen","darkblue"),lwd=2,add=T)
  if (length(colext)==2) {
    mx<-cellStats(xvar,max)
    mn<-cellStats(xvar,min)
    image(xvar1,breaks=c((brk[length(brk)]+0.001),mx),add=T,col=colext[2])
    image(xvar1,breaks=c(mn,(brk[1]-0.001)),add=T,col=colext[1])
  }
#  points(x,y)
  legtxt<-vector()
# good
  n.col1<-length(col1)
  legtxt[1]<-"(No rain)"
  for (b in 1:n.col1) {
    cond<-NULL
    if (b==1) {
      cond<-which(yvar<brk[b] & ydqc>=0)
    }
    if (b==n.col1) {
      cond<-which(yvar>=brk[b-1] & ydqc>=0)
      legtxt[n.col1-b+2]<-paste("[",brk[b-1],",  >mm",sep="")
    }
    if (b>1 & b<n.col1) {
      cond<-which( (yvar>=brk[b-1]) & (yvar<brk[b]) & ydqc>=0)
      legtxt[n.col1-b+2]<-paste("[",brk[b-1],", ",brk[b],"> mm",sep="")
    }
    points(x[cond],y[cond],pch=19,cex=cx[b],col=col1[b])
    points(x[cond],y[cond],cex=(cx[b]+0.1),col=colpoints)
  }
# erroneous 
  points(x[ydqc>0],y[ydqc>0],col="black",pch=4,lwd=3)
# missing
  if (length(y.miss)>0) {
    for (b in 1:length(col1)) {
      cond<-NULL
      if (b==1) cond<-which(yvar[ydqc==-1]<brk[b])
      if (b==length(col1)) cond<-which(yvar[ydqc==-1]>=brk[b-1])
      if (b>1 & b<length(col1)) cond<-which( (yvar[ydqc==-1]>=brk[b-1]) & (yvar[ydqc==-1]<brk[b]))
      if (length(cond)>0) points(x[ydqc==-1][cond],y[ydqc==-1][cond],pch=19,cex=cx[b],col=col1[b])
      if (length(cond)>0) points(x[ydqc==-1][cond],y[ydqc==-1][cond],cex=(cx[b]+0.1),col=colpoints)
    }
  }
#
  if (!is.null(yvartext)) {
    text(x,y,labels=yvartext,cex=yvartextcex)
  }
#
  legend(x=pos,legend=legtxt,cex=legcex,fill=c(col1[1],col1[n.col1:2]))
  if (length(mtxt1)>0) mtext(mtxt1,side=3,cex=1.6)
  contour(xvar.orog,levels=c(0,100,250,500,1500),drawlabels=F,col="black",lwd=0.8,add=T)
  dev.off()
  return()
}

#!
PRECplot<-function(namefileout=NULL,
                   y.data=NULL,
                   r.data=NULL,
#                   scale=NULL,
#                   col.scale=NULL,
                   orog=NULL,
                   bound=NULL,
                   par=NULL) {
#   NA  observation is NA
#   -1  missing DQC info
#   0   good observation
#   100 bad  observation: KDVH flag > 2 | observation not good in external DQC | 
#                         observed value not plausible | station in blacklist/s  
#   200 bad  observation: dry-station surrounded only by wet-stations (close enough)
#   300 bad  observation: wet-stations surrounded only by dry-stations (close enough)
#   400 bad  observation: dry observation is (1) not included in a dry area
#                         (2) is in Norway 
#   500 bad  observation: wet observation is (1) not included in an event (2) in Norway


#x=NULL,y=NULL,yvar=NULL,yvar1=NULL,ydqc=NULL,xvar=NULL,brk=NULL,col=NULL,
#                   namefileout=NULL,mtxt=NULL,mtxt1=NULL,pos="bottomright",
#                   xl=NULL,yl=NULL,bnd=NULL,cx=NULL,colext=NULL,
#                   yvartext=NULL,yvartextcex=NULL,xvar.orog=NULL,
#                   colpoints=NULL,legcex=NULL) {
## ydqc==0 good; ydqc==-1 missing; ydqc==1 erroneous
#----------------------------------------------------------------------------------
  y.data$x<-as.numeric(y.data$x)
  y.data$y<-as.numeric(y.data$y)
  leg.str<-"no rain"
  n.col<-length(par$col.scale)
  #
  png(file=namefileout,width=1200,height=1200)
  plot(y.data$x[!is.na(y.data$yo)],y.data$y[!is.na(y.data$yo)],
       main=par$main,xlab=par$xlab,ylab=par$ylab,xlim=par$xl,ylim=par$yl,cex.main=1.6,col="white")
  image(orog,breaks=c(0,500,1000,1500,2000,2500),col=gray(seq(0.7,1,length=5)),add=T)
  image(r.data,col=par$col.scale,breaks=par$scale,add=T)
  # dry observations
  dry<-!is.na(y.data$yo) & (y.data$yo<0.1)
  aux<-which(dry & y.data$dqcflag<=0)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],col="black",bg="gray",pch=21,cex=1.2)
  aux<-which(dry & y.data$dqcflag==100)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],col="black",bg="gray",pch=24,cex=1.2)
  aux<-which(dry & y.data$dqcflag>100)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],col="black",bg="gray",pch=25,cex=1.2)
  # wet observations
  for (c in 1:n.col) {
    in.break<-(!is.na(y.data$yo)) & (y.data$yo>=par$scale[c]) & (y.data$yo<par$scale[c+1])
    aux<-which(in.break & y.data$dqcflag<=0)
    if (length(aux)>0) 
      points(y.data$x[aux],y.data$y[aux],col="black",bg=par$col.scale[c],pch=21,cex=1.2)
    if (c==1) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=1),", ",
                               formatC(par$scale[c+1],format="f",digits=1),") mm",sep=""))
    } else if (c==2) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=1),", ",
                               formatC(par$scale[c+1],format="f",digits=0),") mm",sep=""))
    } else if (c<n.col) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=0),", ",
                           formatC(par$scale[c+1],format="f",digits=0),") mm",sep=""))
    } else if (c==n.col) {
      leg.str<-c(leg.str,paste(">",formatC(par$scale[c],format="f",digits=0),"mm",sep=""))
    }
    aux<-which(in.break & y.data$dqcflag==100)
    if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],bg=par$col.scale[c],col="black",pch=24,cex=1.2)
    aux<-which(in.break & y.data$dqcflag>100)
    if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],bg=par$col.scale[c],col="black",pch=25,cex=1.2)
  }
  plot(bound,add=T)
  legend(x="bottomright",fill=rev(c("gray",par$col.scale)),legend=rev(leg.str),cex=1.5)
#  contour(orog,levels=c(0,100,250,500,1500),drawlabels=F,col="black",lwd=0.8,add=T)
#  contour(r.data,levels=c(5,10,30,70,100),drawlabels=F,col="black",lwd=0.8,add=T)
  dev.off()
  return()
}

#!
TEMPplot<-function(namefileout=NULL,
                   y.data=NULL,
                   r.data=NULL,
                   orog=NULL,
                   bound=NULL,
                   par=NULL) {
## ydqc==0 good; ydqc==-1 missing; ydqc==1 erroneous
#----------------------------------------------------------------------------------
  y.data$x<-as.numeric(y.data$x)
  y.data$y<-as.numeric(y.data$y)
  leg.str<-vector()
  n.col<-length(par$col.scale)
  #
  png(file=namefileout,width=1200,height=1200)
  plot(y.data$x[!is.na(y.data$yo)],y.data$y[!is.na(y.data$yo)],
       main=par$main,xlab=par$xlab,ylab=par$ylab,xlim=par$xl,ylim=par$yl,cex.main=1.6,col="white")
  image(orog,breaks=c(0,500,1000,1500,2000,2500),col=gray(seq(0.7,1,length=5)),add=T)
  image(r.data,col=par$col.scale,breaks=par$scale,add=T)
  # observations
  for (c in 1:n.col) {
    in.break<-(!is.na(y.data$yo)) & (y.data$yo>=par$scale[c]) & (y.data$yo<par$scale[c+1])
    aux<-which(in.break & y.data$dqcflag<=0)
    if (length(aux)>0) 
      points(y.data$x[aux],y.data$y[aux],col="black",bg=par$col.scale[c],pch=21,cex=1.2)
    if (c==1) {
      leg.str[1]<-paste("<",formatC(par$scale[c+1],format="f",digits=0),sep="")
    } else if (c<n.col) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=0),", ",
                           formatC(par$scale[c+1],format="f",digits=0),")",sep=""))
    } else if (c==n.col) {
      leg.str<-c(leg.str,paste(">",formatC(par$scale[c],format="f",digits=0),sep=""))
    }
    aux<-which(in.break & y.data$dqcflag>0)
    if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],bg=par$col.scale[c],col="black",pch=24,cex=1.2)
  }
  plot(bound,add=T)
  legend(x="bottomright",fill=rev(c(par$col.scale,"white")),legend=rev(c(leg.str,expression(paste("Temp,",degree,"C")))),cex=1.5)
  dev.off()
  return()
}

# Function to plot color bar
color.bar <- function(col, scale, nticks=11, title='') {
    usr<-as.numeric(par("usr"))
    min.y<-usr[3]
    max.y<-usr[4]
    scale.y <- length(col)/(max.y-min.y)
    i.ind<-vector(length=10)
    for (i in 1:10) i.ind[i]<-min(which(scale>(10*i)))-1
    print(usr)
    print(scale.y)
    print(i.ind)

#    dev.new(width=1.75, height=5)
#    par(new=T)
    dx<-(usr[2]-usr[1])/40
    x1<-usr[1]+dx*39
    print(paste(dx,x1,min.y,max.y))
#    plot(c(x1,x1+dx), c(min.y,max.y), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
#    axis(2, ticks, las=1)
    for (i in 1:length(col)) {
     y = (i-1)/scale.y + min.y
     rect(x1,y,x1+dx,y+1/scale.y, col=col[i], border=NA)
    }
    text(x=x1+dx/2,y=((i.ind-1)/scale.y + min.y),labels=round(scale[i.ind],0),cex=1.8,
         col="black")
}


#!
TEMPplot.IDI<-function(namefileout=NULL,
                   y.data=NULL,
                   r.data=NULL,
                   orog=NULL,
                   bound=NULL,
                   par=NULL) {
## ydqc==0 good; ydqc==-1 missing; ydqc==1 erroneous
#----------------------------------------------------------------------------------
  y.data$x<-as.numeric(y.data$x)
  y.data$y<-as.numeric(y.data$y)
  leg.str<-""
  n.col<-length(par$col.scale)
  yidi<-suppressWarnings(as.numeric(y.data$yidi)*100)
  #
  png(file=namefileout,width=1200,height=1200)
  plot(y.data$x[!is.na(y.data$yo)],y.data$y[!is.na(y.data$yo)],
       main=par$main,xlab=par$xlab,ylab=par$ylab,xlim=par$xl,ylim=par$yl,cex.main=1.6,col="white")
  image(orog,breaks=c(0,500,1000,1500,2000,2500),col=gray(seq(0.7,1,length=5)),add=T)
  image(r.data,col=par$col.scale,breaks=par$scale,add=T)
  # observations
  for (c in 1:n.col) {
    in.break<-(!is.na(y.data$yo)) & (yidi>=par$scale[c]) & (yidi<par$scale[c+1])
    aux<-which(in.break & y.data$dqcflag<=0)
    if (length(aux)>0) 
      points(y.data$x[aux],y.data$y[aux],col="black",bg=par$col.scale[c],pch=21,cex=1.2)
    if (c==1) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=1),", ",
                               formatC(par$scale[c+1],format="f",digits=1),") mm",sep=""))
    } else if (c==2) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=1),", ",
                               formatC(par$scale[c+1],format="f",digits=0),") mm",sep=""))
    } else if (c<n.col) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=0),", ",
                           formatC(par$scale[c+1],format="f",digits=0),") mm",sep=""))
    } else if (c==n.col) {
      leg.str<-c(leg.str,paste(">",formatC(par$scale[c],format="f",digits=0),"mm",sep=""))
    }
    aux<-which(in.break & y.data$dqcflag>0)
#    if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],bg=par$col.scale[c],col="black",pch=24,cex=2.2)
    if (length(aux)>0) {
        print("prek")
        points(y.data$x[aux],y.data$y[aux],bg="black",col="black",pch=24,cex=1.2)
    }
  }
  plot(bound,add=T)
#  legend(x="bottomright",fill=rev(par$col.scale),legend=rev(leg.str),cex=1)
  print(paste(length(par$col.scale),length(par$scale)))
  color.bar(par$col.scale,par$scale)
  dev.off()
  return()
}

#!
PRECplot.IDI<-function(namefileout=NULL,
                   y.data=NULL,
                   r.data=NULL,
                   orog=NULL,
                   bound=NULL,
                   par=NULL) {
#   NA  observation is NA
#   -1  missing DQC info
#   0   good observation
#   100 bad  observation: KDVH flag > 2 | observation not good in external DQC | 
#                         observed value not plausible | station in blacklist/s  
#   200 bad  observation: dry-station surrounded only by wet-stations (close enough)
#   300 bad  observation: wet-stations surrounded only by dry-stations (close enough)
#   400 bad  observation: dry observation is (1) not included in a dry area
#                         (2) is in Norway 
#   500 bad  observation: wet observation is (1) not included in an event (2) in Norway
#----------------------------------------------------------------------------------
  y.data$x<-as.numeric(y.data$x)
  y.data$y<-as.numeric(y.data$y)
  leg.str<-""
  n.col<-length(par$col.scale)
  yidi<-suppressWarnings(as.numeric(y.data$yidi))*100
  #
  png(file=namefileout,width=1200,height=1200)
  plot(y.data$x[!is.na(y.data$yo)],y.data$y[!is.na(y.data$yo)],
       main=par$main,xlab=par$xlab,ylab=par$ylab,xlim=par$xl,ylim=par$yl,cex.main=1.6,col="white")
  image(orog,breaks=c(0,500,1000,1500,2000,2500),col=gray(seq(0.7,1,length=5)),add=T)
  image(r.data,col=par$col.scale,breaks=par$scale,add=T)
  # dry observations
  dry<-!is.na(y.data$yo) & (y.data$yo<0.1)
  aux<-which(dry & y.data$dqcflag<=0)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],col="black",bg="gray",pch=21,cex=1.2)
  aux<-which(dry & y.data$dqcflag==100)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],col="black",bg="black",pch=24,cex=1.2)
  aux<-which(dry & y.data$dqcflag>100)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],col="black",bg="black",pch=25,cex=1.2)
  # wet observations
  for (c in 1:n.col) {
    in.break<-(!is.na(y.data$yo)) & (yidi>=par$scale[c]) & (yidi<par$scale[c+1])
    aux<-which(in.break & y.data$dqcflag<=0)
    if (length(aux)>0) 
      points(y.data$x[aux],y.data$y[aux],col="black",bg=par$col.scale[c],pch=21,cex=1.2)
    if (c==1) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=1),", ",
                               formatC(par$scale[c+1],format="f",digits=1),") mm",sep=""))
    } else if (c==2) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=1),", ",
                               formatC(par$scale[c+1],format="f",digits=0),") mm",sep=""))
    } else if (c<n.col) {
      leg.str<-c(leg.str,paste("[",formatC(par$scale[c],format="f",digits=0),", ",
                           formatC(par$scale[c+1],format="f",digits=0),") mm",sep=""))
    } else if (c==n.col) {
      leg.str<-c(leg.str,paste(">",formatC(par$scale[c],format="f",digits=0),"mm",sep=""))
    }
  }
  aux<-which(y.data$dqcflag==100)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],bg="black",col="red",pch=24,cex=1.2)
  aux<-which(y.data$dqcflag>100)
  if (length(aux)>0) points(y.data$x[aux],y.data$y[aux],bg="black",col="red",pch=25,cex=1.2)
  plot(bound,add=T)
#  legend(x="bottomright",fill=rev(par$col.scale),legend=rev(leg.str),cex=1)
  print(paste(length(par$col.scale),length(par$scale)))
  color.bar(par$col.scale,par$scale)
  dev.off()
  return()
}

#!
rainspatplot.cra<-function(x=NULL,y=NULL,yvar=NULL,yvar1=NULL,ydqc=NULL,xvar=NULL,brk=NULL,col=NULL,
                           namefileout=NULL,mtxt=NULL,mtxt1=NULL,legpos="bottomright",
                           xl=NULL,yl=NULL,bnd=NULL,cx=NULL,colext=NULL,
                           yvartext=NULL,yvartextcex=NULL,xvar.orog=NULL,
                           cra.lab=NULL,cra.x=NULL,cra.y=NULL,
                           n.cra=NULL,
                           cra.stn.dens=NULL,
                           mean=NULL,
                           max.x=NULL,
                           max.yo=NULL,
                           max.ya=NULL,
                           max.yav=NULL,
                           cv.rms.rel=NULL,
                           cv.bias.sq=NULL,
                           cv.rmse.sq=NULL,
                           cv.made.sq=NULL,
                           colpoints=NULL,legcex=NULL) {
# ydqc==0 good; ydqc==-1 missing; ydqc==1 erroneous
#----------------------------------------------------------------------------------
  if (is.null(legcex)) legcex<-0.8
  col1<-c(colext[1],col,colext[2])
  if (length(ydqc)==0) ydqc<-rep(0,length=length(yvar))
  if (length(colpoints)==0) colpoints<-"black"
  png(file=namefileout,width=1200,height=1200)
  plot(x[ydqc==0],y[ydqc==0],main=mtxt,xlab="",ylab="",xlim=xl,ylim=yl,cex.main=1.6)
  plot(bnd,add=T)
#  if (!is.null(orog)) contour(orog,levels=c(1,200,700,1500),add=T,col="black",lwd=1)
  image(xvar.orog,breaks=c(0,500,1000,1500,2000,2500),col=gray(seq(0.7,1,length=5)),add=T)
  xvar1<-xvar
  xvar1[xvar==0]<-NA
  image(xvar1,breaks=brk,add=T,col=col)
  contour(xvar.orog,levels=c(0,500,1500),drawlabels=F,col="black",lwd=0.8,add=T)
#  contour(xvar,levels=c(0,-10,-20),drawlabels=F,col=c("red","darkgreen","darkblue"),lwd=2,add=T)
  if (length(colext)==2) {
    mx<-cellStats(xvar,max)
    mn<-cellStats(xvar,min)
    image(xvar1,breaks=c((brk[length(brk)]+0.001),mx),add=T,col=colext[2])
    image(xvar1,breaks=c(mn,(brk[1]-0.001)),add=T,col=colext[1])
  }
#  points(x,y)
  legtxt<-vector()
# good
  n.col1<-length(col1)
  legtxt[1]<-"(No rain)"
  for (b in 1:n.col1) {
    cond<-NULL
    if (b==1) {
      cond<-which(yvar<brk[b] & ydqc>=0)
    }
    if (b==n.col1) {
      cond<-which(yvar>=brk[b-1] & ydqc>=0)
      legtxt[n.col1-b+2]<-paste("[",brk[b-1],",",sep="")
    }
    if (b>1 & b<n.col1) {
      cond<-which( (yvar>=brk[b-1]) & (yvar<brk[b]) & ydqc>=0)
      legtxt[n.col1-b+2]<-paste("[",brk[b-1],", ",brk[b],sep="")
    }
    points(x[cond],y[cond],pch=19,cex=cx[b],col=col1[b])
    points(x[cond],y[cond],cex=(cx[b]+0.1),col=colpoints)
  }
  text(cra.x,cra.y,cra.lab,cex=2.5)
  Lx<-xl[2]-xl[1]
  Lx.3<-Lx/3
  Lx.3.46<-Lx.3/46
  Ly<-yl[2]-yl[1]
  Ly.55<-Ly/55
  seqx<-xl[1]+c(1,3,6,10,14,18,22,26,30,34,38,42)*Lx.3.46
  text(seqx,yl[2],adj=0,cex=1,c("#","n","sdns","avg","mx","mo","ma","mav","cvrel","bia","rms","made"))
  for (i in 1:length(cra.x)) {
    text(seqx,(yl[2]-Ly.55*i),adj=0,cex=1.,
    c(cra.lab[i],
      round(n.cra[i],0),
         round(sqrt(cra.stn.dens[i])/pi,0),
          round(mean[i],1),
          round(max.x[i],1),
          round(max.yo[i],1),
          round(max.ya[i],1),
          round(max.yav[i],1),
          round(cv.rms.rel[i]*100,0),
          round(cv.bias.sq[i],1),
          round(cv.rmse.sq[i],1),
          round(cv.made.sq[i],1)))
  }
# erroneous 
  points(x[ydqc==1],y[ydqc==1],col="black",pch=4,lwd=3)
# missing
  if (length(yvar1)>0) {
    for (b in 1:length(col1)) {
      cond<-NULL
      if (b==1) cond<-which(yvar1<brk[b] & ydqc==-1)
      if (b==length(col1)) cond<-which(yvar1>=brk[b-1] & ydqc==-1)
      if (b>1 & b<length(col1)) cond<-which( (yvar1>=brk[b-1]) & (yvar1<brk[b]) & ydqc==-1)
      points(x[cond],y[cond],pch=19,cex=0.5,col=col1[b])
      points(x[cond],y[cond],cex=0.6,col=colpoints)
    }
  }
#
  if (!is.null(yvartext)) {
    text(x,y,labels=yvartext,cex=yvartextcex)
  }
#
  legend(x=legpos,legend=legtxt,cex=legcex,fill=c(col1[1],col1[n.col1:2]))
  if (length(mtxt1)>0) mtext(mtxt1,side=3,cex=1.6)
  dev.off()
  return()
}


