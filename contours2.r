#Programa para calcular un mapa de contorno a partir de un archivo con 3 columnas x,y,z. 
######################################################################
#DATOS A INTRODUCIR
colpal=rainbow #color
transpuesta=0 #Poner a 1 si la primera columna permancece constante mientras se varía la segunda
datafile=paste('SEP_HCN_10lam.dat',sep="") #titulo del archivo de datos a leer
#pngtit=paste('pot.png',sep="")#titulo del archivo gráfico de salida
xlimite=c(0,180)
ylimite=c(4,12)
zlimites<-c(-30.02657,0)
xetiq=expression(theta)
yetiq='R'
npoint=101
#########################################################################
#png(file=paste(pngtit,sep=""))
z1dat <- as.matrix(read.table(datafile))
z1dat1=as.numeric(z1dat[,2])
z1dat2=as.numeric(z1dat[,1])
z1dat3=as.numeric(z1dat[,3])
zvar<-array(z1dat3,dim=c(npoint,npoint))


if (transpuesta==1) zvar<-t(zvar)

par(xaxs="i",yaxs="i")
plot(1,1,xlab=xetiq,ylab=yetiq,xlim=xlimite,ylim=ylimite,type="n") #type=n es que va a dibujar ejes y va a dibujar th y r
par(new=TRUE,ann=F,xaxs="i",yaxs="i") 

image(y = seq(min(ylimite), max(ylimite), length.out=npoint),x =  seq(min(xlimite), max(xlimite)
     , length.out=npoint),z=zvar,col=colpal(100),zlim=zlimites,axes=FALSE)
par(new=TRUE)
contour(y = seq(min(ylimite), max(ylimite), length.out=npoint),x =  seq(min(xlimite), max(xlimite)
        , length.out=npoint),z=zvar,zlim=zlimites,axes=FALSE,nlevels=10,lwd=1,labcex = 1,vfont=c("sans serif","bold"))
#dev.off()

