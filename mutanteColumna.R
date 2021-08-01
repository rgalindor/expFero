uno <-commandArgs(trailingOnly=T)[1]			
#dos <-commandArgs(trailingOnly=T)[2]			
#tres <-commandArgs(trailingOnly=T)[3]			

####modulo                        #35-(35%/%10)*10 = 5     tipos
####division entera               #35%/%10 = 3             colores


fero<-read.delim(uno,sep=",")
columnas<-length(fero)/10
tiempo<-seq(1,length(fero[,1]),1)
tiempo<-tiempo*(24/length(fero[,1]))
tipos<-c("1 PT1","10 PT1","40 PT1","1 PT2","10 PT2","40 PT2","1 alpha-factor","10 alpha-factor","40 alpha-factor","H2O")
#colores<-rep(0,columnas)
#for(i in 1:columnas){
#	colores[i]<-rgb(1-(i*20/255),1-((1/columnas)*i),(1/columnas)*i) 
#}
colores<-rainbow(columnas)
output_name<-sub(".csv","", uno)
output_name<-paste(output_name, ".pdf", sep="")
pdf(output_name)

k<-0
plot(1,1, type="l", xlim=c(-4,tiempo[length(tiempo)]),ylim=c(0,max(fero)+min(fero)), xlab="Tiempo", ylab="Densidad Optica", main="Celulas que crecieron")
legend(-4, max(fero)+0.1, "[ ](microMol)", bty="n")
legend(4, max(fero)+min(fero)+0.1, "Cepas", bty="n")
for(i in 1:length(fero[0,])){
	if(i<11){
		legend(-4,max(fero)-(i*0.07), tipos[i] , pch=(i-(i%/%10)*10)+1, bty="n")
	}
	if((i-(i%/%10)*10)==0){
		legend(4,max(fero)+min(fero)-(((i%/%10)+1)*0.05), sub("_.+","",names(fero)[i]) ,  lty=1,col=colores[i%/%10], bty="n")
	}
	if(fero[length(fero[,i]),i]>fero[1,i]+.5){
		bandera<-k>0
		lines(tiempo,fero[,i], col=colores[((i+9)%/%10)])
		points(tiempo,fero[,i], col=colores[((i+9)%/%10)], pch=(i-(i%/%10)*10)+1)
		cat(paste(tipos[i-(i%/%10)*10],sub("_.+","",names(fero)[i])), "\n", file="creci.txt", append=bandera)
		k<-k+1
	}
}

k<-0
plot(1,1, type="l",  xlim=c(-4,tiempo[length(tiempo)]),ylim=c(0,max(fero)+min(fero)), xlab="Tiempo", ylab="Densidad Optica", main="Celulas que muestran muerte de algun tipo")
legend(-4, max(fero)+0.1, "[ ](microMol)", bty="n")
legend(4, max(fero)+min(fero)+0.1, "Cepas", bty="n")
for(i in 1:length(fero[0,])){
	if(i<11){
		legend(-4,max(fero)-(i*0.07), tipos[i] , pch=(i-(i%/%10)*10)+1, bty="n")
	}
	if((i-(i%/%10)*10)==0){
		legend(4,max(fero)+min(fero)-(((i%/%10)+1)*0.05), sub("_.+","",names(fero)[i]) ,  lty=1,col=colores[i%/%10], bty="n")
	}
	if(fero[length(fero[,i]),i]<fero[1,i]+.5){
		bandera<-k>0
		lines(tiempo,fero[,i], col=colores[((i+9)%/%10)])
		points(tiempo,fero[,i], col=colores[((i+9)%/%10)], pch=(i-(i%/%10)*10)+1)
		cat(paste(tipos[i-(i%/%10)*10],sub("_.+","",names(fero)[i])), "\n", file="muerte_d.txt", append=bandera)
		k<-k+1
	}
}

k<-0
plot(1,1, type="l",  xlim=c(-4,tiempo[length(tiempo)]),ylim=c(0,max(fero)+min(fero)), xlab="Tiempo", ylab="Densidad Optica", main="Celulas que presentan muerte parcial")
legend(-4, max(fero)+0.1, "[ ](microMol)", bty="n")
legend(4, max(fero)+min(fero)+0.1, "Cepas", bty="n")
for(i in 1:length(fero[0,])){
	if(i<11){
		legend(-4,max(fero)-(i*0.07), tipos[i] , pch=(i-(i%/%10)*10)+1, bty="n")
	}
	if((i-(i%/%10)*10)==0){
		legend(4,max(fero)+min(fero)-(((i%/%10)+1)*0.05), sub("_.+","",names(fero)[i]) ,  lty=1,col=colores[i%/%10], bty="n")
	}
	if((fero[length(fero[,i]),i]<fero[1,i]+.5)&(fero[length(fero[,i]),i]>fero[1,i])){
		bandera<-k>0
		lines(tiempo,fero[,i], col=colores[((i+9)%/%10)])
		points(tiempo,fero[,i], col=colores[((i+9)%/%10)], pch=(i-(i%/%10)*10)+1)
		cat(paste(tipos[i-(i%/%10)*10],sub("_.+","",names(fero)[i])), "\n", file="muerte_parcial.txt", append=bandera)
		k<-k+1
	}
}

k<-0
plot(1,1, type="l",  xlim=c(-4,tiempo[length(tiempo)]),ylim=c(0,max(fero)+min(fero)), xlab="Tiempo", ylab="Densidad Optica", main="Celulas que mueren")
legend(-4, max(fero)+0.1, "[ ](microMol)", bty="n")
legend(4, max(fero)+min(fero)+0.1, "Cepas", bty="n")
for(i in 1:length(fero[0,])){
	if(i<11){
		legend(-4,max(fero)-(i*0.07), tipos[i] , pch=(i-(i%/%10)*10)+1, bty="n")
	}
	if((i-(i%/%10)*10)==0){
		legend(4,max(fero)+min(fero)-(((i%/%10)+1)*0.05), sub("_.+","",names(fero)[i]) ,  lty=1,col=colores[i%/%10], bty="n")
	}
	if(fero[length(fero[,i]),i]<fero[1,i]){
		bandera<-k>0
		lines(tiempo,fero[,i], col=colores[((i+9)%/%10)])
		points(tiempo,fero[,i], col=colores[((i+9)%/%10)], pch=(i-(i%/%10)*10)+1)
		cat(paste(tipos[i-(i%/%10)*10],sub("_.+","",names(fero)[i])), "\n", file="muerte_final.txt", append=bandera)
		k<-k+1
	}
}
#dev.off()
suma_wt<-rep(0,10)
cont<-1
tip_nam<-0
for(j in 1:columnas){
#k<-0
	plot(1,1, type="l",  xlim=c(-4,tiempo[length(tiempo)]),ylim=c(0,max(fero)+min(fero)), xlab="Tiempo", ylab="Densidad Optica", main="Crecimiento de las celulas")
	legend(-4, max(fero)+0.1, "[ ](microMol)", bty="n")
	legend(4, max(fero)+min(fero)+0.1, "Cepas", bty="n")
	for(i in 1:length(fero[0,])){
		if(i<11){
			legend(-4,max(fero)-(i*0.07), tipos[i] , pch=(i-(i%/%10)*10)+1, bty="n")
		}
		if((i-(i%/%10)*10)==0){
			legend(4,max(fero)+min(fero)-(((i%/%10)+1)*0.05), sub("_.+","",names(fero)[i]) ,  lty=1,col=colores[i%/%10], bty="n")
			tip_nam[(i%/%10)]<-sub("_.+","",names(fero)[i])
		}
		if(((i+9)%/%10==j)){
			if(j==columnas){
				for(l in 2:length(fero[,i])){
					suma_wt[cont]<-suma_wt[cont]+(((fero[l,i]+fero[l-1,i])/2)*(24/length(fero[,i])))
				}
				cont<-cont+1
			}
			bandera<-k>0
			lines(tiempo,fero[,i], col=colores[((i+9)%/%10)])
			points(tiempo,fero[,i], col=colores[((i+9)%/%10)], pch=(i-(i%/%10)*10)+1)
	#		cat(paste(tipos[i-(i%/%10)*10],sub("_.+","",names(fero)[i])), "\n", file="muerte_final.txt", append=bandera)
	#		k<-k+1
		}
	}
}
suma<-rep(0,length(fero[0,]))
v_efecto<-suma
for(j in 0:9){
#k<-0
	plot(1,1, type="l",  xlim=c(-4,tiempo[length(tiempo)]),ylim=c(0,max(fero)+min(fero)), xlab="Tiempo", ylab="Densidad Optica", main="Crecimiento de las celulas")
	legend(-4, max(fero)+0.1, "[ ](microMol)", bty="n")
	legend(4, max(fero)+min(fero)+0.1, "Cepas", bty="n")
	for(i in 1:length(fero[0,])){
		if(i<11){
			legend(-4,max(fero)-(i*0.07), tipos[i] , pch=(i-(i%/%10)*10)+1, bty="n")
		}
		if((i-(i%/%10)*10)==0){
			legend(4,max(fero)+min(fero)-(((i%/%10)+1)*0.05), sub("_.+","",names(fero)[i]),lty=1,col=colores[i%/%10], bty="n")
		}
		if(((i+9)-((i+9)%/%10)*10)==j){
			for(l in 2:length(fero[,i])){
				suma[i]<-suma[i]+(((fero[l,i]+fero[l-1,i])/2)*(24/length(fero[,i])))
			}
			v_efecto[i]<-(suma[i]-suma_wt[j+1])/suma_wt[length(suma_wt)]
			bandera<-k>0
			lines(tiempo,fero[,i], col=colores[((i+9)%/%10)])
			points(tiempo,fero[,i], col=colores[((i+9)%/%10)], pch=(i-(i%/%10)*10)+1)
	#		cat(paste(tipos[i-(i%/%10)*10],sub("_.+","",names(fero)[i])), "\n", file="muerte_final.txt", append=bandera)
	#		k<-k+1
		}
	}
}
cosa<-matrix(v_efecto, ncol=columnas)
barplot(cosa, beside=TRUE, col=terrain.colors(10), ylim=c(-1.5,1.5), names=tip_nam, main="Valores estandarizados", ylab="Efecto", xlab="Prueba")
legend(0,1.5, tipos, fill=terrain.colors(10), bty="n", cex=0.75)
dev.off()

