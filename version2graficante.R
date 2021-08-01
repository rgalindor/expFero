uno <-commandArgs(trailingOnly=T)[1]			#Los datos del crecinometro
dos <-commandArgs(trailingOnly=T)[2]			#wt H2O
tres <-commandArgs(trailingOnly=T)[3]			#wt Peptido
cuatro <-commandArgs(trailingOnly=T)[4]		#nombre para grafica 1
cinco <-commandArgs(trailingOnly=T)[5]		#nombre para grafica 2
seis <-commandArgs(trailingOnly=T)[6]			#nombre del archivo output
siete <-commandArgs(trailingOnly=T)[7]

siete<-as.numeric(siete)
opt_dens <- read.delim(uno)
wth <- read.table(dos)
wtp <- read.table(tres)

graph1 <-paste(cuatro,".pdf",sep="")
graph2 <-paste(cinco,".pdf",sep="")

num_pruebas <-length(colnames(opt_dens))-1
num_cepas <-num_pruebas/2

num_mediciones<- length(as.matrix(opt_dens[1]))
area1<-rep(0,num_cepas)
area2<-rep(0,num_cepas)

wth<-as.matrix(wth)
wtp<-as.matrix(wtp)
pta<-0
h2oa<-0

for (i in 2:length(wth)){
	pta<-pta+((wtp[i]+wtp[i-1])/2)
	h2oa<-h2oa+((wth[i]+wth[i-1])/2)
}

#corte<-mean(c(pta,h2oa))
corte<-mean(c(pta,h2oa))*siete
# puedo hacer la de todas las cepas sin que hayan crecido bn
area1<-rep(0,num_cepas)
area2<-area1
nombres<-area1
col_robin<-area1
for(i in 1:num_cepas){
	helper<-as.matrix(opt_dens[(i*2)-1])
	for(j in 2:length(helper)){
		area1[i]<-area1[i]+((helper[j]+helper[j-1])/2)
	}
}
for(i in 1:num_cepas){
	helper<-as.matrix(opt_dens[i*2])
	for(j in 2:length(helper)){
		area2[i]<-area2[i]+((helper[j]+helper[j-1])/2)
	}
}
tiempo<-seq(1,length(helper))
tiempo<-tiempo*(24/length(helper))
pdf(graph1)
k<-0
plot(tiempo, wth, type="lines", col=rgb(0,1,1/4), ylim=c(0,2), ylab="DO", xlab="hours",main="Grafica de cepas que mostraron resistencia", lty=2)#wt+h20
lines(tiempo, wtp, col=rgb(0,1/8,1/10), lty=2)##wt+p
for(i in 1:num_cepas){
	if(mean(c(area1[i],area2[i]))>corte){
####aki poner lo de las graficas
		helper1<-as.matrix(opt_dens[(i*2)-1])
		helper2<-as.matrix(opt_dens[i*2])
		for(j in 1:length(helper1)){
			helper[j]<-mean(c(helper1[j],helper2[j]))
		}
		lines(tiempo, helper, col=rgb(1-(i*5/255),i*5/255,i*10/255))
		points(tiempo, helper, col=rgb(1-(i*5/255),i*5/255,i*10/255), pch=i)
		legend(.2, 1.8-(k*0.1), colnames(opt_dens[(i*2)-1]), col=rgb(1-(i*5/255),i*5/255,i*10/255), lty=1, bty="n", pch=i)
		bandera<-k>0
		cat(colnames(opt_dens[(i*2)-1]),"\n",file=seis, append=bandera)
		k<-k+1
	}
	nombres[i]<-colnames(opt_dens[(i*2)-1])
	col_robin[i]<-rgb(1-(i*5/255),i*5/255,i*10/255)
}
legend(.2,2, "Wild Type H2O", col=rgb(0,1,1/4), lty=2, bty="n")
legend(.2,1.9, "Wild Type PT1", col=rgb(0,1/8,1/10), lty=2, bty="n")
dev.off()
pdf(graph2)
media<-area1*0
for(i in 1:length(area1)){
	media[i]<-mean(c(area1[i],area2[i]))
}
maximo<-max(media)
maximo<-max(maximo,h2oa)
barplot(media, col=col_robin, ylim=c(-1,maximo+3), xlim=c(-4,length(media)+4), main="Area media por cepa", ylab="Area bajo la curva", las=2, names.arg=nombres)
abline(h=h2oa, lty=5, col="gray")
legend(-6,h2oa, "wt+H2O", bty="n")
abline(h=corte, lty=2, col=rgb(0,1,1))
legend(-6,corte, "cut point", bty="n")
abline(h=pta, lty=2, col="gray")
legend(-6,pta, "wt+PT1", bty="n")
dev.off()
