plotGATKStat <- function(gatk, QD=2, FS=60, SOR=3, MQ=40, MQRS=c(-2.5, 2.5), RPRS=c(-4, 4), txt_cex=0.8, exp_id = NULL) {
  mbool <- NULL
  tmp <- gatk$QD < QD
  tmp[is.na(tmp)] <- FALSE
  mbool <- cbind(mbool, tmp)
  tmp <- gatk$FS > FS
  tmp[is.na(tmp)] <- FALSE
  mbool <- cbind(mbool, tmp)
  tmp <- gatk$SOR > SOR
  tmp[is.na(tmp)] <- FALSE
  mbool <- cbind(mbool, tmp)
  tmp <- gatk$MQ < MQ
  tmp[is.na(tmp)] <- FALSE
  mbool <- cbind(mbool, tmp)
  tmp <- gatk$MQRankSum < MQRS[1] | gatk$MQRankSum > MQRS[2]
  tmp[is.na(tmp)] <- FALSE
  mbool <- cbind(mbool, tmp)
  tmp <- gatk$ReadPosRankSum < RPRS[1] | gatk$ReadPosRankSum > RPRS[2]
  tmp[is.na(tmp)] <- FALSE
  mbool <- cbind(mbool, tmp)
  
  cstat <- apply(mbool, 2, sum)
  csnp  <- apply(mbool, 1, sum)
  ntot  <- nrow(gatk)
  nkeep <- sum(csnp == 0)
  ndel  <- ntot - nkeep
  pkeep <- nkeep / ntot
  pdel  <- ndel / ntot
  
  par(mfrow=c(3,2),mar=c(4.1, 4.1, 2.1, 0.1))
  
  col_tot <- grey(0.5)
  col_keep <- rgb(0.2,0.2,0.8,0.5)
  col_del  <- rgb(0.8,0.2,0.2,0.5)
  txt_keep <- rgb(0.2,0.2,0.8)
  txt_del  <- rgb(0.8,0.2,0.2)
  col_sel  <- rgb(0.2,0.8,0.2)
  
  
  # Plot Quality By Depth stat
  #den_tot  <- density(gatk$QD, na.rm = TRUE, from=0, to=40, bw="SJ", kernel="g")
  den_keep <- density(gatk$QD[csnp == 0], na.rm = TRUE, from=0, to=40, bw="SJ", kernel="g")
  den_del  <- density(gatk$QD[csnp != 0], na.rm = TRUE, from=0, to=40, bw="SJ", kernel="g")
  den_tot  <- list(
    x = den_keep$x,
    y = pkeep * den_keep$y + pdel * den_del$y
  ) 
  plot(0, 0, type='n', xlab="QualByDepth (QD)", axes=FALSE, ylab="Density",
       xlim=c(0,40), ylim=c(0, max(den_tot$y)))
  polygon(c(den_tot$x,40,0), c(den_tot$y,0,0), border = NA, col=col_tot)
  polygon(c(den_keep$x,40,0), c(pkeep * den_keep$y, 0,0), border = NA, col=col_keep)
  polygon(c(den_del$x,40,0), c(pdel * den_del$y, 0,0), border = NA, col=col_del)
  axis(1, pos=0)
  axis(2, pos=0)
  polygon(x=c(0,0,QD,QD), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=QD, col.axis="red", col.ticks = col_sel, pos=0)
  mtext(paste("QD kept SNPs: ", ntot - cstat[1], sep=""), col="black", side=1, line = 1.8, at = 10, cex=txt_cex)
  mtext(paste("QD deleted SNPs: ", cstat[1], sep=""), col=col_sel, side=1, line = 1.8, at = 30, cex=txt_cex)
  mtext(paste("Total kept SNPs: ", nkeep, sep=""), col=txt_keep, side = 3, line=0)
  if( !is.null(exp_id) ) {
    text(5, par("usr")[4]*0.9, exp_id, col=grey(0.2,0.5), cex=1.2, adj = 0)
  }
  
  den_keep <- density(log10(gatk$FS[csnp == 0 & gatk$FS != 0]), na.rm = TRUE, from=-1, to=2.8, bw="SJ", kernel="g")
  den_del  <- density(log10(gatk$FS[csnp != 0 & gatk$FS != 0]), na.rm = TRUE, from=-1, to=2.8, bw="SJ", kernel="g")
  den_tot  <- list(
    x = den_keep$x,
    y = pkeep * den_keep$y + pdel * den_del$y
  ) 
  plot(0, 0, type='n', xlab="FisherStrand (FS)", axes=FALSE, ylab="Density",
       xlim=c(-1,2.8), ylim=c(0, max(den_tot$y)))
  polygon(c(den_tot$x,max(den_tot$x),min(den_tot$x)), c(den_tot$y,0,0), border = NA, col=col_tot)
  polygon(c(den_keep$x,max(den_tot$x),min(den_tot$x)), c(pkeep * den_keep$y, 0,0), border = NA, col=col_keep)
  polygon(c(den_del$x,max(den_tot$x),min(den_tot$x)), c(pdel * den_del$y, 0,0), border = NA, col=col_del)
  axis(1, pos=0, at=c(-1,0,1,2),labels = c(0.1, 1, 10, 100))
  axis(2, pos=-1)
  polygon(x=c(log10(FS),log10(FS),2.8,2.8), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=log10(FS), labels = FS, col.axis="red", col.ticks = col_sel, pos=0)
  mtext(paste("FS kept SNPs: ", ntot - cstat[2], sep=""), col="black", side=1, line = 1.8, at = -0.2, cex=txt_cex)
  mtext(paste("FS deleted SNPs: ", cstat[2], sep=""), col=col_sel, side=1, line = 1.8, at = 1.7, cex=txt_cex)
  mtext(paste("Total deleted SNPs: ", ndel, sep=""), col=txt_del, side = 3, line=0)
  
  # Plot SOR
  den_keep <- density(gatk$SOR[csnp == 0], na.rm = TRUE, from=0, to=6, bw="SJ", kernel="g")
  den_del  <- density(gatk$SOR[csnp != 0], na.rm = TRUE, from=0, to=6, bw="SJ", kernel="g")
  den_tot  <- list(
    x = den_keep$x,
    y = pkeep * den_keep$y + pdel * den_del$y
  ) 
  plot(0, 0, type='n', xlab="StrandOddsRatio (SOR)", axes=FALSE, ylab="Density",
       xlim=c(0,6), ylim=c(0, max(den_tot$y)))
  polygon(c(den_tot$x,6,0), c(den_tot$y,0,0), border = NA, col=col_tot)
  polygon(c(den_keep$x,6,0), c(pkeep * den_keep$y, 0,0), border = NA, col=col_keep)
  polygon(c(den_del$x,6,0), c(pdel * den_del$y, 0,0), border = NA, col=col_del)
  axis(1, pos=0)
  axis(2, pos=0)
  #segments(x0=QD, x1=QD, y0=0, y1=max(den_tot$y), col="red", lwd=1.5, lty=2)
  polygon(x=c(SOR,SOR,6,6), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=SOR, col.axis="red", col.ticks = col_sel, pos=0)
  mtext(paste("SOR kept SNPs: ", ntot - cstat[3], sep=""), col="black", side=1, line = 1.8, at = 1.5, cex=txt_cex)
  mtext(paste("SOR deleted SNPs: ", cstat[3], sep=""), col=col_sel, side=1, line = 1.8, at = 4.5, cex=txt_cex)
  
  # Plot MQ
  den_keep <- density(gatk$MQ[csnp == 0], na.rm = TRUE, from=20, to=60, kernel="g")
  den_del  <- density(gatk$MQ[csnp != 0], na.rm = TRUE, from=20, to=60, kernel="g")
  den_tot  <- list(
    x = den_keep$x,
    y = pkeep * den_keep$y + pdel * den_del$y
  ) 
  plot(0, 0, type='n', xlab="RMSMappingQuality (MQ)", axes=FALSE, ylab="Density",
       xlim=c(20,60), ylim=c(0, max(den_tot$y)))
  polygon(c(den_tot$x,60,20), c(den_tot$y,0,0), border = NA, col=col_tot)
  polygon(c(den_keep$x,60,20), c(pkeep * den_keep$y, 0,0), border = NA, col=col_keep)
  polygon(c(den_del$x,60,20), c(pdel * den_del$y, 0,0), border = NA, col=col_del)
  axis(1, pos=0)
  axis(2, pos=20)
  polygon(x=c(20,20,MQ,MQ), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=MQ, col.axis="red", col.ticks = col_sel, pos=0)
  mtext(paste("MQ kept SNPs: ", ntot - cstat[4], sep=""), col="black", side=1, line = 1.8, at = 30, cex=txt_cex)
  mtext(paste("MQ deleted SNPs: ", cstat[4], sep=""), col=col_sel, side=1, line = 1.8, at = 50, cex=txt_cex)
  
  # Plot MQRankSum
  den_keep <- density(gatk$MQRankSum[csnp == 0], na.rm = TRUE, from=-4.5, to=4.5, kernel="g")
  den_del  <- density(gatk$MQRankSum[csnp != 0], na.rm = TRUE, from=-4.5, to=4.5, kernel="g")
  den_tot  <- list(
    x = den_keep$x,
    y = pkeep * den_keep$y + pdel * den_del$y
  ) 
  plot(0, 0, type='n', xlab="MappingQualityRankSumTest (MQRankSum)", axes=FALSE, ylab="Density",
       xlim=c(-4.5,4.5), ylim=c(0, max(den_tot$y)))
  polygon(c(den_tot$x,4.5,-4.5), c(den_tot$y,0,0), border = NA, col=col_tot)
  polygon(c(den_keep$x,4.5,-4.5), c(pkeep * den_keep$y, 0,0), border = NA, col=col_keep)
  polygon(c(den_del$x,4.5,-4.5), c(pdel * den_del$y, 0,0), border = NA, col=col_del)
  axis(1, pos=0)
  axis(2, pos=-4.5)
  polygon(x=c(-4.5,-4.5,MQRS[1],MQRS[1]), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=MQRS[1], col.axis="red", col.ticks = col_sel, pos=0)
  polygon(x=c(MQRS[2],MQRS[2], 4.5, 4.5), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=MQRS[2], col.axis="red", col.ticks = col_sel, pos=0)
  mtext(paste("MQRS kept SNPs: ", ntot - cstat[5], sep=""), col="black", side=1, line = 1.8, at = -2.5, cex=txt_cex)
  mtext(paste("MQRS deleted SNPs: ", cstat[5], sep=""), col=col_sel, side=1, line = 1.8, at = 2.5, cex=txt_cex)
  
  den_keep <- density(gatk$ReadPosRankSum[csnp == 0], na.rm = TRUE, from=-8, to=8, kernel="g")
  den_del  <- density(gatk$ReadPosRankSum[csnp != 0], na.rm = TRUE, from=-8, to=8, kernel="g")
  den_tot  <- list(
    x = den_keep$x,
    y = pkeep * den_keep$y + pdel * den_del$y
  ) 
  plot(0, 0, type='n', xlab="ReadPosRankSumTest (ReadPosRankSum)", axes=FALSE, ylab="Density",
       xlim=c(-8,8), ylim=c(0, max(den_tot$y)))
  polygon(c(den_tot$x,8,-8), c(den_tot$y,0,0), border = NA, col=col_tot)
  polygon(c(den_keep$x,8,-8), c(pkeep * den_keep$y, 0,0), border = NA, col=col_keep)
  polygon(c(den_del$x,8,-8), c(pdel * den_del$y, 0,0), border = NA, col=col_del)
  axis(1, pos=0)
  axis(2, pos=-8)
  #segments(x0=QD, x1=QD, y0=0, y1=max(den_tot$y), col="red", lwd=1.5, lty=2)
  polygon(x=c(-8,-8,RPRS[1],RPRS[1]), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=RPRS[1], col.axis="red", col.ticks = col_sel, pos=0)
  polygon(x=c(RPRS[2],RPRS[2], 8, 8), y=c(0,par("usr")[4], par("usr")[4],0), col=col_sel, density = 5)
  axis(1, at=RPRS[2], col.axis="red", col.ticks = col_sel, pos=0)
  mtext(paste("RPRS kept SNPs: ", ntot - cstat[6], sep=""), col="black", side=1, line = 1.8, at = -4, cex=txt_cex)
  mtext(paste("RPRS deleted SNPs: ", cstat[6], sep=""), col=col_sel, side=1, line = 1.8, at = 4, cex=txt_cex)
  
  return(csnp)
}