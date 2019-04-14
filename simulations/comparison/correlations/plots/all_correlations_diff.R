library(ggplot2)
library(viridis)

    9-md_1000ns_2_A  9-md_1000ns_3_A  corr2dat  do_corrA.sh  do_corr.sh  lllogA  llog   llogB          traj_files
9-md_1000ns_2  9-md_1000ns_3    9-md_1000ns_A

xA = read.table("../100ns/9-md_1000ns_A/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
#
setEPS()
postscript("LMI_100ns.dim1A.eps")
ggplot(data = melt(xa1), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400")) + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

setEPS()
postscript("LMI_100ns.dim1B.eps")
ggplot(data = melt(xb1), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400")) + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

setEPS()
postscript("dif_dim1A.dim1B.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()
######################################################
xA = read.table("../100ns/9-md_1000ns_2_A/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_2/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
#
setEPS()
postscript("LMI_100ns.dim2A.eps")
ggplot(data = melt(xa1), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400")) + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

setEPS()
postscript("LMI_100ns.dim2B.eps")
ggplot(data = melt(xb1), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400")) + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

setEPS()
postscript("dif_dim2A.dim2B.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()
###########################################
xA = read.table("../100ns/9-md_1000ns_3_A/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_3/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
#
setEPS()
postscript("LMI_100ns.dim3A.eps")
ggplot(data = melt(xa1), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400")) + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

setEPS()
postscript("LMI_100ns.dim3B.eps")
ggplot(data = melt(xb1), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400")) + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

setEPS()
postscript("dif_dim3A.dim3B.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

###############################
xA = read.table("../100ns/9-md_1000ns_1_A/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_2_A/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
setEPS()
postscript("dif_dim1A.dim2A.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

###################################
xA = read.table("../100ns/9-md_1000ns_1_A/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_3_A/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
setEPS()
postscript("dif_dim1A.dim3A.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

###################################

xA = read.table("../100ns/9-md_1000ns_2_A/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_3_A/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
setEPS()
postscript("dif_dim2A.dim3A.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

###################################
xA = read.table("../100ns/9-md_1000ns_1/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_2/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
setEPS()
postscript("dif_dim1B.dim2B.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

###################################
xA = read.table("../100ns/9-md_1000ns_1/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_3/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
setEPS()
postscript("dif_dim1B.dim3B.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

###################################
xA = read.table("../100ns/9-md_1000ns_2/corr-average.dat", header = FALSE, fill = TRUE)
xB = read.table("../100ns/9-md_1000ns_3/corr-average.dat", header = FALSE, fill = TRUE)

colnames(xA) <- gsub('V', '', colnames(xA), fixed=TRUE)
xA
xa1 <- as.matrix(xA)
xb1 <- as.matrix(xB)
is.matrix(xa1)
xa1
#
sub = matrix(NA, 413,413)
for (i in 1:nrow(xa1)) #(nrow(z)==nrow(c))
{
  sub [i,]= xa1[i,]-xb1[i,]
}
sub
rng = range(c((xa1), (xb1))) #a range to have the same min and max for both plots
library(reshape2)
setEPS()
postscript("dif_dim2B.dim3B.eps")
ggplot(data = melt(sub), aes(as.factor(Var1), as.factor(Var2), fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "magma", direction = -1, limits=c(-1, 1)) + theme_bw() + scale_y_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + scale_x_discrete(breaks= c("10","50","100","150","200","250","300","350","400"))  + xlab("residue") + ylab("residue") + labs(fill='LMI Correlation Value')
dev.off()

###################################

