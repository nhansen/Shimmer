current_2_datasets <- c("set2_1", "set2_2", "set2_3", "set2_4", "set2_5", "set2_6", "set2_7", "set2_8", "set2_9", "set2_10", "set2_11", "set2_12", "set2_13", "set2_14", "set2_15", "set2_16", "set2_17", "set2_18", "set2_19", "set2_20")
current_10_datasets <- c("set10_1", "set10_2", "set10_3", "set10_4", "set10_5", "set10_6", "set10_7", "set10_8", "set10_9", "set10_10")
jsnvmix_10_datasets <- c("set10_3", "set10_4", "set10_5", "set10_7", "set10_8", "set10_9")
jsnvmix_2_datasets <- c("set2_1", "set2_2", "set2_5", "set2_13", "set2_10", "set2_4", "set2_8", "set2_9", "set2_11", "set2_12", "set2_16", "set2_14", "set2_15", "set2_17", "set2_18", "set2_19")
plotdir <- function(outdir, purity) {
	shimmer <- read.table(paste(outdir, "/shimmer_sens_spec.txt", sep=""), header=FALSE, sep=" ")
	sniper <- read.table(paste(outdir, "/sniper_sens_spec.txt", sep=""), header=FALSE, sep=" ")
	varscan <- read.table(paste(outdir, "/varscan_sens_spec.txt", sep=""), header=FALSE, sep=" ")
	deep <- read.table(paste(outdir, "/deepsnv_sens_spec.txt", sep=""), header=FALSE, sep=" ")
	shimmer_tp <- shimmer[,"V2"]
	shimmer_fp <- shimmer[,"V3"]
	shimmer_fn <- shimmer[,"V4"]
	deep_tp <- deep[,"V2"]
	deep_fp <- deep[,"V3"]
	deep_fn <- deep[,"V4"]
	varscan_tp <- varscan[,"V2"]
	varscan_fp <- varscan[,"V3"]
	varscan_fn <- varscan[,"V4"]
	sniper_tp <- sniper[,"V2"]
	sniper_fp <- sniper[,"V3"]
	sniper_fn <- sniper[,"V4"]
	alldata <- rbind(shimmer, varscan, sniper, deep)
	minmax_x <- range(alldata[,"V3"])	
	minmax_y <- range(alldata[,"V2"]/(alldata[,"V2"]+alldata[,"V4"]))
	plot(shimmer_fp, shimmer_tp/(shimmer_tp+shimmer_fn), col=c("red"), xlim=minmax_x, ylim=minmax_y, main=c(purity, " Percent Tumor"), type="o", lty=1, pch=16)
	points(deep_fp, deep_tp/(deep_tp+deep_fn), col=c("darkblue"), type="o", lty=2, pch=16)
	points(sniper_fp, sniper_tp/(sniper_tp+sniper_fn), col=c("darkgreen"), type="o", lty=3, pch=16)
	points(varscan_fp, varscan_tp/(varscan_tp+varscan_fn), col=c("purple"), type="o", lty=4, pch=16)
}

read_average_frame <- function(files) {
	catpipe <- pipe(paste("cat ", paste(files, collapse = ' ')))
	df <- read.table(catpipe, header=FALSE, sep=" ")
	avgtp <- tapply(df[,2], df["V1"], mean)
	avgfp <- tapply(df[,3], df["V1"], mean)
	avgfn <- tapply(df[,4], df["V1"], mean)
	avgdf <- t(rbind(avgtp, avgfp, avgfn))
	return(avgdf)
}

read_error_frame <- function(files) {
	catpipe <- pipe(paste("cat ", paste(files, collapse = ' ')))
	df <- read.table(catpipe, header=FALSE, sep=" ")
	stderr <- function(x) sqrt(sd(x)/length(x))
	errtp <- tapply(df[,2], df["V1"], stderr)
	errfp <- tapply(df[,3], df["V1"], stderr)
	errfn <- tapply(df[,4], df["V1"], stderr)
	errdf <- t(rbind(errtp, errfp, errfn))
	return(errdf)
}

multiplotdir <- function(outdirs, title='', titlesize=1, text='', textpos=c(0,0), legendpos=c(0,10.0), minmax_x=c(-1,-1), minmax_y=c(-1, -1), no_ylab=0, plotlabel='') {
	shimfiles <- lapply(outdirs, function(x) paste(x, "/shimmer_sens_spec.txt", sep=""))
	avgshim <- read_average_frame(shimfiles)
	shimmer_tp <- avgshim[,1]
	shimmer_fp <- avgshim[,2]
	shimmer_fn <- avgshim[,3]
	deepfiles <- lapply(outdirs, function(x) paste(x, "/deepsnv_sens_spec.txt", sep=""))
	avgdeep <- read_average_frame(deepfiles)
	deep_tp <- avgdeep[,1]
	deep_fp <- avgdeep[,2]
	deep_fn <- avgdeep[,3]
	sniperfiles <- lapply(outdirs, function(x) paste(x, "/sniper_sens_spec.txt", sep=""))
	avgsniper <- read_average_frame(sniperfiles)
	sniper_tp <- avgsniper[,1]
	sniper_fp <- avgsniper[,2]
	sniper_fn <- avgsniper[,3]
	varscanfiles <- lapply(outdirs, function(x) paste(x, "/varscan_sens_spec.txt", sep=""))
	avgvarscan <- read_average_frame(varscanfiles)
	varscan_tp <- avgvarscan[,1]
	varscan_fp <- avgvarscan[,2]
	varscan_fn <- avgvarscan[,3]
	alldata <- rbind(avgshim, avgsniper, avgdeep, avgvarscan)
	if (minmax_x[[1]] < 0) {
		minmax_x <- range(alldata[,"avgfp"])	
	}
	if (minmax_y[[1]] < 0) {
		minmax_y <- range(alldata[,"avgtp"]/(alldata[,"avgtp"]+alldata[,"avgfn"]))
	}
	names <- c("Shimmer", "deepSNV", "SomaticSniper", "VarScan2")
	linetypes <- c(1,2,3,4)
	pointtypes <- c(1,2,3,4)
	pointcolors <- c("black", "blue", "darkgreen", "darkred")
	if (no_ylab == 1) {
		ylabel = c("")
        }
	else {
		ylabel = c("Sensitivity")
	}
	plot(shimmer_fp, shimmer_tp/(shimmer_tp+shimmer_fn), col=pointcolors[[1]], xlab=c("False Positives"), ylab=ylabel, cex.lab=1.4, xlim=minmax_x, ylim=minmax_y, type="o", lty=linetypes[[1]], pch=pointtypes[[1]])
	title(main=title, cex.main=titlesize)
	points(deep_fp, deep_tp/(deep_tp+deep_fn), col=pointcolors[2], type="o", lty=2, pch=2)
	points(sniper_fp, sniper_tp/(sniper_tp+sniper_fn), col=pointcolors[3], type="o", lty=3, pch=3)
	points(varscan_fp, varscan_tp/(varscan_tp+varscan_fn), col=pointcolors[4], type="o", lty=4, pch=4)
	legend(x=legendpos[[1]], y=legendpos[[2]], names, lty=linetypes, pch=pointtypes, col=pointcolors, cex=1.2)
	text_xpos <- minmax_x[[1]] + textpos[[1]]*(minmax_x[[2]] - minmax_x[[1]])
	text_ypos <- minmax_y[[1]] + textpos[[2]]*(minmax_y[[2]] - minmax_y[[1]])
	text(x=text_xpos, y=text_ypos, labels=text, cex=1.3)
	par(xpd=T)
	text(-0.15*(par("usr")[2]-par("usr")[1]),par("usr")[4]+0.08*(par("usr")[4]-par("usr")[3]),plotlabel,cex=1.3)
	par(xpd=F)
}

filelist <- function(set, purity) {
	allfiles <- list.files("/cluster/ifs/projects/seqSim/shimmerpaper/newsimplots")
	txtfiles <- allfiles[grep(glob2rx(paste("set",set,"_*",purity, sep="")), allfiles)]
	return(txtfiles)
}

simfilelist <- function(set, purity) {
	allfiles <- list.files("/cluster/ifs/projects/seqSim/shimmerpaper/newsimplots")
	txtfiles <- allfiles[grep(glob2rx(paste("sim",set,"_*",purity, sep="")), allfiles)]
	return(txtfiles)
}

allfilelist <- function(set, purity, extension='') {
	allfiles <- list.files("/cluster/ifs/projects/seqSim/shimmerpaper/newsimplots")
	globstring <- paste(set,"_",extension,"*", purity, sep="")
	txtfiles <- allfiles[grep(glob2rx(globstring), allfiles)]
	return(txtfiles)
}

mean_frame <- function(set, purity, prog, extension='') {
	dirs <- unlist(lapply(set, function(x) allfilelist(x, purity, extension)))
	topdir <- "/cluster/ifs/projects/seqSim/shimmerpaper/newsimplots"
	progfiles <- lapply(dirs, function(x) paste(topdir,"/",x, "/",prog,"_sens_spec.txt", sep=""))
	meanframe <- read_average_frame(progfiles)
	return(meanframe)
}

error_frame <- function(set, purity, prog, extension='') {
	dirs <- unlist(lapply(set, function(x) allfilelist(x, purity, extension)))
	topdir <- "/cluster/ifs/projects/seqSim/shimmerpaper/newsimplots"
	progfiles <- lapply(dirs, function(x) paste(topdir,"/",x, "/",prog,"_sens_spec.txt", sep=""))
	errorframe <- read_error_frame(progfiles)
	return(errorframe)
}

plotsets <- function(setnames) {
	par(mfrow=c(1,3))
	files100 <- unlist(lapply(setnames, function(x) filelist(x, 100)))
	multiplotdir(files100)
	files60 <- unlist(lapply(setnames, function(x) filelist(x, 60)))
	multiplotdir(files60)
	files20 <- unlist(lapply(setnames, function(x) filelist(x, 20)))
	multiplotdir(files20)
}

errsens <- function(setnames, purity, prog, extension='') {
	meandf <- mean_frame(setnames, purity, prog, extension)
	errdf <- error_frame(setnames, purity, prog, extension)
	senserr <- errdf[,"errtp"]/meandf[,"avgtp"] + (errdf[,"errtp"]+errdf[,"errfn"])/(meandf[,"avgtp"] + meandf[,"avgfn"])
	senserr <- senserr*meandf[,"avgtp"]/(meandf[,"avgtp"] + meandf[,"avgfn"])
	return(senserr)
}

meansens <- function(setnames, purity, prog, extension='') {
	meandf <- mean_frame(setnames, purity, prog, extension)
	sensmean <- meandf[,"avgtp"]/(meandf[,"avgtp"] + meandf[,"avgfn"])
	return(sensmean)
}

errfp <- function(setnames, purity, prog, extension='') {
	errdf <- error_frame(setnames, purity, prog, extension)
	fperr <- errdf[,"errfp"]
	return(fperr)
}

meanfp <- function(setnames, purity, prog, extension='') {
	meandf <- mean_frame(setnames, purity, prog, extension)
	fpmean <- meandf[,"avgfp"]
	return(fpmean)
}

plotsimsets <- function(setnames) {
	par(mfrow=c(1,3))
	files100 <- unlist(lapply(setnames, function(x) simfilelist(x, 100)))
	multiplotdir(files100)
	files60 <- unlist(lapply(setnames, function(x) simfilelist(x, 60)))
	multiplotdir(files60)
	files20 <- unlist(lapply(setnames, function(x) simfilelist(x, 20)))
	multiplotdir(files20)
}

publicationplot <- function(set1names, set2names,extension='') {
	par(mfrow=c(2,3))
	files1_100 <- unlist(lapply(set1names, function(x) allfilelist(x, 100, extension)))
	multiplotdir(files1_100, text='100% purity', textpos=c(0.6, 0.15), minmax_y=c(0.725, 0.765), plotlabel=expression(bold('A')))
	files1_60 <- unlist(lapply(set1names, function(x) allfilelist(x, 60, extension)))
	multiplotdir(files1_60, title='10 Mutations/Mb', titlesize=1.5, text='60% purity', textpos=c(0.6, 0.15), minmax_y=c(0.69,0.75), no_ylab=1, plotlabel=expression(bold('B')))
	files1_20 <- unlist(lapply(set1names, function(x) allfilelist(x, 20, extension)))
	multiplotdir(files1_20, text='20% purity', textpos=c(0.6, 0.15), minmax_y=c(0.0, 0.6), no_ylab=1, plotlabel=expression(bold('C')))
	files2_100 <- unlist(lapply(set2names, function(x) allfilelist(x, 100, extension)))
	multiplotdir(files2_100, text='100% purity', textpos=c(0.6, 0.15), minmax_y=c(0.725, 0.765), plotlabel=expression(bold('D')))
	files2_60 <- unlist(lapply(set2names, function(x) allfilelist(x, 60, extension)))
	multiplotdir(files2_60, title='2 Mutations/Mb', titlesize=1.5, text='60% purity', textpos=c(0.6, 0.15), minmax_y=c(0.69, 0.75), no_ylab=1, plotlabel=expression(bold('E')))
	files2_20 <- unlist(lapply(set2names, function(x) allfilelist(x, 20, extension)))
	multiplotdir(files2_20, text='20% purity', textpos=c(0.6, 0.15), legendpos=c(45, 0.5), minmax_y=c(0.0, 0.6), no_ylab=1, plotlabel=expression(bold('F')))
}

