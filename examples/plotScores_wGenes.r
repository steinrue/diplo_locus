#This sample script plot scan scores with tracks of filters and gene
##Load required packages
library(data.table)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(latex2exp)


##plot the stat scores of the designated region
### use default palette of blue --> red
plotManChrs <- function(dt, left, right, statname, colorstat,
                mask=FALSE, bed=NULL, color_midpoint=0, palname='Greys'){
  # colnames(dt) <- c('Chr','physPos','stat')
  p <- ggplot(dt)
  chr = unique(dt$Chr)[1]
  p <- p  + geom_point(aes(x=physPos, y=stat, color=colorscale), size=3) + xlab(paste('Positions on chromosome', chr))
  ###adding mask (bedfile: chr start end feature)
  if(mask){
    if(nrow(bed)>0){
      yMin=min(dt$stat)
      width = (max(dt$stat)-yMin)/15
      ntype = length(unique(bed$type))
    if(ntype < 3){
      mycols = c("#a3a3a3","#939393")[1:ntype]
    }else{
      mycols = brewer.pal(ntype+1, palname)[2:(ntype+1)]
      names(mycols) = sort(unique(bed$type))
    }
    p <- p + geom_rect(data=bed, aes(xmin=start, xmax=end, ymin=yMin-width, ymax=yMin, fill=type) )+ 
      scale_fill_manual('Filters:', values = mycols)  + 
      theme(legend.position = 'top', legend.direction = 'horizontal', legend.title=element_text(size=15) )
    }
  }
  ### add colorscale
  p <- p + scale_color_gradient2(name=colorstat, high="#a50026", mid="#ffffbf", low="#313695", midpoint=color_midpoint, na.value="#969696")
  ###finishing touches
  ####statname can be customized (instead of a commandline argument) e.g. statname=expression(italic('B')[2])
  p <- p + ylab(statname) + xlim(left, right) + 
    theme(panel.border=element_blank(), panel.grid.minor.y=element_blank(), panel.background=element_blank()
      )
  return(p)
}


##read gene information
getGeneBed <- function(chr, left, right, refSeqFile="hg19_ncbiRefSeq-All_2023May_allFields.txt"){
    genes = data.table(read.table(refSeqFile, sep="\t", comment="#", header=FALSE)) #, header=1
    # header:
    # #bin name  chrom strand  txStart txEnd cdsStart  cdsEnd  exonCount exonStarts  exonEnds  score name2 cdsStartStat  cdsEndStat  exonFrames
    # colnames(genes)=c('Chr','start','end','cdStart','cdEnd','ID','name','exonStarts','exonEnds')
    colnames(genes) <- c('bin', "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames" )
    genes = genes[,.(bin, Chr=chrom, start=txStart, end=txEnd, strand, cdStart=cdsStart, cdEnd=cdsEnd,
                      ID=name, name=name2, exonCount, exonStarts, exonEnds)]
    genes <- genes[Chr == paste("chr", chr, sep="")]
    # starts and ends are in order (starts always < ends) regardless of strands
    print(dim(genes))
    return(genes[ start <= right & end >= left])
}

##plot gene tracks
plotFancyGenes <- function(genes, noRNA=FALSE){
  if(noRNA==TRUE){
    genes <- genes[!grepl("^NR", ID)]
  }
  dups <- genes[,.N, by=name]
  print(dups)
  print(dups[order(-N)])
  ###if any gene in the region has >3 transcript, only show 3
  if (max(dups$N)>3){
    ### make each gene name unique & choose <d>2<\d>1 transcripts by length <-- ok I regret it; only keep 1, otherwise too messy
    genes =  genes[,.(start=start[which.max(end-start)], # c(, which.median(end-start) , which.min(end-start))
                      end=end[which.max(end-start)], # c(, which.min(end-start)), which.median(end-start) 
                      cdStart=cdStart[which.max(end-start)], # c(, which.min(end-start)), which.median(end-start) 
                      cdEnd=cdEnd[which.max(end-start)], # c(, which.min(end-start)), which.median(end-start) 
                      exonStarts=exonStarts[which.max(end-start)], # c(, which.min(end-start)), which.median(end-start) 
                      exonEnds=exonEnds[which.max(end-start)] #c(, which.min(end-start)) , which.median(end-start)
                    ), 
                  by=name]
  }
  
  ###sort the dataframe by starting position
  # genes=genes[order(start)]
  ###initialize parsing
  ends <- c(0)
  Exons = data.frame(layer=integer(),Estart=numeric(),Eend=numeric())
  Tracks = data.frame(layer=integer(), start=numeric(), end=numeric())
  Names = data.frame(layer=integer(), pos=numeric(), name=character())
  ###check to see if two transcripts are too close and parse the exons
  for(i in 1:nrow(genes)){
    start = genes$start[i]
    end = genes$end[i]
    ###determine the layer to lay the gene track
    layer=0
    for(j in 1:length(ends)){
      if(ends[j] + 1000 <= start){
        layer=j
        ends[j] = end
        break
      }
    }
    if (layer == 0){
      ends = c(ends, end)
      layer = length(ends)
    }
    track = data.frame(layer=layer, start=start, end=end)
    Name = data.frame(layer=layer, pos=0.5*(start+end), name=genes$name[i])
    Tracks = rbind(Tracks, track)
    Names = rbind(Names,Name)
    ###parse the exon start and end positions
    exonS = as.integer(unlist(strsplit(as.character(genes$exonStarts[i]),'[,]')))
    exonE = as.integer(unlist(strsplit(as.character(genes$exonEnds[i]),'[,]')))

    exons = data.frame(layer=rep(layer, length(exonS)), Estart=exonS, Eend=exonE)

    Exons = rbind(Exons, exons)
  }
  ###start plotting
  p <- ggplot() + geom_rect(data=Tracks, aes(xmin=start, xmax=end, ymin=-10*layer+5.5, ymax=-10*layer+6.5), fill='#666666')
  p <- p +geom_rect(data=Exons, aes(xmin=Estart, xmax=Eend, ymin=-10*layer+4, ymax=-10*layer+8), fill = '#323232')

  p <- p + geom_text_repel(data=Names, aes(y=-10*layer+4, x=pos, label=name), force=20, size=4, 
                           family='Helvetica', fontface='bold.italic', color='black', 
                           segment.size=0.3, direction='y', ylim=c(NA, -6) ) 

  ###final touch
  p <- p + theme_minimal() + theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(), axis.ticks=element_blank(), 
    panel.grid.major.x=element_line(linewidth=0.3, color='#cccccc'),
    panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
    axis.line.y=element_line(color='black', linewidth=0.5)
    ) 

  res_list=list(nLayer=max(Tracks$layer), plot=p)
  return(res_list)
}

##plot the (left, right) region on chromosome ch, with scores, ancestry track, gene track, and inferred parameter track
plotPeak <- function(ch, inputname, left, right, outputname, statname, uncorrected_cutoff=0.05){
  DT = data.table(read.table(inputname, header=TRUE, sep="\t", comment="#"))
  # colnames(DT) <- c('ID', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood', 'MLR', 'chi2_p')
  DT[, c("Chr", "physPos", "rsID", "placeholder4", "placeholder5")] <- tstrsplit(DT$ID, "_", fixed = TRUE)
  ### format
  DT$physPos <- as.numeric(DT$physPos)

  ###get peak
  peak=DT[physPos > left & physPos < right]

  ###plot the statistics
  if (statname == "MLR"){
    p1 <- plotManChrs(peak[,.(Chr, physPos, stat=MLR, colorscale=s2hat)], 
                      left=left, right=right, statname="Max Likelihood Ratio",
                      colorstat=TeX(r"($\hat{s}_{AA}$)"))
    # label the top SNP
    top <- peak[which.max(MLR)]
    p1 <- p1 + geom_label_repel(data=top, aes(x=physPos, y=MLR, label=rsID), arrow=arrow(length=unit(4, "pt")),
                                  force=20, direction="both", min.segment.length=unit(8, "pt"), segment.linetype=1, #
                                  fill="white", size=6, box.padding=0.5)
    # need to rotate y axis label
    p1 <- p1 + theme_minimal() + theme(axis.title.y=element_text(angle=90))
  }else if (statname == "pval"){
    p1 <- plotManChrs(peak[,.(Chr, physPos, stat=-log10(chi2_p), colorscale=s2hat)],
                      left=left, right=right, statname=TeX(r"($-log_{10}(p_{\chi^2(1)}$))"),
                      colorstat=TeX(r"($\hat{s}_{AA}$)") )
    # need to rotate y axis label
    p1 <- p1 + theme_minimal() #+ theme(axis.title.y=element_text(angle=0))
    # maybe label the top SNPs?
    tops <- peak[order(MLR, decreasing=TRUE)]
    tops$logP <- -log10(tops$chi2_p)
    tops <- tops[logP > 12]
    p1 <- p1 + geom_label_repel(data=tops, aes(x=physPos, y=logP, label=rsID), arrow=arrow(length=unit(4, "pt")),
                                  force=20, direction="both", min.segment.length=unit(8, "pt"), segment.linetype=1, #
                                  fill="white", size=5, box.padding=0.5, ylim=c(9,NA))
    ### count all sites & draw sig. cutoff
    # if (uncorrected_cutoff != FALSE){
    #   n_sites <- nrow(DT)
    #   cat(n_sites, " SNPs on the chromosome\n")
    #   cutoff <- uncorrected_cutoff / n_sites
    #   p1 <- p1 + geom_hline(yintercept=-log10(cutoff), linewidth=0.5, linetype='dashed')
    # }
  } else if (statname == "MLE"){ # TBD: var to be color scale for s2hat
    p1 <- plotManChrs(peak[,.(Chr, physPos, stat=s2hat, colorscale=-log10(chi2_p))],
                      left=left, right=right, statname=TeX(r"($\hat{s}_{AA}$)"),
                      colorstat=TeX(r"($-log_{10}p$)"))
  }

  ## draw y=0
  p1 <- p1 + geom_hline(yintercept=0, color='darkgray', linewidth=0.5)
  
  ###final touches, make sure to remove the x-axis
  p1 <- p1 + 
    theme(panel.grid.major.x=element_line(linewidth=0.3, color='#cccccc'), 
      legend.position = 'top', legend.direction = 'horizontal', 
      legend.key.width = unit(50, "pt"),
      legend.text=element_text(size=16), legend.title=element_text(size=16), 
      axis.title.y=element_text(size=18, hjust=0.5, vjust=0.5, margin=margin(r=0)), #, angle=0
      axis.text.y=element_text(size=15), axis.line.y=element_line(color='black', linewidth=0.5), #
      axis.title.x=element_blank(), axis.text.x=element_blank(), axis.line.x=element_blank(),  #, element_line(color='black', linewidth=0.5)
      axis.ticks.x=element_blank() 
      )

  ###Modify margins
  p1 <- p1 + theme(plot.margin=margin(t=10, r=25, b=-1, l=35, unit='pt'))
  # pa <- pa + theme(plot.margin=margin(t=-1, r=25, b=0, l=35, unit='pt')) 
  # p2 <- p2 + theme(plot.margin=margin(t=-1, r=25, b=0, l=35, unit='pt'))

  ###start plotting gene track
  genes <- getGeneBed(ch, left, right)
  ###adjust dataframe for plotting
  genes = genes[,.(Chr, start=max(start, left), end=min(end, right), cdStart, cdEnd, name, exonStarts, exonEnds), by=ID]
  ###plotting
  if(nrow(genes)>0){
    print(unique(genes$name))
    print('Plotting genes...')
    pgList <- plotFancyGenes(genes, noRNA=TRUE)
    pg <- pgList$plot
    layers <- pgList$nLayer
    pg <- pg + xlim(left, right) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),#axis.line.x=element_blank(),
      axis.ticks=element_blank(), 
      panel.grid.major.x=element_line(linewidth=0.3, color='#cccccc'),
      panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
      axis.line.y=element_line(color='black', linewidth=0.5), plot.margin=margin(t=0, r=25, b=5, l=15, unit='pt'))

    print(layers)
    if (layers > 10){
      n=3
    }else if (layers>5){
      n = 2.5
    }else{
      n=0.5*layers
    }
    ###remove x axis for score plot
    p1 <- p1 + theme( axis.line.x=element_blank(), 
      axis.title.x=element_blank(), axis.text.x=element_blank(),    
      axis.line.y=element_line(color='black', linewidth=0.5),
      plot.margin=margin(t=-1, r=25, b=0, l=35, unit='pt')
      )
    #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.line.x=element_blank() )
    ### keep x axis for gene plot
    pg <- pg  + xlab(paste('Position on chromosome', ch) ) + 
      theme( axis.title.x=element_text(size=18, hjust=0.5), 
        axis.text.x=element_text(size=15, color='black'),
        axis.line=element_line(color='black', linewidth=0.5), plot.margin=margin(t=-1, b=5, unit='pt')
        )

    ###piece tracks together. the "align='v'" argument makes sure the x-axis of all plots align vertically. Labels can be add onto the plots too.
    pics <- ggarrange(p1, pg, heights=c(4, n+1.5), ncol=1, nrow=2, align='v'
      # labels=c('',''), #hjust=-1, vjust=1.1, , 
      # font.label = list(size = 24, color = "black", face = "bold", family = 'Helvetica')
      )
    ###ggsave() will automatically output the plot in the format indicated in the filename
    ggsave(outputname, plot=pics, width=14, height=7+n*.5, units='in', dpi=500, bg = 'white')
  }else{
    ###No gene in the area, dd x axis back to P1
    p1 <- p1 + theme( axis.title.x=element_text(size=18, hjust=0.5), 
      axis.text.x=element_text(size=15, color='black'),
      axis.line=element_line(color='black', linewidth=0.5), plot.margin=margin(b=5, unit='pt')
      )
    ggsave(outputname, plot=p1, width=14, height=7, units='in', dpi=500, bg = 'white') 
  }
  if(interactive())  return(DT)
}


if(!interactive()) {
    args <-commandArgs(trailingOnly=TRUE)
    print(args)
    inputname=args[1]
    ch=as.integer(args[2])
    left=as.numeric(args[3])
    right=as.numeric(args[4])
    statname=args[5]
    outputname=args[6]
    plotPeak(ch, inputname, left=left, right=right, outputname, statname=statname)
}
