source('R/prepare_assays.R',local=TRUE)

library(patchwork)

analytes <- c('lambda1','alpha2a','gamma','ip10')

labels <- o$labels;times <- c(0,3,6,24,48,72,96)

labels <- o$labels;times <- c(0,3,6,24,48,72,96)

dat <- o$msd |> filter(analyte %in% analytes) |> mutate(group=factor(group,levels=c('Asthma','Control')),visit_index=match(time,times),value=msd_analysis,status=msd_status)

summ <- o$msd_summary |> filter(analyte %in% analytes) |> mutate(visit_index=match(time,times),n=n_analysis)

active <- dat |> group_by(id,group) |> summarise(active=any(status!="Missing"),.groups='drop') |>
 filter(active) |> arrange(group,id) |> group_by(group) |>
 mutate(offset=seq(-.105,.105,length.out=n())+ifelse(group=='Asthma',-.12,.12)) |> ungroup()

dat <- dat |> left_join(select(active,id,offset),by='id') |> mutate(x=visit_index+offset)

cols <- c(Asthma='#FF7068',Control='#00B7BF',AsthmaBelow='#FF7068',ControlBelow='#00B7BF',AsthmaND='#FF7068',ControlND='#00B7BF',Above='black')

shapes <- c(Asthma=17,Control=16,AsthmaBelow=2,ControlBelow=1,AsthmaND=4,ControlND=4,Above=8)

legend_labels <- c(Asthma='Asthma: quantifiable',Control='Control: quantifiable',AsthmaBelow='Asthma: detected <LLOQ',ControlBelow='Control: detected <LLOQ',AsthmaND='Asthma: ND',ControlND='Control: ND',Above='Above fit range: source flag')

legend_order <- names(cols)

base <- theme_classic(base_size=23)+theme(axis.text=element_text(size=21,colour='black'),
 axis.title=element_text(size=24,face='bold'),axis.title.y=element_text(margin=margin(r=8)),
 legend.text=element_text(size=18),legend.key.spacing.x=grid::unit(.8,'cm'),legend.key.spacing.y=grid::unit(.2,'cm'),
 plot.tag=element_text(size=27,face='bold'),plot.margin=margin(16,20,8,12))

main <- function(a,tag) {
 d <- filter(dat,analyte==a,!is.na(offset));s <- filter(summ,analyte==a)
 lo <- floor(log10(min(d$value[d$value>0],na.rm=TRUE)));hi <- ceiling(log10(max(d$value,na.rm=TRUE)))
 below <- 10^(lo-.75)
 d <- d |> mutate(y=ifelse(status=='ND',below,value)) |> group_by(id) |> arrange(visit_index,.by_group=TRUE) |>
 mutate(xnext=lead(x),ynext=lead(y),flag_link=status %in% c('ND','Detected <LLOQ') |
 lead(status) %in% c('ND','Detected <LLOQ')) |> ungroup()
 s <- s |> mutate(x=visit_index+ifelse(group=='Asthma',-.12,.12)) |> group_by(group) |>
 arrange(visit_index,.by_group=TRUE) |> mutate(segment=cumsum(is.na(median))) |> ungroup()
 has_below <- any(d$status=='ND')
 breaks <- 10^(lo:hi); labs_y <- vapply(breaks,function(z) format(z,scientific=FALSE,trim=TRUE,big.mark=','),character(1))
 if(has_below) {breaks <- c(below,breaks);labs_y <- c('ND',labs_y)}
 ggplot()+
 {if(has_below) geom_hline(yintercept=below,colour='grey85',linewidth=.4) else NULL}+
 geom_segment(data=filter(d,!flag_link),aes(x=x,xend=xnext,y=y,yend=ynext,colour=group),alpha=.3,linewidth=.6,na.rm=TRUE,show.legend=FALSE)+
 geom_segment(data=filter(d,flag_link),aes(x=x,xend=xnext,y=y,yend=ynext,colour=group),alpha=.35,linewidth=.6,linetype='dashed',na.rm=TRUE,show.legend=FALSE)+
 geom_point(data=filter(d,status=='Quantifiable'),aes(x,y,colour=group,shape=group),size=2.4,alpha=.5,show.legend=FALSE)+
 geom_point(data=filter(d,status=='Detected <LLOQ'),aes(x,y,colour=paste0(group,'Below'),shape=paste0(group,'Below')),size=3,stroke=.9,show.legend=TRUE)+
 geom_point(data=filter(d,status=='ND'),aes(x,y,colour=paste0(group,'ND'),shape=paste0(group,'ND')),size=3,stroke=.9,show.legend=TRUE)+
 {if(unique(d$lloq)>=if(has_below) below/1.5 else 10^lo/1.5) geom_hline(yintercept=unique(d$lloq),linetype='dotted',colour='grey55',linewidth=.5) else NULL}+
 geom_point(data=filter(d,upper_range_flag=='Above Fit Curve Range'),aes(x,y,colour='Above',shape='Above'),size=3.3,stroke=.9,show.legend=TRUE)+
 geom_errorbar(data=s,aes(x,ymin=q1,ymax=q3,colour=group),linewidth=.9,width=.07,na.rm=TRUE,show.legend=FALSE)+
 geom_line(data=s,aes(x,median,colour=group,group=interaction(group,segment)),linewidth=1.7,na.rm=TRUE,show.legend=FALSE)+
 geom_point(data=s,aes(x,median,colour=group,shape=group),size=6,na.rm=TRUE,show.legend=FALSE)+
 scale_fill_manual(values=cols,guide='none')+
 scale_colour_manual(values=cols,limits=legend_order,breaks=legend_order,labels=legend_labels,name=NULL,drop=FALSE)+
 scale_shape_manual(values=shapes,limits=legend_order,breaks=legend_order,labels=legend_labels,name=NULL,drop=FALSE)+
 guides(colour=guide_legend(ncol=2,byrow=TRUE,override.aes=list(alpha=1,size=3.2)),shape=guide_legend(ncol=2,byrow=TRUE))+
 scale_x_continuous(breaks=1:7,labels=times,limits=c(.6,7.4),expand=c(0,0))+
 scale_y_log10(breaks=breaks,labels=labs_y,limits=c(if(has_below) below/1.5 else 10^lo/1.5,10^hi*1.5),expand=c(0,0))+
 labs(x='Time post-inoculation (hours)',y=paste0(labels[[a]],'\n(pg/mL; log scale)'),tag=tag)+base
}

counts <- function(a) ggplot(filter(summ,analyte==a),aes(visit_index,group,label=n,colour=group))+
 geom_text(size=7,show.legend=FALSE)+scale_colour_manual(values=cols)+
 scale_x_continuous(breaks=1:7,labels=times,limits=c(.6,7.4),expand=c(0,0))+
 scale_y_discrete(limits=c('Control','Asthma'),labels=c(Asthma='Asthma (n =)',Control='Control (n =)'))+
 labs(x=NULL,y=NULL)+theme_classic(base_size=21)+theme(axis.text.x=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(colour='black'),plot.margin=margin(0,20,8,12),legend.position='none')

p <- (main('lambda1','a')+main('alpha2a','b')+counts('lambda1')+counts('alpha2a')+
 main('gamma','c')+main('ip10','d')+counts('gamma')+counts('ip10')+guide_area())+
 plot_layout(design='II\nAB\nCD\n##\nEF\nGH',heights=c(1.7,5, .8,.5,5,.8),guides='collect') & theme(legend.position='top')

ggsave(file.path(output_dir,'Figure10_96h_MSD_kinetics.png'),p,device=ragg::agg_png,width=20,height=17.5,dpi=300,bg='white')

