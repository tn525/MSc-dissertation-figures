source('R/prepare_assays.R',local=TRUE)

dat <- o$viral |> mutate(group=factor(group,levels=c('Asthma','Control')),x=time,
 display_x=time+(match(id,sort(unique(id)))-(length(unique(id))+1)/2)/max(1,(length(unique(id))-1)/2)*.36,
 point_key=case_when(viral_status=='Quantifiable'~as.character(group),TRUE~viral_status))

summary <- o$viral_summary |> mutate(x=time,group=factor(group,levels=c('Asthma','Control'))) |>
 group_by(group) |> arrange(time,.by_group=TRUE) |> mutate(segment=cumsum(is.na(median))) |> ungroup()

ng <- o$groups |> count(group)

cols<-c(Asthma='#FF7068',Control='#00B7BF')

keys<-c('Asthma','Control','Detected <LLOQ','Below LLOD','ND')

shapes<-c(Asthma=17,Control=16,'Detected <LLOQ'=5,'Below LLOD'=1,ND=4)

key_labels<-c(paste0('Asthma: quantifiable (n = ',ng$n[ng$group=='Asthma'],')'),paste0('Healthy: quantifiable (n = ',ng$n[ng$group=='Control'],')'),'Detected <LLOQ: 50 copies/µL','Below LLOD','ND')

dat <- dat |> group_by(id) |> mutate(trajectory_nd=any(viral_status=='ND')) |> ungroup()
p <- ggplot(dat,aes(x=display_x,y=viral_plot))+
 geom_line(aes(group=id,colour=group,linetype=trajectory_nd),alpha=.25,linewidth=.45,na.rm=TRUE,show.legend=FALSE)+
 geom_point(aes(shape=point_key,colour=group),size=2.1,alpha=.55,na.rm=TRUE,show.legend=c(colour=FALSE,shape=TRUE))+
 geom_errorbar(data=summary,aes(x=x,ymin=q1,ymax=q3,colour=group),inherit.aes=FALSE,width=1.15,linewidth=.8,na.rm=TRUE,show.legend=FALSE)+
 geom_line(data=summary,aes(x=x,y=median,colour=group,group=interaction(group,segment)),inherit.aes=FALSE,linewidth=1.6,na.rm=TRUE,show.legend=FALSE)+
 geom_point(data=summary,aes(x=x,y=median,colour=group,shape=group),inherit.aes=FALSE,size=5.8,na.rm=TRUE,show.legend=FALSE)+
 scale_fill_manual(values=cols,guide='none')+
 scale_linetype_manual(values=c('FALSE'='solid','TRUE'='dashed'),guide='none')+
 geom_hline(yintercept=c(10,100),linetype=c('dashed','dotted'),colour='grey50',linewidth=.5)+
 annotate('text',x=96,y=c(10,100),label=c('LLOD = 10','LLOQ = 100'),hjust=1,vjust=-.4,size=4)+
 scale_colour_manual(values=cols,name=NULL)+
 scale_shape_manual(values=shapes,breaks=keys,limits=keys,labels=key_labels,name=NULL)+
 scale_y_log10(breaks=c(2,5,10,50,100,1000,10000,100000),labels=c('ND','<LLOD','10','50','100','1,000','10,000','100,000'))+
 scale_x_continuous(breaks=c(0,3,6,24,48,72,96),guide=guide_axis(n.dodge=1))+
 guides(colour='none',shape=guide_legend(nrow=3,byrow=TRUE,override.aes=list(colour=c(cols,'#444444','#444444','#444444'),size=3,alpha=1)))+
 labs(x='Time post-inoculation (hours)',y='Reported viral RNA (copies/µL; logarithmic scale)')+
 theme_classic(base_size=22) +
 theme(strip.background=element_blank(),strip.text=element_text(face="bold",size=15),
       legend.position="top", legend.box="vertical",
       legend.key.spacing.x=grid::unit(.8, "cm"),
       legend.text=element_text(size=16, margin=margin(l=5)),
       axis.text=element_text(size=19,colour="black"),
       axis.title=element_text(size=22),
       plot.caption=element_text(hjust=0,size=10),
       plot.title=element_text(face="bold",size=16))

ggsave(file.path(output_dir,'Figure5_RV_A16_viral_load.png'),p,device=ragg::agg_png,width=14,height=8.5,units='in',dpi=300,bg='white')

