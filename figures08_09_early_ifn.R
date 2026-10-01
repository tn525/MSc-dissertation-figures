source('R/prepare_assays.R',local=TRUE)

library(patchwork)

dat <- o$msd |> filter(time<=6,analyte %in% c('lambda1','alpha2a','beta','gamma')) |>
 mutate(group=factor(group,levels=c('Asthma','Control')),value=msd_analysis)

summary <- o$msd_summary |> filter(time<=6)

tests <- o$change_tests

primary_test <- filter(tests,analyte=='lambda1')

labels <- o$labels

available <- dat |> filter(msd_status!='Missing') |> distinct(id,group) |> count(group)

legend_n <- setNames(available$n,available$group)

colours <- c(Asthma="#FF7068",Control="#00B7BF",AsthmaND="#FF7068",ControlND="#00B7BF",AsthmaLow="#FF7068",ControlLow="#00B7BF")

shapes <- c(Asthma=17,Control=16,AsthmaND=4,ControlND=4,AsthmaLow=2,ControlLow=1)

legend_labels <- c(Asthma=paste0('Asthma: quantifiable (n = ',legend_n['Asthma'],')'),Control=paste0('Control: quantifiable (n = ',legend_n['Control'],')'),
 AsthmaLow='Asthma: detected <LLOQ',ControlLow='Control: detected <LLOQ',AsthmaND='Asthma: ND',ControlND='Control: ND')

legend_order <- names(legend_labels)

format_p <- function(p) if(p<.001) "< 0.001" else paste0("= ",sprintf("%.3f",p))

panel <- function(which,tag=NULL) {
 d <- dat |> filter(analyte==which) |> group_by(id) |> arrange(time,.by_group=TRUE) |>
  mutate(segment=cumsum(is.na(value))) |> ungroup()
 s <- filter(summary,analyte==which)
 tt <- filter(tests,analyte==which)
 test_label <- if(which=="lambda1") paste0("p ",format_p(primary_test$p)) else
  paste0("p ",format_p(tt$p),"; Holm-adjusted p ",format_p(tt$p_holm),"\nChange analysis: asthma n = ",tt$n_asthma,", control n = ",tt$n_control)
 top <- max(d$value,s$q3,na.rm=TRUE)
 nd_y <- -top*.11
 d <- d |> mutate(plot_value=ifelse(msd_status=="ND",nd_y,value)) |>
  group_by(id) |> arrange(time,.by_group=TRUE) |>
  mutate(next_time=lead(time),next_value=lead(plot_value),
         nd_link=msd_status=="ND" | lead(msd_status)=="ND") |> ungroup()
 ticks <- pretty(c(0,top)); ticks <- ticks[ticks>=0 & ticks<=top*1.03]
 ggplot() +
  geom_hline(yintercept=nd_y,colour="grey85",linewidth=.4) +
  geom_segment(data=filter(d,!nd_link),aes(x=time,xend=next_time,y=plot_value,yend=next_value,colour=group),
            linewidth=.65,alpha=.3,na.rm=TRUE,show.legend=FALSE) +
  geom_segment(data=filter(d,nd_link),aes(x=time,xend=next_time,y=plot_value,yend=next_value,colour=group),
            linewidth=.65,alpha=.3,linetype="dashed",na.rm=TRUE,show.legend=FALSE) +
  geom_point(data=filter(d,msd_status=="Quantifiable"),aes(time,value,colour=group,shape=group),
             size=2.4,alpha=.4,show.legend=FALSE) +
  geom_point(data=filter(d,msd_status=="ND"),
             aes(time,plot_value,colour=paste0(group,"ND"),shape=paste0(group,"ND")),
             size=3.4,stroke=1,show.legend=TRUE) +
  geom_point(data=filter(d,msd_status=='Detected <LLOQ'),aes(time,plot_value,colour=paste0(group,'Low'),shape=paste0(group,'Low')),size=3.4,stroke=1,show.legend=TRUE) +
  geom_hline(yintercept=unique(d$lloq),linetype='dotted',colour='grey55',linewidth=.5) +
  geom_errorbar(data=filter(s,!is.na(q1)),aes(time,ymin=q1,ymax=q3,colour=group),width=.13,
                linewidth=1,show.legend=FALSE) +
  geom_segment(data=filter(s,is.na(q1)),aes(x=time,xend=time,y=median,yend=q3,colour=group),
               linewidth=1,show.legend=FALSE) +
  geom_segment(data=filter(s,is.na(q1)),aes(x=time-.065,xend=time+.065,y=q3,yend=q3,colour=group),
               linewidth=1,show.legend=FALSE) +
  geom_line(data=s,aes(time,median,colour=group,group=group),linewidth=1.7,show.legend=FALSE) +
  geom_point(data=s,aes(time,median,colour=group,shape=group),size=6,show.legend=FALSE) +
  scale_fill_manual(values=colours,guide='none') +
  annotate("text",x=3,y=top*1.19,label=test_label,size=6.4,lineheight=1.2) +
  scale_colour_manual(values=colours,breaks=legend_order,limits=legend_order,
                       labels=legend_labels,name=NULL,drop=FALSE) +
  scale_shape_manual(values=shapes,breaks=legend_order,limits=legend_order,
                      labels=legend_labels,name=NULL,drop=FALSE) +
  guides(colour=guide_legend(nrow=3,byrow=TRUE,override.aes=list(size=3.2,alpha=1)),
         shape=guide_legend(nrow=3,byrow=TRUE)) +
  scale_x_continuous(breaks=c(0,3,6),limits=c(-.35,6.35),expand=c(0,0)) +
  scale_y_continuous(limits=c(-top*.17,top*1.32),expand=c(0,0),
                     breaks=c(nd_y,ticks),labels=c("ND",as.character(ticks))) +
  labs(x="Time post-inoculation (hours)",y=paste0(labels[[which]]," (pg/mL)"),tag=tag) +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=22,colour="black"),
        axis.title.x=element_text(size=23,face="bold",margin=margin(t=12)),
        axis.title.y=element_text(size=24,face="bold",margin=margin(r=9)),
        legend.text=element_text(size=16),legend.position="top",
        legend.key.spacing.x=grid::unit(.8,"cm"),
        legend.key.spacing.y=grid::unit(.2,"cm"),
        plot.tag=element_text(size=26,face="bold"),plot.margin=margin(20,20,15,12))
}

p1 <- panel("lambda1")

p2 <- (panel("alpha2a","a") + panel("beta","b") + panel("gamma","c") + guide_area()) +
 plot_layout(design="DDDD\nAABB\n####\n#CC#",heights=c(1.5,6,.35,6),guides="collect") & theme(legend.position="top")

ggsave(file.path(output_dir,"Figure8_Early_IFN_lambda1.png"),p1,
       device=ragg::agg_png,width=10,height=9,units="in",dpi=300,bg="white")

ggsave(file.path(output_dir,"Figure9_Early_IFN_alpha_beta_gamma.png"),p2,
       device=ragg::agg_png,width=20,height=17.5,units="in",dpi=300,bg="white")

