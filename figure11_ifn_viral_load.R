source('R/prepare_assays.R',local=TRUE)

audit <- o$change_pairs |> filter(analyte=='lambda1') |> rename(ifn_change=change) |>
 left_join(o$peaks,by=c('id','group')) |> mutate(included=is.finite(ifn_change)&is.finite(peak)&peak>0)

dat <- audit |> filter(included) |> arrange(id) |> mutate(group=factor(group,levels=c('Asthma','Control')))

gn <- table(dat$group)

glabels <- c(Asthma=paste0('Asthma (n = ',gn['Asthma'],')'),Control=paste0('Control (n = ',gn['Control'],')'))

B <- 99999L

spearman_permutation <- function(d,seed) {
 x<-rank(d$ifn_change,ties.method='average');y<-rank(d$peak,ties.method='average')
 x<-x-mean(x);y<-y-mean(y);den<-sqrt(sum(x*x)*sum(y*y))
 rho<-sum(x*y)/den
 set.seed(seed)
 perm<-replicate(B,sum(x*sample(y))/den)
 extreme<-sum(abs(perm)>=abs(rho)-1e-12)
 pv<-(extreme+1)/(B+1)
 tibble(n=nrow(d),n_asthma=sum(d$group=='Asthma'),n_control=sum(d$group=='Control'),
 rho=rho,p=pv,n_permutations=B,mc_se=sqrt(pv*(1-pv)/(B+1)),seed=seed)
}

main <- spearman_permutation(dat,240926)

cols<-c(Asthma='#FF7068',Control='#00B7BF')

fmtp<-function(p) if(p<.001) '< 0.001' else paste0('= ',sprintf('%.3f',p))

annotation<-sprintf('Exploratory: n = %d\nSpearman ρ = %.2f; permutation p %s',main$n,main$rho,fmtp(main$p))

trend_model <- lm(log10(peak) ~ ifn_change,data=dat)

trend <- data.frame(ifn_change=seq(min(dat$ifn_change),max(dat$ifn_change),length.out=200))

trend$peak <- 10^predict(trend_model,newdata=trend)

p<-ggplot(dat,aes(ifn_change,peak,colour=group,shape=group))+
 geom_vline(xintercept=0,linetype='dashed',colour='grey75',linewidth=.5)+
 geom_line(data=trend,aes(ifn_change,peak),inherit.aes=FALSE,
           colour='grey30',linewidth=1.1,show.legend=FALSE)+
 geom_point(size=4.8)+
 scale_colour_manual(values=cols,labels=glabels,name=NULL)+
 scale_shape_manual(values=c(Asthma=17,Control=16),labels=glabels,name=NULL)+
 scale_y_log10(breaks=c(1,10,100,1000,10000),labels=c('1','10','100','1,000','10,000'),
 limits=c(1,100000),expand=c(0,0))+
 scale_x_continuous(expand=expansion(mult=c(.1,.1)))+
 annotate('text',x=mean(range(dat$ifn_change)),y=60000,label=annotation,size=6.3,lineheight=1.15)+
 labs(x='Early IFN-λ1 change: 6 h − baseline (pg/mL)',
 y='Subsequent peak nasal viral RNA, 24–96 h\n(copies/µL; log scale)')+
 theme_classic(base_size=23)+theme(axis.text=element_text(size=22,colour='black'),
 axis.title=element_text(size=24,face='bold'),axis.title.y=element_text(margin=margin(r=10)),
 axis.title.x=element_text(margin=margin(t=12)),legend.position='top',legend.text=element_text(size=24),
 legend.key.spacing.x=grid::unit(1.3,'cm'),plot.margin=margin(20,24,15,15))

ggsave(file.path(output_dir,'Figure11_IFN_lambda1_vs_peak_viral_load.png'),p,
 device=ragg::agg_png,width=12,height=9,dpi=300,bg='white')

