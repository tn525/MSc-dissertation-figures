source('R/setup.R',local=TRUE)

baseline <- read_excel(input_file,sheet='Baseline characteristics') |>
  transmute(id=trimws(`Subject ID`),group=trimws(Group))
qpcr <- read_excel(input_file,sheet='Viral copies',range=cell_cols('A:F')) |>
  transmute(id=trimws(ID)) |> distinct() |> filter(!is.na(id),nzchar(id)) |>
  left_join(baseline,by='id') |> mutate(infected=id %in% infected_ids)
stopifnot(!anyNA(qpcr$group))
counts <- qpcr |> group_by(group) |>
  summarise(total=n(),infected=sum(infected),.groups='drop') |>
  mutate(group=factor(group,levels=c('Asthma','Control')),
         percent=100*infected/total,
         label=sprintf('%d/%d\n(%.1f%%)',infected,total,percent))
p_value <- fisher.test(with(counts,cbind(infected,total-infected)))$p.value
p <- ggplot(counts,aes(group,percent,fill=group))+
  geom_col(width=.58)+geom_text(aes(label=label),vjust=-.5,size=6)+
  annotate('text',x=1.5,y=94,label=sprintf("Fisher's exact P = %.3f",p_value),size=6)+
  scale_fill_manual(values=c(Asthma='#FF7068',Control='#00B7BF'),guide='none')+
  scale_x_discrete(labels=c(Asthma='Asthma',Control='Healthy'))+
  scale_y_continuous(breaks=seq(0,100,20),labels=function(x) paste0(x,'%'),
                     limits=c(0,100),expand=c(0,0))+
  labs(x=NULL,y='Successfully infected (%)')+
  theme_classic(base_size=22)+theme(axis.text=element_text(colour='black'),
    axis.text.x=element_text(face='bold'),axis.title.y=element_text(face='bold'))
ggsave(file.path(output_dir,'Figure4_infection_proportion.png'),p,
       device=ragg::agg_png,width=10,height=8,dpi=300,bg='white')
