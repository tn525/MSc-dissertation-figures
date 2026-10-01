source('R/setup.R',local=TRUE)

baseline <- read_excel(input_file,sheet="Baseline characteristics") |>
 transmute(id=trimws(`Subject ID`),group=trimws(Group)) |> filter(id %in% infected_ids)

stopifnot(!anyDuplicated(baseline$id),setequal(baseline$id,infected_ids))

raw <- read_excel(input_file,sheet="Cleaned up symptom score") |>
 filter(trimws(`Subject ID`) %in% infected_ids)

stopifnot(!anyDuplicated(raw[c("Subject ID","Day")]),
          all(raw$Day %in% 0:13))

matched_group <- baseline$group[match(trimws(raw$`Subject ID`),baseline$id)]

stopifnot(all(trimws(raw$Group)==matched_group))

parse_item <- function(x) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | x==""
  valid <- grepl("\\([0-3]\\)$",x)
  if(any(!missing & !valid)) stop("Unrecognised symptom response; check source data.")
  result <- rep(NA_real_,length(x))
  result[!missing] <- as.numeric(sub(".*\\(([0-3])\\)$","\\1",x[!missing]))
  result
}

items <- as.data.frame(lapply(raw[6:20],parse_item))

upper <- rowSums(items[1:8],na.rm=FALSE)

lower <- rowSums(items[9:15],na.rm=FALSE)

check_cached <- function(recomputed,cached) {
  ok <- !is.na(recomputed) & !is.na(cached)
  stopifnot(all(recomputed[ok]==cached[ok]))
}

check_cached(upper,raw$`Upper respiratory score`)

check_cached(lower,raw$`Lower respiratory score`)

scores <- data.frame(id=trimws(raw$`Subject ID`),day=raw$Day,upper=upper,lower=lower)

dat <- expand_grid(id=infected_ids,day=0:13) |>
 left_join(scores,by=c("id","day")) |> left_join(baseline,by="id") |>
 pivot_longer(c(upper,lower),names_to="outcome",values_to="score") |>
 mutate(group=factor(group,levels=c("Asthma","Control")))

q <- function(x,p) if(all(is.na(x))) NA_real_ else
 as.numeric(quantile(x,p,na.rm=TRUE,type=7))

summary <- dat |> group_by(outcome,group,day) |>
 summarise(n=sum(!is.na(score)),median=q(score,.5),q1=q(score,.25),q3=q(score,.75),
           .groups="drop")



colours <- c(Asthma="#FF7068",Control="#00B7BF")

shapes <- c(Asthma=17,Control=16)

font_title <- 27

font_label <- 25

font_value <- 22

common_theme <- theme_classic(base_size=font_label) +
 theme(axis.text=element_text(size=font_value,colour="black"),
       axis.title=element_text(size=font_title,face="bold"),
       axis.title.y=element_text(size=22,face="bold",margin=margin(r=8)),
       legend.text=element_text(size=font_label),legend.title=element_blank(),
       legend.key.spacing.x=grid::unit(1.4,"cm"),
       plot.tag=element_text(size=24,face="bold",margin=margin(b=14)),
       plot.margin=margin(22,16,4,12))

main_panel <- function(which,ylab,tag) {
 d <- filter(summary,outcome==which)
 # Both groups use the exact day coordinate, matching the original figure.
 d <- mutate(d,display_day=day)
 ggplot(d,aes(display_day,median,colour=group,group=group,shape=group)) +
  geom_errorbar(aes(ymin=q1,ymax=q3),width=.12,linewidth=.65,na.rm=TRUE,
                show.legend=FALSE) +
  geom_line(linewidth=1,na.rm=TRUE) + geom_point(size=3.3,na.rm=TRUE) +
  scale_colour_manual(values=colours) + scale_shape_manual(values=shapes) +
  scale_x_continuous(breaks=0:13,limits=c(-.45,13.45),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.025,.07))) +
  labs(x="Day post-inoculation",y=ylab,tag=tag) + common_theme
}

count_panel <- function(which) {
 d <- filter(summary,outcome==which)
 ggplot(d,aes(day,group,label=n,colour=group)) +
  geom_text(size=font_value / (72.27 / 25.4),show.legend=FALSE) +
  scale_colour_manual(values=colours) +
  scale_x_continuous(breaks=0:13,limits=c(-.45,13.45),expand=c(0,0)) +
  scale_y_discrete(limits=c("Control","Asthma"),
                   labels=c(Control="Control (n =)",Asthma="Asthma (n =)")) +
  labs(x=NULL,y=NULL) + theme_classic(base_size=font_label) +
  theme(axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.line=element_blank(),axis.text.y=element_text(colour="black",size=20),
        plot.margin=margin(0,16,8,12),legend.position="none")
}

auc_data <- dat |> group_by(id,group,outcome) |>
 summarise(n_days=sum(!is.na(score)),
   auc={
     ord <- order(day); t <- day[ord]; y <- score[ord]
     if(anyNA(y) || !identical(as.integer(t),0:13)) NA_real_ else
       sum(diff(t)*(head(y,-1)+tail(y,-1))/2)
   },.groups="drop")

lower_auc <- auc_data |> filter(outcome=="lower",!is.na(auc))



auc_summary <- lower_auc |> group_by(group) |>
 summarise(n=n(),median=median(auc),q1=as.numeric(quantile(auc,.25)),
           q3=as.numeric(quantile(auc,.75)),.groups="drop")

rank_permutation <- function(value,group) {
 ranks <- rank(value,ties.method="average")
 n_a <- sum(group=="Asthma")
 stopifnot(!anyNA(value),n_a>0,n_a<length(value))
 expected <- n_a*mean(ranks)
 observed <- abs(sum(ranks[group=="Asthma"])-expected)
 permuted <- combn(ranks,n_a,FUN=sum)
 mean(abs(permuted-expected)>=observed-1e-10)
}

tests <- auc_data |> filter(!is.na(auc)) |> group_by(outcome) |>
 summarise(p_raw=rank_permutation(auc,group),.groups="drop") |>
 mutate(p_holm=p.adjust(p_raw,method="holm"),
        method="Exact two-sided permutation rank-sum; average ranks for ties")

auc_panel <- function(which,tag) {
 plot_data <- filter(auc_data,outcome==which,!is.na(auc))
 auc_summary <- plot_data |> group_by(group) |>
  summarise(n=n(),median=median(auc),q1=as.numeric(quantile(auc,.25)),
            q3=as.numeric(quantile(auc,.75)),.groups="drop")
 test_lower <- filter(tests,outcome==which)
label_p <- sprintf("Exact permutation p = %.3f\nHolm-adjusted p = %.3f",test_lower$p_raw,test_lower$p_holm)
top_y <- max(plot_data$auc)*1.15
ggplot(plot_data,aes(group,auc)) +
 geom_point(aes(colour=group,shape=group),
   position=position_jitter(width=.12,height=0,seed=42),size=4,alpha=.85,
   show.legend=FALSE) +
 geom_errorbar(data=auc_summary,aes(x=group,ymin=q1,ymax=q3),
   inherit.aes=FALSE,width=.16,linewidth=1,colour="black") +
 geom_crossbar(data=auc_summary,aes(x=group,y=median,ymin=median,ymax=median),
   inherit.aes=FALSE,width=.36,linewidth=.8,colour="black") +
 annotate("segment",x=1,xend=2,y=top_y,yend=top_y,linewidth=.5) +
 annotate("text",x=1.5,y=top_y*1.045,label=label_p,vjust=0,size=6.5,lineheight=1.1) +
 scale_colour_manual(values=colours) + scale_shape_manual(values=shapes) +
 scale_x_discrete(labels=setNames(paste0(auc_summary$group,"\n(n = ",auc_summary$n,")"),auc_summary$group)) +
 scale_y_continuous(limits=c(0,top_y*1.28),expand=expansion(mult=c(.025,0))) +
 labs(x=NULL,y=paste0(if(which=="upper") "Upper" else "Lower",
      " respiratory symptom AUC\n(score-days, days 0–13)"),tag=tag) +
 common_theme + theme(legend.position="none",axis.title.y=element_text(size=22,margin=margin(r=8)))

}

p <- (free(main_panel("upper","Upper respiratory symptom score","a"),type="label",side="l") +
      free(main_panel("lower","Lower respiratory symptom score","b"),type="label",side="l") +
      count_panel("upper") + count_panel("lower") +
      free(auc_panel("upper","c"),type="label",side="l") +
      free(auc_panel("lower","d"),type="label",side="l") +
      guide_area()) +
 plot_layout(design="GG\nAB\nCD\n##\nEF",heights=c(.65,6,1,.75,5.5),guides="collect") &
 theme(legend.position="top")

ggsave(file.path(output_dir,"Figure6_symptoms_and_AUC.png"),p,
 device=ragg::agg_png,width=20,height=15.8,units="in",dpi=300,bg="white")

