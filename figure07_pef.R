source('R/setup.R',local=TRUE)

groups <- read_excel(clinical_file, sheet="Baseline characteristics") |>
  transmute(id=trimws(`Subject ID`), group=trimws(Group)) |>
  filter(id %in% infected_ids)

stopifnot(!anyDuplicated(groups$id),setequal(groups$id,infected_ids))

spf <- read_excel(pef_file, sheet="Cleaned up Smart Peak Flow", range=cell_cols("A:J"),
                  .name_repair="unique_quiet")

nuvo <- read_excel(pef_file, sheet="Cleaned up Nuvoair", range=cell_cols("A:AK"),
                   .name_repair="unique_quiet")

to_number <- function(x) suppressWarnings(as.numeric(as.character(x)))

spf_data <- tibble(source="Smart Peak Flow", excel_row=seq_len(nrow(spf))+1L,
  id=trimws(as.character(spf[[2]])), timestamp=spf[[5]],
  day_label=as.character(spf[[6]]), period=trimws(as.character(spf[[4]])),
  pef_raw=to_number(spf[[7]]), period_method="Source Day/Night label")

day_start_hour <- 5

night_start_hour <- 17

hour <- as.numeric(format(nuvo[[7]], "%H", tz="UTC")) +
        as.numeric(format(nuvo[[7]], "%M", tz="UTC"))/60

nuvo_data <- tibble(source="NuvoAir", excel_row=seq_len(nrow(nuvo))+1L,
  id=trimws(as.character(nuvo[[2]])), timestamp=nuvo[[7]],
  day_label=as.character(nuvo[[12]]),
  period=ifelse(is.na(hour), NA_character_,
                ifelse(hour>=day_start_hour & hour<night_start_hour, "Day", "Night")),
  pef_raw=to_number(nuvo[[21]]),
  period_method="Inferred: 05:00–16:59 Day; 17:00–04:59 Night; recorded study day retained")

records <- bind_rows(spf_data, nuvo_data) |>
  filter(id %in% infected_ids) |> left_join(groups, by="id") |>
  mutate(day=to_number(day_label),
         valid_pef=is.finite(pef_raw) & pef_raw>0,
         pef=ifelse(valid_pef, pef_raw, NA_real_),
         in_window=!is.na(day) & day %in% 0:13 & period %in% c("Day","Night"))

q <- function(x,p) if(all(is.na(x))) NA_real_ else
  as.numeric(quantile(x,p,na.rm=TRUE,type=7))

devices <- records |> distinct(id,source) |> count(id)

stopifnot(all(devices$n==1))

repeated <- records |> filter(in_window) |> group_by(id,source,day,period) |>
  filter(n()>1) |> ungroup()

daily <- records |> filter(in_window) |> group_by(id,source,day,period) |>
  summarise(n_sessions=sum(!is.na(pef)), pef=q(pef,.5), .groups="drop")

baselines <- daily |> filter(day==0) |> select(id,period,baseline=pef)

followup <- daily |> filter(day>0) |> group_by(id,period) |>
  summarise(n_followup_days=sum(!is.na(pef)), .groups="drop")

eligibility <- expand_grid(id=infected_ids,period=c("Day","Night")) |>
  left_join(groups,by="id") |> left_join(baselines,by=c("id","period")) |>
  left_join(followup,by=c("id","period")) |>
  mutate(n_followup_days=coalesce(n_followup_days,0L), # COUNT only; not a PEF value
         eligible=!is.na(baseline) & baseline>0 & n_followup_days>0,
         reason=case_when(eligible ~ "Included",
            is.na(baseline) & n_followup_days==0 ~ "No usable day-0 baseline or follow-up",
            is.na(baseline) ~ "No usable day-0 baseline for this period",
            TRUE ~ "No usable follow-up for this period"))

dat <- expand_grid(id=infected_ids,period=c("Day","Night"),day=0:13) |>
  left_join(daily,by=c("id","period","day")) |>
  left_join(eligibility,by=c("id","period")) |>
  mutate(change_pct=ifelse(eligible & !is.na(pef),100*(pef-baseline)/baseline,NA_real_),
         group=factor(group,levels=c("Asthma","Control")))

stopifnot(all(is.na(dat$change_pct[is.na(dat$pef)])),
          all(is.na(dat$change_pct[!dat$eligible])),
          all(dat$change_pct[dat$eligible & dat$day==0]==0),
          !anyDuplicated(dat[c("id","period","day")]))

daily_summary <- dat |> group_by(period,group,day) |>
  summarise(n=sum(!is.na(change_pct)), median=q(change_pct,.5),
            q1=q(change_pct,.25),q3=q(change_pct,.75),.groups="drop")

counts <- eligibility |> group_by(period,group) |> summarise(n=sum(eligible),.groups="drop")

colours <- c(Asthma="#FF7068",Control="#00B7BF")

shapes <- c(Asthma=17,Control=16)

font_axis <- 22

font_title <- 25

font_count <- 20

common_theme <- theme_classic(base_size=font_axis) +
  theme(axis.text=element_text(size=font_axis,colour="black"),
        axis.title=element_text(size=font_title,face="bold"),
        axis.title.y=element_text(margin=margin(r=8)),
        axis.title.x=element_text(margin=margin(t=10)),
        legend.title=element_blank(),legend.text=element_text(size=25),
        legend.key.spacing.x=grid::unit(1.3,"cm"),
        plot.tag=element_text(size=27,face="bold",margin=margin(b=12)),
        plot.margin=margin(15,20,5,12))

yr <- range(c(daily_summary$q1,daily_summary$q3,0),na.rm=TRUE)

yr <- yr + c(-1,1)*max(diff(yr)*.1,2)

main_panel <- function(which,tag,label) {
  d <- filter(daily_summary,period==which) |>
    group_by(group) |> arrange(day,.by_group=TRUE) |>
    mutate(segment=cumsum(is.na(median))) |> ungroup()
  ggplot(d,aes(day,median,colour=group,shape=group)) +
    geom_hline(yintercept=0,linetype="dashed",linewidth=.6,colour="grey35") +
    geom_errorbar(aes(ymin=q1,ymax=q3),width=.15,linewidth=.7,na.rm=TRUE,show.legend=FALSE) +
    geom_line(aes(group=interaction(group,segment)),linewidth=1,na.rm=TRUE) +
    geom_point(size=3.5,na.rm=TRUE) +
    scale_colour_manual(values=colours,drop=FALSE) +
    scale_shape_manual(values=shapes,drop=FALSE) +
    scale_x_continuous(breaks=0:13,limits=c(-.45,13.45),expand=c(0,0)) +
    scale_y_continuous(limits=yr,expand=c(0,0)) +
    labs(x="Day post-inoculation",y=label,tag=tag) + common_theme
}

count_panel <- function(which) {
  ggplot(filter(daily_summary,period==which),aes(day,group,label=n,colour=group)) +
    geom_text(size=font_count/(72.27/25.4),show.legend=FALSE) +
    scale_colour_manual(values=colours) +
    scale_x_continuous(limits=c(-.45,13.45),breaks=0:13,expand=c(0,0)) +
    scale_y_discrete(limits=c("Control","Asthma"),
       labels=c(Control="Control (n =)",Asthma="Asthma (n =)")) +
    labs(x=NULL,y=NULL) + theme_classic(base_size=font_count) +
    theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
          axis.text.y=element_text(colour="black",size=font_count),
          plot.margin=margin(0,20,10,12),legend.position="none")
}

p <- (free(main_panel("Day","a","Morning PEF change\nfrom day 0 (%)"),type="label",side="l") +
      free(main_panel("Night","b","Evening PEF change\nfrom day 0 (%)"),type="label",side="l") +
      count_panel("Day") + count_panel("Night") + guide_area()) +
  plot_layout(design="EE\nAB\nCD",heights=c(.65,6.2,1),guides="collect") &
  theme(legend.position="top")

ggsave(file.path(output_dir,"Figure7_PEF.png"),p,device=ragg::agg_png,
       width=20,height=10,units="in",dpi=300,bg="white")

