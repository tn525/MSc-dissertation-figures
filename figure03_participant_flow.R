source('R/prepare_assays.R',local=TRUE)

# Recruitment totals transcribed from the draft, Figure 3; not in the clinical workbook.
recruitment <- list(screened=c(Asthma=58,Control=42),inoculated=c(Asthma=21,Control=15),
                   exclusions=c('Positive RV-A16 antibodies'=40,'Withdrew'=11,
                                'Lost to follow-up'=8,'Unable to schedule'=4,'Medical reasons'=1))
stopifnot(sum(recruitment$screened)-sum(recruitment$exclusions)==sum(recruitment$inoculated))
all_groups <- read_excel(input_file,sheet='Baseline characteristics') |>
  transmute(id=trimws(`Subject ID`),group=trimws(Group))
qpcr <- read_excel(input_file,sheet='Viral copies',range=cell_cols('A:F')) |>
  transmute(id=trimws(ID)) |> distinct() |> filter(!is.na(id),nzchar(id)) |> left_join(all_groups,by='id')
symptoms <- read_excel(input_file,sheet='Cleaned up symptom score') |>
  transmute(id=trimws(`Subject ID`),upper=`Upper respiratory score`,lower=`Lower respiratory score`) |>
  filter(id %in% infected_ids,!is.na(upper)|!is.na(lower)) |> distinct(id) |> left_join(groups,by='id')
proteins <- msd |> filter(!is.na(msd_raw)) |> distinct(id,group)
# Reuse the PEF eligibility rules rather than introducing another definition.
pef_env <- new.env(parent=globalenv())
source('R/figure07_pef.R',local=pef_env)
pef_available <- pef_env$eligibility |> filter(eligible) |> distinct(id,group)
label_counts <- function(d,title) {
 n <- table(factor(d$group,levels=c('Asthma','Control')))
 sprintf('%s (n=%d)\nAsthma (n=%d)\nControl (n=%d)',title,sum(n),n[1],n[2])
}
label_totals <- function(n,title) sprintf('%s (n=%d)\nAsthma (n=%d)\nControl (n=%d)',title,sum(n),n[1],n[2])
box <- function(x,y,w,h,text,fill='white',bold=FALSE) {
 grid::grid.rect(x,y,w,h,gp=grid::gpar(fill=fill,col='grey35'))
 grid::grid.text(text,x,y,gp=grid::gpar(fontfamily='Arial',fontsize=12,
                                     fontface=if(bold) 'bold' else 'plain'))
}
arrow <- function(x0,y0,x1,y1) grid::grid.lines(x=c(x0,x1),y=c(y0,y1),
  arrow=grid::arrow(length=grid::unit(.07,'inches'),type='closed'),gp=grid::gpar(col='grey35'))
ragg::agg_png(file.path(output_dir,'Figure3_participant_flow.png'),width=13,height=11,units='in',res=300)
grid::grid.newpage()
box(.43,.93,.29,.09,label_totals(recruitment$screened,'Participants screened'),bold=TRUE)
box(.80,.82,.34,.15,paste0('Excluded before inoculation (n=',sum(recruitment$exclusions),')\n',
 paste0(names(recruitment$exclusions),' (n=',recruitment$exclusions,')',collapse='\n')),fill='#CDE3F2')
arrow(.43,.885,.43,.755);arrow(.43,.83,.63,.83)
box(.43,.71,.29,.09,label_totals(recruitment$inoculated,'Inoculated with RV-A16'),bold=TRUE)
arrow(.43,.665,.43,.59);arrow(.43,.635,.66,.635)
box(.79,.635,.26,.065,sprintf('Viral qPCR results pending\n(n=%d)',sum(recruitment$inoculated)-nrow(qpcr)),fill='#CDE3F2')
box(.43,.54,.29,.10,label_counts(qpcr,'Viral-load qPCR available'),bold=TRUE)
arrow(.43,.49,.43,.425);arrow(.43,.465,.66,.465)
box(.79,.465,.26,.065,sprintf('No confirmed infection\n(n=%d)',sum(!qpcr$id %in% infected_ids)),fill='#CDE3F2')
box(.43,.38,.29,.09,label_counts(groups,'Confirmed infection'),bold=TRUE)
arrow(.43,.335,.43,.28)
grid::grid.lines(x=c(.18,.76),y=c(.28,.28))
for(i in seq_along(c(.18,.46,.76))) {
 x <- c(.18,.46,.76)[i];arrow(x,.28,x,.24)
 d <- list(symptoms,proteins,pef_available)[[i]]
 title <- c('Symptom data available','MSD protein data available','PEF data available')[i]
 box(x,.19,.26,.10,label_counts(d,title),bold=TRUE)
 reason <- c('Symptom diary not submitted','Samples awaiting MSD analysis','No usable baseline/follow-up PEF')[i]
 box(x,.08,.26,.09,sprintf('Not available: %d/%d\n%s',nrow(groups)-nrow(d),nrow(groups),reason),fill='#D0D3D4')
}
for(i in 1:4) {
 y<-c(.93,.71,.48,.19)[i]
 grid::grid.roundrect(.03,y,.035,c(.13,.14,.24,.13)[i],r=grid::unit(.05,'inches'),gp=grid::gpar(fill='#244F8D',col=NA))
 grid::grid.text(c('Enrolment','Experimental\ninfection','Infection assessment','Analysis')[i],.03,y,rot=90,
                 gp=grid::gpar(col='white',fontfamily='Arial',fontface='bold',fontsize=11))
}
invisible(dev.off())
