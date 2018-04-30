##################################################
# Simulate outcomes of sequence of increasingly
# aggressive screening/treatment policies in Uganda
##################################################
library(bcimodel)
library(reshape)
library(scales)
library(grid)
library(ggplot2)

set.seed(98103)

#datestamp <- '2018-04-13'
datestamp <- '2018-04-19'

plotpath <- file.path('..', '..', '..', '..', 'plots')

##################################################
# Convert k-year survival to exponential rate
##################################################
exp.rate <- function(survival, year=5) -log(survival)/year

##################################################
# Format/initialize treatment distribution model input
##################################################
create_treatment_distribution <- function(nh,
                                          pol,
                                          treat_hrs=c('None'=1, 'Tamoxifen'=0.7),
                                          treat_probs=c('None'=1, 'Tamoxifen'=0),
                                          scale_by_subgroup=FALSE){
    ssids <- with(nh, paste(stage, subgroup, sep='.'))
    tset <- expand.grid(txHR=treat_hrs,
                        SSid=ssids,
                        stringsAsFactors=FALSE)
    tset <- with(tset, data.frame(SSno=as.integer(factor(SSid, levels=unique(SSid))),
                                  SSid,
                                  txSSno=seq_len(nrow(tset)),
                                  txSSid=names(treat_hrs),
                                  txHR,
                                  stringsAsFactors=FALSE))
    tprobs <- data.frame(txSSid=names(treat_probs),
                         treat_probs,
                         row.names=NULL,
                         stringsAsFactors=FALSE)
    tprobs <- reshape::cast(merge(tprobs, data.frame(policy=pol$id)),
                            txSSid~policy,
                            value='treat_probs')
    merged <- merge(tset, tprobs, by='txSSid')
    merged <- transform(merged, subgroup=sub('^[^.]*[.](.*)$', '\\1', SSid))
    if(scale_by_subgroup){
        sgids <- plyr::ddply(nh, 'subgroup', plyr::summarize, prop=sum(prop))
        merged <- merge(merged, sgids, by='subgroup')
        for(policy in as.character(pol$id))
            merged[[policy]] <- merged[[policy]]*merged[['prop']]
    }
    merged <- reshape::sort_df(merged, vars='txSSno')
    merged <- merged[, c('SSno', 'SSid', 'txSSno', 'txSSid', 'txHR', as.character(pol$id))]
    return(merged)
}

##################################################
# Specify sequence of policies
##################################################
polseq <- list(pol=data.frame(num=seq(7),
                              id=c('BC', 'CA', 'CBE', 'M', 'ME', 'MEC', 'MECA'),
                              name=c('Base case',
                                     'Clinical access',
                                     'CBE',
                                     'Mammography',
                                     'Mammography/ET for ER+',
                                     'Mammography/ET for ER+/Chemo for ER-',
                                     'Mammography/ET for ER+/Chemo for ER-/Chemo for adv ER+'),
                              pairnum=c(NA, rep(1, 6)),
                              earlydetHR=c(1, 0.60/0.78, 0.35/0.78, rep(0.30/0.78, 4)),
                              stringsAsFactors=FALSE),
               nh=compile_naturalhist(prop_adv=0.78,
                                      mortrates=c(Early=exp.rate(0.69),
                                                  Advanced=exp.rate(0.35)),
                                      subgroup_probs=c('ER-'=1-0.41, 'ER+'=0.41)))
polseq$map <- create_stageshift_map(polseq$nh)
polseq$tx <- with(polseq, create_treatment_distribution(nh, pol,
                                                        treat_hrs=c('None'=1,
                                                                    'Tamoxifen'=0.7,
                                                                    'Chemo'=0.775,
                                                                    'Tamoxifen+Chemo'=0.5425),
                                                        treat_probs=c('None'=1,
                                                                      'Tamoxifen'=0,
                                                                      'Chemo'=0,
                                                                      'Tamoxifen+Chemo'=0)))
polseq$tx <- within(polseq$tx, {
                         txHR[grepl('ER[-]', SSid) & txSSid == 'Tamoxifen'] <- 1
                         txHR[grepl('ER[-]', SSid) & txSSid == 'Tamoxifen+Chemo'] <- 1
                         ME[grepl('ER[+]', SSid) & txSSid == 'None'] <- 0
                         ME[grepl('ER[+]', SSid) & txSSid == 'Tamoxifen'] <- 1
                         MEC[grepl('ER[+]', SSid) & txSSid == 'None'] <- 0
                         MEC[grepl('ER[+]', SSid) & txSSid == 'Tamoxifen'] <- 1
                         MEC[grepl('ER[-]', SSid) & txSSid == 'None'] <- 0
                         MEC[grepl('ER[-]', SSid) & txSSid == 'Chemo'] <- 1
                         MECA[grepl('ER[+]', SSid) & txSSid == 'None'] <- 0
                         MECA[SSid == 'Early.ER+' & txSSid == 'Tamoxifen'] <- 1
                         MECA[SSid == 'Advanced.ER+' & txSSid == 'Tamoxifen+Chemo'] <- 1
                         MECA[grepl('ER[-]', SSid) & txSSid == 'None'] <- 0
                         MECA[grepl('ER[-]', SSid) & txSSid == 'Chemo'] <- 1
                            })

##################################################
# Confirm treatment distributions sum to 1
##################################################
plyr::d_ply(polseq$tx,
            'SSid',
            with,
            stopifnot(sum(BC) == 1 &
                      sum(CA) == 1 &
                      sum(CBE) == 1 &
                      sum(M) == 1 &
                      sum(ME) == 1 &
                      sum(MEC) == 1 &
                      sum(MECA) == 1))

##################################################
# Confirm treatment distributions sum to 1
##################################################
plyr::d_ply(polseq$tx,
            'SSid',
            with,
            stopifnot(sum(BC) == 1 &
                      sum(CA) == 1 &
                      sum(CBE) == 1 &
                      sum(M) == 1 &
                      sum(ME) == 1 &
                      sum(MEC) == 1 &
                      sum(MECA) == 1))

##################################################
# Extend sequence of policies to include HER2
##################################################
nh <- read.csv('natural_history_inputs_2018-04-18.csv')
nh <- rename(nh, c(Stage='stage', Proportion='prop'))
nh <- transform(nh, subgroup=paste0('ER', ER, 'HER2', HER2))
nh <- subset(nh, select=-c(ER, HER2))
subgroup_dummy <- with(subset(nh, stage == 'Advanced'), prop/sum(prop))
names(subgroup_dummy) <- with(subset(nh, stage == 'Advanced'), subgroup)

polseq <- list(pol=data.frame(num=seq(8),
                              id=c('BC', 'CA', 'CBE', 'M', 'ME', 'MEC', 'MECA', 'MECAH'),
                              name=c('Base case',
                                     'Clinical access',
                                     'CBE',
                                     'Mammography',
                                     'Mammography/ET for ER+',
                                     'Mammography/ET for ER+/Chemo for ER-',
                                     'Mammography/ET for ER+/Chemo for ER-/Chemo for adv ER+',
                                     'Mammography/ET for ER+/Chemo for ER-/Chemo for adv ER+/Trastuzumab for HER2+'),
                              pairnum=c(NA, rep(1, 7)),
                              earlydetHR=c(1, 0.60/0.78, 0.35/0.78, rep(0.30/0.78, 5)),
                              stringsAsFactors=FALSE),
               nh=compile_naturalhist(prop_adv=0.78,
                                      mortrates=c(Early=exp.rate(0.69),
                                                  Advanced=exp.rate(0.35)),
                                      subgroup_probs=subgroup_dummy))
polseq$nh <- merge(subset(polseq$nh, select=-prop),
                   nh,
                   by=c('stage', 'subgroup'),
                   sort=FALSE)[, c('prop', 'stage', 'subgroup', 'mortrate')]
class(polseq$nh) <- c('data.frame', 'naturalhist')
polseq$map <- create_stageshift_map(polseq$nh)
polseq$tx <- with(polseq, create_treatment_distribution(nh, pol,
                                            treat_hrs=c('None'=1,
                                                        'Tamoxifen'=0.7,
                                                        'Chemo'=0.775,
                                                        'Tamoxifen+Chemo'=0.5425,
                                                        'Trastuzumab'=0.7,
                                                        'Chemo+Trastuzumab'=0.5425,
                                                        'Tamoxifen+Trastuzumab'=0.49,
                                                        'Tamoxifen+Chemo+Trastuzumab'=0.37975),
                                            treat_probs=c('None'=1,
                                                          'Tamoxifen'=0,
                                                          'Chemo'=0,
                                                          'Tamoxifen+Chemo'=0,
                                                          'Trastuzumab'=0,
                                                          'Chemo+Trastuzumab'=0,
                                                          'Tamoxifen+Trastuzumab'=0,
                                                          'Tamoxifen+Chemo+Trastuzumab'=0)))
polseq$tx <- within(polseq$tx, {
         txHR[grepl('ER[-]', SSid) & txSSid == 'Tamoxifen'] <- 1
         txHR[grepl('HER2[-]', SSid) & txSSid == 'Trastuzumab'] <- 1
         txHR[grepl('ER[-]', SSid) & txSSid == 'Tamoxifen+Chemo'] <- 0.775
         txHR[grepl('HER2[-]', SSid) & txSSid == 'Chemo+Trastuzumab'] <- 0.775
         txHR[grepl('ER[-]HER2[-]', SSid) & txSSid == 'Tamoxifen+Trastuzumab'] <- 1
         txHR[grepl('ER[-]HER2[+]', SSid) & txSSid == 'Tamoxifen+Trastuzumab'] <- 0.7
         txHR[grepl('ER[+]HER2[-]', SSid) & txSSid == 'Tamoxifen+Trastuzumab'] <- 0.7
         txHR[grepl('ER[-]HER2[-]', SSid) & txSSid == 'Tamoxifen+Chemo+Trastuzumab'] <- 0.775
         txHR[grepl('ER[-]HER2[+]', SSid) & txSSid == 'Tamoxifen+Chemo+Trastuzumab'] <- 0.5425
         txHR[grepl('ER[+]HER2[-]', SSid) & txSSid == 'Tamoxifen+Chemo+Trastuzumab'] <- 0.5425
           })
polseq$tx <- within(polseq$tx, {
                         ME[grepl('ER[+]', SSid) & txSSid == 'None'] <- 0
                         ME[grepl('ER[+]', SSid) & txSSid == 'Tamoxifen'] <- 1
                         MEC[grepl('ER[+]', SSid) & txSSid == 'None'] <- 0
                         MEC[grepl('ER[+]', SSid) & txSSid == 'Tamoxifen'] <- 1
                         MEC[grepl('ER[-]', SSid) & txSSid == 'None'] <- 0
                         MEC[grepl('ER[-]', SSid) & txSSid == 'Chemo'] <- 1
                         MECA[grepl('ER[+]', SSid) & txSSid == 'None'] <- 0
                         MECA[grepl('Early.ER[+]', SSid) & txSSid == 'Tamoxifen'] <- 1
                         MECA[grepl('Advanced.ER[+]', SSid) & txSSid == 'Tamoxifen+Chemo'] <- 1
                         MECA[grepl('ER[-]', SSid) & txSSid == 'None'] <- 0
                         MECA[grepl('ER[-]', SSid) & txSSid == 'Chemo'] <- 1
                         MECAH[grepl('ER[+]', SSid) & txSSid == 'None'] <- 0
                         MECAH[grepl('HER2[+]', SSid) & txSSid == 'None'] <- 0
                         MECAH[SSid == 'Early.ER+HER2-' & txSSid == 'Tamoxifen'] <- 1
                         MECAH[SSid == 'Advanced.ER+HER2-' & txSSid == 'Tamoxifen+Chemo'] <- 1
                         MECAH[SSid == 'Early.ER+HER2+' & txSSid == 'Tamoxifen+Trastuzumab'] <- 1
                         MECAH[SSid == 'Advanced.ER+HER2+' & txSSid == 'Tamoxifen+Chemo+Trastuzumab'] <- 1
                         MECAH[grepl('ER[-]HER2[-]', SSid) & txSSid == 'None'] <- 0
                         MECAH[grepl('ER[-]HER2[-]', SSid) & txSSid == 'Chemo'] <- 1
                         MECAH[grepl('ER[-]HER2[+]', SSid) & txSSid == 'None'] <- 0
                         MECAH[grepl('ER[-]HER2[+]', SSid) & txSSid == 'Chemo+Trastuzumab'] <- 1
                            })

##################################################
# Confirm treatment distributions sum to 1
##################################################
plyr::d_ply(polseq$tx,
            'SSid',
            with,
            stopifnot(sum(BC) == 1 &
                      sum(CA) == 1 &
                      sum(CBE) == 1 &
                      sum(M) == 1 &
                      sum(ME) == 1 &
                      sum(MEC) == 1 &
                      sum(MECA) == 1 &
                      sum(MECAH) == 1))

##################################################
# Tabulate results
##################################################
#knitr::kable(uganda_polseq[['10']])

##################################################
# Visualize results
##################################################
gg_theme <- function(...){
    theme_set(theme_bw())
    theme_update(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.border=element_blank(),
                 axis.line=element_line(colour='black', size=0.5),
                 axis.title=element_text(angle=0, size=16),
                 axis.text=element_text(size=12),
                 axis.text.x=element_text(size=12),
                 strip.background=element_rect(fill=NA, colour=NA),
                 strip.text.x=element_text(angle=0, size=14),
                 strip.text.y=element_text(angle=0, size=14),
                 legend.key=element_rect(colour=NA))
    theme_update(...)
}

seqplot <- function(dset,
                    outcome='Cumulative BC Mortality',
                    fu=10,
                    ext='png',
                    saveit=FALSE){
    ylabel <- switch(outcome,
                     'Cumulative BC Mortality'='Mortality per 100,000 women\nages 30-69',
                     'Years of Life Saved'='Years of life saved per 100,000 women\nages 30-69')
    ylabel <- paste(ylabel, 'after', fu, 'years\n')
    ymax <- switch(outcome,
                   'Cumulative BC Mortality'=450,
                   'Years of Life Saved'=650)
    results <- dset[[as.character(fu)]]
    results <- data.frame(Policy=colnames(results$mean),
                          Mean=results$mean[outcome, ],
                          Lower=results$lower[outcome, ],
                          Upper=results$upper[outcome, ],
                          row.names=NULL)
    results <- transform(results, Policy=gsub('[/]', '\n', Policy))
    results <- transform(results, Policy=factor(Policy,
                                                levels=c('Base case',
                                                         'Clinical access',
                                                         'CBE',
                                                         'Mammography',
                                                         'Mammography\nET for ER+',
                                                         'Mammography\nET for ER+\nChemo for ER-',
                                                         'Mammography\nET for ER+\nChemo for ER-\nChemo for adv ER+',
                                                         'Mammography\nET for ER+\nChemo for ER-\nChemo for adv ER+\nTrastuzumab for HER2+')))
    gg_theme()
    gg <- ggplot(data=results)
    gg <- gg+geom_bar(aes(x=Policy, y=Mean),
                      colour=muted('skyblue'),
                      fill='skyblue',
                      alpha=0.4,
                      position='dodge',
                      stat='identity')
    gg <- gg+geom_errorbar(aes(x=Policy, ymin=Lower, ymax=Upper),
                           colour=muted('skyblue'),
                           alpha=0.6,
                           position='dodge')
    gg <- gg+geom_text(aes(x=Policy,
                           y=Mean/2,
                           label=ifelse(Mean > 0,
                                        sprintf('%3.0f', Mean),
                                        '')),
                       size=5)
    gg <- gg+geom_blank(aes(y=0))
    gg <- gg+scale_x_discrete(name='')
    gg <- gg+scale_y_continuous(name=ylabel,
                                limits=c(0, ymax),
                                breaks=seq(0, ymax, by=100),
                                expand=c(0, 0))
    print(gg)
    if(saveit){
        prefix <- 'bcimodel_uganda'
        midfix <- make.names(tolower(outcome))
        filename <- paste(prefix, midfix, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               height=8,
               width=14)
    }
}

##################################################
# Execute and visualize results
##################################################
execute_sequence <- function(polseq, saveit=FALSE){
    uganda_polseq <- with(polseq, simpolicies(pol, nh, tx,
                                              minage=30,
                                              maxage=69,
                                              futimes=c(10, 20),
                                              popsize=100000,
                                              returnstats=c('mean', 'lower', 'upper'),
                                              sims=100))
    seqplot(uganda_polseq, outcome='Cumulative BC Mortality', saveit=saveit)
    seqplot(uganda_polseq, outcome='Years of Life Saved', saveit=saveit)
}
#execute_sequence(polseq, saveit=TRUE)

##################################################
# Visualize natural history inputs
##################################################
nh <- read.csv('natural_history_inputs_2018-04-18.csv')
nh <- rename(nh, c(Stage='stage', Proportion='prop'))
nh <- transform(nh, merged.label=paste(paste0('ER', ER),
                                       paste0('HER2', HER2),
                                       sep='/'))

nhplot <- function(nh, ext='png', saveit=FALSE){
    nh <- transform(nh, stage=factor(stage,
                                     levels=c('Early', 'Advanced'),
                                     labels=c('\nEarly', '\nAdvanced')))
    gg_theme(legend.position='none',
             axis.ticks=element_blank(),
             axis.text.x=element_text(size=18),
             axis.text.y=element_blank(),
             axis.line=element_blank())
    gg <- ggplot(data=nh)
    gg <- gg+geom_tile(aes(x=stage,
                           y=merged.label,
                           width=prop/max(prop),
                           height=prop/max(prop)),
                       fill='gray85')
    gg <- gg+geom_text(aes(x=stage,
                           y=merged.label,
                           label=merged.label,
                           size=prop/max(prop)),
                       colour='purple')
    gg <- gg+scale_x_discrete(name='', expand=c(0, 0))
    gg <- gg+scale_y_discrete(name='', expand=c(0, 0))
    print(gg)
    if(saveit){
        prefix <- 'bcimodel_uganda_natural_history'
        filename <- paste(prefix, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               height=3,
               width=5)
    }
}
#nhplot(nh, saveit=TRUE)

##################################################
# Visualize survival and early detection inputs
##################################################
survfun <- function(x, rate=1) 1-pexp(x, rate=rate)

survplot <- function(nh, ext='png', saveit=FALSE){
    dset <- plyr::ddply(nh,
                        'stage',
                        plyr::summarize,
                        rate=unique(mortrate))
    dset <- plyr::ddply(dset,
                        'stage',
                        transform,
                        x=seq(0, 10, by=0.1))
    dset <- plyr::ddply(dset,
                        'stage',
                        transform,
                        survival=survfun(x, rate=rate))
    gg_theme(legend.position='none')
    gg <- ggplot(data=dset)
    gg <- gg+geom_line(aes(x=x,
                           y=survival,
                           group=stage,
                           colour=stage),
                       size=0.75)
    gg <- gg+geom_text(data=subset(dset, x == 10),
                       aes(x=x,
                           y=survival,
                           label=stage,
                           colour=stage),
                       hjust=1,
                       vjust=1.2,
                       size=6)
    gg <- gg+scale_x_continuous(name='\nYears since diagnosis',
                                limits=c(0, 11),
                                breaks=seq(0, 10, by=2),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Survival\n',
                                limits=c(0, 1),
                                label=percent_format(),
                                expand=c(0, 0))
    gg <- gg+scale_colour_manual(values=c(Early='darkgreen', Advanced='orange'))
    print(gg)
    if(saveit){
        prefix <- 'bcimodel_uganda_survival'
        filename <- paste(prefix, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               height=3,
               width=5)
    }
}
#survplot(polseq$nh, saveit=TRUE)

##################################################
# Visualize treatment efficacy inputs
##################################################
txplot <- function(tx, ext='png', saveit=FALSE){
    tx <- subset(tx, select=c(SSid, txSSid, txHR))
    tx <- transform(tx,
                    stage=sub('^([^.]*)[.](ER[+-])(HER2[+-])$', '\\1', SSid),
                    ER=sub('^([^.]*)[.](ER[+-])(HER2[+-])$', '\\2', SSid),
                    HER2=sub('^([^.]*)[.](ER[+-])(HER2[+-])$', '\\3', SSid))
    casted <- cast(tx, stage+ER+HER2~txSSid, value='txHR')
    melted <- melt(casted, id.vars=c('stage', 'ER', 'HER2'))
    melted <- rename(melted, c(txSSid='treatment', value='HR'))
    melted <- transform(melted,
                        treatment=factor(treatment, levels=c('None',
                                                             'Tamoxifen',
                                                             'Chemo',
                                                             'Tamoxifen+Chemo',
                                                             'Trastuzumab',
                                                             'Tamoxifen+Trastuzumab',
                                                             'Chemo+Trastuzumab',
                                                             'Tamoxifen+Chemo+Trastuzumab')),
                        stage=factor(stage, levels=c('Early', 'Advanced')))
    gg_theme(legend.position='none',
             panel.spacing=unit(0.01, 'npc'),
             axis.ticks=element_blank(),
             axis.text.x=element_text(size=10, angle=90, hjust=1, vjust=0.5))
    gg <- ggplot(data=melted)
    gg <- gg+geom_rect(aes(xmin=treatment,
                           xmax=treatment,
                           ymin=HR,
                           ymax=1,
                           colour=treatment),
                       size=1)
    gg <- gg+geom_hline(yintercept=1, colour='darkgray')
    gg <- gg+facet_grid(HER2+ER~stage)
    gg <- gg+scale_x_discrete(name='')
    gg <- gg+scale_y_continuous(name='Hazard ratio\n',
                                limits=c(0.2, 1.2),
                                breaks=seq(0.4, 1, by=0.2),
                                expand=c(0, 0))
    print(gg)
    if(saveit){
        prefix <- 'bcimodel_uganda_treatment_efficacy'
        filename <- paste(prefix, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               height=6,
               width=9)
    }
}
#txplot(polseq$tx, saveit=TRUE)

