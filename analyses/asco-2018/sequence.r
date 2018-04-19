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

datestamp <- '2018-04-13'
plotpath <- file.path('..', '..', '..', 'plots')

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
# Collect results
##################################################
if(FALSE)
    uganda_polseq <- with(polseq, simpolicies(pol, nh, tx,
                                              minage=30,
                                              maxage=69,
                                              futimes=c(10, 20),
                                              popsize=100000,
                                              returnstats=c('mean', 'lower', 'upper'),
                                              sims=100))

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
                   'Years of Life Saved'=550)
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
                                                         'Mammography\nET for ER+\nChemo for ER-\nChemo for adv ER+')))
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
               height=5,
               width=12)
    }
}
seqplot(uganda_polseq, outcome='Cumulative BC Mortality', saveit=TRUE)
seqplot(uganda_polseq, outcome='Years of Life Saved', saveit=TRUE)

