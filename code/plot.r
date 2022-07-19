isocode <- 'MWI'
here::i_am("code/plot.r")
library(here)
dir.create(here('fig', isocode))

#+ include = FALSE
knitr::opts_chunk$set(
    include = FALSE,
    message = FALSE,
    warnings = FALSE
)
library(data.table)
devtools::load_all("~/Github/eppasm")
library(tidyverse)
devtools::load_all("~/Code/R/ktools/")
char <- ktools::char
options(ggplot2.discrete.fill = okabe, ggplot2.discrete.colour = okabe)

# Metavars

ago <- c('15-16', '17-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49')
agdb <- c(15:30, "31-34", "35-39", "40-44", "45-49")
agr <- c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49')
young <- c('15-19', '20-24', '25-29')

ktheme <- theme(
    axis.title.y.right = element_text(angle = 0, vjust = 1.05, hjust = 1, margin = margin(l = -20)),
    axis.text.x = element_text(size = 7, angle = 0),
    axis.text.y = element_text(size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey95"),
    strip.background = element_rect(fill = "#FCF4DC", color = NA),
    strip.text = element_text(face = "bold")
)

run_a_model <- function(fit, name = "epp", vs = "R") {
    fp <- fit$fits[[name]]$fp
    pa <- fit$fits[[name]]$par %>% eppasm:::as_parametter_list()
    fp <- modifyList(fp, parameters(pa, fp, "transform"))
    fp$VERSION <- vs
    list(mod = simmod(fp), fp = fp, pa = pa)
}

#' # Outputs:
res <- readRDS(here("fit/res.rds"))

ipd <- readRDS(system.file('extdata', 'inputs_nat.rds', package = 'eppasm'))
iso <- tibble(
    iso2 = names(ipd),
    iso3 = countrycode::countrycode(iso2, "iso2c", "iso3c")
)  %>% 
drop_na()

names(res)  <- iso$iso3

# Data and simulated
likedat <- res[[isocode]]$fits[[paste0('EPP_', isocode)]]$likdat

mwi_epp <- run_a_model(res$MWI, "epp_MWI")
mwi_edb <- run_a_model(res$MWI, "eppdb_MWI")
# Basic r_t
tibble(
    model_EPP = epp$fp$rvec, 
    model_EBD = edb$fp$rvec, 
    time = seq_along(epp$fp$rvec)
) %>% 
pivot_longer(-time, names_sep = '_', names_to = c(NA, 'model')) %>% 
ggplot() + 
    geom_line(aes(time, value, color = model))  +
    scale_y_continuous(trans = 'log', breaks = c(0, .1, .2, .3, .4, .5)) +
    scale_x_continuous(labels = function(x) x / 10 + 1970) +
    ktheme +
    labs(title = 'Transmission rate over time', x = 'Year', y = 'r(t)') +
    theme(legend.position = c(.9, .9))

savePNG(here("fig", isocode, 'r_t_estimate.png'), 7, 4)

#' FRR in EPP vs EDB
frrcd4 <- epp$fp$frr_cd4  %>% 
as.data.table() %>% 
filter(V3 == 35)

edb$fp$frr_cd4  %>% 
as.data.table() %>% 
filter(V3 == 35) %>% 
bind_rows(frrcd4, .id = 'model') %>% 
mutate(model = char(EDB, EPP)[as.numeric(model)]) %>%
ggplot() + 
    geom_line(aes(V2, value, color = factor(V1))) +
    labs(title = 'over agegr') +
    facet_grid(~model, scales = 'free')

#' ## Prevalence general population by age-group
sim_prevagr <- calc_prev_agegr(epp$mod, seq(15, 55, 5)) %>%
    bind_rows(
        calc_prev_agegr(edb$mod, seq(15, 55, 5)),
        .id = "model"
    ) %>%
    mutate(model = char(EPP, EDB)[as.numeric(model)])

#' ## Prevalence in pregnant women over time
pregprev1549 <- agepregprev(
    epp$mod, epp$fp,
    aidx = 1, agspan = 35, yidx = 1:48, expand = TRUE
)
pregprev1549db <- agepregprev(
    edb$mod, edb$fp,
    aidx = 1, agspan = 35, yidx = 1:48, expand = TRUE
)

likedat$ancrtcens.dat  |> 
    mutate(
        ul = prev + 1.96 * sqrt(v.ancrt),
        ll = prev - 1.96 * sqrt(v.ancrt)) |> 
    ggplot() +
    geom_line(aes(year, prev, color = site), alpha = .5, data = likedat$ancsite.dat$df) +
    geom_line(aes(year, prev, color = site, linetype = site), size = 1.2, data = tibble(
        year = 1970:2017, 
        prev = pregprev1549db, 
        site = "National DB"
    )) +
    geom_line(aes(year, prev, color = site, linetype = site), size = 1.2, data = tibble(
        year = 1970:2017, 
        prev = pregprev1549, 
        site = "National EPP"
    )) +
    geom_pointrange(aes(year, prev, ymin = ll, ymax = ul), size = .7) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = "Fitted to HIV prevalence in ANC sites and census", 
        linetype = "Model", y = "Prevalence"
    ) +
    theme(legend.position = c(.1, .9)) +
    guides(color = 'none')

savePNG(here('fig', isocode, 'ANC_National_Prev'), 7, 7)
 
#' ## Prevalence in pregnant women over time by age-group
#'
pregprev1549_agr <- agepregprev(
    epp$mod, epp$fp,
    aidx = seq(1, 35, 5), agspan = 5, yidx = 1:48, expand = TRUE
)

pregprev1549db_agr <- agepregprev(
    edb$mod, edb$fp,
    aidx = seq(1, 35, 5), agspan = 5, yidx = 1:48, expand = TRUE
)

tibble(
    prev_EPP = pregprev1549_agr,
    prev_EDB = pregprev1549db_agr,
    agegr = rep(1:7, times = 48),
    year = rep(1969 + 1:48, each = 7)
) %>%
allot(pregprev_tb)

#' ## Combining general prevalence and pregmant

prevboth <- bind_rows(
    sim_prevagr,
    pregprev_tb  %>% 
    pivot_longer(-c(agegr, year), names_sep = '_', names_to = c(NA, 'model'), values_to = 'prev') %>% 
    mutate(sex = 'female', agegr = agr[agegr]), 
    .id = 'pop'
) %>% 
mutate(pop = char('General', "Pregnant")[as.numeric(pop)]) %>% 
mutate(pop = if_else(pop == 'General', paste(pop, '-', sex), pop))

hhs_dt <- likedat$hhs |> filter(agegr %in% young) %>% 
rename(prev_dt = prev) %>% 
mutate(pop = paste('General -', sex)) %>% 
select(year, prev_dt, ci_u, ci_l, pop, agegr)

prevboth %>%
filter(agegr %in% young) %>%
left_join(hhs_dt) %>% 
mutate(pop = factor(pop, levels = c('General - male', 'General - female', 'Pregnant'))) %>% 
ggplot() +
    geom_line(aes(year, prev, color = pop, linetype = model)) +
    facet_grid(vars(agegr), vars(pop), switch = 'y') +
    scale_y_continuous(labels = scales::percent, position = 'right') +
    geom_pointrange(aes(year, prev_dt, ymin = ci_l, ymax = ci_u, color = pop), size = .3) +
    labs(title = 'Prevalence by age-group and population', y = 'Prevalence') +
    guides(color = 'none', linetype = guide_legend(direction = 'horizontal', )) +
    ktheme +
    theme(legend.position = c(.95, 1.07)) +
    theme(
        axis.title.y.right = element_text(size = 8,
        angle = 0, vjust = 1, hjust = 1, margin = margin(l = -20)))

savePNG(here('fig', isocode, 'combined_preg_general'), 7, 6)

#' ## Prevalence in non-pregnant women over time, by age group 15-19, 20-24
# how to do this?

#' ## Prevalence in sexually active/non-sexually active over time, age 15-19, 20-24
#'
mwi30 <- edb$mod$data[1:16, , , ]
mwi30a <- mwi30 - edb$mod$vpop

mwi30a %<>%
    as.data.table() %>%
    rename_with(~ c("age", "sex", "hiv", "year", "size")) %>%
    mutate(sexual='active')

mwi30n <- edb$mod$vpop %>%
    as.data.table() %>%
    rename_with(~ c("age", "sex", "hiv", "year", "size")) %>%
    mutate(sexual = "nonactive")

mwi30a %>%
    bind_rows(mwi30n) %>%
    pivot_wider(names_from = hiv, values_from = size) %>%
    mutate(
        age = age + 14,
        agr = findInterval2(age, c(15, 20, 25, 30)), 
        sex = c('male', 'female')[sex]
    ) %>%
    filter(age != 30) %>%
    group_by(sex, year, sexual, agr) %>%
    summarise(prev = sum(`2`) / sum(`1` + `2`)) %>%
    allot(mwi30)

#+ prev_active_mwi, fig.cap = 'Prevalence in sexually active and non-active', include=TRUE
#
mwi30 %>%
    ggplot(aes(year, prev, color = sex)) +
    facet_grid(vars(sexual), vars(agr), scales = 'free', switch = 'y') +
    geom_line() +
    scale_y_continuous(labels = scales::percent, position = 'right') +
    scale_x_continuous(labels = \(x) x + 1969) +
    labs(
        title = "HIV prevalence by sexual status and age-group",
        y = "Prevalence", x = 'Year', color = ''
    ) +
    ktheme +
    theme(legend.position = c(1, 1.13)) +
    guides(color = guide_legend(direction = 'horizontal')) +
    coord_cartesian(expand = T) +
    scale_color_manual(values = okabe)

savePNG(here('fig', isocode, 'prev_agr_sex_status'), 7, 4)

# Proportion of sexually active among HIV positive

mwi30a %>%
    bind_rows(mwi30n) %>% 
    filter(hiv == 2) %>% 
    pivot_wider(names_from = sexual, values_from = size) %>% 
    mutate(
        age = age + 14,
        agr = findInterval2(age, c(15, 20, 25, 30)), 
        sex = c('male', 'female')[sex], 
        pa = active / (active + nonactive),
    ) %>% 
    filter(is.finite(pa)) %>% 
    filter(year == max(year)) %>% 
    filter(age < 25) %>%
    allot(pos_active)

pos_active %>% 
    ggplot() + 
        geom_col(aes(age, pa, fill = sex), position = 'dodge') +
        ktheme +
        labs(
            title = 'Proportion of sexually active among HIV positive', 
            y = '', fill = ''
        ) +
        scale_y_continuous(labels = scales::percent) +
        scale_x_continuous(breaks = 15:24) +
        theme(legend.position = c(0.07, .95))

savePNG(here('fig', isocode, 'hivpos_active.png'))

# Proportion newly sexual debutted among HIV positive active
db_long <- edb$fp$db_rate  %>%
as.data.table(1)  %>% 
rename_with(~char(age, sex, year, rate))

mwi30n  %>% 
left_join(db_long)  %>% 
filter(hiv == 2) %>% 
mutate(n_next = size * rate, year = year + 1)  %>% 
select(age, sex, year, n_next) %>% 
allot(newly_debut_pos)

edb$mod$data %>% unclass() %>% 
as.data.table(1) %>% 
rename_with(~char(age, sex, hiv, year, n)) %>% 
filter(hiv == 2, age <= 16) %>% 
left_join(newly_debut_pos) %>% 
mutate(p_db = n_next / n)  %>% 
allot(p_db_pos)

p_db_pos  %>%
mutate(age = 14 + age, year = year + 1969, sex = char(male, female)[sex]) %>% 
filter(age < 20, year %in% c(2005, 2015)) %>%
ggplot() + 
    geom_col(aes(age, p_db, fill = factor(sex)), position = 'dodge') +
    facet_grid(cols = vars(year)) + 
    labs(title = 'Proportion of newly sexual debuted HIV positive', y = '', fill='') +
    ktheme +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = c(0.07, .95))

savePNG(here('fig', isocode, 'hiv_pos_new_debut'), 7, 4)

#' ## HIV+ : HIV- fertility rate/ratio over time
#'
agr_i <- c(15, 17, seq(20, 50, 5))

source(here('code/FRR.r'))
epp_frr <- FRR(epp)
edb_frr <- FRR(edb)

bind_rows(
    epp_frr$frr %>% as.data.table(1) %>% pivot_longer(-rn),
    edb_frr$frr %>% as.data.table(1) %>% pivot_longer(-rn),
    .id = "model"
) %>%
    mutate(model = char(EPPASM, EPPdb)[as.numeric(model)]) %>%
    type.convert() %>%
    filter(rn %in% c("15-19", "20-24", "25-29")) %>%
    ggplot(aes(name, value, color = model)) +
    facet_wrap(~rn) +
    geom_line() +
    labs(color = 'Model', x = "Year", y = "Fertility rate ratio") +
    ktheme +
    theme(legend.position = c(.92, .85))

savePNG(here("fig", isocode, "FRR.png"), 7, 4)

edb_frr$asfr_p  %>% as.data.table(1) %>% 
pivot_longer(-rn) %>% 
mutate(name = as.numeric(name), group = 'Positive', model = 'EDB') %>% 
rename(year = name, FR = value, agegr = rn) %>% 
allot(frr_long)

epp_frr$asfr_p  %>% as.data.table(1) %>% 
pivot_longer(-rn) %>% 
mutate(name = as.numeric(name), group = 'Positive', model = 'EPP') %>% 
rename(year = name, FR = value, agegr = rn) %>% 
bind_rows(frr_long) %>% 
allot(frr_long)

epp_frr$asfr_n  %>% as.data.table(1) %>% 
pivot_longer(-rn) %>% 
mutate(name = as.numeric(name), group = 'Negative', model = 'EPP') %>% 
rename(year = name, FR = value, agegr = rn) %>% 
bind_rows(frr_long) %>% 
allot(frr_long)

edb_frr$asfr_n  %>% as.data.table(1) %>% 
pivot_longer(-rn) %>% 
mutate(name = as.numeric(name), group = 'Negative', model = 'EDB') %>% 
rename(year = name, FR = value, agegr = rn) %>% 
bind_rows(frr_long) %>% 
allot(frr_long)

frr_long  %>% 
filter(agegr %in% young) %>% 
ggplot() + 
    geom_line(aes(year, FR, color = model, linetype = group)) +
    facet_grid(vars(agegr), vars(model), switch = 'y') +
    guides(color = 'none', linetype = guide_legend(direction = 'horizontal')) +
    labs(
        title = 'Estimate fertility rate by age-group',
        linetype = '', 
        y = 'Fertility rate'
    ) +
    ktheme +
    scale_y_continuous(position = 'right') +
    theme(
        legend.position = c(.9, 1.08),
        axis.title.y.right = element_text(angle = 0, vjust = 1, hjust = 1, margin = margin(l = -15), size = 8),
    )

savePNG(here('fig', isocode, 'FR_by_mod_agr_status'), 7, 6)

#' ## HIV incidence rate ratio parameter estimates between the two models

bind_rows(
    epp$fp$incrr_age[, , 1] %>%
        as_tibble() %>%
        rename_with(~ c("male", "female")) %>%
        mutate(model = "EPP", age  = 15:80),
    edb$fp$incrr_age[, , 1] %>%
        as_tibble() %>%
        rename_with(~ c("male", "female")) %>%
        mutate(model = "EDB", age = 15:80)
) %>%
    pivot_longer(-c("model", "age")) %>%
    filter(age < 30) %>%
        ggplot(aes(age, value, color = model)) +
        facet_wrap(~name) +
        geom_line() +
        ktheme +
        labs(
        title = "Estimate incidence rate ratio by age and sex",
        y = "Incidence rate ratio"
        ) +
        theme(
            legend.position = c(.94, 0.1)
        )
savePNG(here('fig', isocode, 'IRR_age'), 7, 4)

#' ## % HIV+ 15-19 and 20-24 who are long-term survivors