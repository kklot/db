#+ include = FALSE
knitr::opts_chunk$set(
    include = FALSE,
    message = FALSE,
    warnings = FALSE
)
here::i_am("code/plot.r")
library(here)
library(data.table)
devtools::load_all("~/Code/R/ktools/")
devtools::load_all("~/Github/eppasm")
library(tidyverse)
char <- ktools::char
options(ggplot2.discrete.fill = okabe, ggplot2.discrete.colour = okabe)

ago <- c('15-16', '17-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49')
agdb <- c(15:30, "31-34", "35-39", "40-44", "45-49")
agr <- c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49')

res <- readRDS(here("fit/res.rds"))

run_a_model <- function(fit, name = "epp", vs = "R") {
    fp <- fit$fits[[name]]$fp
    pa <- fit$fits[[name]]$par %>% eppasm:::as_parametter_list()
    fp <- modifyList(fp, parameters(pa, fp, "transform"))
    fp$VERSION <- vs
    list(mod = simmod(fp), fp = fp, pa = pa)
}

ktheme <- theme(
    axis.title.y.right = element_text(angle = 0, vjust = 1.05, hjust = 1, margin = margin(l = -20)),
    axis.text.x = element_text(size = 7, angle = 0),
    axis.text.y = element_text(size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey95"),
    strip.background = element_rect(fill = "#FCF4DC", color = NA),
    strip.text = element_text(face = "bold")
)

#' # Outputs:
#'
# Data and simulated
likedat <- res$MWI$fits$epp_MWI$likdat

mwi_epp <- run_a_model(res$MWI, "epp_MWI")
mwi_edb <- run_a_model(res$MWI, "eppdb_MWI")

#' FRR in EPP vs EDB
frrcd4 <- mwi_epp$fp$frr_cd4  %>% 
as.data.table() %>% 
filter(V3 == 35)

mwi_edb$fp$frr_cd4  %>% 
as.data.table() %>% 
filter(V3 == 35) %>% 
bind_rows(frrcd4, .id = 'model') %>% 
mutate(model = char(EDB, EPP)[as.numeric(model)]) %>%
ggplot() + 
    geom_line(aes(V2, value, color = factor(V1))) +
    labs(title = 'over agegr') +
    facet_grid(~model, scales = 'free')

#' ## Prevalence general population by age-group
sim_prevagr <- calc_prev_agegr(mwi_epp$mod, seq(15, 55, 5)) %>%
    bind_rows(
        calc_prev_agegr(mwi_edb$mod, seq(15, 55, 5)),
        .id = "model"
    ) %>%
    mutate(model = char(EPP, EDB)[as.numeric(model)])

young <- c('15-19', '20-24', '25-29')

sim_prevagr %>%
    filter(agegr %in% young) %>%
    ggplot(aes(year, prev, color = sex)) +
    facet_grid(vars(sex), vars(agegr), scales = "free_y") +
    geom_line(aes(linetype = model)) +
    geom_pointrange(aes(ymin = ci_l, ymax = ci_u), size = .2,
                    data = likedat$hhs |> filter(agegr %in% young)) +
    scale_colour_manual(values = okabe) +
    theme(axis.text.x = element_text(angle = 45)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Malawi - HIV age-group specific prevalence")

ggsave(here('fig/mw_agr_prev.png'), width = 7, height = 4)

#' ## Prevalence in pregnant women over time
#'
pregprev1549 <- agepregprev(
    mwi_epp$mod, mwi_epp$fp,
    aidx = 1, agspan = 35, yidx = 1:48, expand = TRUE
)
pregprev1549db <- agepregprev(
    mwi_edb$mod, mwi_edb$fp,
    aidx = 1, agspan = 35, yidx = 1:48, expand = TRUE
)

likedat$ancrtcens.dat  |> 
    filter(country == "Malawi")  |> 
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
        linetype = "Model"
    ) +
    guides(color = 'none')

savePNG(here('fig/MW_ANC_National_Prev'), 7, 7)
 
#' ## Prevalence in pregnant women over time by age-group
#'
pregprev1549_agr <- agepregprev(
    mwi_epp$mod, mwi_epp$fp,
    aidx = seq(1, 35, 5), agspan = 5, yidx = 1:48, expand = TRUE
)

pregprev1549db_agr <- agepregprev(
    mwi_edb$mod, mwi_edb$fp,
    aidx = seq(1, 35, 5), agspan = 5, yidx = 1:48, expand = TRUE
)

tibble(
    prev_EPP = pregprev1549_agr,
    prev_EDB = pregprev1549db_agr,
    agegr = rep(1:7, times = 48),
    year = rep(1969 + 1:48, each = 7)
) %>%
allot(pregprev_tb)

pregprev_tb |> 
    filter(agegr < 4) |>
    pivot_longer(
        -c("agegr", "year"),
        names_to = c(NA, "model"), names_sep = "_",
        values_to = "prev"
    ) |>
    mutate(agegr = agr[agegr]) |>
    ggplot() +
    geom_line(aes(year, prev, color = model)) +
    facet_wrap(~agegr) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = okabe) +
    labs(title = "Pregnant prevalence by age-group", y = "Prevalence")

savePNG(here('fig/MW_Pregnant_Age_Group'), 7, 3.5)

#' ## Combining general prevalence and pregmant

prevboth <- bind_rows(
    sim_prevagr %>% mutate(model = if_else(model == "EPPDB", "EDB", "EPP")),
    pregprev_tb  %>% 
    pivot_longer(-c(agegr, year), names_sep = '_', names_to = c(NA, 'model'), values_to = 'prev') %>% 
    mutate(sex = 'female', agegr = agr[agegr]), 
    .id = 'pop'
) %>% 
mutate(pop = char('General', "Pregnant")[as.numeric(pop)]) %>% 
mutate(pop = if_else(pop == 'General', paste(pop, '-', sex), pop))

geom_pointrange(aes(ymin = ci_l, ymax = ci_u), size = .2,
                    data = likedat$hhs |> filter(agegr %in% young))

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

savePNG(here('fig/combined_preg_general'), 7, 6)
resize(7, 7)

#' ## Prevalence in non-pregnant women over time, by age group 15-19, 20-24
# how to do this?

#' ## Prevalence in sexually active/non-sexually active over time, age 15-19, 20-24
#'
mwi30 <- mwi_edb$mod$data[1:16, , , ]
mwi30a <- mwi30 - mwi_edb$mod$vpop

mwi30a %<>%
    as.data.table() %>%
    rename_with(~ c("age", "sex", "hiv", "year", "size")) %>%
    mutate(sexual='active')

mwi30n <- mwi_edb$mod$vpop %>%
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
    theme(
        axis.title.y.right = element_text(angle = 0, vjust = 1.05, hjust = 1, margin = margin(l = -20)),
        axis.text.x = element_text(size = 7, angle = 0), 
        axis.text.y = element_text(size = 7), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'grey95'), 
        strip.background = element_rect(fill = '#FCF4DC', color =NA),
        strip.text = element_text(face = 'bold'),
        legend.position = c(1, 1.13)
    ) +
    guides(color = guide_legend(direction = 'horizontal')) +
    coord_cartesian(expand = T) +
    scale_color_manual(values = okabe)

savePNG(here('fig/prev_agr_sex_status'), 7, 4)

#' ## HIV+ : HIV- fertility rate/ratio over time
#'
agr_i <- c(15, 17, seq(20, 50, 5))

epp_frr <- FRR(mwi_epp)
edb_frr <- FRR(mwi_edb)

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

savePNG(here("fig/FRR_MWI.png"), 7, 4)

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

savePNG(here('fig/FR_by_mod_agr_status'), 7, 6)

# ART Coverage
edb_cov <- Cov(edb_mwi)
bind_rows(edb_cov$cov_a, edb_cov$cov_v, .id = "pop") %>%
    mutate(pop = char("active", "virgin")[as.numeric(pop)]) %>%
    ggplot(aes(year, cov, color = factor(stage), linetype = pop)) +
    geom_line() +
    facet_grid(vars(sex), vars(agr)) + theme_light() +
    labs(title='ART coverage by CD4+ stages and sexual activity statuses')
ggsave(here('fig/ART_COV_ZMB.png'), width = 7, height = 5)


#' ## HIV incidence rate ratio parameter estimates between the two models

#+ IRR, fig.cap = 'age IRR between the models', include=TRUE, eval=FALSE
bind_rows(
    mwi_epp$fp$incrr_age[, , 1] %>%
        as_tibble() %>%
        rename_with(~ c("male", "female")) %>%
        mutate(model = "EPP", age  = 15:80),
    mwi_edb$fp$incrr_age[, , 1] %>%
        as_tibble() %>%
        rename_with(~ c("male", "female")) %>%
        mutate(model = "EDB", age = 15:80)
) %>%
    pivot_longer(-c("model", "age")) %>%
    filter(age < 50)  %>% 
        ggplot(aes(age, value, color = model)) +
        facet_wrap(~name) +
        geom_line() +
        ktheme +
        labs(
            title = 'Estimate incidence rate ratio by age and sex',
            y = 'Incidence rate ratio'
        ) +
        theme(
            legend.position = c(.94, 0.1)
        )

savePNG(here('fig/IRR_age'), 7, 4)
resize(7, 4)

#' ## % HIV+ 15-19 and 20-24 who are long-term survivors