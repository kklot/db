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

#' # Outputs:
#'
# Data and simulated
likedat <- res$MWI$fits$epp_MWI$likdat

mwi_epp <- run_a_model(res$MWI, "epp_MWI")
mwi_edb <- run_a_model(res$MWI, "eppdb_MWI")

#' ## Prevalence 15-49 by age-group
sim_prevagr <- calc_prev_agegr(mwi_epp$mod, seq(15, 55, 5)) %>%
    bind_rows(
        calc_prev_agegr(mwi_edb$mod, seq(15, 55, 5)),
        .id = "model"
    ) %>%
    mutate(model = char(EPP, EPPDB)[as.numeric(model)])

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

options(ggplot2.discrete.color = okabe)

tibble(
    prev_EPP = pregprev1549_agr,
    prev_EDB = pregprev1549db_agr,
    agegr = rep(1:7, times = 48),
    year = rep(1969 + 1:48, each = 7)
) |>
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

mwi_epp %>%
    filter(as.numeric(substr(agr, 1, 2)) <= 25) %>%
    bind_rows(
        mwi_edb %>% filter(as.numeric(substr(agr, 1, 2)) <= 24),
        .id = "model"
    ) %>%
    mutate(model = c("EPP", "EPPDB")[as.numeric(model)]) %>%
    ggplot(aes(year, value, color = agr, linetype = model)) +
        geom_line() +
        coord_cartesian(xlim = c(0, 45)) +
        scale_x_continuous(labels = \(x) x + 1969) +
        scale_y_continuous(labels = scales::percent) +
        labs(title = "Prevalence in pregnant women by age-group", y='')

#' ## Prevalence in sexually active/non-sexually active over time, age 15-19, 20-24
#'
# TODO need to check this
# mwi_edb$vpop[,,2,] %>% sum
# (mwi_edb$vpophiv %>% sum) + (mwi_edb$vpopart %>% sum)
mwi_edb <- run_a_model(res$MWI, "eppdb_MWI")$mod
mwi30 <- mwi_edb$data[1:16, , , ]
mwi30a <- mwi30 - mwi_edb$vpop

mwi30a %<>%
    as.data.table() %>%
    rename_with(~ c("age", "sex", "hiv", "year", "size")) %>%
    mutate(sexual='active')

mwi30n <- mwi_edb$vpop %>%
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
        summarise(
            prev = sum(`2`) / sum(`1` + `2`)
        ) %>%
        allot(mwi30)

#+ prev_active_mwi, fig.cap = 'Prevalence in sexually active and non-active', include=TRUE
mwi30 %>%
    ggplot(aes(year, prev, color = sexual)) +
    facet_grid(vars(sex), vars(agr)) +
    geom_line() +
    scale_y_continuous(labels = scales::percent, trans = "log") +
    scale_x_continuous(labels = \(x) x + 1969) +
    labs(
        title = "Prevalence in sexually active/nonactive by age-group MWI",
        y = "log(Prevalence)"
    )

#' ## HIV+ : HIV- fertility rate/ratio over time
#'
agr_i <- c(15, 17, seq(20, 50, 5))

epp_mwi <- run_a_model(res$ZMB, "epp_ZMB")
edb_mwi <- run_a_model(res$ZMB, "eppdb_ZMB")

epp_frr <- FRR(epp_mwi)
edb_frr <- FRR(edb_mwi)
edb_frr_active <- FRR(edb_mwi, F)

bind_rows(
    epp_frr$frr %>% as.data.table(1) %>% pivot_longer(-rn),
    edb_frr$frr %>% as.data.table(1) %>% pivot_longer(-rn),
    edb_frr_active$frr %>% as.data.table(1) %>% pivot_longer(-rn),
    .id = "model"
) %>%
    type.convert() %>%
        filter(rn %in% c("15-19", "20-24", "25-29")) %>%
        mutate(model = c("epp", "eppdb", "eppdb_active")[model]) %>%
        ggplot(aes(name, value, color = rn, linetype = model)) +
        facet_wrap(~rn) +
        geom_line() +
        labs(subtitle = 'eppdb_active: FRR on sexually active', color = 'Age-group')

ggsave(here("fig/FRR_ZMB.png"), width = 7, height = 4)

#+ FR_2_models, fig.cap = 'Fertility rate between model with and without sexual debut', include=TRUE
joint_frr %>%
    filter(
        # as.numeric(substr(agr, 1, 2)) <= 20,
        year >= 1972
    ) %>%
    ggplot(aes(year, fr_pos, linetype = "EPP", color = "HIV+")) +
        geom_line() +
        geom_line(aes(y = fr_pos_active, linetype = "EPPDB", color = "HIV+")) +
        geom_line(aes(y = fr_neg, linetype = "EPP", color = "HIV-")) +
        geom_line(aes(y = fr_neg_active, linetype = "EPPDB", color = "HIV-")) +
        facet_wrap(~agr, scales = "free_y") +
        coord_cartesian(ylim = c(0, .5)) +
        labs(title = "Fertility rate", y = "FR", color = 'HIV status', linetype = 'Model')
    
ggsave(here('fig/FR_MW.png'), width = 7, height = 4)
ggsave(here('fig/FR_ZMB.png'), width = 7, height = 4)
    

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
    epp_mwi$fp$incrr_age[, , 1] %>%
        as_tibble() %>%
        rename_with(~ c("male", "female")) %>%
        mutate(model = "epp", age  = 15:80),
    eppdb_mwi$fp$incrr_age[, , 1] %>%
        as_tibble() %>%
        rename_with(~ c("male", "female")) %>%
        mutate(model = "eppdb", age = 15:80)
) %>%
    pivot_longer(-c("model", "age")) %>%
        ggplot(aes(age, value, color = model)) +
        facet_wrap(~name) +
        geom_line()

#' ## % HIV+ 15-19 and 20-24 who are long-term survivors