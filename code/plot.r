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

# eppasm::popEPP$debug("infect_spec")
# run_a_model(res$MWI, 'eppdb_MWI')
# devtools::load_all('~/GitHub/eppasm')

#' # Outputs:
#'
#' ## Prevalence in pregnant and non-pregnant women over time, age 15-19, 20-24
#'
mwi_epp <- run_a_model(res$MWI, "epp_MWI")$mod$pregprev_agr %>%
    as_tibble() %>%
    rename_with(~ paste0("y_", 1:52)) %>%
    mutate(agr = ago) %>%
    pivot_longer(-agr, names_sep = "_", names_to = c(NA, "year")) %>%
    mutate(year = as.numeric(year)) 

mwi_edb <- run_a_model(res$MWI, "eppdb_MWI")$mod$pregprev_agr %>%
    as_tibble() %>%
    rename_with(~ paste0("y_", 1:52)) %>%
    mutate(agr = agdb) %>%
    pivot_longer(-agr, names_sep = "_", names_to = c(NA, "year")) %>%
    mutate(year = as.numeric(year)) 

#+ preg_prev_mwi, fig.cap = 'Prevalence in pregnant women', include = TRUE
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