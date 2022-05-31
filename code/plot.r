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
#' Adding output for pregnant women prevalence by age-group as currently only
#' outputting the 15-49 pregnant prevalence.
#'
#' $$\rho_{a, pregnant} = \dfrac{\sum_{stage} FRR_{age, state} H_{female, age, state} +
#' FRR_{age, state, duration} H_{female, age, stage, duration}}{N_{female, age, H^-} +\sum_{stage} FRR_{age, state} H_{female, age, state} +
#' FRR_{age, state, duration} H_{female, age, stage, duration}}
#' $$
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

frrX <- function(x) {
    # AFSR by age
    x$fp$asfr %>%
        as.data.table(T) %>%
        rename(age = rn) %>%
        pivot_longer(-age, names_to = "year", values_to = "asfr") %>%
        mutate(age = as.numeric(age), year = as.numeric(year)) %>%
        allot(asfr_long)
    
    # asfr_long %>%
    #     filter(year %in% sample(1970:2017, 6)) %>%
    #     ggplot(aes(age, asfr, color = factor(year))) +
    #     geom_line()

    # female population size
    x$mod$data %>%
        as.data.table() %>%
        rename_with(~ char(age, sex, hiv, year, n)) %>%
        mutate(
            age = age + 14, year = year + 1969,
            hiv = char(neg, pos)[hiv],
            sex = char(male, female)[sex],
            year = as.numeric(year)
        ) %>%
        allot(pop_long)
    
    pop_long_active <- pop_long

    if (x$fp$ss$MODEL == 2) {
        x$mod$vpop %>%
            as.data.table() %>%
            rename_with(~ char(age, sex, hiv, year, n)) %>%
            mutate(
                age = age + 14, year = year + 1969,
                hiv = char(neg, pos)[hiv],
                sex = char(male, female)[sex],
                year = as.numeric(year)
            ) %>%
            allot(vpop_long)
        
        # vpop_long %>%
        #     filter(year %in% sample(1970:2017, 6)) %>%
        #     ggplot(aes(age, n, color = hiv)) +
        #     facet_wrap(~year) + geom_line()
        
        pop_long_active <- pop_long %>%
            full_join(vpop_long, char(age, sex, hiv, year)) %>%
            mutate(n = if_else(age <= 30, n.x - n.y, n.x))
        
        # pop_long %>%
        #     filter(year %in% sample(1970:2017, 6)) %>%
        #     ggplot(aes(age, n, color = hiv)) +
        #     facet_wrap(~year) + geom_line()
    } 

    # female size - keep full pop so no need to adjust ASFR
    pop_long %>%
        filter(sex == "female") %>%
        select(age, hiv, year, n) %>%
        pivot_wider(names_from = hiv, values_from = n) |>
        mutate(n = neg + pos)  %>%
        group_by(age) %>%
        mutate(n_roll = zoo::rollmean(n, 2, , T, "right")) %>%
        filter(
            year %in% 1971:2017,
            age %in% 15:49
        ) %>%
        allot(female_size)
           
    # female_size %>%
    #     # filter(age %in% sample(15:49, 6)) %>%
    #     ggplot(aes(year, n, color = factor(age))) +
    #     geom_line() 

    # number of births
    female_size %>%
        ungroup() %>%
        left_join(asfr_long, c("age", "year")) %>%
        mutate(
            birth = n_roll * asfr,
            agr = findInterval2(age, agr_i)
        ) %>%
        group_by(year, agr) %>%
        summarise(birth = sum(birth)) %>%
        allot(birth_agr) # this is kept the same between two models

    # birth_agr %>%
    #     ggplot(aes(year, birth, color = agr)) +
    #     geom_line()

    ## preg prevalence

    # hiv - this part nees to take into account virgin
    pop_long_active %>%
        filter(sex == "female", age %in% 15:49, hiv == "neg") %>%
        group_by(age) %>%
        mutate(n_roll = zoo::rollmean(n, 2, , T, "right")) %>%
        filter(year %in% 1971:2017) %>%
        ungroup() %>%
        mutate(agr = findInterval2(age, agr_i)) %>%
        group_by(agr, year) %>%
        summarise(n_roll = sum(n_roll)) %>%
        allot(hivn)
    hivn

    # hivn %>%
    #     # filter(age %in% sample(15:49, 6)) %>%
    #     ggplot(aes(year, n_roll, color = factor(agr))) +
    #     geom_line()

    # hivp - this part nees to take into account virgin but it tracked separatly
    # expect_equal(
    #     x$mod$data[, , 2, ] %>% sum(),
    #     (x$mod$hivpop[] %>% sum()) + (x$mod$artpop[] %>% sum()) +
    #         (x$mod$vpopart %>% sum()) +
    #         (x$mod$vpophiv %>% sum())
    # )
    x$mod$hivpop %>%
        as.data.table() %>%
        rename_with(~ char(stage, agr, sex, year, n)) %>%
        filter(agr < 9, sex == 2) %>%
        mutate(agr = ago[agr], year = 1969 + year) %>%
        group_by(stage, agr) %>%
        mutate(n_roll = zoo::rollmean(n, 2, , T, "right")) %>%
        filter(year %in% 1971:2017) %>%
        ungroup() %>%
        select(stage, agr, year, n_roll) %>%
        allot(hivp)
    hivp

    # hivp %>%
    #     ggplot(aes(year, n_roll, color = agr, linetype = factor(stage))) +
    #     facet_wrap(~stage, scales= 'free') +
    #     geom_line()

    x$fp$frr_cd4 %>%
        as.data.table() %>%
        rename_with(~ char(stage, agr, year, frr)) %>%
        mutate(year = 1969 + year, agr = ago[agr]) %>%
        right_join(hivp) %>%
        mutate(birthx = frr * n_roll) %>%
        group_by(agr, year) %>%
        summarise(birth_hiv = sum(birthx)) %>%
        allot(hivp_birthrr)
    hivp_birthrr

    # hivp_birthrr %>%
    #     ggplot(aes(year, birth_hiv, color = agr)) + geom_line()

    # art
    x$mod$artpop %>%
        as.data.table() %>%
        rename_with(~ char(dur, stage, agr, sex, year, n)) %>%
        filter(agr < 9, sex == 2) %>%
        mutate(agr = ago[agr], year = 1969 + year) %>%
        group_by(dur, stage, agr) %>%
        mutate(n_roll = zoo::rollmean(n, 2, , T, "right")) %>%
        filter(year %in% 1971:2017) %>%
        ungroup() %>%
        select(dur, stage, agr, year, n_roll) %>%
        allot(art)
    art

    x$fp$frr_art %>%
        as.data.table() %>%
        rename_with(~ char(dur, stage, agr, year, frr)) %>%
        mutate(year = 1969 + year, agr = ago[agr]) %>%
        right_join(art, char(dur, stage, agr, year)) %>%
        mutate(birthx = frr * n_roll) %>%
        group_by(agr, year) %>%
        summarise(birth_art = sum(birthx)) %>%
        allot(art_birthrr)
    art_birthrr

    # art_birthrr %>%
    #     ggplot(aes(year, birth_art, color = agr)) + geom_line()    

    # prev
    hivn %>%
        rename(hiv_neg = n_roll) %>%
        left_join(hivp_birthrr, char("agr", "year")) %>%
        left_join(art_birthrr, char("agr", "year")) %>%
        left_join(birth_agr, char("agr", "year")) %>%
        mutate(
            p_factor = birth_hiv + birth_art,
            prev = p_factor / (p_factor + hiv_neg),
            B_pos = birth * prev,
            B_neg = birth - B_pos
        ) %>%
        allot(birth_status)
    
    birth_status

    # birth_status %>%
    # ggplot(aes(year, prev, color = agr)) + geom_line()

    # birth_status %>%
    #     ggplot(aes(year, B_pos, color = agr)) + geom_line()

    pop_long %>%
        left_join(pop_long_active %>%
            select(n_active = n, age, sex, hiv, year)) %>%
        group_by(age, sex, hiv) %>%
        mutate(
            n_roll = zoo::rollmean(n, 2, NA, TRUE, "right"),
            n_active_roll = zoo::rollmean(n_active, 2, NA, TRUE, "right")
        ) %>%
        ungroup() %>%
        filter(sex == "female", year %in% 1971:2017) %>%
        mutate(agr = findInterval2(age, agr_i)) %>%
        group_by(agr, hiv, year) %>%
        summarise(n = sum(n_roll), n_active = sum(n_active_roll)) %>%
        pivot_wider(names_from = hiv, values_from = c(n, n_active)) %>%
        left_join(birth_status, char(agr, year)) %>%
        mutate(
            fr_pos = B_pos / n_pos,
            fr_neg = B_neg / n_neg,
            frr = fr_pos / fr_neg,
            fr_pos_active = B_pos / n_active_pos,
            fr_neg_active = B_neg / n_active_neg,
            frr_active = fr_pos_active / fr_neg_active
        ) %>%
        allot(frr_epp)
    list(
        female_size = female_size,
        birth_agr = birth_agr,
        hiv_neg = hivn,
        hiv_pos_rr = hivp_birthrr,
        art_rr = art_birthrr,
        birth_status = birth_status,
        frr_data = frr_epp
    )
    # frr_epp
}

epp_mwi <- run_a_model(res$ZMB, "epp_ZMB")
frr_epp <- frrX(epp_mwi)

edb_mwi <- run_a_model(res$ZMB, "eppdb_ZMB")
frr_edb <- frrX(edb_mwi)

frr_edb$frr_data %>%
    # select(agr, year, fr_pos, fr_neg, frr) %>%
    select(agr, year, fr_pos_active, fr_neg_active, frr_active) %>%
    left_join(frr_epp$frr_data %>%
        select(agr, year, fr_pos, fr_neg, frr), char(agr, year)) %>%
    allot(joint_frr)

#+ FRR_2_models, fig.cap = 'Fertility rate ratio between model with and without sexual debut', include=TRUE
joint_frr %>%
    filter(
        # as.numeric(substr(agr, 1, 2)) <= 20, 
        year >= 1976
    ) %>%
    ggplot(aes(year, frr, color = "EPP")) +
        geom_line() +
        geom_line(aes(y = frr_active, color = "EPPDB")) +
        facet_wrap(~agr, scales = "free_y") +
        labs(title = "Fertility rate ratio", y = "FRR")

ggsave(here('fig/FRR_MW.png'), width = 7, height = 4)
ggsave(here('fig/FRR_ZMB.png'), width = 7, height = 4)

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