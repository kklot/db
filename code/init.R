rm(list = ls())
here::i_am("code/init.R")
library(here)
library(data.table)
library(dplyr)
library(tidyr)

# remotes::install_github('jeffeaton/anclik')
# remotes::install_github('jeffeaton/epp')

# install.packages(c("abind", "anclik", "binom", "fastmatch", "mvtnorm",
# "RColorBrewer", "reshape2", "BH", "RcppEigen"))
devtools::load_all("~/Github/eppasm", recompile = TRUE)
# install.packages(c('txtplot', 'countrycode', 'snakecase', 'ggplot2'))
devtools::load_all("~/Code/R/ktools")

do_cc <- function(iso2, iso3) {
    fit <- eppasm$new()

    readRDS(here("data/par23.rds")) %>%
        filter(ISO_A3 == iso3) %>%
        #   to fill backward and forward need to separate male feavalmale
        pivot_wider(char(yob, ISO_A3),
            names_from = sex, values_from = char(lambda, skew, shape)
        ) %>%
        #   add needed year, 30 year old in 1970 and 15 2021
        bind_rows(tibble(yob = c(
            (1970 - 30):(1950 - 1),
            (2005 + 1):(2021 - 15)
        ))) %>%
        arrange(yob) %>%
        tidyr::fill(-yob, .direction = "updown") %>%
        #   done filling, generate year cohorts
        mutate(dup = length(1970:2021)) %>%
        uncount(dup, .id = "year") %>%
        mutate(year = year + 1969, age = year - yob) %>%
        filter(age %in% 15:30) %>%
        pivot_longer(
            -char(yob, ISO_A3, year, age),
            names_sep = "_",
            names_to = char(.value, sex)
        ) %>%
        mutate(db = dskewlogis(age, lambda, shape, skew) /
            (1 - pskewlogis(age, lambda, shape, skew))) %>%
        # change age 15 to cumulative as the model does not have 15- age
        mutate(db = if_else(age == 15, pskewlogis(15, lambda, shape, skew), db)) %>%
        mutate(db = if_else(year == 1970, pskewlogis(age, lambda, shape, skew), db)) %>%
        select(year, age, sex, db) %>%
        mutate(sex = if_else(sex == "female", 2, 1)) %>%
        #   sort like epp, year sex age
        arrange(year, sex, age) %>%
        pull(db) %>%
        array(c(16, 2, length(1970:2021))) %>%
        allot(est_db)
    fit$data$est_db[[iso2]] <- est_db
    # No est_senesence
    fit$data$est_senesence[[iso2]][] <- 1
    # no condom use
    fit$data$est_condom[[iso2]] <- array(0, c(66, 2, 52))
    # bw <- prepare_spec_fit("~/Downloads/Malawi_2020_v09.pjnz", proj.end=2022.5) #
    # fp <- attr(bw[[ "Northern Region" ]], "specfp")
    # saveRDS(fp$cd4_prog, here("data/newcd4progress.rds"))
    new_cd4 <- readRDS(here("data/newcd4progress.rds"))
    attr(fit$data$inputs_nat[[iso2]], "specfp")$cd4_prog <- new_cd4

    female_susceptibility <- rep(1, length(15:80))
    control_optim <- list(hessian = FALSE, control = list(maxit = 100, reltol = 1e-4))

    #' # Base model
    fit$fit(
        paste0("EPP", "_", iso3),
        iso2,
        B0 = 1e3,
        eppmod = "rlogistic", control_optim = control_optim
    )

    #' Debut only model
    fit$fit(
        paste0("EDB", "_", iso3),
        iso2,
        B0 = 1e3,
        VERSION = "K",
        eppmod = "rlogistic",
        doParallel = FALSE,
        control_optim = control_optim,
        with_debut = TRUE,
        with_mixing = FALSE
    )
    fit
}

tst <- system.file('extdata', 'inputs_nat.rds', package = 'eppasm')

iso <- tibble(
    iso2 = tst %>% readRDS() %>% names(),
    iso3 = countrycode::countrycode(iso2, "iso2c", "iso3c")
)  %>% 
drop_na()

res <- parallel::mclapply(
    1:nrow(iso),
    function(i) {
        tryCatch(
            do_cc(iso$iso2[i], iso$iso3[i]), 
            error = function(e) e
        )
    },
    mc.cores = min(parallel::detectCores(), nrow(iso))
)

dir.create(here('fit'))
saveRDS(res, here('fit/res.rds'))