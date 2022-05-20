devtools::load_all("~/Code/R/ktools/")
library(tidyverse)

dir.create(here('fit'))
saveRDS(res, here('fit/res.rds'))
res <- readRDS(here('fit/res.rds'))

# Default FRR
ago <- c('15-16', '17-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49')
agdb <- c(15:30, "31-34", "35-39", "40-44", "45-49")

res[["MWI"]]$fits$epp_MWI$fp$frr_cd4[, , ] %>%
    as.data.table() %>%
    rename_with(~ char(CD4_stage, age_group, year, FRR)) %>%
    # mutate(age_group=ago[age_group]) %>%
    filter(year %in% c(1, 10, 20, 30, 40, 50)) %>%
    ggplot(aes(age_group, FRR, color = factor(CD4_stage), linetype = factor(year))) +
    geom_line() +
    scale_x_continuous(labels = ago, breaks = 1:8)

# Fixed FRR
res[["MWI"]]$fits$eppdb_MWI$param$frr_cd4[, , ] |>
    as.data.table() %>%
    rename_with(~ char(CD4_stage, age_group, year, FRR)) %>%
    filter(year %in% c(1, 10, 20, 30, 40, 50)) %>%
    ggplot(aes(age_group, FRR, color = factor(CD4_stage), linetype = factor(year))) +
    geom_line() +
    scale_x_continuous(labels = agdb, breaks = 1:20)

# Age group prevalence
ggsave(here("fig/agp_mwi.png"),
   agp(res[["MWI"]], "MW", "MWI"),
    width = 7, height = 7
)
ggsave(here("fig/agp_uga.png"),
   agp(res[["UGA"]], "UG", "UGA"),
    width = 7, height = 7
)
ggsave(here("fig/agp_zmb.png"),
   agp(res[["ZMB"]], "ZM", "ZMB"),
    width = 7, height = 7
)

ggsave(here("fig/agp_rwa.png"),
   agp(res[["RWA"]], "RW", "RWA"),
    width = 7, height = 7
)
agp <- function(x, iso2, iso3) {
    ddd <- x$data$prev_agesex_nat |>
        filter(ISO == iso2) |>
        mutate(year = as.double(year))
    x$fits[[paste0('epp_', iso3)]]$mod %>% calc_prev_agegr(seq(15, 80, 5)) %>%
        bind_rows(
            x$fits[[paste0('eppdb_', iso3)]]$mod %>% calc_prev_agegr(seq(15, 80, 5)),
            .id = "ver"
        ) %>%
        right_join(ddd, ktools::char(year, sex, agegr)) %>%
        mutate(
            ver = char(epp, eppdb)[as.numeric(ver)],
            agegr = if_else(agegr == "5-9", "05-09", agegr)
        ) %>%
        drop_na(ver) -> merged

    p <- merged |>
        ggplot() +
        geom_line(aes(agegr, prev.x, color = ver, group = ver), size = .8) +
        geom_pointrange(
            aes(agegr, prev.y, ymin = ci_l, ymax = ci_u),
            color = "grey50", shape = 1
        ) +
        facet_wrap(char(year, sex), ncol = 2, scales = "free_y") +
        labs(
            title = "Fitted age-specific",
            x = "Age group", y = "Prevalence", color = "model"
        ) +
        scale_y_continuous(labels = scales::percent) +
            theme(
                legend.position = "bottom",
                axis.text.x = element_text(angle = 90),
                axis.title.x = element_text(vjust = -5)
            )
    p
}

# Rt cross countries
rt <- purrr::map_dfr(iso3_list, \(x) {
    tibble(
        eppdb = res[[x]]$fits[[paste0('eppdb_', x)]]$param$rvec,
        epp = res[[x]]$fits[[paste0('epp_', x)]]$param$rvec,
        iso = x
    ) |> mutate(
        year = 1:n()
    )
})

rt %>%
    pivot_longer(-c(iso, year)) %>%
    ggplot(aes(year, log(value), color = iso, linetype = name)) +
    geom_line() +
    facet_wrap(~iso, scale = "free_y") +
    scale_x_continuous(labels = \(x) 1969 + (x / 10)) +
    labs(x = "Time", y = "Transmission over time")
ggsave(here('fig/rt_6.png'), width = 7, height = 5)

# FRR cross countries
frr <- purrr::map_dfr(iso3_list, \(x) {
    tibble(
        frr = res[[x]]$fits[[paste0("eppdb_", x)]]$param$frr_cd4[1, , 1],
        age_group = agdb,
        type = 'eppdb',
    ) %>% 
    bind_rows(
        tibble(
            frr = res[[x]]$fits[[paste0("epp_", x)]]$param$frr_cd4[1, , 1],
            age_group = ago,
            type = "epp"
        )
    ) %>% mutate(iso = x)
})

frr %>%
    ggplot(aes(age_group, frr, group = type, color = type)) +
    facet_wrap(~iso, scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 60, size = 7)) +
    geom_line() -> p

ggsave(here('fig/frr_6.png'), p, width = 9, height = 5)

# log FRR adjust
purrr::map_dfr(iso3_list, \(x) {
    tibble(
        eppdb = res[[x]]$fits[[paste0("eppdb_", x)]]$param$log_frr,
        epp = res[[x]]$fits[[paste0("epp_", x)]]$param$log_frr,
    ) %>% mutate(iso = x)
}) %>%
    pivot_longer(-iso) %>%
        ggplot(aes(iso, exp(value), fill = name)) +
        geom_col(position = "dodge") +
        labs(x = '', title = "Estimated FRR adjustment")

ggsave(here('fig/log_FRR_6.png'))

# Pregnant prevalence
tibble(
    epp = res[["RWA"]]$fits$epp_RWA$mod$pregprev,
    eppdb = res[["RWA"]]$fits$eppdb_RWA$mod$pregprev,
    year = 1969 + 1:52
) %>%
    filter(year < 2018) |>
        pivot_longer(-year) |>
        ggplot(aes(year, value, color = name)) +
        geom_line() +
        labs(title = "Prevalence in pregnant women", y = "Prevalence") +
        scale_y_continuous(labels = scales::percent)
ggsave(here('fig/preg_MWI.png'), width = 6, height = 4)

# Incidence 
tibble(
    epp = res[["RWA"]]$fits$epp_RWA$mod$incid15to49,
    eppdb = res[["RWA"]]$fits$eppdb_RWA$mod$incid15to49,
    year = 1969 + 1:52
) %>%
    filter(year < 2018) |>
        pivot_longer(-year) |>
        ggplot(aes(year, value, color = name)) +
        geom_line() +
        labs(title = "Incidence rate in 15-49", y = "Incidence rate")
ggsave(here('fig/inc_RWA.png'), width = 6, height = 4)
    
