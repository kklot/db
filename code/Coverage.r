Cov <- function(r) {
    mod <- r$mod

    agreg <- c(rep(c('15-19', '20-24', '25-29'), each = 5), rep('30-34', 2), '35-39', '40-44', '45-49')

    a_long <- mod$artpop %>%
        as.data.table() %>%
        rename_with(~ char(dur, stage, agr, sex, year, n)) %>%
        mutate(agr = agreg[agr]) %>%
        group_by(dur, stage, agr, sex, year) %>%
        summarise(n = sum(n)) %>% ungroup()

    a_long_v <- mod$vpopart %>%
        as.data.table() %>%
        rename_with(~ char(dur, stage, agr, sex, year, n)) %>%
        mutate(agr = agreg[agr]) %>%
        group_by(dur, stage, agr, sex, year) %>%
        summarise(n = sum(n)) %>% ungroup()

    h_long <- mod$hivpop %>%
        as.data.table() %>%
        rename_with(~ char(stage, agr, sex, year, n)) %>%
        mutate(agr = agreg[agr]) %>%
        group_by(stage, agr, sex, year) %>%
        summarise(n = sum(n)) %>% ungroup()

    h_long_v <- mod$vpophiv %>%
        as.data.table() %>%
        rename_with(~ char(stage, agr, sex, year, n)) %>%
        mutate(agr = agreg[agr]) %>%
        group_by(stage, agr, sex, year) %>%
        summarise(n = sum(n)) %>% ungroup()

    cov_a <- a_long %>%
        group_by(stage, agr, sex, year) %>%
        summarise(n = sum(n)) %>%
        full_join(h_long, char(stage, agr, sex, year)) %>%
        mutate(cov = n.x/(n.x + n.y))
    cov_v <- a_long_v %>%
        group_by(stage, agr, sex, year) %>%
        summarise(n = sum(n)) %>%
        full_join(h_long_v, char(stage, agr, sex, year)) %>%
        mutate(cov = n.x / (n.x + n.y))
    
    cov_a %<>%
        filter(year > 30) %>%
        filter(agr %in% char("15-19", "20-24", "25-29")) %>%
        mutate(
            year = year + 1969,
            sex = char("male", "female")[sex]
        )
    
    cov_v %<>%
        filter(year > 30) %>%
        filter(agr %in% char("15-19", "20-24", "25-29")) %>%
        mutate(
            year = year + 1969,
            sex = char("male", "female")[sex]
        )

    list(cov_a = cov_a, cov_v = cov_v)
}