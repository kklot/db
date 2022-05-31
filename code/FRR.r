FRR <- function(r, on_total=T) {
    mod <- r$mod
    fp <- r$fp
    pa <- r$pa

    pop <- mod$data

    births_a <- apply(pop[16:50, 2, , , drop = FALSE], c(1, 4), sum) * fp$asfr

    ## FRR by HIV age groups (coarse age groups)
    ## Note: calculation restricted to only sexually active population
    ha <- fp$ss$hAG - 1
    frr_ha <- (colSums(fp$frr_cd4 * mod$hivpop[, 1:ha, 2, ]) + 
    colSums(fp$frr_art * mod$artpop[, , 1:ha, 2, ], , 2)) /
        (colSums(mod$hivpop[, 1:ha, 2, ]) + colSums(mod$artpop[, , 1:ha, 2, ], , 2))
        
    ## Expand to FRR by single-year age
    frr_a <- frr_ha[fp$ss$ag.idx[1:35], ]

    ## Calculate prevalene among pregnant women by single-year
    ## Note: restricted to sexually active women only
    apop <- pop
    if (fp$ss$MODEL == 2) {
        apop[1:16, , , ] <- apop[1:16, , , ] - mod$vpop
    }

    hivp_apop <- apop[1:35, 2, 2, ]
    hivn_apop <- apop[1:35, 2, 1, ]

    pregprev_a <- 1 - hivn_apop / (hivn_apop + frr_a * hivp_apop)

    ## Births to HIV+ women and HIV- women
    births_hivp <- pregprev_a * births_a
    births_hivn <- births_a - births_hivp

    ## Fertility rate by 5-year age group and HIV status
    ## (births / total population)

    idx <- rep(paste0(3:9 * 5, "-", 3:9 * 5 + 4), each = 5)

    births_hivp_5yr <- apply(births_hivp, 2, tapply, idx, sum)
    births_hivn_5yr <- apply(births_hivn, 2, tapply, idx, sum)

    if (fp$ss$MODEL == 2 & on_total) {
        apop <- pop
    }
    pop_hivp_5yr <- apply(apop[1:35, 2, 2, ], 2, tapply, idx, sum)
    pop_hivn_5yr <- apply(apop[1:35, 2, 1, ], 2, tapply, idx, sum)

    asfr_hivp <- births_hivp_5yr / pop_hivp_5yr
    asfr_hivn <- births_hivn_5yr / pop_hivn_5yr

    frr_5yr <- asfr_hivp / asfr_hivn
    names(dimnames(frr_5yr))[1] <- "agegr"

    list(asfr_n = asfr_hivn, asfr_p = asfr_hivp, frr = frr_5yr)
}