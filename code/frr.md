# Methods

## Sexual debut transition

[copy from uganda draft]

## Force of infection

The transmission probability for each sexual act among homogeneous mixing of
sexually active population is defined as

$$\rho_{t} = \dfrac{N_{t}^{+,*} - A_{t}^{*} (1 - \omega)}{N_{t}^{*}}$$

where $\omega$ is relative reduction among patient who is treated with ART,
$A_{t}^{*}$ is the total HIV positive population on ART treatment, and $N_{t}^{*}$
is the total sexually active population.  

The overall incidence rate over time is then 

$$i_t = r_{t} \delta(\iota, t_0) \rho_t,$$

with $r_{t}$ is a logistic curve, $\delta$ is the Kronecker delta
function of the starting time of the epidemic and the initial prevalence $\iota$.

The incidence rate between male and female, $IRR_{sex}$, is used to calculate
the incidence rate by sex as

$$s_{adj} = \dfrac{N_{15-49}^{+,*}}{N_{15-49, male}^{-,*} +
N_{15-49, female}^{-,*} IRR_{sex}}$$

$$i_{male, t} = i_t s_{adj}$$ 

$$i_{female, t} = i_t s_{adj} IRR_{sex}$$ 

such that the total incidence calculated based on all population or based on
sex-specific population is the same.

Similarly, incidence rate ratio between each single-age against the 20–24
year-old is used to calculate the incidence rate by each single-age as

$$a_{age, male, t} = \rho_{male} \dfrac{N_{15-49, male}^{-,*}}{\sum_{age=15}^{49}N_{age,male}^{-,*}IRR_{age,male}}$$

and similarly for female.  

$$i_{age, sex, t} = IRR_{age, sex} a_{age, sex, t}$$ 

The new number of infections in each age and sex 

$$i_{age, sex} = N_{age, sex}^{-,*} i_{age, sex, t}$$

## HIV+:HIV- fertility rate ratio over time

The fertility rate (FR) for age-group $a = \{15,...,49\}$ and year $t$ is calculated as

$$FR_{a,t}^{+,*} = \dfrac{B_{a,t}^{+}}{F_{a,t}^{+,*}}, \qquad
FR_{a,t}^{-,*} = \dfrac{B_{a,t}^{-}}{F_{a,t}^{-,*}},$$

where $B$ the number of births with HIV positive or negative denoted with plus
(+) or minus (-) sign. $F$ is the number of women sexually active denoted with
star (*) sign. In the model without sexual debut, this number is taken as the
total HIV positive women.

The number of births to HIV positive and sexually active women $B_{a, t}^{+}$ is
estimated as

$$B_{a,t}^{+} = \rho_{a, t}B_{a,t}, \qquad B_{a,t}^{-} = B_{a,t} - B_{a, t}^{+}$$

with $\rho$ the prevalence among pregnant women and $B_{a,t}$ is the estimate
number of births in general population regardless of HIV statuses. In
particular, 

$$B_{a,t} = \frac{1}{2} \left(N_{a, t} + N_{a, t-1}\right) ASFR_{a}$$

with $N$ is the population size including sexually active and non-active women
and ASFR is the age-specific fertility rate input estimated based on general
population. Given that the number of births is fixed and in the model with
sexual debut the number of sexually active in adolescents and young adults is
reduced, the ASFR needs to be adjusted to account for the differences in the
number of births as 

$$ASFR_{a,t}^{*} = \frac{1}{\pi_{a,t}}ASFR_a,\qquad
\pi_{a, t} = \dfrac{N_{a,t}^{*}+N_{a,t-1}^{*}}{N_{a,t}+N_{a,t-1}},$$ 

with $\pi_{a,t}$ is the proportion of sexually active in age $a$ at time $t$.
The number of births is then 

$$B^{*}_{a,t} = \frac{1}{2} \left(N^{*}_{a, t} + N^{*}_{a, t-1}\right) ASFR^{*}_{a}
= B_{a,t}.$$

The prevalence $\rho$ is calculated as

$$b_{a,t} = \sum_{s=1}^{7} FRR_{a,s}^{CD4}F_{a,s,t}^{+,*} +
\sum_{s=1}^{7}\sum_{d=1}^{3}FRR_{a,s,d,t}^{ART} F_{a,s,d,t}^{+,ART},$$
$$\rho_{a,t} = \dfrac{b_{a,t}}{F_{a,t}^{-} + b_{a,t}} =
\dfrac{B_{a,t}^{+}}{B_{a,t}^{+} + B_{a,t}^{-}}.$$

Estimation method of the fertility rate ratio between the different stages of
HIV/AIDS classified by seven stages of CD4+ count ($s$) and three treatment
duration $d$, $FRR_{a,s}^{CD4}$ and $FRR_{a,s}^{ART}$ respectively,
are documented in [cite]. 

## Model fitting

Estimated sexual debut parameters for each countries were obtained from [cite].
Free parameters include the four parameters 'rlogistic' model [cite] to reflect
time trend of the transmission probability, initial prevalence at the start of the
epidemic, age's incidence rate ratio, incidence rate ratio between sex, and five
parameters to incorporate ANC data in the model likelihood [cite] and
non-sampling error in the ANC data [cite].

Assuming that the sexual debut pattern could explain the difference of fertility
rate in adolescent and young adults, in the model with sexual debut transition,
the $FRR_{15-19,s}^{CD4}$ was set to the same level as the
$FRR_{20-24,s}^{CD4}$.

# Results

## Estimated FR among HIV+ and HIV- 

### Malawi

![FR_MW](../fig/FR_MW.png)

![FRR_MW](../fig/FRR_MW.png)

### Zambia

![FR_ZMB](../fig/FR_ZMB.png)

![FRR_ZMB](../fig/FRR_ZMB.png)

both and sexually active

- number HIV+ in each age-group
- prev pregnant age-group
- ASFR by age-group