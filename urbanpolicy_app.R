##### APPLICATION TO URBAN POLICY DATA #####
library(haven)
library(survey)
library(tidyverse)
library(ggpubr)

# Load data
base.menage <- read_sas("menage.sas7bdat")

# Inclusion probabilities of PSUs
pii <- base.menage %>%
  mutate(pii = 1 / SamplingWeight) %>%
  arrange(desc(pii)) %>%
  distinct(pii, .keep_all = TRUE) %>%
  select(str, SamplingWeight, pii)
mean(pii$pii)

# Approximated first-stage sampling rate
NI <- pii %>%
  summarize(sum_samplingweight = sum(SamplingWeight)) %>%
  mutate(rate = 40 / sum_samplingweight) 
# repondant contains data on all the household - each row is a household
repondant <- subset(base.menage, repondant == 1 & !(is.na(f8) | is.na(g3b) | is.na(h1c) | is.na(i6))) 

# Frequencies of variables f8, g3b, h1c, i6
freq_repondant <- table(repondant$f8, repondant$g3b, repondant$h1c, repondant$i6)
print(freq_repondant)

# Sample of respondents
repondant <- base.menage %>%
  filter(repondant == 1) %>%
  mutate(Ni = (poidsmen_init / SamplingWeight) * 40) %>%
  arrange(id_zus)

table_zus <- repondant %>%
  group_by(id_zus) %>%
  summarize(nechi = n(),
             Ni = mean(Ni)) %>%
  mutate(Ni = round(Ni))

# Repondant2 data
repondant2 <- merge(repondant, table_zus, by = "id_zus") %>%
  select(-Ni.x) %>%
  rename(Ni = Ni.y)
repondant2 <- repondant2 %>%
  mutate(
    pi_i = 1 / SamplingWeight,
    d_i = SamplingWeight,
    pi_ksaci = nechi / Ni,
    d_ksaci = 1 / pi_ksaci,
    pi_k = pi_i * pi_ksaci,
    d_k = 1 / pi_k,
    un = 1,
    rep_a = as.numeric(f8 == 1) + as.numeric(f8 == 2),
    rep_b = as.numeric(f8 == 3),
    rep_c = as.numeric(f8 == 4) + as.numeric(f8 == 5) + as.numeric(f8 == 6),
    rep_d = as.numeric(f8 == 7),
    wit_a = as.numeric(g3b == 1),
    wit_b = as.numeric(g3b == 2),
    wit_c = as.numeric(g3b == 3) + as.numeric(g3b == 4),
    wit_d = as.numeric(g3b == 5) + as.numeric(g3b == 6),
    wor_a = as.numeric(h1c == 1) + as.numeric(h1c == 2),
    wor_b = as.numeric(h1c == 3),
    wor_c = as.numeric(h1c == 4) + as.numeric(h1c == 5),
    mov_a = as.numeric(i6 == 1) + as.numeric(i6 == 2),
    mov_b = as.numeric(i6 == 3),
    mov_c = as.numeric(i6 == 4),
    mov_d = as.numeric(i6 == 5) + as.numeric(i6 == 6)
  )

############################
####### REP CATEGORY #######
############################

###### REP_A ######
varia <- data.frame(var = "rep_a")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(rep_a, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_repa <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (rep_a - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_repa %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_repa %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
rep_a_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
rep_a_result <- rep_a_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### REP_B ######
varia <- data.frame(var = "rep_b")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(rep_b, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_repb <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (rep_b - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_repb %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_repb %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
rep_b_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
rep_b_result <- rep_b_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### REP_C ######
varia <- data.frame(var = "rep_c")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(rep_c, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_repc <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (rep_c - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_repc %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_repc %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
rep_c_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
rep_c_result <- rep_c_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### REP_D ######
varia <- data.frame(var = "rep_d")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(rep_d, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_repd <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (rep_d - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_repd %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_repd %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
rep_d_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
rep_d_result <- rep_d_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

############################
####### WIT CATEGORY #######
############################
repondant2[is.na(repondant2)] <- 0

###### WIT_A ######
varia <- data.frame(var = "wit_a")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(wit_a, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_wita <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (wit_a - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_wita %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_wita %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh, by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
wit_a_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
wit_a_result <- wit_a_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### WIT_B ######
varia <- data.frame(var = "wit_b")
pest <- repondant2 %>%
  drop_na() %>%
  summarise(phat = weighted.mean(wit_b, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_witb <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (wit_b - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_witb %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_witb %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
wit_b_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
wit_b_result <- wit_b_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### WIT_C ######
varia <- data.frame(var = "wit_c")
pest <- repondant2 %>%
  drop_na() %>%
  summarise(phat = weighted.mean(wit_c, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_witc <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (wit_c - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_witc %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_witc %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
wit_c_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
wit_c_result <- wit_c_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### WIT_D ######
varia <- data.frame(var = "wit_d")
pest <- repondant2 %>%
  drop_na() %>%
  summarise(phat = weighted.mean(wit_d, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_witd <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (wit_d - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_witd %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_witd %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
wit_d_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
wit_d_result <- wit_d_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

############################
####### WOR CATEGORY #######
############################

###### WOR_A ######
varia <- data.frame(var = "wor_a")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(wor_a, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_wora <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (wor_a - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_wora %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_wora %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
wor_a_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
wor_a_result <- wor_a_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### WOR_B ######
varia <- data.frame(var = "wor_b")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(wor_b, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_worb <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (wor_b - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_worb %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_worb %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
wor_b_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
wor_b_result <- wor_b_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### WOR_C ######
varia <- data.frame(var = "wor_c")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(wor_c, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_worc <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (wor_c - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_worc %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_worc %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
wor_c_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
wor_c_result <- wor_c_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

############################
####### MOV CATEGORY #######
############################

###### MOV_A ######
varia <- data.frame(var = "mov_a")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(mov_a, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_mova <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (mov_a - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_mova %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_mova %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
mov_a_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
mov_a_result <- mov_a_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### MOV_B ######
varia <- data.frame(var = "mov_b")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(mov_b, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_movb <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (mov_b - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_movb %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_movb %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
mov_b_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
mov_b_result <- mov_b_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### MOV_C ######
varia <- data.frame(var = "mov_c")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(mov_c, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_movc <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (mov_c - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_movc %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_movc %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
mov_c_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
mov_c_result <- mov_c_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

###### MOV_D ######
varia <- data.frame(var = "mov_d")
pest <- repondant2 %>%
  summarise(phat = weighted.mean(mov_d, d_k), Nhat = sum(un * d_k)) %>%
  select(phat, Nhat) %>% 
  mutate(Nhat = as.numeric(Nhat), phat = as.numeric(phat))

repondant2_movd <- repondant2 %>%
  mutate(e_ik = (1 / pest$Nhat) * (mov_d - pest$phat)) %>% 
  arrange(id_zus)

# Second variance component
piv1 <- repondant2_movd %>%
  group_by(id_zus) %>%
  summarise(s2ei = var(e_ik), nechi = mean(nechi), 
            Ni = mean(Ni), pi_i = mean(pi_i)) %>%
  mutate(vhtb = ((Ni^2) / pi_i) * (1 / nechi - 1 / Ni) * s2ei) %>%
  ungroup()

vhtb <- piv1 %>% summarise(vhtb = sum(vhtb))

# First variance component
piv2 <- repondant2_movd %>%
  group_by(id_zus) %>%
  summarise(str = mean(str),
            pi_i = weighted.mean(pi_i, d_ksaci),
            Eihat = sum(e_ik * d_ksaci)) %>% 
  ungroup() %>%
  mutate(Eihat = Eihat / pi_i, 
         w = 1 - pi_i) %>%
  arrange(id_zus)

Reh <- piv2 %>%
  mutate(str = as.numeric(str)) %>%
  group_by(str) %>%
  summarise(Reh = weighted.mean(Eihat, w),
            n = n()) %>% 
  ungroup()  %>%
  arrange(str)

piv2 <- merge(piv2, Reh %>% select(Reh, str), by = "str") %>%
  mutate(Eihat = w * ((Eihat - Reh)^2))

vhaja <- piv2 %>% summarise(vhaja = sum(Eihat))

# Gathering
mov_d_result <- data.frame(var = varia$var, phat = pest$phat, vhaja = vhaja$vhaja, vhtb = vhtb$vhtb)
mov_d_result <- mov_d_result %>%
  mutate(
    binf1 = phat - 1.96 * sqrt(vhaja + vhtb),
    bsup1 = phat + 1.96 * sqrt(vhaja + vhtb),
    binf2 = phat - 1.96 * sqrt(vhaja),
    bsup2 = phat + 1.96 * sqrt(vhaja)
  )

############# COMBINING ALL RESULTS ###############
full_res <- rbind(rep_a_result, rep_b_result, rep_c_result, rep_d_result,
                  wit_a_result, wit_b_result, wit_c_result, wit_d_result,
                  wor_a_result, wor_b_result, wor_c_result,
                  mov_a_result, mov_b_result, mov_c_result, mov_d_result) %>%
  as.data.frame() %>%
  select(-c(vhaja, vhtb)) 
full_res <- sapply(full_res, function(x) if(is.numeric(x)) round(x, 3) else x) 
var_type <- c("rep", "rep", "rep", "rep", "wit", "wit", "wit", "wit", 
              "wor", "wor",  "wor", "mov", "mov", "mov", "mov")
full_res <- cbind(full_res, var_type) %>%
  as.data.frame() %>%
  mutate(phat = as.numeric(phat),
         binf1 = as.numeric(binf1),
         bsup1 = as.numeric(bsup1),
         binf2 = as.numeric(binf2),
         bsup2 = as.numeric(bsup2))

############# CREATE FIGURES ###############

# mov variable figure
mov_df <- full_res %>% 
  filter(var_type == "mov") %>%
  mutate(var_meaning = case_when(var == "mov_a" ~ "Certainly or Probably",
                                 var == "mov_b" ~ "Probably Not",
                                 var == "mov_c" ~ "Certainly Not",
                                 var == "mov_d" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "Certainly Not", "Probably Not", "Certainly or Probably"), 
    ordered = TRUE))

mov_fig <- ggplot(data = mov_df, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1, color = "Hajek"), alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2, color = "Hajek-A"), alpha = 0.3) +
  labs(x = "Point Estimate",
       title = "Intention to Leave the District") + 
  scale_color_manual(name='Variance',
                     breaks=c('Hajek', 'Hajek-A'),
                     values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# rep variable figure
rep_df <- full_res %>% 
  filter(var_type == "rep") %>%
  mutate(var_meaning = case_when(var == "rep_a" ~ "Good",
                                 var == "rep_b" ~ "Fair",
                                 var == "rep_c" ~ "Poor",
                                 var == "rep_d" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "Poor", "Fair", "Good"), ordered = TRUE))

rep_fig <- ggplot(data = rep_df, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1), color = "firebrick", alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2), color = "blue", alpha = 0.3) +
  labs(x = "Point Estimate",
       y = "Variable of Interest",
       title = "Perceived District Reputation") + 
  # scale_color_manual(name='Variance',
  #                    breaks=c('Hajek', 'Hajek-A'),
  #                    values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# wit variable figure
wit_df <- full_res %>% 
  filter(var_type == "wit") %>%
  mutate(var_meaning = case_when(var == "wit_a" ~ "Never",
                                 var == "wit_b" ~ "Rarely",
                                 var == "wit_c" ~ "Sometimes",
                                 var == "wit_d" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "Sometimes", "Rarely", "Never"), ordered = TRUE))

wit_fig <- ggplot(data = wit_df, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1, color = "Hajek"), alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2, color = "Hajek-A"), alpha = 0.3) +
  labs(x = "Point Estimate",
       title = "Witnessed Trafficking") + 
  scale_color_manual(name='Variance',
                     breaks=c('Hajek', 'Hajek-A'),
                     values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# wor variable figure
wor_df <- full_res %>% 
  filter(var_type == "wor") %>%
  mutate(var_meaning = case_when(var == "wor_a" ~ "Yes",
                                 var == "wor_b" ~ "No",
                                 var == "wor_c" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "No", "Yes"), ordered = TRUE))

wor_fig <- ggplot(data = wor_df, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1), color = "firebrick", alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2), color = "blue", alpha = 0.3) +
  labs(x = "Point Estimate",
       y = "Variable of Interest",
       title = "Roadworks in Neighborhood") + 
  # scale_color_manual(name='Variance',
  #                    breaks=c('Hajek', 'Hajek-A'),
  #                    values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# combine plots
myres_fullfig <- ggarrange(rep_fig, wit_fig, wor_fig, mov_fig, 
                           nrow = 2, ncol = 2)
ggsave("myres_fullfig.png")

##################################################
######### COMPARISON TO DATA FROM PAPER ##########
##################################################

var_type <- c("rep", "rep", "rep", "rep", "wit", "wit", "wit", "wit", 
              "wor", "wor",  "wor", "mov", "mov", "mov", "mov")
var <- full_res$var
phat <- full_res$phat # same results as in the paper
binf1 <- c(0.182, 0.205, 0.485, 0.018, 0.537, 0.037, 0.135,
           0.036, 0.398, 0.434, 0.022, 0.255, 0.098, 0.531, 
           0.025)
bsup1 <- c(0.253, 0.250, 0.569, 0.038, 0.628, 0.068, 0.192,
           0.063, 0.528, 0.572, 0.045, 0.295, 0.159, 0.594,
           0.043)
binf2 <- c(0.183, 0.206, 0.486, 0.019, 0.538, 0.038, 0.136, 
           0.037, 0.399, 0.435, 0.023, 0.257, 0.099, 0.532, 
           0.036)
bsup2 <- c(0.252, 0.248, 0.568, 0.038, 0.627, 0.068, 0.191, 
           0.062, 0.527, 0.572, 0.044, 0.292, 0.158, 0.593,
           0.042)
real_res <- cbind(var, phat, binf1, bsup1, binf2, bsup2, var_type) %>%
  as.data.frame() %>%
  mutate(phat = as.numeric(phat),
         binf1 = as.numeric(binf1),
         bsup1 = as.numeric(bsup1),
         binf2 = as.numeric(binf2),
         bsup2 = as.numeric(bsup2))


# mov variable figure
mov_df2 <- real_res %>% 
  filter(var_type == "mov") %>%
  mutate(var_meaning = case_when(var == "mov_a" ~ "Certainly or Probably",
                                 var == "mov_b" ~ "Probably Not",
                                 var == "mov_c" ~ "Certainly Not",
                                 var == "mov_d" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "Certainly Not", "Probably Not", "Certainly or Probably"), 
    ordered = TRUE))

mov_fig2 <- ggplot(data = mov_df2, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1, color = "Hajek"), alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2, color = "Hajek-A"), alpha = 0.3) +
  labs(x = "Point Estimate",
       title = "Intention to Leave the District") + 
  scale_color_manual(name='Variance',
                     breaks=c('Hajek', 'Hajek-A'),
                     values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# rep variable figure
rep_df2 <- real_res %>% 
  filter(var_type == "rep") %>%
  mutate(var_meaning = case_when(var == "rep_a" ~ "Good",
                                 var == "rep_b" ~ "Fair",
                                 var == "rep_c" ~ "Poor",
                                 var == "rep_d" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "Poor", "Fair", "Good"), ordered = TRUE))

rep_fig2 <- ggplot(data = rep_df2, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1), color = "firebrick", alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2), color = "blue", alpha = 0.3) +
  labs(x = "Point Estimate",
       y = "Variable of Interest",
       title = "Perceived District Reputation") + 
  # scale_color_manual(name='Variance',
  #                    breaks=c('Hajek', 'Hajek-A'),
  #                    values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# wit variable figure
wit_df2 <- real_res %>% 
  filter(var_type == "wit") %>%
  mutate(var_meaning = case_when(var == "wit_a" ~ "Never",
                                 var == "wit_b" ~ "Rarely",
                                 var == "wit_c" ~ "Sometimes",
                                 var == "wit_d" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "Sometimes", "Rarely", "Never"), ordered = TRUE))

wit_fig2 <- ggplot(data = wit_df2, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1, color = "Hajek"), alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2, color = "Hajek-A"), alpha = 0.3) +
  labs(x = "Point Estimate",
       title = "Witnessed Trafficking") + 
  scale_color_manual(name='Variance',
                     breaks=c('Hajek', 'Hajek-A'),
                     values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# wor variable figure
wor_df2 <- real_res %>% 
  filter(var_type == "wor") %>%
  mutate(var_meaning = case_when(var == "wor_a" ~ "Yes",
                                 var == "wor_b" ~ "No",
                                 var == "wor_c" ~ "No Opinion")) %>%
  mutate(var_meaning = factor(var_meaning, levels = c(
    "No Opinion", "No", "Yes"), ordered = TRUE))

wor_fig2 <- ggplot(data = wor_df2, aes(x = phat, y = var_meaning)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = binf1, xmax = bsup1), color = "firebrick", alpha = 0.4) +
  geom_errorbarh(aes(xmin = binf2, xmax = bsup2), color = "blue", alpha = 0.3) +
  labs(x = "Point Estimate",
       y = "Variable of Interest",
       title = "Roadworks in Neighborhood") + 
  # scale_color_manual(name='Variance',
  #                    breaks=c('Hajek', 'Hajek-A'),
  #                    values=c('Hajek'='firebrick', 'Hajek-A'='blue')) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# combine plots
real_fullfig <- ggarrange(rep_fig2, wit_fig2, wor_fig2, mov_fig2, 
                          nrow = 2, ncol = 2)

ggsave("real_fullfig.png")
