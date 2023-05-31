setwd("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/Graphs")

library(sf)
library(geojsonsf)
library(magic)
library(ggplot2)
library(geofacet)
library(ggnewscale)
library(magic)
library(stringr)
library(ggforce)
library(gridExtra)
library(grid)
library(MetBrewer)
library(ggsci)
library(cowplot)
library(TMB)
library(metR)
library(tidyverse)

dyn.load(dynlib("../new"))

rmvnorm_sparseprec <- function(n, mean = rep(0, nrow(prec)), prec = diag(lenth(mean))) {
  z = matrix(rnorm(n * length(mean)), ncol = n)
  L_inv = Matrix::Cholesky(prec)
  v <- mean + Matrix::solve(as(L_inv, "pMatrix"), Matrix::solve(Matrix::t(as(L_inv, "Matrix")), z))
  as.matrix(Matrix::t(v))
}  
var.sample <- function(fit, nsample){
  r <- fit$obj$env$random
  par_f <- fit$par.full[-r]
  
  par_r <- fit$par.full[r]
  hess_r <- fit$obj$env$spHess(fit$par.full, random = TRUE)
  smp_r <- rmvnorm_sparseprec(nsample, par_r, hess_r)
  
  smp <- matrix(0, nsample, length(fit$par.full))
  smp[ , r] <- smp_r
  smp[ ,-r] <- matrix(par_f, nsample, length(par_f), byrow = TRUE)
  colnames(smp) <- rep("NA", length(fit$par.full))
  colnames(smp)[r] <- names(par_r)
  colnames(smp)[-r] <- names(par_f)
  smp
}

fit_tmb <- readRDS("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/fit.rds")
fit.var.sim <- readRDS("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/data/fit.var.sim.rds")

pyears_data <- readRDS("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/data/pyears_data_smoothed.rds")

age.start <- 10; age.end <- 65
year.start <- 1990; year.end <- 2020

aggr.mat.cohort.0 <- filter(pyears_data, 
                            period >= year.start) %>%
  group_by(country) %>% group_split %>%
  setNames(pyears_data$country %>% levels)

joint.countries <- names(aggr.mat.cohort.0)

m.col <- met.brewer(name = "Juarez", n = 7)[7]
f.col <- met.brewer(name = "Juarez", n = 7)[1]

m.fill.col <- met.brewer(name = "Juarez", n = 7)[7]
f.fill.col <- met.brewer(name = "Lakota", n = 6)[3]

pal <- colorRampPalette(c(f.col, m.col))(1000)

bound.sf <- geojson_sf("C:\\Users\\ktang3\\Desktop\\custom.geo.json")

q4515.sim <- lapply(fit.var.sim, function(i){
  list(m = sapply(seq(length(joint.countries)), function(k){1 - apply(exp(-exp(i$mx_mat_m[15:59 - age.start + 1, , k])), 2, prod)}),
       f = sapply(seq(length(joint.countries)), function(k){1 - apply(exp(-exp(i$mx_mat_f[15:59 - age.start + 1, , k])), 2, prod)})
         )
  })

q4515.CI.list <- list(lower = list(),
                      mean = list(),
                      upper = list())

temp <- c()

for(i in 1:length(q4515.sim[[1]])){
  for(j in 1:length(q4515.sim)){
    temp <- cbind(temp, c(q4515.sim[[j]][[i]]))
  }
  q4515.CI.list[[1]][[ names(q4515.sim[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.025)), 
                                                                            dim(q4515.sim[[1]][[i]]), 
                                                                            dimnames = dimnames(q4515.sim[[1]][[i]]))
  q4515.CI.list[[2]][[ names(q4515.sim[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.5)), 
                                                                            dim(q4515.sim[[1]][[i]]), 
                                                                            dimnames = dimnames(q4515.sim[[1]][[i]]))
  q4515.CI.list[[3]][[ names(q4515.sim[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.975)), 
                                                                            dim(q4515.sim[[1]][[i]]), 
                                                                            dimnames = dimnames(q4515.sim[[1]][[i]]))
  temp <- c()
}

q4515.df <- lapply(seq(dim(fit_tmb$mode$mx_mat_f)[3]), function(i){
    bind_rows(tibble(value = 1 - apply(exp(-exp(fit_tmb$mode$mx_mat_f[15:59 - age.start + 1, , i])), 2, prod),
                       year =  year.start:year.end,
                       sex = "f"),
              
              tibble(value = 1 - apply(exp(-exp(fit_tmb$mode$mx_mat_m[15:59 - age.start + 1, , i])), 2, prod),
                     year =  year.start:year.end,
                     sex = "m"))
  }) %>% setNames(joint.countries) %>%
  map2(names(.), ~ add_column(.x, name = .y)) %>%
  bind_rows() %>%
  full_join(
    pyears_data %>% filter(period >= year.start) %>%
      select(country, mm1, period) %>%
      rename(sex = mm1,
             year = period,
             name = country) %>%
      distinct() %>%
      mutate(est = 1,
             sex = str_replace(sex, "female", "f"),
             sex = str_replace(sex, "male", "m"))
    ) %>%
  mutate(est = ifelse(is.na(est), 0, est),
         sex = factor(sex))

q4515.CI.df <- lapply(q4515.CI.list, function(i){
  bind_rows(
    i$f %>% 
      `rownames<-`(year.start:year.end) %>%
      `colnames<-`(joint.countries) %>%
      reshape2::melt() %>%
      mutate(sex = "f"),
    
    i$m %>%
      `rownames<-`(year.start:year.end) %>%
      `colnames<-`(joint.countries) %>%
      reshape2::melt() %>%
      mutate(sex = "m")
    )}) %>%
  map2(names(.), ~ add_column(.x, CI = .y)) %>%
  bind_rows() %>%
  rename(year = Var1,
         name = Var2) %>%
  pivot_wider(values_from = value,
              names_from = CI)

q4515.df.all <- full_join(q4515.df, q4515.CI.df)

plot.sf <- bound.sf %>% mutate(
  name = str_replace(name, "Swaziland", "Eswatini"),
  name = str_replace(name, "Dem. Rep. Congo", "Congo Democratic Republic"),
  name = str_replace(name, "Côte d'Ivoire", "Cote d'Ivoire")
  ) %>% 
  full_join(crossing(name = .$name %>% unique, year = year.start:year.end)) %>%
  full_join(q4515.df) %>%
  mutate(sex = factor(sex, levels = levels(addNA(sex)), labels = c(levels(sex), "miss"), exclude = NULL))

#45q15 map####
plot.sf %>% filter(!sex == "f" ,
                   year %in% seq(1990, 2015, by = 5),
                   region_wb == "Sub-Saharan Africa") %>%
  ggplot() + 
  geom_sf(aes(fill=value), color = "black", size = 0.1) +
  geom_text(aes(x = -5, y = 0, label = year), size = 8, fontface = "bold", family = "sans") +
  scale_fill_viridis_c(limits = c(min(subset(plot.sf, year %in% seq(1990,2015,by=5))$value, na.rm=T), 
                                  max(subset(plot.sf, year %in% seq(1990,2015,by=5))$value, na.rm=T)), option = "B", na.value = "grey90",
                       name = bquote("Males"~""[45]*q[15]),
                       guide = "none") + 
  scale_y_continuous(expand = c(0.03, 0.03)) +
  scale_x_continuous(expand = c(0.03, 0.03)) +
  ggtitle("Males") + 
  #geom_text(data = data.frame(x = -8, y = 40, year = 1990), aes(x = x, y = y), label = "A", fontface = "bold", size = 20) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        #strip.text = element_text(size = 20, face = "bold", hjust = 0),
        strip.background = element_blank(),
        plot.margin = margin(0, 40, 0, 0, "pt"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold")) +
  facet_wrap( ~ year) -> m.map

plot.sf %>% filter(!sex == "m" ,
                   year %in% seq(1990, 2015, by = 5),
                   region_wb == "Sub-Saharan Africa") %>%
  ggplot() + 
  geom_sf(aes(fill=value), color = "black", size = 0.1) +
  geom_text(aes(x = -5, y = 0, label = year), size = 8, fontface = "bold", family = "sans") +
  scale_fill_viridis_c(limits = c(min(subset(plot.sf,year%in%seq(1990,2015,by=5))$value,na.rm=T), 
                                  max(subset(plot.sf,year%in%seq(1990,2015,by=5))$value,na.rm=T)), option = "B", na.value = "grey90",
                       name = bquote("Females"~""[45]*q[15]),
                       guide = "none") + 
  scale_y_continuous(expand = c(0.03, 0.03)) +
  scale_x_continuous(expand = c(0.03, 0.03)) +
  ggtitle("Females") + 
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        #strip.text = element_text(size = 20, face = "bold", hjust = 0),
        strip.background = element_blank(),
        plot.margin = margin(0, 0, 0, 40, "pt"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))+
  facet_wrap( ~ year) -> f.map

plot.sf %>% filter(!sex == "f" ,
                   year %in% seq(1990, 2015, by = 5)) %>%
  ggplot() + 
  geom_sf(aes(fill=value), color = "black") +
  scale_fill_viridis_c(limits = c(min(subset(plot.sf,year%in%seq(1990,2015,by=5))$value,na.rm=T), 
                                  max(subset(plot.sf,year%in%seq(1990,2015,by=5))$value,na.rm=T)), option = "B", na.value = "grey90",
                       name = bquote(""~""[45]*q[15]),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour ="black",
                                              frame.linewidth = 3, ticks.linewidth = 3)) +
  theme(legend.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom",
        legend.margin = margin(-8 ,0 , 0, 0, "cm"),
        legend.text = element_text(size = 25),
        legend.key.width = unit(2, "cm"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt")) -> legend.plot

map.legend <- cowplot::get_legend(legend.plot)
plot_grid(m.map, f.map, labels = c("A", "B"), label_size = 40, label_y = 0.87) -> top.row
plot_grid(top.row, map.legend, ncol = 1, rel_heights = c(2, 0.25)) -> mf.map.plot

ggsave(filename = "45q15 m and f.png", device = "png", width = 15.4686, height = 8.8946667, unit = "cm", scale = 3.15,
       plot = mf.map.plot)

#45q15 geofacet####
grid.sf <- africa_countries_grid1 %>%
  mutate(name = str_replace(name, "Central African Republic", "Central African Rep."),
         name = str_replace(name, "Democratic Republic of the Congo", "Congo Democratic Republic"),
         name = str_replace(name, "Republic of the Congo", "Congo"),
         name = str_replace(name, "C?te d'Ivoire", "Cote d'Ivoire"),
         name = str_replace(name, "Equatorial Guinea", "Eq. Guinea"),
         name = str_replace(name, "South Sudan", "S. Sudan"),
         name = str_replace(name, "West Sudan", "W. Sudan"),
         #name = str_replace(name, "Somalia", "Somaliland"),
         name = str_replace(name, "S?o Tom? and Principe", "Sao Tome and Principe"),
         name = str_replace(name, "Central African Rep.", "Central African Republic"),
         code = str_replace(code, "NAM", "NA")
  ) %>%
  left_join(data.frame(code = bound.sf$iso_a2, iso_a3 = bound.sf$iso_a3) %>% distinct()) %>%
  mutate(iso_a3 = ifelse(name == "Comoros", "COM", iso_a3),
         iso_a3 = ifelse(name == "Sao Tome and Principe", "STP", iso_a3),
  )

ggplot(q4515.df.all %>%
         full_join(grid.sf) %>%
         filter(sex %in% c("m","f") & est == 1,
                name != "Comoros") %>%
         rowwise() %>%
         mutate(name = paste(strwrap(name, 20), collapse = "\n")) %>%
         ungroup) + 
  geom_line(aes(x = year, y = value, color = sex), size = 1.1) + 
  geom_text(data = grid.sf %>% 
              filter(!name %in% c("Cabo Verde", "Mauritius", "Seychelles", "Comoros"),
                     row > 1) %>%
              rowwise() %>%
              mutate(name = paste(strwrap(name, 20), collapse = "\n")),
            x = 1990, y = 0.9, aes(label = name), size = 8.5, fontface = "bold", family = "sans", hjust = 0) +
  labs(tag = "C") + ggtitle("\n\n\n\n") +
  theme_bw() + scale_y_continuous(trans = "identity", limits = c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, "0.50", 0.75, 1), expand = c(0.05, 0.05)) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = sex), alpha = 0.45) +
  scale_color_manual(values = c("m" = m.col, "f" = f.col), label = c("Males        ", "Females        ")) +
  scale_fill_manual(values = c("m" = m.fill.col, "f" = f.fill.col), guide = "none") +
  scale_linetype_manual(values = c("solid", "solid"), guide = "none") +
  theme(legend.key.width = unit(3, "cm"),
        legend.text = element_text(size = 35),
        title = element_text(size = 20), plot.title = element_text(size = 30, hjust = 0.5),
        #strip.text = element_text(size = 25, face = "bold"),
        strip.text = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = "top", legend.background = element_blank(), legend.box = "vertical",
        aspect.ratio = 0.8,
        strip.background = element_blank(),
        #panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.text.x = element_text(angle = 40, size = 25, hjust = 1),
        axis.text.y = element_text(size = 25),
        plot.tag.position = c(0, 0.95),
        plot.tag = element_text(size = 90, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  guides(color = guide_legend(title = "", override.aes = list(color = c(m.col, f.col), size = 2))) +
  xlab("") + ylab("") +
  #ylab(bquote(""[45]*q[15])) +
  #ggtitle(bquote("Estimated Period"~""[45]*q[15]~"1990-2017")) +
  facet_geo( ~ name, grid = grid.sf %>% 
               filter(!name %in% c("Cabo Verde", "Mauritius", "Seychelles", "Comoros"),
                      row > 1) %>% 
               select(name, iso_a3, row, col) %>%
               rowwise() %>%
               mutate(row = row - 1,
                      name = paste(strwrap(name, 20), collapse = "\n")) %>%
               rename(code = iso_a3) %>%
               ungroup %>%
               rows_update(tibble(code = c("SWZ", "ZAF"), row = c(7, 8), col = c(7, 5)), by = "code")) -> q4515.geofacet

#plot_grid(NULL, q4515.geofacet, labels = c("", "B"), rel_heights = c(0.1, 1), label_size = 40, label_y = 1.01) -> q4515.geofacet.with.B

ggsave(filename="45q15 geofacet.png",device="png",width=15.4686,height=8.898466667,unit = "cm",scale = 7,plot = q4515.geofacet)
#knitr::plot_crop("45q15 geofacet.pdf")

#45q15 csv####
q4515.df %>% 
  filter(est == 1) %>%
  mutate(period = paste(floor(year / 5) * 5, "-", floor(year / 5) * 5 + 4)) %>%
  group_by(sex, name, period) %>%
  filter(!is.na(value),
         year == min(year)) %>%
  ungroup %>%
  bind_rows(crossing(name = joint.countries,
                     value = NA,
                     sex = c("m", "f"),
                     period = paste(seq(1990, 2015, by = 5), "-", seq(1994, 2019, by = 5)))) %>%
  group_by(sex, name, period) %>%
  summarise_at(vars(value, year), first) %>%
  mutate(rate = (value / lag(value)) ^ (1 / (year - lag(year))) - 1) %>%
  select(sex, name, period, rate) %>%
  pivot_wider(names_from = period, values_from = rate) %>%
  ungroup %>%
  mutate(region = NA,
         region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
         region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
         region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
         region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
         region = fct_relevel(region, "Southern", "Eastern", "Central", "Western")
  ) %>% 
  arrange(region, rev(name)) %>%
  select(-`1990 - 1994`) %>%
  rename_with(~paste(as.numeric(str_extract(.x, "\\d{4}")) - 5, "-", as.numeric(str_extract(.x, "\\d{4}"))),
              .cols = -c(1:2, 9)) -> q45.df.rate

write.csv(q45.df.rate, "q45 rate of change.csv", row.names = F)


q4515.df %>% 
  filter(est == 1) %>%
  mutate(period = paste((year + 5) %/% 10 * 10 - 5, "-", (year + 5) %/% 10 * 10 + 4)) %>%
  group_by(sex, name, period) %>%
  filter(!is.na(value),
         year %in% c(min(year), max(year))) %>%
  ungroup-> rate.df.intermediate

rate.df.intermediate %>%
  group_by(sex, name, period) %>%
  summarise_at(vars(year), n_distinct) %>%
  filter(year == 1) %>%
  mutate(year = as.numeric(str_extract(period, "\\d{4}")) + 9,
         value = NA) %>%
  ungroup -> rate.df.dummy

rate.df.intermediate %>%
  bind_rows(rate.df.dummy) %>%
  group_by(sex, name, period) %>%
  arrange(name, year, period) %>%
  summarise(rate = (last(value) / first(value)) ^ (1 / (max(year) - min(year))) - 1) %>%
  ungroup %>%
  arrange(period) %>%
  select(sex, name, period, rate) %>%
  pivot_wider(names_from = period, values_from = rate) %>%
  mutate(region = NA,
         region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
         region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
         region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
         region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
         region = fct_relevel(region, "Southern", "Eastern", "Central", "Western")
  ) %>% 
  arrange(region, rev(name)) -> q45.df.rate.decennial

write.csv(q45.df.rate.decennial, "q45 rate of change decennial.csv", row.names = F)

write.csv(q4515.df %>% filter(sex == "m",
                              est == 1) %>%
            mutate(region = NA,
                   region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
                   region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
                   region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
                   region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
                   region = fct_relevel(region, "Southern", "Eastern", "Central", "Western")
            ) %>%
            arrange(region) %>%
            pivot_wider(names_from = year, values_from = value), "q4515 excel m.csv", row.names = F)

write.csv(q4515.df %>% filter(sex == "f",
                              est == 1) %>%
            mutate(region = NA,
                   region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
                   region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
                   region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
                   region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
                   region = fct_relevel(region, "Southern", "Eastern", "Central", "Western")
            ) %>%
            arrange(region) %>%
            pivot_wider(names_from = year, values_from = value), "q4515 excel f.csv", row.names = F)


#stacked plots####
fit.var.sim.reduced <- lapply(fit.var.sim, function(i){i[!sapply(lapply(i, dim), is.null)]})

fit.var.sim.CI.list <- list(lower = list(),
                            mean = list(),
                            upper = list())
temp <- c()

for(i in 1:length(fit.var.sim.reduced[[1]])){
  for(j in 1:length(fit.var.sim.reduced)){
    temp <- cbind(temp, c(fit.var.sim.reduced[[j]][[i]]))
  }
  fit.var.sim.CI.list[[1]][[ names(fit.var.sim.reduced[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.025)), 
                                                                            dim(fit.var.sim.reduced[[1]][[i]]), 
                                                                            dimnames = dimnames(fit.var.sim.reduced[[1]][[i]]))
  fit.var.sim.CI.list[[2]][[ names(fit.var.sim.reduced[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.5)), 
                                                                            dim(fit.var.sim.reduced[[1]][[i]]), 
                                                                            dimnames = dimnames(fit.var.sim.reduced[[1]][[i]]))
  fit.var.sim.CI.list[[3]][[ names(fit.var.sim.reduced[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.975)),
                                                                            dim(fit.var.sim.reduced[[1]][[i]]), 
                                                                            dimnames = dimnames(fit.var.sim.reduced[[1]][[i]]))
  temp <- c()
}

mort.df <- bind_rows(
  fit_tmb$mode$mx_mat_f %>%
    `dimnames<-`(list(age.start:age.end, year.start:year.end, joint.countries)) %>%
    reshape2::melt() %>%
    mutate(sex = "f"),
  
  fit_tmb$mode$mx_mat_m %>%
    `dimnames<-`(list(age.start:age.end, year.start:year.end, joint.countries)) %>%
    reshape2::melt() %>%
    mutate(sex = "m")
  ) %>%
  rename(age = Var1,
         year = Var2,
         name = Var3) %>%
  full_join(
    pyears_data %>% filter(period >= year.start) %>%
      select(country, mm1, period) %>%
      rename(sex = mm1,
             year = period,
             name = country) %>%
      distinct() %>%
      mutate(est = 1,
             sex = str_replace(sex, "female", "f"),
             sex = str_replace(sex, "male", "m"))
  ) %>%
  mutate(est = ifelse(is.na(est), 0, est),
         sex = factor(sex))
  
mort.CI.df <- lapply(fit.var.sim.CI.list, function(i){
  bind_rows(
    i$mx_mat_f %>% 
      `dimnames<-`(list(age.start:age.end, year.start:year.end, joint.countries)) %>%
      reshape2::melt() %>%
      mutate(sex = "f"),
    
    i$mx_mat_m %>% 
      `dimnames<-`(list(age.start:age.end, year.start:year.end, joint.countries)) %>%
      reshape2::melt() %>%
      mutate(sex = "m")
  )}) %>%
  map2(names(.), ~ add_column(.x, CI = .y)) %>%
  bind_rows() %>%
  rename(age = Var1,
         year = Var2,
         name = Var3) %>%
  pivot_wider(values_from = value,
              names_from = CI)

mort.df.all <- full_join(mort.df, mort.CI.df)

plot.heatmap.coord.flip.vertical <- function(dat.list, agerange, country, log = T, flip.coord = T, region, xlab1 = "Age", ylab1 = "Year",  xlab2 = "Age", ylab2 = "Year"){
   m.dat <- mort.df %>%
     filter(name == country,
            age %in% agerange,
            sex == "m",
            year %in% year.start:year.end)
   
   f.dat <- mort.df %>%
     filter(name == country,
            age %in% agerange,
            sex == "f",
            year %in% year.start:year.end)
   
   if(log == F) {m.dat <- m.dat %>% mutate(value = exp(value)); f.dat <- f.dat %>% mutate(value = exp(value))}
   
   m.mat <- m.dat %>%
     full_join(crossing(age = agerange,
                       year = year.start:year.end))
   
   f.mat <- f.dat %>%
     full_join(crossing(age = agerange,
                        year = year.start:year.end))
   
  a <- ggplot(m.mat) +
    geom_tile(aes(x = age, y = year, fill = value), width = 1) + 
    theme_bw() + ylab(ylab1) + xlab(xlab1) + 
    ggtitle(paste(country, region)) +
    scale_fill_viridis_c(option = "B",na.value = "transparent") +
    geom_contour(aes(x = age, y = year, z = value), color = "white", lwd = 0.8) +
    geom_text_contour(aes(x = age, y = year, z = value), color = "white", size = 7, fontface = "bold") +
    labs(tag = "Male") +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    theme(legend.position = "none",legend.key = element_blank(),legend.background = element_blank(),legend.title=element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size=15),
          legend.key.width = unit(2,"cm"),
          text = element_text(size=20),
          plot.title = element_text(hjust=0.5, face = "bold"),
          axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 25),
          plot.tag.position = c(0.85, 0.15),
          plot.tag = element_text(size = 25, face = "bold", color = "white"),
          plot.margin = margin(0, 0.25, 0, 0.25, "cm"),
          axis.ticks = element_line(size = 1.2))
  
  b <- ggplot(f.mat) +
    geom_tile(aes(x = age, y = year, fill = value), width = 1) + 
    theme_bw() + ylab(ylab2) + xlab(xlab2) +
    #ggtitle(paste(country, region))+
    scale_fill_viridis_c(option = "B",na.value = "transparent") + 
    geom_contour(aes(x = age, y = year, z = value), color = "white", lwd = 0.8) +
    geom_text_contour(aes(x = age, y = year, z = value), color = "white", size = 6, fontface = "bold") +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    labs(tag = "Female") +
    theme(legend.position = "none",legend.key = element_blank(),legend.background = element_blank(),legend.title = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.border=element_blank(),
          axis.line = element_line(colour = "black"),
          legend.text=element_text(size=15),
          legend.key.width=unit(2,"cm"),
          text = element_text(size=20),
          plot.title = element_text(hjust=0.5, face = "bold"),
          axis.text = element_text(size = 20),
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
          axis.title = element_text(size = 25),
          plot.tag.position = c(0.85, 0.27),
          plot.tag = element_text(size = 25, face = "bold", color = "white"),
          plot.margin = margin(0, 0.25, 0, 0.25, "cm"))
  
  if(flip.coord) {a <- a + coord_flip(); b <- b + coord_flip()}
  
  grid.arrange(a,b,nrow=2)
}

a.vert <- plot.heatmap.coord.flip.vertical(fit_tmb, agerange = 15:60, "South Africa", log = F, region = "", xlab1 = "", ylab1 = "", xlab2 = "", ylab2 = "")
b.vert <- plot.heatmap.coord.flip.vertical(fit_tmb, agerange = 15:60, "Kenya", log = F, region = "", xlab1 = "", ylab1 = "", xlab2 = "", ylab2 = "")
c.vert <- plot.heatmap.coord.flip.vertical(fit_tmb, agerange = 15:60, "Gabon", log = F, region = "", xlab1 = "", xlab2 = "", ylab1 = "", ylab2 = "")
d.vert <- plot.heatmap.coord.flip.vertical(fit_tmb, agerange = 15:60, "Sierra Leone", log = F, region = "", ylab1 = "", ylab2 = "")

stacked.heatmaps <- grid.arrange(d.vert, c.vert, b.vert, a.vert, ncol = 4)

plot.cross.section.age.pattern.vertical <- function(country.set, y_scale = "identity", region, ylab = "Mortality rate per 1000"){
  df <- mort.df.all %>% filter(year %in% c(1995, 2005, 2015),
                                  age %in% 15:60,
                                  name %in% country.set) %>%
    mutate(est = factor(est, levels = c("0", "1")),
           value = exp(value),
           lower = exp(lower),
           upper = exp(upper))
  
  ggplot(df) +
    geom_line(data = tibble(est = factor(c(0, 1), levels = c("0", "1")), value = NA, age = 30), aes(x = age, y = value, linetype = est)) +
    geom_line(aes(x = age, y = value, color = sex, linetype = est), size = 1.2) + theme_bw() + 
    scale_y_continuous(trans = y_scale, 
                       limits = c(1, 75)/1000,
                       breaks = c(1, 5, 25) / 1000,
                       #breaks = exp(seq(log(min(df$X2.5.)), log(max(df$X97.5.) + 0.02 * (max(df$X97.5.) - min(df$X2.5.))), length.out = 4)),
                       labels = function(i){signif(i*1000, digits = 3)}, position = "left", expand = c(0.02, 0.02))+
    geom_hline(yintercept = 0) +
    geom_ribbon(aes(x = age, ymin = lower, ymax = upper, fill = sex), alpha = 0.3)+
    geom_text(aes(x = 22, y = 0.8 * 75/1000, label = year), size = 8, fontface = "bold", family = "sans") +
    scale_linetype_manual(values=c("dashed","solid"), guide = "none") +
    scale_color_manual(values=c(f.col, m.col),label=c("Females","Males"), guide = "none") +
    scale_fill_manual(values=c(f.fill.col, m.fill.col),label=c("Females","Males"),guide="none") +
    ggtitle(country.set) +
    theme(legend.key.width = unit(1.5,"cm"),
          legend.text = element_text(size=20),
          plot.title = element_text(size = 25, hjust=0.5, face = "bold"),
          strip.text = element_blank(),
          #strip.text.x = element_text(size=25),
          axis.text = element_text(size=20),
          axis.title = element_text(size = 25),
          strip.background = element_blank(),
          legend.position="bottom",
          legend.background=element_blank(),legend.box="vertical",
          legend.title = element_text(size = 25),
          panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          panel.border = element_blank(),
          panel.spacing.y = unit(0.5, "cm"), 
          axis.line = element_line(colour = "black"))+
    #guides(col = guide_legend(order = 1,reverse=T,title="Sex"))+
    xlab("Age")+ylab(ylab) +
    #ggtitle(paste("Estimated mortality schedules in selected years",region))+
    facet_wrap(~ year, ncol = 1) -> a
  a
}

a.age.vert <- plot.cross.section.age.pattern.vertical(country.set=c("South Africa"),region="on log scale \n(Southern)",y_scale="log", ylab = "")
b.age.vert <- plot.cross.section.age.pattern.vertical(country.set=c("Kenya"),region="on log scale \n(Eastern)",y_scale="log", ylab = "")
c.age.vert <- plot.cross.section.age.pattern.vertical(country.set=c("Gabon"),region="on log scale \n(Central)",y_scale="log", ylab = "")
d.age.vert <- plot.cross.section.age.pattern.vertical(country.set=c("Sierra Leone"),region="on log scale \n(Western)",y_scale="log")

age.legend <- cowplot::get_legend( ggplot(mort.df %>% filter(year %in% c(1995, 2005, 2015),
                                                                    age %in% 15:60,
                                                                    name %in% "Sierra Leone") %>%
                                            mutate(est = factor(est, levels = c("0", "1")))) +
                                     geom_line(data = tibble(est = factor(c(0, 1), levels = c("0", "1")), value = NA, age = 30), 
                                               aes(x = age, y = value, linetype = est)) +
                                     geom_line(aes(x = age, y = value, color = sex, linetype = est), size = 1.2) + theme_bw() + 
                                     scale_linetype_manual(values=c("dashed","solid"), guide = "none")+
                                     scale_color_manual(values=c(f.col, m.col),label=c("Females","Males"), guide = "none")+
                                     theme(legend.key.width = unit(1.5,"cm"),
                                           legend.text = element_text(size=25),
                                           legend.position="bottom",legend.background=element_blank(),legend.box="vertical",
                                           legend.title = element_text(size = 25),
                                           panel.grid = element_blank(),
                                           plot.margin = margin(0, 0, 0, 0, "cm"))+
                                     guides(col = guide_legend(order = 1,reverse=T,title="Sex"))+
                                     facet_wrap(~year, ncol = 1) )

grid.arrange(d.age.vert, c.age.vert, b.age.vert, a.age.vert, ncol = 4) -> top.row.age

plot_grid(top.row.age, age.legend, ncol = 1, rel_heights = c(2, 0.25)) -> age.plot

plot_grid(NULL, stacked.heatmaps, NULL, age.plot, labels = c("", "A", "", "B"), ncol = 1, rel_heights = c(0.1, 0.75, 0.1, 1), 
          label_size = 40, label_y = 1.1) -> stacked.plot.all

ggsave(filename = "stacked plot.png", device = "png", width = 16, height = 24, plot = stacked.plot.all)

#Tips####
tp.df <- fit_tmb$mode$tips_params %>%
  `colnames<-`(joint.countries) %>%
  reshape2::melt() %>%
  rename(tp = Var1, 
         name = Var2) %>%
  mutate(tp = tp - 1,
         tp = as.factor(tp)) %>%
  group_by(name) %>%
  mutate(value = value + fit_tmb$mode$tips_params_common_after)

tp.common <- fit_tmb$mode$tips_params_common_after

tp.df %>%
  ggplot(aes(x = as.factor(tp), y = exp(value))) +
  geom_violin(col = "white", fill = m.fill.col, alpha = 0.5, width = 1.5) +
  geom_boxplot(width = 0.1, outlier.size = 2) +
  geom_line(data = tibble(tp = 1:15, value = tp.common, type = "Common Trend"), aes(x = tp, y = exp(value), col = type),
            size = 1, alpha = 0.9) +
  scale_color_manual(values = rgb(colorRamp(c("black", "orchid"))(0.5)/255), name = "",
                     guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip() +
  theme_bw() +
  xlab("Time Prior to Survey Year") + ylab("") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(size = 15),
        legend.position = c(0.88, 0.08),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 0.75)

tp.df %>%
  mutate(tp = factor(tp)) %>%
  ggplot(aes(x = tp, y = exp(value), group = tp)) +
  geom_violin(trim = FALSE, color = NA, fill = "grey80") +
  geom_dotplot(col = m.col, fill = m.fill.col, binaxis = "y", stackdir ="center", alpha = 0.8, dotsize = 0.5, binpositions = "all", binwidth = 1/40) +
  scale_x_discrete(labels = paste("\n", 0:14, "\n")) +
  #scale_fill_manual(values = rev(met.brewer(name = "Hokusai2", n = 15)), guide = "none") +
  #stat_summary(fun = mean, aes(ymin=..y.., ymax=..y..), color = "black", geom="errorbar", width = 1, size=0.75) +
  stat_summary(data = tibble(tp = as.factor(0:14), value = tp.common, type = "Common Trend"), 
               fun = mean, aes(ymin=..y.., ymax=..y..), geom="errorbar",
               col =  rgb(colorRamp(c("black", "purple"))(0.75)/255),
               size = 2.5, alpha = 0.9, width = 1) +
  geom_line(data = tibble(tp = 1, value = NA, type = "Common Trend"), aes(col = type), size = 2, key_glyph = "vline") +
  scale_color_manual(values = rgb(colorRamp(c("black", "purple"))(0.75)/255), name = "",
                     #guide = guide_legend(override.aes = list(alpha = 1))) +
                     guide = "none") +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip() +
  theme_bw() +
  xlab("Time Prior to Survey Year") + ylab("Odds ratio (log-scale)") +
  theme(axis.text = element_text(size = 25),
        #axis.text.x = element_text(size = 20, margin = margin(t = 0.3, b = 0.3, unit = "cm")),
        axis.title = element_text(size = 25),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(2, "cm"),
        legend.text = element_text(size = 25),
        legend.position = c(0.88, 0.08),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 0.75) -> tp.plot.new


ggsave(filename="tips violin 3.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot = tp.plot.new)
#ggsave(filename="tp geo plot.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=tp.geo.plot)

# library(xtable)
# 
# tp.df %>% 
#   mutate(value = round(exp(value), digits = 3),
#          X2.5. = round(exp(X2.5.), digits = 3),
#          X97.5. = round(exp(X97.5.), digits = 3)) %>%
#   unite("CI", X2.5.:X97.5., sep = ", ") %>%
#   mutate(value.CI = paste0(value, " (", CI, ")")) %>%
#   select(country, tp, value.CI) %>%
#   pivot_wider(names_from = country, values_from = value.CI) -> tp.table
# 
# for(i in split(1:length(tp.all.df$country %>% unique), ceiling(1:length(tp.all.df$country %>% unique) / 7))){
#   print(xtable(tp.table %>% select(c(1, i+1))), include.rownames = F)}
# 
# write.csv(tp.all.df %>% relocate(country, tp) %>%
#             mutate(value = exp(value),
#                    X2.5. = exp(X2.5.),
#                    X97.5. = exp(X97.5.)), "TP excel.csv", row.names = F)
# library(xlsx)
# write.xlsx(tp.all.df %>% relocate(country, tp) %>%
#              mutate(value = exp(value),
#                     X2.5. = exp(X2.5.),
#                     X97.5. = exp(X97.5.)), file = "TP excel")

#45q15 compare regional GBD WPP####
unpdhhd <- readxl::read_excel("C:/Users/ktang3/Desktop/Imperial/SSA_mort/DemoData-LifeTableInput20170206_JE-edit.xlsm", sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0)
unpdhhd <- subset(unpdhhd, LocType=="Country" & DataType=="Household deaths" & LocAreaType== "Whole area" & is.na(AgeGroup)==FALSE )

## The columns AgeGroup, PopulationInput, DeathsInput are the raw input count data (with unknown ages), columns Exposure and Deaths have the unknown ages proportionately redistributed. Mx are the raw mortality rates.
unpdhhd <- unpdhhd[c("Location", "DataCatalogName", "DataProcess", "DataProcessType","FieldWorkMiddle", "StatisticalConcept", 
                     "DataSourceName","DataSourceYear", "DataSourceShortName","Sex", "YearStart", "YearLength", "TimeStart", 
                     "TimeEnd","AgeStart", "AgeEnd", "AgeSpan", "AgeGroup", "PopulationInput","DeathsInput", "ExposureInput", 
                     "Exposure", "Deaths", "Mx", "AgeGroup2","AgeStart2", "AgeEnd2", "AgeSpan2", "mx", "qx", "lx", "dx", "Lx","sx", "Tx", "ex", "ax")]

## Namibia 2001 census appears twice. Use the first 
unpdhhd <- subset(unpdhhd, !(Location == "Namibia" & DataSourceYear == 2015))

## Remove DHS and AIS (analysed below)
unpdhhd <- subset(unpdhhd, !DataProcess %in% c("AIS", "DHS"))

## Remove aggregate 'Both sexes' entries
unpdhhd <- subset(unpdhhd, Sex %in% c("Male", "Female"))

unpdhhd <- unpdhhd %>%
  mutate(
    year = floor(FieldWorkMiddle),
    name = Location,
    source = DataCatalogName,
    sex = tolower(Sex),
    agestart = AgeStart,
    ageend = AgeEnd - 1,
    agespan = AgeSpan,
    type = DataProcess,
    deaths = Deaths,
    pys = Exposure
  ) %>%
  select(c("name", "year", "source", "type", "sex", "agestart", "ageend", "agespan", "pys", "deaths")) %>%
  mutate(
    name = replace(name, name=="United Republic of Tanzania", "Tanzania"),
    name = replace(name, name=="C?te d'Ivoire", "Cote d'Ivoire"),
    name = replace(name, name=="Swaziland", "Eswatini")
  )

## Restrict countries
unique(unpdhhd[c("name", "year", "source", "type")])
unpdhhd <- subset(unpdhhd, name %in% joint.countries & agestart %in% 10:65)

## Calculate q4515
unpdhhd$mx <- unpdhhd$deaths / unpdhhd$pys
unpdhhd_qx <- aggregate(cbind(cummx4515=5*mx, cummx3515=ifelse(agestart >= 50, 0, 5*mx)) ~ name+year+source+type+sex, subset(unpdhhd, agestart %in% 15:59), sum)
unpdhhd_qx$q4515 <- 1-exp(-unpdhhd_qx$cummx4515)
unpdhhd_qx$q3515 <- 1-exp(-unpdhhd_qx$cummx3515)
unpdhhd_qx <- unpdhhd_qx[with(unpdhhd_qx, order(name, sex, year)),]
unpdhhd_qx$sex <- factor(unpdhhd_qx$sex,labels=c("f","m"))
unpdhhd_qx$name<-factor(unpdhhd_qx$name,levels = joint.countries)
unpdhhd_qx$type<-factor(unpdhhd_qx$type,levels=c("Census","GFFS","Annual HH survey","Survey","LSMS","Multiround surv","Panel","WFS"))

unpdhhd$UNPD.qx<-1-exp(-5*unpdhhd$mx)
unpdhhd$age<-0.5*(unpdhhd$agestart+unpdhhd$ageend)

WPP<-read.csv(file="C:/Users/ktang3/Desktop/Imperial/SSA_mort/WPP 45q15.csv",header=T)
WPP$X<-str_trim(WPP$X)
WPP$X<-str_replace(WPP$X,"Democratic Republic of the Congo","Congo Democratic Republic")
#WPP$X<-str_replace(WPP$X,"Cote d'Ivoire","Cote d'Ivoire")
WPP$X<-str_replace(WPP$X,"United Republic of Tanzania","Tanzania")
names(WPP)<-c("ISO code","country","year","note","f","m")
WPP<-subset(WPP,country%in%joint.countries)
WPP$year<-rep(c(1990,1995,2000,2005,2010,2015),length(joint.countries))
WPP<-rbind.fill(WPP,cbind.data.frame(country=rep(joint.countries,each=6),year=rep(c(1992,1997,2002,2007,2012,2017),length(joint.countries))))
WPP<-WPP[order(WPP$country,WPP$year),]

WPP.45q<-subset(WPP[rep(1:nrow(WPP),each=6),c(2,5,6)],country%in%joint.countries)
WPP.45q$year<-rep(c(1990:1995,rep(1995,6),1995:2000,rep(2000,6),2000:2005,rep(2005,6),2005:2010,rep(2010,6),2010:2015,rep(2015,6),2015:2020,rep(2020,6)),length(joint.countries))
#WPP.45q<-rbind(WPP.45q,cbind.data.frame(country=rep(joint.countries,each=6),year=rep(c(1995,2000,2005,2010,2015,2020),length(joint.countries)),f=NA,m=NA))
WPP.45q$q<-"\"\"[45] * q[15]"
WPP.45q$source="WPP"
WPP.45q<-reshape2::melt(WPP.45q,id.vars=c("country","year","q","source"))
WPP.45q$value<-as.numeric(WPP.45q$value)/1000
names(WPP.45q)<-c("name","year","q","source","sex","mle")
WPP.45q$est<-1

GBD<-read.csv(file="C:/Users/ktang3/Desktop/Imperial/SSA_mort/GBD 45q15.csv",header=T)
GBD$location_name<-str_replace(GBD$location_name,"Democratic Republic of the Congo","Congo Democratic Republic")
GBD$location_name<-str_replace(GBD$location_name,"Côte d'Ivoire","Cote d'Ivoire")
GBD$location_name<-str_replace(GBD$location_name,"United Republic of Tanzania","Tanzania")
GBD.45q<-subset(GBD,location_name%in%joint.countries & year_id %in% 1990:2017)[,c(2,4,7,12)]
GBD.45q<-subset(GBD.45q,sex_name!="both")
names(GBD.45q)<-c("name","sex","year","GBD.q")
GBD.45q$sex<-str_replace(GBD.45q$sex,"female","f")
GBD.45q$sex<-str_replace(GBD.45q$sex,"male","m")

q4515.df.all.compare <- q4515.df.all %>%
  full_join(GBD.45q) %>%
  full_join(WPP.45q %>% select(name, year, sex, mle)) %>%
  pivot_longer(cols = c(1, 9, 10),
               names_to = "source") %>%
  mutate(source = factor(source, levels = c("value", "GBD.q", "mle", "UNPD"), labels = c("Spline", "GBD", "WPP", "UNPD")),
         sex = factor(sex),
         sex = fct_relevel(sex, "m"))

plot.45q.compare <- function(country.set, region){
  ggplot(q4515.df.all.compare %>% filter(name %in% country.set,
                           year %in% 1990:2017,
                           source != "UNPD",
                           est == 1)) +
    geom_line(aes(x = year, y = value, color = sex, linetype = source), size = 1.2) + 
    theme_bw() + scale_y_continuous(trans = "identity", limits = c(0, NA)) +
    geom_ribbon(data = q4515.df.all.compare %>% filter(name %in% country.set,
                                         !is.na(lower),
                                         source == "Spline",
                                         est == 1),
                aes(x = year, ymin = lower, ymax = upper, fill = sex), alpha = 0.3) +
    geom_point(data = unpdhhd_qx %>% filter(name %in% country.set,
                                            year %in% 1990:2017),
               aes(x = year, y = q4515, pch = type, color = sex), size = 4) +
    scale_color_manual(values = c(m.col, f.col),label = c("Males", "Females")) +
    scale_fill_manual(values = c(m.fill.col, f.fill.col), guide = "none") +
    scale_linetype_manual(values = c("solid", "dotted", "dashed"), name = "Estimates") +
    scale_shape_manual(values = c(16:18, 7:11), labels = levels(unpdhhd_qx$type), name = "Household \ndeaths") +
    theme(legend.key.width = unit(1.5,"cm"),legend.text = element_text(size=20),
          title = element_text(size = 20),
          plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
          strip.text= element_text(size = 20, face = "bold"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.y = element_blank(), axis.text = element_text(size = 20),
          axis.title = element_text(size = 30),
          legend.position = "right", legend.background = element_blank(), legend.box = "vertical") +
    guides(color = guide_legend(title = "Sex", order = 1, override.aes = list(color = c(m.col, f.col), size = 1.3)),
           linetype = guide_legend(order = 2)) +
    xlab("") + ylab(bquote(""[45]*q[15])) +
    #ggtitle(bquote("Estimated Period"~""[45]*q[15]~"1990-2017"~.(region))) + 
    facet_wrap( ~ name, scales = "free")
}

#plot.45q.compare(country.set=c("Zimbabwe","Zambia","South Africa","Comoros"),region="(Southern)")

ggsave(filename="45q15 compare GBD WPP UNPD southern.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.45q.compare(country.set=c("Zimbabwe","Zambia","South Africa","Namibia"),region="(Southern)"))
ggsave(filename="45q15 compare GBD WPP UNPD eastern.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.45q.compare(country.set=c("Kenya","Uganda","Malawi","Tanzania"),region="(Eastern)"))
ggsave(filename="45q15 compare GBD WPP UNPD central.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.45q.compare(country.set=c("Congo","Cameroon","Gabon","Congo Democratic Republic"),region="(Central)"))
ggsave(filename="45q15 compare GBD WPP UNPD western.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.45q.compare(country.set=c("Sierra Leone","Nigeria","Cote d'Ivoire","Mali"),region="(Western)"))

#35q15 compare####
GBD.age <- read.csv(file="C:/Users/ktang3/Desktop/Imperial/SSA_Mort/GBD age specific.csv",header=T)
GBD.age <- GBD.age %>%
  mutate(name = str_replace(name, "Democratic Republic of the Congo","Congo Democratic Republic"),
         name = str_replace(name, "United Republic of Tanzania","Tanzania"),
  ) %>%
  filter(year %in% 1900:2017,
         sex != "both",
         measure_name == "Probability of death",
         name %in% joint.countries) %>%
  mutate(age = fct_relevel(age_group_name, "<1 year", "1 to 4", "5 to 9", "10 to 14"),
         age = factor(age, labels = c(0.2, 2.5, seq(7, 92, by = 5))),
         age = as.numeric(levels(age))[as.numeric(age)]
  )

WPP.age <- read.csv(file="C:/Users/ktang3/Desktop/Imperial/SSA_Mort/WPP age specific.csv",header=T)
WPP.age <- WPP.age %>%
  mutate(name = str_replace(name, "Democratic Republic of the Congo", "Congo Democratic Republic"),
         name = str_replace(name, "Côte d'Ivoire", "Cote d'Ivoire"),
         name = str_replace(name, "United Republic of Tanzania", "Tanzania"),
         age = factor(AgeGrp, labels = c(0, 2, 7, 100, seq(17, 97, by = 5), 12)),
         age = as.numeric(levels(age))[as.numeric(age)]) %>%
  filter(name %in% joint.countries,
         MidPeriod %in% 1990:2017,
         Sex != "Total",
         age %in% 10:65) %>%
  .[,c(1,3,4,8,9,18)] %>%
  `names<-`(c("name","year","sex","WPP.mx","WPP.qx","age")) %>%
  mutate(sex = str_replace(sex, "Female", "female"),
         sex = str_replace(sex, "Male", "male"))

q3515.df <- lapply(seq(dim(fit_tmb$mode$mx_mat_f)[3]), function(i){
  bind_rows(tibble(value = 1 - apply(exp(-exp(fit_tmb$mode$mx_mat_f[15:49 - age.start + 1, , i])), 2, prod),
                   year =  year.start:year.end,
                   sex = "f"),
            
            tibble(value = 1 - apply(exp(-exp(fit_tmb$mode$mx_mat_m[15:49 - age.start + 1, , i])), 2, prod),
                   year =  year.start:year.end,
                   sex = "m"))
}) %>% setNames(joint.countries) %>%
  map2(names(.), ~ add_column(.x, name = .y)) %>%
  bind_rows() %>%
  full_join(
    pyears_data %>% filter(period >= year.start) %>%
      select(country, mm1, period) %>%
      rename(sex = mm1,
             year = period,
             name = country) %>%
      distinct() %>%
      mutate(est = 1,
             sex = str_replace(sex, "female", "f"),
             sex = str_replace(sex, "male", "m"))
  ) %>%
  mutate(est = ifelse(is.na(est), 0, est),
         sex = factor(sex))

GBD.q.35 <- GBD.age %>%
  filter(age >= 15,
         age < 50) %>%
  group_by(name, sex, year) %>%
  summarise_at(vars(val), function(i){1 - prod(1 - i)}) %>%
  mutate(sex = str_replace(sex, "female", "f"),
         sex = str_replace(sex, "male", "m"))

WPP.q.35 <- WPP.age %>%
  filter(age >= 15,
         age < 50) %>%
  group_by(name, sex, year) %>%
  summarise_at(vars(WPP.qx), function(i){1 - prod(1 - i)}) %>%
  mutate(sex = str_replace(sex, "female", "f"),
         sex = str_replace(sex, "male", "m"))

df.q.35.all <- q3515.df %>%
  full_join(GBD.q.35) %>%
  full_join(WPP.q.35) %>%
  full_join(unpdhhd_qx %>% select(name, sex, year, q3515))

df.q.35.melt <- df.q.35.all %>%
  pivot_longer(cols = val:q3515,
               names_to = "variable",
               values_to = "q") %>%
  mutate(region = NA,
         region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
         region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
         region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
         region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
         region = fct_relevel(region, "Southern", "Eastern", "Central", "Western"),
         variable = str_replace(variable, "val", "GBD"),
         variable = str_replace(variable, "WPP.qx", "WPP"),
         variable = str_replace(variable, "q3515", "Household deaths"),
         variable = fct_relevel(variable, "WPP", "GBD", "Household deaths"),
         sex = str_replace(sex, "m", "male"),
         sex = str_replace(sex, "f", "female"),
         ) 

ggplot(df.q.35.melt %>% filter(variable == "GBD", est == 1)) + 
  geom_point(aes(x = value, y = q, color = sex), size = 2, stroke = 1.5, shape = 1) + 
  geom_abline() + theme_bw() +
  scale_color_manual(name = "", values = c("male" = m.col, "female" = f.col), guide = "none") +
  theme(axis.title = element_text(size = 30), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.background = element_blank(),
        aspect.ratio = 1,
        plot.margin = margin(0, -4, 0, 0, "cm")) +
  xlab("SH Empirical") + ylab(bquote("GBD / WPP")) +
  ggtitle("GBD") +
  scale_y_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  facet_wrap( ~ region) -> q3515.m.plot

ggplot(df.q.35.melt %>% filter(variable == "GBD", est == 1)) + 
  geom_point(aes(x = value, y = q, color = sex), size = 2, stroke = 1.5, shape = 1) + 
  geom_abline() + theme_bw() +
  scale_color_manual(name = "", values = c("male" = m.col, "female" = f.col), guide = "none") +
  theme(axis.title = element_text(size = 30), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1,
        plot.margin = margin(0, 0, 0, -4, "cm")) +
  xlab("SH Empirical") +  ylab("") +
  ggtitle("WPP") +
  scale_y_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  facet_wrap( ~ region) -> q3515.f.plot

cowplot::get_legend(ggplot(df.q.35.melt %>% filter(variable == "WPP", est == 1) %>%
                             mutate(sex = fct_relevel(sex, "male"))) + 
                      geom_point(aes(x = q, y = value, color = sex), size = 5, stroke = 1.5) + 
                      scale_color_manual(name = "", values = c("male" = m.col, "female" = f.col), labels = c("Males", "Females")) +
                      theme(legend.position = "bottom",
                            legend.background = element_blank(),
                            legend.key = element_blank(),
                            legend.text = element_text(size = 25),
                            aspect.ratio = 1) +
                      facet_wrap( ~ region)) -> q3515.legend

plot_grid(q3515.m.plot, q3515.f.plot, nrow = 1) -> btm.row

plot_grid(NULL, btm.row, nrow = 2, rel_heights = c(0.1, 1)) +
  draw_plot_label(c("Sibling~Histories~Empirical~vs~GBD/WPP~\"\"[35]*q[15]", ""), parse = TRUE,
                  size = 40, family = "sans", x = -0.07, y = 1.03, fontface = "bold") -> q3515.toprow

plot_grid(q3515.toprow, q3515.legend, nrow = 2, rel_heights = c(1, 0.1)) -> q3515.both.plot

ggsave(filename="35q15 compare m and f region.png",device="png",width=15.4686,height=8.898466667,unit="cm", scale = 3.15, plot= q3515.both.plot)

#q4515 compare####
GBD.q.45 <- GBD.age %>%
  filter(age >= 15,
         age < 60) %>%
  group_by(name, sex, year) %>%
  summarise_at(vars(val), function(i){1 - prod(1 - i)}) %>%
  mutate(sex = str_replace(sex, "female", "f"),
         sex = str_replace(sex, "male", "m"))

WPP.q.45 <- WPP.age %>%
  filter(age >= 15,
         age < 60) %>%
  group_by(name, sex, year) %>%
  summarise_at(vars(WPP.qx), function(i){1 - prod(1 - i)}) %>%
  mutate(sex = str_replace(sex, "female", "f"),
         sex = str_replace(sex, "male", "m"))

df.q.45.all <- q4515.df %>%
  full_join(GBD.q.45) %>%
  full_join(WPP.q.45) %>%
  full_join(unpdhhd_qx %>% select(name, sex, year, q4515))

df.q.45.melt <- df.q.45.all %>%
  pivot_longer(cols = val:q4515,
               names_to = "variable",
               values_to = "q") %>%
  mutate(region = NA,
         region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
         region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
         region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
         region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
         region = fct_relevel(region, "Southern", "Eastern", "Central", "Western"),
         variable = str_replace(variable, "val", "GBD"),
         variable = str_replace(variable, "WPP.qx", "WPP"),
         variable = str_replace(variable, "q4515", "Household deaths"),
         variable = fct_relevel(variable, "WPP", "GBD", "Household deaths"),
         sex = str_replace(sex, "m", "male"),
         sex = str_replace(sex, "f", "female")
         ) 

ggplot(df.q.45.melt %>% filter(variable == "GBD", est == 1)) + 
  geom_point(aes(x = value, y = q, color = sex), size = 2, stroke = 1.5, shape = 1) + 
  geom_abline() + theme_bw() +
  scale_color_manual(name = "", values = c("male" = m.col, "female" = f.col), guide = "none") +
  theme(axis.title = element_text(size = 30), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        legend.position = "bottom",
        legend.background = element_blank(),
        aspect.ratio = 1,
        plot.margin = margin(0, -4, 0, 0, "cm")) +
  xlab("SH Empirical") + ylab(bquote("GBD / WPP")) +
  ggtitle("GBD") +
  scale_y_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  facet_wrap( ~ region) -> q4515.m.plot

ggplot(df.q.45.melt %>% filter(variable == "GBD", est == 1)) + 
  geom_point(aes(x = value, y = q, color = sex), size = 2, stroke = 1.5, shape = 1) + 
  geom_abline() + theme_bw() +
  scale_color_manual(name = "", values = c("male" = m.col, "female" = f.col), guide = "none") +
  theme(axis.title = element_text(size = 30), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1,
        plot.margin = margin(0, 0, 0, -4, "cm")) +
  xlab("SH Empirical") +  ylab("") +
  ggtitle("WPP") +
  scale_y_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 5, labels = c(0, 0.25, "0.50", 0.75, 1)) +
  facet_wrap( ~ region) -> q4515.f.plot

cowplot::get_legend(ggplot(df.q.45.melt %>% filter(variable == "WPP", est == 1) %>%
                             mutate(sex = fct_relevel(sex, "male"))) + 
                      geom_point(aes(x = value, y = q, color = sex), size = 5, stroke = 1.5) + 
                      scale_color_manual(name = "", values = c("male" = m.col, "female" = f.col), labels = c("Males", "Females")) +
                      theme(legend.position = "bottom",
                            legend.background = element_blank(),
                            legend.key = element_blank(),
                            legend.text = element_text(size = 25),
                            aspect.ratio = 1) +
                      facet_wrap( ~ region)) -> q4515.legend

plot_grid(q4515.m.plot, q4515.f.plot, nrow = 1) -> btm.row

plot_grid(NULL, btm.row, nrow = 2, rel_heights = c(0.1, 1)) +
  draw_plot_label(c("Sibling~Histories~Empirical~vs~GBD/WPP~\"\"[45]*q[15]", ""), parse = TRUE,
                  size = 40, family = "sans", x = -0.07, y = 1.03, fontface = "bold") -> q4515.toprow

plot_grid(q4515.toprow, q4515.legend, nrow = 2, rel_heights = c(1, 0.1)) -> q4515.both.plot

ggsave(filename="45q15 compare m and f region.png",device="png",width=15.4686,height=8.898466667,unit="cm", scale = 3.15, plot= q4515.both.plot)

#nqx compare####
q5x.df <- mort.df %>%
  mutate(age.5 = age %/% 5) %>%
  group_by(year, name, sex, age.5, est) %>%
  summarise_at(vars(value), function(i){1 - prod(exp(-exp(i)))}) %>%
  ungroup %>%
  mutate(age = age.5 * 5 + 2) %>%
  select(-age.5)

q5x.df.all <- q5x.df %>%
  left_join(GBD.age %>% select(name, sex, year, age_group_name, val, age) %>%
              mutate(sex = str_replace(sex, "female", "f"), 
                     sex = str_replace(sex, "male", "m"))) %>%
  left_join(WPP.age %>% select(name, sex, year, WPP.qx, age) %>%
              mutate(sex = str_replace(sex, "female", "f"), 
                     sex = str_replace(sex, "male", "m"))) %>%
  mutate(region = NA,
         region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
         region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
         region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
         region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
         region = fct_relevel(region, "Southern", "Eastern", "Central", "Western"))

q5x.df.all %>%
  mutate(GBD.ratio = val / value,
         age_group_name = gsub(" to ", "-", age_group_name),
         sex = fct_relevel(sex, "m")) %>%
  drop_na(GBD.ratio) %>%
  filter(age %in% 15:60,
         est == 1) %>%
  ggplot() + 
  geom_boxplot(aes(x = age_group_name, y = GBD.ratio, col = sex), lwd = 1) +
  scale_color_manual(values = c("m" = m.col, "f" = f.col), labels = c("Males", "Females")) +
  scale_y_continuous(trans = "log", limits = c(0.25, 4), breaks = c(0.5, 1, 2, 4)) +
  ggtitle("\nGBD") +
  #labs(tag = "A") +
  xlab("") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(legend.text = element_text(size=25), legend.title = element_text(size=25),
        plot.title=element_text(size = 30,hjust=0.5, face = "bold"),
        strip.text= element_text(size = 25, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 18), 
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.title = element_text(size = 25),
        legend.position="bottom",
        legend.background=element_blank(),
        plot.tag.position = c(0.02, 1.002),
        plot.tag = element_text(size = 50, face = "bold"), 
        plot.margin = margin(0, -5, 0, 0, "cm"),
        panel.grid = element_blank(),
        aspect.ratio = 1.1,
        panel.spacing.x = unit(0.5, "cm")) +
  ylab("Ratio of age-specific mortality: GBD or WPP / Spline\n") +
  facet_wrap(~region) -> GBD.ratio.plot

get_legend(GBD.ratio.plot) -> ratio.legend

GBD.ratio.plot <- GBD.ratio.plot + theme(legend.position = "none")

q5x.df.all %>%
  mutate(WPP.ratio = WPP.qx / value,
         age_group_name = gsub(" to ", "-", age_group_name),
         sex = fct_relevel(sex, "m")) %>%
  drop_na(WPP.ratio) %>%
  filter(age %in% 15:60,
         est == 1) %>%
  ggplot() + 
  geom_boxplot(aes(x = age_group_name, y = WPP.ratio, col = sex), lwd = 1) +
  scale_y_continuous(trans = "log", limits = c(0.25, 4), breaks = c(0.5, 1, 2, 4)) +
  scale_color_manual(values = c("m" = m.col, "f" = f.col), guide = "none") +
  ggtitle("\nWPP") +
  #labs(tag = "B") +
  xlab("") + ylab("") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(legend.text = element_text(size=20), legend.title = element_text(size=25),
        plot.title=element_text(size=30,hjust=0.5, face = "bold"),
        strip.text= element_text(size=25, face = "bold"),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18, angle = 50, hjust = 1),
        axis.title = element_text(size = 25),
        legend.position="right", legend.box="vertical",
        legend.background=element_blank(),
        plot.tag.position = c(0.02, 1.002),
        plot.tag = element_text(size = 50, face = "bold"), 
        plot.margin = margin(0, 0, 0, -5, "cm"),
        panel.grid = element_blank(),
        aspect.ratio = 1.1,
        panel.spacing.x = unit(0.5, "cm")) +
  #ylab("Ratio of age-specific mortality: WPP / Spline \n") +
  facet_wrap(~region) -> WPP.ratio.plot

plot_grid(GBD.ratio.plot, WPP.ratio.plot, labels = c("A", "B"), label_size = 45, label_y = 1.01, nrow = 1) -> ratio.top.row
plot_grid(ratio.top.row, ratio.legend, rel_heights = c(2, 0.1), nrow = 2) -> ratio.both.plot
ggsave(filename="GBD WPP ratio.png",device="png",width=15.4686,height=8.898466667,unit="cm", scale = 3.15, plot= ratio.both.plot)


q5x.df.all %>%
  mutate(GBD.ratio = val / value,
         WPP.ratio = WPP.qx / value,
         age_group_name = gsub(" to ", "-", age_group_name)) %>%
  drop_na(GBD.ratio) %>%
  filter(age %in% 15:60,
         est == 1) %>%
  select(-est, -val, -WPP.qx, -value) -> ratio.df

write.csv(ratio.df, "GBD and WPP ratio.csv", row.names = F)

#plot Age patterns####
plot.cross.section.age.pattern <- function(country.set, y_scale = "identity", region){
  ggplot(mort.df.all %>%
           filter(year %in% c(1990,1995,2000,2005,2010,2015),
                  age%in%15:60,
                  country %in% country.set)) +
    geom_line(aes(x = age, y = mle, color = sex, linetype = factor(est)), size = 1.2) + 
    theme_bw() + 
    scale_y_continuous(trans = y_scale, labels = function(i){signif(i, digits = 3)}, position = "right",
                       n.breaks = 4) +
    geom_ribbon(aes(x = age, ymin = X2.5., ymax = X97.5., fill = sex), alpha = 0.3) +
    scale_linetype_manual(values = c("dashed","solid"), guide = "none") +
    scale_color_manual(values = c(f.col, m.col), label = c("Females", "Males")) +
    scale_fill_manual(values = c(f.fill.col, m.fill.col), label = c("Females", "Males"), guide = "none") +
    theme(legend.key.width = unit(1.5,"cm"),
          legend.text = element_text(size=25),
          title = element_text(size=30),
          plot.title = element_text(size = 40, face = "bold"),
          strip.text.x = element_text(size=25),
          strip.text.y = element_text(size=25, face = "bold"),
          strip.background.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          panel.spacing = unit(0, "cm"),
          strip.background = element_blank(),
          legend.position="bottom",legend.background = element_blank(),legend.box="vertical",
          panel.grid = element_blank()) +
    guides(col = guide_legend(order = 1,reverse=T,title="Sex")) +
    xlab("Age") + ylab("Mortality Rates") +
    ggtitle(region) +
    #ggtitle(paste("Estimated mortality schedules in selected years",region)) +
    facet_grid(country~year,switch = "y", scales="free")
}

plot.cross.section.time.pattern <- function(country.set, y_scale = "identity", region){
  age.lab<-c("Age 15","Age 25","Age 35","Age 45","Age 55")
  names(age.lab)<-c(seq(15,60,by=10))
  ggplot(subset(df.age.pattern,age%in%seq(15,60,by=10)&year%in%1990:2017 & country %in% country.set &est==1))+geom_line(aes(x=year,y=mle,color=sex),size=1.2)+theme_bw()+scale_y_continuous(trans=y_scale,labels=NULL)+
    geom_ribbon(aes(x=year,ymin=X2.5.,ymax=X97.5.,fill=sex),alpha=0.3)+
    scale_linetype_manual(values=c("dashed","solid"),guide=FALSE)+
    scale_color_manual(values=c("red","blue"),label=c("Females","Males"))+
    scale_fill_manual(values=c("red","blue"),label=c("Females","Males"),guide=FALSE)+
    theme(legend.key.width = unit(1.5,"cm"),legend.text = element_text(size=15),title = element_text(size=30),plot.title=element_text(size=40,hjust=0.5,face="bold"),
          strip.text.x = element_text(size=25),strip.text.y = element_text(size=15),
          axis.ticks.y=element_blank(),axis.text.x=element_text(size=15,face="bold"),legend.position="right",legend.background=element_blank(),legend.box="vertical")+
    guides(col = guide_legend(order = 1,reverse=T,title="Sex"))+
    xlab("Year")+ylab("Mortality Rates")+
    ggtitle(paste("Estimated mortality schedules at selected ages",region))+
    facet_grid(country~age,labeller=labeller(age=age.lab),switch="y",scales="free")
}

plot.cross.section.age.pattern(country.set=c("South Africa","Namibia","Zimbabwe","Zambia"),region="on log scale \n(Southern)",y_scale="log")

ggsave(filename="Age patterns southern.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("South Africa","Namibia","Zimbabwe","Zambia"),region="(Southern)"))
ggsave(filename="Log Age patterns southern.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("South Africa","Namibia","Zimbabwe","Zambia"),region="A",y_scale="log"))
ggsave(filename="Age patterns easter.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("Kenya","Uganda","Malawi","Tanzania"),region="(Eastern)"))
ggsave(filename="Log Age patterns eastern.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("Kenya","Uganda","Malawi","Tanzania"),region="B",y_scale="log"))
ggsave(filename="Age patterns central.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("Congo","Cameroon","Gabon","Congo Democratic\nRepublic"),region="(Central)"))
ggsave(filename="Log Age patterns central.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("Congo","Cameroon","Gabon","Congo Democratic\nRepublic"),region="C",y_scale="log"))
ggsave(filename="Age patterns western.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("Sierra Leone","Nigeria","Cote d'Ivoire","Benin"),region="(Western)"))
ggsave(filename="Log Age patterns western.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,plot=plot.cross.section.age.pattern(country.set=c("Sierra Leone","Nigeria","Cote d'Ivoire","Benin"),region="D",y_scale="log"))

#plot data####
pyears_data %>%
  filter(period >= year.start) %>%
  group_by(country, survey_type) %>%
  summarise_at(vars(survey_year), unique) %>%
  mutate(survey_year = ifelse(survey_type != "DHS", paste0(survey_year, " (", survey_type, ")"), as.character(survey_year))) %>%
  group_by(country) %>%
  summarise(Survey = paste(survey_year, collapse = ", "),
            n = n()) -> dat.dhs.df
write.csv(dat.dhs.df, "DHS table.csv", row.names = F)

fig.country <- pyears_data %>% 
  filter(period >= year.start) %>%
  rename(name = country) %>%
  mutate(region = NA,
         region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
         region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
         region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
         region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
         region = fct_relevel(region, "Southern", "Eastern", "Central", "Western")) %>%
  arrange(region, name) %>%
  .$name %>% unique
  
dat.dhs.fig.df.both <- pyears_data %>%
  rename(name = country) %>%
  group_by(name, survey_id, survey_type, survey_year) %>%
  summarise_at(vars(pyears2, adjusted, n), sum) %>%
  ungroup %>%
  left_join(tibble(name = fig.country, order = 1:length(joint.countries))) %>%
  mutate(region = NA,
         region = ifelse(name %in% c("Eswatini","Lesotho","Namibia","South Africa","Zimbabwe","Zambia"), "Southern", region),
         region = ifelse(name %in% c("Burundi","Kenya","Malawi","Mozambique","Rwanda","Tanzania","Uganda","Madagascar","Comoros","Ethiopia"), "Eastern", region),
         region = ifelse(name %in% c("Angola","Cameroon","Chad","Congo Democratic Republic","Congo","Gabon","Central African Republic","Sao Tome and Principe"), "Central", region),
         region = ifelse(name %in% c("Benin","Burkina Faso","Cote d'Ivoire","Gambia","Guinea","Liberia","Mali","Niger","Nigeria","Senegal","Sierra Leone","Togo","Mauritania","Guinea-Bissau"), "Western", region),
         order = ifelse(region == "Eastern", order + 1, order),
         order = ifelse(region == "Central", order + 2, order),
         order = ifelse(region == "Western", order + 3, order)) %>%
  full_join(tibble(order = 1:(length(joint.countries) + 3)), by = "order") %>%
  mutate(pyears2 = ifelse(is.na(pyears2), 0, pyears2),
         survey_id = ifelse(is.na(survey_id), 2015, as.character(survey_id)),
         survey_type = ifelse(is.na(survey_type), "DHS", as.character(survey_type)),
         n = ifelse(is.na(n), 0, n),
         region = ifelse(is.na(region), "Southern", as.character(region)),
         region = fct_relevel(region, "Southern", "Eastern", "Central"))

n.countries <- dat.dhs.fig.df.both %>%
  drop_na(name) %>%
  group_by(region) %>%
  summarise_at(vars(name), n_distinct) %>%
  mutate(cumsum = cumsum(name))
  
dat.dhs.fig.df.both %>%
  ggplot(aes(x = survey_year, y = factor(order), color = region, size = n, pch = survey_type)) +
  annotate("text", 1990, c(3.5, 13, 22.5, 34.5), label = c("Southern", "Eastern", "Central", "Western"), angle = 90, fontface = "bold", size = 10) +
  geom_hline(yintercept = dat.dhs.fig.df.both %>% filter(n == 0) %>% .$order, size = 1.2) +
  geom_point(stroke = 1) +
  scale_x_continuous("survey year", breaks = 1992:2019, minor_breaks = NULL) +
  scale_y_discrete(element_blank(), breaks = 1:(length(joint.countries) + 3),
                   labels = c(as.character(fig.country[1:n.countries$cumsum[1]]), " ",
                              as.character(fig.country[(n.countries$cumsum[1] + 1):n.countries$cumsum[2]]), " ", 
                              as.character(fig.country[(n.countries$cumsum[2] + 1):n.countries$cumsum[3]]), " ", 
                              as.character(fig.country[(n.countries$cumsum[3] + 1):n.countries$cumsum[4]])),
                   position = "right", expand = expansion(add = 0.6)) +
  scale_size_area("Sample Size", max_size = 10,  c(2e6, 4e6, 6e6)) +
  scale_color_nejm(guide = "none") +
  scale_shape_manual(values = c("DHS" = 16, "MICS" = 17, "AIS" = 15), name = "Survey Type",
                     guide = guide_legend(override.aes = list(size = 3))) + 
  theme_light(10) +
  guides(#color = guide_legend(order = 1, override.aes = list(size = 5)),
    #shape = guide_legend(order = 2),
    size = guide_legend(order = 3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing = unit(0, "cm"),
        legend.key.size = unit(0, "pt"),
        legend.text.align=0,
        plot.margin = margin(t = 1, l = 50, unit = "pt")) +
  coord_cartesian(xlim = c(1992, 2019), clip = "off")

aaa <- pyears_data %>% group_by(country) %>%
  summarise(year = length(unique(survey_year))) %>%
  filter(year > 1) %>% .$country %>% .[!. %in% c("Togo", "South Africa")]

pyears_data %>% filter(country %in% aaa) %>% .$survey_id %>% unique  %>% length
