
library("odin") 
library(tidyverse)
library(patchwork)

seed_val <- 123
set.seed(seed_val)

# mosquito parameters
kappa <- 0.8     # biting rate
omega <- 0.025   # larval mosquito mortality rate
epsilon <- 0.1   # adult mosquito mortality rate
alpha <- 1 / 19  # one over duration of larval stage
eta <- 1 / 10    # one over duration of incubation period
M_eq <- 1e5      # target number of adult mosquitoes at equilibrium
delta_k <- 0#0.4   # magnitude of fluctuation in carrying capacity
R_mosq <- 2.69   # mosquito reproduction number (average offspring per mosquito in unconstrained environment)

# human parameters
N <- 5e4                # total human population size
p_death <-  0.0136 / 365 # death rate
p_birth <- p_death      # birth rate (set equal to death rate for constant population size)
p_recov <- 1 / 365      # one over duration of cross-protection
extInf <- 1e-6          # additional foi on humans due to external sources
beta_hm <- 1            # per bite transmission probability from humans to mosquitoes
incub <- 6              # duration of intrinsic incubation period (time before a human becomes infectious)
inf_period <- 4         # duration of infectious period in humans
beta_mh_1 <- 0.2          # per bite transmission probability from mosquitoes to humans (serotype 1)
beta_mh_2 <- 0.2          # per bite transmission probability from mosquitoes to humans (serotype 2)
beta_mh_3 <- 0.2          # per bite transmission probability from mosquitoes to humans (serotype 3)
beta_mh_4 <- 0.2          # per bite transmission probability from mosquitoes to humans (serotype 4)

# timing and outbreak parameters
t_burnin <- 100*365
tt <- seq(1, t_burnin + 10*365, by = 1)

# rate at which adult female mosquitoes produce larvae. Calculated based on
# mosquito reproduction number R_mosq
b <- epsilon*R_mosq*(alpha + omega) / alpha

# larval carrying capacity. Calculated based on target equilibrium value M_eq
K0 <- M_eq*epsilon*omega / (alpha^2*b/epsilon - alpha^2 - alpha*omega)

# load and define the odin model
generator <- odin::odin("denv_model.R")
model <- generator$new(kappa = kappa, epsilon = epsilon, omega = omega, alpha = alpha, b = b, eta = eta, 
                       K0 = K0, delta_k = delta_k, N = N, p_death = p_death, p_birth = p_birth, 
                       p_recov = p_recov, extInf = extInf, beta_hm = beta_hm, incub = incub,
                       inf_period = inf_period, beta_mh_1 = beta_mh_1, beta_mh_2 = beta_mh_2,
                       beta_mh_3 = beta_mh_3, beta_mh_4 = beta_mh_4)

# run model and format output
y <- model$transform_variables(model$run(tt)) |>
  as.data.frame() |>
  mutate(t = t - t_burnin) |>
  filter(t > 0)

y_weekly <- y |>
  mutate(week = floor(t/7)) |>  # group by week
  mutate(inc_total = inc_inf1 + inc_inf2 + inc_inf3 + inc_inf4,
         prop1 = inc_inf1 / inc_total,
         prop2 = inc_inf2 / inc_total,
         prop3 = inc_inf3 / inc_total,
         prop4 = inc_inf4 / inc_total) |>
  group_by(week) |>
  summarise(weekly_total=sum(inc_total),
            weekly_inf1=sum(inc_inf1),
            weekly_inf2=sum(inc_inf2),
            weekly_inf3=sum(inc_inf3),
            weekly_inf4=sum(inc_inf4),
            weekly_p1=weekly_inf1 / weekly_total,
            weekly_p2=weekly_inf2 / weekly_total,
            weekly_p3=weekly_inf3 / weekly_total,
            weekly_p4=weekly_inf4 / weekly_total,
            weekly_y1=sum(Y1),
            weekly_y2=sum(Y2),
            weekly_y3=sum(Y3),
            weekly_y4=sum(Y4),
            weekly_rt1=mean(Rt_1),
            weekly_rt2=mean(Rt_2),
            weekly_rt3=mean(Rt_3),
            weekly_rt4=mean(Rt_4),
            mean_rt=mean(c(Rt_1, Rt_2, Rt_3, Rt_4)))

y_weekly <- y_weekly |> 
  mutate(date = seq(from = as.Date("2000-01-01"), by = "week", 
                    length.out = dim(y_weekly)[1]))

y_monthly <- y_weekly |>
  #filter((t %% 7) == 0) |> # reporting frequency
  group_by(month = format(date, "%Y-%m")) |>
  summarise(monthly_total=sum(weekly_total),
            monthly_inf1=sum(weekly_inf1),
            monthly_inf2=sum(weekly_inf2),
            monthly_inf3=sum(weekly_inf3),
            monthly_inf4=sum(weekly_inf4),
            monthly_p1=monthly_inf1 / monthly_total,
            monthly_p2=monthly_inf2 / monthly_total,
            monthly_p3=monthly_inf3 / monthly_total,
            monthly_p4=monthly_inf4 / monthly_total, 
            date = first(date))

#---------------------------
## Plotting

# plot carrying capacity and total adult mosquito population size
y |>
  select(t, K, M) |>
  pivot_longer(cols = -1, names_to = "parameter") |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t / 365, y = value, color = parameter)) +
  geom_hline(yintercept = M_eq, linetype = "dashed") +
  ylim(c(0, NA)) + xlab("Time (years)") + ylab("Value") +
  ggtitle("Mosquito dynamics")

# plot incidence of human infection
y |>
  select(t, sprintf("inc_inf%s", 1:4)) |>
  #mutate(inc_total = inc_inf1 + inc_inf2 + inc_inf3 + inc_inf4) |> # optionally plot total incidence
  pivot_longer(cols = -1, names_to = "parameter") |>
  ggplot() + theme_bw() +
  geom_line(aes(x = t / 365, y = value, color = parameter)) +
  ylim(c(0, NA)) + xlab("Time (years)") + ylab("Value") +
  ggtitle("Incidence of infection")

y_weekly |>
  #filter((t %% 7) == 0) |> # reporting frequency
  select(week, sprintf("weekly_p%s", 1:4)) |>
  rename("Serotype 1" = "weekly_p1",
         "Serotype 2" = "weekly_p2",
         "Serotype 3" = "weekly_p3",
         "Serotype 4" = "weekly_p4") |>
  pivot_longer(cols = -1, names_to = "parameter") |>
  ggplot() + theme_bw(base_size=16) +
  theme(legend.position="inside",
        legend.position.inside = c(0.78,0.85),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = NA),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  geom_step(aes(x = week/52, y = value, color = parameter)) +
  #geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 1) +
  labs(color = NULL) +
  ylim(c(0, 1)) + xlab("Time (years)") + ylab("Serotype proportions") +
  scale_x_continuous(breaks = seq(0, max(y_weekly$week / 52), by = 1))

# plot Rt
y |>
  select(t, sprintf("Rt_%s", 1:4)) |>
  pivot_longer(cols = -1, names_to = "parameter") |>
  ggplot() + theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = t / 365, y = value, color = parameter)) +
  scale_y_continuous(limits = c(0, NA)) +
  ggtitle("Rt for each serotype")


p1 <- y_weekly |>
  select(week, sprintf("weekly_y%s", 1:4)) |>
  pivot_longer(cols = -1, names_to = "parameter") |>
  mutate(parameter = recode(parameter,
                            "weekly_y1" = "Serotype 1",
                            "weekly_y2" = "Serotype 2",
                            "weekly_y3" = "Serotype 3",
                            "weekly_y4" = "Serotype 4")) |>
  ggplot() + theme_bw(base_size=16) +
  theme(legend.position="inside",
        legend.position.inside = c(0.2,0.93),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = NA),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  labs(color = NULL) +
  xlab("Time (years)") + ylab("Number of infected mosquitoes") +
  geom_step(aes(x = week / 52, y = value, color = parameter)) +
  scale_x_continuous(breaks = seq(0, max(y_weekly$week / 52), by = 1)) +
  ggtitle("(a)")

p2 <- y_weekly |>
  select(week, weekly_total) |>
  ggplot() + theme_bw(base_size=16) +
  geom_step(aes(x = week/52, y = weekly_total)) +
  xlab("Time (years)") + ylab("Weekly total incidence") +
  scale_x_continuous(breaks = seq(0, max(y_weekly$week / 52), by = 1)) +
  ggtitle("(d)")

p3 <- y_weekly |>
  #filter((t %% 7) == 0) |> # reporting frequency
  select(week, sprintf("weekly_p%s", 1:4)) |>
  rename("Serotype 1" = "weekly_p1",
         "Serotype 2" = "weekly_p2",
         "Serotype 3" = "weekly_p3",
         "Serotype 4" = "weekly_p4") |>
  pivot_longer(cols = -1, names_to = "parameter") |>
  ggplot() + theme_bw(base_size=16) +
  theme(legend.position="inside",
        legend.position.inside = c(0.2,0.93),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = NA),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  geom_step(aes(x = week/52, y = value, color = parameter)) +
  #geom_vline(xintercept = 4, linetype = "dashed", color = "black", size = 1) +
  labs(color = NULL) +
  ylim(c(0, 1)) + xlab("Time (years)") + ylab("Serotype proportions") +
  scale_x_continuous(breaks = seq(0, max(y_weekly$week / 52), by = 1)) +
  ggtitle("(b)")

# plot Rt
#png("rt_daily_789_3.png", width = 4000, height = 2000, res = 300)
p4 <- y |>
  select(t, sprintf("Rt_%s", 1:4)) |>
  pivot_longer(cols = -1, names_to = "parameter") |>
  mutate(parameter = recode(parameter,
                            "Rt_1" = "Serotype 1",
                            "Rt_2" = "Serotype 2",
                            "Rt_3" = "Serotype 3",
                            "Rt_4" = "Serotype 4")) |>
  ggplot() + theme_bw(base_size = 16) +
  theme(legend.position="inside",
        legend.position.inside = c(0.8,0.85),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = NA),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  geom_step(aes(x = t / 364, y = value, color = parameter)) +
  xlab("Time (years)") + ylab("Value") + ylab(expression(R[t])) +
  labs(color = NULL) +
  scale_x_continuous(breaks = seq(0, max(y$t / 364), by = 1)) + 
  ggtitle("(c)")
#dev.off()

#png("model_plots_123_10.png", width = 4000, height = 6000, res = 300)
plot <- (p1 / p3 / p4 / p2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
plot
#dev.off()
