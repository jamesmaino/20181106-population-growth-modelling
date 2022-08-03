library(tidyverse)

trans <- matrix(
    c(
        0.00, 0.50, 0.00,
        0.00, 0.00, 0.50,
        4.00, 0.00, 0.10
    ),
    ncol = 3, byrow = T
)

rownames(trans) <- c("V1", "V2", "V3")

colnames(trans) <- c("V1", "V2", "V3")

Tmax <- 200

N <- matrix(rep(0, Tmax * 3), nrow = 3)

N[1, 1] <- 1000

for (t in 2:Tmax) N[, t] <- N[, t - 1] %*% trans

ages <- as.data.frame(t(N))

ages$total <- rowSums(ages)

ages$t <- 1:Tmax

agesl <- gather(ages, age.class, value, -t, -total)

ggplot(agesl) +
    geom_area(aes(t, value / total, fill = age.class)) +
    ylab("proportion") +
    theme_minimal() +
    theme(
        plot.background = element_rect(fill = rgb(0.2, 0.21, 0.27)),
        text = element_text(colour = "grey"),
        axis.text = element_text(colour = "grey"),
        panel.grid = element_line(colour = "grey")
    )

ggsave('./plots/plot1.png', width=10, height=6)

df <- data.frame(
    t = 0:(Tmax - 1),
    N = colSums(N)
)

lm1 <- lm(log(N) ~ t, data = df)

coef(lm1)

# plot
df$pred <- predict(lm1)

ggplot(df) +
    geom_line(aes(t, N), colour = "red", alpha = 0.5) +
    geom_line(aes(t, exp(pred)), colour = "grey") +
    theme_minimal() +
    theme(
        plot.background = element_rect(fill = rgb(.2, .21, .27)),
        text = element_text(colour = "grey"),
        axis.text = element_text(colour = "grey"),
        panel.grid = element_line(colour = "grey")
    ) +
    scale_y_log10()

ggsave('./plots/plot2.png', width=10, height=6)

N1_stable <- 1000 * N[, Tmax] / sum(N[, Tmax])

N <- matrix(rep(0, Tmax * 3), nrow = 3)

N[, 1] <- N1_stable

for (t in 2:Tmax) N[, t] <- N[, t - 1] %*% trans

df <- data.frame(
    t = 0:(Tmax - 1),
    N = colSums(N)
)

lm1 <- lm(log(N) ~ t, data = df)

coef(lm1)

exp(coef(lm1)[1])

mite <- read_csv("./data/AtlanticSpiderMite_CareyBradley1982.csv")

mite

# plot

mite$`P(x).m(x)` <- mite$`P(x)` * mite$`m(x)`

dw <- mite %>% gather(variable, value, -x)

ggplot(dw) +
    geom_point(aes(x, value, colour = variable), size = 2) +
    geom_line(aes(x, value, colour = variable)) +
    theme_minimal() +
    theme(
        plot.background = element_rect(fill = rgb(.2, .21, .27)),
        text = element_text(colour = "grey"),
        axis.text = element_text(colour = "grey"),
        panel.grid = element_line(colour = "grey")
    )

ggsave('./plots/plot3.png', width=10, height=6)
