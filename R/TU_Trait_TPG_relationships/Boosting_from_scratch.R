# ______________________________________
# Boosted regression tree from scratch
# ______________________________________
library(data.table)
library(ggplot2)
library(rpart)

# Create some noisy continuous data
set.seed(1234)
fun <- function(x, offset) {
    x * sin(2 * x) + offset
}

dt <- data.table(
    x = seq(0, 10, 0.1),
    y = fun(
        x = seq(0, 10, 0.1),
        offset = 0.5
    )
)
noise <- data.table(
    x = seq(0, 10, 0.01),
    y_noise = fun(
        x = seq(0, 10, 0.01),
        offset = 0.5
    ) + rnorm(1001, 0, .75)
)

ggplot(dt, aes(x = x, y = y)) +
    geom_line() +
    geom_point(
        data = noise,
        aes(x = x, y = y_noise),
        alpha = 0.2
    )

# Boosting Algo: 
target <- noise$y_noise
est_correction <- 0
learning_rate <- 0.2
steps <- 100
resid <- vector(mode = "numeric")
pred <- 0

for (i in 1:steps) {
    correction <- target - learning_rate * est_correction

    regres <- rpart(correction ~ noise$x,
        method = "anova"
    )
    plot(regres)
    text(regres, use.n = TRUE)

    est_correction <- predict(object = regres)

    # Update prediction by adding contirbution from tree
    pred <- pred + learning_rate * est_correction

    # error
    resid[[i]] <- mean((pred - noise$y_noise)^2)

    target <- correction
}
bst_re <- data.table(
    x = seq(0, 10, 0.01),
    pred = pred
)

ggplot(dt, aes(x = x, y = y)) +
    geom_line() +
    geom_point(
        data = noise,
        aes(x = x, y = y_noise),
        alpha = 0.2
    ) +
    geom_line(
        data = bst_re,
        aes(y = pred),
        color = "red"
    )
