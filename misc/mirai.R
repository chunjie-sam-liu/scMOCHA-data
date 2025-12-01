library(tibble)
library(dplyr)
library(mirai)

daemons(5)
df <- tibble(
  id = 1:10000,
  x = rnorm(10000),
  y = runif(10000)
)
df |>
  dplyr::mutate(
    a = mirai::mirai_map(
      .x = x,
      .f = function(.x, .y) {
        Sys.sleep(0.01)
        .x + .y
      },
      .y = y,
    )
  )
daemons(0)
