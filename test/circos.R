set.seed(999)
n = 1000
df = data.frame(
  sectors = sample(letters[1:8], n, replace = TRUE),
  x = rnorm(n), y = runif(n)
)
head(df)
library(circlize)
circos.par("track.height" = 0.1)
circos.initialize(df$sectors, x = df$x)
circos.par

circos.track(df$sectors,
  y = df$y,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$cell.ylim[2] + mm_y(5),
      CELL_META$sector.index
    )
    circos.axis(labels.cex = 0.6)
  }
)
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
circos.clear()



library(yaml)
data = yaml.load_file("https://raw.githubusercontent.com/Templarian/slack-emoji-pokemon/master/pokemon.yaml")
pokemon_list = data$emojis[1:40]
pokemon_name = sapply(pokemon_list, function(x) x$name)
pokemon_src = sapply(pokemon_list, function(x) x$src)

library(EBImage)
circos.par("points.overflow.warning" = FALSE)
circos.initialize(pokemon_name, xlim = c(0, 1))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  image = EBImage::readImage(pokemon_src[CELL_META$sector.numeric.index])
  print(CELL_META)
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(2),
    CELL_META$sector.index,
    facing = "clockwise", niceFacing = TRUE,
    adj = c(1, 0.5), cex = 0.6
  )
  rasterImage(image,
    xleft = pos[1, 1] - 0.05, ybottom = pos[1, 2] - 0.05,
    xright = pos[1, 1] + 0.05, ytop = pos[1, 2] + 0.05
  )
}, bg.border = 1, track.height = 0.15)


circos.clear()


sectors = letters[1:10]
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
circos.initialize(sectors, xlim = cbind(rep(0, 10), runif(10, 0.5, 1.5)))
circos.track(
  ylim = c(0, 1), track.height = mm_h(5),
  panel.fun = function(x, y) {
    circos.lines(c(0, 0 + mm_x(5)), c(0.5, 0.5), col = "blue")
  }
)
circos.track(
  ylim = c(0, 1), track.height = cm_h(1),
  track.margin = c(0, mm_h(2)),
  panel.fun = function(x, y) {
    xcenter = get.cell.meta.data("xcenter")
    circos.lines(c(xcenter, xcenter), c(0, cm_y(1)), col = "red")
  }
)
circos.track(
  ylim = c(0, 1), track.height = inch_h(1),
  track.margin = c(0, mm_h(5)),
  panel.fun = function(x, y) {
    line_length_on_x = cm_x(1 * sqrt(2) / 2)
    line_length_on_y = cm_y(1 * sqrt(2) / 2)
    circos.lines(c(0, line_length_on_x), c(0, line_length_on_y), col = "orange")
  }
)

circos.clear()


par(mar = c(1, 1, 1, 1))
plot(NULL, xlim = c(1, 10), ylim = c(1, 8), axes = FALSE, ann = FALSE)
points(1:10, rep(1, 10),
  pch = 16, cex = 5,
  col = rand_color(10, luminosity = "random")
)
points(1:10, rep(2, 10),
  pch = 16, cex = 5,
  col = rand_color(10, luminosity = "bright")
)
points(1:10, rep(3, 10),
  pch = 16, cex = 5,
  col = rand_color(10, luminosity = "light")
)
points(1:10, rep(4, 10),
  pch = 16, cex = 5,
  col = rand_color(10, luminosity = "dark")
)
points(1:10, rep(5, 10),
  pch = 16, cex = 5,
  col = rand_color(10, hue = "red", luminosity = "bright")
)
points(1:10, rep(6, 10),
  pch = 16, cex = 5,
  col = rand_color(10, hue = "green", luminosity = "bright")
)
points(1:10, rep(7, 10),
  pch = 16, cex = 5,
  col = rand_color(10, hue = "blue", luminosity = "bright")
)
points(1:10, rep(8, 10),
  pch = 16, cex = 5,
  col = rand_color(10, hue = "monochrome", luminosity = "bright")
)


category = paste0("category", "_", 1:9)
percent = sort(sample(40:80, 9))
color = rev(rainbow(length(percent)))

library(circlize)
circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
circos.initialize("a", xlim = c(0, 100)) # 'a` just means there is one sector
circos.track(
  ylim = c(0.5, length(percent) + 0.5), track.height = 0.8,
  bg.border = NA, panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    circos.segments(rep(xlim[1], 9), 1:9,
      rep(xlim[2], 9), 1:9,
      col = "#CCCCCC"
    )
    circos.rect(rep(0, 9), 1:9 - 0.45, percent, 1:9 + 0.45,
      col = color, border = "white"
    )
    circos.text(rep(xlim[1], 9), 1:9,
      paste(category, " - ", percent, "%"),
      facing = "downward", adj = c(1.05, 0.5), cex = 0.8
    )
    breaks = seq(0, 85, by = 5)
    circos.axis(
      h = "top", major.at = breaks, labels = paste0(breaks, "%"),
      labels.cex = 0.6
    )
  }
)

circos.clear()
