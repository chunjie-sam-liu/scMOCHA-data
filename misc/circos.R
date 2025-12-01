set.seed(999)
n = 1000
df = data.frame(
  sectors = sample(letters[1:8], n, replace = TRUE),
  x = rnorm(n), y = runif(n)
)
head(df)

library(circlize)
circos.par(
  "track.height" = 0.1,
  "start.degree" = 90,
  "gap.degree" = 5,
  "gap.after" = 1,
  "track.margin" = c(0.01, 0.01),
)
circos.initialize(df$sectors, x = df$x)
circos.par

circos.track(df$sectors,
  y = df$y,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$cell.ylim[2] + mm_y(5),
      CELL_META$sector.index,
      col = "red",
    )
    circos.axis(labels.cex = 0.6, col = "red")
  }
)
# circos.clear()

col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(df$sectors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)

circos.track(df$sectors,
  x = df$x, y = df$y,
  panel.fun = function(x, y) {
    ind = sample(length(x), 10)
    x2 = x[ind]
    y2 = y[ind]
    od = order(x2)
    circos.lines(x2[od], y2[od])
  }
)
circos.update(
  sector.index = "d", track.index = 2,
  bg.col = "#FF8080", bg.border = "black"
)
circos.points(x = -2:2, y = rep(0.5, 5), col = "white")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "updated", col = "white")


circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 0.1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
    breaks[-1], rep(ylim[2], n_breaks - 1),
    col = rand_color(n_breaks), border = NA
  )
})

circos.link("a", 0, "b", 0, h = 0.4)
circos.link("c", c(-0.5, 0.5), "d", c(-0.5, 0.5),
  col = "red",
  border = "blue", h = 0.2
)
circos.link("e", 0, "g", c(-1, 1), col = "green", border = "black", lwd = 2, lty = 2)

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

library(circlize)

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

library(circlize)

# rand colors
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




library(circlize)
category = paste0("category", "_", 1:9)
percent = sort(sample(40:80, 9))
color = rev(rainbow(length(percent)))
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


col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-5, 5, by = 1))
col2value(col_fun(seq(-5, 5, by = 1)), col_fun = col_fun)


library(circlize)
cell_cycle = data.frame(
  phase = factor(
    c("G1", "S", "G2", "M"),
    levels = c("G1", "S", "G2", "M")
  ),
  hour = c(11, 8, 4, 1)
)
color = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
circos.par(start.degree = 90)
circos.initialize(
  cell_cycle$phase,
  xlim = cbind(rep(0, 4), cell_cycle$hour)
)
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    circos.arrow(
      CELL_META$xlim[1],
      CELL_META$xlim[2],
      arrow.head.width = CELL_META$yrange * 0.8,
      arrow.head.length = cm_x(0.5),
      col = color[CELL_META$sector.numeric.index]
    )
    circos.text(
      CELL_META$xcenter,
      CELL_META$ycenter,
      CELL_META$sector.index,
      facing = "downward"
    )
    circos.axis(
      h = 1,
      major.at = seq(0, round(CELL_META$xlim[2])), minor.ticks = 1,
      labels.cex = 0.6
    )
  },
  bg.border = NA,
  track.height = 0.3
)
circos.clear()



library(png)
image = system.file("extdata", "Rlogo.png", package = "circlize")
image = as.raster(readPNG(image))
circos.par(start.degree = 90)
circos.initialize(letters[1:5], xlim = c(0, 1))
all_facing_options = c(
  "inside",
  "outside",
  "reverse.clockwise",
  "clockwise",
  "downward"
)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.raster(
    image,
    CELL_META$xcenter,
    CELL_META$ycenter,
    width = "1cm",
    facing = all_facing_options[CELL_META$sector.numeric.index]
  )
  circos.text(
    CELL_META$xcenter,
    CELL_META$ycenter,
    all_facing_options[CELL_META$sector.numeric.index],
    facing = "inside",
    niceFacing = TRUE
  )
})
circos.clear()


load(system.file("extdata", "doodle.RData", package = "circlize"))
circos.par("cell.padding" = c(0, 0, 0, 0))
circos.initialize(letters[1:16], xlim = c(0, 1))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  img = img_list[[CELL_META$sector.numeric.index]]
  circos.raster(img, CELL_META$xcenter, CELL_META$ycenter,
    width = CELL_META$xrange, height = CELL_META$yrange,
    facing = "bending.inside"
  )
}, track.height = 0.25, bg.border = NA)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  img = img_list[[CELL_META$sector.numeric.index + 16]]
  circos.raster(img, CELL_META$xcenter, CELL_META$ycenter,
    width = CELL_META$xrange, height = CELL_META$yrange,
    facing = "bending.inside"
  )
}, track.height = 0.25, bg.border = NA)
circos.clear()


library(circlize)

col_fun = colorRamp2(c(-2, 0, 2), c("green", "yellow", "red"))
circlize_plot = function() {
  set.seed(12345)
  sectors = letters[1:10]
  circos.initialize(sectors, xlim = c(0, 1))
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.points(runif(20), runif(20), cex = 0.5, pch = 16, col = 2)
    circos.points(runif(20), runif(20), cex = 0.5, pch = 16, col = 3)
  })
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.lines(sort(runif(20)), runif(20), col = 4)
    circos.lines(sort(runif(20)), runif(20), col = 5)
  })

  for (i in 1:10) {
    circos.link(sample(sectors, 1), sort(runif(10))[1:2],
      sample(sectors, 1), sort(runif(10))[1:2],
      col = add_transparency(col_fun(rnorm(1)))
    )
  }
  circos.clear()
}

library(ComplexHeatmap)
library(circlize)
# discrete
lgd_points = Legend(
  at = c("label1", "label2"), type = "points",
  legend_gp = gpar(col = 2:3), title_position = "topleft",
  title = "Track1"
)
# discrete
lgd_lines = Legend(
  at = c("label3", "label4"), type = "lines",
  legend_gp = gpar(col = 4:5, lwd = 2), title_position = "topleft",
  title = "Track2"
)
# continuous
lgd_links = Legend(
  at = c(-2, -1, 0, 1, 2), col_fun = col_fun,
  title_position = "topleft", title = "Links"
)

lgd_list_vertical = packLegend(lgd_points, lgd_lines, lgd_links)
lgd_list_vertical

circlize_plot()
draw(
  lgd_list_vertical,
  x = unit(4, "mm"),
  y = unit(4, "mm"),
  just = c("left", "bottom")
)
draw(lgd_links,
  x = unit(1, "npc") - unit(2, "mm"),
  y = unit(4, "mm"),
  just = c("right", "bottom")
)


set.seed(123)
mat1 = rbind(
  cbind(
    matrix(rnorm(50 * 5, mean = 1), nr = 50),
    matrix(rnorm(50 * 5, mean = -1), nr = 50)
  ),
  cbind(
    matrix(rnorm(50 * 5, mean = -1), nr = 50),
    matrix(rnorm(50 * 5, mean = 1), nr = 50)
  )
)
rownames(mat1) = paste0("R", 1:100)
colnames(mat1) = paste0("C", 1:10)
mat1 = mat1[sample(100, 100), ] # randomly permute rows
split = sample(letters[1:5], 100, replace = TRUE)
split = factor(split, levels = letters[1:5])


library(ComplexHeatmap)
Heatmap(mat1, row_split = split)


library(circlize) # >= 0.4.10
col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
circos.heatmap(mat1, split = split, col = col_fun1)
circos.clear()

circos.par(start.degree = 90, gap.degree = 10)
circos.heatmap(
  mat1,
  split = split,
  col = col_fun1,
  track.height = 0.5,
  bg.border = "green",
  bg.lwd = 2,
  bg.lty = 2,
  show.sector.labels = TRUE,
  dend.side = "outside",
  rownames.side = "inside"
)
circos.clear()


source("https://gist.githubusercontent.com/jokergoo/0ea5639ee25a7edae3871ed8252924a1/raw/57ca9426c2ed0cebcffd79db27a024033e5b8d52/random_matrices.R")
set.seed(123)
km = kmeans(mat_meth, centers = 5)$cluster
col_meth = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
circos.heatmap(mat_meth, split = km, col = col_meth, track.height = 0.12)

col_direction = c("hyper" = "red", "hypo" = "blue")
circos.heatmap(direction, col = col_direction, track.height = 0.01)

col_expr = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
circos.heatmap(mat_expr, col = col_expr, track.height = 0.12)

col_pvalue = colorRamp2(c(0, 2, 4), c("white", "white", "red"))
circos.heatmap(cor_pvalue, col = col_pvalue, track.height = 0.01)

library(RColorBrewer)
col_gene_type = structure(brewer.pal(length(unique(gene_type)), "Set3"), names = unique(gene_type))
circos.heatmap(gene_type, col = col_gene_type, track.height = 0.01)

col_anno_gene = structure(brewer.pal(length(unique(anno_gene)), "Set1"), names = unique(anno_gene))
circos.heatmap(anno_gene, col = col_anno_gene, track.height = 0.01)

col_dist = colorRamp2(c(0, 10000), c("black", "white"))
circos.heatmap(dist, col = col_dist, track.height = 0.01)

col_enhancer = colorRamp2(c(0, 1), c("white", "orange"))
circos.heatmap(anno_enhancer, col = col_enhancer, track.height = 0.03)

circos.clear()


circos.initializeWithIdeogram()
text(0, 0, "default", cex = 1)
circos.initializeWithIdeogram(species = "hg18")
circos.initializeWithIdeogram(species = "mm10")


cytoband.file = system.file(package = "circlize", "extdata", "cytoBand.txt")
circos.initializeWithIdeogram(cytoband.file)

cytoband.df = read.table(cytoband.file, colClasses = c(
  "character", "numeric",
  "numeric", "character", "character"
), sep = "\t")
circos.initializeWithIdeogram(cytoband.df)
circos.clear()


circos.par("gap.degree" = rep(c(2, 4), 12))
circos.initializeWithIdeogram()
circos.clear()



df = data.frame(
  name  = c("TP53", "TP63", "TP73"),
  start = c(7565097, 189349205, 3569084),
  end   = c(7590856, 189615068, 3652765)
)
circos.genomicInitialize(df)

tp_family = readRDS(system.file(package = "circlize", "extdata", "tp_family_df.rds"))
head(tp_family)
circos.genomicInitialize(tp_family)
circos.track(
  ylim = c(0, 1),
  bg.col = c("#FF000040", "#00FF0040", "#0000FF40"),
  bg.border = NA, track.height = 0.05
)
n = max(tapply(tp_family$transcript, tp_family$gene, function(x) length(unique(x))))
circos.genomicTrack(tp_family,
  ylim = c(0.5, n + 0.5),
  panel.fun = function(region, value, ...) {
    all_tx = unique(value$transcript)
    for (i in seq_along(all_tx)) {
      l = value$transcript == all_tx[i]
      # for each transcript
      # logger.info("Transcript: ", all_tx[i])
      current_tx_start = min(region[l, 1])
      current_tx_end = max(region[l, 2])
      circos.lines(c(current_tx_start, current_tx_end),
        c(n - i + 1, n - i + 1),
        col = "#CCCCCC"
      )
      circos.genomicRect(region[l, , drop = FALSE],
        ytop = n - i + 1 + 0.4,
        ybottom = n - i + 1 - 0.4, col = "orange", border = NA
      )
    }
  }, bg.border = NA, track.height = 0.4
)
circos.clear()
