library(circlize)
set.seed(999)
bed = generateRandomBed()
head(bed)


circos.initializeWithIdeogram()
circos.info()
circos.clear()

{
  set.seed(123)

  circos.par(
    "start.degree" = 90
  )
  circos.initializeWithIdeogram(plotType = NULL)
  circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim
      message("xlim: ", paste(xlim, collapse = ", "))
      circos.rect(
        xlim[1],
        0,
        xlim[2],
        1,
        col = rand_color(1)
      )
      circos.text(
        mean(xlim),
        mean(ylim),
        chr,
        cex = 0.7,
        col = "white",
        facing = "inside",
        niceFacing = TRUE,
      )
    },
    track.height = 0.15,
    bg.border = NA
  )
  circos.clear()
}

{
  df = data.frame(
    name = c("TP53", "TP63", "TP73"),
    start = c(7565097, 189349205, 3569084),
    end = c(7590856, 189615068, 3652765)
  )
  circos.genomicInitialize(df)
  circos.clear()
}


{
  tp_family = readRDS(system.file(
    package = "circlize",
    "extdata",
    "tp_family_df.rds"
  ))
  head(tp_family)
  circos.genomicInitialize(tp_family)
  circos.track(
    ylim = c(0, 1),
    bg.col = c("#FF000040", "#00FF0040", "#0000FF40"),
    bg.border = NA,
    track.height = 0.1
  )

  n = max(tapply(tp_family$transcript, tp_family$gene, function(x) {
    length(unique(x))
  }))
  circos.genomicTrack(
    tp_family,
    ylim = c(0.5, n + 0.5),
    panel.fun = function(region, value, ...) {
      # message("region: ", paste(value, collapse = ", "))
      all_tx = unique(value$transcript)
      for (i in seq_along(all_tx)) {
        l = value$transcript == all_tx[i]
        # for each transcript
        current_tx_start = min(region[l, 1])
        current_tx_end = max(region[l, 2])
        circos.lines(
          c(current_tx_start, current_tx_end),
          c(n - i + 1, n - i + 1),
          col = "#CCCCCC"
        )
        circos.genomicRect(
          region[l, , drop = FALSE],
          ytop = n - i + 1 + 0.4,
          ybottom = n - i + 1 - 0.4,
          col = "orange",
          border = NA
        )
      }
    },
    bg.border = NA,
    track.height = 0.4
  )
  circos.clear()
}
{
  circos.par(
    "track.height" = 0.08,
    start.degree = 90,
    canvas.xlim = c(0, 1),
    canvas.ylim = c(0, 1),
    gap.degree = 270,
    cell.padding = c(0, 0, 0, 0)
  )
  circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
  bed = generateRandomBed(nr = 500)
  circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value)
  })
  circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, area = TRUE)
  })
  circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type = "h")
  })
  bed1 = generateRandomBed(nr = 500)
  bed2 = generateRandomBed(nr = 500)
  bed_list = list(bed1, bed2)
  circos.genomicTrack(
    bed_list,
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicLines(region, value, col = i, ...)
    }
  )
  circos.genomicTrack(
    bed_list,
    stack = TRUE,
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicLines(region, value, col = i, ...)
    }
  )
  bed = generateRandomBed(nr = 500, nc = 4)
  circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = 1:4, ...)
  })
  bed = generateRandomBed(nr = 500, nc = 4)
  circos.genomicTrack(
    bed,
    stack = TRUE,
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicLines(region, value, col = i, ...)
    }
  )
  bed = generateRandomBed(nr = 200)
  circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicLines(
      region,
      value,
      type = "segment",
      lwd = 2,
      col = rand_color(nrow(region)),
      ...
    )
  })
  circos.clear()
}

{
  circos.par(
    "track.height" = 0.15,
    start.degree = 90,
    canvas.xlim = c(0, 1),
    canvas.ylim = c(0, 1),
    gap.degree = 270
  )
  circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
  col_fun = colorRamp2(
    breaks = c(-1, 0, 1),
    colors = c("green", "black", "red")
  )
  bed = generateRandomBed(nr = 100, nc = 4)
  circos.genomicTrack(
    bed,
    stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(
        region,
        value,
        col = col_fun(value[[1]]),
        border = NA,
        ...
      )
    }
  )
  bed1 = generateRandomBed(nr = 100)
  bed2 = generateRandomBed(nr = 100)
  bed_list = list(bed1, bed2)
  circos.genomicTrack(
    bed_list,
    stack = TRUE,
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicRect(
        region,
        value,
        ytop = i + 0.3,
        ybottom = i - 0.3,
        col = col_fun(value[[1]]),
        ...
      )
    }
  )
  circos.genomicTrack(
    bed_list,
    ylim = c(0.5, 2.5),
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicRect(
        region,
        value,
        ytop = i + 0.3,
        ybottom = i - 0.3,
        col = col_fun(value[[1]]),
        ...
      )
    }
  )
  bed = generateRandomBed(nr = 200)
  circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region,
      value,
      ytop.column = 1,
      ybottom = 0,
      col = ifelse(value[[1]] > 0, "red", "green"),
      ...
    )
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
  })
  circos.clear()
}
