# Worked ggplot2 examples

Complete, self-contained figures. Each block runs on its own with synthetic
data, so it can be pasted into a scratch script and executed before being
adapted. Every example here was executed before being written down.

All of them assume this preamble:

```r
suppressMessages({
  library(jutils)          # ggplot2, patchwork, data.table, glue, saveplot()
  load_pkg(ggrepel)
})

fn_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, color = "grey40", size = 10
      )
    )
}

STATUS_COLORS <- c("kept by QC" = "#B09C85", "removed by QC" = "#E64B35")
GROUP_COLORS <- c(Discovery = "#E64B35", Replication = "#F39B7F",
                  Reference = "#374E55")
```

In real code `fn_theme()` lives in the stage config and every color lives in
the track `color.R`. The hexes are spelled out here only so each block runs
standalone; never copy one into a figure script.

---

## 1. Histogram on a log axis, with zeros floored

The pattern to copy: a log axis silently drops zeros, so floor them at half the
smallest non-zero value, keep them visible as their own stripe, and say so in
both the subtitle and the axis label.

```r
set.seed(1)
d <- data.table::data.table(
  maf = c(rep(0, 40), rlnorm(360, log(0.01), 1)),
  status = sample(names(STATUS_COLORS), 400, TRUE, c(0.7, 0.3))
)
THRESH <- 0.01

nz <- d[maf > 0, min(maf)]
floor_maf <- nz / 2
d[, maf_plot := pmax(maf, floor_maf)]
n_zero <- sum(d$maf == 0)

p <- ggplot(d, aes(maf_plot, fill = status)) +
  geom_histogram(bins = 50, color = "white", linewidth = 0.15) +
  geom_vline(xintercept = THRESH, linetype = "dashed",
             color = "grey30", linewidth = 0.4) +
  scale_x_log10() +
  scale_fill_manual(values = STATUS_COLORS, na.value = "grey80") +
  labs(
    title = "Minor allele frequency at the discovery leads",
    subtitle = glue(
      "{nrow(d)} leads \u00b7 {n_zero} monomorphic, floored for display ",
      "\u00b7 dashed line at MAF = {THRESH}"
    ),
    x = "MAF (log scale)",
    y = "leads",
    fill = NULL
  ) +
  fn_theme()

saveplot("hist_log.pdf", p, width = 8, height = 5)
```

---

## 2. Log-log scatter with an identity line

Reference lines go in before the points. The subtitle explains what a position
relative to the dashed line means, so the reader does not have to infer it.

```r
d <- data.table::data.table(
  x = rlnorm(300, log(0.05), 1),
  y = rlnorm(300, log(0.03), 1.2),
  status = sample(names(STATUS_COLORS), 300, TRUE)
)

p <- ggplot(d, aes(x, y, color = status)) +
  geom_abline(slope = 1, intercept = 0, color = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = 0.01, linetype = "dashed",
             color = "grey30", linewidth = 0.35) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = STATUS_COLORS) +
  labs(
    title = "Allele frequency contrast",
    subtitle = glue(
      "Points below the dashed line are cohort-specific \u00b7 {nrow(d)} leads"
    ),
    x = "Discovery MAF (log scale)",
    y = "Reference MAF (log scale)",
    color = NULL
  ) +
  fn_theme()

saveplot("scatter_loglog.pdf", p, width = 7.5, height = 5.5)
```

---

## 3. Forest, one hit per PDF page

`lapply()` over the hits builds a list; `saveplot()` turns a list into a
multi-page PDF, one page per element. Reversing the factor levels puts the
first group at the top, which is what a reader expects in a forest.

```r
hits <- data.table::data.table(
  id = c("chr1:1000:A:T", "chr7:2500:C:G", "chr12:900:G:A"),
  gene = c("GENE1", "GENE2", "GENE3"),
  p_disc = c(2e-8, 4e-7, 9e-6)
)
data.table::setorder(hits, p_disc)

forests <- lapply(seq_len(nrow(hits)), function(i) {
  r <- hits[i]
  d <- data.table::data.table(
    cohort = factor(names(GROUP_COLORS), levels = rev(names(GROUP_COLORS))),
    beta = rnorm(3, 0.3, 0.15),
    se = runif(3, 0.05, 0.2)
  )
  d[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]

  ggplot(d, aes(beta, cohort, color = cohort)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, na.rm = TRUE) +
    geom_point(size = 2.6, na.rm = TRUE) +
    scale_color_manual(values = GROUP_COLORS, guide = "none") +
    labs(
      title = glue("{r$id}  ({r$gene})"),
      subtitle = glue("discovery P = {signif(r$p_disc, 3)}"),
      x = "beta (trait units)",
      y = NULL
    ) +
    fn_theme()
})

saveplot("forest.pdf", forests, width = 7.5, height = 3.5)
```

`na.rm = TRUE` on the geoms matters: an engine that reports no SE (a
score-test-only model) gives `NA` bounds, and without it every such panel
prints a warning.

---

## 4. Bars with value labels, assembled with patchwork

Value labels need axis headroom or they are clipped. Vertical bars expand the
y axis; horizontal bars need an explicit upper limit.

```r
counts <- data.table::data.table(
  group = factor(c("Female", "Male"), levels = c("Female", "Male")),
  n = c(318L, 380L)
)
counts[, pct := n / sum(n) * 100]

p_bar <- ggplot(counts, aes(group, n, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = glue("{n}\n({round(pct, 1)}%)")),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = c(Female = "#F39B7F", Male = "#8491B4"),
                    na.value = "grey80") +
  labs(title = "Sex distribution", x = NULL, y = "Count") +
  fn_theme() +
  theme(legend.position = "none")

miss <- data.table::data.table(
  variable = c("EF", "GLS", "LAVI", "LVMI"),
  n_avail = c(690L, 512L, 470L, 655L),
  n_total = 698L
)
miss[, pct_avail := n_avail / n_total * 100]
miss[, variable := factor(variable, levels = variable)]

p_miss <- ggplot(miss, aes(variable, n_avail)) +
  geom_col(width = 0.7, fill = "#4DBBD5") +
  geom_text(aes(label = glue("{n_avail}/{n_total} ({round(pct_avail, 1)}%)")),
            hjust = -0.05, size = 3) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 698 * 1.35),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Data availability", x = NULL, y = "N with data") +
  fn_theme()

p <- p_bar + p_miss + patchwork::plot_layout(widths = c(1, 2))
saveplot("panel.pdf", p, width = 12, height = 5)
```

Setting `variable` to a factor with `levels = variable` freezes the row order;
without it `coord_flip()` reorders alphabetically.

---

## 5. Manhattan and QQ with point thinning

A genome-wide scatter of millions of points makes an unusable PDF. Keep every
significant point and a fixed percentage of the rest, then say in the subtitle
that the plot is thinned. Save as PNG, not PDF.

```r
set.seed(2)
n <- 200000
gw <- data.table::data.table(
  CHR = sample(1:22, n, TRUE),
  POS = sample(1e6:2.4e8, n, TRUE),
  P = runif(n)
)
gw[sample(.N, 30), P := runif(30, 1e-12, 1e-7)]
gw[, variant_id := paste0("chr", CHR, ":", POS)]

P_GW <- 5e-8
P_SUGG <- 1e-5
THIN_P <- 0.01
THIN_PCT <- 5

# Thin against the full N so the QQ expected quantiles stay correct.
n_total <- nrow(gw)
keep <- gw$P < THIN_P | runif(n_total) < THIN_PCT / 100
plot_dt <- gw[keep]
data.table::setorder(plot_dt, CHR, POS)

chr_len <- plot_dt[, .(mx = max(POS)), by = CHR][order(CHR)]
chr_len[, chr_offset := cumsum(as.numeric(mx)) - mx]
plot_dt <- merge(plot_dt, chr_len[, .(CHR, chr_offset)], by = "CHR")
plot_dt[, `:=`(
  cum_pos = POS + chr_offset,
  logp = -log10(P),
  chr_band = data.table::fifelse(CHR %% 2 == 0, "even", "odd")
)]
data.table::setorder(plot_dt, P)
plot_dt[, expected := -log10((seq_len(.N) - 0.5) / n_total)]

ticks <- plot_dt[, .(center = (min(cum_pos) + max(cum_pos)) / 2),
                 by = CHR][order(CHR)]
hits <- plot_dt[P < P_GW]

p_man <- ggplot(plot_dt, aes(cum_pos, logp)) +
  geom_point(aes(color = chr_band), size = 0.35, alpha = 0.75) +
  scale_color_manual(values = c(odd = "#80cbc4", even = "#b7a6cc"),
                     guide = "none") +
  geom_hline(yintercept = -log10(P_GW), color = "#ff5a66", linewidth = 0.4) +
  geom_hline(yintercept = -log10(P_SUGG), color = "grey50",
             linetype = "dashed", linewidth = 0.35) +
  scale_x_continuous(breaks = ticks$center, labels = ticks$CHR) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(P))),
    title = "Genome-wide association",
    subtitle = glue("points thinned: all P<{THIN_P} + {THIN_PCT}% of the rest")
  ) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7))

if (nrow(hits) > 0) {
  p_man <- p_man +
    geom_point(data = hits, color = "#ff5a66", size = 1.2) +
    ggrepel::geom_text_repel(data = hits, aes(label = variant_id),
                             size = 2.2, min.segment.length = 0,
                             max.overlaps = 20)
}
saveplot("manhattan.png", p_man, width = 10, height = 4)

lambda_gc <- median(qchisq(1 - plot_dt$P, 1)) / qchisq(0.5, 1)
p_qq <- ggplot(plot_dt, aes(expected, logp)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linewidth = 0.4) +
  geom_point(size = 0.35, alpha = 0.75, color = "#80cbc4") +
  labs(
    x = expression(Expected ~ -log[10](italic(P))),
    y = expression(Observed ~ -log[10](italic(P))),
    title = "QQ",
    subtitle = glue("lambda_GC = {round(lambda_gc, 3)}")
  ) +
  fn_theme()
saveplot("qq.png", p_qq, width = 4.5, height = 4.5)
```

Two things that are easy to get wrong:

- The QQ expected quantiles must be computed against the **unthinned** `n_total`,
  or the diagonal shifts and the plot lies.
- Label only the genome-wide subset. Repelling thousands of labels hangs.
