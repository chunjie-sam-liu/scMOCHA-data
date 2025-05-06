my_data <- mtcars
# Print the first 6 rows
head(my_data, 6)
plot(
  x = my_data$wt, y = my_data$mpg,
  pch = 16, frame = FALSE,
  xlab = "wt", ylab = "mpg", col = "#2E9FDF"
)
plot(mtcars$wt, mtcars$mpg,
  main = "Scatterplot in Base R",
  xlab = "Car Weight", ylab = "MPG",
  pch = 4, col = "blue", lwd = 10, cex = 1
)


library(grid)
text_grob <- textGrob("Hello", x = 0.5, y = 0.5)
grid.draw(text_grob) # 显示它
grid.rect(x = unit(0.5, "npc"), y = unit(1, "lines"), width = unit(3, "cm"))


library(grid)

# 创建 grob
my_grob <- grobTree(
  rectGrob(gp = gpar(fill = "lightblue")),
  textGrob("Hello Grid", x = 0.5, y = 0.5, gp = gpar(fontsize = 20))
)

# 绘图到窗口
grid.newpage()
grid.draw(my_grob)
