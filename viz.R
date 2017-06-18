 ## Draw some nice visualizations from the kalman filter
## output file, and compute the state  covariance matrix. 
library(ggplot2)

outfile <- "build/out.txt"
outfile <- read.table(outfile, sep="\t", header=TRUE)
outfile <- data.frame(outfile)

N <- nrow(outfile)

rmse <- sqrt((1/N) * unlist(Map(function(x) { sum(x^2) },
                                cbind(outfile['px_state']-outfile['px_ground_truth'],
                                      outfile['py_state']-outfile['py_ground_truth'],
                                      outfile['v_state']*cos(outfile['yaw_angle_state'])-outfile['vx_ground_truth'],
                                      outfile['v_state']*sin(outfile['yaw_angle_state'])-outfile['vy_ground_truth']))))

pl.data <- data.frame(x = c(outfile[,'px_state'],
                            outfile[,'px_measured'],
                            outfile[,'px_ground_truth']),
                      y = c(outfile[,'py_state'],
                            outfile[,'py_measured'],
                            outfile[,'py_ground_truth']),
                      l = c(rep('e', N),
                            rep('m', N),
                            rep('t', N)))
                            

pl.shapes <- c("e"=16,      "m"=17,     "t"=15)
pl.colors <- c("e"="black", "m"="blue", "t"="firebrick")
pl.alpha  <- c("e"=1,       "m"=0.5,    "t"=0.5)

head(pl.data)

pl <- ggplot(data=pl.data) +
  geom_point(aes(x, y, color=l, shape=l,  alpha=l), size=5) +  
  scale_x_continuous() +
  scale_y_continuous() +
  scale_color_manual(name="uniqueL", values=pl.colors, labels=c("Truth", "Measurement", "Prediction")) +
  scale_shape_manual(name="uniqueL", values=pl.shapes, labels=c("Truth", "Measurement", "Prediction")) +
  scale_alpha_manual(name="uniqueL", values=pl.alpha, labels=c("Truth", "Measurement", "Prediction")) +
  theme(axis.line = element_blank(),
        panel.grid.major = element_line(colour="lightgray"),
        panel.grid.minor = element_line(colour="lightgray"),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        text = element_text(size=30),
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        legend.background = element_blank(),
        legend.key = element_blank()) + 
  ggtitle("Unscented Kalman Filter (position)") +
  annotate("text", label=paste("RMSE: [", paste(round(rmse, 3), collapse=", "),  "]", sep=" "),
           x=10, y=-4, size=8) + 
  labs(x="x position", y="y position")


png("sample-data-output.png", width=1024, height=1024)
pl
dev.off()

## Plot the NIS.
pl.data <- data.frame("x" = 0:(N-1),
                      "y" = outfile[,'NIS'])
pl <- ggplot(data = pl.data) +
  geom_point(aes(x, y), col = "darkgray") +
  theme(axis.line = element_blank(),
        panel.grid.major = element_line(colour="lightgray"),
        panel.grid.minor = element_line(colour="lightgray"),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        text = element_text(size=30),
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        legend.background = element_blank(),
        legend.key = element_blank()) + 
  ggtitle("Unscented Kalman Filter NIS") +
  labs(x = "Timestep (s)", y = "NIS") +
  ylim(c(0,20)) +
  geom_hline(yintercept = 7.815, color = 'red') +
  annotate("segment", x = 75, y = 10, xend = 75, yend = 7.815, col = 'red', arrow = arrow()) + 
  annotate("text", label = "ChiSq NIS 95%", x = 75, y = 10.5 ,size = 5)




png("sample-data-NIS.png", width=1024, height=512)
pl
dev.off()

                      
