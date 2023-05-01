
# max_reso = 2 (say) and step is 0.05
reciprocal_resolution_scale_plot <- function(a, max_reso, reso_step) {

   rr_max_reso = 1/(max_reso * max_reso)
   rr = seq(0.0, rr_max_reso, reso_step)
   print(rr)
   ln = round(sqrt(1/rr), 2)
   x_label = xlab=paste("Resolution", expression(A))
   plot(a$resolution, a$fsc, type='l', xlim=c(0.01,rr_max_reso), xaxt='n', frame.plot=TRUE, xlab=x_label, ylab='FSC', main='Fourier Shell Correlation')
   axis(1, at=rr, labels=ln)
   # resolution_lines(seq(3, 1.5, -0.1))
   abline(h=0, col='lightgrey')
   abline(h=1, col='grey720')
   # points(a$r, a$fsc, cex=0.3, col='black', pch=19, type='l', lwd=2)
}

a = read.table('a.tab', header=TRUE)
a = read.table('0037.tab', header=TRUE)

max_reso = 1.45;

reciprocal_resolution_scale_plot(a, max_reso, 0.05);

