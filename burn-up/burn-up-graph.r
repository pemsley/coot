
ylim=123
xlim=160

prediction_text_x_placement = 42
prediction_text_y_placement = 8

legend_x = 40
legend_y =  7.5


source('arrow.r')

predict = function(x_pos, y_pos) {

   scope_start = a$V3[1]
   scope_end   = a$V3[length(a$V3)]

   done_start = a$V2[1]
   done_end   = a$V2[length(a$V2)]

   s = paste(scope_start, scope_end)
   # print(s)
   s = paste(done_start, done_end)
   # print(s)

   now_day = a$V1[length(a$V1)]

   # print(now_day)

   c1 = done_start
   c2 = scope_start
   m1 = (done_end - done_start)/now_day
   m2 = (scope_end - scope_start)/now_day

   m_diff = m1 - m2
   if (m_diff > 0) {
      X_pred = (c2-c1)/(m1-m2)
      Y_pred_1 = m1 * X_pred + c1
      Y_pred_2 = m2 * X_pred + c2
      r = paste(X_pred+1, Y_pred_1)
      r = paste(r, now_day)
      # print(r)
      days_delta = X_pred - now_day
      today_t = Sys.Date()
      predict_t = today_t + days_delta
      date_s = format(predict_t, format="%d %B %Y")
      t = paste('Projected Release Day:\n', date_s)
      text(x_pos, y_pos, t, pos=3, cex=0.8)
      s = 3 # should depend on xlim, 3 is good when xlim is 200
      # the ratio between s_x and s_y is the ratio of
      # xlim and y_lim
      s_x = 0.008*xlim
      s_y = 0.02 *ylim
      rect(X_pred-s_x, 0.5*(Y_pred_1-s_y), X_pred+s_x, 0.5*(Y_pred_1+s_y), col = 'darkgreen')
      # need list of X values, list of Y values
      # not x,y pairs
      # lines(c(X_pred,X_pred), c(0, Y_pred_1), col='darkgreen', lty=2)
      lines(c(0,X_pred), c(0,      Y_pred_1*0.5), col='darkgreen', lty=2)
      lines(c(0,X_pred), c(c2*0.5, Y_pred_1*0.5), col='darkgreen', lty=2)

      return(r)
   } else {

      nct = paste('Non-converging:', round(m_diff, digits=3))
      print(nct)
      legend(34, 2, nct, cex=0.8, box.lwd=0)

      # (x-values), (y-values)
      lines(c(0, now_day), c(0, done_end*0.5),  col='#202020', lty=2)
      lines(c(0, now_day), c(scope_start*0.5, scope_end*.5), col='#202020', lty=2)

   }
}


a = read.table('burn-up.tab')


# change the plot resolution
# png('burn-up.png', res=92, width=640, height=640)
png('burn-up.png')

do_plot = function() {
   plot(ylim=c(0,ylim), xlim=c(0,xlim), NULL, NULL, t='n',
               main="Coot-3-3-3 Development Progress",
               xlab="Real Days (since development start)",
               ylab="'Dev' Days")

   grid()
   points(a$V1, a$V3*0.5, t='l', lwd=3, lty=1, col='brown')
   points(a$V1, a$V2*0.5, t='l', lwd=3, lty=1, col='black')

   leg.txt <- c("Scope", "Done")
   legend(legend_x, legend_y, legend=leg.txt, col=c("brown", "black"), lty=1:1, lwd=3, cex=0.7)
}

do_plot()

predict(prediction_text_x_placement, prediction_text_y_placement)

# betterArrow(2, 14.5, 2, 11.0, col='grey', code=2)
# text(3, 15.5, labels=" Keele (EM)", col='grey', cex=1.0)

dev.off()

