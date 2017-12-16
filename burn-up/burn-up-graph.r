
source('arrow.r')

predict = function() {

   scope_start = a$V3[1]
   scope_end   = a$V3[length(a$V3)]

   done_start = a$V2[1]
   done_end   = a$V2[length(a$V2)]

   s = paste(scope_start, scope_end)
   print(s)
   s = paste(done_start, done_end)
   print(s)

   now_day = a$V1[length(a$V1)]

   print(now_day)

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
      print(r)
      days_delta = X_pred - now_day
      print(days_delta)
      today_t = Sys.Date()
      predict_t = today_t + days_delta
      date_s = format(predict_t, format="%d %B %Y")
      t = paste('Projected Release Day:\n', date_s)
      text(50, 190, t, pos=3, cex=0.8)
      s = 3
      rect(X_pred-s, Y_pred_1-s, X_pred+s, Y_pred_1+s, col = 'darkgreen')
      print(paste('X_pred', X_pred))
      print(paste('Y_pred_1', Y_pred_1))
      # need list of X values, list of Y values
      # not x,y pairs
      # lines(c(X_pred,X_pred), c(0, Y_pred_1), col='darkgreen', lty=2)
      lines(c(0,X_pred), c(0, Y_pred_1), col='darkgreen', lty=2)
      lines(c(0,X_pred), c(c2, Y_pred_1), col='darkgreen', lty=2)

      return(r)
   } else {
      print('bad m_diff')
      print(m_diff)
   }
}


a = read.table('burn-up.tab')

png('burn-up.png')

plot(ylim=c(0,220), xlim=c(0,280), a$V1, a$V2, t='l', lwd=2,
            main="Coot-0.8.9 Development Progress",
            xlab="Days (since development start)",
            ylab="Dev Points (aka 'Half-Days')")
points(a$V1, a$V3, t='l', lwd=2, lty=2)

leg.txt <- c("Done", "Scope")
legend(180, 30, legend=leg.txt, lty=1:2, lwd=2, cex=0.7)

text(175, 175, labels="CSHL Purge", col='grey', cex=0.7)
# arrows(160, 180, 118, 200, code=2, cex=0.5)
betterArrow(160, 180, 118, 200, col='grey', code=2)

predict()

dev.off()

