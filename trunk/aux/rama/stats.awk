BEGIN { sum = 0; max = 0; min = 999; sumsq = 0; n=0;}

NF == 3 { v = $3+0;
  sum += v;
  sumsq += v*v;
  n++;
  if (v < min)
    min = v;
  if (v > max)
    max = v;
}

END { 
  av = sum/n;
  var=sumsq/n - av*av;
  std=sqrt(var);
  print "mean:", av, "stddev:",std,"max:",max,"min:",min;
}

  