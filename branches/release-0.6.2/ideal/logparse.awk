BEGIN {

  if (iter_no == "") {
    iter_no = 1; # default
  }

  num_grad_on = 0;
  ana_grad_on = 0;

  ana_No = 0;
  num_No = 0;
}

/numerical_gradients/ {
  num_grad_on = 1;
  ana_grad_on = 0;
  num_No++;
  num_count = 0;
  next;
}


/analytical_gradients/ {
  ana_grad_on = 1;
  num_grad_on = 0;
  ana_No++;
  ana_count = 0;
  next;
}

/debug/ {
  ana_grad_on = 0;
  num_grad_on = 0;
}

num_grad_on {  if (NF == 2) num_list[num_No,num_count++] = $1+0; }

ana_grad_on {  if (NF == 2) ana_list[ana_No,ana_count++] = $1+0; }

END {

  if (0) { 
    for (i=0; i< num_count; i++) {
      print num_list[iter_no,i], ana_list[iter_no,i];
    }
  }

  if (1) { 
    lim = 770;
    for (j=0; j<lim; j++) { 
      for (i=0; i< num_count; i++) {
	print num_list[j,i], ana_list[j,i];
      }
    }
  }
}

