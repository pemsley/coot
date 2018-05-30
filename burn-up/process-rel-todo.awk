BEGIN { pre_sum = 0; post_sum = 0; add_to_pre_sum = 1; } 

/^------/ { add_to_pre_sum = 0; } 

{ for (i=1; i<NF; i++) {
	if ($i == ":::") {
	    if (add_to_pre_sum) {
                # print "adding", $(i+1)+0, "to pre_sum"
                if ($(i+2) == "h")
		   pre_sum += ($(i+1)+0) * 0.25
                else
		   pre_sum += $(i+1)+0
	    } else { 
                # print "adding", $(i+1)+0, "to post_sum"
                if ($(i+2) == "h")
		   post_sum += ($(i+1)+0) * 0.25
                else
		   post_sum += $(i+1)+0
	    }
	}
    }
}

END {
    x  = system("date +%s > proc-rel.tmp")
    x2 = system("date > proc-rel-2.tmp")
    getline < "proc-rel.tmp";
    # n_s = $1-1229420879; for 0.6
    # n_s = $1-1260164717; for 0.6.1
    # n_s = $1-1264692452; 0.6.2
    # n_s = $1-1312898037; 0.7
    # n_s = $1 - 1502356115
    # n_s = $1 - 1498000000
    # n_s = $1 - 1519647801 # end of Feb
    # yesterday in guile:  (- (current-time) (* 24 60 60))
    n_s = $1 - 1524348970
    n_days = n_s/(60*60*24)
    getline < "proc-rel-2.tmp";
    date_line = $0
    total = pre_sum + post_sum
    printf("%7.3f %6.2f %6.2f   # %s\n", n_days, post_sum, total, date_line)
}

