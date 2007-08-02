
# coot.log must have standard error included and coot process not terminated by ^C else the
# ending log is not writen.

BEGIN{ print_this = 0}

/^PASS: |^Entered testcase |^UNTESTED: |^FAIL|^UNRESOLVED: / 

/ === greg-tests Summary ===/ {print_this = 1; print ""}

print_this == 1 {print $0}

/# of files abandoned/ {print_this = 0 }

