
# coot.log must have standard error included and coot process not terminated by ^C else the
# ending log is not writen.

BEGIN{ print_this = 0}

/^PASS: |^UNTESTED: |^FAIL:|^UNRESOLVED:|^UPASS: / 

/^Entered testcase/ {print "   ",$0}

/ === greg-tests Summary ===/ {print_this = 1; print ""}

/greg-tests - / {print $0; next }
/Exception: / {print $0; }

print_this == 1 {print $0}

/# of files abandoned/ {print_this = 0 }

