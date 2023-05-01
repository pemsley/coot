r = test_internal()
print "test_interanl returned status", r
if r:
  coot.coot_real_exit(0)
else:
  coot.coot_real_exit(1)
