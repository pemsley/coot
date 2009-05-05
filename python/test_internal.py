r = test_internal()
print "test_interanl returned status", r
if r:
  coot_real_exit(0)
else:
  coot_real_exit(1)
