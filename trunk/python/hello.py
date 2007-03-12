# translated by BL
import time
import os

user = os.getenv('USERNAME')

hour = int(time.strftime("%H", time.localtime()))
if hour < 12: time_str = "Morning"
elif hour < 18: time_str = "Afternoon"
else : time_str = "Evening"
hello_str =  "Good %(a)s %(b)s, Welcome to Coot." % {"a":time_str, "b":user}
# alternative
# hello_str1 =  "%s %s %s%s" % ("Good",time_str, user, ", Welcome to Coot.")

print hello_str
# print hello_str1
set_display_intro_string(hello_str)
