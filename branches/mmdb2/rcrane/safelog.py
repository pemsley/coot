#!/usr/bin/env python
"""A log function that returns negative infinity for log(0) or log(x<0)"""

# Copyright 2010 Kevin Keating
# 
# Licensed under the Educational Community License, Version 2.0 (the
# "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
# 
# http://www.osedu.org/licenses/ECL-2.0
# 
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS"
# BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing
# permissions and limitations under the License.

from math import log

negInf = float("-infinity")   #a stand-in for negative infinity
#negInf = -999999999          #use this if float("-infinity") doesn't work (i.e. on Windows with Pythons earlier than 2.5)


def ln(a):
    """Calculate the log without raising an error
    
    ARGUMENTS:
        a - the value to calculate the log of
    RETURNS:
        log(a) if it is defined
        negative infinity otherwise (when a <= 0)
    """
    
    ln = None
    
    try:
        ln = log(a)
    except (ValueError, OverflowError):
        ln = negInf
    
    return ln
