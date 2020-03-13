#!/usr/bin/env python
"""Simple statistics functions"""

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

def median(values):
    """Calculate the median
    
    ARGUMENTS:
        values - a list of values (does not have to be sorted)
    RETURNS:
        the median
    """
    
    values = sorted(values)
    
    length=len(values)
    
    if length % 2:
        return values[(length-1)//2]
    else:
        lower = values[length//2-1]
        upper = values[length//2]
        return float(lower+upper)/2


def lowerQuartile(values):
    """Calculate the lower quartile 
    
    ARGUMENTS:
        values - a list of values (does not have to be sorted)
    RETURNS:
        the lower quartile
    """
    
    values = sorted(values)
    
    length=len(values)
    
    if length % 2:
        return median(values[0:(length-1)//2])
    else:
        return median(values[0:(length)//2])
