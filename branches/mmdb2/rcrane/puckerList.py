#!/usr/bin/env python
"""Lists of sugar puckers for each conformer

puckerList          - a dictionary of rotamerName => [starting pucker as int, ending pucker as int]
puckerListFullName  - a dictionary of rotamerName => [starting pucker as string (i.e. "C3'"), ending pucker as string]
rotList             - a list of all rotamers (without any pucker information)
rotsByPucker        - a list of all rotamers organized by [startingPucker][endingPucker]
                      so rotsByPucker[3][3] is a list of all C3'-C3' rotamers
"""

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

puckerList = { '1a': [3,3],
               '1m': [3,3],
               '1L': [3,3],
               '&a': [3,3],
               '7a': [3,3],
               '3a': [3,3],
               '9a': [3,3],
               '1g': [3,3],
               '7d': [3,3],
               '3d': [3,3],
               '5d': [3,3],
               '1e': [3,3],
               '1c': [3,3],
               '1f': [3,3],
               '5j': [3,3],
               
               '1b': [3,2],
               '1[': [3,2],
               '3b': [3,2],
               '1z': [3,2],
               '5z': [3,2],
               '7p': [3,2],
               '1t': [3,2],
               '5q': [3,2],
               '1o': [3,2],
               '7r': [3,2],
               
               '2a': [2,3],
               '4a': [2,3],
               '0a': [2,3],
               '#a': [2,3],
               '4g': [2,3],
               '6g': [2,3],
               '8d': [2,3],
               '4d': [2,3],
               '6d': [2,3],
               '2h': [2,3],
               '4n': [2,3],
               '0i': [2,3],
               '6n': [2,3],
               '6j': [2,3],
               
               '2[': [2,2],
               '4b': [2,2],
               '0b': [2,2],
               '4p': [2,2],
               '6p': [2,2],
               '4s': [2,2],
               '2o': [2,2],
               
               '5n': [3,3],   #These conformers are wannabes
               '3g': [3,3],
               #'5p': [3,2],
               '5r': [3,2],
               '2g': [2,3],
               #'0k': [2,3],
               '2z': [2,2],
               '2u': [2,2]
             }

puckerListFullName = dict([(rot, map(lambda x: "C%i'" % x, puckers)) for (rot, puckers) in puckerList.iteritems()])

rotList = puckerList.keys()

#create a list of rotamers searchable by leading and ending puckers
rotsByPucker = {2: {}, 3: {}}
for (leadingPucker, endingPucker) in ((2,2), (2,3), (3,2), (3,3)):
    rotsByPucker[leadingPucker][endingPucker] = [rot for rot in puckerList.keys() if puckerList[rot][0] == leadingPucker and puckerList[rot][1] == endingPucker]