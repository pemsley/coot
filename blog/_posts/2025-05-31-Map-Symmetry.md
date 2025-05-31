---
layout: post
title:  "x"
date: Sat 31 May 18:37:05 BST 2025
---

As of version 1.1.16, Coot and Moorhen use the following "algorithm" to determine if the map is an EM map:

  - If the unit cell angles are 90 degrees
  - If the Space Group is 1 (P1)
  - If the labels do not contain a label with the text "APPLY-SYMMETRY::"
  - If the labels do not contain a label with the text "PANDDA::"
  - Then the map is an "EM map"

EM maps do not display density outside of the "Unit Cell."

So, to make Coot and Moorhen display symmetry-related density add a label "APPLY-SYMMETRY::".


