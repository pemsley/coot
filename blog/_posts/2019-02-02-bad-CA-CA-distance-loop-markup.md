---
layout: post
title:  "'Too Long' CA-CA distance loop markup"
date: Sat 2 Feb 2019 00:02:24 GMT
---

At the request of John Berrisford and Oliver Clarke, I have added changed the
representation of loops where the inter-residue number difference is too short
for the actual distance between the residue - to say it another way, there are
not enough residues to span the missing gap.

Such loops are now represented as a segmented straight line with organge warning
dumbells.

![Bad CA-CA loop distance screenshot]({{"../../../images/screnshot-bad-dist-CA-CA.png"}})

