---
layout: post
title:  "User-defined Unimodal Torsion Restraints"
date:   2018-05-28 04:23:00 +0100
# categories: update
---

``commit 47d61237b3e645087898a0dd9ec8c790e44dcd6f``


Now we can specify unimodal torsion restraints for our own rings for 
our own specified residues, instead of using the types that _Coot_ has built-in:

{% highlight python %}
use_unimodal_ring_torsion_restraints_for_residue("AMV",
                                                 [
                                                     [ " C1 ", " C2 ", " C3 ", " C4 ", -48.2 ],
                                                     [ " C2 ", " C3 ", " C4 ", " C5 ",  49.66],
                                                     [ " C3 ", " C4 ", " C5 ", " O5 ", -56.15],
                                                     [ " C4 ", " C5 ", " O5 ", " C1 ",  63.76],
                                                     [ " C5 ", " O5 ", " C1 ", " C2 ", -62.4],
                                                     [ " O5 ", " C1 ", " C2 ", " C3 ",  53.34]])
{% endhighlight %}
This makes new ring torsion restraints with an esd of 4.0°.


The idea being that you can control the pyranose ring pucking if you wish.

