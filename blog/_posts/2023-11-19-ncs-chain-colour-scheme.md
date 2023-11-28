---
layout: post
title:  "NCS Chain Colour Scheme"
date: Sun 19 Nov 14:52:37 GMT 2023
---

I have recently added the function

`coot.set_colour_by_ncs_chain(imol, goodsell_mode)`

- where `imol` is the molecule index
- `goodsell_mode` is either 0 or 1

The function works as you would expect, given its name - _i.e._ chains that
are related by NCS (or molecular symmetry) have the same colour.

Here, for example, is the `5iv5` assembly (from T4 bacteriophage)

![5iv5 assembly top]({{"../../../images/5iv5-top.png"}})

![5iv5 assembly bottom]({{"../../../images/5iv5-bottom.png"}})

 (I also turned on the outline mode):

I think that there's a bug in the assembly that I got from the RCSB (for whatever that's worth) - you might be able to tell from the occassional overlapping chain colours.

