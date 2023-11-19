---
layout: post
title:  "Setting the Shadow-Box Size"
date: Sun 19 Nov 15:11:22 GMT 2023
---

Most times, the shadow box will just work. However, with large assemblies
one might need to change the shadow box size.

The function to do this is

`coot.set_shadow_box_size(new_size)`

  - where `new_size` is 200 or some such.


Default:

![default shadow box]({{"../../../images/default-shadow_box.png"}})

Resized to 240:

![shadow box at 240]({{"../../../images/shadow-box-240.png"}})
