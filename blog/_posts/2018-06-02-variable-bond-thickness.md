---
layout: post
title:  "Variable Bond Thickness"
date: Sat 2 Jun 2018 14:00:26 BST
---

Now we can render our model molecules with "variable" bond thickness, so that ojects
further away have thinner bond widths (and as we zoom in, they get thicker).

![CootVR Demo]({{"../../../images/static-variable-bond-thickness.png"}})

I am not sure if this should be the default, so I didnt make it so.
You can turn it on with:

```lisp
(set-use-variable-bond-thickness 1)
```

or in python

```python
set_use_variable_bond_thickness(1)
```

Request by Erec Stebbins.
