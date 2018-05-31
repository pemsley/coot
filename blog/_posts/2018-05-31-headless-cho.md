---
layout: post
title:  "Headless N-linked CHO"
date: Thu 31 May 20:43:55 BST 2018
---

Now you can build N-linked carbohydrate without a gui:

``$ coot --no-graphics --script headless-cho.scm``

where _``headless-cho.scm``_ is:

```lisp
(let ((imol (read-pdb "pdb2qc1-no-cho.ent")
   (make-and-draw-map "2qc1_map.mtz" "FWT" "PHWT" "" 0 0)
   (add-linked-residue-tree 0 (list "B" 141 "") oligomannose-tree) ;; add to this ASN
   (write-pdb imol "with-cho-added.pdb"))
```
