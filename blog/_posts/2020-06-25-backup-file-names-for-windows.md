---
layout: post
title:  "Backup file-names for Windows"
date: Thu 25 Jun 10:51:24 BST 2020
---

I have head reports from Windows users that their backup file names fail to be written out.
I think that this is because there are colons in the file name. To prevent this, use
the following setting:

{% highlight python %}
set_decoloned_backup_file_names(1)
{% endhighlight %}

It might work.

(Available in 0.9.1)

