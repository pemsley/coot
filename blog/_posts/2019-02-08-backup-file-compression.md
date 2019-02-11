---
layout: post
title:  "Compression of Backup Files"
date: Fri 8 Feb 2019 03:43:08 PST
---

Sometime Coot reports that it fails to compress the backup files in the coot-backup directory.
This might be, perhaps, because the full path of the directory contains a space. In such a case
you can turn off the backup compression:

{% highlight python %}
set_backup_compress_files(0)
{% endhighlight %}

