---
layout: post
title:  "The Coot Spin Test"
date: Wed 6 Jun 17:22:24 BST 2018
---

1. Open a new coot, with no state script

2.  Edit -> Map Parameters

 * Set the map radius to 12 
 
 * Set the map sampling rate to 1.8

4. Extensions -> Load Tutorial Model and Data

5. Calculate -> FPS - Yes

6. Zoom in until the density touches the edge of the screen

7. Draw -> Spin View on

8. Watch the status bar...

[Note the FPS can top out at 60 fps - synching with the monitor refresh]

For me I had to turn up the radius to 60.0 A to get a frame rate of less
than 60 FPS (it was 50 FPS).  Perhaps that should be the new standard.

