
# 15.3c Arrows =========================================================== 15.3c 
 
# A solid arrowhead function
# 
# This function draws an arrow from (x0,y0) to (x1,y1) with a solid
# arrowhead. 
#
# The default arrowhead length is .25 inches. 
# The default angle of the arrowhead edges is 30 degrees.
# The other parameters are standard for line and polygon formatting.
#
# The basic logic here is to use matrix rotation to translate the arrowhead
# coordinates.
#
# Note that R's trigonometric functions work with radians rather than degrees
#
betterArrow = function(x0, y0, x1, y1,     # Set up arrow function ------------+
  L = .25,                             # Default arrowhead length              |
  angle = 15,                          # Default angle of arrowhead            |
  code = 2,                            # Default arrowhead at x1,y1            |
  col = par("fg"),                     # Default color                         |
  ljoin = par("ljoin"),                # Default line joint style (0)          |
  lty = par("lty"),                    # Default line type                     |
  lwd = par("lwd"),                    # Default line width                    |
  xpd = FALSE){                        # Default stay within plot area         |
                                       # Start function code                   |
    if(code == 1){                     # Reverse arrow direction               |
      tmp = x1; x1 = x0; x0 = tmp      # Switch x values                       |
    }                                  #                                       |
    if(code == 1){                     #                                       |
      tmp = y1; y1 = y0; y0 = tmp      # Switch y values                       |
    }                                  #                                       |
#                                                                              |
# We need to control for the aspect ratio or for different x,y scales          |
#   in setting up the arrow heads. We'll do that by translating the            |
#   usr parameter setting from the original x,y units to units based           |
#   on the plot dimensions measured in inches.  This will allow us to          |
#   adjust the angles to account for different x and y scales. Note,           |
#   however, that rescaling the plot after it is drawn will distort            |
#   the arrowheads.                                                            |
#                                                                              |
    X0 = (x0 - par()$usr[1])/(par()$usr[2] - par()$usr[1]) * par()$fin[1]  #   |
    Y0 = (y0 - par()$usr[3])/(par()$usr[4] - par()$usr[3]) * par()$fin[2]  #   |
    X1 = (x1 - par()$usr[1])/(par()$usr[2] - par()$usr[1]) * par()$fin[1]  #   |
    Y1 = (y1 - par()$usr[3])/(par()$usr[4] - par()$usr[3]) * par()$fin[2]  #   |
#                                                                              |    
    oldusr = par("usr")                # Save original usr settings            |
    par(usr = c(0, par("fin")[1],      # Set up new usr settings based         |
              0, par("fin")[2]))       #  on plot dimensions in inches         |
#                                                                              |
    t = angle * pi/180                 # Convert angle degrees to radians      |
    slope = (Y1 - Y0)/(X1 - X0)        # Calculate slope of line               |
    S = atan(slope)                    # Slope angle in radians                |
                                       #                                       |
    M = ifelse(X1 < X0, -1, 1)         # Set a marker for X1 < X0              |
                                       # Set length of vector XA               |
    XA = sqrt((X1 - X0)^2 + (Y1 - Y0)^2)                                   #   |
#                                                                              |
# Get arrowhead vertices from rotated vectors                                  |
    XC = X0 + M * ((XA - L) * cos(S) + L * tan(t) * sin(S))                #   |
    YC = Y0 + M * ((XA - L) * sin(S) - L * tan(t) * cos(S))                #   |
    XB = X0 + M * ((XA - L) * cos(S) - L * tan(t) * sin(S))                #   |
    YB = Y0 + M * ((XA - L) * sin(S) + L * tan(t) * cos(S))                #   |
#                     |                                                        |
# Draw arrow line stopping at beginning of the arrowhead                       |
    lines(x = c(X0, X1 - M * L * cos(S)),                                  #   |
      y = c(Y0, Y1 - M * L * sin(S)),                                      #   |
      lty = lty, lwd = lwd,            # Apply line format options             |
      col = col, xpd = xpd,            #                                       |
      ljoin = ljoin)                   #                                       |
    polygon(x = c(X1, XB, XC),         # Draw arrow head                       |
      y = c(Y1, YB, YC),               #  at vertices                          |
      col = col,                       # Apply format options                  |
      border = col,                    #                                       |
      xpd = xpd)                       #                                       |
    par(usr = oldusr)                  # Reset to original usr settings        |
}                                      # End of myArrow function --------------+
 
