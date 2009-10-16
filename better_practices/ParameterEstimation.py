#!/usr/bin/env python
# encoding: utf-8
"""
ParameterEstimation.py
Created by Dave Williams on 2009-09-09.
Copyright (c) 2009 Dave Williams. All rights reserved.
"""


from numpy import pi, degrees, radians, cos, sin, hypot, arctan2


def main():
    """Given some 4sXB parameters, calculate 2sXB parameters and other info"""
    # Diagram of the crossbridge head for reference.
    # 
    # =( )==================   Thick filament       
    #   ^ +++                                       
    #   |    +++  <----------- S2 Region            
    #   |       +++                                 
    #   |          ( ) <------ Angle Phi            
    #   |           \\                              
    #   |            \\  <---- LCD                  
    #   |             \\                            
    #   |             ( )  <-- Angle Omega          
    #   |              ||                           
    #   |              ||  <-- CD                   
    #   |              ||                           
    #   |              88  <-- Tip                  
    #   |                                           
    #   Angle Theta                                 
    # 
    # From Liu et al's 2006 paper (which admittedly is using IFM), the most
    # frequently occurring thick fil to s2 angle range is 51-60 degrees. We 
    # assume that this range is being distorted by the compressive radial 
    # force being generated by the rigor crossbridges in the swollen lattice 
    # spacings that are used in this paper. As such, we choose a rest angle 
    # for the s2 domain at the low end of the still common range of 50deg to /
    # 41deg. We do not change this angle between states one, two and three.
    
    theta = radians(40)
    
    # Also from Liu et al's 2006 paper we get a s2 length of 10.5nm that is 
    # invariant with kinetic state.
    
    s2_len = 10.5
    
    # In Taylor's 1999 Cell paper (clearly explained in Davis and Epstein's  
    # 2009 paper): the angle between what Reedy refers to as the light chain 
    # domain (LCD) and the axial axis goes from 125 degrees to 70 degrees. To 
    # determine the rest angle between the s2 region and the LCD we need only 
    # add theta to Taylor's measurements to get 165deg and 110deg as the 
    # smallest of the two angles between S2 and the LCD.
    
    phi = [radians(165), radians(110)]
    
    # From measurements made with MacPyMol on 1dfk from the PDB, the LCD 
    # (about val766-pro835) is about 9.6nm long while the CD (about 
    # r433-val766) is about 5.8nm long.
    
    lcd_len = 9.6
    
    # Now we calculate the rest locations of the s2-LCD junction (hinge) and 
    # the LCD-CD junction (head) before and after the powerstroke. So, head[0] 
    # is pre powerstroke x and y location, while head[1] is the 
    # post-powestroke location. These names a little bit misleading as the tip 
    # of the CD is more rightly the "head" than the LCD-CD junction. To get 
    # the proper angles to use here, there is a little bit of fiddling 
    # required. We are taking the compliment of the angle between the LCD and 
    # the horizontal axis, as opposed to the angle between the LCD and the S2 
    # region.
    
    hinge = [cos(theta)*s2_len, sin(theta)*s2_len] 
    head = [[hinge[0] + cos(pi - (phi[0] - theta)) * lcd_len,
             hinge[1] + sin(pi - (phi[0] - theta)) * lcd_len], 
            [hinge[0] + cos(pi - (phi[1] - theta)) * lcd_len,
             hinge[1] + sin(pi - (phi[1] - theta)) * lcd_len]]
             
    # Display these results
    
    print "Input for the 4sXB:"
    print "==================="
    print "Theta (deg): %5.2f" % degrees(theta)
    print "S2 domain length (nm): %5.2f" % s2_len
    print "Phi (deg), measured from S2 to the LCD: %5.2f pre-powerstroke" \
            % degrees(phi[0])
    print "Phi (deg), measured from S2 to the LCD: %5.2f post-powerstroke" \
            % degrees(phi[1])
    print "LCD length (nm): %5.2f" % lcd_len
    print " "
    
    print "Results for the 4sXB:"
    print "====================="
    print "S2-LCD joined at : [%5.2f, %5.2f]" \
            % (hinge[0], hinge[1])
    print "LCD-CD joined at : [%5.2f, %5.2f] pre-powerstroke" \
            % (head[0][0], head[0][1])
    print "LCD-CD joined at : [%5.2f, %5.2f] post-powerstroke" \
            % (head[1][0], head[1][1])
    print " "
    
    # Now to calculate the parameters for the 2 spring crossbridge we first 
    # find the lengths of the single linear spring finding the distance of the 
    # head from the base of the S2 region before and after the powerstroke.
    
    lin_len = [hypot(head[0][0], head[0][1]), hypot(head[1][0], head[1][1])]
    
    # And find the rest angles of the torsional spring at the base of the 2sXB 
    # before and after the powerstroke. This is the angle between the 2sXB's 
    # linear spring and the horizontal axis.
    
    tor_ang = [arctan2(head[0][1], head[0][0]), \
               arctan2(head[1][1], head[1][0])]
    
    # Check that the 2sXB's pre- and post-powerstroke head locations are close 
    # to those of the 4sXB.
    
    all_match = all([
        round(cos(tor_ang[0]) * lin_len[0], 3) == round(head[0][0], 3), 
        round(sin(tor_ang[0]) * lin_len[0], 3) == round(head[0][1], 3), 
        round(cos(tor_ang[1]) * lin_len[1], 3) == round(head[1][0], 3), 
        round(sin(tor_ang[1]) * lin_len[1], 3) == round(head[1][1], 3)
        ])
    
    # Display the results for the 2sXB
    
    print "Results for the 2sXB:"
    print "====================="
    print "Head locations matches 4sXB: " + str(all_match)
    print "Lin spring lengs (nm): %5.2f pre and %5.2f post" \
            % (lin_len[0], lin_len[1])
    print "Tor spring angs (deg): %5.2f pre and %5.2f post" \
            % (degrees(tor_ang[0]), degrees(tor_ang[1]))
    print " "
    
    # Find the stength of the springs from the known liberated energy of ATP 
    # hydrolysis.
    
    R = 8.31447 # Joule/(Kelvin Mole)
    T = 288 # Kelvin
    Total_Energy = 24 # RT energy liberated over an actomyosin cycle
    eta = 0.68 # fractional free energy drop from states one to three
    # Typically, for a torsional angle U = .5 * k (angle-rest_angle)^2
    # Rearrange to get  k = 2 U/(angle-angle_rest)^2, use this to estimate the 
    # strength of the angular springs
    k4sXB = (2*eta*Total_Energy)/(phi[0]-phi[1])**2
    k2sXB = (2*eta*Total_Energy)/(tor_ang[0]-tor_ang[1])**2
    
    print "Results for the Torsional spring values:"
    print "========================================"
    print "Torsional strength of 4sXB (in RT): " + str(k4sXB)
    print "Torsional strength of 2sXB (in RT): " + str(k2sXB)
    print " "
    
    
    
    # Calculate the lattice spacing conversion formulas
    
    # From Morel's 1997 work: radial forces shift from compressive to 
    # expansive when lattice spacing drops below 34nm d10, making this a 
    # reasonable value for the resting radial offset of our power-generating 
    # crossbridge heads.
    #
    # From Millman's 1998 paper, we derive that the 
    # myosin-center-to-actin-center distance, times 1.5, yields the d10 
    # lattice spacing that is typically used as a measure of lattice spacing. 
    # So if we want to convert our radial coordinates to the d10 lattice 
    # spacing, we need to find the offset that accounts for thickness of the 
    # filaments and the CD. 
    #
    # Note: This is somewhat less than you would expect if you know the thick 
    # filament diameter is 26 to 31nm (from Millman, 1998 on pg378 and Kensler 
    # & Harris, 2008, PMCID: PMC2242758 respectively) and that the thin 
    # filament diameter is 9.5nm (from Millman, 1998 on pg378). We attribute 
    # this to some of the springs we are describing lying within what would 
    # be considered the thick filament. Also, the crossbridge head may 
    # actually be binding to the side of the thin filament, not to the face 
    # oriented to the thick filament, so that the XB is actually reaching past 
    # the surface of the thick filament when viewed in the radial plane.
    
    # d10 = 1.5 * (our_ls + offset)
    # d10 - 1.5 * our_ls = 1.5*offset
    # d10 / 1.5 - our_ls = offset
    ls_offset = 34/1.5 - head[1][1]
    
    # Display the lattice spacing results
    
    print "Results for lattice spacing:"
    print "============================"
    print "Lattice spacing formula: d10 = 1.5 * (ls + offset)"
    print "Offset: %5.2f" % ls_offset
    
    # From Millman's 1998 review: The rest lattice spacing of muscle (type?)
    # is 37nm d10. This is also a reasonable value for the rest d10 spacing 
    # used above and may be substituted if needed.
    

if __name__ == '__main__':
    main()

