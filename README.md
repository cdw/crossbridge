crossbridge
===========

This provides both the mathematical and data generation underpinnings for a single myosin head/crossbridge of the spatially explicit variety. What this translates to is a bunch of relatively simple math and khemical kinetics surrounding a couple of springs (of linear and torsional flavors). 

Types of crossbridges
---------------------

Historically, all models of muscle have treated crossbridges as linear springs aligned with the axis of contraction. This is a big simplification. Look over [here](http://imgur.com/ZBCNxyw) for a glimpse at how myosin actually displaces itself with a rotation of a lever arm over the binding head (image stolen from [here](http://emboj.embopress.org/content/21/11/2517)).

So instead of this paradigm, we created a model of the crossbridge that treats it as two linear and two torsional (angular) springs. This works really well to recapture more of the molecular dynamics we expect from myosin but is too slow to use in larger simulations as it requires optimization to solve for its lowest energy state when it is coupled to a binding site that moves. So we also create a model which treats the crossbridge as a single linear and single torsional spring. This runs a lot faster and is useful in larger scale models of muscle. 

Crossbridge organization
------------------------

The crossbridges used are broken down according to the following schematic: 
```
           H   - Head of the myosin
           |  
           G   - Globular domain
           |  
           C   - Converter region
          /   
         N     - Neck region
        /     
    ===T====== - Thick filament
```
