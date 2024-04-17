undirected-link-breed [redlink redlinks]
breed [cells cell]
cells-own
[
  old-voltage  ;; the voltage of the patch the last time thru go
  voltage  ;; the current voltage of the patch
  old-recov
  GNa
  ENa
  GK
  GM
  EK
  GL
  EL
  m
  mold
  n
  nold
  h
  hold
  istim
  mem1
  CM
]
globals
[
  plate-size  ;; the size of the plate on which heat is diffusing
  ;; Used for scaling the color of the patches
  min-voltage  ;; the minimum voltage at setup time
  max-voltage  ;; the maximum voltage at setup time
  num-turty
  num-turtx
  counter
  dt
  dxx
  dxxnr
  conductance
  outvol1
  outvol2
  T1
  TF
  radius
  internode
  myelin-patches

]

to setup
  clear-all

  set internode 4

  let rho .07
  let fiberdiam 0.001; ;; 0.001 = 10 um

  set dt 0.0001
  set dxx 100 * fiberdiam / internode
  set dxxnr 0.0003


  set T1 20  ;; set the temperature
  set TF  ( 3 ^ (( T1 - 6.3 ) / 10 )) ;; adjust the rate constants

  set radius 0.7 * fiberdiam / 2; %radius
  let resistance rho * dxx / (pi * radius ^ 2)
  set conductance  1.0 / (resistance)
  show conductance


  ask patches [set pcolor 89.9]


  create-cells 100 [set shape "square"  set size 4]

  let kkk 0
  while [ kkk < 100  ] [
  let llc-x kkk * 4 + 2 + 2
  let llc-y 2
  set myelin-patches (patch-set patches with [pxcor >= llc-x and pxcor < (llc-x + (internode * 4) - 2 ) and pycor >= llc-y and pycor < (llc-y + 8)])
  ask myelin-patches [
      set pcolor [green] of patch llc-x llc-y
    ]
    set kkk kkk + internode
    show llc-x
  ]

  if demyelinate? [
    let llc-x1 (site1 - 1 ) * internode * 4 + 2 + 2
    let llc-y1 2
    set myelin-patches (patch-set patches with [pxcor >= llc-x1 and pxcor < (llc-x1 + (internode * 4) - 2 ) and pycor >= llc-y1 and pycor < (llc-y1 + 8)])
    ask myelin-patches [
      set pcolor [white] of patch llc-x1 llc-y1]

    let llc-x2 ( site2 - 1) * internode * 4 + 2 + 2
    let llc-y2 2
    set myelin-patches (patch-set patches with [pxcor >= llc-x2 and pxcor < (llc-x2 + (internode * 4) - 2 ) and pycor >= llc-y2 and pycor < (llc-y2 + 8)])
    ask myelin-patches [
    set pcolor [white] of patch llc-x2 llc-y2]

    let llc-x3 ( site3 - 1) * internode * 4 + 2 + 2
    let llc-y3 2
    set myelin-patches (patch-set patches with [pxcor >= llc-x3 and pxcor < (llc-x3 + (internode * 4) - 2 ) and pycor >= llc-y3 and pycor < (llc-y3 + 8)])
    ask myelin-patches [
    set pcolor [white] of patch llc-x3 llc-y3]

    let llc-x4 ( site4 - 1) * internode * 4 + 2 + 2
    let llc-y4 2
    set myelin-patches (patch-set patches with [pxcor >= llc-x4 and pxcor < (llc-x4 + (internode * 4) - 2 ) and pycor >= llc-y4 and pycor < (llc-y4 + 8)])
    ask myelin-patches [
    set pcolor [white] of patch llc-x4 llc-y4
    ]
]

  ask patch (rec1 * 4 + 2) 9 [set plabel "*" set plabel-color black]
  ask patch (rec2 * 4 + 2) 9 [set plabel "*" set plabel-color black]


  let i 0
  let j 0
  while [ i < 100 ] [
  ask cell i [ setxy j + 2 5 ]
  set i i + 1
  set j j + 4
  ]
  let k 0
  while [ k < 99 ] [

  ask cell k [ create-link-with cell (k + 1) ]
    set k k + 1
  ]


  ask cells ;;corresponds to the internode
  [
   set voltage 0 - 60
    set GNa 0
    set GK 0
    set GL 0
    set GM .0000056
    set ENa 55.0
    set EK -72.0
    set EL -0.05 - 60.0
    set CM 0.001
    set mem1 dxx    ;; mem1 in the internode
    set istim 0

    let am TF * (0.1 * (25.0 - (voltage + 60.0) )) / (e ^ (2.5 - 0.1 * (voltage + 60.0 )) - 1)
    let bm TF * (4 * e ^ ((- (voltage + 60.0))/ 18))
    let ah TF * (0.07 * e ^ ( (- (voltage + 60.0) )/ 20))
    let bh TF * (1 / (e ^ (3 - 0.1 * ((voltage + 60.0))) + 1))
    let an TF * ((0.1 - 0.01 * ((voltage + 60))) / (e ^ (1 - 0.1 * ((voltage + 60.0))) - 1))
    let bn TF * (0.125 * e ^ ((- (voltage + 60.0) )/ 80) )

    set m am / (bm + am)
    set n an / (bn + an)
    set h ah / (bh + ah)

    set old-voltage voltage
    set mold m
    set nold n
    set hold h

    set color gray
  ]

    let iii 0
    while [ iii < 100 ] [


    ask cell iii ;;corresponds to the node

    [

    set GNa 1200.0
    set GK 90.0
    set GL 20
    set GM .0000056
    set mem1 dxxnr
    set CM 1.0


    ]
    set iii iii + internode


  ]


   if demyelinate? [
    let iiii 0
    while [ iiii < 3 ] [



    ask cell ((site1 - 1) * internode + 1 + iiii ) ;;corresponds to the node
      [ set color blue set CM Capacitance_Factor * 0.001]
    ask cell ((site2 - 1) * internode + 1 + iiii ) ;;corresponds to the node
      [ set color blue set CM Capacitance_Factor * 0.001]
    ask cell ((site3 - 1) * internode + 1 + iiii ) ;;corresponds to the node
     [ set color blue set CM Capacitance_Factor * 0.001]
    ask cell ((site4 - 1) * internode + 1 + iiii ) ;;corresponds to the node
     [ set color blue set CM Capacitance_Factor * 0.001]
      set iiii iiii + 1
    ]
  ]



  reset-ticks
end




;;;;;;;;;;;;;;;;;;;;;;;;
;;; Setup Procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Runtime Procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Runs the simulation through a loop
to go
  ask cells [
    ;;

    let am TF * (0.1 * (25.0 - (old-voltage + 60.) )) / (e ^ (2.5 - 0.1 * ((old-voltage + 60.0) )) - 1)
    let bm TF * (4 * e ^ ((- (old-voltage + 60.0))/ 18))
    let ah TF * (0.07 * e ^ ( (- (old-voltage + 60.0) )/ 20))
    let bh TF * (1 / (e ^ (3 - 0.1 * ((old-voltage + 60.0))) + 1))
    let an TF * ((0.1 - 0.01 * ((old-voltage + 60.0))) / (e ^ (1 - 0.1 * ((old-voltage + 60.0))) - 1))
    let bn TF * (0.125 * e ^ ((- (old-voltage + 60.0) )/ 80) )


    set m mold + dt * ( am * (1 - mold) - bm * mold)
    set n nold + dt * ( an * (1 - nold) - bn * nold)
    set h hold + dt * ( ah * (1 - hold) - bh * hold)

    let INa GNa * mold ^ 3 * hold * (old-voltage - ENa)
    let IK  GK *  nold ^ 4 * (old-voltage - EK)
    let IL  GL * (old-voltage - EL)
    let IM  Gm * (old-voltage)

    set voltage old-voltage + (dt / CM) *  ( (conductance / (2 * pi * radius * mem1))  * (sum [old-voltage] of link-neighbors) - ( count link-neighbors * ( conductance / (2 * pi * radius * mem1) )) * old-voltage
       -(INa + IK + IL + IM) + istim)



    set old-voltage voltage
    set mold m
    set nold n
    set hold h

  ]
  set counter counter - 1
  if (counter = 1 )  [ ask cell stimcell [ set istim 0.0 ] ]
  color-turtle
  set outvol1 [voltage] of cell rec1
  set outvol2 [voltage] of cell rec2
  do-plot
  tick-advance dt
end

to stim

 ask turtle stimcell [set istim stim1mag]
  set counter stimdur / dt



end

;; color the patch based on its temperature
to color-turtle  ;; Patch Procedure
  ask turtles[
    set color scale-color red voltage -65 50 ]
end


to do-plot
  if plot-pen-exists? "outvol1" [
    set-current-plot-pen "outvol1"
    plotxy ticks outvol1
  ]
  if plot-pen-exists? "outvol2" [
    set-current-plot-pen "outvol2"
    plotxy ticks outvol2
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
438
510
1296
552
-1
-1
2.1
1
30
1
1
1
0
0
0
1
0
404
0
15
1
1
1
ticks
30.0

BUTTON
270
15
335
48
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
177
15
250
48
Go Once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
450
10
1285
496
plot 1
NIL
NIL
0.0
1.0
-65.0
50.0
true
false
"" ""
PENS
"outvol1" 1.0 0 -16777216 true "" "plot outvol1"
"outvol2" 1.0 0 -5298144 true "" "plot outvol2"

BUTTON
66
68
129
101
NIL
stim
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
47
15
163
48
Initialize Cell
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
44
259
364
292
stim1mag
stim1mag
0
2000
2000.0
1
1
NIL
HORIZONTAL

INPUTBOX
137
68
214
128
stimcell
0.0
1
0
Number

INPUTBOX
230
483
302
543
rec1
37.0
1
0
Number

INPUTBOX
318
485
390
545
rec2
52.0
1
0
Number

INPUTBOX
240
68
320
128
stimdur
0.5
1
0
Number

SWITCH
63
566
199
599
demyelinate?
demyelinate?
0
1
-1000

INPUTBOX
217
564
303
624
site1
10.0
1
0
Number

INPUTBOX
314
564
397
624
site2
11.0
1
0
Number

INPUTBOX
415
565
486
625
site3
12.0
1
0
Number

INPUTBOX
502
565
577
625
site4
13.0
1
0
Number

INPUTBOX
61
482
210
542
Capacitance_Factor
340.0
1
0
Number

@#$#@#$#@
## WHAT IS IT?

This model simulates transient and steady-state temperature distribution of a thin plate.

The View shows a square thin plate as viewed from above.  The plate is thermally isolated on the two faces parallel to the view such that heat can flow only in and out from the perimeter of the plate and not into or out of the world.  Heat is kept constant at the edges.  As the simulation runs, heat is transmitted from warmer parts of the plate to cooler parts of the plate as shown by the varying color of the plate.  Therefore, the temperature of the plate begins to change immediately and possibly differently at different locations, gradually converging to a stable state.  Overall, the temperature distribution over the plate is a function of time and location.  In addition to this simple use of the model, you are encouraged to control various paramaters, such as the temperature of each edge edge of the plate and of the center of the plate before--and even while--the model is running.

Heat diffuses ("spreads") at different rates through different media.  These rates can be determined and are called the Thermal Diffusivity of the material.  The Greek letter alpha is often associated with this value.  The diffusivity of a material does not change based on how much of the material there is.  It is always the same.  Below is a table containing several different materials with different diffusivity rates.  See that wood (bottom row) has a lower heat diffusivity than, say, iron.  This means that it takes a longer for heat to spread through a wooden object than an iron one.  That is one reason why the handles of iron saucepans are wooden, and not the other way round.  Also, think of a marble table with iron legs that has just been put out in the sun in a street-side cafe.  Which material part of the table do you expect will warm up faster?  The model allows you to change thermal diffusivity of the plate in two ways.  You can directly change the value of ALPHA to any value you like, or you can indirectly change ALPHA by selecting a material.

### Thermal diffusivity of selected materials

<table border>
<tr><th>Material<th>Thermal diffusivity<br>(alpha cm*cm/s)
<tr><td>Wood (Maple)<td>0.00128
<tr><td>Stone (Marble)<td>0.0120
<tr><td>Iron<td>0.2034
<tr><td>Aluminum<td>0.8418
<tr><td>Silver<td>1.7004
</table>

## HOW IT WORKS

Initialize the plate and edges to have temperatures that equal their respective slider values.  Each time through the GO procedure, diffuse the heat on each patch in the following way.  Have each patch set its current temperature to the sum of the 4 neighbors' old temperature times a constant based on alpha plus a weighted version of the patch's old temperature.  (For those interested, the updated temperature is calculated by using a Forward Euler Method.)  Then the edges are set back to the specified values and the old temperature is updated to the current temperature.  Then the plate is redrawn.

## HOW TO USE IT

There are five temperature sliders which enable users to set four fixed edge temperatures and one initial plate temperature:
-- TOP-TEMP - Top edge temperature
-- BOTTOM-TEMP - Bottom edge temperature
-- IN-PLATE-TEMP - Initial plate temperature
-- LEFT-TEMP - Left edge temperature
-- RIGHT-TEMP - Right edge temperature

There are two sliders that govern the thermal diffusivity of the plate:
-- MATERIAL-TYPE - The value of the chooser is that of the above chart.  You must press UPDATE ALPHA for this to change the value of ALPHA.
-- ALPHA - The alpha constant of thermal diffusivity

There are four buttons with the following functions:
-- SETUP - Initializes the model
-- GO - Runs the simulation indefinitely
-- GO ONCE - Runs the simulation for 1 time step
-- UPDATE ALPHA - press this if you want to set ALPHA to a preset value based on a material selected by the MATERIAL-TYPE chooser

The TIME monitor shows how many time steps the model has gone through.

## THINGS TO NOTICE

How does the equilibrium temperature distribution vary for different edge temperature settings?

Notice how an equilibrium (the steady-state condition) is reached.

Keep track of the units:

<table border>
<tr><th>Variables<th>Units
<tr><td>time<td>0.1 second
<tr><td>temperature<td>degrees Celsius
<tr><td>length<td>centimeters
<tr><td>diffusivity<td>square centimeters per second
</table>

## THINGS TO TRY

Set the parameters on the temperature sliders.  Pick a value for ALPHA (or pick MATERIAL-TYPE and press UPDATE ALPHA).  After you have changed all the sliders to values you like, press Setup followed by GO or GO ONCE.

Try different materials to observe the heat transfer speed.  How does this compare to physical experiments?

Try the following sample settings:
- Top:100, Bottom:0,   Left:0,   Right:0
- Top:0,   Bottom:100, Left:100, Right:100
- Top:0,   Bottom:66,  Left:99,  Right:33
- Top:25,  Bottom:25,  Left:100, Right:0

## EXTENDING THE MODEL

This model simulates a classic partial differential equation problem (that of heat diffusion). The thin square plate is a typical example, and the simplest model of the behavior.  Try changing the shape or thickness of the plate (e.g. a circular or elliptical plate), or adding a hole in the center (the plate would then be a slice of a torus, a doughnut-shaped geometric object).

Add a slider to alter this thickness.

Try modeling derivative or combined boundary conditions.

## RELATED MODELS

Heat Diffusion - Alternative Gradient

## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Wilensky, U. (1998).  NetLogo Heat Diffusion model.  http://ccl.northwestern.edu/netlogo/models/HeatDiffusion.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 1998 Uri Wilensky.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.

This model was created as part of the project: CONNECTED MATHEMATICS: MAKING SENSE OF COMPLEX PHENOMENA THROUGH BUILDING OBJECT-BASED PARALLEL MODELS (OBPML).  The project gratefully acknowledges the support of the National Science Foundation (Applications of Advanced Technologies Program) -- grant numbers RED #9552950 and REC #9632612.

This model was converted to NetLogo as part of the projects: PARTICIPATORY SIMULATIONS: NETWORK-BASED DESIGN FOR SYSTEMS LEARNING IN CLASSROOMS and/or INTEGRATED SIMULATION AND MODELING ENVIRONMENT. The project gratefully acknowledges the support of the National Science Foundation (REPP & ROLE programs) -- grant numbers REC #9814682 and REC-0126227. Converted from StarLogoT to NetLogo, 2001.

<!-- 1998 2001 -->
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
