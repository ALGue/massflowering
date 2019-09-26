extensions[r profiler] ; r = to compute Beta / exponential distribution with R software ; profiler = check time execution

__includes[
  "landscape.nls"
  "cbcSetup.nls"
  "cbcGo.nls"
  "pollSetup.nls"
  "pollGo.nls"]

globals [

  stockCrops ; tot. number of crop patches
  ; patches groups according to the landcover
  cropPatches
  snhPatches
  ; crop groups according to crop types (defined with cropId)
  cropsA
  cropsB
  ; service requirements according to crop type
  cropsRequiringCBC
  cropsRequiringPollination

  y ; year
  d ; date
  s ; length of growing season
  s1 ; transition date between the 2 successive periods of cbc, and then pollination

  ; K ; carrying capacity / overwintering effect = nb of SNH patches required for the survival of 1 adult
  SNHclusterId_List ; list of unique snh cluster id
  SNHclusterCapacity_List ; list of the number of agents each snh cluster can sustain during overwintering (adults or broods)

  fileName

  ; visitation rate
  visitationRate

  ;delay
  meanDelay

  ; crop loss
  averageCropLoss

  ; final yield
  averageFinalYield
  averageYieldCropsCBC
  averageYieldCropsPollination

  ; total nb of foragers after emergence -> mortality due to dispersal into crops for nesting
  totalForagersAfterEmergence

  ; total nb of foragers end year -> mortalityRate
  totalForagersEndYear

  ; total Broods after conversion of stored pollen -> conversionFactor
  totalBroodsAfterConversion

  ; total Broods after overwintering -> carryingCapacity x SNH
  totalBroodsAfterOverwintering

  ; total pollen resource collected by pollinators and stored into nests
  totalPollenCollected

  ; total initial pollen stock in crops (in the whole landscape) = total pollen before emergence of pollinators -> snh x crop loss
  totalInitPollenStock

  ; total pollen stock in crops (in the whole landscape) at the end of the year -> totalInitPollenStock x totalPollenCollected
  totalFinalPollenStock

  ; total number of flights by pollinator foragers
  totFlights

  ; csq of crop loss on pollen supply for pollinators => 0.7 (1 -> 0.3) ; 0.8 (1 -> 0.2)
  ; deltaPollenSupplied
]

patches-own
[
  landCover ; crops 0 or SNH 1

  cropId ; crop specie 1 2 3

  xClosestSNH
  yClosestSNH
  closestSNHDistance

  clusterSeed
  clusterId

  ; Robert-Addis
  t ; day for beginning of flowering period (beta distribution)
  nBroods ; b ; number of brood developing in cell i
  nForagers ; n ; number of foragers nesting within cell i
  f ; boolean flowering status of the plant in cell i
  p ; boolean pollination status of the plant in cell i

  ;pollenSupplied

  maxPollen
  ; minPollen
  currentStockPollen

  nbVisits

  myPollinationDependence

  ; CBC
  free?
  nextState
  state
  d1
  d2
  delay
  myCropLoss

  ; yield = f°(cbc, poll)
  myFinalYield


]

breed [broods brood]
broods-own [emergenceDate]

breed [foragers forager]
foragers-own [lifespan myNest brood? qPollenCollected myCropDestination-x myCropDestination-y]

breed [adults adult]
adults-own [
  birthyear
  emergenceDate
  move?
  timer_move] ; timer_move = makes possible to limit the mov frequency


;;;;;;;;;;;;;;;;;;;;;;; MAIN

to setup

  clear-all
  reset-ticks

  ; fileName for outputs => see saveES
  set fileName (word
    (random 1000)
    "-"
    snh
    "-"
    (remove "."(word targetForAgregation))
    "-"
    dispersal_behaviour
    "-"
    configurationCrops
    "-"
    ratio_A/B
    "-"
    square_side
    "-"
    choiceCropsRequiringCBC
    "-"
    choiceCropsRequiringPollination
    "-"
    funpollenLoad
    "-"
    (remove "."(word scaleDispersal))
    "-"
    (remove "."(word deltaForagers))
    "-"
    (remove "."(word conversionFactor))
    "-"
    ;(remove "."(word collectionCapacity))
    ;"-"
    (remove "."(word maxPollenSlider))
    "-"
    minPollen
    "-"
    nbVisitsToDepletePollen
    "-"
    (remove "."(word pcr))
    "-"
    (remove "."(word deltaAdults))
    "-"
    (remove "."(word K))
    "-"
    (remove "."(word freq_mov_ad))
    "-"
    (remove "."(word repro_rate))
    "-"
    (remove "."(word perception_range))
    "-"
    %adults
    "-"
    %broods
    "-"
    t2
    "-"
    alpha
    "-"
    (random 1000)".txt")

  ;;; landscape initialization

  ; assign crops and snh
  landscape
  set cropPatches patches with [landCover = 0]
  ask cropPatches[set pcolor 9.9]
  set snhPatches patches with [landCover = 1]
  ask snhPatches [set pcolor green]
  set stockCrops count cropPatches ; var for cbcGo

  ; assign cropId "A" and "B" to cropPatches
  assignCropId
  set cropsA cropPatches with [cropId = "A"]
  set cropsB cropPatches with [cropId = "B"]

  ; assign services requirements for crop types with the help of cropId
  ; écrire condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if choiceCropsRequiringCBC = "cropsA" [set cropsRequiringCBC cropsA]
  if choiceCropsRequiringCBC = "cropsB" [set cropsRequiringCBC cropsB]
  if choiceCropsRequiringCBC = "crops A & B" [set cropsRequiringCBC cropPatches]
  if choiceCropsRequiringCBC = "-" [set cropsRequiringCBC cropPatches with [cropId != "A" and cropId != "B"]]

  if choiceCropsRequiringPollination = "cropsA" [set cropsRequiringPollination cropsA]
  if choiceCropsRequiringPollination = "cropsB" [set cropsRequiringPollination cropsB]
  if choiceCropsRequiringPollination = "crops A & B" [set cropsRequiringPollination cropPatches]
  if choiceCropsRequiringPollination = "-" [set cropsRequiringPollination  cropPatches with [cropId != "A" and cropId != "B"]]

  ; initialise cro
  ask cropsRequiringPollination [set myPollinationDependence 0.8]

  ; list 1/ ids and 2/ capacities for overwintering of snh clusters
  set SNHclusterId_List remove-duplicates [clusterid] of snhPatches
  set SNHclusterCapacity_List map [x -> floor(count patches with [clusterid = x] / K)] SNHclusterId_List


  ; parameter initialization (time)
  set y 1
  set d 0
  set s 180 ; 31 ; 180
  set s1 150 ; 1 ; 150 ; transition date between the 2 successive periods of cbc, and then pollination

  ; pollination initialization
  pollSetup

  ; cbc initialization
  cbcSetup

end

to go
  ifelse d >= s
  ; transition between years
  [
    computeFinalYield
    overwinteringEnemies
    overwinteringBroods
    ;saveDynamics
    if y >= (totYears - 3) [if save? [saveES]] ; save outputs
    transition ; update patches for next year
    if y = totYears [stop]
    set y (y + 1)
    set d 0
    ask cropPatches [set pcolor 9.9]
    resetIndicators
]
  ; processes during the year
  [
    ;saveDynamics
    tick
    set averageFinalYield 0
    set d d + 1
    cbcGo
    if d = s1 [
      ask cropPatches [set pcolor 9.9]
      computePollenSuppAfterCBC
      computeCBCindicators]
    if (d >= s1) and (d <= s) [pollGo]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;; SOME USEFUL METHODS

to profile
  profiler:start         ;; start profiling
  setup
  ; repeat 180 * 10 [ go ]    ;; run something you want to measure => 75 * 180 = 13500 => *10 time for go to have the total time of a complete simulation
  repeat 180 * 5 [ go ]
  profiler:stop          ;; stop profiling
  print profiler:report  ;; view the results
  profiler:reset         ;; clear the data
end

to saveES

  ; directory for outputs
  set-current-directory (word "/home/antoine/Documents/Git/Poll_CBC2/outputs/" folderES)

  file-open fileName

  ask cropPatches
  [file-print (word
    cropId
    " "
    y
    " "
    ;deltaPollenSupplied
    ;" "
    d1
    " "
    d2
    " "
    p
    " "
    myCropLoss
    " "
    myFinalYield
    " "
    ; total quantity of pollen initially in crops (in the whole landscape), before emergence of pollinators -> snh x crop loss
    totalInitPollenStock
    " "
    ; total quantity of pollen in crops (in the whole landscape) at the end of the year -> totalInitPollenStock x totalPollenCollected
    totalFinalPollenStock
    " "
    ; total pollen resource collected by pollinators and stored into nests
    totalPollenCollected
    " "
    ; total nb of foragers after emergence -> mortality due to dispersal into crops for nesting
    totalForagersAfterEmergence
    " "
    ; total nb of foragers end year -> mortalityRate
    totalForagersEndYear
    " "
    ; total Broods after conversion of stored pollen -> conversionFactor
    totalBroodsAfterConversion
    " "
    ; total Broods after overwintering -> carryingCapacity x SNH
    totalBroodsAfterOverwintering
    " "
    ; total nb of flights by pollinators
    totFlights
  )]

  file-close

end

to saveDynamics
  ; directory for outputs
  set-current-directory (word "/home/antoine/Documents/Git/Poll_CBC2/outputs/" folderDynamics)

  file-open fileName

  file-print (word
    y
    " "
    d
    " "
    sum [currentStockPollen] of patches with [nBroods > 0]
    " "
    count foragers
  )

  file-close
end

to computeCBCindicators

  ; visitation rate by NE of crops colonised by pests
  set visitationRate count patches with [landcover = 0 and d1 < 150 and d2 < 180] * 100 / stockCrops
  ; mean delay between pest colonisation and NE arrival on crop patches
  ask cropPatches with [d1 < 150 and d2 < 180][set delay (d2 - d1)]
  if visitationRate > 0 [
    set meanDelay mean [delay] of cropPatches with [d1 < 150 and d2 < 180]]
  ; average crop losses on crop patches after regulation
  ;ask cropPatches[set myCropLoss (maxCL d1) - (reg d1 d2) / 100 * (maxCL d1)]
  ask cropPatches [
    ifelse member? self cropsRequiringCBC
    [set myCropLoss (maxCL d1) - (reg d1 d2) / 100 * (maxCL d1)]
    [set myCropLoss 0]
  ]
  set averageCropLoss mean [myCropLoss] of cropPatches

end

to resetIndicators
  set visitationRate 0
  set meanDelay 0
  set averageCropLoss 0
  set averageYieldCropsCBC 0
  set averageYieldCropsPollination 0
  set totFlights 0

end

to-report BetaDistrib [a b]
  let varX random-gamma a 1
  let varY random-gamma b 1
  report varX / (varX + varY)
end

to transition
  ; cbc
  ask patches [
    set free? TRUE
    set state 0
    set d1 s
    set d2 s
    set myCropLoss 0
    set myFinalYield 0
  ]

  ; pollination

  ask cropPatches[
    set nForagers 0
    set f 0
    ;set currentStockPollen maxPollenSlider
    set p 0
    ; set t round (s1 + (BetaDistrib alpha beta) * (s - s1)) ; new flowering date
    set t s1
    set pcolor 9.9] ; update color of pollinated crops (yellow -> white)
  ask snhPatches[
    set nForagers 0
    set f 0
    set p 0
    set t s1
  ]
    ; set t round (s1 + (BetaDistrib 1 1) * (s - s1))] ; new flowering date
end

to newAdults
  sprout-adults 1 [
    set birthyear y
    set emergenceDate 1
    set size 0.5
    set color 0
    set move? TRUE
    set timer_move 0]
end

to newBroods
  sprout-broods 1 [
    ; set emergenceDate round (s1 + (BetaDistrib mu eta) * (s - s1)) ; emergenceDate between 150 and 180
    set emergenceDate s1
    set shape "square"
    set size 0.5
    set color 0]
end

to overwinteringEnemies

  ;;; cbc

  ; enemies still in crops die
  ask adults [if [landCover] of patch-here = 0 [die]]

  ;show "before"
  ;show count adults

  ; list nb of adult enemies on each snh cluster
  let nbAdultsOnClusters_List map [x -> count adults-on patches with [clusterid = x]] SNHclusterId_List

  ; all enemies die
  ask adults [die]

  ; foreach cluster, sprout randomly the given nb of enemies after comparison between refugees and capacity
  (foreach SNHclusterCapacity_List nbAdultsOnClusters_List SNHclusterId_List [ [capacity nb id] -> ifelse nb >= capacity [ask n-of capacity patches with [clusterId = id][newAdults]][ask n-of nb patches with [clusterId = id][newAdults]]])

  ;show "after"
  ;show count adults

end

to overwinteringBroods
  ;;; pollination

  ; total quantity of pollen in crops at the end of the year -> totalInitPollenStock x totalCollectedPollen
  set totalFinalPollenStock sum [currentStockPollen] of cropPatches

  ; all foragers die, only broods survive
  set totalForagersEndYear count foragers
  ask foragers [die]

  ask patches [set nBroods 0]

  ; compute how many broods emerge from the pollen load in SNH patches
  conversionPollenBroods
  set totalBroodsAfterConversion count broods

  ; list nb of broods on each snh cluster
  let nbBroodsOnClusters_List map [x -> count broods-on patches with [clusterid = x]] SNHclusterId_List

  ; all broods die
  ask broods [die]

  ; foreach cluster, sprout randomly the given nb of broods after comparison between refugees and capacity
  (foreach SNHclusterCapacity_List nbBroodsOnClusters_List SNHclusterId_List [ [capacity nb id] -> ifelse nb >= capacity [ask n-of capacity patches with [clusterId = id][newBroods]][ask n-of nb patches with [clusterId = id][newBroods]]])

  ; total Broods after overwintering -> carryingCapacity x SNH
  set totalBroodsAfterOverwintering count broods

end

to computePollenSuppAfterCBC

  ; calibration:  pollenSupplied = 0.2 - 1, according to crop loss (cl)

  ifelse funPollenLoad = "set" [ask patches with [landCover = 0][
    set maxPollen maxPollenSlider
  set currentStockPollen maxPollen]]
    [if funPollenLoad = "cl_linear"[
      ifelse pcr > 0 [
        ask cropPatches[
        ifelse member? self cropsRequiringCBC[
          let cropLoss ((maxCL d1) - (reg d1 d2) / 100 * (maxCL d1))
          set maxPollen 1 - cropLoss * (1 - minPollen) / 100][set maxPollen 1 set currentStockPollen maxPollen]
        set currentStockPollen maxPollen]] ; deltaPollenSupplied = 0.7
    [ask patches with [landCover = 0][set maxPollen 1 set currentStockPollen maxPollen]]]]

  ;  if funPollenLoad = "cl_convexe"[
;      ifelse pcr > 0 [
;        ask cropPatches[
;            let n 0.5
;            let cropLoss ((maxCL d1) - (reg d1 d2) / 100 * (maxCL d1))
;            set maxPollen (- (1 - minPollen) / 100 ^ (n)) * cropLoss ^ (n) + 1]] ; deltaPollenSupplied = 0.7
;      [ask patches with [landCover = 0][set maxPollen 1]]]
;  if funPollenLoad = "cl_concave"[
;      ifelse pcr > 0 [
;        ask cropPatches[
;            let n 2
;            let cropLoss ((maxCL d1) - (reg d1 d2) / 100 * (maxCL d1))
;            set maxPollen (- (1 - minPollen) / 100 ^ (n)) * cropLoss ^ (n) + 1]] ; deltaPollenSupplied = 0.7
;      [ask cropPatches[set maxPollen 1]]]]

  ; total quantity of pollen initially in the landscape, before emergence of pollinators -> snh x crop loss
  set totalInitPollenStock sum [currentStockPollen] of cropPatches

end

to-report maxCL [a] ; finalCL <- function(date){90 * exp(-(1 / 90) * (date - 1))}
  ifelse (a > 0 and a <= 180) [report 90 * exp(-(1 / 90) * (a - 1))][report 0]
end

to-report reg [x1 x2] ; finalReg3 <- function(date){5 + (95-5) / (1 + exp(-(-0.3)*(date - 15)))^1}
  ifelse  x2 != 180 [report (5 + (95 - 5) / (1 + exp(-(-0.3)*((x2 - x1 - 7) - 15)))^(1))][report 0] ; attention à modifier lorsque l'on passe à 180
  ; condition (x2 - x1) < 40 and
  ; report 95 * exp(- (x2 - x1) / 60)
end

to conversionPollenBroods

  ; total quantity of pollen collected and stored by foragers in nesting cells
  set totalPollenCollected sum [currentStockPollen] of snhPatches

  ask snhPatches with [currentStockPollen > 0][
    sprout-broods round (currentStockPollen * conversionFactor) [
        set emergenceDate round (s1 + (BetaDistrib mu eta) * (s - s1))
        set shape "square"
        set size 0.5
        set color 0]
    set nBroods round (currentStockPollen * conversionFactor)
  set currentStockPollen 0]

end

to computeFinalYield

;  ask cropsRequiringPollination [
;    let V1 (100 - myCropLoss) ; 60% ; 100
;    let V2 (100 - myCropLoss) * myPollinationDependence ; 60% * 0.8 = 48 ; 80
;    let V3 V1 - V2 ; 60% - 48% = 12% ; 20
;     if p = 0 [
;      let V4 0 * V2
;      set myFinalYield V3 + V4]
;    if p = 0.5 [
;      let V4 0.5 * V2
;      set myFinalYield V3 + V4]
;    if p = 1 [
;      let V4 1 * V2
;      set myFinalYield V3 + V4]]

  ; ask cropsRequiringCBC[set myFinalYield (100 - myCropLoss)]

  ask cropPatches[
    if member? self cropsRequiringCbc and member? self cropsRequiringPollination = FALSE [set myFinalYield (100 - myCropLoss)]
    if member? self cropsRequiringCbc = FALSE and member? self cropsRequiringPollination [
      if p = 0 [
        set myFinalYield 100 - (myPollinationDependence * 100) ]
      if p = 0.5 [
        set myFinalYield 100 - (myPollinationDependence * 100) + 0.5 * (myPollinationDependence * 100)]
      if p = 1 [
        set myFinalYield 100]]
    if member? self cropsRequiringCbc and member? self cropsRequiringPollination [
      let V1 (100 - myCropLoss) ; 60% ; 100
      let V2 (100 - myCropLoss) * myPollinationDependence ; 60% * 0.8 = 48 ; 80
      let V3 V1 - V2 ; 60% - 48% = 12% ; 20
      if p = 0 [
        let V4 0 * V2
        set myFinalYield V3 + V4]
      if p = 0.5 [
        let V4 0.5 * V2
        set myFinalYield V3 + V4]
      if p = 1 [
        let V4 1 * V2
        set myFinalYield V3 + V4]] ]



  set averageFinalYield mean [myFinalYield] of cropPatches
  ;set averageYieldCropsCBC mean [myFinalYield] of cropsRequiringCBC
  ;set averageYieldCropsPollination mean [myFinalYield] of cropsRequiringPollination
end
@#$#@#$#@
GRAPHICS-WINDOW
371
10
808
448
-1
-1
13.0
1
10
1
1
1
0
1
1
1
-16
16
-16
16
0
0
1
ticks
30.0

BUTTON
823
12
896
45
NIL
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

BUTTON
827
62
890
95
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
14
10
186
43
snh
snh
0
100
20.0
1
1
NIL
HORIZONTAL

CHOOSER
630
788
768
833
alpha
alpha
1 100
0

SLIDER
617
837
789
870
t2
t2
1
s - s1
13.0
1
1
NIL
HORIZONTAL

INPUTBOX
827
153
886
213
totYears
10.0
1
0
Number

SWITCH
811
220
901
253
save?
save?
1
1
-1000

BUTTON
821
110
898
143
NIL
profile
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
904
10
1425
155
pest control
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"pests" 1.0 0 -1604481 true "" "plot (count patches with [landCover = 0 and state = 1]) / stockCrops"
"adults" 1.0 0 -16777216 true "" "plot (count adults / count patches)"
"juveniles" 1.0 0 -11221820 true "" "plot (count patches with [landCover = 0 and state = 3]) / stockCrops"
"immune" 1.0 0 -955883 true "" "plot (count patches with [landCover = 0 and state = 4]) / stockCrops"

SLIDER
13
50
187
83
targetForAgregation
targetForAgregation
0
10
10.0
1
1
NIL
HORIZONTAL

SLIDER
14
512
186
545
pcr
pcr
0
8
8.0
0.1
1
NIL
HORIZONTAL

PLOT
1440
10
1961
155
pollination 
NIL
NIL
0.0
10.0
0.0
0.3
true
false
"" ""
PENS
"flowers" 1.0 0 -2674135 true "" "plot count patches with [f = 1] / count patches"
"foragers" 1.0 0 -13791810 true "" "plot count foragers / count patches"
"pollCrops" 1.0 0 -955883 true "" "plot count cropPatches with [p = 1] / stockCrops"
"broods" 1.0 0 -16777216 true "" "plot count broods / count patches"
"pollRate" 1.0 0 -5825686 true "" "plot count cropPatches with [p = 1] / stockCrops * 100"

SLIDER
13
384
185
417
%broods
%broods
0
100
100.0
10
1
NIL
HORIZONTAL

SLIDER
14
345
186
378
%adults
%adults
0
100
100.0
10
1
NIL
HORIZONTAL

INPUTBOX
16
783
238
843
folderES
26_09_1
1
0
String

CHOOSER
403
881
585
926
funPollenLoad
funPollenLoad
"set" "cl_linear" "cl_convexe" "cl_concave"
1

SLIDER
404
805
582
838
maxPollenSlider
maxPollenSlider
0
1
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
403
843
583
876
deltaForagers
deltaForagers
0
0.5
0.05
0.05
1
NIL
HORIZONTAL

PLOT
905
422
1427
542
delay
NIL
NIL
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"delay" 1.0 0 -16777216 true "" "plot meanDelay - 7"

SLIDER
405
767
581
800
scaleDispersal
scaleDispersal
1
10
3.0
1
1
NIL
HORIZONTAL

SLIDER
15
605
187
638
freq_mov_ad
freq_mov_ad
0
5
0.0
1
1
NIL
HORIZONTAL

SLIDER
16
642
188
675
K
K
1
5
2.0
1
1
NIL
HORIZONTAL

SLIDER
14
475
186
508
deltaAdults
deltaAdults
0
1
0.1
0.01
1
NIL
HORIZONTAL

PLOT
905
161
1426
287
nb adult NE
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count adults"

SLIDER
15
680
187
713
repro_rate
repro_rate
0
8
1.0
1
1
NIL
HORIZONTAL

CHOOSER
15
722
188
767
perception_range
perception_range
1.5 2.9 4
0

TEXTBOX
404
749
554
767
Pollination
12
0.0
1

TEXTBOX
17
458
167
476
cbc\n
12
0.0
1

PLOT
905
548
1428
668
averageCropLoss
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot averageCropLoss"

PLOT
905
294
1426
414
visitationRate
NIL
NIL
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot visitationRate"

SLIDER
403
932
587
965
minPollen
minPollen
0
1
0.3
0.1
1
NIL
HORIZONTAL

SLIDER
14
94
186
127
ratio_A/B
ratio_A/B
0
100
50.0
1
1
NIL
HORIZONTAL

CHOOSER
13
140
185
185
configurationCrops
configurationCrops
"circle" "square" "random"
1

SLIDER
13
191
185
224
square_side
square_side
1
9
9.0
1
1
NIL
HORIZONTAL

CHOOSER
13
552
186
597
dispersal_behaviour
dispersal_behaviour
"random" "oriented"
1

PLOT
1439
161
1961
286
pollination cropId A/B
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"pollRateA" 1.0 0 -16777216 true "" "plot count cropPatches with [cropId = \"A\" and p = 1] / count cropsA * 100"
"pollRateB" 1.0 0 -5298144 true "" "plot count cropPatches with [cropId = \"B\" and p = 1] / count cropsB * 100"

PLOT
340
451
862
571
average final yield
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"finalYield" 1.0 0 -16777216 true "" "plot averageFinalYield"

SLIDER
608
891
794
924
conversionFactor
conversionFactor
0
1
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
607
930
828
963
nbVisitsToDepletePollen
nbVisitsToDepletePollen
0
10
2.0
1
1
NIL
HORIZONTAL

INPUTBOX
16
851
244
911
folderDynamics
28_08_5_D
1
0
String

PLOT
1486
445
1686
595
nesting resource
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot sum [currentStockPollen] of snhPatches"

CHOOSER
14
230
224
275
choiceCropsRequiringCBC
choiceCropsRequiringCBC
"cropsA" "cropsB" "crops A & B" "-"
2

CHOOSER
14
280
271
325
choiceCropsRequiringPollination
choiceCropsRequiringPollination
"cropsA" "cropsB" "crops A & B" "-"
0

PLOT
392
578
592
728
averageYieldCropsCBC
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"cbc" 1.0 0 -16777216 true "" "plot averageYieldCropsCBC"

PLOT
599
578
799
728
averageYieldCropsPollination
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Poll" 1.0 0 -16777216 true "" "plot averageYieldCropsPollination"

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

; INTERFACE PARAMETERS

;;; landscape + time + save
; snh = proportion of SNH
; targetForAgregation = agregation of SNH => 0 desagregated - 10 agregated
; totYears = nb years simulated
; save? = TRUE / FALSE => records results for ES
;;; pollination
; alpha = flowering concentration parameter (crops) => 1 : uniform, 100 : peak
; t2 = peak emergence date for pollinator foragers (emergence of broods)
;;; cbc
; pcr = pest colonisation rate

; OTHER PARAMETERS

;;; landscape + time + save
; s
; s1

; K = carrying capacity of snh ; nb of cells necessary for survival during overwintering (for adults and broods)

;;; pollination
; l = length flowering period
; t = peak flowering date
; beta
; mu
; eta
; nu
; deltaForagers = daily mortality rate for foragers
; conversionFactor = conversion factor between stored resource and new broods emerging
;;; cbc
; emergenceDate
; deltaAdults = daily mortality rate for adults
; distDispersalAdults = radius dispersal

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

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

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="cbcFloriane" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <steppedValueSet variable="snh" first="0" step="5" last="100"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="10_06" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="pollcbcNointeraction" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="50"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="1"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="20_06" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="10"/>
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="1"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;20_06&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="deltaForagers" first="0.036" step="0.05" last="0.2"/>
    <steppedValueSet variable="qPollenVariability" first="0.2" step="0.4" last="1"/>
  </experiment>
  <experiment name="21_06" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="40"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="1"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;21_06&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.086"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.4"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="21_06_2" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="40"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="1"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;21_06_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="21_06_3" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="40"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;21_06_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="21_06_4" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="40"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;21_06_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="22_06" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="40"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;22_06&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="22_06_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="40"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;22_06_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_06" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="30"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_06&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_06_1" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_06_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_06_4" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_06_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_06_4" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_06_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="29_06" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;29_06&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="29_06_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;29_06_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="01_07" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;01_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
      <value value="0.6"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="01_07_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;01_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.1"/>
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
      <value value="0.6"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="01_07_2" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;01_07_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="02_07" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;02_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
      <value value="0.4"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.3"/>
      <value value="0.5"/>
      <value value="0.7"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="02_07_1" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="10"/>
      <value value="15"/>
      <value value="20"/>
      <value value="25"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;02_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="02_07_3" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="10"/>
      <value value="15"/>
      <value value="20"/>
      <value value="25"/>
      <value value="30"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;02_07_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="05_07" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="40" step="5" last="60"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;05_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="05_07_2" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;05_07_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="05_07_3" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="30"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;05_07_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="05_07_4" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="35" step="5" last="60"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;05_07_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="06_07" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;06_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
      <value value="2"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="07_07" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;07_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
      <value value="2"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="07_07_1" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;07_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="07_07_2" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;07_07_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="08_07_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;08_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="08_07_2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;08_07_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="09_07" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;09_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
      <value value="4"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="10_07" repetitions="3" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;10_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="15_07" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="15"/>
      <value value="30"/>
      <value value="45"/>
      <value value="60"/>
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;15_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="15_07_1" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="55"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;15_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="15_07_2" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="55"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;15_07_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="15_07_3" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;15_07_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.5"/>
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="16_07" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;16_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.25"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;user&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="15_07_5" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="35"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;15_07_5&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="15_07_6" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="35"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;15_07_6&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="15_07_6" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="35"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;15_07_7&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="16_07_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="35"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;16_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="16_07_2" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;16_07_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="16_07_3" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="65"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;16_07_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="16_07_3_bis" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;16_07_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="16_07_4" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="15"/>
      <value value="25"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;16_07_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_07" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="70"/>
      <value value="75"/>
      <value value="80"/>
      <value value="85"/>
      <value value="90"/>
      <value value="95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;18_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_07_1" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="10"/>
      <value value="15"/>
      <value value="25"/>
      <value value="30"/>
      <value value="40"/>
      <value value="45"/>
      <value value="55"/>
      <value value="60"/>
      <value value="70"/>
      <value value="75"/>
      <value value="85"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;18_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_07_2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="15"/>
      <value value="25"/>
      <value value="35"/>
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;18_07_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
      <value value="&quot;cl_convexe&quot;"/>
      <value value="&quot;cl_concave&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_07_3" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;18_07_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;set&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="1"/>
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_07_4" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;18_07_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_concave&quot;"/>
      <value value="&quot;cl_convexe&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.6"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_07_5" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="5"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;18_07_5&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
      <value value="&quot;cl_concave&quot;"/>
      <value value="&quot;cl_convexe&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_07_6" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;18_07_6&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;convexe&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_07" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;19_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;set&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="qPollenVariability">
      <value value="1"/>
      <value value="0.5"/>
      <value value="0.3"/>
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="21_07" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;random&quot;"/>
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;21_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="21_07_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="95"/>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;square&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="square_side">
      <value value="1"/>
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;21_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="22_07" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;random&quot;"/>
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;22_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="22_07_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="90"/>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;square&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="square_side">
      <value value="1"/>
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;22_07_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="23_07" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;square&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="square_side">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;23_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="24_07" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;square&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="square_side">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;24_07&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="26_08_2" repetitions="19" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="10"/>
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;26_08_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;set&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="26_08_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;26_08_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="27_08" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;27_08&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaPollenSupplied">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_08_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_08_1&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_08_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_08_2&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_08_3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_08_3&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_08_4" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folder">
      <value value="&quot;28_08_3&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_08_4" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;28_08_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_4_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_08_5" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;28_08_5&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
      <value value="0.1"/>
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="28_08_6" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;28_08_6&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
      <value value="0.1"/>
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="29_08" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;29_08&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="29_08_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;29_08_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="1"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="29_08_2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;29_08_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="29_08_3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;29_08_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="29_08_4" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;29_08_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="30_08" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;30_08&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="30_08_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;30_08_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="30_08_2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;30_08_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="30_08_3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;30_08_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="3"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="03_09" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;03_09&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="7"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="03_09_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;03_09_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="03_09_2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;03_09_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="04_09_1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;04_09_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="0"/>
      <value value="15"/>
      <value value="45"/>
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="3"/>
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.025"/>
      <value value="0.05"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="1"/>
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.25"/>
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="04_09_2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;04_09_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="04_09_3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;04_09_3&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="04_09_4" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;04_09_4&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderDynamics">
      <value value="&quot;28_08_5_D&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="1"/>
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="06_09" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;06_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="13_09" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;13_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t2">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="17_09" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;17_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_09" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;18_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_09_1_cbconly" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;18_09_1&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_09_2_cbconly_K2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;18_09_2&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_09_3_cbconly_noNE" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;18_09_3&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_09_4_pestcolonisationadapted" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;18_09_4&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="18_09_5_cbconly" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;18_09_5&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_cbconly" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_1_cbconly" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_2_cbconly" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_3_cbcpoll" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_3&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_4_cbcpoll_random" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_4&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_5_cbconly" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_5&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snh">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_6_cbcpoll_random" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_6&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="10" last="90"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_7_cbcpoll_random" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_7&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;random&quot;"/>
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_8_C1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_8&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="19_09_10_C1" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;19_09_10&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="20_09_Circle" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;20_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="0" step="5" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;cropsB&quot;"/>
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
      <value value="&quot;cropsB&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="20_09_1_Circle" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;20_09_1&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;cropsB&quot;"/>
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
      <value value="&quot;cropsB&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="24_09_C1cbcpoll" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;24_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="24_09_3_50cbcAB_50pollB" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;24_09_3&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;crops A &amp; B&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsB&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="24_09_4_50cbcAB_50pollA" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;24_09_4&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="5" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;crops A &amp; B&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="24_09_Acbc_Bpoll" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;24_09_1&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsB&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="24_09_Bcbc_Apoll" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;24_09_2&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;cropsB&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="24_09_C1cbcpoll_scaleD3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;24_09_5&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="26_09_AcbcpollBcbc_scaleD3" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;26_09&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;crops A &amp; B&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="26_09_1_AcbcpollBcbc_paramPoll" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="save?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="folderES">
      <value value="&quot;26_09_1&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="snh" first="5" step="10" last="95"/>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="0"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ratio_A/B">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="configurationCrops">
      <value value="&quot;circle&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringCBC">
      <value value="&quot;crops A &amp; B&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="choiceCropsRequiringPollination">
      <value value="&quot;cropsA&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totYears">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%broods">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scaleDispersal">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaForagers">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="funPollenLoad">
      <value value="&quot;cl_linear&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPollenSlider">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minPollen">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nbVisitsToDepletePollen">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="conversionFactor">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pcr">
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adults">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dispersal_behaviour">
      <value value="&quot;oriented&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="freq_mov_ad">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="repro_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perception_range">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="deltaAdults">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
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
