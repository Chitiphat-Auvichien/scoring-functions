%chk=H2O-MP2-321G.chk
# opt freq=hpmodes mp2/3-21g geom=connectivity

h2o vibration

0 1
 O                  0.00000000    0.00000000   -0.11085125
 H                  0.00000000   -0.78383672    0.44340501
 H                  0.00000000    0.78383672    0.44340501

 1 2 1.0 3 1.0
 2
 3

