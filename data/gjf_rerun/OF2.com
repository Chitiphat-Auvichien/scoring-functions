%chk=OF2_MP2_3-21G.chk
# opt freq=hpmodes mp2/3-21g geom=connectivity

of2 vibration

0 1
 O                 -0.00000000    0.00000000   -0.49563300
 F                  0.00000000    1.01245576    0.22028133
 F                  0.00000000   -1.01245576    0.22028133

 1 2 1.0 3 1.0
 2
 3

