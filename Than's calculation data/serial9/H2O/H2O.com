%chk=H2O.chk
# opt freq=hpmodes mp2/3-21g geom=connectivity

H2O - size1_C2v

0 1
 O        0.000000      0.000000      0.117300
 H        0.000000      0.757200     -0.469200
 H        0.000000     -0.757200     -0.469200

1 2 1.0 3 1.0
2
3

