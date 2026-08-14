%mem=8GB
%nprocshared=8
%chk=XeF2Cl2.chk
# opt freq=hpmodes mp2/3-21g* geom=connectivity

trans-XeF2Cl2 - size1_D2h (hypothetical test compound)

0 1
 Xe       0.000000      0.000000      0.000000
 F        1.900000      0.000000      0.000000
 F       -1.900000      0.000000      0.000000
 Cl       0.000000      2.300000      0.000000
 Cl       0.000000     -2.300000      0.000000

1 2 1.0 3 1.0 4 1.0 5 1.0
2
3
4
5

