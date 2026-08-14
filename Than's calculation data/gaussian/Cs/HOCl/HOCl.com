%mem=8GB
%nprocshared=8
%chk=HOCl.chk
# opt freq=hpmodes mp2/3-21g* geom=connectivity

HOCl - size1_Cs

0 1
 O        0.000000      0.000000      0.000000
 H        0.939300      0.216900      0.000000
 Cl       0.000000     -1.689000      0.000000

1 2 1.0 3 1.0
2
3

