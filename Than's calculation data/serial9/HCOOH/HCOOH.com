%chk=FormicAcid.chk
# opt freq=hpmodes mp2/3-21g geom=connectivity symmetry=loose

Formic acid HCOOH - size2to5_Cs

0 1
 O       -0.776305     -1.638527     -0.703786
 C       -1.014227     -0.406770     -1.179581
 O       -2.068945      0.165751     -0.980208
 H       -1.598267     -1.876082     -0.225482
 H       -0.153985     -0.020238     -1.745299

1 2 1.0 4 1.0
2 3 2.0 5 1.0
3
4
5

