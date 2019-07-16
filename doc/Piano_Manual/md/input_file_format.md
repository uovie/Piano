# Input File Format

## Piano adopts a special input format

```
[job] [description1] [description2] ...

[run time (fs)] [step size (fs)] [data collection period]

[dimension] [volume] [temperature (K)] [pressure]

[model type] [model parameter1] [model parameter2] ...

*****
[number of atoms]

[atom symbol]    [ Rx ]    [ Ry ]    [ Rz ]
[atom symbol]    [ Rx ]    [ Ry ]    [ Rz ]
[atom symbol]    [ Rx ]    [ Ry ]    [ Rz ]
... ...
```

For instance, a practical Piano input file looks like this:

```
pimd nhc middle

1e6 0.1 100

3 0 14 0

LJ 35.6 2.749

*****
13

Ne    0.00000000    0.00000000    0.00000000
Ne    0.00000000    2.18188274    2.18188274
Ne    0.00000000   -2.18188274    2.18188274
Ne    0.00000000   -2.18188274   -2.18188274
Ne    0.00000000    2.18188274   -2.18188274
Ne    2.18188274    0.00000000    2.18188274
Ne    2.18188274    0.00000000   -2.18188274
Ne   -2.18188274    0.00000000   -2.18188274
Ne   -2.18188274    0.00000000    2.18188274
Ne    2.18188274    2.18188274    0.00000000
Ne   -2.18188274    2.18188274    0.00000000
Ne   -2.18188274   -2.18188274    0.00000000
Ne    2.18188274   -2.18188274    0.00000000
```

