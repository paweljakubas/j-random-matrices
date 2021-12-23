NB. --------------------------------
NB. -- Algebra - useful functions --
NB. --------------------------------

NB. Decompose each index in a given array into respective coordinates
toIxs=: 3 : '(#:i.)@$y'
NB. Examples:
   v=: 1 2 3 4 5
   toIxs v
0
1
2
3
4
   ]m=: 3 4 $ i.12
0 1  2  3
4 5  6  7
8 9 10 11
   toIxs m
0 0
0 1
0 2
0 3

1 0
1 1
1 2
1 3

2 0
2 1
2 2
2 3
   ]t=: 3 2 2 $ i.12
 0  1
 2  3

 4  5
 6  7

 8  9
10 11
   toIxs t
0 0 0
0 0 1

0 1 0
0 1 1


1 0 0
1 0 1

1 1 0
1 1 1


2 0 0
2 0 1

2 1 0
2 1 1
