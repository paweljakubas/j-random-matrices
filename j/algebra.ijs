NB. --------------------------------
NB. -- Algebra - useful functions --
NB. --------------------------------

NB. How to load the script - an example
NB. ]scriptdir=: 'PATH-TO-REPO/j-random-matrices/j/'
NB. 0!:1 < scriptdir,'algebra.ijs'


NB. Color palletes
pallete1=: (255 154 162, 255 183 178, 255 218 193, 226 240 203, 191 234 215,: 199 206 234)
NB. load 'viemat'
NB. pallete1 viewmat ?20 20$6

NB. Plot (x,y) pairs in rectangle defined with left down corner x1,y1 and right top corner
NB. x2,y2 using xn,yn number of parts
NB.
NB.    ---------------------- (x2,y2)  --
NB.    |                    |          .
NB.    |                    |          yn intervals
NB.    |                    |          .
NB.    ---------------------           --
NB.  (x1,y1)
NB.    {... xn intervals ...}
NB.
NB. In each axis there will be additional interval before and after added.


getIxs=: 4 : 0
xI=.>0{y
yI=.>1{y
'a b'=.x
(yI I. b), xI I. a
)

toPlotMatrix=: 4 : 0
'x1 y1 x2 y2 xn yn'=.y
assert. ((x2 >: x1) *. (y2 >: y1))
assert. ((xn > 0) *. (yn > 0))
A0=.((yn+2), (xn+2)) $ 0
dx=. (x2 - x1) % xn
xI=. x1 + dx*i.(>:xn)
dy=. (y2 - y1) % yn
yI=. y1 + dy*i.(>:yn)
xyI=.(<xI),(<yI)
1 (getIxs&xyI"1 x) } A0
)
NB.    (2 2 $ _2 10 1 3) toPlotMatrix (_1, 2, 5, 7, 2, 10)
NB. 0 0 0 0
NB. 0 0 0 0
NB. 0 1 0 0
NB. 0 0 0 0
NB. 0 0 0 0
NB. 0 0 0 0
NB. 0 0 0 0
NB. 0 0 0 0
NB. 0 0 0 0
NB. 0 0 0 0
NB. 0 0 0 0
NB. 1 0 0 0

NB. Round y to x decimal places
round=: 40&$: : (4 : 0)
 b %~ <.1r2+y*b=. 10x^x
)
NB. Examples
NB.    3 round 1.2222234555
NB. 611r500
NB.    10 round _4.44089e_16
NB. 0
NB.    16 round _4.44089e_16
NB. _1r2500000000000000
NB.    15 round _4.44089e_16
NB. 0

NB. Round y to 10 decimal places
r10=: 10&round

NB. Transpose of a matrix
transpose=: |:
NB.    transpose (2 3 $ 1 2 3 4 5 6)
NB. 1 4
NB. 2 5
NB. 3 6
NB.    transpose (transpose (2 3 $ 1 2 3 4 5 6))
NB. 1 2 3
NB. 4 5 6


NB. Multiplication of matrices x and y
mult=: +/ .*
NB.    ]a=: 2 2 $ 1 2 3 4
NB. 1 2
NB. 3 4
NB.    ]b=: 2 4 $ 1 2
NB. 1 2 1 2
NB. 1 2 1 2
NB.    ]c=: a mult b
NB. 3  6 3  6
NB. 7 14 7 14
NB.    $c
NB. 2 4

NB. Rotation matrix where x defines the anchor diagonal index
NB. and y is vector
G=: 4 : 0
assert. ( (>:x) < #y)
l=.(< (<0),(<x)) { (|: y)
r=.(< (<0),(<(>:x))) { (|: y)
norm=.%: ( (*:l) + (*:r) )
m=.(2 2 $ l, r, (-r), l) % norm
xs=.x,>:x
sel=. (<(<xs),(<xs))
m sel } =/~ (i.#y)
)
NB. Examples:
NB.    3 G y
NB. |assertion failure: G
NB. |       ((>:x)<#y)
NB.    2 G y
NB. 1 0    0   0
NB. 0 1    0   0
NB. 0 0  0.6 0.8
NB. 0 0 _0.8 0.6
NB.    1 G y
NB. 1        0       0 0
NB. 0   0.5547 0.83205 0
NB. 0 _0.83205  0.5547 0
NB. 0        0       0 1
NB.    0 G y
NB.  0.447214 0.894427 0 0
NB. _0.894427 0.447214 0 0
NB.         0        0 1 0
NB.         0        0 0 1

NB. Givens rotations for a given vector y
NB. Returns G1, G2, ... Givens rotations
givens=: 3 : 0
'r c' =. ,"0 $ y
assert ( c = 1 )
assert ( r > 1 )
ix=.<:<:#y
s=: 1, (#y), #y
>0} ((s$,ix G y);y) ] F.: {{ ( ( ((>0}y)&,) @ ([: (s&$ @ ,)  x&G)) ; ]) (r10 (({:>0}y) mult (>1}y))) }} i.ix
)
NB. Examples:
NB. ]vec=: 4 1 $ 1 2 3 4
NB. 1
NB. 2
NB. 3
NB. 4
NB.    givens vec
NB. ┌────────────────────────────────┬───────────────────────┐
NB. │        1         0        0   0│                      1│
NB. │        0         1        0   0│53851648071r10000000000│
NB. │        0         0      0.6 0.8│                      0│
NB. │        0         0     _0.8 0.6│                      0│
NB. │                                │                       │
NB. │        1         0        0   0│                       │
NB. │        0  0.371391 0.928477   0│                       │
NB. │        0 _0.928477 0.371391   0│                       │
NB. │        0         0        0   1│                       │
NB. │                                │                       │
NB. │ 0.182574  0.983192        0   0│                       │
NB. │_0.983192  0.182574        0   0│                       │
NB. │        0         0        1   0│                       │
NB. │        0         0        0   1│                       │
NB. └────────────────────────────────┴───────────────────────┘
NB.    ]P=:mult/ |. givens vec
NB.  0.182574  0.365148 0.547723 0.730297
NB. _0.983192 0.0678064  0.10171 0.135613
NB.         0 _0.928477 0.222834 0.297113
NB.         0         0     _0.8      0.6
NB.    r10 P mult vec
NB. 54772255751r10000000000
NB.                       0
NB.                       0
NB.                       0

NB. Householder reflection working on y vector
NB. Returns u vector and P matrix
householder=: 3 : 0
'r c' =. ,"0 $ y
assert ( c = 1 )
assert ( r > 1 )
norm=. {{ %: +/ y*y }}
ke=. ($ y) $ (norm y),((<:r) ;@# 0)
v=. y - ke
u=. |: (|: % norm) v
p=. (=/~ (i.r)) - ((2 * u) mult |: u)
u;p
)
NB.    householder (4 1 $ 1 2 3 4)
NB. ┌─────────┬──────────────────────────────────────┐
NB. │_0.639307│0.182574  0.365148  0.547723  0.730297│
NB. │ 0.285582│0.365148  0.836886 _0.244671 _0.326227│
NB. │ 0.428372│0.547723 _0.244671  0.632994 _0.489341│
NB. │ 0.571163│0.730297 _0.326227 _0.489341  0.347545│
NB. └─────────┴──────────────────────────────────────┘
NB.
NB.    ]u=: >0}householder (4 1 $ 1 2 3 4)
NB. _0.639307
NB.  0.285582
NB.  0.428372
NB.  0.571163
NB.    ]P=: >1}householder (4 1 $ 1 2 3 4)
NB. 0.182574  0.365148  0.547723  0.730297
NB. 0.365148  0.836886 _0.244671 _0.326227
NB. 0.547723 _0.244671  0.632994 _0.489341
NB. 0.730297 _0.326227 _0.489341  0.347545
NB.
NB.    NB. X - (2*u) mult (u'mult X)
NB.    r10 y - ((2 * u) mult ((|: u) mult y))
NB. 54772255751r10000000000
NB.                       0
NB.                       0
NB.                       0
NB.
NB.    NB. P mult X
NB.    r10 P mult y
NB. 54772255751r10000000000
NB.                       0
NB.                       0
NB.                       0

NB. QR decomposition, takes matrix to be decomposed as y and number of rounding place as x,
NB. and returns 'Hs R Q R1 Q1'
qr=: 4 : 0
'r c' =: ,"0 $ y
h1=: >1 { householder (r, 1) $ (<(<a:),(<0)) { y
rr=: x&round
r1=: rr h1 mult y
dH=: 1, r, r
dR=: 1, r, c
S0=:(dH$,h1);(dR $, r1)
'Hs R'=:S0 ]F..{{( ((>0}y)&,) ; (rr @ (mult"2&(>1}y))) ) @ (dH$,) (rr>1 { householder (((<:^:x)r), 1) $ (}.^:x) (<(<0),(<a:),(<x)) { (>1}y)) (<(<(}.^:x i.r)),(<(}.^:x i.r))) } =/~ (i.r)}}>:i.<:c
Q=:|: mult/ |. Hs
Hs;((r,c) $ ,R);Q;((<(<i.c),(<a:)) { (r,c) $, R );((<(<a:),(<i.c)) { Q)
)


NB.    Examples:
NB.    ]A=: 4 3 $ 1 1 1 1 2 4 1 3 9 1 4 16
NB. 1 1  1
NB. 1 2  4
NB. 1 3  9
NB. 1 4 16
NB.    10 qr A
NB. ┌────────────────────────────┬────────────────────────────────────────┬────────────────────────────┬────────────────────────────────────────┬──────────────────┐
NB. │0.5       0.5       0.5  0.5│2                   5                 15│0.5  _0.67082  0.5  0.223607│2                   5                 15│0.5  _0.67082  0.5│
NB. │0.5       0.5      _0.5 _0.5│0 894427191r400000000 894427191r80000000│0.5 _0.223607 _0.5  _0.67082│0 894427191r400000000 894427191r80000000│0.5 _0.223607 _0.5│
NB. │0.5      _0.5       0.5 _0.5│0                   0                  2│0.5  0.223607 _0.5   0.67082│0                   0                  2│0.5  0.223607 _0.5│
NB. │0.5      _0.5      _0.5  0.5│0                   0                  0│0.5   0.67082  0.5 _0.223607│                                        │0.5   0.67082  0.5│
NB. │                            │                                        │                            │                                        │                  │
NB. │  1         0         0    0│                                        │                            │                                        │                  │
NB. │  0 _0.894427 _0.447214    0│                                        │                            │                                        │                  │
NB. │  0 _0.447214  0.894427    0│                                        │                            │                                        │                  │
NB. │  0         0         0    1│                                        │                            │                                        │                  │
NB. │                            │                                        │                            │                                        │                  │
NB. │  1         0         0    0│                                        │                            │                                        │                  │
NB. │  0         1         0    0│                                        │                            │                                        │                  │
NB. │  0         0         0    1│                                        │                            │                                        │                  │
NB. │  0         0         1    0│                                        │                            │                                        │                  │
NB. └────────────────────────────┴────────────────────────────────────────┴────────────────────────────┴────────────────────────────────────────┴──────────────────┘
NB.    'Hs R Q R1 Q1'=: 10 qr A
NB.    mult/ |. Hs
NB.      0.5       0.5      0.5       0.5
NB. _0.67082 _0.223607 0.223607   0.67082
NB.      0.5      _0.5     _0.5       0.5
NB. 0.223607  _0.67082  0.67082 _0.223607
NB.    |: mult/ |. Hs
NB. 0.5  _0.67082  0.5  0.223607
NB. 0.5 _0.223607 _0.5  _0.67082
NB. 0.5  0.223607 _0.5   0.67082
NB. 0.5   0.67082  0.5 _0.223607
NB.    (|: mult/ |. Hs) mult R
NB. 1 1  1
NB. 1 2  4
NB. 1 3  9
NB. 1 4 16
NB.    Q mult R
NB. 1 1  1
NB. 1 2  4
NB. 1 3  9
NB. 1 4 16
NB.    Q1 mult R1
NB. 1 1  1
NB. 1 2  4
NB. 1 3  9
NB. 1 4 16


NB. Decompose each index in a given array into respective coordinates
toIxs=: 3 : '(#:i.)@$y'
NB. Examples:
NB.    v=: 1 2 3 4 5
NB.    toIxs v
NB. 0
NB. 1
NB. 2
NB. 3
NB. 4
NB.    ]m=: 3 4 $ i.12
NB. 0 1  2  3
NB. 4 5  6  7
NB. 8 9 10 11
NB.    toIxs m
NB. 0 0
NB. 0 1
NB. 0 2
NB. 0 3
NB.
NB. 1 0
NB. 1 1
NB. 1 2
NB. 1 3
NB.
NB. 2 0
NB. 2 1
NB. 2 2
NB. 2 3
NB.    ]t=: 3 2 2 $ i.12
NB.  0  1
NB.  2  3
NB.
NB.  4  5
NB.  6  7
NB.
NB.  8  9
NB. 10 11
NB.    toIxs t
NB. 0 0 0
NB. 0 0 1
NB.
NB. 0 1 0
NB. 0 1 1
NB.
NB.
NB. 1 0 0
NB. 1 0 1
NB.
NB. 1 1 0
NB. 1 1 1
NB.
NB.
NB. 2 0 0
NB. 2 0 1
NB.
NB. 2 1 0
NB. 2 1 1


NB. interchange of columns in matrix y, x=:col1 col2
interchangeC=: 4 : 0
'r c' =. ,"0 $ y
'c1 c2'=. x
assert. ( (c1 >: 0) *. (c1 < c) *. (c2 >: 0) *. (c2 < c))
col1=:(<(<a:),(<c1))
col2=:(<(<a:),(<c2))
((col1, col2) { y) (col2,col1) } y
)
NB. Examples:
NB.    ]y=: i. 5 6
NB.  0  1  2  3  4  5
NB.  6  7  8  9 10 11
NB. 12 13 14 15 16 17
NB. 18 19 20 21 22 23
NB. 24 25 26 27 28 29
NB.    0 1 interchangeC y
NB.  1  0  2  3  4  5
NB.  7  6  8  9 10 11
NB. 13 12 14 15 16 17
NB. 19 18 20 21 22 23
NB. 25 24 26 27 28 29
NB.    1 4 interchangeC y
NB.  0  4  2  3  1  5
NB.  6 10  8  9  7 11
NB. 12 16 14 15 13 17
NB. 18 22 20 21 19 23
NB. 24 28 26 27 25 29
NB.    1 7 interchangeC y
NB. |assertion failure: interchangeC
NB. |       ((c1>:0)*.(c1<r)*.(c2>:0)*.(c2<r))


NB. interchange of rows in matrix y, x=:row1 row2
interchangeR=: 4 : 0
'r c' =. ,"0 $ y
'r1 r2'=. x
assert. ( (r1 >: 0) *. (r1 < r) *. (r2 >: 0) *. (r2 < r))
row1=:(<(<r1),(<a:))
row2=:(<(<r2),(<a:))
((row1, row2) { y) (row2,row1) } y
)
NB. Examples:
NB.    ]y=: i. 3 5
NB.  0  1  2  3  4
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB.    1 0 interchangeR y
NB.  5  6  7  8  9
NB.  0  1  2  3  4
NB. 10 11 12 13 14
NB.    1 3 interchangeR y
NB. |assertion failure: interchangeR
NB. |       ((r1>:0)*.(r1<r)*.(r2>:0)*.(r2<r))


NB. scaling of a column by a factor in a matrix y, x=:col f
scaleC=: 4 : 0
'r c' =. ,"0 $ y
'f col'=. x
assert. ( (col >: 0) *. (col < c) )
col0=:(<(<a:),(<col))
(f * (col0 { y)) col0 } y
)
NB. Examples:
NB.    y
NB.  0  1  2  3  4
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB.    10 0 scaleC y
NB.   0  1  2  3  4
NB.  50  6  7  8  9
NB. 100 11 12 13 14
NB.    10 10 scaleC y
NB. |assertion failure: scaleC
NB. |       ((col>:0)*.(col<c))


NB. scaling of a row by a factor in a matrix y, x=:row f
scaleR=: 4 : 0
'r c' =. ,"0 $ y
'f row'=. x
assert. ( (row >: 0) *. (row < r) )
row0=:(<(<row),(<a:))
(f * (row0 { y)) row0 } y
)
NB. Examples:
NB.    y
NB.  0  1  2  3  4
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB.    10 0 scaleR y
NB.  0 10 20 30 40
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB.    10 3 scaleR y
NB. |assertion failure: scaleR
NB. |       ((row>:0)*.(row<r))


NB. addition of a scaled column col1 to an another column col2 in a matrix y, x=: f col2 col1
additionC=: 4 : 0
'r c' =. ,"0 $ y
'f c1 c2'=. x
assert. ( (c1 >: 0) *. (c1 < c) *. (c2 >: 0) *. (c2 < c))
col1=:(<(<a:),(<c1))
col2=:(<(<a:),(<c2))
((col1 { y) + (f *(col2 { y))) col1 } y
)
NB. Examples:
NB.    y
NB.  0  1  2  3  4
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB.    1 0 1 additionC y
NB.  1  1  2  3  4
NB. 11  6  7  8  9
NB. 21 11 12 13 14
NB.    1 1 2 additionC y
NB.  0  3  2  3  4
NB.  5 13  7  8  9
NB. 10 23 12 13 14
NB.    10 1 2 additionC y
NB.  0  21  2  3  4
NB.  5  76  7  8  9
NB. 10 131 12 13 14
NB.    10 1 9 additionC y
NB. |assertion failure: additionC
NB. |       ((c1>:0)*.(c1<c)*.(c2>:0)*.(c2<c))


NB. addition of a scaled row row1 to an another row row2 in a matrix y, x=: f row2 row1
additionR=: 4 : 0
'r c' =. ,"0 $ y
'f r1 r2'=. x
assert. ( (r1 >: 0) *. (r1 < r) *. (r2 >: 0) *. (r2 < r))
row1=:(<(<r1),(<a:))
row2=:(<(<r2),(<a:))
((row1 { y) + (f *(row2 { y))) row1 } y
)
NB. Examples:
NB.    y
NB.  0  1  2  3  4
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB.    1 0 1 additionR y
NB.  5  7  9 11 13
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB.    1 1 0 additionR y
NB.  0  1  2  3  4
NB.  5  7  9 11 13
NB. 10 11 12 13 14
NB.    1 1 4 additionR y
NB. |assertion failure: additionR
NB. |       ((r1>:0)*.(r1<r)*.(r2>:0)*.(r2<r))


NB. Unique elements of vector
nub=: 3 : '{./.~ y'
NB. Examples:
NB.    nub 1 2 3 1 4 5 5
NB. 1 2 3 4 5


NB. Sort rows of matrix `y` by ascending column `x`
sortRowsByColumnAsc=: ]/:{"1
NB. Examples:
NB.    a=: 2 1 3 ,. 15 6 7
NB.    a
NB. 2 15
NB. 1  6
NB. 3  7
NB.    0 sortRowsByColumnAsc a
NB. 1  6
NB. 2 15
NB. 3  7
NB.    1  a
NB. 1  6
NB. 3  7
NB. 2 15


NB. Sort rows of matrix `y` by descending column `x`
sortRowsByColumnDesc=: ]\:{"1
NB. Examples:
NB.    a=: 2 1 3 ,. 15 6 7
NB.    a
NB. 2 15
NB. 1  6
NB. 3  7
NB.    0 sortRowsByColumnDesc a
NB. 3  7
NB. 2 15
NB. 1  6
NB.    1 sortRowsByColumnDesc a
NB. 2 15
NB. 3  7
NB. 1  6


NB. Get all principal submatrices of a square matrix, consecutive planes contain row-wise principal submatrices
principalSubmatrices =: 3 : 0
'r c' =. ,"0 $ y
assert. (r = c)
1 |:\."2^:2 ] y
)
NB. Examples:
NB.    ]m=: i. 3 3
NB. 0 1 2
NB. 3 4 5
NB. 6 7 8
NB.    principalSubmatrices (i. 3 3)
NB. 4 5
NB. 7 8
NB.
NB. 3 5
NB. 6 8
NB.
NB. 3 4
NB. 6 7
NB.
NB.
NB. 1 2
NB. 7 8
NB.
NB. 0 2
NB. 6 8
NB.
NB. 0 1
NB. 6 7
NB.
NB.
NB. 1 2
NB. 4 5
NB.
NB. 0 2
NB. 3 5
NB.
NB. 0 1
NB. 3 4


NB. Get a principal submatrix specified by x of a square matrix y
principalSubmatrix=: 4 : 0
i=. 0 { x
j=. 1 { x
'r c' =. ,"0 $ y
assert. (r = c)
assert. ( (i >: 0) *. (i < r) )
assert. ( (j >: 0) *. (j < c) )
(<(<<i),(<<j)) { y
)
NB. Examples:
NB.    ]m=: i. 3 3
NB. 0 1 2
NB. 3 4 5
NB. 6 7 8
NB.    0 1 principalSubmatrix m
NB. 3 5
NB. 6 8
NB.    3 3 principalSubmatrix m
NB. |assertion failure: principalSubmatrix
NB. |       ((i>:0)*.(i<r))
NB.    2 2 principalSubmatrix m
NB. 0 1
NB. 3 4
NB.    ]m=: i. 2 3
NB. 0 1 2
NB. 3 4 5
NB.    0 1 principalSubmatrix m
NB. |assertion failure: principalSubmatrix
NB. |       (r=c)


NB. Determinant of a square matrix y
det=: 3 : 0
'r c' =. ,"0 $ y
assert. (r = c)
-/ .* y
)
NB. Examples:
NB.    ]m=: i. 3 3
NB. 0 1 2
NB. 3 4 5
NB. 6 7 8
NB.    det m
NB. 0
NB.    ]m=: i. 2 3
NB. 0 1 2
NB. 3 4 5
NB.    det m
NB. |assertion failure: det
NB. |       (r=c)


NB. Adjoint of matrix y
adjoint=: 3 : 0
'r c' =. ,"0 $ y
assert. (r = c)
det=.-/ .*
([: |: */~@($&1 _1)@# * det@principalSubmatrices) y
)
NB. Examples:
NB.    ]m=: i. 3 3
NB. 0 1 2
NB. 3 4 5
NB. 6 7 8
NB.    adjoint m
NB. _3   6 _3
NB.  6 _12  6
NB. _3   6 _3


NB. Trace of a square matrix
trace=: 3 : 0
'r c' =. ,"0 $ y
assert. (r = c)
idM=. =/~ (i.#y)
length=.*/$y
+/ ( (idM #&,i.length) { ,y)
)
NB. Examples:
NB.    ]m=: i. 3 3
NB. 0 1 2
NB. 3 4 5
NB. 6 7 8
NB.    trace m
NB. 12
NB.    ]m=: i. 4 4
NB.  0  1  2  3
NB.  4  5  6  7
NB.  8  9 10 11
NB. 12 13 14 15
NB.    trace m
NB. 30

NB. Check if a square matrix is lower triangular L
isL=: 3 : 0
'r c' =. ,"0 $ y
assert. (r = c)
uppertriang=. </~ (i.r)
elems=. uppertriang #&,y
*/ (0 = elems)
)
NB. Examples
NB.    ]A=: i. 3 3
NB. 0 1 2
NB. 3 4 5
NB. 6 7 8
NB.    isL A
NB. 0
NB.    ]lowertriang=. (>:)/~ (i.4)
NB. 1 0 0 0
NB. 1 1 0 0
NB. 1 1 1 0
NB. 1 1 1 1
NB.    isL lowertriang
NB. 1
NB.    ]nonsquare=: i.3 2
NB. 0 1
NB. 2 3
NB. 4 5
NB.    isL nonsquare
NB. |assertion failure: isL
NB. |       (r=c)

NB. Check if a square matrix is upper triangular U
isU=: 3 : 0
'r c' =. ,"0 $ y
assert. (r = c)
lowertriang=. >/~ (i.r)
elems=. lowertriang #&,y
*/ (0 = elems)
)
NB. Examples
NB.    ]A=: i. 3 3
NB. 0 1 2
NB. 3 4 5
NB. 6 7 8
NB.    isU A
NB. 0
NB.    ]uppertriang=. (<:)/~ (i.4)
NB. 1 1 1 1
NB. 0 1 1 1
NB. 0 0 1 1
NB. 0 0 0 1
NB.    isU uppertriang
NB. 1
NB.    ]nonsquare=: i.3 2
NB. 0 1
NB. 2 3
NB. 4 5
NB.    isU nonsquare
NB. |assertion failure: isU
NB. |       (r=c)

NB. Check if a square matrix is triangular, ie. either L or U
isT=: 3 : '(isL y) +. (isU y)'

NB. Represent number y > 0 as 1 and y - 1 0s
expressNum=: 3 : 0
assert (y > 0)
assert (y = (>. y))
if. y = 1 do. 1 else. (1, (y - 1) $ 0) end.
)
NB.    expressNum 0
NB. |assertion failure: assert
NB. |       assert(y>0)
NB.    expressNum 1
NB. 1
NB.    expressNum 10
NB. 1 0 0 0 0 0 0 0 0 0
NB.    expressNum 10.1
NB. |assertion failure: assert
NB. |       assert(y=(>.y))
NB.

NB. Partition matrix y by specifying consecutive partition block lenghts rows and columns
partitionMatrix=: 4 : 0
rows=.>0{x
cols=.>1{x
rowDim=.0{$y
colDim=.1{$y
assert ( (+/ rows) = rowDim)
assert ( (+/ cols) = colDim)
((;<@:expressNum"0 rows); (;<@:expressNum"0 cols)) <;.1 y
)
NB.    ]m=: i. 10 10
NB.  0  1  2  3  4  5  6  7  8  9
NB. 10 11 12 13 14 15 16 17 18 19
NB. 20 21 22 23 24 25 26 27 28 29
NB. 30 31 32 33 34 35 36 37 38 39
NB. 40 41 42 43 44 45 46 47 48 49
NB. 50 51 52 53 54 55 56 57 58 59
NB. 60 61 62 63 64 65 66 67 68 69
NB. 70 71 72 73 74 75 76 77 78 79
NB. 80 81 82 83 84 85 86 87 88 89
NB. 90 91 92 93 94 95 96 97 98 99
NB.    (2 4 1 3; 3 3 4) partitionMatrix m
NB. ┌────────┬────────┬───────────┐
NB. │ 0  1  2│ 3  4  5│ 6  7  8  9│
NB. │10 11 12│13 14 15│16 17 18 19│
NB. ├────────┼────────┼───────────┤
NB. │20 21 22│23 24 25│26 27 28 29│
NB. │30 31 32│33 34 35│36 37 38 39│
NB. │40 41 42│43 44 45│46 47 48 49│
NB. │50 51 52│53 54 55│56 57 58 59│
NB. ├────────┼────────┼───────────┤
NB. │60 61 62│63 64 65│66 67 68 69│
NB. ├────────┼────────┼───────────┤
NB. │70 71 72│73 74 75│76 77 78 79│
NB. │80 81 82│83 84 85│86 87 88 89│
NB. │90 91 92│93 94 95│96 97 98 99│
NB. └────────┴────────┴───────────┘
NB.    (2 4 1 3; 3 3 4 1) partitionMatrix m
NB. |assertion failure: assert
NB. |       assert((+/cols)=colDim)
NB.    ]m=: i. 10 5
NB.  0  1  2  3  4
NB.  5  6  7  8  9
NB. 10 11 12 13 14
NB. 15 16 17 18 19
NB. 20 21 22 23 24
NB. 25 26 27 28 29
NB. 30 31 32 33 34
NB. 35 36 37 38 39
NB. 40 41 42 43 44
NB. 45 46 47 48 49
NB.    (2 4 1 3; 3 2) partitionMatrix m
NB. ┌────────┬─────┐
NB. │0 1 2   │3 4  │
NB. │5 6 7   │8 9  │
NB. ├────────┼─────┤
NB. │10 11 12│13 14│
NB. │15 16 17│18 19│
NB. │20 21 22│23 24│
NB. │25 26 27│28 29│
NB. ├────────┼─────┤
NB. │30 31 32│33 34│
NB. ├────────┼─────┤
NB. │35 36 37│38 39│
NB. │40 41 42│43 44│
NB. │45 46 47│48 49│
NB. └────────┴─────┘


NB. Prepare frequency report of vector 'y' ordered in ascending order by elements
discreteHist=: 3 : 0
  freq=. #/.~ y
  freqRel=.(] % +/) freq
  elems=. nub y
  t=. 0 sortRowsByColumnAsc (elems ,. freqRel)
  2 2 $ 'elem' ; 'freq' ; (,.((<(<a:),(<0)) { t)) ; (,.((<(<a:),(<1)) { t))
)
NB. Examples:
NB.    discreteHist 1 1 2 3 4
NB. ┌────┬────┐
NB. │elem│freq│
NB. ├────┼────┤
NB. │1   │0.4 │
NB. │2   │0.2 │
NB. │3   │0.2 │
NB. │4   │0.2 │
NB. └────┴────┘


NB. Prepare count report of element occurance of vector 'y' in ordered bins that right boundaries are specified by consequitive elements of x
intervalHist=: 4 : 0
  hist=. <: @ (#/.~) @ (i.@#@[ , I.)
  count=. x hist y
  freq=. (% +/) count
  2 3 $ 'interval' ; 'count' ; 'freq' ; (,.x) ; (,.count) ; (,.freq)
)
NB. Examples:
NB.    ]bins=: 0.2*i:15
NB. _3 _2.8 _2.6 _2.4 _2.2 _2 _1.8 _1.6 _1.4 _1.2 _1 _0.8 _0.6 _0.4 _0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3
NB.    bins intervalHist (0 1 rnorm 100)
NB. ┌────────┬─────┬────┐
NB. │interval│count│freq│
NB. ├────────┼─────┼────┤
NB. │  _3    │ 0   │   0│
NB. │_2.8    │ 0   │   0│
NB. │_2.6    │ 0   │   0│
NB. │_2.4    │ 0   │   0│
NB. │_2.2    │ 1   │0.01│
NB. │  _2    │ 2   │0.02│
NB. │_1.8    │ 0   │   0│
NB. │_1.6    │ 1   │0.01│
NB. │_1.4    │ 2   │0.02│
NB. │_1.2    │ 3   │0.03│
NB. │  _1    │ 6   │0.06│
NB. │_0.8    │ 5   │0.05│
NB. │_0.6    │ 5   │0.05│
NB. │_0.4    │ 9   │0.09│
NB. │_0.2    │13   │0.13│
NB. │   0    │12   │0.12│
NB. │ 0.2    │10   │ 0.1│
NB. │ 0.4    │ 5   │0.05│
NB. │ 0.6    │ 6   │0.06│
NB. │ 0.8    │ 6   │0.06│
NB. │   1    │ 5   │0.05│
NB. │ 1.2    │ 5   │0.05│
NB. │ 1.4    │ 2   │0.02│
NB. │ 1.6    │ 0   │   0│
NB. │ 1.8    │ 1   │0.01│
NB. │   2    │ 0   │   0│
NB. │ 2.2    │ 1   │0.01│
NB. │ 2.4    │ 0   │   0│
NB. │ 2.6    │ 0   │   0│
NB. │ 2.8    │ 0   │   0│
NB. │   3    │ 0   │   0│
NB. └────────┴─────┴────┘


NB. mean of a vector
mean=: +/ % #
NB. Examples:
NB.    mean 1 1 1 1
NB. 1
NB.    mean 1 1 1 2 3
NB. 1.6


NB. variance of a sample
var=: (+/@(*:@(] - +/ % #)) % #)"1
NB. Examples:
NB.    var 1 1 1 1
NB. 0
NB.    var 1 2 1 2
NB. 0.25

rnd=. ?@#&0
rne=. -@* ^.@rnd
rno=: 3 : '(%:2 rne y)*2 o.(rnd y)*o.2'


NB. Generation of y samples of normal distribution having 'mean var'=:x
rnorm=: 4 : 0
  (rno y) #.every <|.x
)
NB. Examples:
NB. 10 samples of N(0,1)
NB.      0 1 rnorm 10
NB. 0.583949 0.151037 1.44553 1.18409 _0.53704 1.49066 0.649399 0.569303 0.855299 0.738213


NB. Generation of y samples of uniform distribution having bounds 'from to'=:x
runiform=: 4 : 0
  rnd=. ?@#&0
  relist=. (-/ , {:)@|.
  > (rnd y) #.each <relist x
)
NB. Examples:
NB. 10 samples of U(0,1)
NB.      0 1 runiform 10
NB. 0.183411 0.0968962 0.587723 0.165308 0.68218 0.0916652 0.00554653 0.149567 0.340257 0.370271

NB. --------------------
NB. - Property testing -
NB. --------------------

NB. Check equality for the property requiring matrices and scalars
NB. x is gerund relation having leftR`rightR and both
NB. S leftR (<M)
NB. S rightR (<M)
NB. where S=:S1,S2,..,SI and M=:M1;M2;..;MI
NB. and y is data=:S;<M
NB. Please notice that M is also boxed as it can assume a sequence of
NB. matrices with different shapes
checkEqOfMatricesScalarsRel=: {{
  ( (>0{y) x@.0 (>1{y) ) -: ( (>0{y) x@.1 (>1{y) )
}}
NB. Examples:
NB.    scalars=: _100 100 runiform 1
NB.    matrices=: (genUniformMatrix 4 2);(genUniformMatrix 2 6)
NB.    ]data=:scalars;<matrices
NB. ┌───────┬──────────────────────────────────────────────────────────────────────┐
NB. │24.4943│┌─────────────────┬──────────────────────────────────────────────────┐│
NB. │       ││  892.71  850.185│_619.827 464.316 _601.967 309.364 851.114 _181.237││
NB. │       ││  636.44  248.583│ 238.782 783.326  777.188 433.257 992.401 _44.5578││
NB. │       ││_714.152 _577.337│                                                  ││
NB. │       ││_556.986 _412.428│                                                  ││
NB. │       │└─────────────────┴──────────────────────────────────────────────────┘│
NB. └───────┴──────────────────────────────────────────────────────────────────────┘
NB.
NB.    rightR=: {{
NB. s=.0{x
NB. A=.>(0{y)
NB. B=.>(1{y)
NB. A mult (s*B)
NB. }}
NB.       leftR=: {{
NB. s=:0{x
NB. A=:>(0{y)
NB. B=:>(1{y)
NB. s*(A mult B)
NB. }}
NB.    relation=: leftR`rightR
NB.    relation checkEqOfMatricesScalarsRel data
NB. 1
