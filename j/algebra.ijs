NB. --------------------------------
NB. -- Algebra - useful functions --
NB. --------------------------------

NB. How to load the script - an example
NB. ]scriptdir=: 'PATH-TO-REPO/j-random-matrices/j/'
NB. 0!:1 < scriptdir,'algebra.ijs'

NB. Transpose of a matrix
tr=: |:
NB.    tr (2 3 $ 1 2 3 4 5 6)
NB. 1 4
NB. 2 5
NB. 3 6
NB.    tr (tr (2 3 $ 1 2 3 4 5 6))
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
NB.

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

NB. Check equality for the property requiring two matrices
NB. x is gerund having leftR`rightR and both
NB.  M1 leftR M2
NB. M1 rightR M2
NB. are expected to work
NB. y is M1,:M2
checkEqTwoMatrices=: 4 : '( (0{y) x@.0 (1{y) ) -: ( (0{y) x@.1 (1{y) )'
NB. Examples:
NB.    leftR=: 4 : 'x + y'
NB.    rightR=: 4 : 'y + x'
NB.    relation=: leftR`rightR
NB.       relation@.0
NB. 4 : 'x + y'
NB.       relation@.1
NB. 4 : 'y + x'
NB.    ]matrices=: (genUniformMatrix 4 2),:(genUniformMatrix 4 2)
NB.  783.326  777.188
NB.  433.257  992.401
NB. _44.5578   892.71
NB.  850.185   636.44
NB.
NB.  211.161 _619.827
NB.  464.316 _601.967
NB.  309.364  851.114
NB. _181.237  238.782
NB. relation checkEqTwoMatrices matrices
NB. 1

NB. Check equality for the property requiring more than 1 matrices
NB. x is gerund having leftR`rightR and both
NB. leftR M
NB. rightR M
NB. are expected to work with
NB. y is M=,M1,M2,:M3 and so on
checkEqOfMatrices=: 4 : '( x@.0 y ) -: ( x@.1 y )'
NB. Examples:
NB.    leftR=: 3 : '( (0{y) + (1{y) ) + (2{y)'
NB.    rightR=: 3 : '(0{y) + ( (1{y) + (2{y) )'
NB.    relation=: leftR`rightR
NB.    genUniformMatrix=: 3 : 'y $ _1000 1000 runiform ((0{y) * (1{y))'
NB.    ]matrices=: (genUniformMatrix 4 2),(genUniformMatrix 4 2),:(genUniformMatrix 4 2)
NB.  211.161 _619.827
NB.  464.316 _601.967
NB.  309.364  851.114
NB. _181.237  238.782
NB.
NB. _428.685  853.433
NB. _400.652 _164.792
NB.  375.372  675.547
NB.  584.175 _69.8546
NB.
NB.  244.943 _350.586
NB.   815.65 _873.687
NB.  _226.76 _322.805
NB. _808.466  202.957
NB.    relation checkEqOfMatrices matrices
NB. 1

NB. Check equality for the property requiring matrices and scalars
NB. x is gerund relation having leftR`rightR and both
NB. S leftR M
NB. S rightR M
NB. are expected to work with
NB. y is data=:S;M
checkEqOfMatricesScalarsRel=: 4 : '( (>0{y) x@.0 (>1{y) ) -: ( (>0{y) x@.1 (>1{y) )'
NB. Examples:
NB.    ]scalars=: _100 100 runiform 1
NB. 24.4943
NB.    ]matrices=: (genUniformMatrix 4 2),:(genUniformMatrix 4 2)
NB.  853.433 _400.652
NB. _164.792  375.372
NB.  675.547  584.175
NB. _69.8546  211.161
NB.
NB. _350.586   815.65
NB. _873.687  _226.76
NB. _322.805 _808.466
NB.  202.957 _428.685
NB.    data=:scalars;matrices
NB.    >0{data
NB. 24.4943
NB.    >1{data
NB.  853.433 _400.652
NB. _164.792  375.372
NB.  675.547  584.175
NB. _69.8546  211.161
NB.
NB. _350.586   815.65
NB. _873.687  _226.76
NB. _322.805 _808.466
NB.  202.957 _428.685
NB.
NB.    leftR=: 4 : '(0{x) * ( (0{y) + (1{y) )'
NB.    rightR=: 4 : '( (0{x) * (0{y) ) + ( (0{x) * (1{y) )'
NB.    relation=: leftR`rightR
NB.    relation checkEqOfMatricesScalarsRel data
NB. 1