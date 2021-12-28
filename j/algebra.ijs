NB. --------------------------------
NB. -- Algebra - useful functions --
NB. --------------------------------

NB. How to load the script - an example
NB. ] scriptdir =: 'PATH-TO-REPO/j-random-matrices/j/'
NB. 0!:1 < scriptdir,'algebra.ijs'

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
NB.    nub 1 2 3 1 4 5 5
NB. 1 2 3 4 5

NB. Sort rows of matrix `y` by ascending column `x`
sortRowsByColumnAsc=: ]/:{"1
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
  rcount=. (% +/) count
  2 3 $ 'interval' ; 'count' ; 'relative count' ; (,.x) ; (,.count) ; (,.rcount)
)
NB.    ]bins=: 0.2*i:15
NB. _3 _2.8 _2.6 _2.4 _2.2 _2 _1.8 _1.6 _1.4 _1.2 _1 _0.8 _0.6 _0.4 _0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3
NB.    bins intervalHist (0 1 rnorm 100)
NB. ┌────────┬─────┬──────────────┐
NB. │interval│count│relative count│
NB. ├────────┼─────┼──────────────┤
NB. │  _3    │ 0   │        0     │
NB. │_2.8    │ 0   │        0     │
NB. │_2.6    │ 0   │        0     │
NB. │_2.4    │ 1   │ 0.010101     │
NB. │_2.2    │ 2   │ 0.020202     │
NB. │  _2    │ 0   │        0     │
NB. │_1.8    │ 0   │        0     │
NB. │_1.6    │ 1   │ 0.010101     │
NB. │_1.4    │ 1   │ 0.010101     │
NB. │_1.2    │ 5   │0.0505051     │
NB. │  _1    │ 4   │ 0.040404     │
NB. │_0.8    │ 6   │0.0606061     │
NB. │_0.6    │ 3   │ 0.030303     │
NB. │_0.4    │ 5   │0.0505051     │
NB. │_0.2    │ 5   │0.0505051     │
NB. │   0    │ 4   │ 0.040404     │
NB. │ 0.2    │ 7   │0.0707071     │
NB. │ 0.4    │12   │ 0.121212     │
NB. │ 0.6    │13   │ 0.131313     │
NB. │ 0.8    │ 7   │0.0707071     │
NB. │   1    │10   │  0.10101     │
NB. │ 1.2    │ 4   │ 0.040404     │
NB. │ 1.4    │ 6   │0.0606061     │
NB. │ 1.6    │ 1   │ 0.010101     │
NB. │ 1.8    │ 1   │ 0.010101     │
NB. │   2    │ 0   │        0     │
NB. │ 2.2    │ 0   │        0     │
NB. │ 2.4    │ 0   │        0     │
NB. │ 2.6    │ 0   │        0     │
NB. │ 2.8    │ 0   │        0     │
NB. │   3    │ 1   │ 0.010101     │
NB. │        │ 0   │        0     │
NB. └────────┴─────┴──────────────┘

NB. mean of a vector
mean=: +/ % #
NB.    mean 1 1 1 1
NB. 1
NB.    mean 1 1 1 2 3
NB. 1.6

NB. variance of a sample
var=: (+/@(*:@(] - +/ % #)) % #)"1
NB.    var 1 1 1 1
NB. 0
NB.    var 1 2 1 2
NB. 0.25

NB. Generation of y samples of normal distribution having 'mean var'=:x
rnorm=: 4 : 0
  rnd=. ?@#&0
  rne=. -@* ^.@rnd
  rno=. 3 : '(%:2 rne y)*2 o.(rnd y)*o.2'
  (rno y) #.every <|.x
)