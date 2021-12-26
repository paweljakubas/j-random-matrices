# Linear algebra

## Contents
1. [Selecting from a matrix](#selecting-from-matrix)
2. [Updating a matrix](#updating-of-matrix)
3. [Generate a random matrix](#generate-random-matrix) - TODO
4. [Elementary operations of a matrix](#elementary-operations-in-matrix)
5. [Transpose of a matrix](#transpose-of-matrix) - TODO
6. [Inverse of a matrix](#inverse-of-matrix) - TODO
7. [Determinant of a matrix](#determinant-of-matrix) - TODO
8. [Trace of a matrix](#trace-of-matrix) - TODO
9. [A partitioned matrix](#partitioned-matrix) - TODO
10. [A matrix decompositions](#matrix-decompositions) - TODO


[Solutions to exercices](#linear-algebra-solutions-to-exercises)

# Linear algebra
Let's start with fundamental results and operations on matrices
that are prerequisite for going further. There are topics that I dive into quite deeply,
others I just stop investigating when accomplishing basic result. This is the reflection
of my subjective thinking what I regard will be especially important to master the later steps.

## Selecting from matrix
Generic selector can be defined by specifying both indices, row's or column's. The way the selection
works can also be extended to higher dimensional arrays, for example, in case of tensor the
plane index should also be additionally specified.
In the first example below the selector allows cutting submatrix by taking first two rows and columns,
ie., via enforcing 0 and 1 indices in the selector:
```
   ]m=: 100 + 4 4 $ i.16
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115

   sel=: (<(<0 1),(<0 1))
   sel { m
100 101
104 105
```
There is a way to choose either whole columns or rows by using *a:* :
```
   sel=: (<(<a:),(<0 1))
   sel { m
100 101
104 105
108 109
112 113

   sel=: (<(<2),(<a:))
   sel { m
108 109 110 111
```
Also one can define a selection by saying which indices to omit:
```
sel=: (<(<<2 3 1),(<a:))
   sel { m
100 101 102 103

   sel=: (<(<a:),(<<0 1))
   sel { m
102 103
106 107
110 111
114 115
```
There is also a way to specify indices starting the count from the end via prefixing negative sign
to the index. In that case *_1* is the last index:
```
   sel=: (<(<a:),(<_1))
   sel { x
103 107 111 115
```
It is also possible to combine selections with each other:
```
   sel=:((< (<0 1), (<0 1)),(<(<2 3),(<2 3)))
   sel { m
100 101
104 105

110 111
114 115
```
Let's now see how we can negate a selector and select all elements except those specified by the selector.
Above we were omitting indices in a given axis, now we want to learn how to treat a selection as mask
and take everything except what the selection pinpoints. As the result can be non-rectangular, we need to
realize the operation in a linearized form to make sure we have a general solution.
After we get the result we can reshape it as we want:
```
   sel=:(< (<1 2), (<1 2))
   sel { m
105 106
109 110

   sel { i.$m
5  6
9 10
   ,sel { i.$m
5 6 9 10

   ]selOmittedIxs =: (< (<< (,sel { i.$m)) ) { ,i.$m
0 1 2 3 4 7 8 11 12 13 14 15

   (< (<< (,sel { i.$m)) ) { ,m
100 101 102 103 104 107 108 111 112 113 114 115

   2 6 $(< (<< (,sel { i.$m)) ) { ,m
100 101 102 103 104 107
108 111 112 113 114 115
```

**Exercise 1**
How to get *m1* that has only elements of *m* with odd indices?
Try two approaches with a selection and a negated selection.
```
100 101 102 103
104 105 106 107  -->   105 107
108 109 110 111        113 115
112 113 114 115
```
[Solution to exercise 1](#solution-to-exercise-1)

**Exercise 2**
Define the tensor that has three planes, the first plane is *m* and the next ones has
corresponding elements increased by 100 and 200 with respect to the first plane.
  (a) cut the tensor using a selection in such a way that only edges containing 100, 112, 103 and 115 (plus 100 and 200) are maintained
  (b) cut the tensor using a selection that only inner (non-surface) elements are maintained

[Solution to exercise 2](#solution-to-exercise-2)

We can also select from a given array by specifying a predicate that filters the array's values.
Let's say we want to select only those values that divide without remainder by 3
```
   3 | (,m)
1 2 0 1 2 0 1 2 0 1 2 0 1 2 0 1
   (0 = 3 | (,m))
0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0
   pred =: 3 : '(0 = 3 | y) # (i. # y)'
   pred (,m)
2 5 8 11 14
   (pred (,m)) { (,m)
102 105 108 111 114
```

**Exercise 3**
Select all values from *m* that does not divide without remainder by either 3 or 5

[Solution to exercise 3](#solution-to-exercise-3)

Finally, we will see how to introduce functions that act on indices of elements of arrays.
First let's show how to get handle on them:
```
      ]m=: 100 + 4 4 $ i.16
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115
   toIxs=: 3 : '(#:i.)@$y'
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

3 0
3 1
3 2
3 3
```
In order to get diagonal elements one can do the following:
```
   nub=: 3 : '{./.~ y'
   nub 1 2 3 4 5 6
   nub 1 2 3 4 5 6
1 2 3 4 5 6
   nub 1 2 3 4 5 6 1 2 3 4
1 2 3 4 5 6
   < "1 (toIxs m)
┌───┬───┬───┬───┐
│0 0│0 1│0 2│0 3│
├───┼───┼───┼───┤
│1 0│1 1│1 2│1 3│
├───┼───┼───┼───┤
│2 0│2 1│2 2│2 3│
├───┼───┼───┼───┤
│3 0│3 1│3 2│3 3│
└───┴───┴───┴───┘
   diag=: 3 : '1 = # nub y'
   diag "1 (toIxs m)
1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1

   (diag "1 (toIxs m)) #&,m
100 105 110 115
```
One can also utilize the following scheme to generate diagonal or triangular selections.
*(x u/ y)* returns a table having entries *(a u b)* for every *a* in *x* and *b* in *y*.
```
   ]diag=: =/~ (i.#m)
1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1
   diag #&,m
100 105 110 115

   ]uppertriang=: (<:)/~ (i.#m)
1 1 1 1
0 1 1 1
0 0 1 1
0 0 0 1
   uppertriang #&,m
100 101 102 103 105 106 107 110 111 115
```

**Exercise 4**
Select two diagonals (diagonal and cross diagonal) of *m* using index function.
Wilkinson diagrams for diagonal and cross diagonal matrices:
```
diagonal             cross diagonal
X 0 0 0 0              0 0 0 0 X
0 X 0 0 0              0 0 0 X 0
0 0 X 0 0              0 0 X 0 0
0 0 0 X 0              0 X 0 0 0
0 0 0 0 X              X 0 0 0 0
```


[Solution to exercise 4](#solution-to-exercise-4)

**Exercise 5**
Select tridiagonal values of *m* that are odd.  Wilkinson diagram of band matrix is
below (a tridiagonal matrix is a special case of a band matrix):
```
X X 0 0 0
X X X 0 0
0 X X X 0
0 0 X X X
0 0 0 X X
```
Do the same for lower bidiagonal matrix (Wilkinson diagram below)
```
X 0 0 0 0
X X 0 0 0
0 X X 0 0
0 0 X X 0
0 0 0 X X
```

[Solution to exercise 5](#solution-to-exercise-5)

**Exercise 6**
Use the tensor from Exercise 2
  (a) cut the diagonal plane from the tensor
  (b) cut the any two ortogonal planes inside the tensor and produce a new two plane tensor out of them

[Solution to exercise 6](#solution-to-exercise-6)

**Summary**: in order to select from a matrix we have quite extensive arsenal at disposal:
(a) we can use J selections,
(b) negated J selections,
(c) deliver functions filtering values,
(d) write index functions,
and finally (e) use combination of those approaches

## Updating of matrix
We can update a matrix with new values and a selection:
```
   sel=: (<(<0 1),(<0 1))
   sel { m
100 101
104 105
   1 sel } m
  1   1 102 103
  1   1 106 107
108 109 110 111
112 113 114 115
```
We need to be sure that the shape of the delivered new values is compatible with what selection determines:
```
   (2 2 $ 1 2 3 4) sel } m
  1   2 102 103
  3   4 106 107
108 109 110 111
112 113 114 115

   sel=: ((<(<0 1),(<0 1)),(<(<2 3),(<2 3)))
   sel { m
100 101
104 105

110 111
114 115

   (1 2 2 $ 1 2 3 4), (1 2 2 $ 5 6 7 9)
1 2
3 4

5 6
7 9
   ((1 2 2 $ 1 2 3 4), (1 2 2 $ 5 6 7 9)) sel } m
  1   2 102 103
  3   4 106 107
108 109   5   6
112 113   7   9

   (2 2 2 $ 1 2 3 4 5 6 7 9) sel } m
  1   2 102 103
  3   4 106 107
108 109   5   6
112 113   7   9
```
The question of ovelapping selection arises. The overwriting is deterministic in a sense that
the selection on the right overwrites the values of the left selection:
```
   sel=: ((<(<0 1),(<0 1)),(<(<1 2),(<1 2)))
   sel { m
100 101
104 105

105 106
109 110
   ((1 2 2 $ 0 0 0 0), (1 2 2 $ 1 1 1 1)) sel } m
  0   0 102 103
  0   1   1 107
108   1   1 111
112 113 114 115
   ((1 2 2 $ 1 1 1 1), (1 2 2 $ 0 0 0 0)) sel } m
  1   1 102 103
  1   0   0 107
108   0   0 111
112 113 114 115
```
There is also a way to update with predicate function:
```
pred =: 4 : '(y > x) # (i. # y)'
   (110 pred (,m)) { (,m)
111 112 113 114 115
   ($m) $ 100 (110 pred (,m)) } ,m
100 101 102 103
104 105 106 107
108 109 110 100
100 100 100 100
```
For matrices, as seen above, one can go through linearized form of matrix (*,*) and
after applying the operation retrieve the shape using (*$*)
We can also use two predicates and apply different updates for each predicate in one go:
```
   pred1 =: 4 : '((y > (0 { x)) *. (y < (1 { x))) # (i. # y)'
   (101 104 pred1 (,m)) { (,m)
102 103

   100 (110 pred (,m)) } ,m
100 101 102 103 104 105 106 107 108 109 110 100 100 100 100 100
   200 (101 104 pred1 (,m)) } ,m
100 101 200 200 104 105 106 107 108 109 110 111 112 113 114 115

   ((5$100),(2$200)) ((110 pred (,m)),(101 104 pred1 (,m))) } ,m
100 101 200 200 104 105 106 107 108 109 110 100 100 100 100 100

   ($m) $ ((5$100),(2$200)) ((110 pred (,m)),(101 104 pred1 (,m))) } ,m
100 101 200 200
104 105 106 107
108 109 110 100
100 100 100 100
```
One can notice that the updating using selection is rectangular in shape and one can
operate in matrix shape doing this. If the updating is done using predicate function
there is no guarantee that the updated slice of a matrix is rectangular, and one needs to
use linearized indices. The same is, in general, true for the negated selection,
ie. updating everything in an array except what a given selection determines. Example below:
```
   sel=:(< (<1 2), (<1 2))
   sel { m
105 106
109 110
   sel { i.$m
5  6
9 10
   ]toUpdateIxs=:(< (<< (,sel { i.$m)) ) { ,i.$m
0 1 2 3 4 7 8 11 12 13 14 15
   ($m) $ 0 toUpdateIxs } ,m
0   0   0 0
0 105 106 0
0 109 110 0
0   0   0 0
```
Now rather than to update with new values let's update the selected elements using
the specified function that uses the old values.
```
   f =: 3 : '(y + 100)'
   (f (sel { m))
205 206
209 210
   (f (sel { m)) sel } m
100 101 102 103
104 205 206 107
108 209 210 111
112 113 114 115
```

**Exercise 7**
Work with the m matrix and using a negated selection update all elements except what below two selections point to:
one selecting 0-1 row and 0-1 column, the second one selecting 3,3 element.
  (a) Set all to be updated elements to 0
  (b) Increment all to be updated elements by 1
  (c) Set all to be updated elements to value that is average of selected elements
  (d) Set all to be updated elements to value that is average of unselected elements

[Solution to exercise 7](#solution-to-exercise-7)

**Exercise 8**
Turn *m* to (a) upper triangular matrix, (b) upper Hessenberg matrix, and (c) upper cross triangular matrix.
Corresponding Wilkinson diagrams:
```
upper Hessenberg     upper triangular    upper cross triangular
X X X X X            X X X X X           X X X X X
X X X X X            0 X X X X           X X X X 0
0 X X X X            0 0 X X X           X X X 0 0
0 0 X X X            0 0 0 X X           X X 0 0 0
0 0 0 X X            0 0 0 0 X           X 0 0 0 0
```

[Solution to exercise 8](#solution-to-exercise-8)

Let's finally investigate how to update some neighborhood of a given point in a matrix. In an example below we want to
increment the nearest neighbors of the point. If the point is going to be defined by *(x,y)* then we will
update 4 points: *(x-1,y)*,*(x+1,y)*,*(x,y-1)*,*(x,y+1)*. The complication that arises here regards points chosen on the
boundary of a matrix. Let's assume that in that case we want to omit the points tresspassing the boundary (in other version
we could also contemplate including them wrapped up on the oppostite side).
```
   nn=: 3 : '( ((0{y)-1), (1{y) ),( ((0{y)+1), (1{y) ),( (0{y), ((1{y)-1) ),:( (0{y), ((1{y)+1) )'

   ]nn11=: nn 1 1
0 1
2 1
1 0
1 2
   ]nn00=: nn 0 0
_1  0
 1  0
 0 _1
 0  1

   NB. let's filter out neighbors that are outside boundary of the matrix
   ]maxR=. (0{$m) - 1
3
   ]maxC=. (1{$m) - 1
3
   validateNeighbors=: 3 : '((0{y) >: 0) *. ((0{y) <: maxR) *. ((1{y) >: 0) *. ((1{y) <: maxC)'
   validateNeighbors (3 1)
1
   validateNeighbors (4 1)
0
   validateNeighbors (3 _1)
0
   validateNeighbors"1 nn00
0 1 0 1
   ]nn00Filtered=: (validateNeighbors"1 nn00) # (i.(0{$nn00)) { nn00
1 0
0 1
   ]nn11Filtered=: (validateNeighbors"1 nn11) # (i.(0{$nn11)) { nn11
0 1
2 1
1 0
1 2

  NB. let's create 0-1 matrix of what to update
   calcIxs=: 3 : '(1{y) + ((0{y)*(maxC+1))'
   calcIxs 0 3
3
   calcIxs 1 3
7
   ($m) $ 1 (calcIxs"1 nn00Filtered) } (,($m) $ 0)
0 1 0 0
1 0 0 0
0 0 0 0
0 0 0 0

   ($m) $ 1 (calcIxs"1 nn11Filtered) } (,($m) $ 0)
0 1 0 0
1 0 1 0
0 1 0 0
0 0 0 0
   ]nn11Ixs=: calcIxs"1 nn11Filtered
1 9 4 6
   ]nn11Vals=: (calcIxs"1 nn11Filtered) { ,m
101 109 104 106
   ]m11=:($m) $ (>: nn11Vals) nn11Ixs } ,m
100 102 102 103
105 105 107 107
108 110 110 111
112 113 114 115

   NB.checking the correctness
   m11 - m
0 1 0 0
1 0 1 0
0 1 0 0
0 0 0 0

   NB. the same for (0,0) point
   ]nn00Ixs=: calcIxs"1 nn00Filtered
4 1
   ]nn00Vals=: (calcIxs"1 nn00Filtered) { ,m
104 101
   ]m00=:($m) $ (>: nn00Vals) nn00Ixs } ,m
100 102 102 103
105 105 106 107
108 109 110 111
112 113 114 115
   m00 - m
0 1 0 0
1 0 0 0
0 0 0 0
0 0 0 0
```

**Exercise 9**
Show capability to update (decrement) the nearest neighbors in a tensor

[Solution to exercise 9](#solution-to-exercise-9)

**Exercise 10**
(a) Show capability to update (increment) the second nearest neighbors in a matrix.
(b) Also show updating simultaneuously the nearest neighbors of two points at the same time
(both incrementing, then one incrementing second decrementing)
(c) Also show updating simultaneuously the nearest neighbors (incrementing) and
second nearest neighbors (decrementing)

[Solution to exercise 10](#solution-to-exercise-10)

**Exercise 11**
Show capability to update (make 0) both diagonals passing through the point.

[Solution to exercise 11](#solution-to-exercise-11)

**Summary**: As in the case of selecting the updating of a matrix can be realized in a number of ways:
(a) via J selections which requires rectangular updating values mimicking the shape of the selection,
(b) negated selections but then we need to work with linearized indices,
(c) functions acting on both values or indices to filtering out indices to be updated
(d) on top of indices to be updated we can provide new values independent or dependent
on the current values
(e) finally we can extend updating spacially, ie. beyond pointwise updating, and come up with neighborhood updating

## Generating random matrix

## Elementary operations in matrix
There are three elementary operations we are going to cover here, all three in the context of both rows and columns.
Let's start with **interchange** elementary operations.
For column case we define two selection for each column we want to interchange, then use the pair selection in 'view'
and the interchanged pair in the same update:
```
   col0=:(<(<a:),(<0))
   col0 { m
100 104 108 112
   col2=:(<(<a:),(<2))
   col2 { m
102 106 110 114

   (col0, col2) { m
100 104 108 112
102 106 110 114

   ((col0, col2) { m) (col2,col0) } m
102 101 100 103
106 105 104 107
110 109 108 111
114 113 112 115
```
For a row case we do analogically but we choose corresponding row selections:
```
   row0=:(<(<0),(<a:))
   row0 { m
100 101 102 103
   row2=:(<(<2),(<a:))
   row2 { m
108 109 110 111
   (row0, row2) { m
100 101 102 103
108 109 110 111
   ((row0, row2) { m) (row2,row0) } m
108 109 110 111
104 105 106 107
100 101 102 103
112 113 114 115
```

In **scaling** elementary operation each element in a give column (row) is multiplied
by the provided factor. For both column and row cases this could be realized as follows:
```
f =: 4 : '(y * x)'
   (5 f (col0 { m))
500 520 540 560
   (5 f (col0 { m)) col0 } m
500 101 102 103
520 105 106 107
540 109 110 111
560 113 114 115

   (5 f (row0 { m)) row0 } m
500 505 510 515
104 105 106 107
108 109 110 111
112 113 114 115
```

The **addition** elementary operation in a column case entails adding element-wise to some column
another column scaled by some factor. The case for the row varies with the choice of the row selector
rather than column selector. Following we are replacing column 0 with the result of the
addition of column 0 and column 2 that was scaled by factor 5. Once again we use update operation
with selectors. Then we do the same for the row 0:
```
   (col0 { m) + (5 f (col2 { m))
610 634 658 682

   ((col0 { m) + (5 f (col2 { m))) col0 } m
610 101 102 103
634 105 106 107
658 109 110 111
682 113 114 115

   (row0 { m) + (5 f (row2 { m))
640 646 652 658
   (row0 { m) + (5 f (row2 { m)) row0 } m
640 645 650 655
205 206 207 208
210 211 212 213
215 216 217 218
```

**Exercise 12**
Show that the three basic operations can be realized by multiplication of transformed identity matrices.

[Solution to exercise 12](#solution-to-exercise-12)

## Linear algebra. Solutions to exercises
### Solution to exercise 1
```
   m
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115
   sel=:(< (<1 3), (<1 3))
   sel { m
105 107
113 115

   ]omittedVals=: (,(<(<0 2),(<a:)) { m),(,(<(<a:),(<0 2)) { m)
100 101 102 103 108 109 110 111 100 102 104 106 108 110 112 114
   ]omittedIxs=: (,(<(<0 2),(<a:)) { i.$m),(,(<(<a:),(<0 2)) { i.$m)
0 1 2 3 8 9 10 11 0 2 4 6 8 10 12 14
   2 2 $ (<(<<omittedIxs))    {,m
105 107
113 115
```

### Solution to exercise 2
```
   ]t=: (m ,: m+100) , m+200
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115

200 201 202 203
204 205 206 207
208 209 210 211
212 213 214 215

300 301 302 303
304 305 306 307
308 309 310 311
312 313 314 315
   $t
3 4 4

   NB. (a) cut the tensor using a selection in such a way that only edges containing 100, 112, 103 and 115 (plus 100 and 200) are maintained
      e1=:(<(<a:),(<0),(<0))
   e1 { t
100 200 300
   e2=:(<(<a:),(<0),(<3))
   e2 { t
103 203 303
   e3=:(<(<a:),(<3),(<0))
   e3 { t
112 212 312
   e4=:(<(<a:),(<3),(<3))
   e4 { t
115 215 315
   edges =: 3 2 2 $ , (e1 { i.$t) ,. (e2 { i.$t) ,. (e3 { i.$t) ,. (e4 { i.$t)
   edges { ,t
100 103
112 115

200 203
212 215

300 303
312 315

   NB. (b) cut the tensor using a selection that only inner (non-surface) elements are maintained
   sel=:(<(<1),(<1 2),(<1 2))
   ]t1=:sel { t
205 206
209 210
   $t1
2 2
   sel=:(<(<<0 2),(<<0 3),(<<0 3))
   ]t1=:sel { t
205 206
209 210
   $t1
1 2 2
```

### Solution to exercise 3
```
   m
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115

   (0 = 3 | (,m)) +. (0 = 5 | (,m))
1 0 1 0 0 1 0 0 1 0 1 1 0 0 1 1
   -. ((0 = 3 | (,m)) +. (0 = 5 | (,m)))
0 1 0 1 1 0 1 1 0 1 0 0 1 1 0 0
   pred =: 3 : '(-. ((0 = 3 | y) +. (0 = 5 | y))) # (i. # y)'
   pred (,m)
1 3 4 6 7 9 12 13
   (pred (,m)) { (,m)
101 103 104 106 107 109 112 113
```

### Solution to exercise 4
```
   m
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115
   crossdiag=: 4 : '(#x) = (0 { y) + (1 { y) +1'
   m crossdiag "1 (toIxs m)
0 0 0 1
0 0 1 0
0 1 0 0
1 0 0 0
   diags=: (m crossdiag "1 (toIxs m)) + (diag "1 (toIxs m))
1 0 0 1
0 1 1 0
0 1 1 0
1 0 0 1
   diags #&,m
100 103 105 106 109 110 112 115
```

### Solution to exercise 5
```
   m
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115
   f =: 4 : '((x - y) <: 1) *. ((x - y) >: _1)'
   ]tridiag=: (f"0)/~ (i.#m)
1 1 0 0
1 1 1 0
0 1 1 1
0 0 1 1
   tridiag #&,m
100 101 104 105 106 109 110 111 114 115

   f =: 4 : '((x - y) <: 1) *. (x >: y)'
   ]lowerbidiag=: (f"0)/~ (i.#m)
1 0 0 0
1 1 0 0
0 1 1 0
0 0 1 1
   lowerbidiag #&,m
100 104 105 109 110 114 115
```

### Solution to exercise 6
```
   ]t=: (m ,: m+100) , m+200
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115

200 201 202 203
204 205 206 207
208 209 210 211
212 213 214 215

300 301 302 303
304 305 306 307
308 309 310 311
312 313 314 315
   $t
3 4 4

   NB.(a) cut the diagonal plane from the tensor
   ]diag=: =/~ (i.#m)
1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1
   process=:3 : 'diag#&,y'
   process"2 t
100 105 110 115
200 205 210 215
300 305 310 315

   NB. (b) cut the any two ortogonal planes inside the tensor and produce a new two plane tensor out of them
   ((<(<a:),(<1),(<a:)),.(<(<a:),(<a:),(<1))) { t
104 105 106 107
204 205 206 207
304 305 306 307

101 105 109 113
201 205 209 213
301 305 309 313
   $ ((<(<a:),(<1),(<a:)),.(<(<a:),(<a:),(<1))) { t
2 3 4
```

### Solution to exercise 7
```
   m
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115

   sel1=:(< (<0 1), (<0 1))
   sel1 { m
100 101
104 105
   sel2=:(< (<3), (<3))
   sel2 { m
115
   ((,sel1 { m), (,sel2 { m))
100 101 104 105 115
   sels=: 3 : '((,sel1 { y), (,sel2 { y))'
      ]toUpdateIxs=:(< (<< (sels i.$m)) ) { ,i.$m
2 3 6 7 8 9 10 11 12 13 14

   NB. (a) Set all to be updated elements to 0
   ($m) $ 0 toUpdateIxs } ,m
100 101 0   0
104 105 0   0
  0   0 0   0
  0   0 0 115

  NB. (b) Increment all to be updated elements by 1
   ]toUpdateVals=:((< (<< (sels i.$m)) ) { ,i.$m) { ,m
102 103 106 107 108 109 110 111 112 113 114
   ($m) $ (>: toUpdateVals) toUpdateIxs } ,m
100 101 103 104
104 105 107 108
109 110 111 112
113 114 115 115

  NB. (c) Set all to be updated elements to value that is average of selected elements
  mean=: +/ % #
   ($m) $ (mean (,sels m)) toUpdateIxs } ,m
100 101 105 105
104 105 105 105
105 105 105 105
105 105 105 115

  NB. (d) Set all to be updated elements to value that is average of unselected elements
   ($m) $ (mean toUpdateVals) toUpdateIxs } ,m
    100     101 108.636 108.636
    104     105 108.636 108.636
108.636 108.636 108.636 108.636
108.636 108.636 108.636     115
```

### Solution to exercise 8
```
   m
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115

  NB. (a) upper triangular matrix
   ]uppertriang=: <:/~ (i.#m)
1 1 1 1
0 1 1 1
0 0 1 1
0 0 0 1
   ]ixsIntact=: (,uppertriang) # (i.#,m)
0 1 2 3 5 6 7 10 11 15
   ]ixsToZero=: (<(<<ixsIntact)) { (i.#,m)
4 8 9 12 13 14
   ($m) $ 0 ixsToZero } ,m
100 101 102 103
  0 105 106 107
  0   0 110 111
  0   0   0 115

  NB. (b) upper Hessenberg matrix
   upperHessenberg=: 3 : '(0 { y) < (1 { y) + 2'
   upperHessenberg "1 (toIxs m)
1 1 1 1
1 1 1 1
0 1 1 1
0 0 1 1
   ]ixsIntact=: (,(upperHessenberg "1 (toIxs m))) # (i.#,m)
0 1 2 3 4 5 6 7 9 10 11 14 15
   ]ixsToZero=: (<(<<ixsIntact)) { (i.#,m)
8 12 13
   ($m) $ 0 ixsToZero } ,m
100 101 102 103
104 105 106 107
  0 109 110 111
  0   0 114 115

  NB. (b) upper cross triangular matrix
   upperCrossTriang=: 4 : '(#x) > (0 { y) + (1 { y)'
   m upperCrossTriang "1 (toIxs m)
1 1 1 1
1 1 1 0
1 1 0 0
1 0 0 0
   ]ixsIntact=: (,(m upperCrossTriang "1 (toIxs m))) # (i.#,m)
0 1 2 3 4 5 6 8 9 12
   ]ixsToZero=: (<(<<ixsIntact)) { (i.#,m)
7 10 11 13 14 15
   ($m) $ 0 ixsToZero } ,m
100 101 102 103
104 105 106   0
108 109   0   0
112   0   0   0
```

### Solution to exercise 9
```
   ]t=: (m ,: m+100) , m+200
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115

200 201 202 203
204 205 206 207
208 209 210 211
212 213 214 215

300 301 302 303
304 305 306 307
308 309 310 311
312 313 314 315
   $t
3 4 4
   neighbors3D=: 3 : '( ((0{y),(1{y)-1),(2{y) ),((0{y),((1{y)+1),(2{y) ),( (0{y),(1{y),((2{y)-1) ),( (0{y),(1{y),((2{y)+1) ), ( ((0{y)+1),(1{y),(2{y) ),:( ((0{y)-1),(1{y),(2{y) )'
   ]n111=: neighbors3D 1 1 1
1 0 1
1 2 1
1 1 0
1 1 2
2 1 1
0 1 1
   ]n011=: neighbors3D 0 1 1
 0 0 1
 0 2 1
 0 1 0
 0 1 2
 1 1 1
_1 1 1
   ]maxP=. (0{$t) - 1
2
   ]maxR=. (1{$t) - 1
3
   ]maxC=. (2{$t) - 1
3
   validateNeighbors=: 3 : '((0{y) >: 0) *. ((0{y) <: maxP) *. ((1{y) >: 0) *. ((1{y) <: maxR) *. ((2{y) >: 0) *. ((2{y) <: maxC)'
   validateNeighbors (1 1 1)
1
   validateNeighbors (_1 1 1)
0
   validateNeighbors"1 n111
1 1 1 1 1 1
   validateNeighbors"1 n011
1 1 1 1 1 0
   ]n111Filtered=: (validateNeighbors"1 n111) # (i.(0{$n111)) { n111
1 0 1
1 2 1
1 1 0
1 1 2
2 1 1
0 1 1
   ]n011Filtered=: (validateNeighbors"1 n011) # (i.(0{$n011)) { n011
0 0 1
0 2 1
0 1 0
0 1 2
1 1 1
   calcIxs=: 3 : '(2{y) + ((1{y)*(maxC+1)) + ((0{y)*(maxR+1)*(maxC+1))'
   calcIxs 0 0 3
3
   calcIxs 0 2 3
11
   calcIxs 1 2 3
27
   ($t) $ 1 (calcIxs"1 n111Filtered) } (,($t) $ 0)
0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0

0 1 0 0
1 0 1 0
0 1 0 0
0 0 0 0

0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0
   ($t) $ 1 (calcIxs"1 n011Filtered) } (,($t) $ 0)
0 1 0 0
1 0 1 0
0 1 0 0
0 0 0 0

0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0

0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
   ]n111Ixs=: calcIxs"1 n111Filtered
17 25 20 22 37 5
   ]n111Vals=: (calcIxs"1 n111Filtered) { ,t
201 209 204 206 305 105
   ]t111=:($t) $ (<: n111Vals) n111Ixs } ,t
100 101 102 103
104 104 106 107
108 109 110 111
112 113 114 115

200 200 202 203
203 205 205 207
208 208 210 211
212 213 214 215

300 301 302 303
304 304 306 307
308 309 310 311
312 313 314 315
   t - t111
0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0

0 1 0 0
1 0 1 0
0 1 0 0
0 0 0 0

0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0

   ]n011Ixs=: calcIxs"1 n011Filtered
1 9 4 6 21
   ]n011Vals=: (calcIxs"1 n011Filtered) { ,t
101 109 104 106 205
   ]t011=:($t) $ (<: n011Vals) n011Ixs } ,t
100 100 102 103
103 105 105 107
108 108 110 111
112 113 114 115

200 201 202 203
204 204 206 207
208 209 210 211
212 213 214 215

300 301 302 303
304 305 306 307
308 309 310 311
312 313 314 315
   t - t011
0 1 0 0
1 0 1 0
0 1 0 0
0 0 0 0

0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0

0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
```

### Solution to exercise 10
```
   NB. (a) Show capability to update (increment) the second nearest neighbors in a matrix.
   m
100 101 102 103
104 105 106 107
108 109 110 111
112 113 114 115
  NB. For (x,y) we will update 4 points: (x-1,y-1),(x-1,y+1),(x+1,y-1),(x+1,y+1)
   nnn=: 3 : '( ((0{y)-1), ((1{y)-1) ),( ((0{y)-1), ((1{y)+1) ),( ((0{y)+1), ((1{y)-1) ),:( ((0{y)+1), ((1{y)+1) )'
   ]nnn22=: nnn 2 2
1 1
1 3
3 1
3 3
   ]nnn00=: nnn 0 0
_1 _1
_1  1
 1 _1
 1  1
   validateNeighbors=: 3 : '((0{y) >: 0) *. ((0{y) <: maxR) *. ((1{y) >: 0) *. ((1{y) <: maxC)'
   validateNeighbors"1 nnn00
0 0 0 1
   validateNeighbors"1 nnn22
1 1 1 1
   ]nnn00Filtered=: (validateNeighbors"1 nnn00) # (i.(0{$nnn00)) { nnn00
1 1
   ]nnn22Filtered=: (validateNeighbors"1 nnn22) # (i.(0{$nnn22)) { nnn22
1 1
1 3
3 1
3 3
   ($m) $ 1 (calcIxs"1 nnn00Filtered) } (,($m) $ 0)
0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0
   ]nnn00Ixs=: calcIxs"1 nnn00Filtered
5
   ]nnn00Vals=: (calcIxs"1 nnn00Filtered) { ,m
105
   ]m00=:($m) $ (>: nnn00Vals) nnn00Ixs } ,m
100 101 102 103
104 106 106 107
108 109 110 111
112 113 114 115
   m00 - m
0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0

   ($m) $ 1 (calcIxs"1 nnn22Filtered) } (,($m) $ 0)
0 0 0 0
0 1 0 1
0 0 0 0
0 1 0 1
   ]nnn22Ixs=: calcIxs"1 nnn22Filtered
5 7 13 15
   ]nnn22Vals=: (calcIxs"1 nnn22Filtered) { ,m
105 107 113 115
   ]m22=:($m) $ (>: nnn22Vals) nnn22Ixs } ,m
100 101 102 103
104 106 106 108
108 109 110 111
112 114 114 116
   m22 - m
0 0 0 0
0 1 0 1
0 0 0 0
0 1 0 1


   NB. (b) Also show updating simultaneuously the nearest neighbors of two points at the same time (both incrementing, then one incrementing second decrementing)
   ]twopointsA=:($m) $ (>: nnn22Vals, >: nnn00Vals) (nnn22Ixs, nnn00Ixs) } ,m
100 101 102 103
104 107 106 108
108 109 110 111
112 114 114 116
   twopointsA - m
0 0 0 0
0 2 0 1
0 0 0 0
0 1 0 1
   ]twopointsB=:($m) $ (>: nnn22Vals, <: nnn00Vals) (nnn22Ixs, nnn00Ixs) } ,m
100 101 102 103
104 105 106 108
108 109 110 111
112 114 114 116
   twopointsB - m
0 0 0 0
0 0 0 1
0 0 0 0
0 1 0 1

   NB. (c) Also show updating simultaneuously the nearest neighbors (incrementing) and second nearest neighbors (decrementing)
   NB. We are going to reuse the (a) case and the one from the main text, and then combine them.
   NB. the nearest neighbors
   ($m) $ 1 (calcIxs"1 nn22Filtered) } (,($m) $ 0)
0 0 0 0
0 0 1 0
0 1 0 1
0 0 1 0
   ($m) $ 1 (calcIxs"1 nn00Filtered) } (,($m) $ 0)
0 1 0 0
1 0 0 0
0 0 0 0
0 0 0 0
   ]nn22Ixs=: calcIxs"1 nn22Filtered
6 14 9 11
   ]nn00Ixs=: calcIxs"1 nn00Filtered
4 1
   ]nn22Vals=: (calcIxs"1 nn22Filtered) { ,m
106 114 109 111
   ]nn00Vals=: (calcIxs"1 nn00Filtered) { ,m
104 101

   NB. the second nearest neighbors
   ($m) $ 1 (calcIxs"1 nnn00Filtered) } (,($m) $ 0)
0 0 0 0
0 1 0 0
0 0 0 0
0 0 0 0
   ($m) $ 1 (calcIxs"1 nnn22Filtered) } (,($m) $ 0)
0 0 0 0
0 1 0 1
0 0 0 0
0 1 0 1
   ]nnn00Ixs=: calcIxs"1 nnn00Filtered
5
   ]nnn22Ixs=: calcIxs"1 nnn22Filtered
5 7 13 15
   ]nnn00Vals=: (calcIxs"1 nnn00Filtered) { ,m
105
   ]nnn22Vals=: (calcIxs"1 nnn22Filtered) { ,m
105 107 113 115

   ]twoneighbors00=:($m) $ ((>: nn00Vals), (<: nnn00Vals)) (nn00Ixs, nnn00Ixs) } ,m
100 102 102 103
105 104 106 107
108 109 110 111
112 113 114 115
   twoneighbors00 - m
0  1 0 0
1 _1 0 0
0  0 0 0
0  0 0 0

   ]twoneighbors22=:($m) $ ((>: nn22Vals), (<: nnn22Vals)) (nn22Ixs, nnn22Ixs) } ,m
100 101 102 103
104 104 107 106
108 110 110 112
112 112 115 114
   twoneighbors22 - m
0  0 0  0
0 _1 1 _1
0  1 0  1
0 _1 1 _1
```

### Solution to exercise 11
```
   NB. x is point where diagonals cross, y is matrix where diagonal generation acts
   diags=: 4 : '( ((0{y) + (1{x)) = ((1{y) + (0{x)) ) +. ( ((0{x)+(1{x)) = ((0{y) + (1{y)) )'
   2 3 diags"1(toIxs m)
0 1 0 0
0 0 1 0
0 0 0 1
0 0 1 0
   1 1 diags"1(toIxs m)
1 0 1 0
0 1 0 0
1 0 1 0
0 0 0 1
   0 0 diags"1(toIxs m)
1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1
   1 2 diags"1(toIxs m)
0 1 0 1
0 0 1 0
0 1 0 1
1 0 0 0

   ]ixsToZero=: (,(1 2 diags"1 (toIxs m))) # (i.#,m)
1 3 6 9 11 12
   ($m) $ 0 ixsToZero } ,m
100   0 102   0
104 105   0 107
108   0 110   0
  0 113 114 115
```
