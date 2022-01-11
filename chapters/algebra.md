# Linear algebra

## Contents
1. [Selecting from a matrix](#selecting-from-matrix)
2. [Updating a matrix](#updating-of-matrix)
3. [Generate a random matrix](#generating-random-matrix)
4. [Testing matrix properties](#testing-matrix-properties)
5. [Elementary operations of a matrix](#elementary-operations-in-matrix)
6. [Transpose of a matrix](#transpose-of-matrix)
7. [Inverse of a matrix](#inverse-of-matrix) - TODO
8. [Determinant of a matrix](#determinant-of-matrix) - TODO
9. [Trace of a matrix](#trace-of-matrix) - TODO
10. [A partitioned matrix](#partitioned-matrix) - TODO
11. [Matrix decompositions](#matrix-decompositions) - TODO


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
Let's start to revisit what are basic functionalities in J when it comes to vector random generation.
When we want to pick N natural numbers from 0 up to M-1, then it can be achieved via `?N#M`:
```
   ?6#6
1 3 3 2 0 2
   ?6#6
5 2 1 5 2 0
   ?6#6
2 5 5 0 1 5
   ?10#6
2 3 3 1 1 5 4 3 4 4
   ?10#6
0 0 3 1 1 1 2 3 5 2
   1+?6#6
3 6 6 3 1 4
   1+?6#6
3 5 3 6 1 3
```
In the first three examples we picked 6 numbers from the set {0,1,2,3,4,5}, the next two 10 numbers from the same set.
In the last example, we picked 6 numbers from the set {1,2,3,4,5,6}.
But what if we want to choose the set in completely arbitrary way? Well, we can exploit selecting capability we have already mastered:
```
   domain=: _1 1
   ( 10 ?@$ #domain) { domain
_1 1 1 1 1 _1 1 1 _1 _1
   ( 10 ?@$ #domain) { domain
1 1 1 _1 _1 _1 1 _1 _1 1

   domain=: 1 10 100 1000
   ( 10 ?@$ #domain) { domain
100 100 10 1 10 1000 100 10 100 10
   ( 10 ?@$ #domain) { domain
1 1000 1 100 100 1 1 1000 1000 100
```

These were **random sampling with replacement**. If we want to pick a sequence of elements from a given domain **without replacement** we can use `?` with
the remark that the number of picked numbers cannot exceed the cardinality of the domain:
```
   domain=: 1 10 100 1000
   ( 1 ? #domain) { domain
100
   ( 2 ? #domain) { domain
1000 1
   ( 3 ? #domain) { domain
10 1000 100
   ( 4 ? #domain) { domain
100 10 1 1000
   ( 5 ? #domain) { domain
|domain error
|   (5    ?#domain){domain

```

Now what if we want to set relative weights to the elements of domain? On the one hand we can replicate accordingly elements in the domain.
So if we want to have 100 picked 3 times more frequently than any other number in the domain we could set domain the following:
```
   domain=: 1 10 100 100 100 1000
```

The following reference [https://code.jsoftware.com/wiki/Fifty_Shades_of_J/Chapter_28] suggests:
```
   cumulativeWeights=: +/\ % +/
   cumulativeWeights 1 1 3 1
0.166667 0.333333 0.833333 1
   rnd=: ?@#&0
   rnd 5
0.797844 0.357451 0.817211 0.397167 0.160679
   rndWeighted=: cumulativeWeights@[ I. rnd@]
   1 1 3 1 rndWeighted 10
2 2 2 3 1 2 3 2 3 1
```
`rnd` picks random number from (0,1) interval, the last line picks 10 random numbers from a set {0,1,2,3} where
element 3 has relative weight 3 with respect to other elements of the domain that all have relative weight equal to 1.

Now, let's see how we can use arbitrary domain using both approaches.
```
   domain=: 1 10 100 1000
   (1 1 3 1 rndWeighted 10) { domain
1000 1 100 100 1000 100 10 1 100 10
   domain=: 1 10 100 100 100 1000
   ( 10 ?@$ #domain) { domain
1000 1000 100 1000 100 100 1000 100 100 100
```
**Exercise 12**
Show that the both above-mentioned approaches give approximately the same result when the number of picks is substantial.
In the both approaches element 100 has 3 times bigger frequency than the rest elements of the domain, ie., element 100 should occur
approximately 50% time, the rest three elements should be equally frequent.

[Solution to exercise 12](#solution-to-exercise-12)

The exercise 12 is of great importance as it elucidates one of the possible powerful strategies that we will use later to show that our assumptions are valid
or to confirm experimentally the mathematical relations. Basically the strategy bogs down to repeating million or so times the toss, utilizing the randomness of number generation,
proper counting of the resultant observations, correct aggregation of them and drawing the proper conclusion.
It is great asset of J that such massive experiments are up for grabs for us. Moreover, we will encounter numerous situations that the simulation act not just as
a proxy for lemma or theorem or just some finding, which is very reassuring, but sometimes is the only quick way to get to the result as analytical solution
is very nontrivial or only intricate approximation can be provided.

Now we are empowered to replicate **binomial distribution**. We need to set binary domain and weights being probabilities,
ie. two positive numbers that add up to 1. Traditionally, domain reflects **Bernoulli trial** which is experiment with two outcomes
possible, 1 representing success with probability p, 0 representing failure with probability (1-p) [6, page 89-91].
We can generalize domain though.
```
   cumulativeWeights=: +/\ % +/
   rnd=: ?@#&0
   rndWeighted=: cumulativeWeights@[ I. rnd@]
   domain=: 0 1
   (0.3 0.7 rndWeighted 10) { domain
1 1 1 0 1 1 0 1 0 1
   domain=: _1 1
   (0.15 0.85 rndWeighted 10) { domain
1 1 1 _1 1 1 1 1 1 1
```
**Exercise 13**
Explore two generalized binomial distributions:
(a) (1,p) and (0,1-p) where p is (0,1),
(b) (1, p) and (22, 1-p) where p is (0,1).
Demonstrate that if the number of tials is big enough the weights (ie. probabilities) are replicated.
Calculate experimental mean and variance of Bernoulli random variables and compare to the theoretical results.

[Solution to exercise 13](#solution-to-exercise-13)

Moreover, random numbers 0 and 1 in binomial distribution can be obtained via `binomialrand` - see [https://code.jsoftware.com/wiki/Addons/stats/base/random]
```
   load 'stats/base/random'
   NB. probability of success=0.2, number of trials 10
   binomialrand 0.2 10
0 0 0 0 1 0 0 0 0 0
```

Let's now investigate two continuous distributions: normal and uniform. A **normal distribution** example is below:
```
   NB. `rnorm` is defined in j/algebra.ijs and takes as x mean and variance, and number of samples as y
   0 1 rnorm 10
_0.22246 0.565404 _0.81757 _1.44307 1.37019 1.32798 _0.325787 0.85836 _0.586362 0.751552

   10 2 rnorm 10
7.79113 7.16799 12.1581 8.78351 10.3067 9.38921 11.7871 11.0787 8.18868 10.8973
```
We can now see using `intervalHist` from j/algebra.ijs how the generated samples are distributed:
```
   ]bins=: 0.2*i:15
_3 _2.8 _2.6 _2.4 _2.2 _2 _1.8 _1.6 _1.4 _1.2 _1 _0.8 _0.6 _0.4 _0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3
   bins intervalHist (0 1 rnorm 100)
┌────────┬─────┬────┐
│interval│count│freq│
├────────┼─────┼────┤
│  _3    │ 0   │   0│
│_2.8    │ 0   │   0│
│_2.6    │ 0   │   0│
│_2.4    │ 0   │   0│
│_2.2    │ 1   │0.01│
│  _2    │ 2   │0.02│
│_1.8    │ 0   │   0│
│_1.6    │ 1   │0.01│
│_1.4    │ 2   │0.02│
│_1.2    │ 3   │0.03│
│  _1    │ 6   │0.06│
│_0.8    │ 5   │0.05│
│_0.6    │ 5   │0.05│
│_0.4    │ 9   │0.09│
│_0.2    │13   │0.13│
│   0    │12   │0.12│
│ 0.2    │10   │ 0.1│
│ 0.4    │ 5   │0.05│
│ 0.6    │ 6   │0.06│
│ 0.8    │ 6   │0.06│
│   1    │ 5   │0.05│
│ 1.2    │ 5   │0.05│
│ 1.4    │ 2   │0.02│
│ 1.6    │ 0   │   0│
│ 1.8    │ 1   │0.01│
│   2    │ 0   │   0│
│ 2.2    │ 1   │0.01│
│ 2.4    │ 0   │   0│
│ 2.6    │ 0   │   0│
│ 2.8    │ 0   │   0│
│   3    │ 0   │   0│
└────────┴─────┴────┘

   bins intervalHist (0 1 rnorm 1e6)
┌────────┬─────┬──────────┐
│interval│count│freq      │
├────────┼─────┼──────────┤
│  _3    │ 1337│  0.001337│
│_2.8    │ 1191│  0.001191│
│_2.6    │ 2124│  0.002124│
│_2.4    │ 3577│  0.003577│
│_2.2    │ 5733│0.00573301│
│  _2    │ 8883│0.00888301│
│_1.8    │13485│  0.013485│
│_1.6    │18618│  0.018618│
│_1.4    │26150│   0.02615│
│_1.2    │34150│   0.03415│
│  _1    │43414│  0.043414│
│_0.8    │53318│ 0.0533181│
│_0.6    │62517│ 0.0625171│
│_0.4    │70226│ 0.0702261│
│_0.2    │76195│ 0.0761951│
│   0    │78716│ 0.0787161│
│ 0.2    │79239│ 0.0792391│
│ 0.4    │76277│ 0.0762771│
│ 0.6    │70208│ 0.0702081│
│ 0.8    │62095│ 0.0620951│
│   1    │53518│ 0.0535181│
│ 1.2    │43695│  0.043695│
│ 1.4    │34461│  0.034461│
│ 1.6    │26138│  0.026138│
│ 1.8    │18898│  0.018898│
│   2    │12913│  0.012913│
│ 2.2    │ 8858│0.00885801│
│ 2.4    │ 5760│0.00576001│
│ 2.6    │ 3556│  0.003556│
│ 2.8    │ 2120│   0.00212│
│   3    │ 1175│  0.001175│
│        │ 1454│  0.001454│
└────────┴─────┴──────────┘
```

**Exercise 14**
Show that 1-,2-, 3- sigma interval probabilities can be reasonable assessed using `rnorm`.

[Solution to exercise 14](#solution-to-exercise-14)


Finally, we can look at **uniform distribution** sample generation (`runiform` is in j/algebra.ijs ):
```
   NB. 10 samples of U(0,1)
   0 1 runiform 10
0.183411 0.0968962 0.587723 0.165308 0.68218 0.0916652 0.00554653 0.149567 0.340257 0.370271

   ]bins=.(0.2&*) <: i.8
_0.2 0 0.2 0.4 0.6 0.8 1 1.2
   bins intervalHist (0 1 runiform 100)
┌────────┬─────┬────┐
│interval│count│freq│
├────────┼─────┼────┤
│_0.2    │ 0   │   0│
│   0    │ 0   │   0│
│ 0.2    │23   │0.23│
│ 0.4    │17   │0.17│
│ 0.6    │17   │0.17│
│ 0.8    │21   │0.21│
│   1    │22   │0.22│
│ 1.2    │ 0   │   0│
└────────┴─────┴────┘

   bins intervalHist (0 1 runiform 1e6)
┌────────┬──────┬────────┐
│interval│count │freq    │
├────────┼──────┼────────┤
│_0.2    │     0│       0│
│   0    │     0│       0│
│ 0.2    │199601│0.199601│
│ 0.4    │200025│0.200025│
│ 0.6    │199935│0.199935│
│ 0.8    │199972│0.199972│
│   1    │200467│0.200467│
│ 1.2    │     0│       0│
└────────┴──────┴────────┘
```

**Exercise 15**
Generate 10x10 random upper triangular matrix where elements are N(2,3).

[Solution to exercise 15](#solution-to-exercise-15)

**Exercise 16**
Generate 10x10 random diagonal matrix where elements are U(20,30).

[Solution to exercise 16](#solution-to-exercise-16)

**Exercise 17**
Generate 8x5 matrix where consecutive columns consist of random and evenly probable pairs: (1,2), (3,4), (5,6),...

[Solution to exercise 17](#solution-to-exercise-17)

In-depth coverage of many both discrete and continuous distribution families will be included in statistics inference chapter.

The last topics to cover when introducing basic random generation are **random generators** and **seeds**.
We have the following random generators at our disposal (adapted from [https://code.jsoftware.com/wiki/Vocabulary/query])

|     code      |      rng           | relative cost  |
| ------------- |:------------------:| --------------:|
|       1       |    GB_Flip         |        1       |
|       2       | Mersenne Twister   |        1       |
|       3       |    DX-1597-4d      |        3       |
|       4       |     MRG32k3a       |        8       |
|       0       | combination of all |       12       |

Below is self-explanatory code snippet.
```
   NB. show current rng (Mersenne Twister is default)
   9!:42 ''
2
   NB. set new rng
   9!:43 ]1

   9!:42 ''
1

   rng=.9!:42
   rng ''
1
   NB. set seed=2000 to current rng
   rngWithSeed2000=. 2{.2000,rng ''
   9!:43 {:rngWithSeed2000

   NB. reset the state of rng with seed
   9!:1  {.rngWithSeed2000

   NB. generating 5 samples of N(0,1)
   0 1 rnorm 5
_1.42716 _0.353494 0.0569464 _0.300366 0.752976
   0 1 rnorm 5
0.58307 _1.25819 1.05434 1.09644 _0.294877

   NB. reset the state of rng with seed
   9!:1  {.rngWithSeed2000

   NB. now 5 sample of N(0,1) should be the same as immediately after previous rng's state reset
   NB. the rng's state reset functionality will be important later for experiment's reproducibility
   0 1 rnorm 5
_1.42716 _0.353494 0.0569464 _0.300366 0.752976
```

The above is functionality is enough to have a basic control of random generation.
More information can be found here [https://code.jsoftware.com/wiki/Essays/RNG].

**Summary**: The basic random generation capabilities of J were covered. We know how to
(a) toss with repetition and without repetition,
(b) set arbitrary domain and give weights to elements of the domain,
(c) pick random samples from normal and uniform distributions,
(d) get frequency report of the generated samples,
(e) use substantial number generation to reason about properties of distributions, like mean or variance,
(f) set random number generator with default or arbitrary seed,
(g) reset the state of the generator.

## Testing matrix properties
In the coming chapters we will develop many techniques and recipies, and to have reasonable confidence the proposed
solution is correct we will adapt **property testing**. The scheme I will adopt is following:
1. Implement a concept **C** (eg. transpose, SVD, ...)
2. Refer to the facts, formulas, lemmas and proofs of mathematics and construct **leftR R rightR**
Here both leftR and rightR can contain the concept **C** (and possibly others) and establish relation **R** (eg. =, <=, ...)
3. As both **leftR** and **rightR**, in general, act on sequence of arrays we will need to deliver it. Very often they will need to be special
arrays, due to the shape or the array type constraints
4. Rather than handcrafting the arrays we will rely on the generated array instances. The developments of previous section will be very useful indeed
5. We will repeat the experiment many times, with the expectation that in every experiment the property we are verifying
holds

This is a standard procedure, for example in Haskell development, were we construct a property, implement `Arbitrary` instances, and then
upon property testing, proper array instances are generated and the property is tried with them. As I am convinced this is the proper and
the required approach I will adopt it as well here.

## Elementary operations in matrix
We have the following basic results.

### Matrix addition
Let's have matrices of the same order: A, B, C and scalars s<sub>1</sub> and s<sub>2</sub>. Then we have [2, pages 5]:
- <img src="https://latex.codecogs.com/svg.image?A&space;&plus;&space;B&space;=&space;B&space;&plus;&space;A" title="A + B = B + A" />
- <img src="https://latex.codecogs.com/svg.image?(A&space;&plus;&space;B)&space;&plus;&space;C&space;=&space;A&space;&plus;&space;(B&space;&plus;&space;C)" title="(A + B) + C = A + (B + C)" />
- <img src="https://latex.codecogs.com/svg.image?(s_1&space;&plus;&space;s_2)A&space;=&space;s_1&space;A&space;&plus;&space;s_2&space;A" title="(s_1 + s_2)A = s_1 A + s_2 A" />
- <img src="https://latex.codecogs.com/svg.image?s_1&space;(A&space;&plus;&space;B)&space;=&space;s_1&space;A&space;&plus;&space;s_1&space;B" title="s_1 (A + B) = s_1 A + s_1 B" />
- <img src="https://latex.codecogs.com/svg.image?s_1&space;(s_2&space;A)&space;=&space;(s_1&space;s_2)A" title="s_1 (s_2 A) = (s_1 s_2)A" />
- <img src="https://latex.codecogs.com/svg.image?A&space;&plus;&space;(-1)A&space;=&space;0&space;" title="A + (-1)A = 0 " />

Let's develop how we can perform property test of the first equality. According to the scheme proposed above we need to:
1. have a way to generate two arbitrary matrices of the same shape
```
genUniformMatrix=: 3 : 'y $ _1000 1000 runiform ((0{y) * (1{y))'
   genUniformMatrix 2 2
 _226.76 _322.805
_808.466  202.957
   genUniformMatrix 4 2
_428.685  853.433
_400.652 _164.792
 375.372  675.547
 584.175 _69.8546
```
2. have a way to check **leftR R rightR** with the generated matrices
```
   leftR=: 4 : 'x + y'
   rightR=: 4 : 'y + x'
   relation=: leftR`rightR
      relation@.0
4 : 'x + y'
      relation@.1
4 : 'y + x'

   checkEqTwoMatrices=: 4 : '( (0{y) x@.0 (1{y) ) = ( (0{y) x@.1 (1{y) )'
   ]matrices=: (genUniformMatrix 4 2),:(genUniformMatrix 4 2)
 783.326  777.188
 433.257  992.401
_44.5578   892.71
 850.185   636.44

 211.161 _619.827
 464.316 _601.967
 309.364  851.114
_181.237  238.782
   NB. all elements the same
   relation checkEqTwoMatrices matrices
1 1
1 1
1 1
1 1
   (0{matrices) (relation@.0) (1{matrices)
994.487  157.36
897.573 390.434
264.806 1743.82
668.948 875.221
   (0{matrices) (relation@.1) (1{matrices)
994.487  157.36
897.573 390.434
264.806 1743.82
668.948 875.221
   NB. let's redefine check in such a way that it returns 1 only if all elements are the same,
   NB. ie. we have perfect matching
   checkEqTwoMatrices=: 4 : '( (0{y) x@.0 (1{y) ) -: ( (0{y) x@.1 (1{y) )'
   relation checkEqTwoMatrices matrices
1
```
3. finally, repeat the check many times for different shapes (as an example dimensions are independently picked from 1 ... 100 domain)
```
   run=: 3 : 0
shape=.1+?2#100
m=.(genUniformMatrix shape),:(genUniformMatrix shape)
relation checkEqTwoMatrices m
)
   NB. now let's run it 100 times
   (+/)(run"0)100#0
100
```
Now, I am reasonably confident that matrix addition as implemented in J is a commutative operator, ie. <img src="https://latex.codecogs.com/svg.image?A&space;&plus;&space;B&space;=&space;B&space;&plus;&space;A" title="A + B = B + A" /> holds. I intend to verify each property like that onwards.

**Exercise 18**
Perform property testing for the rest addition properties specified above

[Solution to exercise 18](#solution-to-exercise-18)

#### Note on floating point arithmetics.

Let's revisit the following property <img src="https://latex.codecogs.com/svg.image?s_1&space;(A&space;&plus;&space;B)&space;=&space;s_1&space;A&space;&plus;&space;s_1&space;B" title="s_1 (A + B) = s_1 A + s_1 B" />


```
   leftR=: 4 : '(0{x) * ( (0{y) + (1{y) )'
   rightR=: 4 : '( (0{x) * (0{y) ) + ( (0{x) * (1{y) )'
   relation=: leftR`rightR

   checkEqOfMatricesScalarsRel=: 4 : '( (>0{y) x@.0 (>1{y) ) -: ( (>0{y) x@.1 (>1{y) )'

   run=: 3 : 0
shape=.1+?2#100
m=.(genUniformMatrix shape),:(genUniformMatrix shape)
s=. _100 100 runiform 1
data=.s;m
relation checkEqOfMatricesScalarsRel data
)

   (+/)(run"0)100#0
43
```
This means in this batch run 57 sample cases failed to validate this property. The reason is that floating-point addition is not associative or distributive.
When testing the properties we would like to have a way to (1) detect the failing cases, (2) avoid the floating point deficiencies.

TODO detect failing cases

In order to harness floating point we will try two approaches. In the first we will make decrease strictness of tolerance.
```
   9!:18 ''
5.68434e_14
   9!:19 ]1e_11

   9!:18 ''
1e_11
   (+/)(run"0)100#0
100
   (+/)(run"0)1000#0
991
```
Making equal tolerance less strict helped a lot, but we are still not perfect.

### Matrix multiplication

We have the following basic properties of matrix multiplication:
- <img src="https://latex.codecogs.com/svg.image?(AB)C=A(BC)" title="(AB)C=A(BC)" />
- <img src="https://latex.codecogs.com/svg.image?A(B&plus;C)=AB&plus;AC" title="A(B+C)=AB+AC" />
- <img src="https://latex.codecogs.com/svg.image?(A&plus;B)C=AC&plus;BC" title="(A+B)C=AC+BC" />

The matrix multiplcation is defined as
```
   mult=: +/ .*
   ]a=: 2 2 $ 1 2 3 4
1 2
3 4
   ]b=: 2 4 $ 1 2
1 2 1 2
1 2 1 2
   ]c=: a mult b
3  6 3  6
7 14 7 14
   $c
2 4
```

**Exercise 19**
Compute the n-th Fibonacci number by using the matrix form [https://en.wikipedia.org/wiki/Fibonacci_number#Matrix_form]

[Solution to exercise 19](#solution-to-exercise-19)

**Exercise 20**
Perform property testing for the multiplication properties specified above

[Solution to exercise 20](#solution-to-exercise-20)

### Matrix elementary row and column operations

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

**Exercise 21**
Show that the three basic operations can be realized by multiplication of transformed identity matrices.

[Solution to exercise 21](#solution-to-exercise-21)

### Transpose of matrix
The transpose is defined as follows:
```
   tr=: |:
   tr (2 3 $ 1 2 3 4 5 6)
1 4
2 5
3 6
   tr (tr (2 3 $ 1 2 3 4 5 6))
1 2 3
4 5 6
```
We have also the following properties
- <img src="https://latex.codecogs.com/svg.image?(A^T)^T=A" title="(A^T)^T=A" />
- <img src="https://latex.codecogs.com/svg.image?(A&plus;B)^T=A^T&plus;B^T" title="(A+B)^T=A^T+B^T" />
- <img src="https://latex.codecogs.com/svg.image?(AB)^T=B^TA^T" title="(AB)^T=B^TA^T" />

**Exercise 22**
Perform property testing for transpose properties.

[Solution to exercise 22](#solution-to-exercise-22)

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

### Solution to exercise 12
```
   domain=: 1 10 100 100 100 1000
   runs=: (10000000 ?@$ #domain) { domain
   discreteHist runs
┌────┬────────┐
│elem│freq    │
├────┼────────┤
│   1│0.166741│
│  10│0.166661│
│ 100│0.499965│
│1000│0.166633│
└────┴────────┘

   cumulativeWeights=: +/\ % +/
   rnd=: ?@#&0
   rndWeighted=: cumulativeWeights@[ I. rnd@]
   domain=: 1 10 100 1000
   runs=: (1 1 3 1 rndWeighted 10000000) { domain
   discreteHist runs
┌────┬────────┐
│elem│freq    │
├────┼────────┤
│   1│ 0.16648│
│  10│0.166672│
│ 100│0.500081│
│1000│0.166768│
└────┴────────┘

  NB. discreteHist is function defined in j/algebra.ijs
```

### Solution to exercise 13
```
   domain=: 0 1
   discreteHist ((0.15 0.85 rndWeighted 1e6) { domain)
┌────┬────────┐
│elem│freq    │
├────┼────────┤
│0   │0.149227│
│1   │0.850773│
└────┴────────┘

   domain=: 1 22
   discreteHist ((0.35 0.65 rndWeighted 1e6) { domain)
┌────┬────────┐
│elem│freq    │
├────┼────────┤
│ 1  │0.349786│
│22  │0.650214│
└────┴────────┘

   NB. mean of Bernoulli trial for {(val1, p),(val2,p)} is (val1*p + val2(1-p))
   NB. When {(1,p),(0,1-p)} then mean is p
   NB. mean is defined in j/algebra.ijs

   domain=: 0 1
   mean ((0.15 0.85 rndWeighted 1e6) { domain)
0.850207

   domain=: 1 22
   mean ((0.35 0.65 rndWeighted 1e6) { domain)
14.6405
   (0.35*1) + (0.65*22)
14.65

   NB. variance of Bernoulli trial for {(1, p),(0,p)} is p(1-p)
   NB. var is defined in j/algebra.ijs

   domain=: 0 1
   var ((0.15 0.85 rndWeighted 1e6) { domain)
0.127708
   0.15*(1-0.15)
0.1275

   domain=: 1 22
   var ((0.35 0.65 rndWeighted 1e6) { domain)
100.414
```

### Solution to exercise 14
```
   NB. 1-sigma, 2-sigma and 3-sigma for N(0,1)
   bins=: _1 1 _
   bins intervalHist (0 1 rnorm 1e6)
┌────────┬──────┬────────┐
│interval│count │freq    │
├────────┼──────┼────────┤
│_1      │158765│0.158765│
│ 1      │682634│0.682635│
│ _      │158600│  0.1586│
└────────┴──────┴────────┘
   bins=: _2 2 _
   bins intervalHist (0 1 rnorm 1e6)
┌────────┬──────┬────────┐
│interval│count │freq    │
├────────┼──────┼────────┤
│_2      │ 23235│0.023235│
│ 2      │954106│0.954107│
│ _      │ 22658│0.022658│
└────────┴──────┴────────┘
   bins=: _3 3 _
   bins intervalHist (0 1 rnorm 1e6)
┌────────┬──────┬────────┐
│interval│count │freq    │
├────────┼──────┼────────┤
│_3      │  1316│0.001316│
│ 3      │997328│0.997329│
│ _      │  1355│0.001355│
└────────┴──────┴────────┘

    NB. 1-sigma, 2-sigma and 3-sigma for N(10,4)
    bins=: 6 14 _
┌────────┬──────┬────────┐
│interval│count │freq    │
├────────┼──────┼────────┤
│ 6      │159015│0.159015│
│14      │682149│ 0.68215│
│_       │158835│0.158835│
└────────┴──────┴────────┘
   bins=: 2 18 _
   bins intervalHist (10 4 rnorm 1e6)
┌────────┬──────┬────────┐
│interval│count │freq    │
├────────┼──────┼────────┤
│ 2      │ 22878│0.022878│
│18      │954398│0.954398│
│ _      │ 22724│0.022724│
└────────┴──────┴────────┘
   bins=: _2 22 _
   bins intervalHist (10 4 rnorm 1e6)
┌────────┬──────┬────────┐
│interval│count │freq    │
├────────┼──────┼────────┤
│_2      │  1318│0.001318│
│22      │997293│0.997293│
│ _      │  1389│0.001389│
└────────┴──────┴────────┘
```

### Solution to exercise 15
```
   ]uppertriang=: (<:)/~ (i.10)
1 1 1 1 1 1 1 1 1 1
0 1 1 1 1 1 1 1 1 1
0 0 1 1 1 1 1 1 1 1
0 0 0 1 1 1 1 1 1 1
0 0 0 0 1 1 1 1 1 1
0 0 0 0 0 1 1 1 1 1
0 0 0 0 0 0 1 1 1 1
0 0 0 0 0 0 0 1 1 1
0 0 0 0 0 0 0 0 1 1
0 0 0 0 0 0 0 0 0 1
   ]rsquare=: 10 10 $ 2 3 rnorm 100
     4.0108  3.72568   4.08369   5.47384  0.501336   2.13195   1.68755 _0.885756  2.23265 _2.12981
   _1.22991 _1.65395   0.12277   1.72872  _2.41367   1.16439     2.086  _4.52652  6.96851    2.817
0.000702556  4.14621   3.75542   2.08605  0.637916    5.0991 _0.118078 _0.397692  1.60716 0.213515
    4.55878 _2.63454   _0.4625   4.41528 _0.301713 _0.931771   5.11125  _8.55276  2.81113  0.46742
    4.13462  1.85034  _3.68801  _1.34194   2.32792   6.77218   5.74401  _4.34345  3.68939  1.26638
    5.09811  4.95374 _0.755182   1.31864   3.74399  _3.20322  0.356599   1.77643  2.46348  3.91575
 _0.0623435 0.536959  0.621922 0.0958049  _4.42459   3.59502    1.1154   5.66886 0.428384  2.67899
    4.24616   4.5828  0.819118   5.94708   5.75729  _2.91884  _1.49414   3.10644   2.4001  4.69477
    3.58613 0.022629   1.10225  _4.06354   3.86263  _2.09862   6.52759   6.45744 _3.94655  8.84773
     1.5702 _1.26945   4.18487 _0.629048  _1.81789  _2.27127   6.02709   2.48604  6.65651 _2.52068

   rsquare * uppertriang
4.0108  3.72568 4.08369 5.47384  0.501336   2.13195   1.68755 _0.885756  2.23265 _2.12981
     0 _1.65395 0.12277 1.72872  _2.41367   1.16439     2.086  _4.52652  6.96851    2.817
     0        0 3.75542 2.08605  0.637916    5.0991 _0.118078 _0.397692  1.60716 0.213515
     0        0       0 4.41528 _0.301713 _0.931771   5.11125  _8.55276  2.81113  0.46742
     0        0       0       0   2.32792   6.77218   5.74401  _4.34345  3.68939  1.26638
     0        0       0       0         0  _3.20322  0.356599   1.77643  2.46348  3.91575
     0        0       0       0         0         0    1.1154   5.66886 0.428384  2.67899
     0        0       0       0         0         0         0   3.10644   2.4001  4.69477
     0        0       0       0         0         0         0         0 _3.94655  8.84773
     0        0       0       0         0         0         0         0        0 _2.52068
```

### Solution to exercise 16
```
   ]diag=: =/~ (i.10)
1 0 0 0 0 0 0 0 0 0
0 1 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 0 0 0
0 0 0 1 0 0 0 0 0 0
0 0 0 0 1 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0
0 0 0 0 0 0 1 0 0 0
0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 1 0
0 0 0 0 0 0 0 0 0 1
   ]runif=: 20 30 runiform 10
26.2247 23.2471 29.0782 20.6316 23.8662 23.386 20.9577 26.0148 22.8566 29.2672
   diag * runif
26.2247       0       0       0       0      0       0       0       0       0
      0 23.2471       0       0       0      0       0       0       0       0
      0       0 29.0782       0       0      0       0       0       0       0
      0       0       0 20.6316       0      0       0       0       0       0
      0       0       0       0 23.8662      0       0       0       0       0
      0       0       0       0       0 23.386       0       0       0       0
      0       0       0       0       0      0 20.9577       0       0       0
      0       0       0       0       0      0       0 26.0148       0       0
      0       0       0       0       0      0       0       0 22.8566       0
      0       0       0       0       0      0       0       0       0 29.2672
```

### Solution to exercise 17
```
   domain=: 1 2
   ]c1=:(8 ?@$ #domain) { domain
2 1 2 1 1 1 2 2
   domain=: 3 4
   ]c2=:(8 ?@$ #domain) { domain
3 4 4 3 3 3 4 4
   domain=: 5 6
   ]c3=:(8 ?@$ #domain) { domain
5 5 5 5 6 6 6 6
   domain=: 7 8
   ]c4=:(8 ?@$ #domain) { domain
8 8 7 7 7 7 7 7
   domain=: 9 10
   ]c5=:(8 ?@$ #domain) { domain
9 9 9 9 10 10 9 9
   (,.c1),.(,.c2),.(,.c3),.(,.c4),.(,.c5)
2 3 5 8  9
1 4 5 8  9
2 4 5 7  9
1 3 5 7  9
1 3 6 7 10
1 3 6 7 10
2 4 6 7  9
2 4 6 7  9
```

### Solution to exercise 18
- <img src="https://latex.codecogs.com/svg.image?(A&space;&plus;&space;B)&space;&plus;&space;C&space;=&space;A&space;&plus;&space;(B&space;&plus;&space;C)" title="(A + B) + C = A + (B + C)" />
```
   leftR=: 3 : '( (0{y) + (1{y) ) + (2{y)'
   rightR=: 3 : '(0{y) + ( (1{y) + (2{y) )'
   relation=: leftR`rightR
   checkEqOfMatrices=: 4 : '( x@.0 y ) -: ( x@.1 y )'
   genUniformMatrix=: 3 : 'y $ _1000 1000 runiform ((0{y) * (1{y))'
   ]matrices=: (genUniformMatrix 4 2),(genUniformMatrix 4 2),:(genUniformMatrix 4 2)
 211.161 _619.827
 464.316 _601.967
 309.364  851.114
_181.237  238.782

_428.685  853.433
_400.652 _164.792
 375.372  675.547
 584.175 _69.8546

 244.943 _350.586
  815.65 _873.687
 _226.76 _322.805
_808.466  202.957
   relation checkEqOfMatrices matrices
1
   run=: 3 : 0
shape=.1+?2#100
m=.(genUniformMatrix shape),(genUniformMatrix shape),:(genUniformMatrix shape)
relation checkEqOfMatrices m
)
   (+/)(run"0)100#0
100
```
- <img src="https://latex.codecogs.com/svg.image?(s_1&space;&plus;&space;s_2)A&space;=&space;s_1&space;A&space;&plus;&space;s_2&space;A" title="(s_1 + s_2)A = s_1 A + s_2 A" />
```
   leftR=: 4 : '( (0{x) + (1{x) ) * y'
   rightR=: 4 : '( (0{x) * y ) + ( (1{x) * y )'
   relation=: leftR`rightR
   ]scalars=: _100 100 runiform 2
61.2493 16.9455
   ]matrix=: genUniformMatrix 4 2
 472.662  94.5818
 _267.82 _744.562
 181.815 _345.643
_141.425  45.4723
   checkEqOfMatricesWithScalars=: 4 : '(x relation@.0 y) -: (x relation@.1 y)'
   scalars checkEqOfMatricesWithScalars matrix
1
      run=: 3 : 0
shape=.1+?2#100
m=.genUniformMatrix shape
s=. _100 100 runiform 2
s checkEqOfMatricesWithScalars m
)
   (+/)(run"0)100#0
100
```
- <img src="https://latex.codecogs.com/svg.image?s_1&space;(A&space;&plus;&space;B)&space;=&space;s_1&space;A&space;&plus;&space;s_1&space;B" title="s_1 (A + B) = s_1 A + s_1 B" />
```
   ]scalars=: _100 100 runiform 1
24.4943
   ]matrices=: (genUniformMatrix 4 2),:(genUniformMatrix 4 2)
 853.433 _400.652
_164.792  375.372
 675.547  584.175
_69.8546  211.161

_350.586   815.65
_873.687  _226.76
_322.805 _808.466
 202.957 _428.685
   data=:scalars;matrices
   >0{data
24.4943
   >1{data
 853.433 _400.652
_164.792  375.372
 675.547  584.175
_69.8546  211.161

_350.586   815.65
_873.687  _226.76
_322.805 _808.466
 202.957 _428.685

   leftR=: 4 : '(0{x) * ( (0{y) + (1{y) )'
   rightR=: 4 : '( (0{x) * (0{y) ) + ( (0{x) * (1{y) )'
   relation=: leftR`rightR

   NB. it is important here to set tolerance to the least strict value - floating point arithmetic
   9!:18 ''
1e_11

   checkEqOfMatricesScalarsRel=: 4 : '( (>0{y) x@.0 (>1{y) ) -: ( (>0{y) x@.1 (>1{y) )'
   relation checkEqOfMatricesScalarsRel data
1

   run=: 3 : 0
shape=.1+?2#100
m=.(genUniformMatrix shape),:(genUniformMatrix shape)
s=. _100 100 runiform 1
data=.s;m
relation checkEqOfMatricesScalarsRel data
)
   (+/)(run"0)100#0
100
```
- <img src="https://latex.codecogs.com/svg.image?s_1&space;(s_2&space;A)&space;=&space;(s_1&space;s_2)A" title="s_1 (s_2 A) = (s_1 s_2)A" />
```
   leftR=: 4 : '(0{x) * ( (1{x) * (0{y) )'
   rightR=: 4 : '( (0{x) * (1{x) ) * (0{y)'
   relation=: leftR`rightR
   run=: 3 : 0
shape=.1+?2#100
m=.genUniformMatrix shape
s=. _100 100 runiform 2
data=.s;m
relation checkEqOfMatricesScalarsRel data
)
   (+/)(run"0)100#0
100
```
- <img src="https://latex.codecogs.com/svg.image?A&space;&plus;&space;(-1)A&space;=&space;0&space;" title="A + (-1)A = 0 " />
```
   leftR=: 4 : '(0{y) + ( (0{x) * (0{y) )'
   rightR=: 4 : '1{y'
   relation=: leftR`rightR
   run=: 3 : 0
shape=.1+?2#100
m=.(genUniformMatrix shape),:(shape $ 0)
s=. _1
data=.s;m
relation checkEqOfMatricesScalarsRel data
)
   (+/)(run"0)1000#0
1000
```

### Solution to exercise 19
```
    m=: 2 2 $ 0 1 1 1
    mp=: +/ .*
    f=: 3 : 'y mp m'
    fibb=: 3 : '(<(<0),(<1)) { (f^:y m)'
    fibb"0 i.30
1 1 2 3 5 8 13 21 34 55 89 144 233 377 610 987 1597 2584 4181 6765 10946 17711 28657 46368 75025 121393 196418 317811 514229 832040
    NB. this is very efficient implementation
```
