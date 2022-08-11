# Intermediate linear algebra

## Contents
1. [A partitioned matrix](#partitioned-matrix) - IN PROGRESS
2. [SVD](#svd)
3. [Jordan decomposition](#jordan-decomposition)
4. [Kronecker product](#kronecker-product)
5. [Moore-Penrose inverse](#moore-penrose-inverse)
6. [Adjoint matrix](#adjoint-matrix)
7. [Hadamard product](#hadamard-product)
8. [Commutation matrix](#commutation-matrix)
9. [Duplication matrix](#duplication-matrix)
10. [Gramian matrix](#gramian-matrix)
11. [FFT](#fft)

[Project3](#project3)

[Solutions to exercices](#intermediate-linear-algebra-solutions-to-exercises)


## Partitioned matrix

A partitioned matrix is obtained when an underlying matrix is carved out both row-wise and column-wise in such a way that inner blocks are created.
Those blocks cover fully the underlying matrix and do not mutually overlap. Moreover, we usually restrict the partitioning in such a way that
each block row has the same number of rows within the same block row, and each block column has the same number of columns within the same block column.
At first we will develop a technique that partitions a given underlying matrix into the same dimensional blocks. We will use `x(u;._3)y` that
applies verb `u` to each tile of a regular tiling of `y` specified by `x`. Basically we will define
```j
   NB. We will define tile as 2x2 which is translated throughout the underlying matrix by vector (2,2)
   NB. So we will not have overlap here upon tile translation
   translationVector=: 2 2
   tileSize=: 2 2
   ] tile=: translationVector ,: tileSize
2 2
2 2

   NB. boxing as verb
   u=: <
   NB. discarding uncompleted tiles
   tile < ;._3 input
   NB. including uncompleted tiles
   tile < ;.3 input

   toR=:x:
   genUniformMatrix=: 3 : 'y $ toR <. ( _10 10 runiform ((0{y) * (1{y)))'
   ]m=:genUniformMatrix 10 10
_6  _3  _2  5 _1  5  _2 _4  _3  _8
 6   6 _10 _8  9 _6 _10  4  _5  _6
 3   9   4  2  9  8   5  5  _1  _1
 4  _6   7  8  9  6  _1  1  _7   0
_7 _10  _9  7  9 _6   0 _6  _9  _6
_5  _9  _3  5 _9  9   7  7  _7   3
 6  _1  _8 _4  8 _3  _6  8   2 _10
_9   6  _8  4 _2 _7  _6 _9   2   9
 4  _4   1 _7 _3  3  _8  0 _10   0
_2  _2   1  1 _9 _8   9 _2   4   3
   ]mblock=: (2 2,: 2 2) <;._3 m
┌──────┬──────┬─────┬──────┬─────┐
│_6 _3 │ _2  5│_1  5│ _2 _4│_3 _8│
│ 6  6 │_10 _8│ 9 _6│_10  4│_5 _6│
├──────┼──────┼─────┼──────┼─────┤
│3  9  │4 2   │9 8  │ 5 5  │_1 _1│
│4 _6  │7 8   │9 6  │_1 1  │_7  0│
├──────┼──────┼─────┼──────┼─────┤
│_7 _10│_9 7  │ 9 _6│0 _6  │_9 _6│
│_5  _9│_3 5  │_9  9│7  7  │_7  3│
├──────┼──────┼─────┼──────┼─────┤
│ 6 _1 │_8 _4 │ 8 _3│_6  8 │2 _10│
│_9  6 │_8  4 │_2 _7│_6 _9 │2   9│
├──────┼──────┼─────┼──────┼─────┤
│ 4 _4 │1 _7  │_3  3│_8  0 │_10 0│
│_2 _2 │1  1  │_9 _8│ 9 _2 │  4 3│
└──────┴──────┴─────┴──────┴─────┘
  NB. selecting operators applies here
   (<(<1),(<1)) { mblock
┌───┐
│4 2│
│7 8│
└───┘
   (<(<1,2),(<1,2,3)) { mblock
┌────┬─────┬────┐
│4 2 │9 8  │ 5 5│
│7 8 │9 6  │_1 1│
├────┼─────┼────┤
│_9 7│ 9 _6│0 _6│
│_3 5│_9  9│7  7│
└────┴─────┴────┘
   ]mblock=: (5 5,: 5 5) <;._3 m
┌────────────────┬────────────────┐
│_6  _3  _2  5 _1│ 5  _2 _4 _3 _8 │
│ 6   6 _10 _8  9│_6 _10  4 _5 _6 │
│ 3   9   4  2  9│ 8   5  5 _1 _1 │
│ 4  _6   7  8  9│ 6  _1  1 _7  0 │
│_7 _10  _9  7  9│_6   0 _6 _9 _6 │
├────────────────┼────────────────┤
│_5 _9 _3  5 _9  │ 9  7  7  _7   3│
│ 6 _1 _8 _4  8  │_3 _6  8   2 _10│
│_9  6 _8  4 _2  │_7 _6 _9   2   9│
│ 4 _4  1 _7 _3  │ 3 _8  0 _10   0│
│_2 _2  1  1 _9  │_8  9 _2   4   3│
└────────────────┴────────────────┘
   ]mblock=: (2 3,: 2 3) <;._3 m
┌─────────┬────────┬─────────┐
│_6 _3  _2│ 5 _1  5│ _2 _4 _3│
│ 6  6 _10│_8  9 _6│_10  4 _5│
├─────────┼────────┼─────────┤
│3  9 4   │2 9 8   │ 5 5 _1  │
│4 _6 7   │8 9 6   │_1 1 _7  │
├─────────┼────────┼─────────┤
│_7 _10 _9│7  9 _6 │0 _6 _9  │
│_5  _9 _3│5 _9  9 │7  7 _7  │
├─────────┼────────┼─────────┤
│ 6 _1 _8 │_4  8 _3│_6  8 2  │
│_9  6 _8 │ 4 _2 _7│_6 _9 2  │
├─────────┼────────┼─────────┤
│ 4 _4 1  │_7 _3  3│_8  0 _10│
│_2 _2 1  │ 1 _9 _8│ 9 _2   4│
└─────────┴────────┴─────────┘
   ]mblock=: (5 2,: 5 2) <;._3 m
┌──────┬──────┬─────┬──────┬───────┐
│_6  _3│ _2  5│_1  5│ _2 _4│_3 _8  │
│ 6   6│_10 _8│ 9 _6│_10  4│_5 _6  │
│ 3   9│  4  2│ 9  8│  5  5│_1 _1  │
│ 4  _6│  7  8│ 9  6│ _1  1│_7  0  │
│_7 _10│ _9  7│ 9 _6│  0 _6│_9 _6  │
├──────┼──────┼─────┼──────┼───────┤
│_5 _9 │_3  5 │_9  9│ 7  7 │ _7   3│
│ 6 _1 │_8 _4 │ 8 _3│_6  8 │  2 _10│
│_9  6 │_8  4 │_2 _7│_6 _9 │  2   9│
│ 4 _4 │ 1 _7 │_3  3│_8  0 │_10   0│
│_2 _2 │ 1  1 │_9 _8│ 9 _2 │  4   3│
└──────┴──────┴─────┴──────┴───────┘

   ]m=:genUniformMatrix 10 10
 2 _4  8 _9 _3 _4 _9  2 _5  8
_5 _2  3  6  5 _1  2 _7  4 _7
 3  8 _2  2  7  7  4  9 _1  8
 8  6  2 _8 _6 _6 _5  3  6  4
 1 _1  7 _2 _6  4  9 _8 _1 _7
_1 _3 _4  5  4 _9  3  4  3  3
 9  3 _9 _6  6  1  2  9 _6 _8
_9 _7  8  0  5 _4 _5  0  4 _1
_7 _2  9  4 _5  3 _2 _9  9 _8
 0  3 _7 _8  7  8  2  1 _1  6
   (2 2,: 2 2) <;._3 m
┌─────┬─────┬─────┬─────┬─────┐
│ 2 _4│8 _9 │_3 _4│_9  2│_5  8│
│_5 _2│3  6 │ 5 _1│ 2 _7│ 4 _7│
├─────┼─────┼─────┼─────┼─────┤
│3 8  │_2  2│ 7  7│ 4 9 │_1 8 │
│8 6  │ 2 _8│_6 _6│_5 3 │ 6 4 │
├─────┼─────┼─────┼─────┼─────┤
│ 1 _1│ 7 _2│_6  4│9 _8 │_1 _7│
│_1 _3│_4  5│ 4 _9│3  4 │ 3  3│
├─────┼─────┼─────┼─────┼─────┤
│ 9  3│_9 _6│6  1 │ 2 9 │_6 _8│
│_9 _7│ 8  0│5 _4 │_5 0 │ 4 _1│
├─────┼─────┼─────┼─────┼─────┤
│_7 _2│ 9  4│_5 3 │_2 _9│ 9 _8│
│ 0  3│_7 _8│ 7 8 │ 2  1│_1  6│
└─────┴─────┴─────┴─────┴─────┘
   (2 2,: 2 2) <@mean@,;._3 m
┌────┬────┬────┬────┬─────┐
│_9r4│2   │_3r4│_3  │0    │
├────┼────┼────┼────┼─────┤
│25r4│_3r2│1r2 │11r4│17r4 │
├────┼────┼────┼────┼─────┤
│_1  │3r2 │_7r4│2   │_1r2 │
├────┼────┼────┼────┼─────┤
│_1  │_7r4│2   │3r2 │_11r4│
├────┼────┼────┼────┼─────┤
│_3r2│_1r2│13r4│_2  │3r2  │
└────┴────┴────┴────┴─────┘
```
Let's assume now that we want to adopt variable tile sizes that satisfy:
- tiles do not overlap and cover the underlying matrix completely
- each block row has the same number of rows within the same block row
- each block column has the same number of columns within the same block column

It is achievable using `x(u;.1)y` when x marks both rows (column) block starting boundaries with 1s and their height (lenght) using 0s.
There is also an option of specifying the ending of partitions with 1s - `x(u;.2)y` when defining blocks in the underlying matrix.
```
   ]m=:genUniformMatrix 10 10
 2 _4  8 _9 _3 _4 _9  2 _5  8
_5 _2  3  6  5 _1  2 _7  4 _7
 3  8 _2  2  7  7  4  9 _1  8
 8  6  2 _8 _6 _6 _5  3  6  4
 1 _1  7 _2 _6  4  9 _8 _1 _7
_1 _3 _4  5  4 _9  3  4  3  3
 9  3 _9 _6  6  1  2  9 _6 _8
_9 _7  8  0  5 _4 _5  0  4 _1
_7 _2  9  4 _5  3 _2 _9  9 _8
 0  3 _7 _8  7  8  2  1 _1  6
   (1 0 1 0 0 0 1 1 0 0; 1 0 0 1 0 0 1 0 0 0) <;.1 m
┌────────┬────────┬───────────┐
│ 2 _4 8 │_9 _3 _4│_9  2 _5  8│
│_5 _2 3 │ 6  5 _1│ 2 _7  4 _7│
├────────┼────────┼───────────┤
│ 3  8 _2│ 2  7  7│ 4  9 _1  8│
│ 8  6  2│_8 _6 _6│_5  3  6  4│
│ 1 _1  7│_2 _6  4│ 9 _8 _1 _7│
│_1 _3 _4│ 5  4 _9│ 3  4  3  3│
├────────┼────────┼───────────┤
│9 3 _9  │_6 6 1  │2 9 _6 _8  │
├────────┼────────┼───────────┤
│_9 _7  8│ 0  5 _4│_5  0  4 _1│
│_7 _2  9│ 4 _5  3│_2 _9  9 _8│
│ 0  3 _7│_8  7  8│ 2  1 _1  6│
└────────┴────────┴───────────┘
   NB. It would be more user friendly to be able to specify partitioning the following (2 4 1 3; 3 3 4)
   expressNum=: 3 : 'if. y = 1 do. 1 else. (1, (y - 1) $ 0) end.'
   ;<@:expressNum"0 (2 4 1 3)
1 0 1 0 0 0 1 1 0 0
   partitionMatrix=: 4 : 0
rows=.>0{x
cols=.>1{x
((;<@:expressNum"0 rows); (;<@:expressNum"0 cols)) <;.1 y
)
   (2 4 1 3; 3 3 4) partitionMatrix m
┌────────┬────────┬───────────┐
│ 2 _4 8 │_9 _3 _4│_9  2 _5  8│
│_5 _2 3 │ 6  5 _1│ 2 _7  4 _7│
├────────┼────────┼───────────┤
│ 3  8 _2│ 2  7  7│ 4  9 _1  8│
│ 8  6  2│_8 _6 _6│_5  3  6  4│
│ 1 _1  7│_2 _6  4│ 9 _8 _1 _7│
│_1 _3 _4│ 5  4 _9│ 3  4  3  3│
├────────┼────────┼───────────┤
│9 3 _9  │_6 6 1  │2 9 _6 _8  │
├────────┼────────┼───────────┤
│_9 _7  8│ 0  5 _4│_5  0  4 _1│
│_7 _2  9│ 4 _5  3│_2 _9  9 _8│
│ 0  3 _7│_8  7  8│ 2  1 _1  6│
└────────┴────────┴───────────┘
```

## SVD

## Jordan decomposition

## Kronecker product

## Moore-Penrose inverse

## Adjoint matrix

## Hadamard product

## Commutation matrix

## Duplication matrix

## Gramian matrix

## FFT

## Intermediate linear algebra solutions to exercises
