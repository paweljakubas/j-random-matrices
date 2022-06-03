# Linear algebra and random matrices in J

by Pawel Jakubas, PhD

These are notes that cover a number of topics from linear algebra that I have
found fundamental to master random matrices using J language. The prerequisite
for fully comprehending the examples below is *Learning J* which is the recommended first
introductory material when learning J. The great example how beautifully and effciently J can be used in a specific domain is
a wonderful *Fractals, Visualization and J*. Besides that a list of high quality book references is specified.
This tutorial is supposed to be very practical and, I hope, has an educational flavor. I do not enforce a tacit representation as I think,
although undoubtly very elegant, can be difficult to appreciate for a beginner J enthusiasta and could steepen his/her learning curve.
I also omit discussing the performance optimisations that could be adopted, once again for the sake of clarity and for being user-friendly
as much as possible. The aim of this tutorial is to document my adventures and experiments with J that I found interesting when exploring a number of topics
culminating in random matrices, which although used in physics for decades with many successes, only relatively recently have started to play
an important role in other areas of science like biology, computational statistics or machine learning.
The aim of these notes is to document my efforts for the future myself, but there is a chance they could also be useful for someone else.
And I mean here anyone quantitative and interested to learn about J and witness the potential of this powerful array programming language.
For J experts it will be rather very far from what they would write down, maybe not interesting at all,
but for beginning J programmers, like myself, presumably the proof that this array language, as in the current form, is very
useful, powerful and intellectually attractive. I attach a number of examples, along with the solutions,
as a documentation of my hacking, also for anyone that would be eager to learn through doing. I strongly believe
in correctness of Richard Feynman's famous quote:
*It's the way I study - to understand something by trying to work it out or, in other words, to understand something by creating it.*

The sequels are planned and will cover the author's J experiments in:
- probability and statistics
- data mining
- quantum mechanics
- network science and graph spectra
- statistical physics

I find this hacking so enjoyable that I hope this will become a multi-year project for me.
I bet that investing in J will prove very fruitful as I will not only get extremely powerful calculator at my disposal,
but also quantitative toolbox massively outperfoming other popular approaches. For sure, when it comes to
quantitative prototyping.

Following is the list of books that I found especially useful during the first installement:
### Linear algebra, optimization and machine learning
1. Linear Algebra and Optimization for Machine Learning. A Textbook, Charu C. Aggarwal, Springer 2020
2. Matrix Differential Calculus with Applications in Statistics and Econometrics, 3rd Edition, Jan R. Magnus, Heinz Neudecker, Wiley 2019
3. Scalar, Vector, and Matrix Mathematics: Theory, Facts, and Formulas - Revised and Expanded Edition, Dennis S. Bernstein, Princeton University Press 2018
4. Matrix Computations, 4th Edition, Gene H. Golub and Charles F. Van Loan, John Hopkins University Press
### Random matrices
5. A First Course in Random Matrix Theory for Physicists, Engineers and Data Scientists, Marc Potters, Jean-Philippe Bouchaud, Cambridge University Press 2021
### J language
6. Learning J. An Introduction to the J Programming Language, Roger Stokes, [https://www.jsoftware.com/help/learning/contents.htm#toc]
7. Fractals, Visualization & J, 4th ed. (2 Parts), Clifford A. Reiter 2016
8. Fifty Shades of J, Norman Thomson, [https://code.jsoftware.com/wiki/Fifty_Shades_of_J]

   *Throughout the code J903 version of J lang was used.*


## Contents
### [Basic linear algebra](chapters/basic-algebra.md)
### [Intermediate linear algebra](chapters/intermediate-algebra.md)
### [Advanced linear algebra I]
### [Advanced linear algebra II]
### [Random matrix theory]
### [Random matrix theory applications]
