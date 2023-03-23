# Stan User guide

## arrays, vectors and indexes

matrices just store the data and the dimensions
in contrast with arrays, that store arrays of arrays

you need matrices to do linear algebra

matrices store in column major order, so you should loop over them like that

matrices store contiguously

if you want to do stuff by rows or columns, do arrays
if you want to do linear algebra, do matrices

no automatic conversion between vectors, arrays or matrices
but plenty of functions to do it

expressions are evaluated before assignment happens
so there is no risk of aliasing

versatile implementation of multiple indexing

slicing using the same conventions as R

there is a subtlety when you pass multiple indexes to an array of vectors
you may get a vector or an array of reals
the thing to know is that when you pass a single index
you reduce the corresponding dimension

when you slice a matrix, there is an automatic conversion to vectors

same rules for arrays

## functions

defined in their own block, which has to be at the top
used in the `generated quantities` block

reject is like a throw error
but it's a little tricky because they affect the running in different ways
depending on where the function is called

arguments of functions don't require size
but also cannot restrict by putting bounds or anything like that
however for arrays you do need to specify dimension, with proper syntax

`data` qualifier for arguments

you can have `void` functions
they still need a `return` but without argument
they can be used as statements, rather than assigned to something

functions ending in `_lp` can use distribution statements

it's different to try to define the distribution of parameters
than to generate random numbers
for the latter, you can use predefined prgn
but your function has to end in `_rng`
