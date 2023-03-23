name ends in \_lpdf
when you use them in a distributional assignment
you drop the suffix
and move the first argument to the left of the ~ sign

similar for \_lpmf

you can overload
Stan resolves by looking at the arguments
there is some ambiguity because of promotions

it allows for recursive functions

how to generate truncated distributions by hand

assignment statements: alternatice to distributional statements
allow for custom probability statements, by incrementing the total log prob

oh, so actually, the distributional statements are syntatic sugar
they are shorthand for augmenting the log-probability
but we augment by unnormalized densities to save computation

user defined distros
can you define custom unnormalized distributions

## Examples

vectorized statements
if you don't specify a prior, it will be an improper one
