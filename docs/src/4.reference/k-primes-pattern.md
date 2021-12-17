# K primes pattern

I was looking for a way to iterate through prime number k's.
I wanted the method that would avoid [twin primes](https://en.wikipedia.org/wiki/Twin_prime).
Twin primes are only 2 apart, e.g.:
(3, 5), (5, 7), (11, 13), (17, 19), (29, 31), (41, 43), (59, 61), (71, 73), (101, 103), (107, 109).
Re-assessing data using such closely related window sizes seems like it would be for minimal gain.
With the same principle of trying to avoid reprocessing data for minimal gain, an ideal iteration would steadily increase in distance between the current k-length and the next.

```julia
julia> step = 2
2

julia> current_prime = 11
11

julia> while Primes.isprime(current_prime)
           current_prime = current_prime + step
           step += 2
           @show current_prime
       end
current_prime = 13
current_prime = 17
current_prime = 23
current_prime = 31
current_prime = 41
current_prime = 53
current_prime = 67
current_prime = 83
current_prime = 101
current_prime = 121
```