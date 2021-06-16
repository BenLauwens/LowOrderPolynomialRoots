using Polynomials: fromroots, Polynomial
import Polynomials
using LowOrderPolynomialRoots: roots
using BenchmarkTools


p = fromroots([1.0])
@show p roots(p);

p = fromroots([1.0, 1.0])
@show p roots(p);

p = Polynomial([1.0, 0.0, 1.0])
@show p roots(p);

p = Polynomial([-1.0, 1.0, -1.0, 1.0])
@show p roots(p);

p = fromroots([1.0, 1.0, 1.0])
@show p roots(p);

p = fromroots([1.0, 1.0, -1.0])
@show p roots(p);

p = fromroots([2.0, 1.0, -1.0])
@show p roots(p);

p = Polynomial([-2.0e9, 3.0e9, -1.0e9, 1.0])
@show p roots(p);

p = Polynomial([1.0, 1.0, 1.0, 1.0])
@show p roots(p);

p = fromroots([1.0, 1.0, 1.0, 1.0])
@show p roots(p);

p = fromroots([1.0, -1.0, 2.0, 1.0])
@show p roots(p);

p = Polynomial([-3.0, -10.0, -5.0, 6.0, 1.0])
@show p roots(p);

p = Polynomial([-12.0, -7.0, 7.0, 6.0, 1.0])
@show p roots(p);

p = Polynomial([4.0, 0.0, 4.0, 0.0, 1.0])
@show p roots(p);

p = Polynomial([-1.0, 0.0, 0.0, 0.0, 1.0])
@show p roots(p);

p = Polynomial([24.0, -50.0, 35.0, -10.0, 1.0])
@show p roots(p);

p = fromroots([2.35, 2.36, 2.37])
@show p roots(p);

p = fromroots([2.35, 2.35, 2.35, 2.56])
@show p roots(p);

p = fromroots([1.0, -1.0, 2.0, 1.0, 3.0])
@show p roots(p);

a = 1.0e3
p = fromroots([a, a, 1/a])
@show p roots(p);

a = 1.0e6
p = fromroots([a, a, 1/a])
@show p roots(p);

a = 1.0e9
p = fromroots([a, a, 1/a])
@show p roots(p);

a = 1.0e10
p = fromroots([a, -a, 1])
@show p roots(p);

a = 1.0e-10
p = fromroots([a, -a, 1])
@show p roots(p);

@benchmark filter(x->isreal(x) && x >= 0, roots(p)) setup=(p=Polynomial(vec(rand(1, 3)) .- 0.5))
@benchmark filter(x->isreal(x) && real(x) >= 0, Polynomials.roots(p)) setup=(p=Polynomial(vec(rand(1, 3)) .- 0.5))

@benchmark filter(x->isreal(x) && x >= 0, roots(p)) setup=(p=Polynomial(vec(rand(1, 4)) .- 0.5))
@benchmark filter(x->isreal(x) && real(x) >= 0, Polynomials.roots(p)) setup=(p=Polynomial(vec(rand(1, 4)) .- 0.5))

@benchmark filter(x->isreal(x) && real(x) >= 0, roots(p)) setup=(p=Polynomial(vec(rand(1, 5)) .- 0.5))
@benchmark filter(x->isreal(x) && real(x) >= 0, Polynomials.roots(p)) setup=(p=Polynomial(vec(rand(1, 5)) .- 0.5))

@benchmark filter(x->isreal(x) && x >= 0, roots(p)) setup=(p=fromroots([0.1, 0.2, 0.3, 0.4]))
@benchmark filter(x->isreal(x) && real(x) >= 0, Polynomials.roots(p)) setup=(p=fromroots([0.1, 0.2, 0.3, 0.4]))