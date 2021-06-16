### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ d12cd828-5a8c-11eb-09b0-6b8f09352ec4
using PlutoUI

# ╔═╡ 46716430-5a8f-11eb-0840-97b849c24784
using MacroTools: flatten

# ╔═╡ b859cf6c-5a92-11eb-0f1d-fffd383f3bc9
using BenchmarkTools

# ╔═╡ 95a04d4e-5c37-11eb-3dbd-5d60a3ab9da8
using Random

# ╔═╡ 4cb69cda-5c49-11eb-11f9-8d28ac5252bf
using Polynomials

# ╔═╡ 19a96fd2-bbe7-4a12-b83e-4745307bd756
using PolynomialRoots

# ╔═╡ ca19a1bb-9721-40c5-8806-53b2c6fe70ae
using StaticArrays

# ╔═╡ a9d53ab8-5a8c-11eb-2606-4f346de6853d
md"# Interval Newton method"

# ╔═╡ 16aa3ff2-5dc7-11eb-3034-fdc0827b9088
@inline function quadraticsmallestpositiveroot(a::F, b::F, c::F) where F <: AbstractFloat
	_a = one(F) / a
	b, c = -0.5b * _a, c * _a
	Δ = muladd(b, b, -c) # b * b - c
	if Δ < -4eps(F) # Complex roots
		typemax(F) 
	elseif Δ > 4eps(F) # Real roots
		if b > eps(F)
			c > eps(F) ? c / (b + sqrt(Δ)) : b + sqrt(Δ)
		elseif b < -eps(F)
			c < eps(F) ? c / (b - sqrt(Δ)) : typemax(F)
		else
			sqrt(-c)
		end
	else # Double real root
		b > -eps(F) ? b : typemax(F) 
	end
 end

# ╔═╡ bf557e01-2f78-4715-b13f-9113fc68b62a
with_terminal() do
	for _ in 1:1000000
		a = vec(2(rand(1,3) .- 0.5))
		pol = Polynomials.roots(Polynomial(a))
		p = filter(r->isreal(r) && real(r) > -eps(), pol)
		n = quadraticsmallestpositiveroot(a[end:-1:1]...)
		if length(p) !== 0 && n !== nothing 
			if abs(minimum(real.(p)) - n) > 1.0e-7
				@show a p n pol
				break
			end
		end
	end
end

# ╔═╡ 5811c086-5dc8-11eb-3cf8-df6b37715bb3
@inline function cubicsmallestpositiveroot(a::F, b::F, c::F, d::F) where F <: AbstractFloat
    _a = one(F) / a
    b, c, d = b * _a, c * _a, d * _a
    m = b < c ? b : c
    m = d < m ? d : m
    m > eps(F) && return typemax(F) # Cauchy bound
	_3 = one(F) / 3
	_9 = one(F) / 9
	SQ3 = sqrt(3one(F))
    xₙ = -b * _3
    b²_9 = b * b * _9
    yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)
    δ² = muladd(-_3, c, b²_9)
    h² = 4δ² * δ² * δ²
    Δ = muladd(yₙ, yₙ, -h²)
    if Δ > 4eps(F) # one real root and two complex roots
		p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
		q = δ² / p
		z = xₙ + p + q
		z > -eps(F) ? z : typemax(F)
    elseif Δ < -4eps(F) # three real roots
		θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
		δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
		z₁ = 2δ * cos(θ)
		z₂ = muladd(-0.5, z₁, xₙ)
		z₃ = SQ3 * δ * sin(θ)
		x₁ = xₙ + z₁
		x₂ = z₂ + z₃
		x₃ = z₂ - z₃
		x = x₁ > -eps(F) ? x₁ : typemax(F)
		x = x₂ > -eps(F) && x₂ < x ? x₂ : x
		x₃ > -eps(F) && x₃ < x ? x₃ : x
    else # double or triple real roots
		δ = cbrt(0.5yₙ)
		x₁ = xₙ + δ
		x₂ = xₙ - 2δ
		x = x₁ > -eps(F) ? x₁ : typemax(F)
		x₂ > -eps(F) && x₂ < x ? x₂ : x
    end
 end

# ╔═╡ 20078858-600a-11eb-0953-4569a4a1eb8a
with_terminal() do
	for _ in 1:1000000
		a = vec(2(rand(1,4) .- 0.5))
		pol = Polynomials.roots(Polynomial(a))
		p = filter(r->isreal(r) && real(r) > -eps(), pol)
		n = cubicsmallestpositiveroot(a[end:-1:1]...)
		if length(p) !== 0 && n !== nothing 
			if abs(minimum(real.(p)) - n) > 1.0e-7
				@show a p n pol
				#break
			end
		end
	end
end

# ╔═╡ aa2e9fd4-5dcc-11eb-2fe2-417d13773cf6
begin
	@inline function quadraticrealroots(b::F, c::F) where F <: AbstractFloat
		b, c = -0.5b, c
		Δ = muladd(b, b, -c) #b * b - c
		if Δ < -4eps(F) # Complex roots
			tuple()
		elseif Δ > 4eps(F) # Real roots
			q = b > 0 ? b + sqrt(Δ) : b - sqrt(Δ)
			(q, c / q)
		else # Double real root
			tuple(b)
		end
	end

	@inline function cubicmaxroot(b::F, c::F, d::F) where F <: AbstractFloat
		_3 = one(F) / 3
		_9 = one(F) / 9
		SQ3 = sqrt(3one(F))
		xₙ = -b * _3
		b²_9 = b * b * _9
		yₙ = muladd(muladd(-2, b²_9, c), xₙ, d) #d + xₙ * (c - 2b²_9)
		δ² = muladd(-_3, c, b²_9) #b²_9 - c *_3
		h² = 4δ² * δ² * δ²
		Δ = muladd(yₙ, yₙ, -h²) #yₙ * yₙ - h²
		if Δ > 4eps(F) # one real root and two complex roots
			p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
			q = δ² / p # cbrt(0.5 * (-yₙ - √Δ)) : cbrt(0.5 * (-yₙ + √Δ))
			xₙ + p + q
		elseif Δ < -4eps(F) # three real roots
			θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
			δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
			z₁ = 2δ * cos(θ)
			z₂ = muladd(-0.5, z₁, xₙ) #xₙ - 0.5z₁   # x₂ = xₙ + 2δ * cos(2π / 3 - θ)
			z₃ = SQ3 * δ * sin(θ) # x₃ = xₙ + 2δ * cos(4π / 3 - θ)
			x₁ = xₙ + z₁
			x₂ = z₂ + z₃
			x₃ = z₂ - z₃
			x = x₁ > x₂ ? x₁ : x₂
			x₃ > x ? x₃ : x
		else # double or triple real roots
			δ = cbrt(0.5yₙ)
			x₁ = xₙ + δ
			x₂ = xₙ - 2δ
			x₁ > x₂ ? x₁ : x₂
		end
	end

	@inline function quarticsmallestpositiveroot(a::F, b::F, c::F, d::F, e::F) where F <: AbstractFloat
		_a = inv(a)
		b, c, d, e = b * _a, c * _a, d * _a, e * _a
		m = b < c ? b : c
		m = d < m ? d : m
		m = e < m ? e : m
		m > eps(F) && return typemax(F) # Cauchy bound
		xₙ = -0.25b
		b² = b * b
		p = muladd(-0.375, b², c)
		q = muladd(muladd(-0.5, b², 2c), xₙ, d)
		r = muladd(muladd(muladd(-0.1875, b², c), xₙ, d), xₙ, e)
		if abs(r) < 4eps(F)
			xₙ
		elseif abs(q) < 4eps(F)
			y² = quadraticrealroots(p, r)
			l = length(y²)
			l === 0 && return typemax(F)
			y₁ = y²[1] < -eps(F) ? typemax(F) : √y²[1]
			z₁ = xₙ + y₁
			z₂ = xₙ - y₁
			x = z₁ > -eps(F) ? z₁ : typemax(F)
			x = z₂ > -eps(F) && z₂ < x ? z₂ : x
			l === 1 && return x
			y₂ = y²[2] < -eps(F) ? typemax(F) : √y²[2]
			z₃ = xₙ + y₂
			z₄ = xₙ - y₂
			x = z₃ > -eps(F) && z₃ < x ? z₃ : x
			z₄ > -eps(F) && z₄ < x ? z₄ : x
		else
			h² = cubicmaxroot(2p, muladd(p, p, -4r), -q * q)
			h = √h²
			j = 0.5(p + h² - q / h)
			y = quadraticrealroots(h, j)
			l = length(y)
			x = if l === 0
				typemax(F)
			elseif l === 1
				z₁ = xₙ + y[1]
				z₁ > -eps(F) ? z₁ : typemax(F)
			else
				z₁ = xₙ + y[1]
				z₁ = z₁ > -eps(F) ? z₁ : typemax(F)
				z₂ = xₙ + y[2]
				z₂ > -eps(F) && z₂ < z₁ ? z₂ : z₁
			end
			y = quadraticrealroots(-h, r / j)
			l = length(y)
			if l === 0
				x
			elseif l === 1
				z₁ = xₙ + y[1]
				#z₁ > -eps(F) ? z₁ : x
				z₁ > -eps(F) && z₁ < x ? z₁ : x
			else
				z₁ = xₙ + y[1]
				#x = z₁ > -eps(F) ? z₁ : x
				x = z₁ > -eps(F) && z₁ < x ? z₁ : x
				z₂ = xₙ + y[2]
				z₂ > -eps(F) && z₂ < x ? z₂ : x
			end
		end
	end
end

# ╔═╡ a9381b20-6009-11eb-35ca-9fc500eb33d8
with_terminal() do
	for _ in 1:1000000
		a = vec(2(rand(1,5) .- 0.5))
		p = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		n = quarticsmallestpositiveroot(a[end:-1:1]...)
		if length(p) !== 0 && n !== nothing 
			if abs(minimum(real.(p)) - n) > 1.0e-3
				@show a p n
				break
			end
		end
	end
end

# ╔═╡ 53ec20a4-5a95-11eb-0c42-3335e92fa572
@generated function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	ex = quote m = arg1 < arg2 ? arg1 : arg2 end
	for i in 1:length(args)
		ex = quote
			$ex
			@inbounds m = args[$i] < m ? args[$i] : m
		end
	end
	ex |> flatten
end

# ╔═╡ 7a3b1d28-5a95-11eb-1651-f92916f5435e
with_terminal() do
	@code_warntype smallest(1.0, 2.0, -3.0)
end

# ╔═╡ c9a1943e-5c3a-11eb-18be-6fdf094d0329
with_terminal() do
	@btime smallest(1.0, 2.0, -3.0)
end

# ╔═╡ 112bcf98-5a91-11eb-32c8-65076d5ddb87
@generated function horner(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
 	l = length(coeffs)
	l === 0 && return quote coeff1 end
	ex = :(coeff1)
	for i in 1:l
		ex = :(@inbounds muladd(x, $ex, coeffs[$i]))
	end
	ex |> flatten
end

# ╔═╡ 439e5aea-5a91-11eb-30e9-63c970d2b687
with_terminal() do
	@code_warntype horner(2.0, -1.0, 2.0, 3.0, 4.0)
end

# ╔═╡ a4f7e650-5c3a-11eb-0cbd-2567ba93a120
with_terminal() do
	@btime horner(2.0, -1.0, 2.0, 3.0, 4.0)
end

# ╔═╡ 1d794e94-5a8f-11eb-2403-97e74ea215d6
@generated function posintervalhorner(low::F, high::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	l === 0 && return quote coeff1, coeff1 end
	ex = quote colow, cohigh = if coeff1 > eps(F)
				@inbounds muladd(low, coeff1, coeffs[1]), muladd(high, coeff1, coeffs[1])
			elseif coeff1 < -eps(F)
				@inbounds muladd(high, coeff1, coeffs[1]), muladd(low, coeff1, coeffs[1])
			else
				@inbounds muladd(high, coeff1, coeffs[1]), muladd(high, coeff1, coeffs[1])
			end
		 end
	for i in 2:l
	ex = quote
		$ex
		colow, cohigh = if colow > eps(F)
				@inbounds muladd(low, colow, coeffs[$i]), muladd(high, cohigh, coeffs[$i])
			elseif cohigh < -eps(F)
				@inbounds muladd(high, colow, coeffs[$i]), muladd(low, cohigh, coeffs[$i])
			else
				@inbounds muladd(high, colow, coeffs[$i]), muladd(high, cohigh, coeffs[$i])
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 52548908-5a8f-11eb-3a7e-cbb7173b7d4a
with_terminal() do
	@code_warntype posintervalhorner(1.0, 2.0, 2.0, -3.0, 4.0)
end

# ╔═╡ a00080ea-5c23-11eb-3cf8-9b04bb41f04b
struct List{F<:AbstractFloat}
	prev :: Union{Nothing, List{F}}
	low :: F
	high :: F
end

# ╔═╡ 8897ddc0-5c20-11eb-1d0c-85b3daa6cc4a
@generated function gensmallestpositiverootintervalnewton(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		mm > eps(F) && return typemax(F)
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				low, high = if comid > eps(F)
					mid - comid / codom′low, mid - comid  / codom′high
				elseif comid < -eps(F)
					mid - comid  / codom′high, mid - comid / codom′low
				else
					return mid
				end
				low = domlow < low ? low : domlow
				high = domhigh > high ? high : domhigh
				#@show low high
				domlow, domhigh = if high - low > -eps(F)
					if high - low < 1e-8mid
						return 0.5(low + high)
					end
					low, high
				elseif list === nothing
					return typemax(F)
				else
					list, low, high = list.prev, list.low, list.high
					low, high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 620642ce-5c33-11eb-345f-f5dbc336f066
with_terminal() do
	@code_warntype gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)
end

# ╔═╡ 364c6578-5c38-11eb-35a9-f7aae15de7dc
with_terminal() do
	@show gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)
end

# ╔═╡ 1efb1472-5c43-11eb-053f-4fa0a141fcb7
with_terminal() do
	for _ in 1:1
		a = vec(2(rand(1,5) .- 0.5))
		p = sort(real.(filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))))
		n = gensmallestpositiverootintervalnewton(a[end:-1:1]...)
		@show a p n
	end
end

# ╔═╡ 617dd1c2-5c48-11eb-2794-87688a920328
with_terminal() do
	Random.seed!(150)
	for _ in 1:100000
		a = vec(2(rand(1,6) .- 0.5))
		p = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		n = gensmallestpositiverootintervalnewton(a[end:-1:1]...)
		if length(p) !== 0 && n !== nothing 
			if abs(minimum(real.(p)) - n) > 1.0e-6
				@show a p n
				break
			end
		end
	end
end

# ╔═╡ a7d6ada0-5c49-11eb-2e0d-7b4a1e00da32
@generated function smallestpositiverootintervalnewtonpure(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						if abs(delta) < 1.0e-8x return x - delta end
						x -= delta
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 63da7c4a-5d01-11eb-2d97-09fa3949188b
@generated function smallestpositiverootintervalnewton(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						if abs(delta) < 1.0e-8x return x - delta end
						x -= delta
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 3d142600-5d02-11eb-3b01-efb06298021e
smallestpositiverootintervalnewton(1.0, -2.0, 1.0)

# ╔═╡ f556e97a-5dc3-11eb-0ea3-75dca8552301
with_terminal() do
	@code_warntype smallestpositiverootintervalnewton(1.0, -2.0, -1.0)
end

# ╔═╡ 200b051e-5d0b-11eb-1a48-65c65a49cfd4
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show Polynomials.roots(Polynomial(a))
	@show PolynomialRoots.roots(a)
	@show smallestpositiverootintervalnewton(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
end

# ╔═╡ 6ed0c3d8-5d07-11eb-0db8-e94277a78e89
with_terminal() do
	for _ in 1:1
		a = vec(2(rand(1,5) .- 0.5))
		p = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(p) === 0 ? Inf : minimum(real.(p)) 
		n = smallestpositiverootintervalnewton(a[end:-1:1]...)
		@show a p n
	end
end

# ╔═╡ 253115b6-5d02-11eb-25cb-29e5e60cd302
with_terminal() do
	Random.seed!(rand(1:150))
	i = 7
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol)) 
		n = gensmallestpositiverootintervalnewton(a[end:-1:1]...)
		o = if i === 3
			quadraticsmallestpositiveroot(a[end:-1:1]...)
		elseif i === 4
			cubicsmallestpositiveroot(a[end:-1:1]...)
		elseif i === 5
			p#quarticsmallestpositiveroot(a[end:-1:1]...)
		else
			p
		end
		if abs(p - n) > 1e-5p || abs(p - o) > 1e-5p
			@show a p n o pol j
			break
		end
	end
end

# ╔═╡ 452cfa54-5dca-11eb-25c8-a11fc7991956
@generated function smallestpositiverootintervalnewtonhalley(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	syms′′ = [gensym() for _ in 1:l-1]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	for i in 1:l-1
		ex = quote
			$ex
			$(syms′′[i]) = $(l-i) * $(syms′[i])
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						f′′ = horner(x, $(syms′′...))
						delta = f / f′
						delta = delta / (one(F) - 0.5delta * f′′ / f′)
						if abs(delta) < 1.0e-8x return x - delta end
						x -= delta
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 64eee2fe-5dcb-11eb-22d9-57ef47971aae
with_terminal() do
	@btime smallestpositiverootintervalnewtonhalley(1.0, -2.0, 1.0)
end

# ╔═╡ fec63166-5dcb-11eb-0943-710e42a3d0b4
with_terminal() do
	Random.seed!(rand(1:150))
	i = 5
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = smallestpositiverootintervalnewton(a[end:-1:1]...)
		o = smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
		if abs(p - o) > 1e-7 || abs(n - o) > 1e-7
			@show a pol p n o j
			break
		end
	end
end

# ╔═╡ fee1bc36-5e90-11eb-1099-51a99093e50a
@generated function smallestpositiverootintervalnewtonregulafalsi(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				if codomlow * codomhigh < eps(F)
					side = 0
					while domhigh - domlow > 0.5e-8(domhigh + domlow)
						#@show domlow domhigh codomlow codomhigh
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow)
						comid = horner(mid, $(syms...))
						if comid * codomlow < -eps(F)
							domhigh = mid
							codomhigh = comid
							if side === -1
								codomlow *= 0.5
							end
							side = -1
						elseif comid * codomhigh < -eps(F)
							domlow = mid
							codomlow = comid
							if side === 1
								codomhigh *= 0.5
							end
							side = 1
						else
							break
						end
					end
					return 0.5(domhigh + domlow)
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 1dffade6-5e9b-11eb-0b59-f752ee12b637
smallestpositiverootintervalnewtonregulafalsi(1.0, -2.0, 1.0)

# ╔═╡ 610fb1be-5e9c-11eb-209b-a3ce2af1d734
with_terminal() do
	Random.seed!(rand(1:150))
	i = 6
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = smallestpositiverootintervalnewton(a[end:-1:1]...)
		o = smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
		q = smallestpositiverootintervalnewtonregulafalsi(a[end:-1:1]...)
		if abs(p - o) > 1e-7p || abs(p - n) > 1e-8p || abs(p - q) > 1e-5p
			@show a pol p n o q j
			break
		end
	end
end

# ╔═╡ 443e4dc2-5e9f-11eb-0b9b-75858fadca7b
with_terminal() do
	@btime smallestpositiverootintervalnewtonregulafalsi(1.0, -2.0, 1.0)
end

# ╔═╡ 81730ab8-5f51-11eb-0aaa-637c3c953a7c
@generated function smallestpositiverootintervalnewtonridders(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				if codomlow * codomhigh < eps(F)
					ans = typemax(F)
					while true
						s = sqrt(comid*comid - codomlow*codomhigh)
						#if s === 0.0 return ans end
						next = mid + (codomlow > codomhigh ? 1 : -1) * (mid-domlow) * comid / s
						if abs(ans - next) < 1e-7ans return ans end
						ans = next
						conext = horner(next, $(syms...))
						if comid * conext < 0
							domlow = mid
							codomlow = comid
							domhigh = next
							codomhigh = conext
						elseif codomlow * conext < 0
							domhigh = next
							codomhigh = conext
						else
							domlow = next
							codomlow = conext
						end
						if abs(domhigh - domlow) < 1e-8ans return ans end
						mid = 0.5(domlow + domhigh)
						comid = horner(mid, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ c9876212-5f52-11eb-3d6c-cf78970daf25
smallestpositiverootintervalnewtonridders(1.0, -2.0, -1.0)

# ╔═╡ eaf11bf0-5f52-11eb-1c71-e37dd834d1a0
with_terminal() do
	Random.seed!(rand(1:150))
	i = 8
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = gensmallestpositiverootintervalnewton(a[end:-1:1]...)
		o = smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
		q = smallestpositiverootintervalnewtonridders(a[end:-1:1]...)
		if abs(p - o) > 1e-6p || abs(p - q) > 1e-6p || abs(p - n) > 1e-6p
			@show a pol p n o q j
			break
		end
	end
end

# ╔═╡ 3f5526c9-3b74-47a2-b3bf-60ed9257b5be
@generated function smallestpositiverootintervalnewtonrobust(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f*codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f*codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							x = 0.5(domlow + domhigh)
						end
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 7606e46a-f984-4f41-9cf2-f0129912d2a4
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show Polynomials.roots(Polynomial(a))
	@show smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
end

# ╔═╡ 16dc8a56-7354-4dd7-bcc1-e728c367abd0
with_terminal() do
	Random.seed!(rand(1:150))
	i = 10
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = gensmallestpositiverootintervalnewton(a[end:-1:1]...)
		q = smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
		if abs(p - q) > 1e-6p || abs(p - n) > 1e-6p
			@show a pol p n q j
			break
		end
	end
end

# ╔═╡ 4da291a7-d6d9-4078-a72d-fb50b6212504
@generated function smallestpositiverootintervalnewtonrobustregulafalsi(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f*codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f*codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = if domlow < newx < domhigh
								newx
							else
								0.5(domlow + domhigh) # bisection
							end
						end
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 09641ff0-40c3-4ce7-a932-6e65e9b656ef
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
	@show PolynomialRoots.roots(a)
	@show smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
end

# ╔═╡ 8a991fb6-67fb-4b3d-a14c-4637b19a8018
with_terminal() do
	Random.seed!(rand(1:150))
	i = 10
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = gensmallestpositiverootintervalnewton(a[end:-1:1]...)
		q = smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
		if abs(p - q) > 1e-6p || abs(p - n) > 1e-6p
			@show a pol p n q j
			break
		end
	end
end

# ╔═╡ 46820ceb-a88b-4494-9bc3-2de11ecd49e7
@generated function smallestpositiverootintervalnewtonrobustregulafalsistatic(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{$(2l), NTuple{2, F}}(undef)
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index] = tuple(domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f * codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f * codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = domlow < newx < domhigh ? newx : 0.5(domlow + domhigh) # bisection
						end
						f = horner(x, $(syms...))
					end
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ f66e1ade-056f-4fa1-be44-d16d896d009a
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
	@show PolynomialRoots.roots(a)
	@show smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsistatic(a[end:-1:1]...)
end

# ╔═╡ 4e6ca50e-a966-40c5-8ce5-5a066dcdf534
with_terminal() do
	Random.seed!(rand(1:150))
	i = 9
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		q = smallestpositiverootintervalnewtonrobustregulafalsistatic(a[end:-1:1]...)
		if abs(p - q) > 1e-6p
			@show a pol p q j
			break
		end
	end
end

# ╔═╡ d5e70a4f-a10e-4f42-be76-5d91b275737b
@generated function smallestpositiverootintervalnewtonrobustregulafalsiarray(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{40, F}(undef)
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index] = domlow < rightlow ? rightlow : domlow
						index += 1
					   	@inbounds list[index] = domhigh < righthigh ? domhigh : righthigh
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domhigh = list[index]
					index -= 1
					@inbounds domlow = list[index]
					index -= 1
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f*codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f*codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = if domlow < newx < domhigh
								newx
							else
								0.5(domlow + domhigh) # bisection
							end
						end
						f = horner(x, $(syms...))
					end
				elseif index === 0
					return typemax(F)
				else
					@inbounds domhigh = list[index]
					index -= 1
					@inbounds domlow = list[index]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ ce041037-a8c2-43f7-be16-7fe704dcb55e
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
	@show PolynomialRoots.roots(a)
	@show smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsistatic(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsiarray(a[end:-1:1]...)
end

# ╔═╡ baecb3c3-8aa4-4b2e-9c42-6460f67915f2
with_terminal() do
	Random.seed!(rand(1:150))
	i = 10
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = gensmallestpositiverootintervalnewton(a[end:-1:1]...)
		q = smallestpositiverootintervalnewtonrobustregulafalsiarray(a[end:-1:1]...)
		if abs(p - q) > 1e-6p || abs(p - n) > 1e-6p
			@show a pol p n q j
			break
		end
	end
end

# ╔═╡ d44bca78-d35e-4c1c-8fb7-9502660b0c36
@generated function gensmallestpositiverootintervalnewtonstatic(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		mm > eps(F) && return typemax(F)
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{20, NTuple{2, F}}(undef)
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index] = tuple(domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			else
				low, high = if comid > eps(F)
					mid - comid / codom′low, mid - comid  / codom′high
				elseif comid < -eps(F)
					mid - comid  / codom′high, mid - comid / codom′low
				else
					return mid
				end
				low = domlow < low ? low : domlow
				high = domhigh > high ? high : domhigh
				#@show low high
				if high - low > -eps(F)
					if high - low < 1e-8mid
						return 0.5(low + high)
					end
					domlow, domhigh = low, high
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ b29579f9-954d-4950-a9b3-c40d0f89e831
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
	@show PolynomialRoots.roots(a)
	@show smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewtonstatic(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsistatic(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsiarray(a[end:-1:1]...)
end

# ╔═╡ 71726ca0-c7ff-4292-9f97-f5dff052cead
with_terminal() do
	Random.seed!(rand(1:150))
	i = 10
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = gensmallestpositiverootintervalnewtonstatic(a[end:-1:1]...)
		q = smallestpositiverootintervalnewtonrobustregulafalsiarray(a[end:-1:1]...)
		if abs(p - q) > 1e-6p || abs(p - n) > 1e-6p
			@show a pol p n q j
			break
		end
	end
end

# ╔═╡ 6bfbc7a8-e793-4489-a135-2d5b4d95c196
@generated function smallestpositiverootintervalnewtonrobustregulafalsistaticref(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = (Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}())
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index][] = tuple(domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index][]
					index -= 1
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f * codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f * codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = domlow < newx < domhigh ? newx : 0.5(domlow + domhigh) # bisection
						end
						f = horner(x, $(syms...))
					end
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index][]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 4f8c1390-b2f4-4b7d-a98c-2ed6ad93b8e0
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
	@show PolynomialRoots.roots(a)
	@show smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewtonstatic(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsistatic(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsistaticref(a[end:-1:1]...)
end

# ╔═╡ 37dc92fc-f7de-4569-b4c3-25ce858add1d
mutable struct Tupple2FloatList{N, F<:AbstractFloat}
	data :: NTuple{N, NTuple{2, F}}
	index :: Int
	function Tupple2FloatList{N, F}() where {N, F<:AbstractFloat}
		t2fl = new{N, F}()
		t2fl.index = 0
		t2fl
	end
end

# ╔═╡ ca5fde85-5b99-4b55-a520-07d0a7e670a7
@inline function Base.push!(t2fl::Tupple2FloatList{N, F}, val::NTuple{2, F}) where {N, F<:AbstractFloat}
	i = t2fl.index + 1
	GC.@preserve t2fl unsafe_store!(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(t2fl)), convert(NTuple{2, F}, val), i)
	t2fl.index = i
	t2fl
end

# ╔═╡ 093a262a-3345-4180-b042-def2120b7fb6
@inline function Base.pop!(t2fl::Tupple2FloatList{N, F} ) where {N, F<:AbstractFloat}
	i = t2fl.index
	val = GC.@preserve t2fl unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(t2fl)), i)
	t2fl.index = i - 1
	val
end

# ╔═╡ 10f6769e-5386-476d-87eb-ac3ec97d0cf4
@inline function Base.length(t2fl::Tupple2FloatList{N, F} ) where {N, F<:AbstractFloat}
	t2fl.index
end

# ╔═╡ 198031bc-48c4-48aa-8bc2-18ff732e5b2a
@generated function smallestpositiverootintervalnewtonrobustregulafalsit2fl(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = Tupple2FloatList{$(2l), F}()
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   	push!(list, (domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh))
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif length(list) === 0
					return typemax(F)
				else
					domlow, domhigh = pop!(list)
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f * codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f * codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = domlow < newx < domhigh ? newx : 0.5(domlow + domhigh) # bisection
						end
						f = horner(x, $(syms...))
					end
				elseif length(list) === 0
					return typemax(F)
				else
					domlow, domhigh = pop!(list)
				end
			end
		end
	end
	ex |> flatten
end

# ╔═╡ 9a8f96a9-9ed5-41f2-b541-5764e6f92f98
with_terminal() do
	a = [0.749287907532135, -0.09641722621619575, 0.021665680177560098, -0.7380526416931921, -0.6881261131859229, 0.7199512046406764]
	@show filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
	@show PolynomialRoots.roots(a)
	@show smallestpositiverootintervalnewtonhalley(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobust(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewton(a[end:-1:1]...)
	@show gensmallestpositiverootintervalnewtonstatic(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsistatic(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsistaticref(a[end:-1:1]...)
	@show smallestpositiverootintervalnewtonrobustregulafalsit2fl(a[end:-1:1]...)
end

# ╔═╡ 48b330b1-ba27-43b9-a19c-2fc14de99333
with_terminal() do
	Random.seed!(rand(1:150))
	i = 10
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		n = gensmallestpositiverootintervalnewtonstatic(a[end:-1:1]...)
		q = smallestpositiverootintervalnewtonrobustregulafalsit2fl(a[end:-1:1]...)
		if abs(p - q) > 1e-6p || abs(p - n) > 1e-6p
			@show a pol p n q j
			break
		end
	end
end

# ╔═╡ ad7967ed-ce1b-4fe8-b348-ee92fab63856
with_terminal() do
	n = 4
	r = rand(1:150)
	Random.seed!(r)
	bench = @benchmark cubicsmallestpositiveroot(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewton(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonhalley(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobust(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsistatic(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsistaticref(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsit2fl(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsiarray(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark gensmallestpositiverootintervalnewton(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark gensmallestpositiverootintervalnewtonstatic(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark begin pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(p)); length(pol) === 0 ? Inf : minimum(real.(pol)) end setup=(p=Polynomial(vec(2(rand(1,$n) .- 0.5))))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark begin pol = filter(r->isreal(r) && real(r) > -eps(), PolynomialRoots.roots(p)); length(pol) === 0 ? Inf : minimum(real.(pol)) end setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
end

# ╔═╡ 118011bf-9f51-42ad-8c4a-383668d85031
with_terminal() do
	r = rand(1:150)
	n = 9
	Random.seed!(r)
	bench = @benchmark gensmallestpositiverootintervalnewton(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobust(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsistatic(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsit2fl(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsiarray(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark gensmallestpositiverootintervalnewton(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark gensmallestpositiverootintervalnewtonstatic(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark begin pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(p)); length(pol) === 0 ? Inf : minimum(real.(pol)) end setup=(p=Polynomial(vec(2(rand(1,$n) .- 0.5))))
	@show mean(bench.times) maximum(bench.times)
	Random.seed!(r)
	bench = @benchmark begin pol = filter(r->isreal(r) && real(r) > -eps(), PolynomialRoots.roots(p)); length(pol) === 0 ? Inf : minimum(real.(pol)) end setup=(p=vec(2(rand(1,$n) .- 0.5)))
	@show mean(bench.times) maximum(bench.times)
end

# ╔═╡ 42e249b8-7cce-42db-b63d-53072c85bbac
let
	Random.seed!(150)
	n = 9
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsistatic(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ d329bb22-7af2-4f00-bc70-f903d22ff33d
let
	Random.seed!(150)
	n = 9
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsit2fl(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 9e126207-9cbf-479f-bdf1-837c69aadca6
let
	Random.seed!(150)
	n = 9
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsiarray(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ fe6d4fdc-420c-4035-93b8-31846bf3f405
let
	Random.seed!(150)
	n = 4
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsistatic(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 139c4080-7f59-497a-a255-b5c802fc362c
let
	Random.seed!(150)
	n = 4
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsit2fl(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ e75f2df8-b497-4fba-96fd-0f0773a9c27f
let
	Random.seed!(150)
	n = 4
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsiarray(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ baac569e-f403-4026-b541-14f6e9af80c5
let
	Random.seed!(150)
	n = 5
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsistatic(p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ Cell order:
# ╟─a9d53ab8-5a8c-11eb-2606-4f346de6853d
# ╠═d12cd828-5a8c-11eb-09b0-6b8f09352ec4
# ╠═46716430-5a8f-11eb-0840-97b849c24784
# ╠═b859cf6c-5a92-11eb-0f1d-fffd383f3bc9
# ╠═95a04d4e-5c37-11eb-3dbd-5d60a3ab9da8
# ╠═4cb69cda-5c49-11eb-11f9-8d28ac5252bf
# ╠═19a96fd2-bbe7-4a12-b83e-4745307bd756
# ╠═ca19a1bb-9721-40c5-8806-53b2c6fe70ae
# ╠═16aa3ff2-5dc7-11eb-3034-fdc0827b9088
# ╠═bf557e01-2f78-4715-b13f-9113fc68b62a
# ╠═5811c086-5dc8-11eb-3cf8-df6b37715bb3
# ╠═20078858-600a-11eb-0953-4569a4a1eb8a
# ╠═aa2e9fd4-5dcc-11eb-2fe2-417d13773cf6
# ╠═a9381b20-6009-11eb-35ca-9fc500eb33d8
# ╠═53ec20a4-5a95-11eb-0c42-3335e92fa572
# ╠═7a3b1d28-5a95-11eb-1651-f92916f5435e
# ╠═c9a1943e-5c3a-11eb-18be-6fdf094d0329
# ╠═112bcf98-5a91-11eb-32c8-65076d5ddb87
# ╠═439e5aea-5a91-11eb-30e9-63c970d2b687
# ╠═a4f7e650-5c3a-11eb-0cbd-2567ba93a120
# ╠═1d794e94-5a8f-11eb-2403-97e74ea215d6
# ╠═52548908-5a8f-11eb-3a7e-cbb7173b7d4a
# ╠═a00080ea-5c23-11eb-3cf8-9b04bb41f04b
# ╠═8897ddc0-5c20-11eb-1d0c-85b3daa6cc4a
# ╠═620642ce-5c33-11eb-345f-f5dbc336f066
# ╠═364c6578-5c38-11eb-35a9-f7aae15de7dc
# ╠═1efb1472-5c43-11eb-053f-4fa0a141fcb7
# ╠═617dd1c2-5c48-11eb-2794-87688a920328
# ╠═a7d6ada0-5c49-11eb-2e0d-7b4a1e00da32
# ╠═63da7c4a-5d01-11eb-2d97-09fa3949188b
# ╠═3d142600-5d02-11eb-3b01-efb06298021e
# ╠═f556e97a-5dc3-11eb-0ea3-75dca8552301
# ╠═200b051e-5d0b-11eb-1a48-65c65a49cfd4
# ╠═6ed0c3d8-5d07-11eb-0db8-e94277a78e89
# ╠═253115b6-5d02-11eb-25cb-29e5e60cd302
# ╠═452cfa54-5dca-11eb-25c8-a11fc7991956
# ╠═64eee2fe-5dcb-11eb-22d9-57ef47971aae
# ╠═fec63166-5dcb-11eb-0943-710e42a3d0b4
# ╠═fee1bc36-5e90-11eb-1099-51a99093e50a
# ╠═1dffade6-5e9b-11eb-0b59-f752ee12b637
# ╠═610fb1be-5e9c-11eb-209b-a3ce2af1d734
# ╠═443e4dc2-5e9f-11eb-0b9b-75858fadca7b
# ╠═81730ab8-5f51-11eb-0aaa-637c3c953a7c
# ╠═c9876212-5f52-11eb-3d6c-cf78970daf25
# ╠═eaf11bf0-5f52-11eb-1c71-e37dd834d1a0
# ╠═3f5526c9-3b74-47a2-b3bf-60ed9257b5be
# ╠═7606e46a-f984-4f41-9cf2-f0129912d2a4
# ╠═16dc8a56-7354-4dd7-bcc1-e728c367abd0
# ╠═4da291a7-d6d9-4078-a72d-fb50b6212504
# ╠═09641ff0-40c3-4ce7-a932-6e65e9b656ef
# ╠═8a991fb6-67fb-4b3d-a14c-4637b19a8018
# ╠═46820ceb-a88b-4494-9bc3-2de11ecd49e7
# ╠═f66e1ade-056f-4fa1-be44-d16d896d009a
# ╠═4e6ca50e-a966-40c5-8ce5-5a066dcdf534
# ╠═d5e70a4f-a10e-4f42-be76-5d91b275737b
# ╠═ce041037-a8c2-43f7-be16-7fe704dcb55e
# ╠═baecb3c3-8aa4-4b2e-9c42-6460f67915f2
# ╠═d44bca78-d35e-4c1c-8fb7-9502660b0c36
# ╠═b29579f9-954d-4950-a9b3-c40d0f89e831
# ╠═71726ca0-c7ff-4292-9f97-f5dff052cead
# ╠═6bfbc7a8-e793-4489-a135-2d5b4d95c196
# ╠═4f8c1390-b2f4-4b7d-a98c-2ed6ad93b8e0
# ╠═37dc92fc-f7de-4569-b4c3-25ce858add1d
# ╠═ca5fde85-5b99-4b55-a520-07d0a7e670a7
# ╠═093a262a-3345-4180-b042-def2120b7fb6
# ╠═10f6769e-5386-476d-87eb-ac3ec97d0cf4
# ╠═198031bc-48c4-48aa-8bc2-18ff732e5b2a
# ╠═9a8f96a9-9ed5-41f2-b541-5764e6f92f98
# ╠═48b330b1-ba27-43b9-a19c-2fc14de99333
# ╠═ad7967ed-ce1b-4fe8-b348-ee92fab63856
# ╠═118011bf-9f51-42ad-8c4a-383668d85031
# ╠═42e249b8-7cce-42db-b63d-53072c85bbac
# ╠═d329bb22-7af2-4f00-bc70-f903d22ff33d
# ╠═9e126207-9cbf-479f-bdf1-837c69aadca6
# ╠═fe6d4fdc-420c-4035-93b8-31846bf3f405
# ╠═139c4080-7f59-497a-a255-b5c802fc362c
# ╠═e75f2df8-b497-4fba-96fd-0f0773a9c27f
# ╠═baac569e-f403-4026-b541-14f6e9af80c5
