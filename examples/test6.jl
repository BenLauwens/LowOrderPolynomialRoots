using Random
using BenchmarkTools
using StaticArrays

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

@inline function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	m = arg1 < arg2 ? arg1 : arg2
	for arg in args
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	v = coeff1
	for i in eachindex(coeffs)
		@inbounds v = muladd(x, v, coeffs[i])
	end
	v
end

@inline function horner2(x::F, y::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	v = coeff1
	w = coeff1
	for coeff in coeffs
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	v = coeff1
	d = zero(F)
	for coeff in coeffs
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function posintervalhorner(low::F, high::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	colow, cohigh = coeff1, coeff1
	for coeff in coeffs
		colow, cohigh = if colow > eps(F)
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh < -eps(F)
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end

@generated function smallestpositiverootintervalnewtonrobustregulafalsi(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	if l === 0
		return quote typemax(F) end
	elseif l === 1
		return quote 
			res = -coeffs[1] / coeff1
			if res < -zero(F)
				typemax(F)
			else
				res
			end
		end
	end
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
	quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{$(2l+3), NTuple{2,F}}(undef)
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
				#@show leftlow lefthigh rightlow righthigh
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
						@inbounds list[index] = (domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
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
				codomlow, codomhigh = horner2(domlow, domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					while true
						comid, comid′ = hornerd(mid, $(syms...))
						delta = comid / comid′
						newmid = mid - delta
						if abs(delta) < 1.0e-8mid 
							return newmid
						elseif domlow < newmid < domhigh
							mid = newmid
						else
							if comid * codomlow > -eps(F)
								domlow = mid
								codomlow = comid
							elseif comid * codomhigh > -eps(F)
								domhigh = mid
								codomhigh = comid
							else
								return mid
							end
							newmid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							mid = if domlow < newmid < domhigh
								newmid
							else 
								0.5(domlow + domhigh) # bisection
							end
						end
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
end

let
    res = quadraticsmallestpositiveroot(1.0, -2.0, -1.0)
    bench = @benchmark quadraticsmallestpositiveroot(1.0, -2.0, -1.0) samples=10000 evals=1000
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times)) 
    res = smallestpositiverootintervalnewtonrobustregulafalsi(1.0, -2.0, -1.0)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(1.0, -2.0, -1.0) samples=10000 evals=1000
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    res = smallestpositiverootintervalnewtonrobustregulafalsi(1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    Random.seed!(150)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[1], p[2], p[3]) setup=(p = tuple((2(rand(3) .- 0.5))...)) samples=10000 evals=1000
    println(3, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    Random.seed!(150)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[1], p[2], p[3], p[4]) setup=(p= tuple((2(rand(4) .- 0.5))...)) samples=10000 evals=1000
    println(4, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    Random.seed!(150)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[1], p[2], p[3], p[4], p[5]) setup=(p= tuple((2(rand(5) .- 0.5))...)) samples=10000 evals=1000
    println(5, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    Random.seed!(150)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[1], p[2], p[3], p[4], p[5], p[6]) setup=(p= tuple((2(rand(6) .- 0.5))...)) samples=10000 evals=1000
    println(6, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    Random.seed!(150)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[1], p[2], p[3], p[4], p[5], p[6], p[7]) setup=(p= tuple((2(rand(7) .- 0.5))...)) samples=10000 evals=1000
    println(7, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    Random.seed!(150)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]) setup=(p= tuple((2(rand(8) .- 0.5))...)) samples=10000 evals=1000
    println(8, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    Random.seed!(150)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]) setup=(p= tuple((2(rand(9) .- 0.5))...)) samples=10000 evals=1000
    println(9, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times)) 
	Random.seed!(150)
	p = tuple((2(rand(9) .- 0.5))...)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($(p[1]), $(p[2]), $(p[3]), $(p[4]), $(p[5]), $(p[6]), $(p[7]), $(p[8]), $(p[9])) samples=10000 evals=1000
    println(9, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
	N = 100
	Random.seed!(150)
	times = Vector{Float64}(undef, N)
	for i = 1:N
		p = tuple((2(rand(9) .- 0.5))...)
		times[i] = 1000000000 * @belapsed smallestpositiverootintervalnewtonrobustregulafalsi($(p[1]), $(p[2]), $(p[3]), $(p[4]), $(p[5]), $(p[6]), $(p[7]), $(p[8]), $(p[9]))
	end
	println(9, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times))
end