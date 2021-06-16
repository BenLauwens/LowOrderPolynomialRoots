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

#= @inline function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	m = arg1 < arg2 ? arg1 : arg2
	for arg in args
		m = arg < m ? arg : m
	end
	m
end =#

@generated function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	ex = quote m = arg1 < arg2 ? arg1 : arg2 end
	for i in 1:length(args)
		ex = quote
			$ex
			@inbounds m = args[$i] < m ? args[$i] : m
		end
	end
	ex
end

#= @inline function horner(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	v = coeff1
	for coeff in coeffs
		v = muladd(x, v, coeff)
	end
	v
end =#

@generated function horner(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	l === 0 && return quote coeff1 end
	ex = :(coeff1)
	for i in 1:l
		ex = :(@inbounds muladd(x, $ex, coeffs[$i]))
	end
	ex
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

@generated function smallestpositiverootintervalnewtonrobustregulafalsi(coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	if N === 1
		return quote typemax(F) end
	elseif N === 2
		return quote
			@inbounds ret = -coeffs[2] / coeffs[1]
			if ret < 0
				typemax(F)
			else
				ret
			end
		end
	end
	syms = [gensym() for _ in 1:N]
  	syms′ = [gensym() for _ in 1:N-1]
	ex = quote
		@inbounds _coeff1 = inv(coeffs[1])
		$(syms[1]) = one(F)
	end
	for i in 1:N-1
		ex = quote
			$ex
			$(syms′[i]) = $(N-i) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$(i+1)]
		end
	end
	quote
		$ex
		mm = smallest($(syms[2:N]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{$(2N+1), NTuple{2,F}}(undef)
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
    p = (1.0, -2.0, -1.0)
    res = smallestpositiverootintervalnewtonrobustregulafalsi(p)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($p)
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    p = (1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)
    res = smallestpositiverootintervalnewtonrobustregulafalsi(p)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($p)
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    N = 1000
	times = Vector{Float64}(undef, N)
    for n in 3:9
        Random.seed!(150)
        for i = 1:N
            p = tuple((2(rand(n) .- 0.5))...)
            times[i] = 1000000000 * @belapsed smallestpositiverootintervalnewtonrobustregulafalsi($p)
        end
        println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
    end
end