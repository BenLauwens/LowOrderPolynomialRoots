using BenchmarkTools
using Random

@inline function smallest(args::NTuple{N,F}) where {N, F <: AbstractFloat}
	m = typemax(F)
	for arg in args
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	v = zero(F)
	for coeff in coeffs
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::F, y::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	v = zero(F)
	w = zero(F)
	for coeff in coeffs
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	v = zero(F)
	d = zero(F)
	for coeff in coeffs
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	colow = zero(F)
	cohigh = zero(F)
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

function smallestpositiverootintervalnewtonrobustregulafalsi(coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	if N === 1
		return typemax(F)
	elseif N === 2
		ret = -coeffs[2] / coeffs[1]
		if ret < 0
			return typemax(F)
		else
			return ret
		end
	end
    _coeff1 = inv(coeffs[1])
    syms = ntuple(N) do i
        _coeff1 * coeffs[i]
    end
    syms′ = ntuple(N - 1) do i
        (N-i) * coeffs[i]
    end
    mm = smallest(syms[2:end])
	if mm > eps(F) return typemax(F) end
	domlow, domhigh = zero(F), one(F) - mm
	list = Vector{NTuple{2,F}}(undef, 2N+1)
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, syms)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, syms′)
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
			codomlow, codomhigh = horner2(domlow, domhigh, syms)
			if codomlow * codomhigh < eps(F)
				while true
					comid, comid′ = hornerd(mid, syms)
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
let
    p = (1.0, -2.0, -1.0)
    res = smallestpositiverootintervalnewtonrobustregulafalsi(p)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($p) samples=10000 evals=1000
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    p = (1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)
    res = smallestpositiverootintervalnewtonrobustregulafalsi(p)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($p) samples=10000 evals=1000
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    for n in 3:9
        Random.seed!(150)
        bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p) setup=(p=tuple((2(rand($n) .- 0.5))...)) samples=10000 evals=1000 
        println(n, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    end
    Random.seed!(150)
	p = Tuple(2(rand(9) .- 0.5))
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($p) samples=10000 evals=1000
    println(9, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    N = 100
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