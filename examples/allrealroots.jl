using StaticArrays
using Base.Threads
using BenchmarkTools
using Random

@inline function horner(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::F, y::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	@inbounds w = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	d = zero(F)
	for i in 2:N
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function intervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	@inbounds cohigh = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > eps(F)
            if low > -eps(F)
			    muladd(low, colow, coeff), muladd(high, cohigh, coeff)
            elseif high < eps(F)
                muladd(low, cohigh, coeff), muladd(high, colow, coeff)
            else
                muladd(low, cohigh, coeff), muladd(high, cohigh, coeff)
            end
		elseif cohigh < eps(F)
            if low > -eps(F)
			    muladd(high, colow, coeff), muladd(low, cohigh, coeff)
            elseif high < eps(F)
                muladd(high, cohigh, coeff), muladd(low, colow, coeff)
            else
                muladd(high, colow, coeff), muladd(low, colow, coeff)
            end
		else
            if low > -eps(F)
			    muladd(high, colow, coeff), muladd(high, cohigh, coeff)
            elseif high < eps(F)
                muladd(low, cohigh, coeff), muladd(low, colow, coeff)
            else
                min(muladd(low, cohigh, coeff), muladd(high, colow, coeff)), 
                max(muladd(low, colow, coeff), muladd(high, cohigh, coeff))
            end
        end
	end
    colow, cohigh
end

function allrealrootintervalnewtonrobustregulafalsi(coeffs::NTuple{N,F}, res::MVector{M, F}) where {N, M, NM, F <: AbstractFloat}
	if N === 1
		return 0
	elseif N === 2
		@inbounds res[1] = -coeffs[2] / coeffs[1]
		return 1
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	mm = one(F)
	MM = one(F)
	s = -one(F)
	for i in 2:N-1
		@inbounds coeff = poly[i]
		if coeff < MM; MM = coeff end
		_coeff = s * coeff
		if _coeff < mm; mm = _coeff end
		s = -s
	end
	if mm > zero(F); mm = one(F) end
	if MM > zero(F); MM = one(F) end
	if mm === one(F) && MM === one(F); return 0 end
	domlow = mm - one(F)
	domhigh = one(F) - MM
	list = MVector{2N+1, NTuple{2,F}}(undef)
    counter = 0
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = intervalhorner(domlow, domhigh, poly′)
		#@show mid comid codom′low codom′high
		if codom′low < -eps(F) && codom′high > eps(F)
			leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
				typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
			elseif comid < -eps(F)
				typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
			else
                counter += 1
				@inbounds res[counter] = mid
                if index === 0
                    return counter
                else
                    @inbounds domlow, domhigh = list[index]
                    index -= 1
                end
                continue
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
				return counter
			else
				@inbounds domlow, domhigh = list[index]
				index -= 1
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
            #@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < eps(F)
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
                    #@show mid comid comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-8abs(mid) 
						counter += 1
                        @inbounds res[counter] = newmid
                        if index === 0
                            return counter
                        else
                            @inbounds domlow, domhigh = list[index]
                            index -= 1
			            end
                        break
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
							counter += 1
                            @inbounds res[counter] = mid
                            if index === 0
                                return counter
                            else
                                @inbounds domlow, domhigh = list[index]
                                index -= 1
                            end
                            break
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
				return counter
			else
				@inbounds domlow, domhigh = list[index]
				index -= 1
			end
		end
	end
end

function allrealrootintervalnewtonrobustregulafalsithreads(coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
    if N === 1
		return NTuple{0, F}()
	elseif N === 2
		return NTuple{1, F}(-coeffs[2] / coeffs[1])
	end
    _coeff1 = inv(coeffs[1])
    coeffs = ntuple(N) do i
        _coeff1 * coeffs[i]
    end
    coeffs′ = ntuple(N - 1) do i
        (N-i) * coeffs[i]
    end
	mm = typemax(F)
	MM = typemax(F)
	s = one(F)
	for coeff in coeffs
		if coeff < MM; MM = coeff end
		_coeff = s * coeff
		if _coeff < mm; mm = _coeff end
		s *= -one(F)
	end
	domlow = mm - one(F)
	domhigh = one(F) - MM
    #domlow = min(zero(F), smallest(syms2[2:end]) - one(F))
	#domhigh = max(zero(F), one(F) - smallest(coeffs[2:end]))
	#@show domlow domhigh
    ch = Channel{F}(N-1)
    intervalcheck(ch, coeffs, coeffs′, domlow, domhigh)
    tuple(ch.data...)
end

function intervalcheck(ch::Channel{F}, coeffs::NTuple{N,F}, coeffs′::NTuple{M,F}, domlow::F, domhigh::F) where {N, M, F <: AbstractFloat}
	#println("start: ", Threads.threadid())
    @sync while true
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, coeffs)
		codom′low, codom′high = intervalhorner(domlow, domhigh, coeffs′)
		if codom′low < -eps(F) && codom′high > eps(F)
			leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
				typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
			elseif comid < -eps(F)
				typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
			else
                put!(ch, mid)
                break
			end
			if !(domhigh < leftlow || lefthigh < domlow)
				if !(domhigh < rightlow || righthigh < domlow)
                    newlow = domlow < rightlow ? rightlow : domlow
                    newhigh = domhigh < righthigh ? domhigh : righthigh
					Threads.@spawn intervalcheck(ch, coeffs, coeffs′, $newlow, $newhigh)
				end
				domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
			elseif !(domhigh < rightlow || righthigh < domlow)
				domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
			else
                break
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, coeffs)
			if codomlow * codomhigh < eps(F)
				while true
					comid, comid′ = hornerd(mid, coeffs)
					delta = comid / comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-8abs(mid)
						put!(ch, newmid)
                        break
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
							put!(ch, mid)
                            break
						end
						newmid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
						mid = if domlow < newmid < domhigh
							newmid
						else 
							0.5(domlow + domhigh) # bisection
						end
					end
                end
			end
            break
		end
	end
	#println("stop: ", Threads.threadid())
end


using Polynomials
let N = parse(Int64, ARGS[1])
	p = Tuple(2(rand(N) .- 0.5))
	println(p)
	res = MVector{N-1, Float64}(undef)
	count = allrealrootintervalnewtonrobustregulafalsi(p, res)
	println(count)
	println(res[1:count])
	res = real.(filter(x->isreal(x), roots(Polynomial(p[end:-1:1]))))
	println(res)
	res = allrealrootintervalnewtonrobustregulafalsithreads(p)
	println(res)
	res = MVector{N-1, Float64}(undef)
	bench = @benchmark allrealrootintervalnewtonrobustregulafalsi($p, $res)
	println(minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
	bench = @benchmark roots(Polynomial($(p[end:-1:1])))
	println(minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
	bench = @benchmark allrealrootintervalnewtonrobustregulafalsithreads($p)
	println(minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
	N = 100
	times = Vector{Float64}(undef, N)
	for n in 3:9
		Random.seed!(150)
		res = MVector{n-1, Float64}(undef)
		for i = 1:N
			p = tuple((2(rand(n) .- 0.5))...)
			times[i] = 1000000000 * @belapsed allrealrootintervalnewtonrobustregulafalsi($p, $res)
		end
		println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
	end
end