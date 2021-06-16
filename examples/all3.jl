using StaticArrays

@inline function smallest(args::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds m = args[1]
	for i in 2:N
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end

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

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	@inbounds cohigh = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
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
		@inbounds ret = -coeffs[2] / coeffs[1]
		if ret < zero(F)
			return typemax(F)
		else
			return ret
		end
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	MM = smallest(poly)
	if MM > eps(F) return typemax(F) end
	domlow, domhigh = zero(F), one(F) - MM
	list = MVector{2N+1, NTuple{2,F}}(undef)
	#list = MMatrix{2, 2N+1, F}(undef)
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, poly′)
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
					#@inbounds list[index] = (domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
                    GC.@preserve list unsafe_store!(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), convert(NTuple{2, F}, (domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)), index)
					#@inbounds list[1, index] = domlow < rightlow ? rightlow : domlow
					#@inbounds list[2, index] = domhigh < righthigh ? domhigh : righthigh
				end
				domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
			elseif !(domhigh < rightlow || righthigh < domlow)
				domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
			elseif index === 0
				return typemax(F)
			else
				#@inbounds domlow, domhigh = list[index]
                domlow, domhigh = GC.@preserve list unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), index)
                #domlow = list[1, index]
				#domhigh = list[2, index]
				index -= 1
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
			#@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < eps(F)
				while true
					comid, comid′ = hornerd(mid, poly)
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
				#@inbounds domlow, domhigh = list[index]
                domlow, domhigh = GC.@preserve list unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), index)
                #domlow = list[1, index]
				#domhigh = list[2, index]
				index -= 1
			end
		end
	end
end

function allrealrootintervalnewtonrobustregulafalsi(coeffs::NTuple{N,F}, res::MVector{M, F}) where {N, M, F <: AbstractFloat}
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
				#@inbounds res[counter] = mid
                GC.@preserve res unsafe_store!(Base.unsafe_convert(Ptr{F}, pointer_from_objref(res)), convert(F, mid), counter)
                if index === 0
                    return counter
                else
                    #@inbounds domlow, domhigh = list[index]
                    domlow, domhigh = GC.@preserve list unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), index)
                    index -= 1
                end
                continue
			end
			#@show leftlow lefthigh rightlow righthigh
			if !(domhigh < leftlow || lefthigh < domlow)
				if !(domhigh < rightlow || righthigh < domlow)
					index += 1
					#@inbounds list[index] = (domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
                    GC.@preserve list unsafe_store!(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), convert(NTuple{2, F}, (domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)), index)
				end
				domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
			elseif !(domhigh < rightlow || righthigh < domlow)
				domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
			elseif index === 0
				return counter
			else
				#@inbounds domlow, domhigh = list[index]
                domlow, domhigh = GC.@preserve list unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), index)
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
                        #@inbounds res[counter] = newmid
                        GC.@preserve res unsafe_store!(Base.unsafe_convert(Ptr{F}, pointer_from_objref(res)), convert(F, newmid), counter)
                        if index === 0
                            return counter
                        else
                            #@inbounds domlow, domhigh = list[index]
                            domlow, domhigh = GC.@preserve list unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), index)
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
                            #@inbounds res[counter] = mid
                            GC.@preserve res unsafe_store!(Base.unsafe_convert(Ptr{F}, pointer_from_objref(res)), convert(F, mid), counter)
                            if index === 0
                                return counter
                            else
                                #@inbounds domlow, domhigh = list[index]
                                domlow, domhigh = GC.@preserve list unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), index)
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
				#@inbounds domlow, domhigh = list[index]
                domlow, domhigh = GC.@preserve list unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(list)), index)
                index -= 1
			end
		end
	end
end

using Random
using BenchmarkTools
using ProgressMeter
using UnicodePlots

let N = 1000
	for n in 3:15
		#Random.seed!(150)
		times = Vector{Float64}(undef, N)
		res = MVector{n-1, Float64}(undef)
		@showprogress for i = 1:N
			p = tuple((2(rand(n) .- 0.5))...)
			times[i] = 1000000000 * @belapsed allrealrootintervalnewtonrobustregulafalsi($p, $res)
		end
		show(histogram(times))
		println()
		println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
	end
	for n in 3:15
        #Random.seed!(150)
		times = Vector{Float64}(undef, N)
        @showprogress for i = 1:N
            p = tuple((2(rand(n) .- 0.5))...)
            times[i] = 1000000000 * @belapsed smallestpositiverootintervalnewtonrobustregulafalsi($p)
        end
		show(histogram(times))
		println()
        println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
    end
end