using StaticArrays
using BenchmarkTools

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
    ex = quote
        $ex
        mm = $(syms[1])
    end
    for i in 2:N
        ex = quote
            $ex
            if $(syms[i]) < mm; mm = $(syms[i]) end
        end
    end
    horner = quote
        comid = $(syms[1])
    end
    for i in 2:N
        horner = quote
            $horner
            comid = muladd(mid, comid, $(syms[i]))
        end
    end
    posintervalhorner = quote
        codom′low, codom′high = $(syms′[1]), $(syms′[1])
    end
    for i in 2:N-1
        posintervalhorner = quote
            $posintervalhorner
            codom′low, codom′high = if codom′low > eps(F)
                muladd(domlow, codom′low, $(syms′[i])), muladd(domhigh,  codom′high, $(syms′[i]))
            elseif codom′high < -eps(F)
                muladd(domhigh, codom′low, $(syms′[i])), muladd(domlow,  codom′high, $(syms′[i]))
            else
                muladd(domhigh, codom′low, $(syms′[i])), muladd(domhigh,  codom′high, $(syms′[i]))
            end
        end
    end
    horner2 = quote
        codomlow, codomhigh = $(syms[1]), $(syms[1])
    end
    for i in 2:N
        horner2 = quote
            $horner2
            codomlow, codomhigh = muladd(domlow, codomlow, $(syms[i])), muladd(domhigh, codomhigh, $(syms[i]))
        end
    end
    hornerd = quote
        comid′, comid = zero(F), $(syms[1])
    end
    for i in 2:N
        hornerd = quote
            $hornerd
            comid′, comid = muladd(mid, comid′, comid), muladd(mid, comid, $(syms[i]))
        end
    end
	quote
		$ex
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{$(2N+1), NTuple{2,F}}(undef)
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			$horner
			$posintervalhorner
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
					break
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			else
				$horner2
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					while true
						$hornerd
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
					break
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			end
		end
        return typemax(F)
	end
end

let
    using Random
    using Polynomials
    p = (1.0, -2.0, -1.0)
    res = smallestpositiverootintervalnewtonrobustregulafalsi(p)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($p)
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    p = (1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)
    res = smallestpositiverootintervalnewtonrobustregulafalsi(p)
    bench = @benchmark smallestpositiverootintervalnewtonrobustregulafalsi($p)
    println(res, '\t', minimum(bench.times), '\t', sum(bench.times)/ length(bench.times), '\t', maximum(bench.times))
    N = 100
	times = Vector{Float64}(undef, N)
    for n in 3:9
        #Random.seed!(150)
        for i = 1:N
            p = tuple((2(rand(n) .- 0.5))...)
            times[i] = 1000000000 * @belapsed smallestpositiverootintervalnewtonrobustregulafalsi($p)
        end
        println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
    end
end

