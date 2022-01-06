using StaticArrays

@generated function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	if N == 1
		return quote typemax(F) end
	elseif N == 2
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
	lows = [gensym() for _ in 1:2N+1]
	highs = [gensym() for _ in 1:2N+1]
	init_expr = quote 
		index = 0
	end
	for i in 1:2N+1
		init_expr = quote
			$init_expr
			$(lows[i]) = zero(F)
			$(highs[i]) = zero(F)
		end
	end
	push_expr = quote 
		index += 1
	end
	for i in 1:2N+1
		push_expr = quote
			$push_expr
			if index == $i
				$(lows[i]) = rightlow
				$(highs[i]) = righthigh
				@goto push_end
			end
		end
	end
	push_expr = quote
		$push_expr
		@label push_end
	end
	pop_expr = quote end
	for i in 1:2N+1
		pop_expr = quote
			$pop_expr
			if index == $i
				domlow = $(lows[i])
				domhigh = $(highs[i])
				@goto pop_end
			end
		end
	end
	pop_expr = quote
		$pop_expr
		@label pop_end
		index -= 1
	end
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
		$(Expr(:meta, :inline))
		$ex
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		$init_expr 
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			$horner
			$posintervalhorner
			#@show comid codom′low codom′high
			if codom′low < zero(F) < codom′high
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
				elseif comid < -eps(F)
					domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
				else
					return mid
				end
				#@show leftlow lefthigh rightlow righthigh
				if leftlow < lefthigh
					if rightlow < righthigh
						$push_expr
					end
					domlow, domhigh = leftlow, lefthigh
					continue
				elseif rightlow < righthigh
					domlow, domhigh = rightlow, righthigh
					continue
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
							if comid * codomlow > eps(F)
								domlow, codomlow = mid, comid
							elseif comid * codomhigh > eps(F)
								domhigh, codomhigh = mid, comid
							else
								return mid
							end
							mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
						end
					end
				end
			end
			if index == 0
				return typemax(F)
			else
				$pop_expr
			end
		end
	end
end

@generated function allrealrootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, res::MVector{M, F}) where {N, M, F <: AbstractFloat}
	if N == 1
		return quote return 0 end
	elseif N == 2
		return quote
			@inbounds res[1] = -coeffs[2] / coeffs[1]
            return 1
		end
	end
	syms = [gensym() for _ in 1:N]
  	syms′ = [gensym() for _ in 1:N-1]
	lows = [gensym() for _ in 1:2N+1]
	highs = [gensym() for _ in 1:2N+1]
	init_expr = quote 
		index = 0
	end
	for i in 1:2N+1
		init_expr = quote
			$init_expr
			$(lows[i]) = zero(F)
			$(highs[i]) = zero(F)
		end
	end
	push_expr = quote 
		index += 1
	end
	for i in 1:2N+1
		push_expr = quote
			$push_expr
			if index == $i
				$(lows[i]) = rightlow
				$(highs[i]) = righthigh
				@goto push_end
			end
		end
	end
	push_expr = quote
		$push_expr
		@label push_end
	end
	pop_expr = quote end
	for i in 1:2N+1
		pop_expr = quote
			$pop_expr
			if index == $i
				domlow = $(lows[i])
				domhigh = $(highs[i])
				@goto pop_end
			end
		end
	end
	pop_expr = quote
		$pop_expr
		@label pop_end
		index -= 1
	end
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
        MM = $(syms[1])
        s = -one(F)
    end
    for i in 2:N
        ex = quote
            $ex
            if $(syms[i]) < MM; MM = $(syms[i]) end
            coeff = s * $(syms[i])
            if coeff < mm; mm = coeff end
            s = -s
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
    intervalhorner = quote
        codom′low, codom′high = $(syms′[1]), $(syms′[1])
    end
    for i in 2:N-1
        intervalhorner = quote
            $intervalhorner
            codom′low, codom′high = if codom′low > eps(F)
                if domlow > -eps(F)
                    muladd(domlow, codom′low, $(syms′[i])), muladd(domhigh,  codom′high, $(syms′[i]))
                elseif domhigh < eps(F)
                    muladd(domlow, codom′high, $(syms′[i])), muladd(domhigh, codom′low, $(syms′[i]))
                else
                    muladd(domlow, codom′high, $(syms′[i])), muladd(domhigh, codom′high, $(syms′[i]))
                end
            elseif codom′high < -eps(F)
                if domlow > -eps(F)
                    muladd(domhigh, codom′low, $(syms′[i])), muladd(domlow,  codom′high, $(syms′[i]))
                elseif domhigh < eps(F)
                    muladd(domhigh, codom′high, $(syms′[i])), muladd(domlow, codom′low, $(syms′[i]))
                else
                    muladd(domhigh, codom′low, $(syms′[i])), muladd(domlow, codom′low, $(syms′[i]))
                end
            else
                if domlow > -eps(F)
                    muladd(domhigh, codom′low, $(syms′[i])), muladd(domhigh,  codom′high, $(syms′[i]))
                elseif domhigh < eps(F)
                    muladd(domlow, codom′high, $(syms′[i])), muladd(domlow, codom′low, $(syms′[i]))
                else
                    min(muladd(domlow, codom′high, $(syms′[i])), muladd(domhigh, codom′low, $(syms′[i]))), 
                    max(muladd(domlow, codom′low, $(syms′[i])), muladd(domhigh, codom′high, $(syms′[i])))
                end
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
		$(Expr(:meta, :inline))
		$ex
		if mm == zero(F) && MM == zero(F) return 0 end
		$init_expr
		domlow, domhigh = if mm < zero(F)
			if MM < zero(F)
				index = 1
				$(lows[1]) = zero(F)
				$(highs[1]) = one(F) - MM
			end
			mm - one(F), zero(F)
		else
			zero(F), one(F) - MM
		end
		counter = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			$horner
			$intervalhorner
			#@show mid comid codom′low codom′high
			if codom′low < zero(F) < codom′high
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
				elseif comid < -eps(F)
					domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
				else
					counter += 1
                    @inbounds res[counter] = mid
					zero(F), zero(F), zero(F), zero(F)
				end
				#@show leftlow lefthigh rightlow righthigh
				if leftlow < lefthigh
					if rightlow < righthigh
						$push_expr
					end
					domlow, domhigh = leftlow, lefthigh
					continue
				elseif rightlow < righthigh
					domlow, domhigh = rightlow, righthigh
					continue
				end
			else
				$horner2
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					while true
						$hornerd
                        #@show mid comid comid′
						delta = comid / comid′
						newmid = mid - delta
						if abs(delta) < 1.0e-8abs(mid) 
							counter += 1
                            @inbounds res[counter] = newmid
                            break
						elseif domlow < newmid < domhigh
							mid = newmid
						else
							if comid * codomlow > eps(F)
								domlow, codomlow = mid, comid
							elseif comid * codomhigh > eps(F)
								domhigh, codomhigh = mid, comid
							else
								counter += 1
                                @inbounds res[counter] = mid
                                break
							end
							mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
						end
					end
				end
			end
			if index == 0
				return counter
			else
				$pop_expr
			end
		end
	end
end

using Random
using BenchmarkTools
using ProgressMeter
using UnicodePlots

let N = 1000
	for n in 3:9
		Random.seed!(150)
		times = Vector{Float64}(undef, N)
		res = MVector{n-1, Float64}(undef)
		@showprogress for i = 1:N
			p = tuple((2(rand(n) .- 0.5))...)
			times[i] = 1000000000 * @belapsed allrealrootintervalnewtonregulafalsi($p, $res)
		end
		show(histogram(times))
		println()
		println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		open("out7.txt", "a") do io
			println(io, n, "\t", minimum(times), "\t", sum(times)/ N, "\t", maximum(times), "\t", sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		end
	end
	for n in 3:9
        Random.seed!(150)
		times = Vector{Float64}(undef, N)
        @showprogress for i = 1:N
            p = tuple((2(rand(n) .- 0.5))...)
            times[i] = 1000000000 * @belapsed smallestpositiverootintervalnewtonregulafalsi($p)
        end
		show(histogram(times))
		println()
        println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		open("out7.txt", "a") do io
			println(io, n, "\t", minimum(times), "\t", sum(times)/ N, "\t", maximum(times), "\t", sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		end
    end
end
