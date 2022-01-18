@generated function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, lows::Ptr{F}, highs::Ptr{F}) where {N, F <: AbstractFloat}
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
            codom′low, codom′high = if codom′low > zero(F)
                muladd(domlow, codom′low, $(syms′[i])), muladd(domhigh,  codom′high, $(syms′[i]))
            elseif codom′high < zero(F)
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
		if mm > zero(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			$horner
			$posintervalhorner
			#@show comid codom′low codom′high
			if codom′low < zero(F) < codom′high
				leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
					domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
				else
					domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
				end
				#@show leftlow lefthigh rightlow righthigh
				if leftlow < lefthigh
					if rightlow < righthigh
						index += 1
						unsafe_store!(lows, rightlow, index)
						unsafe_store!(highs, righthigh, index)
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
				if codomlow * codomhigh < zero(F)
					while true
						$hornerd
						delta = comid / comid′
						newmid = mid - delta
						if abs(delta) < 1.0e-8mid 
							return newmid
						elseif domlow < newmid < domhigh
							mid = newmid
						else
							if comid * codomlow < zero(F)
								domhigh, codomhigh = mid, comid
							else
								domlow, codomlow = mid, comid	
							end
							mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
						end
					end
				end
			end
			if index == 0 break end
			domlow = unsafe_load(lows, index)
			domhigh = unsafe_load(highs, index)
			index -= 1
		end
		return typemax(F)
	end
end

@generated function allrealrootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, res::Ptr{F}, lows::Ptr{F}, highs::Ptr{F}) where {N, F <: AbstractFloat}
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
    posintervalhorner = quote end
	negintervalhorner = quote end
    for i in 2:N-1
        posintervalhorner = quote
            $posintervalhorner
            codom′low, codom′high = if codom′low > zero(F)
                muladd(domlow, codom′low, $(syms′[i])), muladd(domhigh,  codom′high, $(syms′[i]))
            elseif codom′high < zero(F)
                muladd(domhigh, codom′low, $(syms′[i])), muladd(domlow,  codom′high, $(syms′[i]))
            else
                muladd(domhigh, codom′low, $(syms′[i])), muladd(domhigh,  codom′high, $(syms′[i]))
            end
        end
		negintervalhorner = quote
            $negintervalhorner
            codom′low, codom′high = if codom′low > zero(F)
                muladd(domlow, codom′high, $(syms′[i])), muladd(domhigh, codom′low, $(syms′[i]))
            elseif codom′high < zero(F)
                muladd(domhigh, codom′high, $(syms′[i])), muladd(domlow, codom′low, $(syms′[i]))
            else
                muladd(domlow, codom′high, $(syms′[i])), muladd(domlow, codom′low, $(syms′[i]))
            end
        end
    end
	intervalhorner = quote
        codom′low, codom′high = $(syms′[1]), $(syms′[1])
		if domlow < zero(F)
			$negintervalhorner
		else
			$posintervalhorner
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
		index = 0
		domlow, domhigh = if mm < zero(F)
			if MM < zero(F)
				index = 1
				unsafe_store!(lows, zero(F), 1)
				unsafe_store!(highs, one(F) - MM, 1)
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
				leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
					domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
				else
					domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
				end
				#@show leftlow lefthigh rightlow righthigh
				if leftlow < lefthigh
					if rightlow < righthigh
						index += 1
						unsafe_store!(lows, rightlow, index)
						unsafe_store!(highs, righthigh, index)
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
				if codomlow * codomhigh < zero(F)
					while true
						$hornerd
                        #@show mid comid comid′
						delta = comid / comid′
						newmid = mid - delta
						if abs(delta) < 1.0e-8abs(mid) 
							counter += 1
                            unsafe_store!(res, newmid, counter)
                            break
						elseif domlow < newmid < domhigh
							mid = newmid
						else
							if comid * codomlow < zero(F)
								domhigh, codomhigh = mid, comid
							else
								domlow, codomlow = mid, comid
							end
							mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
						end
					end
				end
			end
			if index == 0 || counter == $(N-1) break end
			domlow = unsafe_load(lows, index)
			domhigh = unsafe_load(highs, index)
			index -= 1
		end
		return counter
	end
end

using BenchmarkTools
using ProgressMeter
using UnicodePlots

let N = 1000
	times = Vector{Float64}(undef, N)
	for n in 3:9
		res = pointer(Vector{Float64}(undef, n-1))
		lows = pointer(Vector{Float64}(undef, 2n+1))
		highs = pointer(Vector{Float64}(undef, 2n+1))
		open("polynomials$n.txt", "r") do io
			@showprogress for i = 1:N
				str = readline(io)
				p = tuple(parse.(Float64, split(str))...)
				times[i] = 1000000000 * @belapsed allrealrootintervalnewtonregulafalsi($p, $res, $lows, $highs) samples=10000 evals=1000
			end
		end
		show(histogram(times))
		println()
		println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		open("out5s.txt", "a") do io
			println(io, n, "\t", minimum(times), "\t", sum(times)/ N, "\t", maximum(times), "\t", sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		end
	end
	for n in 3:9
		lows = pointer(Vector{Float64}(undef, 2n+1))
		highs = pointer(Vector{Float64}(undef, 2n+1))
		open("polynomials$n.txt", "r") do io
			@showprogress for i = 1:N
				str = readline(io)
				p = tuple(parse.(Float64, split(str))...)
            	times[i] = 1000000000 * @belapsed smallestpositiverootintervalnewtonregulafalsi($p, $lows, $highs) samples=10000 evals=1000
			end
        end
		show(histogram(times))
		println()
        println(n, '\t', minimum(times), '\t', sum(times)/ N, '\t', maximum(times), '\t', sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		open("out5s.txt", "a") do io
			println(io, n, "\t", minimum(times), "\t", sum(times)/ N, "\t", maximum(times), "\t", sqrt((sum(times.^2)-sum(times)^2/N)/(N-1)))
		end
    end
end

