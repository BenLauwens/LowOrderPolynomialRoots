### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 2c2f2e31-3886-4782-adda-6563dc355905
begin
	using Polynomials
	using Random
	using PlutoUI
	using BenchmarkTools
end

# ╔═╡ 3f521521-cbfa-4305-8880-e40f3a37ceca
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

# ╔═╡ dcbb99d9-55e5-4f0a-9819-a79cf166fc8f
quadraticsmallestpositiveroot(1.0, -2.0, -1.0)

# ╔═╡ 9c97ce34-66f5-4786-b9b9-91e9aa15eef5
@benchmark quadraticsmallestpositiveroot(1.0, -2.0, -1.0)

# ╔═╡ 7dabb617-af0b-4b69-b363-f79846213d91
mutable struct Tupple2FloatList{N, F<:AbstractFloat}
	data :: NTuple{N, NTuple{2, F}}
	index :: Int
	function Tupple2FloatList{N, F}() where {N, F<:AbstractFloat}
		t2fl = new{N, F}()
		t2fl.index = 0
		t2fl
	end
end

# ╔═╡ 5f42573b-7aee-43a6-9e51-aa6e55de85a6
@inline function Base.push!(t2fl::Tupple2FloatList{N, F}, val::F...) where {N, F<:AbstractFloat}
	i = t2fl.index + 1
	GC.@preserve t2fl unsafe_store!(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(t2fl)), convert(NTuple{2, F}, val), i)
	t2fl.index = i
	t2fl
end

# ╔═╡ 4d7f0ada-784e-47d0-8661-3c8765fe7366
@inline function Base.pop!(t2fl::Tupple2FloatList{N, F} ) where {N, F<:AbstractFloat}
	i = t2fl.index
	val = GC.@preserve t2fl unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(t2fl)), i)
	#val = t2fl.data[i]
	t2fl.index = i - 1
	val
end

# ╔═╡ 5061a440-2689-4805-b17e-de20fe41a0e6
@inline function Base.length(t2fl::Tupple2FloatList{N, F} ) where {N, F<:AbstractFloat}
	t2fl.index
end

# ╔═╡ 6b270c32-d475-4543-bef1-f4338713258a
@inline function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	m = arg1 < arg2 ? arg1 : arg2
	for arg in args
		m = arg < m ? arg : m
	end
	m
end

# ╔═╡ 53759139-bd86-4ca8-b406-c5c684ad7a65
@inline function horner(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	v = coeff1
	for i in eachindex(coeffs)
		@inbounds v = muladd(x, v, coeffs[i])
	end
	v
end

# ╔═╡ 1e3f662b-c76e-43ee-a134-060bf8d3aadc
@inline function horner2(x::F, y::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	v = coeff1
	w = coeff1
	for coeff in coeffs
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

# ╔═╡ 212617f3-4748-450c-9e50-8501704f509a
@inline function hornerd(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
	v = coeff1
	d = zero(F)
	for coeff in coeffs
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

# ╔═╡ ee4acb92-0f90-4c12-b3c3-f0c32fab6683
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

# ╔═╡ c0f50886-5b2c-4c63-a2ec-3b8248779a07
@generated function smallestpositiverootintervalnewtonrobustregulafalsigetset(coeff1::F, coeffs::F...) where F <: AbstractFloat
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
			$(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		index = 0
	end
	bounds = [gensym() for _ in 1:4l+6]
	for i in 1:4l+6
		ex = quote
			$ex
			$(bounds[i]) = zero(F)
		end
	end
	set = quote index += 1 end
	for i in 1:2l+3
		set = quote
			$set
			if $i === index
				$(bounds[i]) = domlow < rightlow ? rightlow : domlow
				$(bounds[i+2l+3]) = domhigh < righthigh ? domhigh : righthigh
				@goto endset 
			end
		end
	end
	set = quote
		$set
		@label endset
	end
	get = quote end
	for i in 1:2l+3
		get = quote
			$get
			if $i === index
				domlow = $(bounds[i])
				domhigh = $(bounds[i+2l+3])
				index -= 1
				@goto endget
			end
		end
	end
	quote
		$ex
		while true
			@label endget
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
						$set
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					$get
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
					$get
				end
			end
		end
	end
end

# ╔═╡ 3349b746-ed88-450b-9bfb-ebf0602e64db
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
			$(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = Tupple2FloatList{$(2l+3), F}()
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
						push!(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
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
				elseif length(list) === 0
					return typemax(F)
				else
					domlow, domhigh = pop!(list)
				end
			end
		end
	end
end

# ╔═╡ 1ca4a294-d9df-4c89-b54d-f45045624a58
smallestpositiverootintervalnewtonrobustregulafalsigetset(1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)

# ╔═╡ d7b71da5-2b8d-4a6a-807a-f0bc9c8b80ad
@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)

# ╔═╡ 3e5299d2-1e60-47fd-9d6f-19be531cbab3
smallestpositiverootintervalnewtonrobustregulafalsi(1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)

# ╔═╡ da011543-f010-491d-b145-1095b59b7a49
@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(1.0, -2.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)

# ╔═╡ 0c8c2a93-bd57-4b0e-b9c4-68efabfd006e
with_terminal() do
	Random.seed!(rand(1:150))
	i = 5
	for j in 1:1000000
		a = vec(2(rand(1,i) .- 0.5))
		pol = filter(r->isreal(r) && real(r) > -eps(), Polynomials.roots(Polynomial(a)))
		p = length(pol) === 0 ? Inf : minimum(real.(pol))
		q = smallestpositiverootintervalnewtonrobustregulafalsi(a[end:-1:1]...)
		if abs(p - q) > 1e-7p
			@show a pol p q
			break
		end
	end
end

# ╔═╡ a76ea3db-84dd-449c-a9b8-842ef6a5ea46
let
	Random.seed!(150)
	n = 3
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ f83f1da6-649f-4b20-9806-57ad4d6b2b1d
let
	Random.seed!(150)
	n = 3
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 316c5ebe-dae0-4cc8-b4a7-b6e3b3f9b469
let
	Random.seed!(150)
	n = 4
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 7acbe692-db30-4e57-bb2e-476bf6b7bb73
let
	Random.seed!(150)
	n = 4
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 18c8d847-5240-42ba-8c6b-7f46251f156d
let
	Random.seed!(150)
	n = 5
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 4c461d39-e086-4b42-a2e5-864d9d618567
let
	Random.seed!(150)
	n = 5
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ ae097f56-99f3-49f2-aaa7-32a88e456d66
let
	Random.seed!(150)
	n = 6
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 6d3594bd-ffb5-4fb1-b02c-169d6b825902
let
	Random.seed!(150)
	n = 6
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ e2f61960-c3cb-4ad6-9c78-2c018df8af24
let
	Random.seed!(150)
	n = 7
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 66641ad7-f507-47f4-8743-7ec2912251c4
let
	Random.seed!(150)
	n = 7
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 1f57fc65-3567-4768-a780-e15d52c04dc5
let
	Random.seed!(150)
	n = 8
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 660be593-97c6-4a1b-877e-3acd8608ebfe
let
	Random.seed!(150)
	n = 8
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ fdba44ce-2e73-4b97-bc08-f897dc4e306f
let
	Random.seed!(150)
	n = 9
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 8d633551-04bf-45a7-8ba0-4071d1ee65f0
let
	Random.seed!(150)
	n = 9
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ 23d49206-e78d-48c5-a2b8-514c0540fb62
let
	Random.seed!(150)
	n = 10
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsi(p[10], p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ f71ac031-0167-43f9-8dd9-783e2d11b0ba
let
	Random.seed!(150)
	n = 10
	@benchmark smallestpositiverootintervalnewtonrobustregulafalsigetset(p[10], p[9], p[8], p[7], p[6], p[5], p[4], p[3], p[2], p[1]) setup=(p=vec(2(rand(1,$n) .- 0.5)))
end

# ╔═╡ Cell order:
# ╠═3f521521-cbfa-4305-8880-e40f3a37ceca
# ╠═dcbb99d9-55e5-4f0a-9819-a79cf166fc8f
# ╠═9c97ce34-66f5-4786-b9b9-91e9aa15eef5
# ╠═7dabb617-af0b-4b69-b363-f79846213d91
# ╠═5f42573b-7aee-43a6-9e51-aa6e55de85a6
# ╠═4d7f0ada-784e-47d0-8661-3c8765fe7366
# ╠═5061a440-2689-4805-b17e-de20fe41a0e6
# ╠═6b270c32-d475-4543-bef1-f4338713258a
# ╠═53759139-bd86-4ca8-b406-c5c684ad7a65
# ╠═1e3f662b-c76e-43ee-a134-060bf8d3aadc
# ╠═212617f3-4748-450c-9e50-8501704f509a
# ╠═ee4acb92-0f90-4c12-b3c3-f0c32fab6683
# ╠═c0f50886-5b2c-4c63-a2ec-3b8248779a07
# ╠═3349b746-ed88-450b-9bfb-ebf0602e64db
# ╠═1ca4a294-d9df-4c89-b54d-f45045624a58
# ╠═d7b71da5-2b8d-4a6a-807a-f0bc9c8b80ad
# ╠═3e5299d2-1e60-47fd-9d6f-19be531cbab3
# ╠═da011543-f010-491d-b145-1095b59b7a49
# ╠═2c2f2e31-3886-4782-adda-6563dc355905
# ╠═0c8c2a93-bd57-4b0e-b9c4-68efabfd006e
# ╠═a76ea3db-84dd-449c-a9b8-842ef6a5ea46
# ╠═f83f1da6-649f-4b20-9806-57ad4d6b2b1d
# ╠═316c5ebe-dae0-4cc8-b4a7-b6e3b3f9b469
# ╠═7acbe692-db30-4e57-bb2e-476bf6b7bb73
# ╠═18c8d847-5240-42ba-8c6b-7f46251f156d
# ╠═4c461d39-e086-4b42-a2e5-864d9d618567
# ╠═ae097f56-99f3-49f2-aaa7-32a88e456d66
# ╠═6d3594bd-ffb5-4fb1-b02c-169d6b825902
# ╠═e2f61960-c3cb-4ad6-9c78-2c018df8af24
# ╠═66641ad7-f507-47f4-8743-7ec2912251c4
# ╠═1f57fc65-3567-4768-a780-e15d52c04dc5
# ╠═660be593-97c6-4a1b-877e-3acd8608ebfe
# ╠═fdba44ce-2e73-4b97-bc08-f897dc4e306f
# ╠═8d633551-04bf-45a7-8ba0-4071d1ee65f0
# ╠═23d49206-e78d-48c5-a2b8-514c0540fb62
# ╠═f71ac031-0167-43f9-8dd9-783e2d11b0ba
