### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 003643f0-a1df-11ef-02bc-a3bb1aac6d2e
using DelimitedFiles, LsqFit, Plots, PlutoUI,  FFTW, Loess, Glob, NLsolve

# ╔═╡ 3f2957da-93d6-43bc-b111-0c110ea5647a
plotly()

# ╔═╡ f03eb307-16cc-428f-924e-e057181428d4
TableOfContents()

# ╔═╡ 5290a6ce-89ae-4ed9-b02f-f1e6033e07dd
md"""# Analyze __*I(V) vs. T*__ en bloc
to do:
- hi-res surface image!
- parametric dif using ``\Delta (T)`` vector? ...parameters relative to ``\Delta``
- separate en bloc fitting that reads .txt files from individual porcessing
- add Specs filetype
"""

# ╔═╡ 833c981f-d250-4194-8454-ad10be888005
md"""#### Define functions"""

# ╔═╡ 432f4d98-109c-49e6-8be7-eba090783fc7
#linear fit function of equidistant values
function linear_fit(ydata)
	#x coordinate is the index
	xdata = collect(eachindex(ydata)) #more general than 1:length*(ydata)
    # Define the model function
    model(x, p) = p[1] .* x .+ p[2]
    # Initial guess for the parameters
    p0 = [1.0, 0.0]
    # Perform the curve fitting
    fit = curve_fit(model, xdata, ydata, p0)
    # Extract the fitting parameters
    slope = fit.param[1]
    intercept = fit.param[2]
    return slope, intercept
end

# ╔═╡ f4d390b1-4da9-41a7-bfb1-1c1b3d0d89a1
#differentiation using linear regression over a predifined number of points
function pts_dif(current, pts)
    n = length(current)
    dIdV = similar(current)  # creates an array of the same type and size

	for i in 1:n
		# range adjusts for left and right boundaries
		# div outputs the quotient from Euclidean division. Computes x/y, truncated to an integer.
		range = max(1, i - pts[i]) : min(n, i + pts[i])
        dIdV[i] = linear_fit(current[range])[1]  # store the slope
    end
    return dIdV
end

# ╔═╡ 93924403-fbae-41a1-addd-0ef35c399ab4
@bind reset Button("Uncheck all boxes below") 

# ╔═╡ bebf269e-fd93-4602-8c26-3b708e8c9c99
md"""## Import raw *I(V)* spectra in *mV* & *pA*
*e.g.* from: \
C:\Users\PC\OneDrive - UPJŠ\Dokumenty\data\DQ01\1Q1H\Ondro 1Q1H-stm24-teplotne+poloveZavislosti\Vzorka zhora\Teplotne zavislosti\tempdep\_2020-09-09\_4 \
or \
C:\Users\tomas\OneDrive - UPJŠ\Dokumenty\data\DQ01\1Q1H\Ondro 1Q1H-stm24-teplotne+poloveZavislosti\Vzorka zhora\Teplotne zavislosti\tempdep\_2020-09-09\_4 \
\
(Check filename to match. Watch for path and underscores!)
"""

# ╔═╡ cd219bb0-425e-41e4-bbf4-cf0bd4ed5914
begin
	folder = @bind folder TextField(70)
	md"""Paste folder path: $folder"""
end

# ╔═╡ 9944e204-ff10-4d32-bcce-d25ba099908b
println(folder)

# ╔═╡ 6eee3717-9994-4069-8063-57fa70617e24
begin
	controller = @bind controller Select([1 => "WSxM", 2 => "Nanonis", 3 => "raw data"])
	md"Select file type: $(controller)"
end

# ╔═╡ eb63dc7f-3b64-430f-82e4-f9715091acbd
md"""Select among:
1. WSxM .cur file

under construction:
2. Nanonis .dat file
3. raw data file - specify header, delimiter, etc. below
Data is automatically sorted from lowest to highest bias. \
Decimal comma is automatically replaced with dot for "raw data file".
"""

# ╔═╡ 16523253-6d02-4f0e-9985-14e6423ced42
begin
	if controller == 1
		bwd = @bind bwd Select([0 => "forward", 2 => "backward"])
		md"Select sweep direction: $(bwd)"
		else
		bwd = 0
		md"loading first two columns of the selected file as bias & current"
	end
end

# ╔═╡ 1e9f4924-374f-41a9-8921-1f1906d0172f
begin
	# Initialize results arrays
	temps_raw = Float64[] #vector of temperatures in K
	nums_raw = Int[] #vector of file numbers
	currents = []
	voltages = []

	# WSxM files
	if controller == 1
	# Get all matching files
	files = glob("*Avg.iv.cur", folder)

	
		for file in files
		    # Extract temperature and number from filename
				
			m = match(r"_(\d+)mK_.*?_(\d+)_Current", basename(file)) 
			#m = match(r"-(\d+)mK_(\d+)_Current", basename(file)) 
			# basename(file) returns just the filename part of a full file path — it removes the directory path
		    
				if m === nothing
		        @warn "Filename doesn't match expected pattern: $file"
		        continue
		    	end
			push!(temps_raw, parse(Float64, m.captures[1]) / 1000.0)  # Convert mK to K
		    push!(nums_raw, parse(Int, m.captures[2]))
			
			# read data from file
			# Read file lines
		    lines = readlines(file)
		
		    # Find header end
		    hendIV = findfirst(contains("[Header end]"), lines)
			    if hendIV === nothing
			        @warn "Header end not found in file: $file"
			        continue
			    end
		
		    # Read data after header
		    data = readdlm(file, ' ', skipstart=hendIV)
		
		    unsorted_bias = data[:, 1+bwd].*10^3  # Voltage in mV
		    unsorted_cur = data[:, 2+bwd].*10^3  # Current in pA
			
			# sort bias from negative to positive 	
			# Get the permutation of indices that sorts the bias
			perm = sortperm(unsorted_bias)
			# sort by the bias and remove NaN elements
			sort_bias = unsorted_bias[perm]
			sort_cur = unsorted_cur[perm]
		
			#remove NaN
			# find all indexes of Numbers in current
			# !!! might change the voltage step and mess with differentiation!!!
			ind_N = findall(!isnan, sort_cur)
			bias = sort_bias[ind_N]
			cur = sort_cur[ind_N]
			
		    # Store result
		    push!(currents, [cur])
			push!(voltages, [bias])
		end
	end
	if controller == 2
		@warn "⚠️ Nanonis files import not ready, yet."
	end
	if controller == 3
		@warn "⚠️ Raw files import not ready, yet."
	end
	
	#order data by temperature
	order = sortperm(temps_raw)
	temps = temps_raw[order]
	nums = nums_raw[order]
	
	# Convert to matrix if all rows are same length
	I_matrix = hcat(currents...)[order]
	U_matrix = hcat(voltages...)[order]
	
	# check each column of voltages from the second onward against the first
	if all(vec -> vec == U_matrix[1], U_matrix)
	    @info "✅ All bias vectors are identical."
	else
	    @warn "⚠️ Bias vectors are not the same!"
	end
end

# ╔═╡ 15fce7af-32fe-418f-9f94-185e56777faa
plot(U_matrix[1], I_matrix[:], xlabel = "Bias voltage [mV]", ylabel = "Tunneling current [pA]", title = "$(length(temps)) curves", label = string.(temps'))

# ╔═╡ 3fdc633d-f0e6-47cb-b35c-3337710fadbd
md"""## Process lowest temperature spectrum"""

# ╔═╡ 5d830960-1b18-4343-a52f-8e5eec9233b7
begin
	cur = I_matrix[1]
	raw_bias = U_matrix[1]
	plot(raw_bias, cur, xlabel = "Bias voltage [mV]", ylabel = "Tunneling current [pA]", title = "$(temps[1]) K", legend = false)
end

# ╔═╡ 5c9d391a-c80d-4dd1-961d-3a0d309664f7
md"### Differentiate tunneling current
no bias voltage! \
[Julia documentation](https://www.jlhub.com/julia/manual/en/function/diff)"

# ╔═╡ 54dbb287-d5bf-474f-874e-3e384a98aa93
md"""#### Interval derivative
calculates the slope of linear fit to an interval around each point
- 0 corresponds no derivative
- 1 corresponds to *diff()*
"""

# ╔═╡ 491ce86b-0cfd-406c-bcac-d0c25a6ea917
begin
	points = @bind points Slider(0:500, 50, true)
	md"Points for linear regression: $(points)"
end

# ╔═╡ d91352cd-c72b-4863-99a0-5f48a0672eda
begin
	#create a vector of constant intervals for each points
	pts = ones(Int, length(cur)).*points
	#plot derivative
	plot(pts_dif(cur, pts), legend = false)
end

# ╔═╡ 94ec2390-2c02-42b9-978f-5374716dba34
md"### Bias offset"

# ╔═╡ f321657b-cfbb-46fd-b679-c5e9a2bdfd57
md"""#### Find bias offset
Either move to center using slider or calculate average value of raw peak positions.
"""

# ╔═╡ 12235fc5-3ffa-44fe-86f1-d3d4ac867734
begin
	left_peak = @bind left_peak NumberField(-100.0:0.01:100.0, default = 0)
	right_peak = @bind right_peak NumberField(-100.0:0.01:100.0, default = 0)
	md"left peak: - $(left_peak) mV; right peak: + $(right_peak) mV"
end

# ╔═╡ 635a1125-6df8-4d34-b112-d6969a415652
md"Bias offset from peak positions is $(round((-left_peak + right_peak)/2, digits = 4)) mV."

# ╔═╡ ad5b9aa8-7117-44aa-a768-f3106a1da4b7
begin
	bioff = @bind bioff Slider(-1.0:0.001:1.0, 0.0, true)
	md"Bias offset: $(bioff)"
end

# ╔═╡ 9b7ddf69-97d1-4ccf-bb86-f2aac6aaf4c5
begin
	plot(raw_bias, pts_dif(cur, pts), xlabel = "Bias voltage [mV]", ylabel = "dI/dV [a.u.]", title = "Bias offset: $(bioff) mV.", label = false)
	vline!([-left_peak], label = "left peak")
	vline!([right_peak], label = "right peak")
	vline!([bioff], label = "bias offset")
end

# ╔═╡ 1a09db4e-bdb3-445e-ad11-da2a1db93aa7
md"""#### Remove bias offset
Input bias offset: $(@bind bias_offset NumberField(-100.0:0.01:100.0, default = 0.03)) mV	
"""

# ╔═╡ 2aedf711-745e-4330-b954-0a3592583eb2
begin
	bias = raw_bias .- bias_offset
	plot(bias, pts_dif(cur, pts), xlabel = "Bias voltage [mV]", ylabel = "dI/dV [a.u.]", title = "bias offset removed", label = false)
end

# ╔═╡ 6b7510e8-f34e-43b9-9b1f-71b9d43a0fd5
begin
	#zero bias index...abs.(vector) computes the absolute value of each element in the vector, argmin returns the index of the smallest value in the resulting array, which corresponds to the value closest to 0 in the original vector
	zbi_raw = argmin(abs.(bias))
	
	# lower boundary index and value from center
	lim_ind = min(zbi_raw, length(bias)-zbi_raw)
	lim_val = bias[zbi_raw+lim_ind]
	
	md" Bias voltage closest to 0 is $(bias[zbi_raw]) \
	@ index $zbi_raw"

end

# ╔═╡ 7d91f4cb-512e-41c4-848d-56c85d8a3cfb
md"""#### Magnificent 7
- calculates the slope of linear fit for each point using 3 regions of different values connected by [Fermi-Dirac distribution](https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics#Distribution_of_particles_over_energy) of bins (*i.e.* 7 parameters):``f=\left(\frac{max_{in} - min}{\exp\left(\frac{\left|index - step_{in}\right|}{width_{in}} + 1\right)} \right) + \left( (\text{max} - \text{min}) \left( 1 - \frac{1}{\exp\left(\frac{\left|index - step_{out}\right|}{width_{out}}\right) + 1} \right) + min \right)``
- points can be switched to mV (comment in code)
"""

# ╔═╡ 569dd4e4-ea5b-4ccb-8841-0c530a5173f4
begin
	step = @bind step Slider(0:lim_ind, 900, true)
	width = @bind width Slider(1:500, 10, true)
	max_points = @bind max_points Slider(1:1000, 200, true)
	min_points = @bind min_points Slider(1:1000, 3, true)
	stepin = @bind stepin Slider(0:lim_ind, 100, true)
	widthin = @bind widthin Slider(1:500, 10, true)
	maxin = @bind maxin Slider(1:1000, 3, true)
	md"Outside gap: \
	- ``step_{out}`` position: $(step) points \
	- step ``width_{out}``: $(width) points \
	- ``max``imum: $(max_points) points \
	- ``min``imum: $(min_points) points \
	\
	Inside gap: \
	- ``step_{in}`` position: $(stepin) points \
	- ``width_{in}``: $(widthin) points \
	- ``max_{in}``: $(maxin) points \
	"
end

# ╔═╡ efc0e0c1-af66-4593-a97d-1df83af452a2
begin
	# calculates points for pts_dif for 3 regions connected by Fermi-Dirac distribution of bins
	
	# can be switched to voltage parameters - change index to bias, and in console change values of step, stepin, width & widthin to mV, and channge lim_ind to lim_val
	
	index = (1:length(bias)).-zbi_raw # center value
	
	pts7 = @. round(Int, 
    ((maxin - min_points) / (exp((abs(index) - stepin) / widthin) + 1)) +
    ((max_points - min_points) * (1 - 1 / (exp((abs(index) - step) / width) + 1)) + 	min_points))

	plot(bias, pts_dif(cur, pts7).*2000, label = "", xlabel = "Bias voltage [mV]", ylabel = "Differential conductance [a.u.]", legend=:topright)
	plot!(bias, pts7,  label = "")
	vline!([bias[zbi_raw + step], bias[zbi_raw - step+1]], label = "step outside gap")
	vline!([bias[zbi_raw + stepin], bias[zbi_raw - stepin+1]], label = "step inside gap")
end

# ╔═╡ 200ad45d-8923-46e3-abaa-f72b644fc7cd
md"""#### Select data for further processing
"""

# ╔═╡ 7246d556-39a9-4bc8-aec0-a4da6562ee3f
begin
	plot(bias, pts_dif(cur, pts), xlabel = "Bias voltage [mV]", ylabel = "Differential conductance [a.u.]", label = "Interval derivative")
	plot!(bias, pts_dif(cur, pts7), label = "Magnificent 7", alpha = 0.9)
end

# ╔═╡ 77a12ae9-a85d-4432-9263-9006fb107bba
begin
	dif_met = @bind dif_met Select([1 => "Interval ", 2 => "Magnificent 7"])
	md"Select differentiation method: $(dif_met)"
end

# ╔═╡ a23ba2a6-6213-4837-8179-df6498a299f0
begin
	if dif_met == 1
		didv = pts_dif(cur, pts)
	end
	if dif_met == 2
		didv = pts_dif(cur, pts7)
	end
	if dif_met == 3
		didv = finder
	end
	plot(bias, didv, xlabel = "Bias voltage [mV]", ylabel = "Differential conductance [a.u.]", title = "Selected data", label = false)
end

# ╔═╡ 69c16299-6b9d-449d-ac19-4ea71487ed43
md"""### Symmetric bias crop
- __*zbi*__ - zero bias index - is the index of the bias closest to zero
- we want *zbi* to be in the middle, *i.e.* the vector of non-positive bias values should be one element longer than the vector of positive values
- length of non-positive bias vector *__bias[1:zbi]__* can be i) larger than, ii) smaller than, or iii) equal to vector of positive bias *__bias[zbi+1:end]__*
- i)  ``zbi > length(bias)-zbi \implies bias[2\times zbi- end:end]``
- ii) `` zbi < length(bias)-zbi \implies bias[1:2\times zbi - 1]``
- iii) ``2 \times zbi = length(bias) \implies bias[1:end-1] = bias[1:2\times zbi - 1]``
"""

# ╔═╡ 80a4da83-9bf3-4dca-b520-e35be37ac50d
begin
	#calculate index range symmetric around zero bias
	lastind = lastindex(bias) # index of "end"
	if 2*zbi_raw > length(bias)
		cut_sym = 2*zbi_raw-lastind:lastind
	else
		cut_sym = 1:2*zbi_raw-1
	end

	mid_ind = 	fld(length(cut_sym) + 1, 2) #computes the floor division to get the middle index
	# add points to cut
	cut_add = @bind cut_add Slider(0:mid_ind-1, 0, true)
	md"Additional points to cut: $(cut_add)"
end

# ╔═╡ 9d70f2f1-7188-453b-968e-0c6ca886bbda
begin
	# crop data symmetrically around zbi
	cut = cut_sym[cut_add+1:end - cut_add] #vector of chosen indexes
	
	#cut data
	cut_bias = bias[cut]
	cut_cur = cur[cut]
	cut_cond = didv[cut]
	lo_cut = cut[1]
	hi_cut = cut[end]
	
	plot(bias, didv, label = "processed data", xlabel = "Bias voltage [mV]", ylabel = "Normalized dI/dV [a.u.]")
	#plot range
	vline!([bias[lo_cut],bias[hi_cut]], legend=:bottomright, label = "bias cutoff")
end

# ╔═╡ e226ec3d-8e8b-4243-9bad-ab0d092a7c6b
md"""### Smooth
__LOESS__ [locally estimated scatterplot smoothing](https://en.wikipedia.org/wiki/Local_regression)
"""

# ╔═╡ beb37c65-e341-417b-b5b8-091142bcf475
begin
	span = @bind span Slider(0.0:0.01:1.0, 0.1, true)
	md"Span: $(span)"
end

# ╔═╡ 786d617f-4d8b-4adc-a079-5382508ca353
begin
	model = loess(cut_bias, cut_cond, span = span)
	smooth_cond = predict(model, cut_bias)
	#us = range(extrema(bias)...; step = 0.001)
	#vs = predict(model, us)
	
	plot(cut_bias, cut_cond, label = "processed data", xlabel = "Bias voltage [mV]", ylabel = "Differential conductance [a.u.]", title = "Selected data")
	plot!(cut_bias, smooth_cond, label = "LOESS", lw = 4, alpha = 0.5)
end

# ╔═╡ 8c7dc91b-c6d3-4963-992d-7cf351c514e7
md"### Normalize
to normal metal conductance using linear fit to points out of gap "

# ╔═╡ b76a33c3-818f-4e76-8170-c2195ed63fe0
#normalize to normal metal conductance using linear fit to points out of gap
function normetal(data, percent)
    length(data) * percent / 100 |>
    x -> round(Int, x) |>
    pts -> vcat(data[1:pts], data[end-pts+1:end]) |>
    linear_fit |>
    fit -> fit[2] .+ fit[1] .* collect(eachindex(data))
end

# ╔═╡ 319e762d-6917-49b1-b6a9-8c29e24ab339
begin
	# set percentage of points from each side to fit a line
	metal_perc = @bind metal_perc Slider(0:50, 25, true)
	md"""Normal conductance points from each side: $(metal_perc) %"""
end

# ╔═╡ c6333d05-7d51-497d-9bf8-5557769737b5
begin
	#zero bias is in the middle
	#zbi = fld(length(cut_bias) + 1, 2)
	
	fit_metal = normetal(smooth_cond, metal_perc)

	plot(smooth_cond, ylabel = "Differential conductance [a.u.]", label = "data", legend = :bottomright)
	plot!(fit_metal, label = "linear fit")
	metal_pts = round(Int, length(cut_bias) * metal_perc / 100)
	vline!([metal_pts length(cut_bias)-metal_pts], label = "$(metal_perc)%")
end

# ╔═╡ c90af035-9f38-4928-86d1-8963755631a9
begin
	norm_cond = smooth_cond./fit_metal
	plot(cut_bias, norm_cond, label = "processed data", xlabel = "Bias voltage [mV]", ylabel = "Normalized dI/dV [a.u.]", title = "final figure")
	hline!([1], ls = :dash, label = "")
end

# ╔═╡ c1ac85d0-29d6-45dc-8ff9-6e15d06e5ee8
md"#### Save processed spectrum"

# ╔═╡ 55ab325f-6fb1-4494-af20-0c4edb63dd1e

begin
	DownloadButton([cut_bias norm_cond], "$(split(basename(files[1]), ".")[1]).txt") # download generated curve to a file with the name of input file but extension txt
end

# ╔═╡ 9c919d2b-8278-49ad-8700-7dc845b9794b
md"#### Save processessing parameters"

# ╔═╡ 1833ebc2-df5b-4f4a-bd4f-a2524dab6274
begin
	DownloadButton([points, bias_offset, span, metal_perc, lo_cut, hi_cut], "$(split(basename(files[1]), ".")[1])_par.txt") # download generated curve to a file with the name of input file but extension txt
end

# ╔═╡ 45eabd6f-2030-4323-9527-0dbba24f6fe2
md"### Dynes DOS fit
[Dynes et al., PRL '78](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.41.1509)\
\
``N(E,\Gamma)=\Re\left[\frac{E - i\Gamma}{\sqrt{(E-i\Gamma)^2-\Delta^2}}\right]``"

# ╔═╡ 77c591ea-2c3c-4951-b0bf-97c7828a531c
md"#### Fermi smearing at finite temperatures
[wiki] (https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics#Distribution_of_particles_over_energy) & code from Jose Gabriel\
using the [exponential form] (https://en.wikipedia.org/wiki/Hyperbolic_functions#Exponential_definitions) of [the Fermi function derivative](https://lampz.tugraz.at/~hadley/ss1/materials/thermo/gp/gp/Fermi-function.html)
#### ``f(E)=\frac{1}{e^{\frac{(E-E_F)}{k_BT}}+1};~f'(E)\approx\frac{1}{4k_BTcosh^2\left(\frac{E-E_F}{2k_BT}\right)}``"

# ╔═╡ 8828be87-4df7-42b7-95ad-f84713faa69b
begin
	# Code from Jose Gabriel Rodrigo
	# check "C:\Users\PC\OneDrive - UPJŠ\Dokumenty\code\Matlab\IV\BCS & Dynes\Fermi fft"
	function fermifilter(E,t)
		E=E.*0.001
		kb=8.617343*10^(-5)
		#find center
		mid = div(length(E), 2)
		1.0 ./ (4 * kb * t  .* (cosh.((E .- E[mid]) ./ (2 * kb * t))).^2)
	end
	# convolution of Fermi and DOS usinf fft instead of integral
	# check "C:\Users\PC\OneDrive - UPJŠ\Dokumenty\code\Matlab\IV\BCS & Dynes\Fermi fft"
	function convol(energy,N,t)
		du=energy[2]-energy[1]
		#find center by dividing and discarding decimal values
		mid = div(length(energy), 2)
		# Calculate the array to be shifted and transformed
		array = fermifilter(energy,t)
		# Perform circular shift
		shifted_array = circshift(array, mid + 1)
		# Perform FFT
		g = fft(shifted_array)
		# Perform FFT on N
		fft_N = fft(N)
		# Element-wise multiplication of f and fft_N
		multiplied = g .* fft_N
		# Perform IFFT on the result
		ifft_result = ifft(multiplied)
		# Take the real part and multiply by du
		G = du * real(ifft_result)
		#normalize to leftmost value
		GN = G./G[1]
	end
end

# ╔═╡ 80701395-09d1-4eb6-a788-79eb18715097
begin
	function dynes(u, del, gam)
		N = abs.(real.((u .- 1im*gam)./sqrt.(Complex.((u.- 1im*gam).^2 .-del^2))))
		#normalize to leftmost value
		NN=N./N[1]
	end
	function model_dos(bias, p)
		 #calculate DOS with Dynes
	 	 zero = dynes(bias, p[1], p[2])
		 # temperature
		 result = convol(bias, zero, p[3])
	 end
end

# ╔═╡ c15bafa3-921b-4b64-805c-362ae9d4b8de
begin
	t = @bind t NumberField(0.01:0.01:10, default = temps[1])
	terror = @bind terror NumberField(0:1:100, default = 0)
	del_dynes = @bind del_dynes Slider(0.00:0.01:10.00, 0.12, true)
	gam_dynes = @bind gam_dynes Slider(0.001:0.01:0.5, 0.011, true)
	md"Fitting parameters: \
	``\Delta``: $(del_dynes) meV \
	``\Gamma``: $(gam_dynes) meV \
	``T=`` $(t)K ``\pm`` $(terror)%"
end

# ╔═╡ 8c8a3cd1-7fc3-4700-9039-80e737d87b33
p0_dos=[del_dynes, gam_dynes, t]

# ╔═╡ cf7a6448-ad8e-4ec7-b3bc-25b261fa4af1
lo_dos = [0.0, 0.0, t - t*terror/100]

# ╔═╡ eab02cdc-8f97-4c1c-b346-ce52c6e4df1f
up_dos = [2.0, 1.0, t + t*terror/100]

# ╔═╡ 6f79f9ef-d61c-4a7d-ac27-40db625bc1b5
begin
	reset
	@bind fit1 CheckBox()
	md"#### Fit the spectrum? $(@bind fit1 CheckBox())"
end

# ╔═╡ 2c255601-3a0d-4e1a-9118-9f92c97dc4ce
begin
	if fit1 == true
		plot(cut_bias, norm_cond,
				title = "Dynes DOS guess",
				label="experiment", 	
				ylabel = "Normalized dI/dV [a.u.]", 
				xlabel = "Bias [meV]")
		plot!(cut_bias, model_dos(cut_bias, p0_dos), label="initial guess")
		else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ 85892a58-8847-4b07-9c24-fc65755c21d2
begin
	if fit1 == true
		fitdos = curve_fit(model_dos, cut_bias, norm_cond, p0_dos, lower = lo_dos, upper = up_dos)
		plot(cut_bias, norm_cond,
				title = "Dynes DOS fit",
				label="experiment", 	
				ylabel = "Normalized dI/dV [a.u.]", 
				xlabel = "Bias [meV]")
		plot!(cut_bias, model_dos(cut_bias, fitdos.param),label="LsqFit")
		else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ d7e755bb-b254-4759-b2bc-e8285337e705
begin
	if fit1 == true
		md"### Fit result: ``\Delta \approx`` $(round(fitdos.param[1], digits = 2)) meV, ``\Gamma \approx`` $(fitdos.param[2])  meV, ``T \approx`` $(round(fitdos.param[3], digits = 2)) K"
		
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ 9e13f528-8c03-4af2-8255-bdb73396f195
md"""## Process spectra en bloc"""

# ╔═╡ ef71f4aa-1bb9-43dc-aaa9-359ee6f8573d
begin
	reset
	@bind proc CheckBox()
	md"#### Process all spectra? $(@bind proc CheckBox())"
end

# ╔═╡ 60d10afa-edc7-485a-be14-9a1c9142e4db
begin
	if proc == true
		norm_t = I_matrix |>
		# differentiate each vector of matrix
			x -> map(v -> pts_dif(v, pts), x) |>
		# cut each vector of matrix
			x -> map(v -> v[cut], x) |>
		# smooth spectra		
			x -> map(v -> predict(loess(cut_bias, v, span=span), cut_bias), x) |>
		# normalize each curve to its linear fit
			x -> map(v -> v ./ normetal(v, metal_perc), x)
				
		plot(cut_bias, norm_t[:], label = string.(temps'), xlabel = "Bias voltage [mV]", ylabel = "Normalized dI/dV [a.u.]", title = "Processed spectra", legend=:bottomleft)
		
		else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ ea44c5f7-b55b-4b09-892a-d34659d4619d
begin
	if proc == true
		#vec(mat) flattens the 1×59 matrix into a vector of vectors, reduce(vcat, ...) stacks them into a 2D matrix.
		red_mat = reduce(vcat, [v' for v in vec(norm_t)])
		
		# Plot as a surface
		surface(cut_bias, temps, red_mat, title="Surface Plot", xlabel = "Bias voltage [mV]", ylabel = "temperature [K]", zlabel="Normalized dI/dV [a.u.]", color=:turbo, colorbar=false, zlims = (0, :auto))
		
		else
		println("Waiting for checkbox...")
	end
end


# ╔═╡ 648ad72d-ff6e-4390-8d97-c95b020c9d26
md"""## Fit Dynes DOS to spectra en bloc"""

# ╔═╡ 2ab86e14-ae79-4996-a751-e3112d2bae1a
# BCS delta(t)
function deltanh(t, par)
    result = nlsolve(d -> tanh.(d / (t./par[1])) - d, [0.5])
    return result.zero[1].*par[2]
end

# ╔═╡ b40516b6-70f3-4dd6-8150-e0c0f18cfed9
function fitanh(temps, par)
    return map(t -> deltanh(t, par), temps)
end

# ╔═╡ 125e2aea-9f87-4e57-b2a5-8191770c258e
begin
	del0_guess = @bind del0_guess NumberField(0.01:0.01:10, default = del_dynes)
	tc_guess = @bind tc_guess NumberField(0:0.01:10, default = del_dynes*6.6)
	md"Initial estimate of \
	``\Delta(0) = `` $(del0_guess)meV \
	``T_c`` = $(tc_guess)K \
	for processing."
end

# ╔═╡ ff5cd9f7-6e4d-4e34-8386-2c09e4015572
begin
	reset
	@bind fitall CheckBox()
	md"#### Fit all spectra? $(@bind fitall CheckBox())"
end

# ╔═╡ 21284cd3-54d6-4468-818d-c029c29a9909
begin
	if fitall == true
	# initial guess of delta(T)
	del_ini = fitanh(temps, [tc_guess, del0_guess])
	
	# fit all curves in norm_t with Dynes and correpsonding delta and temperature
	fit_all = [curve_fit(model_dos, cut_bias, norm_t[i], [del_ini[i], gam_dynes, temps[i]], lower=[lo_dos[1], lo_dos[2], temps[i]-temps[i]*terror/100], upper=[up_dos[1], up_dos[2], temps[i]+temps[i]*terror/100]).param for i in eachindex(norm_t)]
	
	# extract deltas
	deltas_fit = [v[1] for v in fit_all]

	# extract gammas
	gammas_fit = [v[2] for v in fit_all]
	
	# fit delta vs T
	delta_t = curve_fit(fitanh, temps, deltas_fit, [0.8, 0.1]).param

	# plot Dynes fit results
	scatter(temps, deltas_fit, label = "Δ from Dynes fit ", xlabel = "T [K]", ylabel = "Δ [meV]")

	# plot delta(T) fit
	tempfit = range(0.0, delta_t[1], 1001)
	plot!(tempfit, fitanh(tempfit, delta_t), label = "Δ(T) fit")
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ abb50450-3960-4a88-b36c-a9399244e1b3
begin
	if fitall == true
	md"#### Fit result: ``T_c \approx`` $(round(delta_t[1], digits = 2)) K, ``\Delta_0 \approx`` $(round(delta_t[2], digits = 2)) meV, ``\frac{2\Delta}{k_BT_c} \approx`` $(round(2*11.6*delta_t[2]/delta_t[1], digits = 2))"
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ 30612e72-404e-4def-baff-d4ed522b328e
md"""### Individual Dynes fits"""

# ╔═╡ c2280603-a1f6-4465-90d0-4103ef2447c9
begin
	if fitall == true
	scatter(deltas_fit, label = "Δ from Dynes fit ", xlabel = "spectrum nr.", ylabel = "Δ [meV]")
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ ab015edf-6deb-4538-9dc0-f6829af06ae0
begin
	if fitall == true
	scatter(gammas_fit, label = "Γ from Dynes fit ", xlabel = "spectrum nr.", ylabel = "Γ [meV]")
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ aef936df-9d8c-4f0c-bab3-447a55c24841
begin
	if fitall == true
	dynes_nr = @bind dynes_nr NumberField(1:length(temps), default = 1)
	md"Plot curve nr. $(dynes_nr)"
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ 9c0eaa78-3368-4d2f-9995-228ae33eeba5
begin
	if fitall == true
	plot(cut_bias, norm_t[dynes_nr],
				title = "Spectrum measured at $(temps[dynes_nr]) K",
				label="experiment", 	
				ylabel = "Normalized dI/dV [a.u.]", 
				xlabel = "Bias [meV]")
		plot!(cut_bias, model_dos(cut_bias, fit_all[dynes_nr]),label="LsqFit")
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ dabf93e1-8d86-4c01-bc55-5c13537ead30
begin
	if fitall == true
	md"#### Fit #$(dynes_nr): ``\Delta \approx`` $(round(fit_all[dynes_nr][1], digits = 2)) meV, ``\Gamma \approx`` $(fit_all[dynes_nr][2])  meV, ``T \approx`` $(round(fit_all[dynes_nr][3], digits = 2)) K"
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ 2085e01a-6dcd-4300-8b1e-abd50de0d3bd
begin
	reset
	@bind remove_out CheckBox()
	md"### Remove outliers & zeros? $(@bind remove_out CheckBox())
	input indexes in cell below..."
end

# ╔═╡ e1d48916-0931-4487-aca5-744ce230cc3c
begin
	if remove_out == true
		
	# Indices to remove
	last = length(deltas_fit)
	remove_indices = [33, 38, 40:42..., 44:last...]
	
	# Keep only the elements not in remove_indices
	deltas_rem = deltas_fit[setdiff(1:length(deltas_fit), remove_indices)]
	temps_rem = temps[setdiff(1:length(deltas_fit), remove_indices)]
	red_mat_rem = red_mat[setdiff(1:size(red_mat, 1), remove_indices),:]

	# plot remaining spectra
	surface(cut_bias, temps_rem, red_mat_rem, title="Surface Plot", xlabel = "Bias voltage [mV]", ylabel = "temperature [K]", zlabel="Normalized dI/dV [a.u.]", color=:turbo, colorbar=false, zlims = (0, :auto))
	else
		println("Waiting for checkbox...")
	end
end


# ╔═╡ 2da46c7f-2b61-428e-8a52-9c952cee9078
println("Indices removed:", remove_indices)

# ╔═╡ b85701ab-8bf3-42f3-a918-778c610ddc30
begin
	if remove_out == true
	# fit delta vs T
	delta_t_rem = curve_fit(fitanh, temps_rem, deltas_rem, [0.8, 0.1]).param

	# plot Dynes fit results
	scatter(temps_rem, deltas_rem, label = "Δ from Dynes fit ", xlabel = "T [K]", ylabel = "Δ [meV]")

	# plot delta(T) fit
	tempfit_rem = range(0.0, delta_t_rem[1], 1001)
	plot!(tempfit_rem, fitanh(tempfit_rem, delta_t_rem), label = "Δ(T) fit")
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ d3bd25a2-f249-4d1b-9415-6d51f051fdbb
begin
	if remove_out == true
	md"#### Fit result: ``T_c \approx`` $(round(delta_t_rem[1], digits = 2)) K, ``\Delta_0 \approx`` $(round(delta_t_rem[2], digits = 2)) meV, ``\frac{2\Delta}{k_BT_c} \approx`` $(round(2*11.6*delta_t_rem[2]/delta_t_rem[1], digits = 2))"
	else
		println("Waiting for checkbox...")
	end
end

# ╔═╡ 4a689046-9125-4b8a-bee4-b16029367d10
md"""#### My estimate: ``T_c \approx`` 1.25 K, ``\Delta_0 \approx`` 0.23 meV, ``\frac{2\Delta}{k_BT_c} \approx`` $(round(2*11.6*0.23/1.25, digits = 2))
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
Loess = "4345ca2d-374a-55d4-8d30-97f9976e7612"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
FFTW = "~1.8.0"
Glob = "~1.3.1"
Loess = "~0.6.4"
LsqFit = "~0.15.0"
NLsolve = "~4.5.1"
Plots = "~1.40.9"
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "122b4d1fb6196efff5d8d460d0bfc584dbcb33f3"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "d80af0733c99ea80575f612813fa6aa71022d33a"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.0"

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "c785dfb1b3bfddd1da557e861b919819b82bbe5b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.27.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

    [deps.Distances.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "d7477ecdafb813ddee2ae727afa94e9dcb5f3fb0"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.112"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc5231d52eb1771251fbd37171dbc408bcc8a1b6"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "84e3a47db33be7248daa6274b287507dd6ff84e8"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.26.2"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "fa8e19f44de37e225aa0f1695bc223b05ed51fb4"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.3+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee28ddcd5517d54e417182fec3886e7412d3926f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.8"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f31929b9e67066bee48eec8b03c0df47d31a74b3"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.8+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b36c7e110080ae48fdef61b0c31e6b17ada23b33"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.2+0"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "01979f9b37367603e2848ea225918a3b3861b606"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ae350b8225575cc3ea385d4131c81594f86dfe4f"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.12"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "71b48d857e86bf7a1838c4736545699974ce79a2"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.9"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6ce1e19f3aec9b59186bdf06cdf3c4fc5f5f3e6"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.50.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "61dfdba58e585066d8bce214c5a51eaa0539f269"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "84eef7acd508ee5b3e956a2ae51b05024181dee0"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "edbf5309f9ddf1cab25afc344b1e8150b7c832f9"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.2+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "f749e7351f120b3566e5923fefdf8e52ba5ec7f9"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.4"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "Printf", "StatsAPI"]
git-tree-sha1 = "40acc20cfb253cf061c1a2a2ea28de85235eeee1"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.15.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "dae01f8c2e069a683d3a6e17bbae5070ab94786f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.9"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "a2fccc6559132927d4c5dc183e3e01048c6dcbd6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7d1671acbe47ac88e981868a078bd6b4e27c5191"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.42+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "15e637a697345f6743674f1322beefbc5dcd5cfc"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2b0e27d52ec9d8d483e2ca0b72b3cb1a8df5c27a"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+1"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "02054ee01980c90297412e4c809c8694d7323af3"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+1"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee57a273563e273f0f53275101cd41a8153517a"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+1"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b9ead2d2bdb27330545eb14234a2e300da61232e"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6e50f145003024df4f5cb96c7fce79466741d601"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.56.3+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─003643f0-a1df-11ef-02bc-a3bb1aac6d2e
# ╟─3f2957da-93d6-43bc-b111-0c110ea5647a
# ╟─f03eb307-16cc-428f-924e-e057181428d4
# ╟─5290a6ce-89ae-4ed9-b02f-f1e6033e07dd
# ╟─833c981f-d250-4194-8454-ad10be888005
# ╟─432f4d98-109c-49e6-8be7-eba090783fc7
# ╟─f4d390b1-4da9-41a7-bfb1-1c1b3d0d89a1
# ╟─93924403-fbae-41a1-addd-0ef35c399ab4
# ╟─bebf269e-fd93-4602-8c26-3b708e8c9c99
# ╟─cd219bb0-425e-41e4-bbf4-cf0bd4ed5914
# ╟─9944e204-ff10-4d32-bcce-d25ba099908b
# ╟─6eee3717-9994-4069-8063-57fa70617e24
# ╟─eb63dc7f-3b64-430f-82e4-f9715091acbd
# ╟─16523253-6d02-4f0e-9985-14e6423ced42
# ╟─1e9f4924-374f-41a9-8921-1f1906d0172f
# ╟─15fce7af-32fe-418f-9f94-185e56777faa
# ╟─3fdc633d-f0e6-47cb-b35c-3337710fadbd
# ╟─5d830960-1b18-4343-a52f-8e5eec9233b7
# ╟─5c9d391a-c80d-4dd1-961d-3a0d309664f7
# ╟─54dbb287-d5bf-474f-874e-3e384a98aa93
# ╟─491ce86b-0cfd-406c-bcac-d0c25a6ea917
# ╟─d91352cd-c72b-4863-99a0-5f48a0672eda
# ╟─94ec2390-2c02-42b9-978f-5374716dba34
# ╟─f321657b-cfbb-46fd-b679-c5e9a2bdfd57
# ╟─12235fc5-3ffa-44fe-86f1-d3d4ac867734
# ╟─635a1125-6df8-4d34-b112-d6969a415652
# ╟─ad5b9aa8-7117-44aa-a768-f3106a1da4b7
# ╟─9b7ddf69-97d1-4ccf-bb86-f2aac6aaf4c5
# ╟─1a09db4e-bdb3-445e-ad11-da2a1db93aa7
# ╟─2aedf711-745e-4330-b954-0a3592583eb2
# ╟─6b7510e8-f34e-43b9-9b1f-71b9d43a0fd5
# ╟─7d91f4cb-512e-41c4-848d-56c85d8a3cfb
# ╟─569dd4e4-ea5b-4ccb-8841-0c530a5173f4
# ╟─efc0e0c1-af66-4593-a97d-1df83af452a2
# ╟─200ad45d-8923-46e3-abaa-f72b644fc7cd
# ╟─7246d556-39a9-4bc8-aec0-a4da6562ee3f
# ╟─77a12ae9-a85d-4432-9263-9006fb107bba
# ╟─a23ba2a6-6213-4837-8179-df6498a299f0
# ╟─69c16299-6b9d-449d-ac19-4ea71487ed43
# ╟─80a4da83-9bf3-4dca-b520-e35be37ac50d
# ╟─9d70f2f1-7188-453b-968e-0c6ca886bbda
# ╟─e226ec3d-8e8b-4243-9bad-ab0d092a7c6b
# ╟─beb37c65-e341-417b-b5b8-091142bcf475
# ╟─786d617f-4d8b-4adc-a079-5382508ca353
# ╟─8c7dc91b-c6d3-4963-992d-7cf351c514e7
# ╟─b76a33c3-818f-4e76-8170-c2195ed63fe0
# ╟─319e762d-6917-49b1-b6a9-8c29e24ab339
# ╟─c6333d05-7d51-497d-9bf8-5557769737b5
# ╟─c90af035-9f38-4928-86d1-8963755631a9
# ╟─c1ac85d0-29d6-45dc-8ff9-6e15d06e5ee8
# ╟─55ab325f-6fb1-4494-af20-0c4edb63dd1e
# ╟─9c919d2b-8278-49ad-8700-7dc845b9794b
# ╟─1833ebc2-df5b-4f4a-bd4f-a2524dab6274
# ╟─45eabd6f-2030-4323-9527-0dbba24f6fe2
# ╟─77c591ea-2c3c-4951-b0bf-97c7828a531c
# ╟─8828be87-4df7-42b7-95ad-f84713faa69b
# ╟─80701395-09d1-4eb6-a788-79eb18715097
# ╟─c15bafa3-921b-4b64-805c-362ae9d4b8de
# ╟─8c8a3cd1-7fc3-4700-9039-80e737d87b33
# ╟─cf7a6448-ad8e-4ec7-b3bc-25b261fa4af1
# ╟─eab02cdc-8f97-4c1c-b346-ce52c6e4df1f
# ╟─6f79f9ef-d61c-4a7d-ac27-40db625bc1b5
# ╟─2c255601-3a0d-4e1a-9118-9f92c97dc4ce
# ╟─85892a58-8847-4b07-9c24-fc65755c21d2
# ╟─d7e755bb-b254-4759-b2bc-e8285337e705
# ╟─9e13f528-8c03-4af2-8255-bdb73396f195
# ╟─ef71f4aa-1bb9-43dc-aaa9-359ee6f8573d
# ╟─60d10afa-edc7-485a-be14-9a1c9142e4db
# ╟─ea44c5f7-b55b-4b09-892a-d34659d4619d
# ╟─648ad72d-ff6e-4390-8d97-c95b020c9d26
# ╟─2ab86e14-ae79-4996-a751-e3112d2bae1a
# ╟─b40516b6-70f3-4dd6-8150-e0c0f18cfed9
# ╟─125e2aea-9f87-4e57-b2a5-8191770c258e
# ╟─ff5cd9f7-6e4d-4e34-8386-2c09e4015572
# ╟─21284cd3-54d6-4468-818d-c029c29a9909
# ╟─abb50450-3960-4a88-b36c-a9399244e1b3
# ╟─30612e72-404e-4def-baff-d4ed522b328e
# ╟─c2280603-a1f6-4465-90d0-4103ef2447c9
# ╟─ab015edf-6deb-4538-9dc0-f6829af06ae0
# ╟─aef936df-9d8c-4f0c-bab3-447a55c24841
# ╟─9c0eaa78-3368-4d2f-9995-228ae33eeba5
# ╟─dabf93e1-8d86-4c01-bc55-5c13537ead30
# ╟─2085e01a-6dcd-4300-8b1e-abd50de0d3bd
# ╟─e1d48916-0931-4487-aca5-744ce230cc3c
# ╟─2da46c7f-2b61-428e-8a52-9c952cee9078
# ╟─b85701ab-8bf3-42f3-a918-778c610ddc30
# ╟─d3bd25a2-f249-4d1b-9415-6d51f051fdbb
# ╟─4a689046-9125-4b8a-bee4-b16029367d10
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
