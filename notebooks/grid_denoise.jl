### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ a5540242-d382-4ac5-87b9-c0c45f74cef3
begin
	using PlutoUI, PlutoPlotly, JLD2,  Lowess, NoiseRobustDifferentiation, Statistics, FFTW
	import Plots as gr
	md"`packages`"
end

# ╔═╡ 6a6ea7ad-0a03-464d-8b84-4f642db7cfdd
TableOfContents()

# ╔═╡ 6df8df70-3346-11f1-a692-31d015c5335e
begin
	source = @bind source Select([1 => "WSxM", 2 => "Nanonis"])
md"# Spectral grid
using [PlutoPlotly](https://github.com/JuliaPluto/PlutoPlotly.jl)\
Select file source: $(source) ...under construction...

process:
1. smooth_iv & dif -- switch to Savitzky-Golay
1. compare with lowess & nrd
1. save spectrum to .tsv [bias, di/dv, current]

plots:
- dI/dV <switch> (dI/dV)/(I/V)
- peak plots (renormalized option) -- one range (position, value), two ranges (gap map)
- Matlab plots (...Matlab\Gmaps\\)
- FFT (QPI, ±E antisymmetrization)
- profile & spectra along multiple clicks
- save .jld2 using MultiCheckBox

misc:
- check tunneling current values
- add Nanonis
"
end

# ╔═╡ 89ec1e0d-3811-4ef7-9bbc-4a0f9d8d41b4
md"""## Import data
**drag & drop** grid file: $(@bind gsi FilePicker()) (e.g."...\Pluto\STM\datatest\test.gsi")
"""

# ╔═╡ 76884c57-b7e3-46b2-b095-5fa39a0de81f
md"""
!!! warning "Check units!!!"
	Load data in mV & pA. \
	Process topography in *e.g.* [WSxM](https://wsxm.eu/) or [Gwyddion](https://gwyddion.net/).\
	Data is automatically sorted from lowest to highest bias. \
"""

# ╔═╡ a52ec94e-a22e-4c7f-86af-c9eaf0cffdca
@bind en_bloc Button("Disable en bloc processing")

# ╔═╡ c0f60871-fde9-429a-8389-f014844cd2a2
md"# Process *I(V)* map"

# ╔═╡ 557407a7-fcda-48a6-9a60-2f7626e579dc
md"## Topography"

# ╔═╡ dd53ab84-5693-4c87-bb41-710841cc02aa
@bind unclick Button("unclick")

# ╔═╡ 181f2fab-cc28-4016-b9a9-ce0db79ef2ab
md"## *I(V)* spectrum @ clixel
[1,1] if no click"

# ╔═╡ 94bf84cc-78ab-4612-b673-d12baf9b74e0
md"## Process spectrum"

# ╔═╡ 6d5bd5bb-5e4b-41d4-ab70-9066ce13bddc
md"""### Interval derivative of raw *I(V)*
Calculates the slope of linear fit to an interval around each point.
"""

# ╔═╡ ae073f35-b095-491c-8347-33e9b69764f9
begin
	md"Points for linear regression: $(@bind points Slider(1:50, default = 9, show_value = true))"
end

# ╔═╡ 08296e70-9f4b-4aca-834d-21c9a4e7d1d7
md"""### Fourier transform of ``\frac{dI}{dV}``"""

# ╔═╡ 22107084-f3f7-4d1e-acbb-6ab4f8a87ab3
md"### Bias offset"

# ╔═╡ 1b05e239-7419-4aeb-81ce-9c8e66ad4b51
md"""#### Find bias offset
Either move to center using slider or calculate average value of raw peak positions.
"""

# ╔═╡ 1feba41d-4dfd-4757-abf4-f3397ac48b65
begin
	left_peak = @bind left_peak NumberField(-50.0:0.01:50.0, default = 0)
	right_peak = @bind right_peak NumberField(-50.0:0.01:50.0, default = 0)
	md"left peak: - $(left_peak) mV; right peak: + $(right_peak) mV"
end

# ╔═╡ 1e60e654-0a68-4efd-a870-0d56b65e9731
md"Bias offset from peak positions is $(round((-left_peak + right_peak)/2, digits = 4)) mV."

# ╔═╡ 92e25728-fa53-4910-a779-c739c3eee693
md"Bias offset: $(@bind bioff Slider(-1.0:0.01:1.0, default = 0.0, show_value = true))"

# ╔═╡ 73f16b28-e715-4352-9698-e2679fd664ce
md"""#### Remove bias offset
Input bias offset: $(@bind bias_offset NumberField(-100.0:0.01:100.0, default = 0)) mV	
"""

# ╔═╡ c2547995-2493-4e1e-a956-7ab960ccf5c5
md"""### Smooth *I(V)*
[Lowess.jl](https://github.com/xKDR/Lowess.jl)
[wiki](https://en.wikipedia.org/wiki/Local_regression)
[R](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/lowess.html)

!!! warning "Apply carefully!!!"
"""

# ╔═╡ 9bcc955e-6fcf-49e2-a548-09368e2c6df8
begin
	md"""``f =`` $(@bind lwf Slider(0.0:0.001:0.2, default =  0.03, show_value = true)) ``\boxed{n_{steps} = 3, \delta = 0.0}``"""
end

# ╔═╡ fd48b9a7-f51b-4ee8-b3dd-f3982154414f
md"""### Current offset
Points around the zero bias conductance: $(@bind curoff_pts NumberField(1:500, default = 1))
"""

# ╔═╡ b092b5d4-66ef-48e5-b218-856b4b108d49
md"""### Noise robust differentiation of *I(V)*
[noise robust differentiation](https://adrianhill.de/NoiseRobustDifferentiation.jl/), [Rick Chartrand](https://onlinelibrary.wiley.com/doi/10.5402/2011/164564)
- [Heteroscedastic noise](https://en.wikipedia.org/wiki/Homoscedasticity_and_heteroscedasticity) is suppressed by segment-wise TV differentiation.
"""

# ╔═╡ 20513a32-5bb6-44cd-9f68-869ca54a01ff
begin
	reset
	iter = @bind iter NumberField(1:10000, default = 50)
	md""" iterations = $(iter)\
	`diff_kernel:` $(@bind kernel Select(["abs", "square"], default = "square"))\
	"""
end

# ╔═╡ e2375ff2-7095-4751-80ea-eb3a3d1a6f5b
md"### Normalize
to normal metal conductance using linear fit to points out of gap."

# ╔═╡ 208a0062-fc5b-4494-9ae2-29b02c8c3935
begin
	# set percentage of points from each side to fit a line
	metal_perc = @bind metal_perc Slider(0:50, default = 25, show_value = true)
	md"""Normal conductance points from each side: $(metal_perc) %"""
end

# ╔═╡ 0c0bca6d-1364-4913-b8d1-6f55a26c8a69
md""" Apply normalization? $(@bind normalize CheckBox())"""

# ╔═╡ c9af2faa-1fb2-4bd3-9fce-d00cfe5d4aa1
md"""## *I(V)* & *dI/dV* en bloc"""

# ╔═╡ 9f2f3928-00d2-42a2-b53e-5499a84527f4
begin
	en_bloc
	md"#### Process all spectra? $(@bind proc CheckBox())"
end

# ╔═╡ c3f3be9e-65a2-47f2-9aff-abbd1c991f45
md"# Plot maps"

# ╔═╡ 1006d8dd-738c-42ff-b275-7521d2ebde33
md"###  *dI/dV* slice heatmap"

# ╔═╡ f3ba921b-3b47-4fab-b616-c478a0e91566
md"## *I(V)* spectrum @ clixel"

# ╔═╡ 049c19cd-91a6-4459-aa2e-55c1089e1b1d
md"""# Hic sunt functiones programmatoriae"""

# ╔═╡ 45e0787a-c377-4c84-b639-aeabf66e946f
md"## import .gsi file"

# ╔═╡ 20c7f0e7-9be5-4922-9209-f5bbd03eae98
"""
    read_gsi_header(io::IO)

Read a WSxM/GSI text header up to `[Header end]` and extract key metadata using regex.

Extracted scalar fields (if present):
- `begbytes`    : Image header size (byte offset where binary payload begins)
- `dattype`     : Image Data Type (e.g. "short", "simple", "double")
- `ncol`/`nrow` : Number of columns / rows
- `ramp`        : Number of points per ramp
- `f1`          : Conversion factor 0 for input channel (numeric part only)
- `x_amp`       : X Amplitude numeric value (from `[Control]`)
- `x_amp_unit`  : X Amplitude unit (e.g. "nm")
- `y_amp`       : Y Amplitude numeric value (from `[Control]`)
- `y_amp_unit`  : Y Amplitude unit (e.g. "nm")
- `z_amp`       : Z Amplitude numeric value (from `[General Info]`)
- `z_amp_unit`  : Z Amplitude unit (e.g. "nm")

Additionally parses the section:
- `[Spectroscopy images ramp value list]`

Returns:
- `ramp_unit`     : unit token from ramp list (e.g. "mV") or `missing`
- `ramp_values`   : ramp values in that unit (Float64[])
- `ramp_indices`  : indices (Int[])
- `ramp_values_V` : ramp values converted to volts if unit is V/mV/uV/µV/nV/pV, else `missing`

Notes:
- Parses the header in a single pass and preserves stream position at the end of header.
- Validates ramp-unit consistency and (optionally) ramp length vs `Number of points per ramp`.
"""
function read_gsi_header(io::IO; validate_ramp_length::Bool=true)
    begbytes   = missing
    dattype    = missing
    ncol       = missing
    nrow       = missing
    ramp       = missing
    f1         = missing

    x_amp      = missing
    x_amp_unit = missing
    y_amp      = missing
    y_amp_unit = missing
	z_amp      = missing
    z_amp_unit = missing

    ramp_unit     = missing
    ramp_values   = Float64[]
    ramp_indices  = Int[]
    ramp_values_V = missing

    pending_line::Union{Nothing,String} = nothing

    # Scalar field regexes
    rx_begbytes = r"^Image header size:\s*(\d+)\b"
    rx_dattype  = r"^Image Data Type:\s*([A-Za-z]+)\b"
    rx_ncol     = r"^Number of columns:\s*(\d+)\b"
    rx_nrow     = r"^Number of rows:\s*(\d+)\b"
    rx_ramp     = r"^Number of points per ramp:\s*(\d+)\b"
    rx_f1       = r"^Conversion factor 0 for input channel:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\b"

    # Control amplitudes
    rx_xamp = r"^X Amplitude:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*([^\s]+)\b"
    rx_yamp = r"^Y Amplitude:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*([^\s]+)\b"
    rx_zamp = r"^Z Amplitude:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*([^\s]+)\b"

    # Ramp list entry
    rx_ramp_line = r"^Image\s+(\d+)\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*([^\s]+)\s*$"

    function volt_factor(unit::AbstractString)
        u = strip(unit)
        return u == "V"  ? 1.0   :
               u == "mV" ? 1e-3  :
               (u == "uV" || u == "µV") ? 1e-6 :
               u == "nV" ? 1e-9  :
               u == "pV" ? 1e-12 :
               nothing
    end

    seekstart(io)

    while !eof(io)
        line = pending_line === nothing ? readline(io) : pending_line
        pending_line = nothing

        s = strip(line)
        s == "[Header end]" && break

        # ---- Ramp section ----
        if s == "[Spectroscopy images ramp value list]"
            while !eof(io)
                l2 = readline(io)
                s2 = strip(l2)

                # end section (includes "[Header end]")
                if startswith(s2, "[")
                    pending_line = l2
                    break
                end

                isempty(s2) && continue

                m = match(rx_ramp_line, s2)
                m === nothing && continue

                idx = parse(Int, m.captures[1])
                val = parse(Float64, m.captures[2])
                u   = m.captures[3]

                if ramp_unit === missing
                    ramp_unit = u
                elseif ramp_unit != u
                    error("Mixed units in ramp list: saw $(ramp_unit) then $u")
                end

                push!(ramp_indices, idx)
                push!(ramp_values, val)
            end

            if ramp_unit !== missing
                fac = volt_factor(ramp_unit)
                fac !== nothing && (ramp_values_V = ramp_values .* fac)
            end

            continue
        end

        # ---- X/Y/Z amplitude extraction ----
        if (m = match(rx_xamp, s)) !== nothing
            x_amp = parse(Float64, m.captures[1])
            x_amp_unit = m.captures[2]
            continue
        end

        if (m = match(rx_yamp, s)) !== nothing
            y_amp = parse(Float64, m.captures[1])
            y_amp_unit = m.captures[2]
            continue
        end

		if (m = match(rx_zamp, s)) !== nothing
            z_amp = parse(Float64, m.captures[1])
            z_amp_unit = m.captures[2]
            continue
        end

        # ---- Scalar fields ----
        if (m = match(rx_begbytes, s)) !== nothing
            begbytes = parse(Int, m.captures[1]); continue
        end

        if (m = match(rx_dattype, s)) !== nothing
            dattype = m.captures[1]; continue
        end

        if (m = match(rx_ncol, s)) !== nothing
            ncol = parse(Int, m.captures[1]); continue
        end

        if (m = match(rx_nrow, s)) !== nothing
            nrow = parse(Int, m.captures[1]); continue
        end

        if (m = match(rx_ramp, s)) !== nothing
            ramp = parse(Int, m.captures[1]); continue
        end

        if (m = match(rx_f1, s)) !== nothing
            f1 = parse(Float64, m.captures[1]); continue
        end
    end

    if validate_ramp_length && ramp !== missing && !isempty(ramp_values)
        length(ramp_values) == ramp ||
            @warn "Ramp list length ($(length(ramp_values))) != header ramp ($ramp)."
    end

    return (
        begbytes      = begbytes,
        dattype       = dattype,
        ncol          = ncol,
        nrow          = nrow,
        ramp          = ramp,
        f1            = f1,

        x_amp         = x_amp,
        x_amp_unit    = x_amp_unit,
        y_amp         = y_amp,
        y_amp_unit    = y_amp_unit,
		z_amp         = z_amp,
        z_amp_unit    = z_amp_unit,

        ramp_unit     = ramp_unit,
        ramp_values   = ramp_values,
        ramp_indices  = ramp_indices,
        ramp_values_V = ramp_values_V,
    )
end

# ╔═╡ 7b9eb7f5-cfe6-40c6-a093-0e9f20a8a19c
"""
Read topo + IV cube from a GSI file.

Assumes payload layout:
  payload starts at begbytes
  then topo: nrow x ncol values
  then iv:   nrow x ncol x ramp values

Returns (header, topo, iv).
"""
function read_gsi(io::IO; file_nbytes::Union{Nothing,Int}=nothing, max_elems::Int=600_000_000)
    header = read_gsi_header(io)

    # Validate required fields
    header.begbytes === missing && error("Missing 'Image header size' (begbytes) in header.")
    header.dattype  === missing && error("Missing 'Image Data Type' in header.")
    header.nrow     === missing && error("Missing 'Number of rows' in header.")
    header.ncol     === missing && error("Missing 'Number of columns' in header.")
    header.ramp     === missing && error("Missing 'Number of points per ramp' in header.")

    dt = header.dattype
    T = dt == "double" ? Float64 :
        dt == "simple" ? Float32 :
        dt == "short"  ? Int16  :
        error("Unsupported Image Data Type: $(repr(dt))")

    nrow = header.nrow
    ncol = header.ncol
    ramp = header.ramp

    # Use Int64 for safe size math
    nxy = Int64(nrow) * Int64(ncol)
    niv = nxy * Int64(ramp)
    total_vals = nxy + niv

    total_vals ≤ max_elems || error("Dataset too large ($total_vals values). Increase max_elems only if you really intend this.")

    expected_payload_bytes = total_vals * sizeof(T)
    expected_total_bytes = Int64(header.begbytes) + expected_payload_bytes

    if file_nbytes !== nothing
        if Int64(file_nbytes) < expected_total_bytes
            error("File too small: have $(file_nbytes) bytes, need at least $(expected_total_bytes).")
        elseif Int64(file_nbytes) != expected_total_bytes
            @warn "File size ($(file_nbytes)) differs from expected ($(expected_total_bytes)). Continuing (possible extra channels/trailer)."
        end
    end

    # Jump to binary payload start and read topo + iv
    seek(io, header.begbytes)

    topo_vec = Vector{T}(undef, nxy)
    read!(io, topo_vec)
    topo = reshape(topo_vec, nrow, ncol)

    iv_vec = Vector{T}(undef, niv)
    read!(io, iv_vec)
    iv = reshape(iv_vec, nrow, ncol, ramp)

    return header, topo, iv
end

# ╔═╡ c9856d06-9463-45eb-a686-f0b51b5d3629
if gsi === nothing
    "Upload a .gsi file above."
else
	if source == 2
		"TBA..."
	end
	if source == 1
	bytes = gsi["data"]          # <- String key
    name  = gsi["name"]
    mime  = get(gsi, "type", missing)

    io = IOBuffer(bytes)  # safe and seekable
    header, topo, unsorted_iv = read_gsi(io; file_nbytes=length(bytes))
	unsorted_bias = header.ramp_values_V.*1000 #bias in mV
	x_range = range(0, header.x_amp, length = header.ncol)
	y_range = range(0, header.y_amp, length = header.nrow)
	topo_unit = header.x_amp_unit
	z_unit = header.z_amp_unit

	# sort bias voltage from lo to hi
	perm = sortperm(unsorted_bias)
	raw_volts = unsorted_bias[perm]
	raw_iv = unsorted_iv[:, :, perm]
		
	header
    end
end

# ╔═╡ 05594a29-1279-4428-80cb-940268e60e6b
println("grid file: ",name)

# ╔═╡ b1aa9714-624e-4940-a57a-cc2dacf8a6ff
begin
unclick
	@bind click let
    p = plot(
        heatmap(x=x_range, y=y_range, z=topo, colorscale="Viridis",
                colorbar=attr(title="z [$(z_unit)]")),
        Layout(
            title="Topography",
            xaxis=attr(title="x [$(topo_unit)]"),
            yaxis=attr(title="y [$(topo_unit)]", scaleanchor="x", scaleratio=1),
            width=650, height=650
        )
    )

    add_plotly_listener!(p, "plotly_click", """
    (e) => {
        const dt = e.points[0];
        const r = dt.pointNumber[0] + 1;
        const c = dt.pointNumber[1] + 1;
        PLOT.value = [r, c];
        PLOT.dispatchEvent(new CustomEvent('input'));
    }
    """)

    # When cell reruns because `reset` changed, immediately clear the value:
    add_plotly_listener!(p, "plotly_afterplot", """
    (e) => {
        // `reset` is captured by Pluto reactivity: rerun = reset pressed
        PLOT.value = null;
        PLOT.dispatchEvent(new CustomEvent('input'));
    }
    """)

    p
end
end

# ╔═╡ 00d9b5f3-1d87-4443-aaca-7a20feefc016
println("click index: ",click)

# ╔═╡ 50474186-101b-45de-9540-cd0bedc75323
begin
    if click === nothing
		col = 1
		row = 1
	else
		row = click[1]
		col = click[2]
	end
	raw_curs = vec(raw_iv[row, col, :])
	gr.plot(raw_volts, raw_curs, label = "row=$(row), col=$(col)",xlabel ="U[mv]",ylabel="I[pA]]")
end

# ╔═╡ 8f73478c-1edf-4a44-8145-71e303d5d95e
println("position: x = ", x_range[col], " y = ", y_range[row])

# ╔═╡ c5ee2976-cdfc-4f41-8608-3fb46cc8e1c2
begin
	click_surf = gr.heatmap(x_range, y_range, topo, xlabel = "x [$(topo_unit)]", ylabel = "y [$(topo_unit)]", c = :viridis, colorbar_title = "z [$(topo_unit)]")
	if !(ismissing(click) || click === nothing)
        # If your click is a Dict with "row"/"col" keys:
        gr.scatter!([x_range[col]],[y_range[row]], marker = :x, msw = 5, label = false)
    end
	click_surf
end

# ╔═╡ 46efaead-594d-49ec-abc4-97b4db9c1eb3
@bind bias Slider(eachindex(raw_volts); default = 1, show_value = i -> "$(raw_volts[i]) mV")

# ╔═╡ 4a8c9f62-c65d-449b-b9b1-7ab41471f576
begin
unclick
	@bind click2 let
    p2 = plot(
        heatmap(x=x_range, y=y_range, z=topo, colorscale="Viridis",
                colorbar=attr(title="z [$(z_unit)]")),
        Layout(
            title="Topography",
            xaxis=attr(title="x [$(topo_unit)]"),
            yaxis=attr(title="y [$(topo_unit)]", scaleanchor="x", scaleratio=1),
            width=650, height=650
        )
    )

    add_plotly_listener!(p2, "plotly_click", """
    (e) => {
        const dt = e.points[0];
        const r = dt.pointNumber[0] + 1;
        const c = dt.pointNumber[1] + 1;
        PLOT.value = [r, c];
        PLOT.dispatchEvent(new CustomEvent('input'));
    }
    """)

    # When cell reruns because `reset` changed, immediately clear the value:
    add_plotly_listener!(p2, "plotly_afterplot", """
    (e) => {
        // `reset` is captured by Pluto reactivity: rerun = reset pressed
        PLOT.value = null;
        PLOT.dispatchEvent(new CustomEvent('input'));
    }
    """)

    p2
end
end

# ╔═╡ 3452ba8e-eadb-4827-a6d0-f18a238ac807
md"## differentiation"

# ╔═╡ a120eb8e-0010-4da7-ba19-b95317480846
"""
    linfit(ydata; xdata=nothing)

Fits y ~ a*x + b with linear least squares (no iterative solver).
If `xdata` is not provided, uses the axis of `ydata`.
Returns a NamedTuple: (slope, intercept, r2, residuals, predict).
"""
function linfit(ydata::AbstractVector{<:Real}; xdata=nothing)
	# use linear_fit(curs; xdata = volts)
   	#=  =#
	
	n = length(ydata)
    n ≥ 2 || throw(ArgumentError("Need at least two points"))

    xdata === nothing && (xdata = collect(axes(ydata, 1)))
    x = collect(float.(xdata))
    y = collect(float.(ydata))

    X = hcat(x, ones(n))
    β = X \ y
    slope, intercept = β[1], β[2]

	#complete output
    ŷ = X * β
    resid = y .- ŷ
    ss_res = sum(abs2, resid)
    ss_tot = sum(abs2, y .- mean(y))
    r2 = 1 - ss_res / ss_tot

    predict = (xnew -> slope .* xnew .+ intercept)
    return (slope=slope, intercept=intercept, r2=r2, residuals=resid, predict=predict)
end

# ╔═╡ 7e481848-564d-4270-80d5-c10f07883f4e
"""
    pts_dif(volts, currents; pts::Integer = 1)

Differentiates using linfit().slope over a predifined number of points.
"""
function pts_dif(volts, currents; pts::Integer = 1)
    n = length(volts)
    n == length(currents) || throw(ArgumentError("current and volts must have same length"))
    pts ≥ 0 || throw(ArgumentError("pts (window radius) must be nonnegative"))

    didv = zeros(n)  # or choose a type based on inputs, see below

    for i in 1:n
        start = max(1, i - pts)
        stop  = min(n, i + pts)

        # Fit I (y) as a function of V (x) to get dI/dV
        coefs = linfit(currents[start:stop]; xdata = volts[start:stop])
        didv[i] = coefs.slope
    end

    return didv
end

# ╔═╡ d23920d6-a4b1-423d-8f90-710dcf905db3
begin
	# near-uniform grid (your data satisfy this at ~0.01%)
    @assert all(isapprox.(diff(raw_volts), diff(raw_volts); rtol=1e-3)) "Use uniform V or resample first"
	#create a vector of constant intervals for each points
	pts = ones(length(raw_volts)).*points
	raw_dif = pts_dif(raw_volts, raw_curs; pts = points)
	#plot derivative
	gr.plot(raw_volts, raw_dif, label = "$(points) points diff of raw I(V)", xlabel = "voltage [mV]") 
end

# ╔═╡ 9b1d9de3-20af-4f4d-8403-2d5e22a7bb12
let
    ymin, ymax = extrema(raw_dif)

    tr_data = scatter(
        x = raw_volts,
        y = raw_dif,
        mode = "lines",
        name = "",
        showlegend = false
    )

    vline_trace(x; name, color, dash="dash") = scatter(
        x = [x, x],
        y = [ymin, ymax],
        mode = "lines",
        name = name,
        line = attr(color=color, width=2, dash=dash),
        hoverinfo = "skip"
    )

    traces = [
        tr_data,
        vline_trace(-left_peak;  name="left peak",  color="royalblue"),
        vline_trace( right_peak; name="right peak", color="crimson"),
        vline_trace( bioff;      name="bias offset", color="black", dash="dot"),
    ]

    lay = Layout(
        title = "Bias offset: $(bioff) mV.",
        xaxis = attr(title = "voltage [mV]"),
        yaxis = attr(title = "dI/dV [a.u.]"),
        legend = attr(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        uirevision = 1
    )

    plot(traces, lay)
end


# ╔═╡ bc524503-7132-46e2-bb80-e69d0d10d3dd
begin
	volts = raw_volts .- bias_offset
	gr.plot(volts, raw_dif, xlabel = "bias [mV]", ylabel = "dI/dV [a.u.]", title = "bias offset removed", label = false)
end

# ╔═╡ ce788855-e6ee-402a-a2eb-2c06b0c752b1
begin
	#zero bias index...abs.(vector) computes the absolute value of each element in the vector, argmin returns the index of the smallest value in the resulting array, which corresponds to the value closest to 0 in the original vector
	zbi = argmin(abs.(volts))
	
	md" Bias voltage closest to 0 is $(volts[zbi])
	@ index $zbi"
end

# ╔═╡ cf585ea8-e949-4573-8700-13d394d16b81
begin
	curs_lw = lowess(volts, raw_curs, lwf, 3, 0.0)
	gr.plot(volts, raw_curs, label = "raw data", title = "f = $(lwf)", xlabel = "bias [mV]", ylabel = "tunneling current [pA]", legend = :bottomright)
	gr.plot!(volts, curs_lw, label = "LOWESS", lw = 2)
end

# ╔═╡ 1269ecf1-82b1-426f-95c8-5dd73d18040e
begin
	curs_offset = linfit(curs_lw[zbi-curoff_pts:zbi+curoff_pts]; xdata = volts[zbi-curoff_pts:zbi+curoff_pts]).intercept
	md"""Current offset calculated as the intercept of the linear fit of currents around zero bias with indexes $([zbi-curoff_pts:zbi+curoff_pts]) is $(curs_offset) *pA*.
	"""
end

# ╔═╡ e89da3e2-9505-4c3e-9c47-67dcac520b20
begin
	curs = curs_lw.-curs_offset
	gr.plot(volts, curs, label = "LOWESS - offset", xlabel = "bias [mV]", ylabel = "tunneling current [pA]", legend = :bottomright)
end

# ╔═╡ 3b24d2af-6ac5-43e0-8ac1-280e6a1852dc
println("current @ zbi = ", curs[zbi])

# ╔═╡ af2c4054-4393-4540-8677-ebba7756ca99
begin
md"""
**Smoothing parameters:**\
``\alpha_{gap} =`` $(@bind al_gap Slider(-10:0.1:1.0, default = -4.7, show_value = true)) 
``\alpha_{out} =`` $(@bind al_out Slider(-10:0.1:1.0, default = -0.7, show_value = true))\
**Segment voltages:**\
``U_{gap} =`` $(@bind u_gap Slider(0:0.01:maximum(volts), default = 0.5, show_value = true))\
**Segment transition half-widths:**\
``w_{gap} =`` $(@bind wg Slider(1:500, default = 40, show_value = true)) ``dV``
"""
end

# ╔═╡ 4caebd47-2069-4d18-be70-5d7166a7f60f
begin
	dv = mean(diff(raw_volts))
	
	# near-uniform grid (your data satisfy this at ~0.01%)
    @assert all(isapprox.(diff(raw_volts), dv; rtol=1e-3)) "Use uniform V or resample first"
	
	fs = 1/dv # Sample rate
	freqs =  fftshift(fftfreq(length(raw_volts), fs))
	F2 = fftshift(fft(pts_dif(raw_volts, raw_curs; pts = points)))
	gr.plot(freqs, abs.(F2), xlabel="frequency in 1/mV", label="fft", xlims = (0, :auto))
end

# ╔═╡ 6dfe666b-3b7e-47fa-ac54-c1666cce044f
"""
    cross_nrd(volts, curs; iter=200, kernel="square", al_gap, al_out, u_gap, wg=5,
                  dv=nothing, scale="small", eps=1e-6) -> derivative & crossfade

Compute a blended derivative from two total-variation (TV) regularized derivatives:
one tuned for the gap region (`al_gap`) and one for the outer region (`al_out`),
then **raised-cosine crossfade** them around ±`u_gap` to ensure continuity.
"""
function cross_nrd(volts::AbstractVector, curs::AbstractVector;
                       iter::Integer=50,
                       kernel="square",
                       al_gap::Real,
                       al_out::Real,
                       u_gap::Real,
                       wg::Real,
                       dv::Union{Nothing,Real}=nothing,
                       scale::AbstractString="small",
                       eps::Real=1e-6)

    n = length(curs)
    length(volts) == n || throw(ArgumentError("volts and curs must have the same length"))

    # Infer dv if not provided; check uniform spacing
    if dv === nothing
        if volts isa AbstractRange
            dv = step(volts)
        else
            dv = mean(diff(volts))
            # sanity check for uniform grid
            maxerr = maximum(abs, diff(volts) .- dv)
            maxerr ≤ 1e-8*max(1.0, abs(dv)) ||
                @warn "volts not perfectly uniform (max deviation = $maxerr)"
        end
    end

    wid = wg * dv
    wid > 0 || throw(ArgumentError("Crossfade half-width wid = wg*dv must be > 0"))

    # Two TV-regularized derivatives (α given as log10 exponents)
    α_gap = 10.0^al_gap
    α_out = 10.0^al_out
    d_gap = tvdiff(curs, iter, α_gap; dx=dv, scale=scale, ε=eps, diff_kernel=kernel)
    d_out = tvdiff(curs, iter, α_out; dx=dv, scale=scale, ε=eps, diff_kernel=kernel)

    # Raised-cosine crossfade around ±u_gap
    absV = abs.(volts)
    crossfade = 0.5 .* (1 .- cospi.(clamp.((absV .- (u_gap - wid)) ./ (2*wid), 0.0, 1.0)))  # 0→1
    wout = crossfade
    wgap = 1 .- crossfade

    # Blended derivative (continuous across ±u_gap)
    seg = wgap .* d_gap .+ wout .* d_out

    return (seg = seg, crossfade = crossfade)

	#= Steps
1) `d_gap = tvdiff(curs, iter, 10^al_gap; dx=dv, scale, ε=eps, diff_kernel=kernel)`
2) `d_out = tvdiff(curs, iter, 10^al_out; dx=dv, scale, ε=eps, diff_kernel=kernel)`
3) Build a Hahn/Hann-style crossfade around ±`u_gap` with half-width `wid = wg*dv`:
       λ_g(V) = 0.5 * (1 - cos(π * clamp((|V| - (u_gap - wid)) / (2*wid), 0, 1)))
   Then blend:
       seg = (1 - λ_g) .* d_gap .+ λ_g .* d_out

Arguments
- `volts`::AbstractVector  — x-axis (bias), uniformly spaced.
- `curs` ::AbstractVector  — y-values (trace), same length as `volts`.
- `iter` ::Integer         — iteration budget for `tvdiff` (default 200).
- `kernel`                  — optional derivative kernel passed to `tvdiff` as `diff_kernel`. Default "square"
- `al_gap`::Real           — log10(α) used near the gap (α = 10^al_gap).
- `al_out`::Real           — log10(α) used in the outer region (α = 10^al_out).
- `u_gap` ::Real           — center of the crossfade in |V| units (same units as `volts`).
- `wg`    ::Real           — half-width in **samples** for the crossfade (default 5); `wid = wg*dv`.
- `dv`    ::Union{Nothing,Real} — sample spacing in `volts` (if `nothing`, it is inferred).
- `scale` ::AbstractString  — passed to `tvdiff` (default "small").
- `eps`   ::Real            — ε tolerance passed to `tvdiff` (default 1e-6).
- `return_all`::Bool        — if `true`, returns a NamedTuple with internals.

Notes
- `wg` is in **samples**, so the crossfade **energy half-width** is `wid = wg * dv`.
- The raised-cosine gives λ_g = 0 below `u_gap - wid` and λ_g = 1 above `u_gap + wid`,
  smoothly transitioning in between.
- Assumes `volts` is uniformly spaced; if `dv` is not given, it is inferred and checked. =#
end

# ╔═╡ f0ecf5f4-7552-4b3c-9ab5-4afc93c14ead
begin
    seg_der = cross_nrd(volts, curs; al_gap = al_gap, al_out = al_out, u_gap = u_gap, wg = wg, kernel = kernel, iter = iter)
	seg = seg_der.seg
    # 4) plot
	crossfade = seg_der.crossfade
	gr.plot(volts, raw_dif, label = "$(points) points diff of raw I(V)")
	gr.plot!(volts, seg, xlabel = "bias [mV]", label="dI/dV (crossfaded)")
	gr.plot!(volts, crossfade.*maximum(seg), ls=:dash, xlabel = "bias [mV]", ylabel = "dI/dV [a.u.]", label="crossfade", legend=:bottomright, lc = :black)
end

# ╔═╡ f5960c88-4a4f-4f2b-b5d6-e3dc3d8b1f50
begin
    gr.plot(volts, crossfade.*maximum(seg), ls=:dash, xlabel = "bias [mV]", ylabel = "dI/dV [a.u.]", label="crossfade", lc = :black, legend=:bottomleft)
    gr.plot!(volts, seg, label="dI/dV (crossfaded)")
	gr.plot!(gr.twinx(), volts, pts_dif(volts, seg), label = "d²I/dV²", lc = 3, legend = :bottomright)	
end

# ╔═╡ dd801718-4a17-4eb2-969b-8cc291142f6d
md"## normalization"

# ╔═╡ 84c8f90b-0871-4086-91cf-973842daec47
"""
    normetal(x, y, percent)

Normalize by fitting a linear baseline to the outer `percent` of points (from both ends) and dividing the entire series by that baseline.

- `percent` is the percentage of total length taken from each end (so 2×`percent`% total points).
"""
function normetal(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, percent::Real)
    n = length(y)
    n ≥ 4 || throw(ArgumentError("Need at least 4 points"))
 	# full x axis
    length(x) == n || throw(ArgumentError("xdata must have the same length as y"))
    
	# number of points from EACH side
    pts = round(Int, n * percent / 100)

	# indices used for the baseline fit (outer segments)
    idx = vcat(1:pts, n-pts+1:n)

    # fit baseline to outer points, with correct xdata
    fit = linfit(y[idx]; xdata = x[idx])

    # evaluate baseline across the whole axis
    baseline = fit.predict(x)
    return baseline
end

# ╔═╡ 6f9e6764-9f30-428e-a6e7-2289704a56f6
begin
	bkg = normetal(volts, seg, metal_perc)
	metal_pts = round(Int, length(volts) * metal_perc / 100)

	gr.plot(volts, seg, xlabel = "bias [mV]", ylabel = "dI/dV [a.u.]", label="dI/dV (crossfaded)", legend = :bottomright)
	gr.plot!(volts, bkg, label = "linear fit")
	gr.vline!([volts[metal_pts] volts[end-metal_pts]], label = "$(metal_perc)%")
end

# ╔═╡ 90a8da52-2d64-4cd6-9566-c0ba13e70c66
begin
	norm_conds = seg ./ bkg
	gr.plot(volts, norm_conds, xlabel = "bias [mV]", ylabel = "dI/dV [a.u.]", label = "normalized dI/dV")
	gr.hline!([1], ls = :dash, label = "1")
	gr.hline!([maximum(norm_conds)], ls = :dash,label = "maximum")
end

# ╔═╡ ee57e7db-0875-4df3-9f3b-db8e9157da31
md"""## *en bloc* processing"""

# ╔═╡ c752f559-a0ea-4d54-8e4c-e5e9ee7708f1
"""
Process a 3D spectroscopy cube `I` (nx×ny×nE) by applying, for each spectrum v = I[i,j,:]:

1) smoothing: lowess(volts, curs, lwf, 3, 0.0)
2) offset removal: subtract intercept of linfit on a window around zbi

Input: grid, volts, lowess factor, zbi, current offset points

Returns: Array{Float64,3} with same size as `I`.
"""
function iv_en_bloc(I::AbstractArray{<:Real,3},
                                 volts::AbstractVector,
                                 lwf,
                                 zbi::Integer,
                                 curoff_pts::Integer)

    nx, ny, nE = size(I)
    @assert length(volts) == nE "length(volts) must equal size(I,3)"

    out = Array{Float64}(undef, nx, ny, nE)
    buf = Vector{Float64}(undef, nE)   # reusable buffer

    @views for i in 1:nx, j in 1:ny
        # fill Float64 buffer
        @inbounds for k in 1:nE
            buf[k] = Float64(I[i, j, k])
        end

        # smooth
        s = lowess(volts, buf, lwf, 3, 0.0)

        # fit window (clamped)
        lo = max(1,  zbi - curoff_pts)
        hi = min(nE, zbi + curoff_pts)

        # IMPORTANT: use view(...) not @view ...
        b0 = linfit(view(s, lo:hi)).intercept

        # write result
        @inbounds for k in 1:nE
            out[i, j, k] = s[k] - b0
        end
    end

    return out
end

# ╔═╡ 7ac91382-2102-4cb2-b2c4-a2934591968b
"""
Process a 3D spectroscopy cube `I` (nx×ny×nE) by applying, for each spectrum v = I[i,j,:]:

1) crossfaded Noise Robust Differentiation
2) optional normalization

Input: grid, volts, al_gap, al_out, u_gap, wg, dv, kernel, iter, normalize::Bool=false, metal_perc=nothing

Returns: Array{Float64,3} with same size as `I`.
"""
function didv_en_bloc(I::AbstractArray{<:Real,3},
                                volts::AbstractVector;
                                # cross_nrd kwargs
                                al_gap, al_out, u_gap, wg, dv, kernel, iter,
                                # normalization control
                                normalize::Bool=false,
                                metal_perc=nothing)

    nx, ny, nE = size(I)
    @assert length(volts) == nE "length(volts) must equal size(I,3)"

    out = Array{Float64}(undef, nx, ny, nE)

    # Reusable Float64 buffer for one spectrum (avoids per-spectrum allocations)
    buf = Vector{Float64}(undef, nE)

    @views for i in 1:nx, j in 1:ny
        # Load spectrum into Float64 buffer (works for Int16 or float inputs)
        @inbounds for k in 1:nE
            buf[k] = Float64(I[i, j, k])
        end

        # 1) cross_nrd step
        seg = cross_nrd(volts, buf;
            al_gap=al_gap, al_out=al_out, u_gap=u_gap,
            wg=wg, dv=dv, kernel=kernel, iter=iter
        ).seg

        # 2) optional normalization
        if normalize
            @assert metal_perc !== nothing "metal_perc must be provided when normalize=true"
            nrm = normetal(volts, seg, metal_perc)
            # nrm might be scalar or vector; broadcasting handles both
            @inbounds for k in 1:nE
                out[i, j, k] = Float64(seg[k] / nrm[k])  # works if nrm is a vector
            end
        else
            @inbounds for k in 1:nE
                out[i, j, k] = Float64(seg[k])
            end
        end
    end

    return out
end

# ╔═╡ 681fe248-8d90-4ae0-bbaa-0f37ef614f7b
begin
	if proc == true
		iv = iv_en_bloc(raw_iv, volts, lwf, zbi,curoff_pts)
		didv =  didv_en_bloc(iv, volts; al_gap = al_gap, al_out = al_out, u_gap = u_gap, wg = wg, dv = dv, kernel = kernel, iter = iter, normalize = normalize, metal_perc= metal_perc)
				
		else
		println("Check please!")
	end
end

# ╔═╡ 384426f4-fea3-41a1-9328-d486eb7ba3f1
gr.heatmap(x_range, y_range, didv[:, :, bias], c = gr.cgrad(:RdBu, rev = true), colorbar_title="I[pA]", xlabel = "x [$(topo_unit)]", ylabel = "y [$(topo_unit)]")

# ╔═╡ 60993909-da30-49de-8ee6-3875a2df2174
begin
	tr = surface(x = x_range, y = y_range, z = topo,
	    surfacecolor = didv[:, :, bias],     # <-- heatmap painted on surface
	    colorscale = "Viridis",
	    cmin = minimum(didv[:, :, bias]),
	    cmax = maximum(didv[:, :, bias]),
	    colorbar = attr(title="IV value"),
	    showscale = true
	)
	
	lay = Layout(
	    title = "Surface colored by heatmap (surfacecolor) @ $(raw_volts[bias]) mV",
	    scene = attr(
	        xaxis = attr(title="x [$(topo_unit)]"),
	        yaxis = attr(title="y [$(topo_unit)]"),
	        zaxis = attr(title="z [$(z_unit)]")
	    ),
	    width = 600, height = 400, uirevision = "keep"
	)
	
	plot(tr, lay)
end

# ╔═╡ c8a5e136-cc38-4138-99a5-de264cfcd15c
begin
    if click2 === nothing
		col2 = 1
		row2 = 1
	else
		row2 = click2[1]
		col2 = click2[2]
	end
	check_iv = vec(iv[row, col, :])
	check_didv = vec(didv[row2, col2, :])
	gr.plot(volts, check_didv, label = "row=$(row2), col=$(col2)",xlabel ="U[mv]",ylabel="I[pA]]")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
Lowess = "a2847a04-7795-4581-bff1-964344ea125d"
NoiseRobustDifferentiation = "470638dc-0858-4731-a73a-678bdc45695b"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
FFTW = "~1.10.0"
JLD2 = "~0.6.3"
Lowess = "~0.1.0"
NoiseRobustDifferentiation = "~0.2.4"
Plots = "~1.41.6"
PlutoPlotly = "~0.6.5"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "13b18bd2b70baf6b54e889b68244ac59b77d0e07"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AMD]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "45a1272e3f809d36431e57ab22703c6896b8908f"
uuid = "14f7f29c-3bd6-536c-9a0b-7339e30b5a3e"
version = "0.5.3"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "2eeb2c9bef11013efc6f8f97f32ee59b146b09fb"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.44"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "35ea197a51ce46fcd01c4a44befce0578a1aaeca"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AlgebraicMultigrid]]
deps = ["CommonSolve", "LinearAlgebra", "LinearSolve", "Printf", "Reexport", "SparseArrays", "TimerOutputs"]
git-tree-sha1 = "223bc2d88e05f471586fbb79bde602fc372f4add"
uuid = "2169fc97-5a83-5252-b627-83903c6c433c"
version = "1.2.0"

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
git-tree-sha1 = "78b3a7a536b4b0a747a0f296ea77091ca0a9f9a3"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.23.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChunkCodecCore]]
git-tree-sha1 = "1a3ad7e16a321667698a19e77362b35a1e94c544"
uuid = "0b6fb165-00bc-4d37-ab8b-79f91016dbe1"
version = "1.0.1"

[[deps.ChunkCodecLibZlib]]
deps = ["ChunkCodecCore", "Zlib_jll"]
git-tree-sha1 = "cee8104904c53d39eb94fd06cbe60cb5acde7177"
uuid = "4c0bbee4-addc-4d73-81a0-b6caacae83c8"
version = "1.0.0"

[[deps.ChunkCodecLibZstd]]
deps = ["ChunkCodecCore", "Zstd_jll"]
git-tree-sha1 = "34d9873079e4cb3d0c62926a225136824677073f"
uuid = "55437552-ac27-4d47-9aa3-63184e8fd398"
version = "1.0.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

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
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6522cfb3b8fe97bec632252263057996cbd3de20"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.18.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers", "PrecompileTools", "TruncatedStacktraces"]
git-tree-sha1 = "1b14f69c661abc9d21264eeddb3ec00d517134df"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "1.3.0"

    [deps.FunctionWrappersWrappers.extensions]
    FunctionWrappersWrappersEnzymeExt = ["Enzyme", "EnzymeCore"]
    FunctionWrappersWrappersMooncakeExt = "Mooncake"

    [deps.FunctionWrappersWrappers.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

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
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "59545b0a2b27208b0650df0a46b8e3019f85055b"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["ChunkCodecLibZlib", "ChunkCodecLibZstd", "FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues"]
git-tree-sha1 = "8f8ff711442d1f4cfc0d86133e7ee03d62ec9b98"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.6.3"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "c4d19f51afc7ba2afbe32031b8f2d21b11c9e26e"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LimitedLDLFactorizations]]
deps = ["AMD", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "0a79513936f19e9829c15f6c5ff68e571a4b1ba1"
uuid = "f5a24dde-3ab7-510b-b81b-6a72c6098d3b"
version = "0.5.2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LinearMaps]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7f6be2e4cdaaf558623d93113d6ddade7b916209"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.4"
weakdeps = ["ChainRulesCore", "SparseArrays", "Statistics"]

    [deps.LinearMaps.extensions]
    LinearMapsChainRulesCoreExt = "ChainRulesCore"
    LinearMapsSparseArraysExt = "SparseArrays"
    LinearMapsStatisticsExt = "Statistics"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "35b6d3b3cfa50f97a7e13bfdd3b82d9f9fcb68af"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.72.0"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveAlgebraicMultigridExt = "AlgebraicMultigrid"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveChainRulesCoreExt = "ChainRulesCore"
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveElementalExt = "Elemental"
    LinearSolveEnzymeExt = ["EnzymeCore", "SparseArrays"]
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveGinkgoExt = ["Ginkgo", "SparseArrays"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolveMooncakeExt = "Mooncake"
    LinearSolvePETScExt = ["PETSc", "SparseArrays"]
    LinearSolveParUExt = ["ParU_jll", "SparseArrays"]
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    AlgebraicMultigrid = "2169fc97-5a83-5252-b627-83903c6c433c"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Elemental = "902c3f28-d1ec-5e7e-8399-a24c3845ee38"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Ginkgo = "4c8bd3c9-ead9-4b5e-a625-08f1338ba0ec"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    ParU_jll = "9e0b026c-e8ce-559c-a2c4-6a3d5c955bc9"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

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
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.Lowess]]
deps = ["Interpolations"]
git-tree-sha1 = "e661d3cd1e4ab170b29e8391f66b5216479234b3"
uuid = "a2847a04-7795-4581-bff1-964344ea125d"
version = "0.1.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.NoiseRobustDifferentiation]]
deps = ["IterativeSolvers", "LinearAlgebra", "LinearMaps", "Preconditioners", "SparseArrays"]
git-tree-sha1 = "108fc03c3419164898540545988d0df13fa5239b"
uuid = "470638dc-0858-4731-a73a-678bdc45695b"
version = "0.2.4"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
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
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "6256ab3ee24ef079b3afa310593817e069925eeb"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.23"

    [deps.PlotlyBase.extensions]
    DataFramesExt = "DataFrames"
    DistributionsExt = "Distributions"
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyBase.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

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

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "8acd04abc9a636ef57004f4c2e6f3f6ed4611099"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.5"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "e16b73bf892c55d16d53c9c0dbd0fb31cb7e25da"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "1.2.0"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preconditioners]]
deps = ["AlgebraicMultigrid", "LimitedLDLFactorizations", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "85c25d9f5a9e2941129aa9a314c613b5fbbeb269"
uuid = "af69fa37-3177-5a40-98ee-561f696e4fcd"
version = "0.6.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

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

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "PrecompileTools", "RecipesBase", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "d0282d612f22dcad7b81cf487b746e63aa2a6709"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.54.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsFastBroadcastPolyesterExt = ["FastBroadcast", "Polyester"]
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStatisticsExt = "Statistics"
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

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
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cfcdc949c4660544ab0fdeed169561cb22f835f4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.18"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLLogging", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "908c0bf271604d09393a21c142116ab26f66f67c"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.154.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDifferentiationInterfaceExt = "DifferentiationInterface"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLLogging]]
deps = ["Logging", "LoggingExtras", "Preferences"]
git-tree-sha1 = "0161be062570af4042cf6f69e3d5d0b0555b6927"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "1.9.1"

    [deps.SciMLLogging.extensions]
    SciMLLoggingTracyExt = "Tracy"

    [deps.SciMLLogging.weakdeps]
    Tracy = "e689c965-62c8-4b79-b2c5-8359227902fd"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "234869cf9fee9258a95464b7a7065cc7be84db00"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.16.0"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.SciMLStructures]]
deps = ["ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "607f6867d0b0553e98fc7f725c9f9f13b4d01a32"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.10.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "ac4b837d89a58c848e85e698e2a2514e9d59d8f6"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.6.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

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
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

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
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "28145feabf717c5d65c1d5e09747ee7b1ff3ed13"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.3"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

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

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

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
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "5f24e158cf4cee437052371455fe361f526da062"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.6"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e2a7072fc0cdd7949528c1455a3e5da4122e1153"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.56+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─a5540242-d382-4ac5-87b9-c0c45f74cef3
# ╟─6a6ea7ad-0a03-464d-8b84-4f642db7cfdd
# ╟─6df8df70-3346-11f1-a692-31d015c5335e
# ╟─89ec1e0d-3811-4ef7-9bbc-4a0f9d8d41b4
# ╟─c9856d06-9463-45eb-a686-f0b51b5d3629
# ╟─05594a29-1279-4428-80cb-940268e60e6b
# ╟─76884c57-b7e3-46b2-b095-5fa39a0de81f
# ╟─a52ec94e-a22e-4c7f-86af-c9eaf0cffdca
# ╟─c0f60871-fde9-429a-8389-f014844cd2a2
# ╟─557407a7-fcda-48a6-9a60-2f7626e579dc
# ╟─dd53ab84-5693-4c87-bb41-710841cc02aa
# ╟─b1aa9714-624e-4940-a57a-cc2dacf8a6ff
# ╟─181f2fab-cc28-4016-b9a9-ce0db79ef2ab
# ╟─00d9b5f3-1d87-4443-aaca-7a20feefc016
# ╟─50474186-101b-45de-9540-cd0bedc75323
# ╟─8f73478c-1edf-4a44-8145-71e303d5d95e
# ╟─c5ee2976-cdfc-4f41-8608-3fb46cc8e1c2
# ╟─94bf84cc-78ab-4612-b673-d12baf9b74e0
# ╟─6d5bd5bb-5e4b-41d4-ab70-9066ce13bddc
# ╟─ae073f35-b095-491c-8347-33e9b69764f9
# ╟─d23920d6-a4b1-423d-8f90-710dcf905db3
# ╟─08296e70-9f4b-4aca-834d-21c9a4e7d1d7
# ╟─4caebd47-2069-4d18-be70-5d7166a7f60f
# ╟─22107084-f3f7-4d1e-acbb-6ab4f8a87ab3
# ╟─1b05e239-7419-4aeb-81ce-9c8e66ad4b51
# ╟─1feba41d-4dfd-4757-abf4-f3397ac48b65
# ╟─1e60e654-0a68-4efd-a870-0d56b65e9731
# ╟─92e25728-fa53-4910-a779-c739c3eee693
# ╟─9b1d9de3-20af-4f4d-8403-2d5e22a7bb12
# ╟─73f16b28-e715-4352-9698-e2679fd664ce
# ╟─bc524503-7132-46e2-bb80-e69d0d10d3dd
# ╟─ce788855-e6ee-402a-a2eb-2c06b0c752b1
# ╟─c2547995-2493-4e1e-a956-7ab960ccf5c5
# ╟─9bcc955e-6fcf-49e2-a548-09368e2c6df8
# ╟─cf585ea8-e949-4573-8700-13d394d16b81
# ╟─fd48b9a7-f51b-4ee8-b3dd-f3982154414f
# ╟─1269ecf1-82b1-426f-95c8-5dd73d18040e
# ╟─e89da3e2-9505-4c3e-9c47-67dcac520b20
# ╟─3b24d2af-6ac5-43e0-8ac1-280e6a1852dc
# ╟─b092b5d4-66ef-48e5-b218-856b4b108d49
# ╟─20513a32-5bb6-44cd-9f68-869ca54a01ff
# ╟─af2c4054-4393-4540-8677-ebba7756ca99
# ╟─f0ecf5f4-7552-4b3c-9ab5-4afc93c14ead
# ╟─f5960c88-4a4f-4f2b-b5d6-e3dc3d8b1f50
# ╟─e2375ff2-7095-4751-80ea-eb3a3d1a6f5b
# ╟─208a0062-fc5b-4494-9ae2-29b02c8c3935
# ╟─6f9e6764-9f30-428e-a6e7-2289704a56f6
# ╟─90a8da52-2d64-4cd6-9566-c0ba13e70c66
# ╟─0c0bca6d-1364-4913-b8d1-6f55a26c8a69
# ╟─c9af2faa-1fb2-4bd3-9fce-d00cfe5d4aa1
# ╟─9f2f3928-00d2-42a2-b53e-5499a84527f4
# ╟─681fe248-8d90-4ae0-bbaa-0f37ef614f7b
# ╟─c3f3be9e-65a2-47f2-9aff-abbd1c991f45
# ╟─1006d8dd-738c-42ff-b275-7521d2ebde33
# ╟─46efaead-594d-49ec-abc4-97b4db9c1eb3
# ╟─384426f4-fea3-41a1-9328-d486eb7ba3f1
# ╟─60993909-da30-49de-8ee6-3875a2df2174
# ╟─4a8c9f62-c65d-449b-b9b1-7ab41471f576
# ╟─f3ba921b-3b47-4fab-b616-c478a0e91566
# ╟─c8a5e136-cc38-4138-99a5-de264cfcd15c
# ╟─049c19cd-91a6-4459-aa2e-55c1089e1b1d
# ╟─45e0787a-c377-4c84-b639-aeabf66e946f
# ╟─20c7f0e7-9be5-4922-9209-f5bbd03eae98
# ╟─7b9eb7f5-cfe6-40c6-a093-0e9f20a8a19c
# ╟─3452ba8e-eadb-4827-a6d0-f18a238ac807
# ╟─a120eb8e-0010-4da7-ba19-b95317480846
# ╟─7e481848-564d-4270-80d5-c10f07883f4e
# ╟─6dfe666b-3b7e-47fa-ac54-c1666cce044f
# ╟─dd801718-4a17-4eb2-969b-8cc291142f6d
# ╟─84c8f90b-0871-4086-91cf-973842daec47
# ╟─ee57e7db-0875-4df3-9f3b-db8e9157da31
# ╟─c752f559-a0ea-4d54-8e4c-e5e9ee7708f1
# ╟─7ac91382-2102-4cb2-b2c4-a2934591968b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
