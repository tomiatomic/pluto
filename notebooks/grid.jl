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
	using PlutoUI, PlutoPlotly, JLD2
	md"`packages`"
end

# ╔═╡ 6a6ea7ad-0a03-464d-8b84-4f642db7cfdd
TableOfContents()

# ╔═╡ 6df8df70-3346-11f1-a692-31d015c5335e
begin
	source = @bind source Select([1 => "WSxM", 2 => "Nanonis"])
md"# Process & plot spectral maps
using [PlutoPlotly](https://github.com/JuliaPluto/PlutoPlotly.jl)\
Select file source: $(source) ...under construction...

process:
- lowess switch (or zero) - smooth_iv
- NRD (test on clicked spectra and then apply to all witch checkbox)
- add dIdV/I/V plot
- add Matlab plots
- add profile along multiple clicks

misc:
- check current values
- add Nanonis
"
end

# ╔═╡ 89ec1e0d-3811-4ef7-9bbc-4a0f9d8d41b4
md"""## Import data
**drag & drop** grid file: $(@bind gsi FilePicker()) (e.g."...\Pluto\STM\datatest\test.gsi")
"""

# ╔═╡ c0f60871-fde9-429a-8389-f014844cd2a2
md"## Plot data"

# ╔═╡ 557407a7-fcda-48a6-9a60-2f7626e579dc
md"### Topography
process in *e.g.* [WSxM](https://wsxm.eu/) or [Gwyddion](https://gwyddion.net/)"

# ╔═╡ dd53ab84-5693-4c87-bb41-710841cc02aa
@bind unclick Button("unclick")

# ╔═╡ 181f2fab-cc28-4016-b9a9-ce0db79ef2ab
md"### *I(V)* spectrum @ clicked pixel
[1,1] if no click"

# ╔═╡ 1006d8dd-738c-42ff-b275-7521d2ebde33
md"###  *I(V)* slice heatmap"

# ╔═╡ 049c19cd-91a6-4459-aa2e-55c1089e1b1d
md"""# Hic sunt functiones programmatoriae"""

# ╔═╡ 45e0787a-c377-4c84-b639-aeabf66e946f
md"## Import .gsi file"

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
	volts = unsorted_bias[perm]
	iv = unsorted_iv[:, :, perm]
		
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
	curve = vec(iv[row, col, :])
	Plot(
	    scatter(x=volts, y=curve, mode="lines", line=attr(width=2)),
	    Layout(
	        title="IV curve at (row=$(row), col=$(col))",
	        xaxis=attr(title="bias [mv]"),
	        yaxis=attr(title="Signal"),
	        width=600, height=400
	    	)
		)
end

# ╔═╡ 8f73478c-1edf-4a44-8145-71e303d5d95e
println("position: x = ", x_range[col], " y = ", y_range[row])

# ╔═╡ 110c7446-f4ce-47d5-8cab-02aa293de31b
let
    # IMPORTANT: typed vector so Plot(traces, ...) works
    traces = PlotlyBase.AbstractTrace[]

    push!(traces,
        heatmap(
            z = topo,
            x = x_range,
            y = y_range,
            colorscale = "Viridis",
            colorbar = attr(title="z [$(z_unit)]")
        )
    )

    if !(ismissing(click) || click === nothing)
        # If your click is a Dict with "row"/"col" keys:
        push!(traces,
            scatter(
                x = [x_range[col]],                 # must be arrays
                y = [y_range[row]],
                mode = "markers",
                marker = attr(size=12, color="red", symbol="x"),
                name = "selected"
            )
        )
    end

    plot(
        traces,
        Layout(
            title = "Topography (selected pixel marked)",
            xaxis = attr(title="x [$(topo_unit)]"),
            yaxis = attr(title="y [$(topo_unit)]", scaleanchor="x", scaleratio=1),
            width = 650, height = 650,
            uirevision = "keep"
        )
    )
    end

# ╔═╡ 46efaead-594d-49ec-abc4-97b4db9c1eb3
@bind bias Slider(eachindex(volts); default = 1, show_value = i -> "$(volts[i]) mV")

# ╔═╡ 384426f4-fea3-41a1-9328-d486eb7ba3f1
plot(
        heatmap(x=x_range, y=y_range, z=iv[:, :, bias], colorscale="RdBu", reversescale = true, colorbar=attr(title="current")),
        Layout(
            title="I(V) map",
            xaxis=attr(title="x [$(topo_unit)]"),
            yaxis=attr(title="y [$(topo_unit)]", scaleanchor="x", scaleratio=1),
            width=650, height=650
        )
    )

# ╔═╡ 60993909-da30-49de-8ee6-3875a2df2174
begin
	tr = surface(x = x_range, y = y_range, z = topo,
	    surfacecolor = iv[:, :, bias],     # <-- heatmap painted on surface
	    colorscale = "Viridis",
	    cmin = minimum(iv[:, :, bias]),
	    cmax = maximum(iv[:, :, bias]),
	    colorbar = attr(title="IV value"),
	    showscale = true
	)
	
	lay = Layout(
	    title = "Surface colored by heatmap (surfacecolor) @ $(volts[bias]) mV",
	    scene = attr(
	        xaxis = attr(title="x [$(topo_unit)]"),
	        yaxis = attr(title="y [$(topo_unit)]"),
	        zaxis = attr(title="z [$(z_unit)]")
	    ),
	    width = 600, height = 400, uirevision = "keep"
	)
	
	plot(tr, lay)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
JLD2 = "~0.6.3"
PlutoPlotly = "~0.6.5"
PlutoUI = "~0.7.80"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "b9af1ee716e1c23f3cb15502bd55f3a91d54ba8f"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

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

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6522cfb3b8fe97bec632252263057996cbd3de20"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.18.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLD2]]
deps = ["ChunkCodecLibZlib", "ChunkCodecLibZstd", "FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues"]
git-tree-sha1 = "8f8ff711442d1f4cfc0d86133e7ee03d62ec9b98"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.6.3"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

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

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

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

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

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
git-tree-sha1 = "fbc875044d82c113a9dee6fc14e16cf01fd48872"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.80"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

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

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "fa95b3b097bcef5845c142ea2e085f1b2591e92c"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.7.1"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

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

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─a5540242-d382-4ac5-87b9-c0c45f74cef3
# ╟─6a6ea7ad-0a03-464d-8b84-4f642db7cfdd
# ╟─6df8df70-3346-11f1-a692-31d015c5335e
# ╟─89ec1e0d-3811-4ef7-9bbc-4a0f9d8d41b4
# ╟─c9856d06-9463-45eb-a686-f0b51b5d3629
# ╟─05594a29-1279-4428-80cb-940268e60e6b
# ╟─c0f60871-fde9-429a-8389-f014844cd2a2
# ╟─557407a7-fcda-48a6-9a60-2f7626e579dc
# ╟─dd53ab84-5693-4c87-bb41-710841cc02aa
# ╟─b1aa9714-624e-4940-a57a-cc2dacf8a6ff
# ╟─181f2fab-cc28-4016-b9a9-ce0db79ef2ab
# ╟─00d9b5f3-1d87-4443-aaca-7a20feefc016
# ╟─50474186-101b-45de-9540-cd0bedc75323
# ╟─8f73478c-1edf-4a44-8145-71e303d5d95e
# ╟─110c7446-f4ce-47d5-8cab-02aa293de31b
# ╟─1006d8dd-738c-42ff-b275-7521d2ebde33
# ╟─46efaead-594d-49ec-abc4-97b4db9c1eb3
# ╟─384426f4-fea3-41a1-9328-d486eb7ba3f1
# ╟─60993909-da30-49de-8ee6-3875a2df2174
# ╟─049c19cd-91a6-4459-aa2e-55c1089e1b1d
# ╟─45e0787a-c377-4c84-b639-aeabf66e946f
# ╟─20c7f0e7-9be5-4922-9209-f5bbd03eae98
# ╟─7b9eb7f5-cfe6-40c6-a093-0e9f20a8a19c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
