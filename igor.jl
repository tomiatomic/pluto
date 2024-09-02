### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ cdbe0240-4285-11ef-0a0b-2d0cfaf6901f
using PlutoUI

# ╔═╡ d7d289f3-b855-4b06-90da-8b1d5f940173
TableOfContents()

# ╔═╡ 5c4fa047-8d68-4255-9c43-4a8b546735b3
md"""## Pauli  (or [Clogston](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.9.266)–[Chandrasekhar](https://pubs.aip.org/aip/apl/article-abstract/1/1/7/976920/A-NOTE-ON-THE-MAXIMUM-CRITICAL-FIELD-OF-HIGH-FIELD)) paramagnetic limit:
(Darshana's preliminary intro to our [PRB '23 paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.L220501) from 22.7.2022.) \
In superconductors, an external magnetic field can break the Cooper pairs via two main mechanisms. The *orbital pair breaking* occurs when magnetic field acts on the Cooper pair’s electrons’ charge. The Cooper pair breaks when the induced Lorentz force exceeds its binding force. In 2D superconductors, this mechanism is absent for when the magnetic field is oriented in-plane. Nevertheless, in spin-singlet superconductors, the second mechanism – *spin pair breaking* – breaks the Cooper pair by aligning both its antiparallel spins with the magnetic field via the Zeeman effect. The minimum field necessary for spin pair breaking is known as the Pauli paramagnetic limit ``B_P``. For conventional superconductors with critical temperature ``T_c``, the critical magnetic field cannot exceed
#### ``B_P \approx 1.86 \times T_c ~~ [T/K]``. 
"""

# ╔═╡ 72bad491-2ad6-4440-9de7-f42b734e714b
md""" ## Qualitative reasoning
It can be understood from the simple energy balance: in the normal state, an external field ``B`` is Pauli-screened, with an energy gain of ``\frac{\chi_{Pauli}B^2}{2}``, compared to an absence of such screening, where ``\chi_{Pauli}`` is the Pauli susceptibility, which in atomic units, is the density of states at the Fermi level ``N(0)``. In the Bardeen-Cooper-Schrieffer theory the pairing energy is ``\frac{N(0)\Delta^2}{2}``, so above the thermodynamic critical field ``B_t=\Delta`` it is energetically favorable to break the superconducting state. In reality both relations are renormalized, the former via Stoner enhancement, the latter via strong coupling effects, but the qualitative physics remains. \
In SI units:
* ``\frac{2\Delta}{k_B T_c} \approx 3.53`` ...from [BCS theory](https://en.wikipedia.org/wiki/BCS_theory#Implications)
* ``\mu_B \approx 5.8 \times 10^{-5}`` is the [Bohr magneton](https://en.wikipedia.org/wiki/Bohr_magneton)
* ``k_B \approx 8.62 \times 10^{-5}`` is the [Boltzmann constant](https://en.wikipedia.org/wiki/Boltzmann_constant)
#### ``B_t \approx \frac{3.53}{2}\frac{k_B}{\mu_B}T_c \approx2.6 × T_c ~~[T/K]``, which is close to ``B_P``.
"""

# ╔═╡ fdb79c73-3777-423f-8f10-4e9766f87652
begin
critemp = @bind critemp NumberField(0.00:0.01:100.00, default = 6.0)
md"### ``T_c=`` $(critemp) K"
end

# ╔═╡ fe9d9b05-e71d-4815-8031-c34d3d4e0483
@bind reset Button("Reset")

# ╔═╡ 5f3634fc-2829-44cf-99dc-605d22474979
begin
	reset
	coupling = @bind coupling Slider(3.00:0.01:5.00, 3.53, true)
	md"### BCS coupling: $(coupling)"
end

# ╔═╡ 9d0d7817-01f9-41e2-aba2-056be056bce5
md"## Pauli limit is $(round(1.86*coupling*critemp/3.53; digits = 2)) T
 ``B_P \approx 1.86\times \frac{coupling}{3.53} \times T_c ~~[T/K]``"

# ╔═╡ 035aa9e5-e159-4b83-be1e-de7427772a76
md"# ``4H-NbSe_2``:
- ``T_c \approx 6 K`` (SS2 crystal4)
- splitting ca. 150 meV between spin-up and spin-down (4 pairs, ``\Gamma - K - M``)``\implies \zeta = 75~meV``
- bandwidth in ``k_z`` direction for ``k`` in-plane at Fermi level crossing is ca. 17 meV ``\implies t \approx 4.25~ meV``
- offset (interlayer potential) is ca. 56 meV ``\implies T_0 \approx 28~meV``
"

# ╔═╡ 37801fb2-ac6b-484c-803f-25a278d56a0b
md"""
![Hc2 from Heat capacity](https://github.com/tomiatomic/pics/blob/main/C_P.jpg?raw=true)
"""

# ╔═╡ 7b630861-efab-42a3-969d-7084d07e28f2
md"""
![DFT](https://github.com/tomiatomic/pics/blob/main/dft.png?raw=true)
"""

# ╔═╡ 64566669-b988-4c20-9c0b-663ecbb70357
md"# Decoupled bilayer"

# ╔═╡ 8ff8de6a-7896-4351-8ec9-6473e138a143
md"## ``B_c = B_P\frac{\sqrt{\zeta^2+t^2}}{t}``
- ``\zeta`` is half of the spin-orbit splitting \
- ``t`` is the interlayer hopping: 1/4 of the bandwidth along the out-of-plane direction \
- and ``1.86[T/K]*T_c`` is the Pauli limit"

# ╔═╡ f69f4293-0d1b-4379-9354-4f1d1370527c
function bilayer(z,t,temp)
	@. (sqrt(z^2+t^2)/t)*1.86*temp
end

# ╔═╡ 33fe0602-e8ca-4a9c-a723-d523b4a61799
md"1T1H values: ``\zeta = 75 meV``; ``t \approx 25 meV``; ``T_c = 1.23 K``"

# ╔═╡ afead579-8f54-434d-b81a-9bb9932572a7
begin
	zeta = @bind zeta NumberField(0.00:100.00, default = 75.00)
	hop = @bind hop NumberField(0.00:0.05:100.00, default = 4.25)
	md"### ``\zeta=`` $(zeta) meV; ``t=`` $(hop) meV"
end

# ╔═╡ f0c35575-1f53-48e7-b8e5-b9cc59ed5fd5
md"## ``B_c =`` $(round(bilayer(zeta,hop,critemp); digits = 2)) ``T``; ``\frac{B_c}{B_P} \approx`` $(round(bilayer(zeta,hop,critemp)/(critemp*1.86); digits = 2))"

# ╔═╡ b570423d-8afb-44cf-b654-11342985c872
md"# Bilayer in asymmetric environment"

# ╔═╡ 1bb2d710-39b4-47e8-be01-bd1bb7cc446a
md"## ``\mu_BB_t = \sqrt{t_0^2+(T_0+\zeta)^2}-\sqrt{t_0^2+(T_0-\zeta)^2}``
- ``\zeta`` is half of the spin-orbit splitting \
- ``t_0 = 2\zeta`` is the interlayer coupling as in bulk ``NbSe_2`` \
- ``T_0`` is half of the mutual shift of the spin-split bands\
- ``\mu_B \approx 5.8×10^{−2}meV/T``"

# ╔═╡ af20d126-ce8e-4a46-a3e6-8304e227eafa
function asym(z,t)
	bohr = 5.7883818060e-2
	@. (sqrt((2*z)^2+(t+z)^2)-sqrt((2*z)^2+(t-z)^2))/bohr
end

# ╔═╡ 18ad4fdd-89e7-445d-ae37-c269033211fb
md"1T2H values: ``\zeta = 55~meV``; ``T_0 \approx 10~meV``"

# ╔═╡ d05fb8ca-0a99-4bab-b44d-8e285fe514ec
begin
	zetasym = @bind zetasym NumberField(0.00:100.00, default = 75.00)
	shift = @bind shift NumberField(0.00:0.01:100.00, default = 28.00)
	md"### ``\zeta=`` $(zetasym) meV; ``T_0=`` $(shift) meV"
end

# ╔═╡ 27c4de0b-47c6-46a3-8f17-2f8582cd1872
md"## ``B_c =`` $(round(asym(zetasym,shift); digits = 2)) ``T``; ``\frac{B_c}{B_P} \approx`` $(round(asym(zetasym,shift)/(critemp*1.86); digits = 2))"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "6e7bcec4be6e95d1f85627422d78f10c0391f199"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─cdbe0240-4285-11ef-0a0b-2d0cfaf6901f
# ╟─d7d289f3-b855-4b06-90da-8b1d5f940173
# ╟─5c4fa047-8d68-4255-9c43-4a8b546735b3
# ╟─72bad491-2ad6-4440-9de7-f42b734e714b
# ╟─fdb79c73-3777-423f-8f10-4e9766f87652
# ╟─5f3634fc-2829-44cf-99dc-605d22474979
# ╟─fe9d9b05-e71d-4815-8031-c34d3d4e0483
# ╟─9d0d7817-01f9-41e2-aba2-056be056bce5
# ╟─035aa9e5-e159-4b83-be1e-de7427772a76
# ╟─37801fb2-ac6b-484c-803f-25a278d56a0b
# ╟─7b630861-efab-42a3-969d-7084d07e28f2
# ╟─64566669-b988-4c20-9c0b-663ecbb70357
# ╟─8ff8de6a-7896-4351-8ec9-6473e138a143
# ╟─f69f4293-0d1b-4379-9354-4f1d1370527c
# ╟─33fe0602-e8ca-4a9c-a723-d523b4a61799
# ╟─afead579-8f54-434d-b81a-9bb9932572a7
# ╟─f0c35575-1f53-48e7-b8e5-b9cc59ed5fd5
# ╟─b570423d-8afb-44cf-b654-11342985c872
# ╟─1bb2d710-39b4-47e8-be01-bd1bb7cc446a
# ╟─af20d126-ce8e-4a46-a3e6-8304e227eafa
# ╟─18ad4fdd-89e7-445d-ae37-c269033211fb
# ╟─d05fb8ca-0a99-4bab-b44d-8e285fe514ec
# ╟─27c4de0b-47c6-46a3-8f17-2f8582cd1872
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
