__precompile__()
module Scallop

using Rayden, CoordinateTransformations, StaticArrays, OffsetArrays, Rotations, Optim, Cuba#ture

export Volume, Plane, getmembranes, scallop, Light, Near, Far, getLight, getRay!, getRay, getpsf, get2dpsf, getfwhm, coordinates2d, coordinates3d, retinas, retinai

const retinas = ["proximal_retina", "distal_retina"]
const retinai = Dict(k => i for (i, k) in enumerate(retinas))

immutable Volume
    ri::Float64
    thickness::Float64
end

immutable Plane
    radii::Vec
    down::Bool # is the dome facing down?
end

getz(x::Float64, rx::Float64, rz::Float64) = rz*sqrt(1 - (x/rx)^2)
function getcosα(x::Float64, rx::Float64, rz::Float64)
    X = getz(x, rx, rz)/x
    return X/sqrt(1 + X^2)
end

function getmembranes(volumes::Dict{String, Volume}, planes::Dict{String, Plane}, aperture::Float64)
    z_last = cumsum([-volumes[k].thickness for k in ["water", "cornea", "lens", "lens2distal_retina", "distal_retina2proximal_retina", "gap"]])
    # return Dict(k => Ellipsoid(Vec(0, 0, z_last_i - (-1)^planes[k].down*planes[k].radii[3]), planes[k].radii, Vec(0, 0, (-1)^planes[k].down), k == "ksjfhksdljfhlskfhlka" ? getcosα(aperture/2, planes[k].radii[1], planes[k].radii[3]) : 0.0) for (z_last_i, k) in zip(z_last, ["cornea_distal", "lens_distal", "lens_proximal", "distal_retina", "proximal_retina", "mirror"]))
    return Dict(k => Ellipsoid(Vec(0, 0, z_last_i - (-1)^planes[k].down*planes[k].radii[3]), planes[k].radii, Vec(0, 0, (-1)^planes[k].down), k == "cornea_distal" ? getcosα(aperture/2, planes[k].radii[1], planes[k].radii[3]) : 0.0) for (z_last_i, k) in zip(z_last, ["cornea_distal", "lens_distal", "lens_proximal", "distal_retina", "proximal_retina", "mirror"]))
end


function scallop(planes::Dict{String, Ellipsoid}, volumes::Dict{String, Volume}, register::String)
    mind = ["cornea_distal", "lens_distal", "lens_proximal", "distal_retina", "proximal_retina", "mirror", "proximal_retina", "distal_retina"]
    vind1 = ["water", "cornea", "lens", "lens2distal_retina", "distal_retina2proximal_retina", "gap", "gap", "distal_retina2proximal_retina"]
    vind2 = ["cornea", "lens", "lens2distal_retina", "distal_retina2proximal_retina", "gap", "behind_mirror", "distal_retina2proximal_retina", "lens2distal_retina"]
    return [OpticUnit(planes[k], i > 6 ? planes[k].dir[3] > 0 : planes[k].dir[3] < 0, volumes[v1].ri/volumes[v2].ri, contains(k, register) && i > 6, k) for (k, v1, v2, i) in zip(mind, vind1, vind2, 1:8)]
end

abstract type Light end

immutable Near <: Light
    cosa::Float64
    orig::Vec
    rotm::LinearMap{RotY{Float64}}
end

immutable Far <: Light
    R::Float64
    z::Float64
    dir::Vec
    rotm::AffineMap{RotY{Float64}, Vec}
end

function getLight(l::Float64, a::Float64, s::Ellipsoid, θ::Float64)
    r = a/2
    h = sqrt((1 - r^2/s.r[1]^2)*s.r[3]^2) + s.c[3]
    rotm = LinearMap(RotY(θ))
    tran = recenter(rotm, Vec(0,0,h))
    if l < 1e7
        orig = tran(Vec(0,0,l))
        p2 = Vec(r,0,h) - orig 
        cosalpha1 = dot(normalize(-orig), normalize(p2))
        p2 = Vec(-r,0,h) - orig
        cosalpha2 = dot(normalize(-orig), normalize(p2))
        cosalpha3 = cos(atan(r/(l - h)))
        cosa = minimum([cosalpha1, cosalpha2, cosalpha3])
        Near(cosa, orig, rotm)
    else
        Far(r, 100.0, rotm(Vec(0,0,-1)), tran)
    end
end

#=function Ray(ray::Ray, l::NearLight, b1::Float64, b2::Float64)
    z = b1*(l.cosa - 1) - l.cosa
    k = sqrt(1 - z^2)
    theta = 2π*b2
    ray.orig = l.orig
    ray.dir = Vec(k*cos(theta), k*sin(theta), z)
end=#

function getRay!(r::Ray, l::Near, b1::Float64, b2::Float64)
    z = b1*(l.cosa - 1) - l.cosa
    k = sqrt(1 - z^2)
    theta = 2b2
    r.orig = l.orig
    r.dir = l.rotm(Vec(k*cospi(theta), k*sinpi(theta), z))
end

function getRay!(r::Ray, l::Far, b1::Float64, b2::Float64)
    k = l.R*sqrt(b1)
    theta = 2b2
    r.orig = l.rotm(Vec(k*cospi(theta), k*sinpi(theta), l.z))
    r.dir = l.dir
end

function getRay{T <: Light}(l::T, b1::Float64, b2::Float64)
    r = Ray()
    getRay!(r, l, b1, b2)
    return r
end

#=function Ray(ray::Ray, l::FarLight, b1::Float64, b2::Float64)
    k = l.R*sqrt(b1)
    theta = 2π*b2
    ray.orig = Vec(k*cos(theta), k*sin(theta), l.z)
    ray.dir = l.dir
end=#

function getpsf{L <: Light}(ous::Vector{OpticUnit}, l::L, n::Int)
    psf = Dict{String, OffsetArray{Int, 2, Matrix{Int}}}()
    for k in retinas
        i = findfirst(x -> x.name == k, ous)
        r = ceil(Int, ous[i].body.r[1])
        psf[k] = OffsetArray(zeros(Int, 2r + 1, 2r + 1), -(r + 1), -(r + 1))
    end
    i = 0
    r = Ray()
    while i < n 
        getRay!(r, l, rand(), rand())
        complete = true
        for ou in ous
            failiure = raytrace!(r, ou)
            if failiure
                complete = false
                break
            end
            if ou.register 
                x = round(Int, r.orig[1])
                y = round(Int, r.orig[2])
                psf[ou.name][x, y] += 1
            end
        end
        if complete
            i += 1
        end
    end
    return psf
end

function get2dpsf{L <: Light}(ous::Vector{OpticUnit}, l::L, n::Int)
    psf = Dict{String, OffsetArray{Int,1,Array{Int,1}}}()
    for k in retinas
        i = findfirst(x -> x.name == k, ous)
        r = ceil(Int, ous[i].body.r[1])
        psf[k] = OffsetArray(zeros(Int, r + 1), -1)
    end
    i = 0
    r = Ray()
    while i < n 
        getRay!(r, l, rand(), rand())
        complete = true
        for ou in ous
            failiure = raytrace!(r, ou)
            if failiure
                complete = false
                break
            end
            if ou.register 
                col = abs(round(Int, r.orig[1]))
                psf[ou.name][col] += 1
            end
        end
        if complete
            i += 1
        end
    end
    return psf
end

#=function getfwhm(ous::Vector{OpticUnit}, viewdistance::Float64, aperture::Float64, n::Int, δ::Float64)
    θ = Dict(k => Nullable{Float64}() for k in retinas)
    ratio = Dict(k => 0.0 for k in retinas)
    θ1 = 0.0
    while any(v < 2 for v in values(ratio))
        θ1 += δ
        l = getLight(viewdistance, aperture, ous[1].body, θ1)
        psf = get2dpsf(ous, l, n)
        for k in keys(ratio)
            # ratio[k] = maximum(psf[k])/psf[k][0]
            ratio[k] = quantile(collect(psf[k]), .999)/psf[k][0]
            # ratio[k] = quantile(collect(psf[k]), .999)/mean(collect(psf[k][-1:1]))
            if isnull(θ[k]) && ratio[k] ≥ 2
                θ[k] = Nullable(rad2deg(θ1 - δ/2))
            end
        end
    end
    return Dict(k => get(v) for (k,v) in θ)
end=#

function centerpixel{L <: Light}(ous::Vector{OpticUnit}, l::L, θ::Float64, ψ::Float64, r::Ray, v::Vector{Float64})
    v[1] = 0.0
    getRay!(r, l, θ, ψ)
    for ou in ous
        raytrace!(r, ou) && break
        if ou.register
            if r.orig[1]^2 + r.orig[2]^2 < 30
                v[1] += 0.5
            end
        end
    end
    return v
end

function shortest_distance2OA{L <: Light}(ous::Vector{OpticUnit}, l::L, θ::Float64, ψ::Float64, r::Ray)
    v = Inf
    getRay!(r, l, θ, ψ)
    for ou in ous
        raytrace!(r, ou) && break
        if ou.register
            d = r.orig[1]^2 + r.orig[2]^2
            if d < v
                v = d
            end
        end
    end
    return sqrt(v)
end

function centerpixel{L <: Light}(ous::Vector{OpticUnit}, l::L, θ::Float64, ψ::Float64, r::Ray)
    v = 0.0
    getRay!(r, l, θ, ψ)
    for ou in ous
        raytrace!(r, ou) && break
        if ou.register
            if r.orig[1]^2 + r.orig[2]^2 < 6
                v += 0.5
            end
        end
    end
    return v
end


function sumcenterpixel{L <: Light}(ous::Vector{OpticUnit}, l::L)
    r = Ray()
    # val,err = hcubature(x -> centerpixel(ous, l, x[1], x[2], r), [0.,0.], [1.,1.], abstol = 1e-5)
    val,_ = divonne((x, y) -> centerpixel(ous, l, x[1], x[2], r, y), 2, 1, abstol = 1e-5)
    # val = sum(centerpixel(ous, l, θ, ψ, r) for θ in linspace(0, 1, 100) for ψ in linspace(0, 1, 100))
    # return val
    return val[1]
end

function sumcenterpixel_theta(distance::Float64, aperture::Float64, ous::Vector{OpticUnit}, theta::Float64)
    l = getLight(distance, aperture, ous[1].body, theta)
    return sumcenterpixel(ous, l)
end

function getfwhm(distance::Float64, aperture::Float64, volumes::Dict{String, Volume}, planes::Dict{String,Plane}, retina::String, alpha_guess::Float64)
    ellipsoids = getmembranes(volumes, planes, aperture)
    ous = scallop(ellipsoids, volumes, retina)
    target = sumcenterpixel_theta(distance, aperture, ous, 0.0)/2
    max_α = alpha_guess
    δ = .2
    max_α += δ
    min_α = alpha_guess
    y = sumcenterpixel_theta(distance, aperture, ous, max_α)
    while y > target
        min_α += δ
        max_α += δ
        y = sumcenterpixel_theta(distance, aperture, ous, max_α)
    end
    result = optimize(theta -> abs(target - sumcenterpixel_theta(distance, aperture, ous, theta)), min_α, max_α)
    return 2rad2deg(Optim.minimizer(result))
end

# plots:
function coordinates3d(s::Ellipsoid, n::Int)
    x = NaN*zeros(n,n)
    y = NaN*zeros(n,n)
    z = NaN*zeros(n,n)
    for (i, θ) in enumerate(linspace(0,π,n)), (j, ψ) in enumerate(linspace(0,π,n))
        p = Vec(cos(θ)*sin(ψ), sin(θ)*sin(ψ), cos(ψ))
        p = s.unscale(p)
        cosα = s.dir⋅normalize(p)
        if cosα > s.open
            p = s.uncenter(p)
            x[i, j] = p[1]
            y[i, j] = p[2]
            z[i, j] = p[3]
        end
    end
    return (x, y, z)
end

function coordinates2d(s::Ellipsoid, n::Int)
    α = acos(s.open)
    θs = linspace(-α + π/2, α + π/2, n)
    x = [s.c[1] + s.dir[3]*s.r[1]*cos(θ) for θ in θs]
    y = [s.c[3] + s.dir[3]*s.r[3]*sin(θ) for θ in θs]
    return (x, y)
end
function coordinates2d(planes::Dict{String, Ellipsoid})
    n = 25
    m = length(planes)
    x = zeros(n, m)
    y = zeros(n, m)
    for (i, v) in enumerate(values(planes))
        x[:,i], y[:,i] = coordinates2d(v, n)
    end
    return (x, y)
end

function fwhmstats(distance::Float64, aperture::Float64, volumes::Dict{String, Volume}, planes::Dict{String,Plane}, n::Int)
    ellipsoids = getmembranes(volumes, planes, aperture)
    ous = scallop(ellipsoids, volumes, "retina")
    l = getLight(distance, aperture, ous[1].body, 0.0)
    r = Ray()
    sigma = Dict(retina => 0.0 for retina in retinas)
    for i = 1:n
        getRay!(r, l, rand(), rand())
        # first = Dict(retina => true for retina in retinas)
        for ou in ous
            raytrace!(r, ou) && break
            if ou.register
                #=if first[ou.name]
                    first[ou.name] = false
                    continue
                end=#
                sigma[ou.name] += r.orig[1]^2 + r.orig[2]^2
            end
        end
    end
    fun = length2angle(distance, aperture, volumes, planes)
    return Dict(k => rad2deg(2*sqrt(2*log(2))*fun[k]*sqrt(v/n)) for (k, v) in sigma)
end

function length2angle(distance::Float64, aperture::Float64, volumes::Dict{String, Volume}, planes::Dict{String,Plane})
    r = Ray()
    fun = Dict(retina => 0.0 for retina in retinas)
    for retina in retinas
        ellipsoids = getmembranes(volumes, planes, aperture)
        ous = scallop(ellipsoids, volumes, retina)
        n = 100
        x = linspace(1e-6,.5,n)
        y = zeros(n)
        for (i, xi) in enumerate(x)
            l = getLight(distance, aperture, ous[1].body, xi)
            y[i] = shortest_distance2OA(ous, l, 1.,.0,r)
        end
        a = (y\collect(x))[1]
        fun[retina] = a
    end
    return fun
end






end # module
