__precompile__()
module Scallop

using Rayden, CoordinateTransformations, StaticArrays, OffsetArrays, Rotations

export Volume, Plane, getmembranes, scallop, Light, Near, Far, getLight, getRay!, getRay, getpsf, get2dpsf, getfwhm, coordinates2d, retinas

const retinas = ["distal_retina", "proximal_retina"]

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

function getmembranes(volumes::Dict{String, Volume}, planes::Dict{String, Plane}, morphz::Float64, aperture::Float64)
    morphxy = sqrt(1/morphz)
    squeeze = LinearMap(@SMatrix([morphxy 0 0; 0 morphxy 0; 0 0 morphz]))
    z_last = cumsum([-volumes[k].thickness for k in ["water", "cornea", "lens", "lens2distal_retina", "distal_retina2proximal_retina", "gap"]])
    return Dict(k => Ellipsoid(squeeze(Vec(0, 0, z_last_i - (-1)^planes[k].down*planes[k].radii[3])), squeeze(planes[k].radii), Vec(0, 0, (-1)^planes[k].down), k == "cornea_distal" ? getcosα(aperture/2, planes[k].radii[1], planes[k].radii[3]) : cospi(.3)) for (z_last_i, k) in zip(z_last, ["cornea_distal", "lens_distal", "lens_proximal", "distal_retina", "proximal_retina", "mirror"]))
end


function scallop(planes::Dict{String, Ellipsoid}, volumes::Dict{String, Volume})
    mind = ["cornea_distal", "lens_distal", "lens_proximal", "distal_retina", "proximal_retina", "mirror", "proximal_retina", "distal_retina"]
    vind1 = ["water", "cornea", "lens", "lens2distal_retina", "distal_retina2proximal_retina", "gap", "gap", "distal_retina2proximal_retina"]
    vind2 = ["cornea", "lens", "lens2distal_retina", "distal_retina2proximal_retina", "gap", "behind_mirror", "distal_retina2proximal_retina", "lens2distal_retina"]
    return [OpticUnit(planes[k], i > 6 ? planes[k].dir[3] > 0 : planes[k].dir[3] < 0, volumes[v1].ri/volumes[v2].ri, r"retina"(k), k) for (k, v1, v2, i) in zip(mind, vind1, vind2, 1:8)]
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
    if l < 1e6
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
        psf[k] = OffsetArray(zeros(Int, 2r + 1), -(r + 1))
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
                col = round(Int, r.orig[1])
                psf[ou.name][col] += 1
                psf[ou.name][-col] += 1
            end
        end
        if complete
            i += 1
        end
    end
    return psf
end

function getfwhm(ous::Vector{OpticUnit}, viewdistance::Float64, aperture::Float64, n::Int, δ::Float64)
    θ = Dict(k => Nullable{Float64}() for k in retinas)
    ratio = Dict(k => 0.0 for k in retinas)
    θ1 = 0.0
    while any(v < 2 for v in values(ratio))
        θ1 += δ
        l = getLight(viewdistance, aperture, ous[1].body, θ1)
        psf = get2dpsf(ous, l, n)
        for k in keys(ratio)
            ratio[k] = maximum(psf[k])/psf[k][0]
            if isnull(θ[k]) && ratio[k] ≥ 2
                θ[k] = Nullable(rad2deg(θ1 - δ/2))
            end
        end
    end
    return θ
end

# plots:

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



end # module
