using Rayden, Scallop, Unitful
import Unitful: nm, μm, mm, cm, m, km, inch, ft, mi, rad, °
include("plottingFunctions.jl")

mutable struct RI 
    water
    cornea
    lens
    lens2distal_retina
    distal_retina2proximal_retina
    gap
    behind_mirror
end

mutable struct Thick
    water
    cornea
    lens
    lens2distal_retina
    distal_retina2proximal_retina
    gap
    behind_mirror
end

mutable struct Radii
    cornea_distal
    lens_distal
    lens_proximal
    distal_retina
    proximal_retina
    mirror
end

length_conv{T <: Number}(x::T) = Float64(ustrip(uconvert(μm, x)))
angle_conv{T <: Number}(x::T) = Float64(ustrip(uconvert(rad, x)))

function format_radii{T <: Number}(x::Tuple{T})
    xx1, = x
    x1 = length_conv(xx1)
    return Vec(x1, x1, x1)
end
function format_radii{T <: Number}(x::Tuple{T, T})
    xx1, xx2 = x
    x1 = length_conv(xx1)
    x2 = length_conv(xx2)
    return Vec(x1, x1, x2)
end
function format_radii{T <: Number}(xx::T) 
    x = length_conv(xx)
    return Vec(x, x, x)
end

ri = RI(zeros(7)...)
thick = Thick(μm*zeros(7)...)
radii = Radii(μm*zeros(6)...)

include("variables.jl")
thick.behind_mirror = Inf*μm
ri.behind_mirror = 0.0
thick.water = 0μm
morphing_factor = Float64(morphing_factor)
source_distance = length_conv(source_distance)

volumes = Dict(String(k) => Volume(Float64(getfield(ri, k)), length_conv(getfield(thick, k))) for k in fieldnames(ri))
planes = Dict(String(k) => Plane(format_radii(getfield(radii, k)), !r"distal$"(String(k))) for k in fieldnames(radii))
aperture = length_conv(aperture)
morphing_factor = morphing_factor

ellipsoids = getmembranes(volumes, planes, morphing_factor, aperture)
opticunits = scallop(ellipsoids, volumes)
l = getLight(source_distance, aperture, ellipsoids["cornea_distal"], angle_conv(source_angle))

min_aperture = length_conv(min_aperture)
max_aperture = length_conv(max_aperture)
min_distance = length_conv(min_distance)
max_distance = length_conv(max_distance)


ploteye(ellipsoids, l, opticunits)

# plotpixels(opticunits, l)
plotfwhm_aperture(opticunits, source_distance, min_aperture, max_aperture)
# plotfwhm_distance(opticunits, aperture, min_distance, max_distance)
# plotfwhm_aperture_distance(opticunits, min_aperture, max_aperture, min_distance, max_distance, n_rays = 1000)
