using Rayden, Scallop, Plots
gr()

volumes = Dict(
"water"                         => Volume(1.34, 0),
"cornea"                        => Volume(1.37, 23),
"lens"                          => Volume(1.42, 215),
"lens2distal_retina"            => Volume(1.35, 6),
"distal_retina2proximal_retina" => Volume(1.35, 81),
"gap"                           => Volume(1.34, 49),
"behind_mirror"                 => Volume(0.0,  Inf)
)
planes = Dict(
"cornea_distal"   => Plane(cospi(.3), Vec(244, 244, 227), false),
"lens_distal"     => Plane(cospi(1/3),  Vec(195, 195, 213), false),
"lens_proximal"   => Plane(cospi(.3), Vec(337, 337, 337), true),
"distal_retina"   => Plane(cospi(.3), Vec(337, 337, 337), true),
"proximal_retina" => Plane(cospi(.3), Vec(337, 337, 337), true),
"mirror"          => Plane(cospi(.3), Vec(337, 337, 337), true)
)

ellipsoids = getmembranes(volumes, planes, 1.04)# 1.157)
opticunits = scallop(ellipsoids, volumes)
aperture = 240.

l = getLight(1e9, aperture, ellipsoids["cornea_distal"], 0.)
x, y = coordinates2d(ellipsoids)
ph = plot(x,y, aspect_ratio = :equal, leg=false, color = :black)
r = Ray()
# plot!([-aperture/2, aperture/2], [radius2z_shift(aperture/2, planes["cornea_distal"]), radius2z_shift(aperture/2, planes["cornea_distal"])])#, marker = :circle, line = :red)
for b in linspace(0, 1, 25), b2 in [0., 0.5]
    getRay!(r, l, b, b2)
    x = Float64[]
    y = Float64[]
    push!(x, r.orig[1])
    push!(y, r.orig[3])
    for i in opticunits
        raytrace!(r, i)
        push!(x, r.orig[1])
        push!(y, r.orig[3])
        # p = i.body.c + i.body.dir[3]*i.body.r
        # annotate!(0, p[3], text(i.name, 10, :black, :center))
    end
    plot!(x, y, color = :grey)
end
plot!(ylims = (-500,100), xlims = (-300,300))

psf = getpsf(opticunits, l, 100000)
plot(heatmap(collect(psf["distal_retina"][-5:5,-5:5]), aspect_ratio=:equal), heatmap(collect(psf["proximal_retina"][-110:110, -110:110]), aspect_ratio=:equal))


l = getLight(1e9, aperture, ellipsoids["cornea_distal"], .05)
psf = get2dpsf(opticunits, l, 100000)
plot()
for (k, v) in psf
    plot!(collect(v), label = k)
end
gui()


nt = 20
thetas = linspace(0, .5, nt)
ratios = Dict(k => zeros(nt) for k in retinas)
for (i, θ) in enumerate(thetas)
    l = getLight(1e9, aperture, ellipsoids["cornea_distal"], θ)
    psf = get2dpsf(opticunits, l, 100000)
    for k in keys(ratios)
        ratios[k][i] = maximum(psf[k])/psf[k][0]
    end
end
plot(ylims = (0,5))
for (k, v) in ratios
    plot!(thetas, v, label = k)
end
gui()

n = 10
apertures = linspace(100, 300, 10)
fwhms = Dict(k => zeros(n) for k in retinas)
for (i, aperture) in enumerate(apertures), (k, v) in getfwhm(opticunits, 1e9, aperture, 10000, deg2rad(1/4))
    fwhms[k][i] = get(v)
end
plot(ylims = (0,20))
for (k, v) in fwhms
    plot!(apertures, v, label = k)
end
gui()

aperture = 240.
n = 10
viewdistances = logspace(3, 9, 10)
fwhms = Dict(k => zeros(n) for k in retinas)
for (i, viewdistance) in enumerate(viewdistances), (k, v) in getfwhm(opticunits, viewdistance, aperture, 10000, deg2rad(1/4))
    fwhms[k][i] = get(v)
end
plot()
for (k, v) in fwhms
    plot!(viewdistances, v, label = k)
end
plot!(ylims = (0,22), xscale = :log10)

n = 10
viewdistances = logspace(3, 9, 10)
apertures = linspace(100, 300, 10)
fwhms = Dict(k => zeros(n,n) for k in retinas)
for (i, viewdistance) in enumerate(viewdistances), (j, aperture) in enumerate(apertures), (k, v) in getfwhm(opticunits, viewdistance, aperture, 10000, deg2rad(1/4))
    fwhms[k][i,j] = get(v)
end

ph = []
for (k, v) in fwhms
    push!(ph, heatmap(apertures, viewdistances, v, label = k))
end
plot(ph..., zlims = (0,10))
