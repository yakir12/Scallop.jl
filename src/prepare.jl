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

ri = RI(zeros(7)...)
thick = Thick(μm*zeros(7)...)
radii = Radii(μm*zeros(6)...)
