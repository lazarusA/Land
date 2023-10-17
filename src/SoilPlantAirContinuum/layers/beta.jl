"""
    spac_beta_max(spac::SPACMono{FT}, beta::BetaGLinearPsoil{FT}) where {FT<:AbstractFloat}

Compute the beta tuning factor for SPAC by taking the maximum, given
- `spac` SPAC
- `beta` `BetaGLinearPsoil` type scheme

"""
function spac_beta_max(spac::SPACMono{FT}, beta::BetaGLinearPsoil{FT}, return_swc::Bool = false) where {FT<:AbstractFloat}
    _βm::FT = 0;
    swc_value::FT = zero(FT)

    for _i in eachindex(spac.plant_hs.roots)
        _β = β_factor(spac.plant_hs.leaves[1], spac.plant_hs.roots[_i].sh, beta, FT(0), spac.plant_hs.roots[_i].p_ups, spac.swc[_i]);
        _βm = max(_β, _βm);
        swc_value = spac.swc[_i]
    end

    if return_swc
        return _βm, swc_value
    else
        return _βm
    end
end

