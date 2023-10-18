"""
    GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the GPP of the SPAC per ground area
"""
function GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _gpp::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _gpp += numerical∫(_iPS.Ag, _iPS.LAIx) * _iPS.LA;
    end;

    return _gpp / spac.ga
end


"""
    CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the canopy NPP of the SPAC per ground area
"""
function CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _cnpp::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _cnpp += numerical∫(_iPS.An, _iPS.LAIx) * _iPS.LA;
    end;

    return _cnpp / spac.ga
end


"""
    T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the transpiration of the SPAC per ground area
"""
function T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _t_veg::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        _iPS = spac.plant_ps[_i_can];
        _t_veg += numerical∫(_iPS.g_lw, _iPS.LAIx) * (_iPS.p_sat - _iEN.p_H₂O) / _iEN.p_atm * _iPS.LA;
    end;

    return _t_veg / spac.ga  
end


function Canopy_cond(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _g_lw::FT = 0;

    for _i_can in 1:spac.n_canopy
        # @show _i_can;
        _iPS = spac.plant_ps[_i_can];
        _g_lw += numerical∫(_iPS.g_lw, _iPS.LAIx);
    end;

    return  _g_lw / spac.ga
 
end

function gsw_ss_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    gsw_ss::FT = 0;

    for _i_can in 1:spac.n_canopy
        # @show _i_can;
        _iPS = spac.plant_ps[_i_can];
        gsw_ss= _iPS.gsw_ss
    end;

    return  gsw_ss
 
end


function g_sw_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    g_sw::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        g_sw += numerical∫(_iPS.g_sw, _iPS.LAIx);
    end;

    return   g_sw / spac.ga
 
 
end

function g_sw0_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    g_sw0::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        g_sw0 += numerical∫(_iPS.g_sw0, _iPS.LAIx);
    end;

    return   g_sw0 / spac.ga
 
 
end


function g_bw_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    g_bw::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        g_bw += numerical∫(_iPS.g_bw, _iPS.LAIx);
    end;

    return   g_bw / spac.ga
 
 
end


function An_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _An::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _An += numerical∫(_iPS.An, _iPS.LAIx);
    end;

    return _An/ spac.ga
end

function LAIx_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    LAIx::FT = 0.0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        LAIx += numerical∫(_iPS.LAIx, _iPS.LAIx)
    end;

    return LAIx
end


function LA_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    LA::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        LA =  _iPS.LA;
    end;

    return LA  
end

function tao_esm_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    τ_esm::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        τ_esm=  _iPS.τ_esm;
    end;

    return τ_esm  
end




function p_sat_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    p_sat::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        p_sat= _iPS.p_sat
    end;

    return p_sat
end

function vpd_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    vpd::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        _iPS = spac.plant_ps[_i_can];
        vpd = max(FT(0.001), _iPS.p_sat - _iEN.p_H₂O)
    end;

    return vpd
end



function p_H₂O_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    p_H₂O::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        p_H₂O = _iEN.p_H₂O
    end;

    return p_H₂O  
end

function p_atm_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    p_atm::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        p_atm=  _iEN.p_atm 
    end;

    return p_atm  
end

function p_a_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    p_a::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        p_a=  _iEN.p_a 
    end;

    return p_a 
end


function gamma_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    Γ_star::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        Γ_star=  _iPS.ps.Γ_star 
    end;

    return Γ_star 
end

function p_s_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    p_s::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        p_s+= numerical∫(_iPS.p_s, _iPS.LAIx);
    end;

    return p_s
end

function p_i_out(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    p_i::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        p_i+= numerical∫(_iPS.p_i, _iPS.LAIx);
    end;

    return p_i
end




"""
    PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the cumulative PPAR of the SPAC per ground area
"""
function PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _ppar::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _ppar += numerical∫(_iPS.APAR, _iPS.LAIx) * FT(1e-6) * spac.canopy_rt.LAI / spac.canopy_rt.nLayer;
    end;

    return _ppar
end