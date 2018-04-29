
namespace Foam
{

//- Return vapour pressure [Pa]
scalar soilVegetationModel::calc_pv(const scalar& p, const scalar& w)
{
    // vapor pressure pv
    return (p*w)/(0.6219438209 + w);
}

//- Return saturated water vapour pressure [Pa]
scalar soilVegetationModel::calc_pvsat(const scalar& T)
{
    // saturated vapor pressure pws - ASHRAE 1.2
    return exp( - 5.8002206e3/T
                + 1.3914993
                - 4.8640239e-2*T
                + 4.1764768e-5*pow(T,2)
                - 1.4452093e-8*pow(T,3)
                + 6.5459673*log(T) );
}

//- Return vapour pressure deficit [Pa/Pa]
scalar soilVegetationModel::calc_VPD(const scalar& T, const scalar& p, const scalar& w)
{
    // saturation vapour pressure [Pa]
    scalar pvs = calc_pvsat(T);
    // vapour pressure [Pa]
    scalar pv = calc_pv(p, w);

    return (pvs - pv)/Pstd_.value();
}


}