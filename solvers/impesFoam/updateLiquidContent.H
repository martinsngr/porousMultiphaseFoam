//- Update readable liquid-content diagnostics (theta, mLiquid) from the
//  current saturation field, and log the domain totals at write times.
theta = eps*Sb;
mLiquid.primitiveFieldRef() = rhob.value()*theta.primitiveField()*mesh.V();

if (runTime.writeTime() && Pstream::master())
{
    scalar totalLiquidVolume = fvc::domainIntegrate(theta).value(); // m3
    scalar totalLiquidMass = rhob.value()*totalLiquidVolume;        // kg
    liquidBalanceFile()
        << runTime.timeName() << ","
        << totalLiquidVolume*1e6 << ","
        << totalLiquidMass*1e3 << endl;
}
