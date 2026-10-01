//- Update readable liquid-content diagnostics (theta, mLiquid) from the
//  current saturation field, and log the domain totals at write times.
theta = eps*Sb;
mLiquid.primitiveFieldRef() = rhob.value()*theta.primitiveField()*mesh.V();

if (runTime.writeTime())
{
    //- domainIntegrate is a global reduction: every rank must call it,
    //  otherwise the master deadlocks waiting for the other processors.
    scalar totalLiquidVolume = fvc::domainIntegrate(theta).value(); // m3
    scalar totalLiquidMass = rhob.value()*totalLiquidVolume;        // kg
    scalar meanSaturation =
        fvc::domainIntegrate(Sb).value()/gSum(mesh.V());

    if (Pstream::master())
    {
        liquidBalanceFile()
            << runTime.timeName() << ","
            << totalLiquidVolume*1e6 << ","
            << totalLiquidMass*1e3 << ","
            << meanSaturation << endl;
    }
}
