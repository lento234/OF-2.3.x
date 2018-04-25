pc=pc.prevIter(); 
forAll(pc.boundaryField(), patchi)
{ 
    forAll(pc.boundaryField()[patchi], facei)
    {
        pc.boundaryField()[patchi][facei] = pc.prevIter().boundaryField()[patchi][facei];
    }
}
Ts=Ts.prevIter();
forAll(Ts.boundaryField(), patchi)
{ 
    forAll(Ts.boundaryField()[patchi], facei)
    {
        Ts.boundaryField()[patchi][facei] = Ts.prevIter().boundaryField()[patchi][facei];
    }
}
