#include "buildingMaterialModel.H"

forAll(Materials, MaterialsI)
{
	const dictionary& dict = Materials[MaterialsI];

    const word cellZoneName(dict.lookup("name"));
    const word cellZoneModel(dict.lookup("buildingMaterialModel"));

    //create the buildingMaterial model
    autoPtr<buildingMaterialModel> buildingMaterial
    (
       buildingMaterialModel::New("buildingMaterial", solidTransportProperties[i], cellZoneModel)
    );

    label cellZoneID = mesh.cellZones().findZoneID(cellZoneName);
    if (cellZoneID == -1)
    {
		Info << "Something is wrong, cannot find at least one of the necessary cellZones! Exiting!" << endl;
    	Foam::FatalError();
    }
    const labelList& cells = mesh.cellZones()[cellZoneID];
    //const polyBoundaryMesh& patches = mesh.boundaryMesh();
    boolList zoneCell(mesh.nCells(), false);


    forAll(cells, cellsI)
    {
        label celli = cells[cellsI];
        cellType[celli] = MaterialsI;
        zoneCell[celli] = true;

        if (update_sequence == 1)
        {
            // buildingMaterial->update_w_C_cell(pc,ws,Crel,celli);
            // buildingMaterial->update_Krel_cell(pc,ws,Krel,celli);
            // buildingMaterial->update_Kv_cell(pc,ws,Ts,K_v,celli);
            // buildingMaterial->update_Kpt_cell(pc,ws,Ts,K_pt,celli);
						Info << "Not implemtented !!" << endl;
        }
        else if (update_sequence == 2)
        {
            rho_m.internalField()[celli] = buildingMaterial->rho().value();
            cap_m.internalField()[celli] = buildingMaterial->cap().value();
            //buildingMaterial->update_lambda_cell(ws,lambda_m,celli);
        }
        else if (update_sequence == 3)
        {
            //buildingMaterial->update_w_C_cell(pc,ws,Crel,celli);
						Info << "Not implemtented !!" << endl;
        }
    }


}
