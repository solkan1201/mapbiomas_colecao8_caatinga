
var vis = {
    visMask : {
        min: 0,
        max: 1
    },
    incidente: {
      min: 1,
      max: 5,
      palette :["a6a6a6","d9d9d9","dbed55","ff5050","990033"]
    }
};

var param = {
    'classMapB': [3, 4, 5, 9, 12, 13, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26, 29, 30, 31, 32, 33,
                  36, 39, 40, 41, 46, 47, 48, 49, 50, 62],
    'classNew':  [3, 4, 3, 3, 12, 12, 15, 18, 18, 18, 18, 22, 22, 22, 22, 33, 29, 22, 33, 12, 33,
                  18, 18, 18, 18, 18, 18, 18,  4,  4, 21],
    'assetMapbiomas71': 'projects/mapbiomas-workspace/public/collection7_1/mapbiomas_collection71_integration_v1',
    'asset_fire': 'projects/mapbiomas-workspace/FOGO_COL2/SUBPRODUTOS/mapbiomas-fire-collection2-annual-burned-v1',
    'asset_befFilters': 'projects/mapbiomas-workspace/AMOSTRAS/col7/CAATINGA/classification_Col71_S1v18',
    'asset_filtered': 'projects/mapbiomas-workspace/AMOSTRAS/col7/CAATINGA/class_filtered_Tp',
    'asset_alerts': 'users/data_sets_solkan/Alertas/layersClassTP',
    'asset_output_mask' : 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/masks/maks_layers',
    'assetrecorteCaatCerrMA' : 'projects/mapbiomas-workspace/AMOSTRAS/col7/CAATINGA/recorteCaatCeMA',
    "anoIntInit": 1985,
    "anoIntFin": 2022
};

var imMask_alerts = ee.Image(ee.ImageCollection(param['asset_alerts'])
                        .filter(ee.Filter.eq('yearDep', 2020))
                        .reduce(ee.Reducer.max()).eq(0)).rename('mask_alerta');
print("show image mask alert ", imMask_alerts);

var imMask_fire = ee.Image(param['asset_fire']).select('burned_area_2019')
                          .unmask(0).eq(0).rename('mask_fire');
print("show image mask fire ", imMask_fire);
// var imMaskInc= ee.ImageCollection(param['asset_output_mask'])
var imMaskInc = ee.Image(param['asset_output_mask'] + '/masks_pixels_incidentes_2019');              
print("show imagen mask incidentes ", imMaskInc);

var imMaskpixelEst = ee.Image(param['asset_output_mask'] + '/masks_estatic_pixels_2019');
print("show imagen mask pixeis estaveis ", imMaskInc);

Map.addLayer(imMask_alerts, vis.visMask, 'alerta');
Map.addLayer(imMask_fire, vis.visMask, 'fire');
Map.addLayer(imMaskpixelEst, vis.visMask, 'pixel estav');
Map.addLayer(imMaskInc, vis.incidente, 'Incidentre');