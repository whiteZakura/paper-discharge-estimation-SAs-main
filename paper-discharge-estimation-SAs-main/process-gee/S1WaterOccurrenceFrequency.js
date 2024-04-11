//Please load the imports for S1 first.
         
//Load JRC maximum extent data for reference
Map.addLayer(JRC.select('max_extent'),{min:0,max:1}, 'JRC max_extent');        


var region = region_sittang
var description = 'water_fre_sittang'
var assetId = 'water_fre_VH_sittang'
    


//OTSU
function otsu(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  var indices = ee.List.sequence(1, size);
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
          bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  return means.sort(bss).get([-1]);
}


var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
          .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
          .filter(ee.Filter.eq('instrumentMode', 'IW'))
          .filterDate('2014-01-01','2022-12-31')
          .filterBounds(region)
          //.filter(ee.Filter.calendarRange(6, 8, 'month'))//only for Arctic Rivers

//Extract water bodies
function calWater(n){

  var imgVH = n.select('VH').rename('VH')
  var edge = imgVH.lt(-30.0);
  var maskedImage = imgVH.mask().and(edge.not())
  imgVH = imgVH.updateMask(maskedImage)
                   .clip(region)
                   .rename('VH')

  var histogram =  imgVH.reduceRegion({
      reducer: ee.Reducer.histogram(), 
      geometry: region, 
      scale: 10,
      maxPixels: 1e13,
      tileScale: 8
    })
    

  var threshold = otsu(histogram.get("VH"))
  var mask =imgVH.lte(threshold)
  var water = mask.rename("water")
  return water.clip(region).copyProperties(n, ["system:time_start"])

}


var water = s1.map(calWater)


//obtain water occurrence frequency
var water_fre = water.reduce(ee.Reducer.mean()).rename("water_fre")


Export.image.toAsset({
  image: water_fre,
  description: description,
  assetId: assetId,
  scale: 10,
  region: region,
  maxPixels: 1e13,
})



