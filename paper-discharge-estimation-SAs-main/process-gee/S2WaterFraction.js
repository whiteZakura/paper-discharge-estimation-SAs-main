//Please load the imports for S2 first.

//-------------------------1. Run one river at a time--------------------------------------//
//Please load the water occurrence frequency of the river into the Imports and name it water_fre_X, where X is the name of the river.

// set folder, water_fre, region 
var folder = 'sittang-cpp0.6-rcpp0.2';
var water_fre = water_fre_sittang;
var region = region_sittang;

var image_cloud_pixel_percentage = 60
var region_cloud_pixel_percentage = 0.2


var Tmin_list = [0.1, 0.2, 0.3, 0.4];
var Tmax_list = [0.6, 0.7, 0.8, 0.9];



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


//-----------------------2. Calculate the total number of pixels within SAs----------------------//

// Select a representative image
var representativeImage = s2_raw
         .filterBounds(region)
         .filterDate('2017-01-01','2017-05-01')
         .median()
         .clip(region)
         .select('B3');

// Calculate the total number of pixels within the region
var totalPixels = representativeImage.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: region,
  scale: 10, 
  maxPixels: 1e9
}).values().get(0);

// Threshold for image coverage within the region
var thresholdTotalPixels = ee.Number(totalPixels).divide(120);



//-----------------------3. Extract water pixels and total pixels within the region----------  //

//Double main loop
Tmin_list.forEach(function(Tmin) {
  Tmax_list.forEach(function(Tmax) {

    

    var SAs = water_fre.gt(Tmin).and(water_fre.lt(Tmax));
    
    
    var filteredImages = s2_raw
      .filterBounds(region) 
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', image_cloud_pixel_percentage)) 
      .filterDate('2015-01-01','2021-01-01')
      //.filter(ee.Filter.calendarRange(6, 8, 'month'))//only for Arctic rivers
      .map(function(image) {
        
        var clipped = image.clip(region);
  

        var pixelCount = clipped.reduceRegion({
          reducer: ee.Reducer.count(),
          geometry: region,
          scale: 10,  
          maxPixels: 1e9
        }).values().get(0);
        
        

        var cloudBitMask = 1 << 10;
        var cirrusBitMask = 1 << 11;
        
        var qa = image.select('QA60');
        var qamask = qa.bitwiseAnd(cloudBitMask).gt(0)
                      .or(qa.bitwiseAnd(cirrusBitMask).gt(0))
                      .updateMask(SAs)
                          
                               
        var cloudPixels  = qamask.reduceRegion({
          reducer: ee.Reducer.sum(),
          geometry: SAs.geometry(),
          scale: 10,
          maxPixels: 1e9
        }).get('QA60');
     
                                  
        var SAsPixels = qamask.reduceRegion({
          reducer: ee.Reducer.count(),
          geometry: SAs.geometry(),
          scale: 10,
          maxPixels: 1e9
        }).get('QA60');
        
        var rcpp = ee.Number(cloudPixels).divide(SAsPixels);
      
        
        return image.set({'rcpp':rcpp,
                          'pixelCount':pixelCount
        }); 
      })
      .filter(ee.Filter.gt('pixelCount', thresholdTotalPixels))
      .filter(ee.Filter.lt('rcpp', region_cloud_pixel_percentage)); 
    
    

    var results = filteredImages.map(function(image) {

      var mndwi = image.normalizedDifference(["B3", "B11"]).rename('mndwi');
      var histogram = mndwi.reduceRegion({
        reducer: ee.Reducer.histogram(), 
        geometry: region, 
        scale: 10,
        maxPixels: 1e13,
        tileScale: 8
      });
      var threshold = otsu(histogram.get("mndwi"));
      var waterMask = mndwi.gte(threshold);
      

      var waterPixel = waterMask.updateMask(SAs).rename("water");
    

      var waterPixelCount = waterPixel.reduceRegion({
        reducer: ee.Reducer.sum(), 
        geometry: SAs.geometry(), 
        scale: 10,
        maxPixels: 1e13
      });
    
      var totalPixelCount = waterPixel.reduceRegion({
        reducer: ee.Reducer.count(), 
        geometry: SAs.geometry(), 
        scale: 10,
        maxPixels: 1e13
      });
    

      return image.set({
        'waterPixelCount': waterPixelCount.get('water'),
        'totalPixelCount': totalPixelCount.get('water'),
        'date': image.date().format('YYYY-MM-dd')
      });
    });
    

    var featureCollection = ee.FeatureCollection(results.map(function(image) {
      return ee.Feature(null, {
        'date': image.get('date'),
        'waterPixelCount': image.get('waterPixelCount'),
        'totalPixelCount': image.get('totalPixelCount')
      });
    }));
    

    var Tmin_int = parseInt(Tmin * 10);
    var Tmax_int = parseInt(Tmax * 10);


    var filename = 'water_' + Tmin_int.toString() + '_' + Tmax_int.toString();
    

    Export.table.toDrive({
      collection: featureCollection,
      description: filename, 
      folder: folder,
      fileFormat: 'CSV'
    });
  });
});













