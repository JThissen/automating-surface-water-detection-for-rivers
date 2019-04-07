/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var SRTM = ee.Image("USGS/SRTMGL1_003"),
    NED = ee.Image("USGS/NED"),
    HAND1000 = ee.Image("users/gena/GlobalHAND/30m/hand-1000"),
    fa = ee.Image("users/gena/GlobalHAND/90m-global/fa"),
    HAND5000 = ee.Image("users/gena/GlobalHAND/30m/hand-5000"),
    Rhine = /* color: #d63000 */ee.Geometry.Polygon(
        [[[5.9157120243488635, 51.93205527887385],
          [5.916741992610582, 51.83796589417415],
          [6.178697253841051, 51.81122896914191],
          [6.222299243587145, 51.92760954201599]]]),
    Exclude = /* color: #d63000 */ee.Geometry({
      "type": "GeometryCollection",
      "geometries": [
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                -94.81714085772768,
                30.010325095887058
              ],
              [
                -94.81656150058046,
                30.007705146572096
              ],
              [
                -94.81396512225405,
                30.009061582869844
              ],
              [
                -94.81516675189272,
                30.011421366200835
              ]
            ]
          ],
          "evenOdd": true
        },
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                5.963434233133739,
                51.68904256365271
              ],
              [
                5.94918657333892,
                51.68712703278048
              ],
              [
                5.946783353382102,
                51.66679612784789
              ],
              [
                5.963605891628731,
                51.651781592721505
              ],
              [
                5.992616186922646,
                51.68702061206681
              ],
              [
                5.981458381023458,
                51.69117083628353
              ]
            ]
          ],
          "geodesic": true,
          "evenOdd": true
        },
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                5.995977169599087,
                51.869084747141045
              ],
              [
                5.977437740888149,
                51.852758693287385
              ],
              [
                5.99288726481393,
                51.84766884198093
              ],
              [
                6.005933529462368,
                51.85975629897288
              ],
              [
                6.005933529462368,
                51.867600805240244
              ]
            ]
          ],
          "evenOdd": true
        },
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                5.982579420079446,
                51.92046759533146
              ],
              [
                5.984210203160501,
                51.91866771885133
              ],
              [
                5.985154340733743,
                51.91707953261963
              ],
              [
                5.99073333548472,
                51.91538540538169
              ],
              [
                5.996913145055032,
                51.916920710906524
              ],
              [
                5.991334150304056,
                51.92258500468874
              ],
              [
                5.982751081456399,
                51.92480817702574
              ],
              [
                5.981635282506204,
                51.92274380636441
              ]
            ]
          ],
          "evenOdd": true
        },
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                5.858497809206256,
                51.7607736816366
              ],
              [
                5.862338904044577,
                51.760175999406634
              ],
              [
                5.864699353379251,
                51.759511898760834
              ],
              [
                5.866330209281159,
                51.75911343368237
              ],
              [
                5.865772284928198,
                51.7617166752547
              ],
              [
                5.8527683548064715,
                51.772284256092256
              ],
              [
                5.845386585830738,
                51.7711423032311
              ],
              [
                5.846888690053788,
                51.76691948256114
              ],
              [
                5.8513091679857325,
                51.76431654098316
              ]
            ]
          ],
          "geodesic": true,
          "evenOdd": true
        },
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                5.834667688266677,
                51.864031722234216
              ],
              [
                5.835697656528396,
                51.86336917979553
              ],
              [
                5.837156778232497,
                51.8624946088347
              ],
              [
                5.858185296909255,
                51.85332384108857
              ],
              [
                5.858957773105544,
                51.86450874674787
              ],
              [
                5.82960367764656,
                51.87611478412917
              ],
              [
                5.8229947146338645,
                51.870391633161596
              ]
            ]
          ],
          "evenOdd": true
        },
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                5.828959947482986,
                51.863081726242974
              ],
              [
                5.825140481845779,
                51.86506933698577
              ],
              [
                5.823810106174392,
                51.864857329357484
              ],
              [
                5.819346910373611,
                51.85405220754775
              ],
              [
                5.833079820529861,
                51.85362808744885
              ],
              [
                5.834367280857009,
                51.85877027405492
              ],
              [
                5.827929979221267,
                51.861897708387005
              ]
            ]
          ],
          "evenOdd": true
        },
        {
          "type": "Polygon",
          "coordinates": [
            [
              [
                6.029231297652132,
                51.94403257571618
              ],
              [
                6.042280628756771,
                51.91650517309937
              ],
              [
                6.088983452067737,
                51.900192522047966
              ],
              [
                6.105724349273714,
                51.88670892753365
              ],
              [
                6.10791354336493,
                51.885543161563746
              ],
              [
                6.111873409104646,
                51.88371496762613
              ],
              [
                6.1140411409892295,
                51.88489402833141
              ],
              [
                6.096394280772984,
                51.916628154545535
              ],
              [
                6.056815256387608,
                51.944394183961855
              ]
            ]
          ],
          "geodesic": true,
          "evenOdd": true
        }
      ],
      "coordinates": []
    }),
    Trinity = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-94.82052328104442, 29.9736458377149],
          [-94.81992246622508, 29.92768611853083],
          [-94.77511884684031, 29.92738857543992],
          [-94.77511884684031, 29.973274080577227],
          [-94.82052328104442, 29.9736458377149]]]),
    Budapest = /* color: #d63000 */ee.Geometry.Polygon(
        [[[18.993118161230427, 47.65377145030989],
          [18.998611325292927, 47.37365543516635],
          [19.155166501074177, 47.38574421033091],
          [19.155166501074177, 47.651921294998765]]]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//Thesis script Joost Thissen February, 2019.

//unsolved issue (centerline) [bug in GEE]: https://groups.google.com/forum/#!searchin/google-earth-engine-developers/geometry$20bug%7Csort:date/google-earth-engine-developers/z65gS-IfaqQ/sQmpLXNDEQAJ
//MAIN VARIABLES
//*--------------------------------*//
//Satellite of choice: Landsat 8 [0] or Sentinel-2[1]
var preferredSatellite = 1;

//Enter a date of interest [startDate]. In case an image is not available, the next most recent following image is selected (yyyy-mm-dd).
var dateOfInterest = {startDate: "2016-04-01", endDate: "2018-10-01"}; //2017-01-26

//Bounds of the image. Create a geometry under "Geometry Imports" and draw a shape of the area of interest.
var imageBounds = Budapest;//Trinity;//Rhine;

//Water occurrence composite days relative to the date of interest.
var samplingRange = {daysPositive: 120, daysNegative: 120};

//Projection scale.
var scale = 10;

//Opening radius.
var smoothingAmount = 20;

//Water occurrence composite threshold value.
var cutOffThreshold = {min: 0.2, max: 1};

//Prune undesirable branches modifier.
var pruneRadiusMultiplier = 3;

//NDWI or MNDWI.
var waterBands = ['green', 'nir']; //ndwi
//waterBands = ["green", "swir1"] //mndwi

//Resampling method ("bilinear" or "bicubic");
var resample = "bilinear";

//Simplify polygon by a certain factor.
var simplificationFactor = 2;

//The pruning strength. Define how many iterations you'd like to run.
var pruneIterations = 3;

//Sample river width every [sampleStepSize] meters.
var sampleStepSize = 10;

//Weigh every image of the historical composite based on the time between the point of reference and the sampled image. If false, the water occurrence composite represents an average.
var useInverseDistanceWeight = true;

//weight distribution equals 1/abs(d)^x, where x = weightExponent. 0 equals no weight.
var weightExponent = 0.1;

//kernel radius r (erosion)
var focalMin = {lower: ee.Number(0), upper: ee.Number(25)};

//K - proportionality between kernel radius r and Aratio (posterior/i)
var scaleFactor = {lower: ee.Number(1), upper: ee.Number(1.5)};

//Filtering based on a river width threshold is resource intensive, as PrunedCenterline is called for every detected polygon. 
//For performance reasons, filter by a polygon that contains the most amount of vertices within an image [false].
var useRiverWidthThreshold = false;

//River width thershold range.
var riverWidthThreshold = {minimum: 0, maximum: 300};

//Zoom level of the map.
var zoomScale = 13;

//Maximum number of vertices (Eucledean distance map).
var maxVerticesPerPolygon = 500;

//Simplify centerline by a certain factor.
var simplificationFactorCenterline = 3;
//*--------------------------------*//

var waterOccurrenceComposite = require("users/JThissen/ThesisScripts:WaterOccurrenceComposite");
var debugScript = false;
var debug = false;
var handOnly = false;
var handMask = true;
var usePrior = true;
var fillGaps = true;
var prior = null;
var projection = ee.Projection("EPSG:3857").atScale(scale);
var error = ee.ErrorMargin(scale, "meters");
var dem = SRTM;
var bounds = imageBounds.bounds();
var availableSatellites = ["Landsat8", "Sentinel2"];
var handTh = 15;
var hand = HAND1000;
Map.centerObject(imageBounds, zoomScale); 

if(preferredSatellite == 0)
  var images = Landsat8(bounds, dateOfInterest);
else if(preferredSatellite == 1)
  var images = Sentinel2(bounds, dateOfInterest);

GenerateRiverPolygon(images, 'system:start_time', true, waterBands, handOnly, handMask, usePrior, fillGaps);

//FUNCTIONS
//*--------------------------------*//
// rescales to given ranges
function rescale(img, exp, thresholds) {
    return img.expression(exp, { img: img }).subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
}
// used as aside function for debugging
function show(image, name, vis) {
    if (debug) {
        Map.addLayer(image, vis || {}, '  ' + name, false);
    }
    return image;
}

//deg2Rad image
function radians(img) { return img.toFloat().multiply(3.1415927).divide(180); }

//deg2Rad
function radians(deg) {
  return deg.multiply(3.1415927).divide(180);
}

//canny edge detector
function getEdge(mask) 
{
    var canny = ee.Algorithms.CannyEdgeDetector(mask, 0.99, 0);
    return canny.mask(canny);
}

//return the DN that maximizes interclass variance in B5 (in the region).
function otsu(histogram) {
    histogram = ee.Dictionary(histogram);

    var counts = ee.Array(histogram.get('histogram'));
    var means = ee.Array(histogram.get('bucketMeans'));
    var size = means.length().get([0]);
    var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
    var mean = sum.divide(total);

    var indices = ee.List.sequence(1, size);

    // Compute between sum of squares, where each mean partitions the data.
    var bss = indices.map(function (i) {
        var aCounts = counts.slice(0, 0, i);
        var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
        var aMeans = means.slice(0, 0, i);
        var aMean = aMeans.multiply(aCounts).reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount);
        var bCount = total.subtract(aCount);
        var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
        return aCount.multiply(aMean.subtract(mean).pow(2)).add(bCount.multiply(bMean.subtract(mean).pow(2)));
    });

    // Return the mean value corresponding to the maximum BSS.
    return means.sort(bss).get([-1]);
};

//Otsu thresholding
function computeThresholdUsingOtsu(image, scale, bounds, th, g, skipShort, weightGradient, minValue, usePrior) 
{
    // clip image edges
    var mask = image.mask().gt(0).focal_min(ee.Number(scale).multiply(3), 'circle', 'meters');

    // detect sharp changes
    var edge = ee.Algorithms.CannyEdgeDetector(image, th, g);
    edge = edge.multiply(mask);
    //print(edge, "edge")

    if(usePrior) {
      edge = edge.multiply(prior).updateMask(prior);
       //print(edge, "edge updated")
    }

    // take the largest changes, estimate gradient around edge and use that as a weight
    if (weightGradient) {
        var gradient = image.gradient().abs();
        var edgeGradient = gradient.select(0).max(gradient.select(1)).mask(edge.gt(th)).reproject(image.projection().scale(2, 2));
        // take the upper percentiles only
        var p85 = ee.Number(ee.Dictionary(edgeGradient.reduceRegion(
          {
            reducer: ee.Reducer.percentile([85]),
            geometry: bounds,
            scale: scale,
            maxPixels: 1e9
            })).values().get(0));
        var significantEdgesMask = edgeGradient.gt(p85);
        
        // take the mode
        /*
        var mode = ee.Number(ee.Dictionary(edgeGradient.reduceRegion(ee.Reducer.mode(), bounds, scale)).values().get(0));
        var σ = ee.Number(ee.Dictionary(edgeGradient.reduceRegion(ee.Reducer.stdDev(), bounds, scale)).values().get(0));
        var _buckets = 50;
        var significantEdgesMask = edgeGradient.gt(mode);
        */

        edge = edge.updateMask(significantEdgesMask);
        //print(edge)

        if (debug) {
            // gradient around edges
            if (edgeGradient) {
                print(ui.Chart.image.histogram(edgeGradient, bounds, scale, _buckets));
                Map.addLayer(edgeGradient, {}, 'edge gradient', false);
                Map.addLayer(significantEdgesMask.mask(significantEdgesMask), {}, 'significant edges', false);

                print('Mode: ', mode);
                print('Sigma: ', σ);
                //Map.addLayer(edgeGradient.updateMask(significantEdgesMask), {min:0, max:mode.add(σ.multiply(2)), palette:['ffffff', 'ff0000']}, 'edge gradient, upper percentiles', false)
            }
        }
    }

    // advanced, detect edge lengths
    var coonnectedVis = void 0;
    if (skipShort) {
        var connected = edge.mask(edge).lt(0.8).connectedPixelCount(50, true);

        var edgeLong = connected.gte(50);

        edge = edgeLong;

        coonnectedVis = connected.updateMask(edgeLong).visualize({ palette: ['ffffff', 'ff0000'], min: 0, max: 50 });
    }

    // buffer around NDWI edges
    var edgeBuffer = edge.focal_max(ee.Number(scale).multiply(1), 'square', 'meters');
    var imageEdge = image.mask(edgeBuffer);

    // compute threshold using Otsu thresholding
    var buckets = 100;
    var hist = ee.Dictionary(ee.Dictionary(imageEdge.reduceRegion(
      {
        reducer: ee.Reducer.histogram(buckets),
        geometry: bounds,
        scale: scale,
        maxPixels: 1e9
      })).values().get(0));

    var threshold = ee.Algorithms.If(hist.contains('bucketMeans'), otsu(hist), 0);
    threshold = ee.Number(threshold); //.add(0.05)

    if (debug) {
        Map.addLayer(edge.mask(edge), { palette: ['ff0000'] }, 'edges', false);

        if (skipShort) {
            Map.addLayer(coonnectedVis, {}, 'edges (connected)', false);
        }

        print('Threshold: ', threshold);

        print(ui.Chart.image.histogram(image, bounds, scale, buckets));
        print(ui.Chart.image.histogram(imageEdge, bounds, scale, buckets));
        Map.addLayer(mask.mask(mask), { palette: ['000000'] }, 'image mask', false);
    }

    return {
      threshold: minValue ? threshold.max(minValue) : threshold,
      edge: edge
    };
}

//Sample water occurrence composite and merge
function fillWater(i) 
{
  var i1 = i.focal_min(30, 'square', 'meters').focal_max(15, 'square', 'meters');
  var edge = getEdge(i);
  var occurrenceExpected = prior.mask(edge)
    .reduceRegion(
      {
        reducer: ee.Reducer.intervalMean(15, 50),
        geometry: bounds,
        scale: 10,
        maxPixels: 1e9
      }).values().get(0);
      
  occurrenceExpected = ee.Algorithms.If(ee.Algorithms.IsEqual(occurrenceExpected, null), 1, occurrenceExpected);
  var guess = prior.gt(ee.Image.constant(occurrenceExpected)).or(i);
  
  var guessArea = guess.reduceRegion(
    {
      reducer: ee.Reducer.sum(),
      geometry: bounds,
      scale: scale,
      maxPixels: 1e9
    }).values().get(0);
    
  var totalArea = i.reduceRegion(
    {
      reducer: ee.Reducer.sum(),
      geometry: bounds,
      scale: scale,
      maxPixels: 1e9
    }).values().get(0);
    
  var updatedAreaFraction = ee.Number(guessArea).divide(totalArea);
  var slope = (focalMin.lower.subtract(focalMin.upper)).divide((scaleFactor.upper.subtract(scaleFactor.lower)));
  var y = slope.multiply(updatedAreaFraction.subtract(scaleFactor.lower)).add(focalMin.upper);
  y = y.min(focalMin.upper);
  y = y.max(focalMin.lower);
  var posterior = ee.Algorithms.If(updatedAreaFraction.lt(scaleFactor.lower),
    prior.gt(ee.Image.constant(occurrenceExpected)).focal_min(focalMin.upper, 'circle', 'meters').or(i),
      ee.Algorithms.If(updatedAreaFraction.gt(scaleFactor.upper),
        prior.gt(ee.Image.constant(occurrenceExpected)).focal_min(focalMin.lower, 'circle', 'meters').or(i),
           prior.gt(ee.Image.constant(occurrenceExpected)).focal_min(y, 'circle', 'meters').or(i)));
  
  print("Kernel radius [pixels]: ", y, "", "A-ratio [-]: ", updatedAreaFraction, "");
  posterior = ee.Image(posterior)
    .copyProperties(i)
    .set('occurrence', occurrenceExpected)  
    .set('system:time_start', i.get('system:time_start'));
    
  return {area: ee.Algorithms.If(ee.Algorithms.IsEqual(occurrenceExpected, null), null, posterior), fraction: updatedAreaFraction};
}

//Acquire Landsat 8 data
function Landsat8(bounds, dateOfInterest)
{
  var bands = ['swir1', 'nir', 'green', 'red', 'blue', 'cloud'];
  var l8 = new ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').map(function(i){return ee.Algorithms.Landsat.simpleCloudScore(i);})
                                                            .select(['B6', 'B5', 'B3', 'B4', 'B2', 'cloud'], bands);
  var images = l8.filterBounds(bounds.centroid(100))
                  .filterDate(dateOfInterest.startDate, dateOfInterest.endDate)
                  //.filterMetadata("CLOUD_COVER", "less_than", 15)
                  .sort("system:time_start");
                  
  print("Landsat 8 images available:", images.toList(5000, 0).size());
  return images;
}

//Acquire Sentinel 2 data
function Sentinel2(bounds, dateOfInterest)
{
  var s2BandsNative = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12', 'QA10', 'QA20', 'QA60'];
  var s2BandsReadable = ['coastal', 'blue', 'green', 'red', 'red2', 'red3', 'red4', 'nir', 'nir2', 'water_vapour', 'cirrus', 'swir1', 'swir2', 'QA10', 'QA20', 'QA60'];
  var bands = ['swir1', 'nir', 'green', 'red', 'blue', 'cloud'];
  var bandsSentinel2 = ['nir', 'green', 'red', 'blue'];
  var filter = ee.Filter.and(ee.Filter.bounds(bounds.centroid(100)), ee.Filter.date(dateOfInterest.startDate, dateOfInterest.endDate));
  var images = ee.ImageCollection('COPERNICUS/S2').select(s2BandsNative, s2BandsReadable).filter(filter);
  images = (images.select(bandsSentinel2, bandsSentinel2)).sort('system:time_start')
  //.filterMetadata("CLOUDY_PIXEL_PERCENTAGE", "less_than", 15).sort("system:time_start");
  print("Sentinel-2 images available:", images.toList(5000, 0).size());
  return images;
}

//Generate the river polygon
function GenerateRiverPolygon(images, sortProperty, waterOtsu, waterBands, handOnly, handMask, usePrior, fillGaps) 
{
  var boundsImage = ee.Image().toInt().paint(bounds, 1).reproject(ee.Image(images.first()).projection()) ;
  var listOfImages = images.toList(images.size());
  var dayToMs = ee.Number(86400000);
  var value = ee.Number(0);
  var i = ee.Image(ee.List(listOfImages).get(value)).resample(resample);
  var systemTime = i.get("system:time_start");
  var minimumDate = ee.Date(ee.Number(systemTime).subtract(dayToMs.multiply(ee.Number(samplingRange.daysNegative))));
  var maximumDate = ee.Date(ee.Number(systemTime).add(dayToMs.multiply(ee.Number(samplingRange.daysPositive))));
  prior = waterOccurrenceComposite.app(minimumDate, maximumDate, imageBounds, cutOffThreshold.min, cutOffThreshold.max, systemTime, weightExponent, useInverseDistanceWeight);
  prior = prior.unmask();
  Map.addLayer(prior.clip(bounds), {}, "Unmasked historical composite (bounded)", false)/
  print("Date satellite image: ", ee.Date(systemTime), "", "Historical imagery min. threshold", minimumDate, "", 
        "Historical imagery max. threshold", maximumDate);
  
  if(preferredSatellite == 0)
    Map.addLayer(i.clip(bounds), {}, "SelectedImage", true);
  else if(preferredSatellite == 1)
    Map.addLayer(i.clip(bounds), { bands: ["red", "green", "blue"], min: 250, max: 3000}, "SelectedImage", true);
  
  Map.addLayer(ee.Image(0).clip(imageBounds), {}, "BlackBackground", false);
  i = i.updateMask(boundsImage.multiply(i.mask()));
  if(typeof waterBands === 'undefined') {waterBands = ['green', 'nir']}
  var waterScore = i.normalizedDifference(waterBands);
  var dateString = ee.Date(i.get('system:time_start')).format("dd-MM-YYYY");
  
  if(handMask) 
    waterScore = waterScore.mask(waterScore.mask().multiply(hand.lt(handTh)));
  else 
    waterScore = waterScore.mask(waterScore.mask());
    
  var th = 0;
  var edge = ee.Image();
  
  if(waterOtsu) 
  {
    var o = computeThresholdUsingOtsu(waterScore, 30, bounds, 0.4, 0.5, false, false, -0.1, usePrior);
    th = o.threshold;
    edge = o.edge.reproject(i.projection().scale(0.5, 0.5));
  }
  
  var waterMask = waterScore.gte(th);
  
  if(usePrior) 
  {
    //Map.addLayer(waterMask, {}, "watermask before multiplication", false);
    waterMask = waterMask.multiply(prior).gt(0);
    
    if(debugScript)
      Map.addLayer(waterMask, {}, "watermask after multiplication", false);
    
    if(fillGaps) 
    {
      var fill = fillWater(waterMask);
      waterMask = fill.area;
    }
  }
  
  var waterEdge = getEdge(waterMask).focal_max(1);
  waterMask = ee.Image(waterMask);

  if(handOnly) 
  {
    return hand.mask(hand.gt(handTh)).visualize({min:handTh, max:1000})
      .updateMask(boundsImage);
  }
    
  waterMask = ee.Image(waterMask);
  waterMask = waterMask.focal_min({radius: ee.Number(smoothingAmount), units: "meters"}).focal_max({radius: ee.Number(smoothingAmount), units: "meters"});
  Map.addLayer(waterMask, {}, "Surface water mask (SWM)", false)
  var waterOnly = waterMask.updateMask(waterMask);
  var sampleImage = ee.Image(images.first());
  var projection = sampleImage.select(0).projection();
  var waterPolygon = waterOnly.reduceToVectors(
  {
    reducer: ee.Reducer.countEvery(),
    geometry: imageBounds,
    scale: scale,
    crs: projection,
    geometryType: "polygon", 
    eightConnected: false,
    bestEffort: true,
    tileScale: 4
  });
  
  var simplifiedPolygon = waterPolygon.map(function(item)
  {
    item = item.simplify(scale * simplificationFactor);
    item = ee.Feature((item.geometry(error)).difference(Exclude, error));
    var type = item.geometry().type();
    return item.set({size: item.geometry().coordinates().flatten().length().divide(2), geometryType: type});
  }, true).filter(ee.Filter.eq("geometryType", "Polygon"));
  
  if(useRiverWidthThreshold)
  {
    //Filter based on river width. Resource-intensive since PrunedCenterline is called for every feature.
    var distancePerPolygon = simplifiedPolygon.map(function(feature)
    {
      var geom = feature.geometry(error);
      var prunedCenterline = PrunedCenterline(geom, pruneRadiusMultiplier, maxVerticesPerPolygon);
      return feature.set({meanWidth: prunedCenterline.width, centerline: prunedCenterline.centerline});
    }, true).filter(ee.Filter.greaterThan("meanWidth", ee.Number(riverWidthThreshold.minimum)))
            .filter(ee.Filter.lessThan("meanWidth", ee.Number(riverWidthThreshold.maximum)));
    
    var dppList = distancePerPolygon.toList(distancePerPolygon.size());
    var dppSize = distancePerPolygon.size().getInfo();
    
    for(var j = 0; j < dppSize; j++)
      Map.addLayer(ee.FeatureCollection(ee.Feature(dppList.get(j)).get("centerline")), {color: "f71eff"}, "Centerline", false);
  }
  else
  {
   
    //Filter based on maximum vertices within a polygon. Fast since PrunedCenterline is called only once.
    var largestValue = simplifiedPolygon.aggregate_max("size");
    var largestPolygon = simplifiedPolygon.filter(ee.Filter.eq("size", largestValue));
    var filteredGeo = (ee.Feature(largestPolygon.first()).geometry());//.difference(Exclude, error);
    var result = PrunedCenterline(filteredGeo, pruneRadiusMultiplier, maxVerticesPerPolygon);
    Map.addLayer(filteredGeo, {color: "00ffff"}, "RiverPolygon", true);
    Map.addLayer(result.centerline, {color: "f71eff"}, "FinalCenterline", true);
    print("Mean width [m]", result.width);
    print("Geometry area [m2]", filteredGeo.area());
  }
  
  var paintPoints = ee.Geometry.MultiPoint(filteredGeo.coordinates().flatten());
  Map.addLayer(ee.Image().paint(paintPoints.buffer(9)), {palette:['ffffff']}, 'PolygonPoints', false)
  
  if(debugScript)
  {
    Map.addLayer(waterScore, {}, "normalizedDifference", false);
    print(ee.Date(i.get('system:time_start')));
    print(currentImage, "selected image");
    Map.addLayer(waterScore, {}, dateString.getInfo(), false);
    Map.addLayer(hand, {}, "hand", false);
    Map.addLayer(hand.lt(handTh), {}, "hand lower than handTh", false);
    Map.addLayer(waterScore.mask(), {}, "waterScore mask only", false);
    Map.addLayer(waterScore, {}, "waterScore w/ handMask", false);
    Map.addLayer(waterMask, {}, "waterMask after Otsu", false);
    Map.addLayer(prior, {}, "priorLayer", false);
    Map.addLayer(waterEdge, {}, "waterEdge focal_max(1)", false);
    Map.addLayer(waterPolygon, {color: "00ffff"}, "waterPolygon", false);
    Map.addLayer(simplifiedPolygon, {color: "00ffff"}, "simplifiedPolygon", false);
    Map.addLayer(points, {color: "ffffff"}, "points", false);
    Map.addLayer(ee.Image().paint(verticesPolygon.buffer(5)), {palette:["ffffff"]}, 'verticesPolygon', false);
  }
 
  return null;
}

//Generate the pruned centerline
function PrunedCenterline(geometryInput, pruneRadiusMultiplier, verticesCap)
{
  var error = ee.ErrorMargin(scale, "meters");
  var geom = geometryInput;
  var amount =  ee.List(geom.coordinates().get(0)).size().multiply(ee.Number(scale)).min(verticesCap);
  var simplifyFactor = simplificationFactorCenterline;
  var proj = ee.Projection("EPSG:3857").atScale(scale); 
  var geometryBuffer = geom.buffer(scale * 2.5, error);
  var perimeterGeometry = geometryBuffer.difference(geometryBuffer.buffer(-scale * 3.5, error), error);
  var points = GeneratePoints(geometryBuffer, amount, error);
  var generatedPolygons = GeneratePolygons(geometryBuffer, points, scale, error, proj);
  var polygons = generatedPolygons.resultingPolygons;
  //var geometryBuffer2 = geom.buffer(scale, error);
  var points2 = GeneratePoints(geom, amount, error);
  print("Amount of points generated (EDM): ", points2.size());
  Map.addLayer(points2, {color: "ffffff"}, "EDM Points", false);
  var distance = ee.Image(0).byte().paint(points2, 1).fastDistanceTransform()
                           .sqrt().multiply(ee.Image.pixelArea().sqrt()).multiply(2).clip(geom).reproject(proj);
                           
  Map.addLayer(distance, {min: 0, max: 400}, "Euclidean distance map (EDM)", false);
  var distFilter = ee.Filter.and(ee.Filter.intersects({leftField: '.geo', rightField: '.geo', maxError: error}),
                                  ee.Filter.notEquals({leftField: 'labels', rightField: 'labels',}));
  var distSaveAll = ee.Join.saveAll({matchesKey: 'matches',});
  var featureCollection = distSaveAll.apply(polygons, polygons, distFilter);
  var f = featureCollection.map(function(feature) 
  { 
    var matches = ee.FeatureCollection(ee.List(feature.get('matches'))) ;
    return matches.map(function(matchFeature) 
    {
      var line = matchFeature.intersection(feature, error, proj);
      return line
            .set({ intersectsPerimeter: line.intersects(perimeterGeometry, error, proj) })
            .set({ intersectsPolygon: line.intersects(geometryBuffer, error, proj) });
    });
  }).flatten();
  
  var centerline = f.filter(ee.Filter.eq('intersectsPerimeter', false))
                    .filter(ee.Filter.eq('intersectsPolygon', true));
  
  //Prune
  var setGeometryType = ee.FeatureCollection(centerline.geometry().geometries().map(function(item)
  {
    item = ee.Geometry(item);
    var result = item.type();
    return ee.Feature(item).set({geometryType: result});
  }));

  centerline = setGeometryType.filter(ee.Filter.neq("geometryType", "Polygon")).filter(ee.Filter.neq("geometryType", "Point"));
  centerline = centerline.geometry().dissolve(error, proj).simplify(scale * simplifyFactor, proj);
  var originalCenterline = centerline;
  Map.addLayer(originalCenterline, {color: "f71eff"}, "InitialCenterline", false);
  
  var multipoints = ee.FeatureCollection(centerline.coordinates().map(function(i)
  {
    return ee.Feature(ee.Geometry.MultiPoint(i, proj));
  }));
  
  var requiredPoints = multipoints.map(function(feature)
  {
    var first = feature.geometry().coordinates().get(0);
    var last = feature.geometry().coordinates().get(-1);
    var combined = ee.Geometry.MultiPoint([first, last], proj);
    return ee.Feature(combined);
  }).geometry();
  
  var allRequiredPoints = ee.FeatureCollection(requiredPoints.coordinates().map(function(item)
  {
    return ee.Feature(ee.Geometry.Point(item, proj));
  }));
  
  var centerlineFC = ee.FeatureCollection(centerline.geometries().map(function(i)
  {
    return ee.Feature(ee.Geometry(i));
  })).map(function(j)
  {
    var id = j.id();
    return j.set({featureID: id});
  });
  
  for(var i = 0; i < pruneIterations; i++)
  {
    var pointTypes = allRequiredPoints.map(function(point)
    {
      var checkIntersections = centerlineFC.map(function(feature)
      {
        var doIntersect = ee.Algorithms.If(point.intersects(feature, ee.ErrorMargin(scale, "meters"), proj), ee.Number(1), ee.Number(0));
        return feature.set({intersects: doIntersect});
      });
      var intersectCount = checkIntersections.filter(ee.Filter.gte("intersects", 1)).size();
      var pointDefinition = ee.Algorithms.If(intersectCount.lte(ee.Number(1)), ee.String("endPoint"), ee.String("branchPoint"));
      return ee.Feature(point).set({pointType: pointDefinition, intersections: intersectCount});
    });
    
    var endPoints = pointTypes.filter(ee.Filter.eq("pointType", "endPoint"));
    var branchPoints = pointTypes.filter(ee.Filter.eq("pointType", "branchPoint"));
    
    branchPoints = branchPoints.map(function(feature)
    {
      var distanceAmount = distance.reduceRegion(ee.Reducer.first(), feature.geometry(), scale);
      feature = ee.Feature(feature.geometry(error).buffer(ee.Number(distanceAmount.get("distance")).multiply(0.5).multiply(ee.Number(pruneRadiusMultiplier))));
      return feature;
    });
    
    var intersectFilter = ee.Filter.intersects({leftField: ".geo", rightField: ".geo"});
    var intersectJoin = ee.Join.saveAll({matchesKey: "intersects", measureKey: "bar"});
    var intersectJoined = intersectJoin.apply(branchPoints, endPoints, intersectFilter);//.aside(print, "intersectJoined");
    var toRemove = ee.FeatureCollection(ee.List(intersectJoined.aggregate_array("intersects")).flatten());
    
    var filter = ee.Filter.intersects({leftField: ".geo", rightField: ".geo"});
    var simpleJoin = ee.Join.saveAll({matchesKey: "intersects", measureKey: "bar"});
    var simpleJoined = simpleJoin.apply(centerlineFC, toRemove, filter);//.aside(print, "simpleJoined");
    var listToRemove = simpleJoined.aggregate_array("featureID");//.aside(print);
    centerlineFC  = centerlineFC.filter(ee.Filter.inList("featureID", listToRemove).not());

    Map.addLayer(centerlineFC, {color: "f71eff"}, "PrunedCenterline iteration: " + i.toString(), false);
  }

  var paintCenterline = ee.Image(1).int().paint(centerlineFC, 0, 5).not();
  var test = centerlineFC.map(function(feature)
  {
    var pointsSampled = sampleLinePoints(feature.geometry(ee.ErrorMargin(scale, "meters")), sampleStepSize);
    
    return pointsSampled;
  }).flatten();
  
  var k = ee.Number(test.size()).subtract(2);
  var centerlineList = test.toList(test.size());
  var setOffset = ee.List.sequence(0, k).iterate(function(n, list) 
   {
    var index = ee.Number(n);
    n = ee.Number(n).multiply(sampleStepSize);
    var previous = ee.Feature(centerlineList.get(index)).set("offset", n);
    return ee.List(list).add(previous);
   }, ee.List([]));
  
  var collection = ee.FeatureCollection(ee.List(setOffset));

  var riverWidthCollection = collection.map(function(feature)
  {
    var sample = distance.reduceRegion(ee.Reducer.mean(), feature.geometry(error), scale);
    return feature.set({riverWidth: sample.values().get(0)});
  }).aside(print, riverWidthCollection);
  
  var meanRiverWidth = ee.Number(riverWidthCollection.aggregate_mean("riverWidth")).int();
  print(ui.Chart.feature.byFeature(riverWidthCollection, "offset", "riverWidth").setOptions({ lineWidth: 1, pointSize: 0 }));
  
  //Export river width and distance data for use.
  Export.table.toDrive(
    {
      description: "CenterlineData",
      collection: riverWidthCollection,
      fileFormat: "CSV"
    });
    
  //Color river width based on XML properties (rainbow);
  // var distanceColors = '\
  // <RasterSymbolizer>\
  //   <ColorMap extended="false" >\
  //     <ColorMapEntry color="#000099" quantity="10" label="lte 10" />\
  //     <ColorMapEntry color="#0000CC" quantity="20" label="10-20" />\
  //     <ColorMapEntry color="#0000FF" quantity="30" label="20-30" />\
  //     <ColorMapEntry color="#0033FF" quantity="40" label="30-40" />\
  //     <ColorMapEntry color="#0066FF" quantity="50" label="40-50" />\
  //     <ColorMapEntry color="#0099FF" quantity="60" label="50-60" />\
  //     <ColorMapEntry color="#00CCFF" quantity="70" label="60-70" />\
  //     <ColorMapEntry color="#00FFFF" quantity="80" label="70-80" />\
  //     <ColorMapEntry color="#33FFCC" quantity="90" label="80-90" />\
  //     <ColorMapEntry color="#66FF99" quantity="100" label="90-100" />\
  //     <ColorMapEntry color="#99FF66" quantity="110" label="100-110" />\
  //     <ColorMapEntry color="#CCFF33" quantity="120" label="110-120" />\
  //     <ColorMapEntry color="#FFFF00" quantity="130" label="120-130" />\
  //     <ColorMapEntry color="#FFCC00" quantity="140" label="130-140" />\
  //     <ColorMapEntry color="#FF9900" quantity="150" label="140-150" />\
  //     <ColorMapEntry color="#FF6600" quantity="160" label="150-160" />\
  //     <ColorMapEntry color="#FF3300" quantity="170" label="160-170" />\
  //     <ColorMapEntry color="#FF0000" quantity="180" label="170-180" />\
  //     <ColorMapEntry color="#CC0000" quantity="190" label="180-190" />\
  //     <ColorMapEntry color="#990000" quantity="200" label="190-200" />\
  //   </ColorMap>\
  // </RasterSymbolizer>';

  // var centerlineColored = centerlineMasked.sldStyle(distanceColors);
  // //Map.addLayer(centerlineColored, {}, "centerlineColored", false);
  // //Map.addLayer(centerlineColored.focal_max(1),{}, "centerlineColoredDilated", false);
  
  return {
    width: meanRiverWidth, 
    centerline: centerlineFC
  };
}

function Clamp(num, min, max)
{
  return num <= min ? min : num >= max ? max : num;
}

//EDM points
function GeneratePoints(geom, stepSize, errorMargin)
{
  var peri = geom.perimeter(errorMargin);
  var lineString = ee.Algorithms.GeometryConstructors.MultiLineString(geom.coordinates());
  var distances = ee.List.sequence(0, peri, peri.divide(stepSize));
  var feature = ee.Feature(lineString).set({distanceList: distances}).set({distanceLength: distances.length()});
  var distances = feature.get("distanceList");
  var segments = feature.geometry().cutLines(distances, errorMargin).geometries();
  var result = segments.map(function(i)
  {
    return ee.Feature(ee.Geometry(i).centroid(errorMargin));
  });
  return ee.FeatureCollection(result);
}

//Compute Voronoi-based polygons.
function GeneratePolygons(geom, vertices, scale, error, projection)
{
  var radiusMultiplier = 3;
  var distanceTransform = ee.Image(0).byte().paint(vertices, 1).fastDistanceTransform().sqrt().multiply(ee.Image.pixelArea().sqrt()).multiply(2).clip(geom).reproject(projection);
  var convolution = distanceTransform.convolve(ee.Kernel.laplacian8());
  
  var polygons = ee.Image(1).mask(convolution.gt(0))
  .reduceToVectors({
    scale: scale,
    geometry: geom
  });
  
  convolution = convolution.multiply(distanceTransform);
  var edges = convolution.lt(0);
  
  var connected = edges.not().connectedComponents(ee.Kernel.circle(1), 256).clip(geom)
    .focal_max(scale * radiusMultiplier, 'circle', 'meters')
    .focal_min(scale * radiusMultiplier, 'circle', 'meters') 
    .focal_mode(scale * radiusMultiplier, 'circle', 'meters')
    .reproject(projection.atScale(scale));
   
  var polygons = connected.select('labels').reduceToVectors(
    {
      scale: scale,
      crs: projection,
      geometry: geom,
      eightConnected: true,
      labelProperty: 'labels',
      tileScale: 4
    });
  
  //Map.addLayer(polygons.style({color: 'black', fillColor: "lightblue", width: 1}), {}, 'polygons', false, 1.0);
  Map.addLayer(polygons.map(function(f) { return f.simplify(scale * 2)}).style({color: 'black', fillColor: "lightblue"}), {}, 'Voronoi polygons (smoothed)', false);
  return {resultingPolygons: polygons, distanceTransform: distanceTransform};
}

function sampleLinePoints(lineString, step) 
{
  var length = lineString.length(ee.ErrorMargin(scale, "meters"));
  step = ee.Algorithms.If(ee.Number(length).lte(ee.Number(step)), ee.Number(length), ee.Number(step));
  var distances = ee.List.sequence(0, length, step);

  function makePointFeature(coord, offset) {
    var pt = ee.Algorithms.GeometryConstructors.Point(coord).buffer(5, ee.ErrorMargin(scale, "meters"));
    return new ee.Feature(pt).set('offset', offset);
  }
  
  var lines = lineString.cutLines(distances, ee.ErrorMargin(scale, "meters")).geometries();
  var points =   lines.zip(distances).map(function(s) {
    var line = ee.List(s).get(0);
    var offset = ee.List(s).get(1);
    return makePointFeature(ee.Geometry(line).coordinates().get(0), offset);
  });
  points = points.add(makePointFeature(lineString.coordinates().get(-1), length));

  return new ee.FeatureCollection(points);
}
//*--------------------------------*//